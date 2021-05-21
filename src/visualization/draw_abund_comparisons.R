
library(tidyverse)
library(ape)


bact_abund_tips <-
  read.tree(snakemake@input[[1]]) %>%
  .$tip.label %>%
  tibble %>%
  separate(".", c("Taxonomy", "Abundance"), sep="____", remove=FALSE) %>%
  mutate(Abundance = as.numeric(Abundance),
         Abundance = ifelse(is.na(Abundance), 51, Abundance)) %>%
  filter(Abundance > 50) %>%
  pull(".")


bact_tre <-
  read.tree(snakemake@input[[1]]) %>%
  keep.tip(bact_abund_tips) %>%
  root(outgroup = "JQ837894.1.1415",
       resolve.root = TRUE)


euk_abund_tips <-
  read.tree(snakemake@input[[2]]) %>%
  .$tip.label %>%
  tibble %>%
  separate(".", c("Taxonomy", "Abundance"), sep="____", remove=FALSE) %>%
  mutate(Abundance = as.numeric(Abundance),
         Abundance = ifelse(is.na(Abundance), 51, Abundance)) %>%
  filter(Abundance > 50) %>%
  pull(".")


euk_tre <-
  read.tree(snakemake@input[[2]]) %>%
    keep.tip(euk_abund_tips) %>%
    root(outgroup = "Human_18S_rRNA_gene",
         resolve.root = TRUE)


is_tip <- function(tree) tree$edge[,2] <= length(tree$tip.label)
ordered_tips <- function(tree) tree$edge[is_tip(tree), 2]
tip_coord <- function(tree) {
    tre_len <- length(tree$tip.label)
    tibble(Tip = tree$tip.label[ordered_tips(tree)]) %>%
        mutate(Ix = seq(0, 1, 1 / (tre_len - 1)))
}


plot_hist <- function(tip_percs, sample, color) {
  max_val <- max(tip_percs$Count)
  plot(NULL, xlim = c(0, max_val),
       ylim = c(0, 1),
       type='n',
       axes=FALSE,
       ann=FALSE)

  tip_percs %>%
    filter(Sample == sample) %>%
    with(walk2(Ix, Count, ~ lines(c(0, .y), c(.x, .x), col=color)))
}


plot_histograms <- function(type,
                            samples,
                            colors,
                            no_labels=10) {

  if (type == "bacteria") {
    tree <- bact_tre
    abunds <- read_tsv(snakemake@input[[3]])
  } else {
    tree <- euk_tre
    abunds <- read_tsv(snakemake@input[[4]])
  }

  tip_coords <-
    tip_coord(tree) %>%
    separate(Tip, c("Taxonomy", "Abundance"), sep="____")

  tip_labels <-
    tip_coord(tree) %>%
    separate(Tip, c("Taxonomy", "Abundance"), sep="____") %>%
    mutate(Abundance = as.numeric(Abundance)) %>%
    arrange(desc(Abundance)) %>%
    head(no_labels)

  tip_percs <-
    abunds %>%
    mutate(Taxonomy = str_replace_all(Taxonomy, ":", "_"),
           Taxonomy = str_replace_all(Taxonomy, "\\(", "__"),
           Taxonomy = str_replace_all(Taxonomy, "\\),", "__"),
           Taxonomy = str_replace_all(Taxonomy, "\\)", ""),
           Taxonomy = str_replace_all(Taxonomy, " ", "_")) %>%
    inner_join(tip_coords, by="Taxonomy") %>%
    group_by(Sample) %>%
    mutate(Perc = Count / sum(Count))

  sample_table <-
    data.frame(Sample = samples,
               Color = colors)

  widths <- c(2, rep(1, length(samples)), 8)
  lmat <- matrix(1:length(widths), ncol = length(widths))
  layout(lmat, widths = widths, heights = 1)
  par(mar=c(1, 1, 1, 1))

  plot(ladderize(tree), cex = 0.4,
       align.tip.label = TRUE, show.tip.label = FALSE)

  sample_table %>%
    with(walk2(Sample, Color,
               ~plot_hist(tip_percs=tip_percs, sample=.x, color=.y)))

  plot(NULL, xlim = c(0, 5e3),
       ylim = c(0, 1), type='n', axes=FALSE, ann=FALSE)

  tip_labels %>%
    with(walk2(Ix, Taxonomy, ~ text(0, .x, .y, cex=0.6, adj=0)))
}


for (fname_ix in seq_along(snakemake@output)) {
  png(snakemake@output[[fname_ix]], units="in", width=8, height=5, res=300)
  plot_histograms(snakemake@params[[fname_ix]]$Type,
                  snakemake@params[[fname_ix]]$Samples,
                  snakemake@params[[fname_ix]]$Colors,
                  snakemake@params[[fname_ix]]$N_Labels)
  dev.off()
}

