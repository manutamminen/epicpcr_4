
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


bact_tip_coords <-
  tip_coord(bact_tre) %>%
  separate(Tip, c("Taxonomy", "Abundance"), sep="____") %>%
  select(-Abundance)


euk_tip_coords <-
  tip_coord(euk_tre) %>%
  separate(Tip, c("Taxonomy", "Abundance"), sep="____") %>%
  select(-Abundance)


connections <-
    read_tsv(snakemake@input[[3]]) %>%
    mutate(Eukaryotic_taxonomy = str_replace_all(Eukaryotic_taxonomy, ":", "_"),
           Eukaryotic_taxonomy = str_replace_all(Eukaryotic_taxonomy, "\\(", "__"),
           Eukaryotic_taxonomy = str_replace_all(Eukaryotic_taxonomy, "\\),", "__"),
           Eukaryotic_taxonomy = str_replace_all(Eukaryotic_taxonomy, "\\)", ""),
           Eukaryotic_taxonomy = str_replace_all(Eukaryotic_taxonomy, " ", "_"),
           Bacterial_taxonomy = str_replace_all(Bacterial_taxonomy, ":", "_"),
           Bacterial_taxonomy = str_replace_all(Bacterial_taxonomy, "\\(", "__"),
           Bacterial_taxonomy = str_replace_all(Bacterial_taxonomy, "\\),", "__"),
           Bacterial_taxonomy = str_replace_all(Bacterial_taxonomy, "\\)", ""),
           Bacterial_taxonomy = str_replace_all(Bacterial_taxonomy, " ", "_")) %>%
    inner_join(bact_tip_coords, by = c("Bacterial_taxonomy" = "Taxonomy")) %>%
    inner_join(euk_tip_coords, by = c("Eukaryotic_taxonomy" = "Taxonomy")) %>%
    group_by(Sample) %>%
    mutate(Perc = Count/sum(Count))


bact_tip_percs <-
  read_tsv(snakemake@input[[4]]) %>%
  mutate(Taxonomy = str_replace_all(Taxonomy, ":", "_"),
         Taxonomy = str_replace_all(Taxonomy, "\\(", "__"),
         Taxonomy = str_replace_all(Taxonomy, "\\),", "__"),
         Taxonomy = str_replace_all(Taxonomy, "\\)", ""),
         Taxonomy = str_replace_all(Taxonomy, " ", "_")) %>%
  inner_join(bact_tip_coords, by="Taxonomy") %>%
  group_by(Sample) %>%
  mutate(Perc = Count / sum(Count))

bact_tip_labels <-
  tip_coord(bact_tre) %>%
  separate(Tip, c("Taxonomy", "Abundance"), sep="____") %>%
  mutate(Abundance = as.numeric(Abundance)) %>%
  arrange(desc(Abundance))


euk_tip_percs <-
  read_tsv(snakemake@input[[5]]) %>%
  mutate(Taxonomy = str_replace_all(Taxonomy, ":", "_"),
         Taxonomy = str_replace_all(Taxonomy, "\\(", "__"),
         Taxonomy = str_replace_all(Taxonomy, "\\),", "__"),
         Taxonomy = str_replace_all(Taxonomy, "\\)", ""),
         Taxonomy = str_replace_all(Taxonomy, " ", "_")) %>%
  inner_join(euk_tip_coords, by="Taxonomy") %>%
  group_by(Sample) %>%
  mutate(Perc = Count / sum(Count))


euk_tip_labels <-
  tip_coord(euk_tre) %>%
  separate(Tip, c("Taxonomy", "Abundance"), sep="____") %>%
  mutate(Abundance = as.numeric(Abundance)) %>%
  arrange(desc(Abundance))


plot_bact_hist <- function(sample, color) {
  plot(NULL, xlim = c(0, 1e4), ylim = c(0, 1), type='n', axes=FALSE, ann=FALSE)
  bact_tip_percs %>%
    filter(Sample == sample) %>%
    with(walk2(Ix, Count, ~ lines(c(0, .y), c(.x, .x), col=color)))
}


plot_euk_hist <- function(sample, color) {
  plot(NULL, xlim = c(-1e4, 0), ylim = c(0, 1), type='n', axes=FALSE, ann=FALSE)
  euk_tip_percs %>%
    filter(Sample == sample) %>%
    with(walk2(Ix, Count, ~ lines(c(1, (1 - .y)), c(.x, .x), col=color)))
}


plot_connections <- function(sample, color) {
  connections %>%
    filter(Sample == sample) %>%
    with(pwalk(list(Ix.x, Ix.y),
               ~ lines(c(0, 1), c(.x, .y), col = alpha(color, 0.05))))
}


plot_tanglegram <- function(samples,
                            color_list,
                            normalize_connections=TRUE,
                            n_labels) {
  sample_table <-
    data.frame(Sample = samples,
               Color = color_list)

  if (!normalize_connections)
    connections$Perc <- 0.05

  widths <- c(5, 2, rep(1, length(samples)), 3, rep(1, length(samples)), 2, 5)
  lmat <- matrix(1:length(widths), ncol = length(widths))
  layout(lmat, widths = widths, heights = 1)
  par(mar=c(1, 1, 1, 1))

  plot(NULL, xlim = c(0, 5e3),
       ylim = c(0, 1), type='n', axes=FALSE, ann=FALSE)

  bact_tip_labels %>%
    head(n_labels) %>%
    with(walk2(Ix, Taxonomy, ~ text(5e3, .x, .y, cex=0.4, adj=1)))

  plot(ladderize(bact_tre), cex = 0.4,
       align.tip.label = TRUE, show.tip.label = FALSE)

  sample_table %>%
    with(walk2(Sample, Color, ~plot_bact_hist(.x, .y)))

  plot(NULL, xlim = c(0, 1), ylim = c(0, 1), type='n', axes=FALSE, ann=FALSE)

  connections %>%
    filter(Sample %in% samples) %>%
    with(pwalk(list(Ix.x, Ix.y, Sample, Perc),
               function(Ixx, Ixy, Sample, Perc)
                 lines(c(0, 1), c(Ixx, Ixy),
                       col = alpha(sample_table[sample_table == Sample, 2], Perc))))

  sample_table %>%
    map_dfr(rev) %>%
    with(walk2(Sample, Color, ~plot_euk_hist(.x, .y)))

  plot(ladderize(euk_tre), align.tip.label = TRUE,
       direction = "leftwards", show.tip.label = FALSE)

  plot(NULL, xlim = c(0, 5e3),
       ylim = c(0, 1), type='n', axes=FALSE, ann=FALSE)

  euk_tip_labels %>%
    head(n_labels) %>%
    with(walk2(Ix, Taxonomy, ~ text(0, .x, .y, cex=0.4, adj=0)))
}


for (fname_ix in seq_along(snakemake@output)) {
  png(snakemake@output[[fname_ix]], units="in", width=10, height=5, res=300)
  plot_tanglegram(snakemake@params[[fname_ix]]$Samples,
                  snakemake@params[[fname_ix]]$Colors,
                  snakemake@params[[fname_ix]]$Normalize_connections,
                  snakemake@params[[fname_ix]]$N_Labels)
  dev.off()
}

