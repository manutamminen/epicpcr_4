
library(tidyverse)


bact_bc <- read_tsv(snakemake@input[[1]])
euk_bc <- read_tsv(snakemake@input[[2]])


bact_bc_counts <-
  bact_bc %>%
  group_by(Sample, BC) %>%
  summarise(n = sum(Count)) %>%
  arrange(Sample, desc(n)) %>%
  group_by(Sample) %>%
  mutate(id = row_number()) %>%
  separate(Sample, c("Sample", "Junk"), sep="_") %>%
  select(-Junk) %>%
  mutate(Ribo = "16S",
         Sample = str_replace(Sample, "16S", ""))


euk_bc_counts <-
  euk_bc %>%
  group_by(Sample, BC) %>%
  summarise(n = sum(Count)) %>%
  arrange(Sample, desc(n)) %>%
  group_by(Sample) %>%
  mutate(id = row_number()) %>%
  separate(Sample, c("Sample", "Junk"), sep="_") %>%
  select(-Junk) %>%
  mutate(Ribo = "18S",
         Sample = str_replace(Sample, "18S", ""))


bc_counts <-
  bind_rows(bact_bc_counts, euk_bc_counts)


png(snakemake@output[[1]], units="in", width=5, height=5, res=300)
bc_counts %>%
  ggplot(aes(x = id, y = n)) +
  geom_density(stat = "identity") +
  facet_grid(Sample ~ Ribo, scales = "free") +
  scale_x_log10() +
  theme(strip.text.y = element_text(angle = 0),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size = 5))
dev.off()
# ggsave(snakemake@output[[1]], last_plot())


bact_bc_taxa <-
  bact_bc %>%
  select(-Count) %>%
  count(Sample, BC) %>%
  arrange(Sample, desc(n)) %>%
  group_by(Sample) %>%
  mutate(id = row_number()) %>%
  separate(Sample, c("Sample", "Junk"), sep="_") %>%
  select(-Junk) %>%
  mutate(Ribo = "16S",
         Sample = str_replace(Sample, "16S", ""))


euk_bc_taxa <-
  euk_bc %>%
  select(-Count) %>%
  count(Sample, BC) %>%
  arrange(Sample, desc(n)) %>%
  group_by(Sample) %>%
  mutate(id = row_number()) %>%
  separate(Sample, c("Sample", "Junk"), sep="_") %>%
  select(-Junk) %>%
  mutate(Ribo = "18S",
         Sample = str_replace(Sample, "18S", ""))


bc_taxa <-
  bind_rows(bact_bc_taxa, euk_bc_taxa)

png(snakemake@output[[2]], units="in", width=5, height=5, res=300)
bc_taxa %>%
  ggplot(aes(x = id, y = n)) +
  geom_density(stat = "identity") +
  facet_grid(Sample ~ Ribo, scales = "free") +
  scale_x_log10() +
  theme(strip.text.y = element_text(angle = 0),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size = 5))
dev.off()
# ggsave(snakemake@output[[2]], last_plot())

