
library(tidyverse)


euk_bc <-
  read_tsv(snakemake@input[[2]]) %>%
  separate(Sample, c("Sample", "Junk"), sep="_") %>%
  mutate(Sample = str_replace(Sample, "18S", "")) %>%
  select(-Junk, -Count)


bact_bc <-
  read_tsv(snakemake@input[[1]]) %>%
  separate(Sample, c("Sample", "Junk"), sep="_") %>%
  mutate(Sample = str_replace(Sample, "16S", "")) %>%
  select(-Junk, -Count)


joined_bcs <-
  inner_join(euk_bc, bact_bc, by=c("Sample", "BC")) %>%
  count(Sample, Taxonomy.x, Taxonomy.y) %>%
  arrange(Sample, desc(n)) %>%
  setNames(c("Sample", "Eukaryotic_taxonomy", "Bacterial_taxonomy", "Count"))


joined_bcs %>%
  write_tsv(snakemake@output[[1]])


joined_bcs %>%
  filter(!str_detect(Eukaryotic_taxonomy, "d:Mock"),
         !str_detect(Bacterial_taxonomy, "d:Mock")) %>%
  write_tsv(snakemake@output[[2]])


joined_bcs %>%
  filter(!str_detect(Eukaryotic_taxonomy, "d:Mock"),
         !str_detect(Bacterial_taxonomy, "d:Mock")) %>%
  filter(str_detect(Bacterial_taxonomy, "loroplas")) %>%
  count(Sample, Eukaryotic_taxonomy) %>%
  arrange(Sample, desc(n)) %>%
  write_tsv(snakemake@output[[1]])
