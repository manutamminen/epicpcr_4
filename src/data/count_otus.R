
library(tidyverse)


read_tsv(snakemake@input[[1]]) %>%
  separate(Sample, c("Sample", "Junk"), sep="_") %>%
  mutate(Sample = str_replace(Sample, "16S", "")) %>%
  filter(!str_detect(Taxonomy, "Mock")) %>%
  select(-Junk, -Count) %>%
  count(Sample, Taxonomy) %>%
  setNames(c("Sample", "Taxonomy", "Count")) %>%
  arrange(Sample, desc(Count)) %>%
  write_tsv(snakemake@output[[1]])


read_tsv(snakemake@input[[2]]) %>%
  separate(Sample, c("Sample", "Junk"), sep="_") %>%
  mutate(Sample = str_replace(Sample, "18S", "")) %>%
  filter(!str_detect(Taxonomy, "Mock")) %>%
  select(-Junk, -Count) %>%
  count(Sample, Taxonomy) %>%
  setNames(c("Sample", "Taxonomy", "Count")) %>%
  arrange(Sample, desc(Count)) %>%
  write_tsv(snakemake@output[[2]])

