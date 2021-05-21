library(tidyverse)


euk_bcs <-
  read_tsv(snakemake@input[[1]])


bact_bcs <-
  read_tsv(snakemake@input[[2]])


euk_hybrid_bcs <-
  euk_bcs %>%
  filter(str_detect(Sample, "Mock")) %>%
  mutate(Mock = str_detect(Taxonomy, "Mock")) %>%
  group_by(BC) %>%
  mutate(Has_mock = any(Mock)) %>%
  filter(!Mock && Has_mock) %>%
  arrange(BC) %>%
  select(-Mock, -Has_mock)


bact_hybrid_bcs <-
  bact_bcs %>%
  filter(str_detect(Sample, "Mock")) %>%
  mutate(Mock = str_detect(Taxonomy, "Mock")) %>%
  group_by(BC) %>%
  mutate(Has_mock = any(Mock)) %>%
  filter(!Mock && Has_mock) %>%
  arrange(BC) %>%
  select(-Mock, -Has_mock)


bind_rows(euk_hybrid_bcs, bact_hybrid_bcs) %>%
  write_tsv(snakemake@output[[1]])




