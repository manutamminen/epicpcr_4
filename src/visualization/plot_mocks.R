
library(tidyverse)

bact_bcs <-
  read_tsv(snakemake@input[[1]])


euk_bcs <-
  read_tsv(snakemake@input[[2]])

png(snakemake@output[[1]], units="in", width=5, height=5, res=300)
bact_bcs %>%
  count(Sample, Taxonomy) %>%
  filter(str_detect(Taxonomy, "Mock")) %>%
  arrange(Sample, desc(n)) %>%
  separate(Taxonomy, into=c(NA, "Ix"), sep="ock", remove=FALSE) %>%
  separate(Ix, into=c("Ix", NA), sep="_") %>%
  separate(Sample, into=c("Sample", NA), sep="_") %>%
  mutate(Ix = as.numeric(Ix)) %>%
  mutate(Fus = ifelse(Ix < 4, "No_fusion", "Fusion")) %>%
  select(-Taxonomy) %>%
  group_by(Sample, Ix, Fus) %>%
  summarise(n = sum(n)) %>%
  ggplot(aes(x=Ix, y=n)) +
  geom_point() +
  geom_smooth(method = "glm",
              family = gaussian(link="log")) +
  facet_grid(Sample ~ Fus, scales="free") +
  theme(strip.text.y = element_text(angle = 0),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size = 5))
dev.off()


png(snakemake@output[[2]], units="in", width=5, height=5, res=300)
euk_bcs %>%
  count(Sample, Taxonomy) %>%
  filter(str_detect(Taxonomy, "Mock")) %>%
  arrange(Sample, desc(n)) %>%
  separate(Taxonomy, into=c(NA, "Ix"), sep="ock", remove=FALSE) %>%
  separate(Ix, into=c("Ix", NA), sep="_") %>%
  separate(Sample, into=c("Sample", NA), sep="_") %>%
  mutate(Ix = as.numeric(Ix)) %>%
  mutate(Fus = ifelse(Ix < 4, "No_fusion", "Fusion")) %>%
  select(-Taxonomy) %>%
  group_by(Sample, Ix, Fus) %>%
  summarise(n = sum(n)) %>%
  ggplot(aes(x=Ix, y=n)) +
  geom_point() +
  geom_smooth(method = "glm",
              family = gaussian(link="log")) +
  facet_grid(Sample ~ Fus, scales="free") +
  theme(strip.text.y = element_text(angle = 0),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size = 5))
dev.off()

