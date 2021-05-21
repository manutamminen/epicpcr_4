
library(tidyverse)

png(snakemake@output[[2]], units="in", width=5, height=5, res=300)
read_tsv(snakemake@input[[1]]) %>%
  arrange(Sample, desc(Count)) %>%
  group_by(Sample) %>%
  mutate(Ix = row_number()) %>%
  ggplot(aes(x=Ix, y=Count)) +
  geom_density(stat = "identity") +
  facet_grid(Sample~., scales="free") +
  scale_x_log10() +
  theme(strip.text.y = element_text(angle = 0),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size = 5))
dev.off()
# ggsave(snakemake@output[[1]])
