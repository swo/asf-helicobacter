library(tidyverse)

short_tax <- function(x) {
  x %>%
    # remove species
    str_replace('s__\\S*$', '') %>%
    # get only non-empty levels
    str_extract_all('[kpcofg]__\\S*') %>%
    first() %>%
    last()
}

raw <- read_tsv(snakemake@input[[1]])

tax <- raw %>%
  select(Taxon) %>%
  distinct() %>%
  mutate(tax = map_chr(Taxon, short_tax))

plot_data <- raw %>%
  left_join(tax, by = 'Taxon') %>%
  group_by(sample) %>%
  mutate(ra = counts / sum(counts)) %>%
  mutate(
    tax = case_when(
      ra <= 0.015 ~ 'Other',
      is.na(tax) & str_length(otu) == 32 ~ 'New open-ref OTU',
      TRUE ~ tax
    )
  )

plot <- plot_data %>%
  ggplot(aes(sample, ra, fill = fct_reorder(tax, -ra))) +
  geom_col() +
  coord_flip()

ggsave(snakemake@output[[1]])
