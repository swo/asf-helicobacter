library(tidyverse)

read_otu <- function(fn) read_tsv(fn, col_types = cols(otu = 'c'))

data <- tibble(
  fn = snakemake@input,
  fields = map(fn, ~ str_split_fixed(., '-', 4)),
  direction = map_chr(fields, ~ .[[1]]),
  trim = map_chr(fields, ~ .[[2]]),
  picking = map_chr(fields, ~ .[[3]]),
  data = map(fn, read_otu)
) %>%
  select(direction, trim, picking, data) %>%
  unnest(cols = 'data') %>%
  mutate_at('otu', as.character) %>%
  separate(sample, c('direction2', 'trim2', 'sample')) %>%
  separate(Taxon, c('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'), '; ', fill = 'right', remove = FALSE)

with(data, {
  stopifnot(all(direction == direction2))
  stopifnot(all(trim == trim2))
})

discarded <- data %>%
  group_by(direction, trim, sample, picking, sample) %>%
  summarize_at('counts', sum) %>%
  mutate(
    Taxon = 'Discarded',
    counts = max(counts) - counts
  ) %>%
  ungroup() %>%
  filter(counts > 0)

plot <- data %>%
  bind_rows(discarded) %>%
  group_by(direction, trim, picking, sample) %>%
  mutate(ra = counts / sum(counts)) %>%
  ungroup() %>%
  mutate(
    tax = case_when(
      Taxon == 'Discarded' ~ 'Discarded',
      ra <= 0.015 ~ 'Other',
      picking == 'openref' & is.na(Taxon) & str_length(otu) == 32 ~ 'New open-ref OTU',
      genus != 'g__' ~ genus,
      family != 'f__' ~ family,
      order != 'o__' ~ order,
      class != 'c__' ~ class,
      phylum != 'p__' ~ phylum,
      kingdom != 'k__' ~ kingdom
    )
  ) %>%
  group_by(direction, trim, picking, sample, tax) %>%
  summarize_at('ra', sum) %>%
  ungroup() %>%
  ggplot(aes(interaction(direction, trim, picking), ra, fill = fct_reorder(tax, -ra))) +
  facet_wrap(~ sample) +
  geom_col() +
  coord_flip() +
  theme_minimal()

ggsave(snakemake@output[[1]])
