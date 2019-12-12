library(tidyverse)

data <- tibble(
  fn = snakemake@input,
  fields = map(fn, ~ str_split_fixed(., '-', 4)),
  direction = map_chr(fields, ~ .[[1]]),
  pick = map_chr(fields, ~ .[[2]]),
  data = map(fn, read_tsv)
) %>%
  select(direction, pick, data) %>%
  unnest(cols = 'data') %>%
  separate(sample, c('direction2', 'trim2', 'sample'))

with(data, {
  stopifnot(all(direction == direction2))
  stopifnot(all(trim2 == "trim"))
})

# discarded <- data %>%
#   group_by(direction, sample) %>%
#   summarize_at('counts', sum) %>%
#   mutate(
#     short_tax = 'discarded',
#     counts = max(counts) - counts
#   ) %>%
#   ungroup() %>%
#   filter(counts > 0)

plot_data <- data %>%
  # bind_rows(discarded) %>%
  # group up relative abundances
  group_by(direction, pick, sample) %>%
  mutate(
    ra = counts / sum(counts),
    short_tax = if_else(ra <= 0.05, 'Other', short_tax)
  ) %>%
  group_by(direction, pick, sample, short_tax) %>%
  summarize_at('ra', sum) %>%
  ungroup()

plot <- plot_data %>%
  filter(sample %in% c("heli5", "noheli5")) %>%
  ggplot(aes(interaction(direction, pick), ra, fill = fct_reorder(short_tax, -ra))) +
  facet_wrap(~ sample) +
  geom_col() +
  coord_flip() +
  theme_minimal()

ggsave(snakemake@output[["plot"]])
write_tsv(plot_data, snakemake@output[["data"]])
