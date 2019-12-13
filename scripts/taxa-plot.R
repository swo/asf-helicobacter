library(tidyverse)
library(RColorBrewer)

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

discarded <- data %>%
  group_by(direction, sample, pick) %>%
  summarize_at("counts", sum) %>%
  mutate(
    short_tax = "discarded",
    counts = max(counts) - counts
  ) %>%
  ungroup() %>%
  filter(counts > 0)

plot_data <- data %>%
  bind_rows(discarded) %>%
  # group up relative abundances
  group_by(direction, pick, sample) %>%
  mutate(ra = counts / sum(counts)) %>%
  ungroup()

others <- plot_data %>%
  group_by(short_tax) %>%
  summarize(is_other = max(ra) < 0.10)

plot_data <- plot_data %>%
  left_join(others, by = "short_tax") %>%
  mutate(short_tax = if_else(is_other, "Other", short_tax)) %>%
  select(-is_other) %>%
  group_by(direction, pick, sample, short_tax) %>%
  summarize_at(c("counts", "ra"), sum) %>%
  ungroup() %>%
  mutate(short_tax = case_when(
    .$pick %in% c("q1closed", "q1open") & .$short_tax == "no_taxonomy" ~ "discarded",
    TRUE ~ .$short_tax
  ))

display_levels <- c(
  "g__Parabacteroides",
  "g__Mucispirillum",
  "discarded",
  "no_taxonomy",
  "g__Helicobacter",
  "g__Flexispira",
  "g__Turicibacter",
  "f__Lachnospiraceae",
  "f__Bacillaceae",
  "f__[Mogibacteriaceae]",
  "Other"
)

tax_names <- Filter(function(x) str_detect(x, "__"), display_levels)
colors <- c(
  setNames(brewer.pal(length(tax_names), "Paired"), tax_names),
  "discarded" = "black", "no_taxonomy" = "gray50", "Other" = "white"
)


plot <- plot_data %>%
  filter(sample %in% c("heli5", "noheli5")) %>%
  mutate(short_tax = factor(short_tax, levels = display_levels)) %>%
  ggplot(aes(interaction(direction, pick), ra, fill = short_tax)) +
  facet_wrap(~ sample) +
  geom_col(color = "black") +
  coord_flip() +
  scale_fill_manual(values = colors) +
  theme_minimal()

ggsave(snakemake@output[["plot"]])
write_tsv(plot_data, snakemake@output[["data"]])
