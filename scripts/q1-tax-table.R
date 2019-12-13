#!/usr/bin/env Rscript --vanilla

library(tidyverse)
library(optparse)

opt_list <- list(
  make_option("--otu_map"),
  make_option("--tax", help = "Taxonomy map"),
  make_option("--table", help = "Feature table"),
  make_option(c("-o", "--output"))
)

opt <- parse_args(OptionParser(option_list = opt_list))

first_line <- function(fn) {
  read_lines(fn, n_max = 1) %>%
    first() %>%
    trimws()
}

read_table <- function(fn) {
  stopifnot(first_line(fn) == '# Constructed from biom file')

  read_tsv(fn, skip = 1) %>%
    rename(feature = `#OTU ID`) %>%
    mutate_at("feature", as.character)
}

otu_map <- melt_tsv(opt$otu_map) %>%
  group_by(row) %>%
  mutate(otu = value[col == 1]) %>%
  ungroup() %>%
  filter(col > 1) %>%
  select(feature = value, otu)

tax <- read_tsv(opt$tax, col_names = c("otu", "short_tax")) %>%
  mutate_at("otu", as.character)

table <- read_table(opt$table)

table %>%
  gather("sample", "counts", -feature) %>%
  filter(counts > 0) %>%
  left_join(otu_map, by = "feature") %>%
  left_join(tax, by = "otu") %>%
  replace_na(list(short_tax = "no_taxonomy")) %>%
  group_by(sample, short_tax) %>%
  summarize_at("counts", sum) %>%
  ungroup() %>%
  write_tsv(opt$output)
