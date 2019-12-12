#!/usr/bin/env Rscript

library(tidyverse)
library(optparse)

opt_list <- list(
  make_option(c("-i", "--input")),
  make_option(c("-o", "--output"))
)

opt <- parse_args(OptionParser(option_list = opt_list))

last_tax <- function(x) {
  if_else(
    x == "Unassigned",
    "Unassigned",
    Find(function(x) str_detect(x, "[kpcofg]__\\S+"), str_split(x, "; ")[[1]], right = TRUE)
  )
}

read_tab <- function(fn) {
  read_tsv(fn, col_names = FALSE) %>%
    select(feature = X1, tax = X2)
}

read_tab(opt$input) %>%
  mutate(tax = map_chr(tax, last_tax)) %>%
  write_tsv(opt$output, col_names = FALSE)
