library(tidyverse)

remove_tails <- function(x) {
  for (i in 1:7) {
    x <- str_replace(x, '; [kpcofgs]__$', '')
  }
  x
}

first_line <- function(fn) {
  read_lines(fn, n_max = 1) %>%
    first() %>%
    trimws()
}

read_tax <- function(fn) {
  line <- first_line(fn)

  if (line == 'Feature ID\tTaxon\tConfidence') {
    tax <- read_tsv(fn) %>%
      rename(otu = `Feature ID`) %>%
      mutate_at('otu', as.character)
  } else if (line == '367523\tk__Bacteria; p__Bacteroidetes; c__Flavobacteriia; o__Flavobacteriales; f__Flavobacteriaceae; g__Flavobacterium; s__') {
    tax <- read_tsv(fn, col_names = c('otu', 'Taxon')) %>%
      mutate_at('Taxon', remove_tails) %>%
      mutate_at('otu', as.character)
  } else {
    stop('Unrecognized file format')
  }
}

read_table <- function(fn) {
  stopifnot(first_line(fn) == '# Constructed from biom file')

  read_tsv(fn, skip = 1) %>%
    rename(otu = `#OTU ID`) %>%
    mutate_at('otu', as.character)
}

tax <- read_tax(snakemake@input[['tax']])
table <- read_table(snakemake@input[['table']])

table %>%
  gather('sample', 'counts', -otu) %>%
  left_join(tax, by = 'otu') %>%
  replace_na(list(Taxon = 'no_taxonomy_found')) %>%
  write_tsv(snakemake@output[[1]])
