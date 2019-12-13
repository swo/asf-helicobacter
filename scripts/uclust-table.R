library(tidyverse)

first_line <- function(fn) {
  read_lines(fn, n_max = 1) %>%
    first() %>%
    trimws()
}

read_table <- function(fn) {
  stopifnot(first_line(fn) == '# Constructed from biom file')

  read_tsv(fn, skip = 1) %>%
    rename(otu = `#OTU ID`) %>%
    mutate_at('otu', as.character)
}

read_otu_assignments <- function(fn) {
  raw <- melt_tsv(fn)

  otus <- raw %>%
    filter(col == 1) %>%
    select(row, otu = value)

  seqs <- raw %>%
    filter(col != 1) %>%
    select(row, seq = value)

  left_join(otus, seqs, by = 'row') %>%
    select(-row)
}

sequence_table <- read_table(snakemake@input[['sequence_table']]) %>%
  rename(seq = otu)

otu_assignments <- read_otu_assignments(snakemake@input[['otu_assignments']])

output_table <- sequence_table %>%
  left_join(otu_assignments, by = 'seq') %>%
  select(-seq) %>%
  group_by(otu) %>%
  summarize_at(vars(-group_cols()), sum) %>%
  rename(`#OTU ID` = otu) %>%
  format_tsv()

output_table <- str_c('# Constructed from biom file', '\n', output_table, sep = '')

write_file(output_table, snakemake@output[[1]])
