source("~/.Rprofile")

library(tidyverse)

usearch <- function(input, output, accepts, rejects, id = 0.90) {
  system(str_glue("usearch -usearch_global {input} -db {snakemake@input[['otus']]} -id {id} -strand both -blast6out {output} -maxaccepts {accepts} -maxrejects {rejects}"))
}

read_b6 <- function(fn) {
  raw <- read_tsv(fn, col_names = FALSE)

  if (nrow(raw) == 0) {
    tibble(otu = character(0), id = numeric(0))
  } else {
    raw %>%
      select(otu = X2, id = X3) %>%
      mutate_at('otu', as.character)
  }
}

table <- read_tsv(snakemake@input[['table']], skip = 1) %>%
  rename(otu = `#OTU ID`) %>%
  mutate_at('otu', as.character)

top_otu <- table %>%
  gather('sample', 'counts', -otu) %>%
  group_by(otu) %>%
  summarize_at('counts', sum) %>%
  # keep only the "new" otus (with 32-character OTU IDs)
  filter(str_length(otu) == 32) %>%
  arrange(desc(counts)) %>%
  pull(otu) %>%
  first()

fasta <- read_lines(snakemake@input[['seqs']])
i <- Position(function(x) x == str_c('>', top_otu), fasta)
top_fasta_lines <- fasta[i: (i + 1)]

top_fn <- tempfile()
write_lines(top_fasta_lines, top_fn)

b6_fn <- tempfile()
usearch(top_fn, b6_fn, 0, 0)
exhaustive <- read_b6(b6_fn)

usearch(top_fn, b6_fn, 1, 8)
default <- read_b6(b6_fn)

unlink(top_fn)
unlink(b6_fn)

tax <- read_tsv(
    snakemake@input[['taxonomy']],
    col_names = c('otu', 'taxonomy')
  ) %>%
  mutate_at('otu', as.character)

bind_rows(
  default = default,
  exhaustive = exhaustive,
  .id = 'search_type'
) %>%
  left_join(tax, by = 'otu') %>%
  write_tsv(snakemake@output[[1]])
