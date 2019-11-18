library(tidyverse)

b6 <- read_tsv(snakemake@input[[1]], col_names = FALSE) %>%
  select(sequence = X1, otu = X2, identity = X3)

print(b6)

tax <- read_tsv(snakemake@params$taxonomy, col_names = c('otu', 'taxonomy'))

print(tax)

b6 %>%
  left_join(tax, by = 'otu') %>%
  write_tsv(snakemake@output[[1]])
