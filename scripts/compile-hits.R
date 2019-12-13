library(tidyverse)

# Load blast6 hits ----------------------------------------------------

read_b6 <- function(fn) {
  read_tsv(fn, col_names = FALSE) %>%
    select(feature = X1, otu = X2, id = X3) %>%
    mutate_at("otu", as.character)
}

hits <- tibble(
  fn = snakemake@input[["hits"]],
  direction = str_match(fn, "^[a-z]+")[, 1],
  data = map(fn, read_b6)
) %>%
  unnest() %>%
  select(-fn)

# Load deblur feature tables ------------------------------------------

read_otu <- function(fn) {
  read_tsv(fn, skip = 1) %>%
    rename(feature = `#OTU ID`)
}

features <- tibble(
  fn = snakemake@input[["tables"]],
  direction = str_match(fn, "^[a-z]+")[, 1],
  data = map(fn, read_otu)
) %>%
  unnest() %>%
  select(-fn) %>%
  gather("sample", "counts", -direction, -feature) %>%
  filter(!is.na(counts)) %>%
  separate(sample, c("direction2", "trim", "sample"))

stopifnot(all(features$trim == "trim"))
stopifnot(all(features$direction == features$direction2))

features <- features %>% select(-direction2, -trim)

# Load taxonomies -----------------------------------------------------

tax <- read_tsv(snakemake@input[["tax"]], col_names = c("otu", "tax")) %>%
  mutate_at("otu", as.character)

# Combine tables ------------------------------------------------------

features %>%
  left_join(hits, by = c("direction", "feature")) %>%
  left_join(tax, by = "otu") %>%
  write_tsv(snakemake@output[[1]])
