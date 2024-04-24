#!/usr/bin/env Rscript

library(dplyr)
library(arrow)

if (exists("snakemake")) {
  input_files <- snakemake@input[["files"]]
  station <- snakemake@wildcards[["station"]]
  output_filename <- snakemake@output[[1]]
} else {
  input_files <- list.files("results/preprocessing/C3S", pattern = "^SEAS5-v20171101_[0-9]{8}_.*.parquet", full.names = TRUE)
  station <- 10002
  output_filename <- "results/preprocessing/ERA5/10002.parquet"
}

## Open files as a dataset and filter on ID
ds <- open_dataset(input_files, unify_schemas = FALSE)
ds <- ds %>% filter(id %in% station)
ds <- ds %>% collect()

write_parquet(ds, output_filename)
