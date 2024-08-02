#!/usr/bin/env Rscript

if (exists("snakemake")) {
  input_files <- snakemake@input[["files"]]
  station <- snakemake@wildcards[["efas_station"]]
  output_filename <- snakemake@output[[1]]
  snakemake@source("extract.R")
} else {
  input_files <- list.files("results/preprocessing/EFAS", pattern = "[0-9]{8}.parquet", full.names = TRUE)
  # input_files <- list.files("results/preprocessing/HadUK", pattern = "[a-z]+_[0-9]{8}-[0-9]{8}.parquet", full.names = TRUE)
  station <- "10002"
  output_filename <- "results/preprocessing/HadUK/10002.parquet"
  source("workflow/scripts/extract.R")
}

## Open files as a dataset and filter on ID
ds <- open_dataset(input_files, unify_schemas = FALSE)
ds <- ds %>% dplyr::filter(id %in% station)
ds <- ds %>% collect()
ds <- ds %>% arrange(member, init_time, time)

## Write output
write_parquet(ds, output_filename)
