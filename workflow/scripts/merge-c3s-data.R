#!/usr/bin/env Rscript

library(dplyr)
library(arrow)

if (exists("snakemake")) {
  config <- snakemake@config
  flux_input_files <- snakemake@input[["flux_files"]]
  inst_input_files <- snakemake@input[["inst_files"]]
  station <- snakemake@wildcards[["station"]]
  output_filename <- snakemake@output[[1]]
  snakemake@source("utils.R")
} else {
  config <- read_yaml("config/config.yml")
  flux_input_files <- list.files("results/preprocessing/C3S", pattern = "tp_.*.parquet", full.names = TRUE)
  inst_input_files <- list.files("results/preprocessing/C3S", pattern = "t2m_.*.parquet", full.names = TRUE)
  station <- "10002"
  output_filename <- "results/preprocessing/C3S/10002.parquet"
  source("workflow/scripts/utils.R")
}

input_files <- c(flux_input_files, inst_input_files)

## Open files as a dataset and filter on ID
ds <- open_dataset(input_files, unify_schemas = FALSE)

## Restrict members
selected_members <- config$hindcast$members |> parse_range()
ds <- ds %>% filter(member %in% selected_members)
ds <- ds %>% filter(id %in% station)
ds <- ds %>% collect()

write_parquet(ds, output_filename)
