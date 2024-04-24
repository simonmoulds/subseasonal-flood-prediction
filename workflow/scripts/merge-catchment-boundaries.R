#!/usr/bin/env Rscript

library(tidyverse)
library(lubridate)
library(arrow)
library(sf)

sf::sf_use_s2(FALSE)
options(dplyr.summarise.inform = FALSE)

if (exists("snakemake")) {
  input_filenames <- snakemake@input
  output_filename <- snakemake@output[[1]]
} else {
  input_filenames <- list.files("results/preprocessing/streamflow/catchment_boundaries", pattern = "[0-9]+.gpkg", full.names = TRUE)
  output_filename <- "results/preprocessing/streamflow/catchment_boundaries/merged.gpkg"
}

merged <- dplyr::bind_rows(lapply(input_filenames, FUN=function(x) st_read(x)))
st_write(merged, output_filename, quiet = TRUE, delete_layer = TRUE)
