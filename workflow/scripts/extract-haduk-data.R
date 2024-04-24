#!/usr/bin/env Rscript

if (exists("snakemake")) {
  config <- snakemake@config
  haduk_input_file <- snakemake@input[[1]]
  poly_input_file <- snakemake@input[[2]]
  output_file <- snakemake@output[[1]]
  variable <- snakemake@wildcards[["haduk_variable"]]
  snakemake@source("extract.R")
  snakemake@source("utils.R")
} else {
  config <- read_yaml("config/config.yml")
  station <- 28083
  variable <- "rainfall"
  output_directory <- 'results/input/haduk-field'
  streamflow_directory <- 'results/input/streamflow'
  source("workflow/scripts/extract.R")
  source("workflow/scripts/utils.R")
}

catchment_poly <- st_read(poly_input_file, quiet = TRUE)

x <- extract_haduk_daily_catchment_data(
  haduk_input_file, variable, catchment_poly,
  "ID_STRING", progress = FALSE
)

## ## Assign consistent variable names
## x <- x %>% rename_variables()

write_parquet(x, output_file)
