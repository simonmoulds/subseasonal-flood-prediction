#!/usr/bin/env Rscript

if (exists("snakemake")) {
  config <- snakemake@config
  c3s_input_file <- snakemake@input[[1]]
  poly_input_file <- snakemake@input[[2]]
  output_file <- snakemake@output[[1]]
  variable <- snakemake@wildcards[["c3s_variable"]]
  snakemake@source("extract.R")
  snakemake@source("utils.R")
} else {
  config <- read_yaml("config/config.yml")
  poly_input_file <- "results/preprocessing/streamflow/catchment_boundaries/merged.gpkg"
  output_file <- ""
  source("workflow/scripts/extract.R")
  source("workflow/scripts/utils.R")
}

catchment_poly <- st_read(poly_input_file, quiet = TRUE)

x <- extract_c3s_daily_catchment_data(
  c3s_input_file, variable, catchment_poly,
  "ID_STRING", progress = FALSE
)
x <- x %>% rename(system = source_id)

## Restrict members
selected_members <- config$hindcast$members |> parse_range()
x <- x %>% filter(member %in% selected_members)

## ## Assign consistent variable names
## x <- x %>% rename_variables()
## ## Add product type to variable name
## x <- x %>% unite("variable", "variable", "product_type")

## Units
if (str_detect(variable, "tp")) {
  ## Convert to mm
  x <- x %>% mutate(value = value * 1000)
}

write_parquet(x, output_file)
