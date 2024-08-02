#!/usr/bin/env Rscript

## if (exists("snakemake")) {
##   config <- snakemake@config
##   c3s_input_file <- snakemake@input[[1]]
##   poly_input_file <- snakemake@input[[2]]
##   output_file <- snakemake@output[[1]]
##   variable <- snakemake@wildcards[["c3s_variable"]]
##   snakemake@source("extract.R")
##   snakemake@source("utils.R")
## } else {
##   config <- read_yaml("config/config.yml")
##   poly_input_file <- "results/preprocessing/streamflow/catchment_boundaries/merged.gpkg"
##   output_file <- ""
##   source("workflow/scripts/extract.R")
##   source("workflow/scripts/utils.R")
## }

library(argparse)

source("workflow/scripts/extract.R")
source("workflow/scripts/utils.R")

parser <- ArgumentParser(description = "Extract C3S data")

parser$add_argument("--var", help = "Variable name")
parser$add_argument("--poly", help = "Polygons")
parser$add_argument("--input", help = "Input file")
parser$add_argument("--output", help = "Output file")

xargs <- parser$parse_args()

variable <- xargs$var
polygons <- xargs$poly
input_file <- xargs$input
output_file <- xargs$output

catchment_poly <- st_read(polygons, quiet = TRUE)

# input_file <- 'data/C3S_hindcast_daily_SEAS5/ecmf_SEAS5-v20171101_hindcast_S19840801_tp_daily.nc'
# input_file <- 'results/preprocessing/C3S/precip_disagg/ecmf_SEAS5-v20171101_hindcast_S19840801_tp_daily_disagg.nc'
# catchment_poly <- st_read("results/preprocessing/streamflow/catchment_boundaries/merged.gpkg")
x <- extract_c3s_daily_catchment_data(
  input_file, variable, catchment_poly,
  "ID_STRING", progress = FALSE
)
x <- x %>% rename(system = source_id)

## ## Restrict members
## selected_members <- config$hindcast$members |> parse_range()
## x <- x %>% filter(member %in% selected_members)

## ## Assign consistent variable names
## x <- x %>% rename_variables()
## ## Add product type to variable name
## x <- x %>% unite("variable", "variable", "product_type")

## Units
if (str_detect(variable, "tp")) {
  ## Convert to mm
  x <- x |> mutate(value = value * 1000)
}

# We focus on subseasonal timescales, so no need to keep all 215 days 
x <- x |> filter(lead_time <= 90)

write_parquet(x, output_file)

# ## TESTING
# library(terra)
# r <- terra::rast('data/C3S_hindcast_daily_SEAS5/ecmf_SEAS5-v20171101_hindcast_S19840801_tp_daily.nc')
# ix <- grep("number=10", names(r))
# r <- r[[ix]]
# r <- rotate(r)
# v <- terra::extract(r, data.frame(x=-1.83, y=57.5))

# ## ## r <- terra::rast('data/C3S_hindcast_daily_SEAS5/ecmf_SEAS5-v20171101_1993_8_daily.nc', 2)
# ## ## ix <- grep("number=0", names(r))
# ## ## r <- r[[ix]]
# ## ## r <- rotate(r)
# ## ## v0 <- terra::extract(r, data.frame(x=-1.83, y=57.5)) |> as.numeric()

# r <- terra::rast("results/preprocessing/C3S/precip_disagg/ecmf_SEAS5-v20171101_hindcast_S19840801_tp_daily_disagg.nc")
# ix <- grep("number=0", names(r))
# r <- r[[ix]]
# r <- rotate(r)
# v0 <- terra::extract(r, data.frame(x=-1.83, y=57.5))
