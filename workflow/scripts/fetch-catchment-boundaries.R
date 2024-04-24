#!/usr/bin/env Rscript

library(tidyverse)
library(lubridate)
library(arrow)
library(sf)

sf::sf_use_s2(FALSE)
options(dplyr.summarise.inform = FALSE)

if (exists("snakemake")) {
  station <- snakemake@wildcards[["station"]]
  output_filename <- snakemake@output[[1]]
} else {
  ## TESTING
  station <- 24001
  output_filename <- "results/preprocessing/streamflow/catchment_boundaries/24001.gpkg"
}

nrfa_basins_filename <- "data/NRFA_shp/nrfa_all_catchments.shp"
nrfa_catchment_poly <- st_read(nrfa_basins_filename) %>%
  mutate(
    ID = as.integer(ID),
    ID_STRING = as.character(ID),
    SOURCE = "NRFA",
    VERSION = NA,
    EXPORTED = NA
  ) %>%
  dplyr::select(ID_STRING, ID, SOURCE, VERSION, EXPORTED, geometry)
nrfa_catchment_poly <- st_transform(nrfa_catchment_poly, 4326)
poly <- nrfa_catchment_poly %>% filter(ID %in% station)
st_write(poly, output_filename, quiet = TRUE, delete_layer = TRUE)


