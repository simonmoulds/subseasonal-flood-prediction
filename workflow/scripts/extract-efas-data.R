#!/usr/bin/env Rscript

library(tidyverse)
library(arrow)
library(terra)
library(sf)
library(ncdf4)

if (exists("snakemake")) {
  sites_filename <- snakemake@input[[1]]
  input_filename <- snakemake@input[[2]]
  output_filename <- snakemake@output[[1]]
  init_time <- snakemake@wildcards[["efas_time"]]
  snakemake@source("utils.R")
} else {
  sites_filename <- "resources/efas_sites.parquet"
  input_filename <- 'data/efas_stations/ftp.ecmwf.int/for_louise_uni_oxford/dis_stations_SR1991010100.nc'
  output_filename <- 'results/preprocessing/EFAS/SR1991010100.parquet'
  init_time <- "19910101"
  source("workflow/scripts/utils.R")
}

sites <- read_parquet(sites_filename)
efas_ids <- sites$efas_id
nrfa_ids <- sites$nrfa_id

catchment_area <- sites |>
  dplyr::select(nrfa_id, efas_catchment_area_ldd) |>
  rename(id = nrfa_id)

init_time <- as.Date(init_time, format="%Y%m%d")

nc <- nc_open(input_filename)

tm <- ncvar_get(nc, "time")
tm_units <- ncatt_get(nc, "time", "units")$value
origin <- sub("(days since )([0-9]{4}-[0-9]{2}-[0-9]{2})", "\\2", tm_units)
origin <- as.Date(origin)
tm <- origin + days(tm)

stn <- ncvar_get(nc, "station")
idx <- match(efas_ids, stn)

dis <- ncvar_get(nc, "dis")
dis <- dis[idx, ,1:25]
dis <- aperm(dis, c(2, 3, 1)) # time, ensemble, station

data <- expand.grid(time = tm, member = 0:24, id = nrfa_ids) |> as_tibble()
data <- data |>
  mutate(member = as.character(member), id = as.character(id)) |>
  mutate(Qd = as.numeric(dis), init_time = as.Date(init_time)) |>
  pivot_longer(Qd, names_to = "variable", values_to = "value") |>
  dplyr::select(id, member, time, init_time, variable, value)

## Normalize by catchment area
data <- data |>
  left_join(catchment_area, by = "id") |>
  mutate(value = value / efas_catchment_area_ldd) |>
  mutate(lead_time = interval(init_time, time) %/% days(1))

nc_close(nc)


data <- data |>
  dplyr::select(
           id, member, init_time, lead_time, time,
           variable, value, efas_catchment_area_ldd
         )

write_parquet(data, output_filename)
