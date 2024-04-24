#!/usr/bin/env Rscript

library(tidyverse)
library(lubridate)
library(arrow)

options(dplyr.summarise.inform = FALSE)

if (exists("snakemake")) {
  station <- snakemake@wildcards[["station"]]
  normalize <- snakemake@params[[1]]
  daily_output_filename <- snakemake@output[[1]]
  snakemake@source("utils.R")
} else {
  ## TESTING
  station <- 10002
  normalize <- TRUE
  daily_output_filename <- "results/preprocessing/streamflow/timeseries/daily/24001.parquet"
  source("workflow/scripts/utils.R")
}

## 1 - Catchment metadata
metadata <- read_parquet("resources/stations_selected.parquet")
metadata <- metadata %>% filter(id %in% station)
catchment_area <- metadata$area[1]

## 2 - Read raw data
x <- read_csv(file.path("data/NRFA", paste0(station, ".csv")))
x <- x %>%
  rename(Qd = gdf) %>%
  mutate(time = as.Date(time)) %>%
  mutate(id = as.character(station), .before = time)

if (normalize) {
  ## We convert to mm/day
  x <- x %>% mutate(Qd = (Qd * 60 * 60 * 24) / (catchment_area * 1000 * 1000) * 1000)
}

## ## 3 - Get monthly data
## x_month <- x %>%
##   pivot_longer(-(id:time), names_to = "variable", values_to = "value") %>%
##   mutate(year = lubridate::year(time) %>% unname(), month = lubridate::month(time) %>% unname()) %>%
##   summary_stats_fun("value", c("id", "year", "month", "variable"), probs = NULL) %>%
##   unite(variable, variable, summary_statistic) %>%
##   mutate(time = lubridate::make_date(year, month, 1) %>% unname(), .before = year)

## x_month <- x_month %>%
##   mutate(
##     id = as.character(id),
##     time = as.Date(time),
##     year = as.integer(year),
##     month = as.integer(month),
##     variable = as.character(variable),
##     value = as.numeric(value)
##   )

write_parquet(x, daily_output_filename)
## write_parquet(x_month, monthly_output_filename)
