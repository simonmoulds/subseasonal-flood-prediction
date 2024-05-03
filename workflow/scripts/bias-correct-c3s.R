#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(lubridate)
library(arrow)
library(CSTools)


if (exists("snakemake")) {
  c3s_input_file <- snakemake@input[[1]]
  obs_input_file <- snakemake@input[[2]]
  output_filename <- snakemake@output[[1]]
} else {
  c3s_input_file <- "results/preprocessing/C3S/10002.parquet"
  obs_input_file <- "results/preprocessing/HadUK/10002.parquet"
  output_filename <- "results/preprocessing/ERA5/10002.parquet"
}


exp <- read_parquet(c3s_input_file)
obs <- read_parquet(obs_input_file)

## ############################## ##
## Algorithm 1: Use qmap directly ##
## ############################## ##

month_difference <- function(init_tm, tm) {
  diff <- tm - init_tm
  neg_index <- diff < 0
  diff[neg_index] = diff[neg_index] + 12
  diff
}

## Add lead time month
exp <- exp |> mutate(lead_time_month = month_difference(month(init_time), month(time)))

exp_train <- exp |> filter(init_time < as.Date("1994-01-01"))
exp_test <- exp |> filter(init_time >= as.Date("1994-01-01"))
obs_train <- obs |> filter(time >= as.Date("1981-01-01") & time < as.Date("1994-01-01"))
obs_test <- obs |> filter(time >= as.Date("1994-01-01"))

## Approach 1: All days of the year within the calibration data set are used
lead_time_months <- c(0, 1, 2, 3, 4, 5, 6)

## Bias-correct each lead time month in turn
exp_train <- exp_train |> filter(lead_time_month %in% lead_time_months[1])

## Precipitation
exp_train_pr <- exp_train |> filter(variable %in% "tp") |> arrange(lead_time, init_time, member)
obs_train_pr <- obs_train |> filter(variable %in% "rainfall")

## Per member
members <- exp_train_pr |> pull(member) |> unique()
for (i in 1:length(members)) {
  exp_train_pr0 <- exp_train_pr |> filter(member %in% members[i])
  obs_train_pr0 <- exp_train_pr0 |> dplyr::select(init_time, time) |> left_join(obs_train_pr, by = "time") |> dplyr::select(time, value)
  qm_fit <- fitQmapDIST(obs_train_pr0$value, exp_train_pr0$value, distr = "berngamma", qstep = 0.001)
  qm <- doQmapDIST(exp_train_pr0$value, qm_fit)
}
## Convert to arrays with the format expected by CSTools

## ######################################## ##
## Algorithm 2: Use CSTools wrapper to qmap ##
## ######################################## ##
exp_train <- exp |> filter(init_time < as.Date("1994-01-01"))
exp_test <- exp |> filter(init_time >= as.Date("1994-01-01"))
obs_train <- obs |> filter(time >= as.Date("1981-01-01") & time < as.Date("1994-01-01"))
obs_test <- obs |> filter(time >= as.Date("1994-01-01"))

exp_train <- exp_train |> filter(lead_time %in% 1:30)

## Precipitation
exp_train_pr <- exp_train |> filter(variable %in% "tp") |> arrange(lead_time, init_time, member)
obs_train_pr <- obs_train |> filter(variable %in% "rainfall")

members <- exp_train_pr |> pull(member) |> unique()
sdates <- exp_train_pr |> pull(init_time) |> unique()
fdates <- exp_train_pr |> pull(lead_time) |> unique()

## Complete specification
grd <- expand_grid(lead_time = fdates, init_time = sdates, member = members) |> dplyr::select(member, init_time, lead_time)
exp_train_pr <- grd |> left_join(exp_train_pr, by = c("member", "init_time", "lead_time"))

exp_arr <- exp_train_pr$value
dim(exp_arr) <- c(member = length(members), sdate = length(sdates), ftime = length(fdates))
## ## Check array assignment - OK
## member_ix = 25
## sdate_ix = 136
## all.equal(exp_arr[member_ix,sdate_ix,], exp_train_pr |> filter(member %in% members[member_ix] & init_time %in% sdates[sdate_ix]) |> dplyr::pull(value))

grd <- expand_grid(lead_time = fdates, init_time = sdates, member = 0) |> mutate(time = (init_time + days(lead_time)))
obs_train_pr <- grd |> left_join(obs_train_pr, by = "time")

obs_arr <- obs_train_pr$value
dim(obs_arr) <- c(member = 1, sdate = length(sdates), ftime = length(fdates))
## ## Check array assignment
## sdate_ix = 1
## all.equal(obs_arr[1,sdate_ix,], obs_train_pr |> filter(init_time %in% sdates[sdate_ix]) |> dplyr::pull(value))

res <- QuantileMapping(exp_arr, obs_arr)


## ## Open files as a dataset and filter on ID
## ds <- open_dataset(input_files, unify_schemas = FALSE)
## ds <- ds %>% filter(id %in% station)
## ds <- ds %>% collect()

## write_parquet(ds, output_filename)

## rootdir <- "/exports/geos.ed.ac.uk/moulds_hydro/c3s-data-download/C3S_hindcast_daily_SEAS5"
## fpath <- file.path(rootdir, "ecmf_SEAS5-v20171101_hindcast_S$sdate$_$var$_daily.nc")
## sdates <- c("20100101", "20100201", "20100301")

## x <- CST_Start(
##   dataset = fpath,
##   sdate = sdates,
##   var = "tp",
##   latitude = 'all',
##   longitude = 'all',
##   number = 'all',
##   time = 'all',
##   ## return_vars = list(latitude = NULL, longitude = NULL, number = NULL, time = 'sdate'),
##   retrieve = TRUE
## )

## ## Use synthetic data
## exp <- 1 : c(1 * 3 * 5 * 4 * 3 * 2)
## dim(exp) <- c(dataset = 1, member = 3, sdate = 5, ftime = 4, lat = 3, lon = 2)
## obs <- 101 : c(100 + 1 * 1 * 5 * 4 * 3 * 2)
## dim(obs) <- c(dataset = 1, member = 1, sdate = 5, ftime = 4, lat = 3, lon = 2)
## res <- QuantileMapping(exp, obs)

## ## OK, so do this using ncdf4 to read data?
## nc <- nc_open(file.path(rootdir, "ecmf_SEAS5-v20171101_hindcast_S20100101_tp_daily.nc"))
## tp <- ncvar_get(nc, "tp", count = c(10, 10, -1, -1))
## tp <- aperm(tp, c(4, 3, 2, 1))

## ## How to do this practically?
