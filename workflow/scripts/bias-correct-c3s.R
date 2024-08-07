#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(lubridate)
library(arrow)
library(qmap)


if (exists("snakemake")) {
  c3s_input_file <- snakemake@input[[1]]
  obs_input_file <- snakemake@input[[2]]
  output_filename <- snakemake@output[[1]]
} else {
  c3s_input_file <- "results/preprocessing/C3S/10002.parquet"
  obs_input_file <- "results/preprocessing/HadUK/10002.parquet"
  output_filename <- "results/preprocessing/C3S/{station}_bc.parquet"
}

exp <- read_parquet(c3s_input_file)
obs <- read_parquet(obs_input_file)

## ############################## ##
## Algorithm 1: Use qmap directly ##
## ############################## ##

month_difference <- function(init_tm, tm) {
  diff <- tm - init_tm
  neg_index <- diff < 0
  diff[neg_index] <- diff[neg_index] + 12
  diff
}

## Add lead time month
exp <- exp |> mutate(lead_time_month = month_difference(month(init_time), month(time)))

# Add init time year 
exp <- exp |> mutate(init_year = lubridate::year(init_time))
obs <- obs |> mutate(year = lubridate::year(time))

# # Compute tmean 
# obs <- obs |> pivot_wider(names_from = variable, values_from = value)
# obs$tmean <- rowMeans(obs |> dplyr::select(tasmin, tasmax))
# obs <- obs |> pivot_longer(-(id:year), names_to = "variable", values_to = "value")

## Approach 1: All days of the year within the calibration data set are used
test_years <- seq(1994, 2016)

lead_time_months <- c(0, 1, 2) #, 3, 4, 5, 6)

members <- exp |> pull(member) |> unique() |> sort()

bc_list <- list()
for (i in 1:length(test_years)) {
  yr <- test_years[i]

  for (j in 1:length(lead_time_months)) { 
    mo <- lead_time_months[j]

    exp_train <- exp |> filter(init_year < yr & lead_time_month == mo)
    exp_test <- exp |> filter(init_year == yr & lead_time_month == mo)
    obs_train <- obs |> filter(year < yr)

    ## Temperature - no bias correction for now
    exp_train_temp <- exp_train |> 
      filter(variable %in% "t2m") |> 
      arrange(lead_time, init_time, member)
    exp_test_temp <- exp_test |> 
      filter(variable %in% "t2m")
    exp_test_temp <- exp_test_temp |> mutate(value_bc = value)
    bc_list[[length(bc_list) + 1]] <- exp_test_temp

    ## Precipitation
    exp_train_pr <- exp_train |> 
      filter(variable %in% "tp") |> 
      arrange(lead_time, init_time, member)
    exp_test_pr <- exp_test |> 
      filter(variable %in% "tp") |> 
      arrange(lead_time, init_time, member)
    obs_train_pr <- obs_train |> 
      filter(variable %in% "rainfall")

    ## Per member
    for (k in 1:length(members)) {
      num <- members[k]

      exp_train_pr0 <- exp_train_pr |> filter(member %in% num)
      exp_test_pr0 <- exp_test_pr |> filter(member %in% num)
      obs_train_pr0 <- exp_train_pr0 |> 
        dplyr::select(init_time, time) |> 
        left_join(obs_train_pr, by = "time") |> 
        dplyr::select(time, value)
      qm_fit <- fitQmapDIST(
        obs_train_pr0$value, 
        exp_train_pr0$value, 
        distr = "berngamma", 
        qstep = 0.001
      )
      qm <- doQmapDIST(exp_test_pr0$value, qm_fit)
      exp_test_pr0 <- exp_test_pr0 |> 
        mutate(value_bc = qm)
      bc_list[[length(bc_list) + 1]] <- exp_test_pr0 
    }
  }
}

exp_test_bc <- do.call("rbind", bc_list)

exp_test_bc <- exp_test_bc |> 
  arrange(lead_time, init_time, member, variable)

write_parquet(exp_test_bc, output_filename)