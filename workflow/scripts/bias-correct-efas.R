#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(lubridate)
library(arrow)
library(qmap)

if (exists("snakemake")) {
  efas_input_file <- snakemake@input[[1]]
  obs_input_file <- snakemake@input[[2]]
  output_filename <- snakemake@output[[1]]
} else {
  efas_input_file <- "results/preprocessing/EFAS/10002.parquet"
  obs_input_file <- "results/preprocessing/streamflow/timeseries/daily/10002.parquet"
  output_filename <- "results/preprocessing/EFAS/10002_bc.parquet"
}

exp <- read_parquet(efas_input_file)
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
exp <- exp |> mutate(lead_time = as.integer(time - init_time), lead_time_month = month_difference(month(init_time), month(time)))

# Add init time year 
exp <- exp |> mutate(init_year = lubridate::year(init_time))
obs <- obs |> mutate(year = lubridate::year(time)) |> pivot_longer(Qd, names_to = "variable", values_to = "value")

# # Compute tmean 
# obs <- obs |> pivot_wider(names_from = variable, values_from = value)
# obs$tmean <- rowMeans(obs |> dplyr::select(tasmin, tasmax))
# obs <- obs |> pivot_longer(-(id:year), names_to = "variable", values_to = "value")

## Approach 1: All days of the year within the calibration data set are used
test_years <- seq(2004, 2016)

lead_time_months <- c(0, 1, 2, 3) #, 4, 5, 6)

members <- exp |> pull(member) |> unique() |> sort()

bc_list <- list()
for (i in 1:length(test_years)) {
  yr <- test_years[i]

  for (j in 1:length(lead_time_months)) { 
    mo <- lead_time_months[j]

    exp_train <- exp |> filter(init_year < yr & lead_time_month == mo)
    exp_test <- exp |> filter(init_year == yr & lead_time_month == mo)
    obs_train <- obs |> filter(year < yr)

    # Streamflow 
    exp_train_q <- exp_train |> 
        filter(variable %in% "Qd") |> 
        arrange(lead_time, init_time, member) 

    exp_test_q <- exp_test |> filter(variable %in% "Qd")

    obs_train_q <- obs_train |> filter(variable %in% "Qd") 

    ## Per member
    for (k in 1:length(members)) {
      num <- members[k]

      exp_train_q0 <- exp_train_q |> filter(member %in% num)
      exp_test_q0 <- exp_test_q |> filter(member %in% num)
      obs_train_q0 <- exp_train_q0 |> 
        dplyr::select(init_time, time) |> 
        left_join(obs_train_q, by = "time") |> 
        dplyr::select(time, value)
      
      qm_fit <- fitQmapQUANT(
        obs_train_q0$value, 
        exp_train_q0$value,
        qstep = 0.001
      )

      # Handle missing values
      test_values <- exp_test_q0$value 
      qm <- rep(NA, length(test_values))
      na_index <- is.na(test_values)
      test_values <- test_values[!na_index]
      qm_nonmissing <- doQmapQUANT(test_values, qm_fit)
      qm[!na_index] <- qm_nonmissing 

      exp_test_q0 <- exp_test_q0 |> 
        mutate(value_bc = qm)
      bc_list[[length(bc_list) + 1]] <- exp_test_q0 
    }
  }
}

exp_test_bc <- do.call("rbind", bc_list)

exp_test_bc <- exp_test_bc |> 
  arrange(lead_time, init_time, member, variable)

write_parquet(exp_test_bc, output_filename)