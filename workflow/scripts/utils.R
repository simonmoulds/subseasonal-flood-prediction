## Author : Simon Moulds
## Date   : Nov-Dec 2021

library(tidyverse)
library(arrow)
library(zoo)
library(terra)
library(exactextractr)
library(lubridate)
library(yaml)


options(readr.show_col_types = FALSE)
options(dplyr.summarise.inform = FALSE)
options(readr.show_progress = FALSE)


## check_data_availability <- function(datadir,
##                                     stations,
##                                     max_start_date,
##                                     min_end_date,
##                                     min_data_availability = 0.95, ...) {
##   ## Identify stations which do not meet a minimum data
##   ## availability per month
##   exclude <- c()
##   n_stations <- length(stations)
##   pb <- txtProgressBar(min = 0, max = n_stations, initial = 0)
##   for (i in 1:n_stations) {
##     stn <- stations[i]
##     fn <- file.path(datadir, paste0(stn, ".csv"))
##     if (!file.exists(fn)) {
##       exclude <- c(exclude, stn)
##       next
##     }
##     x <- read_csv(fn, show_col_types = FALSE, progress = FALSE)
##     if (nrow(x) == 0) {
##       exclude <- c(exclude, stn)
##       next
##     }
##     ts <- tibble(
##       time = seq.Date(max_start_date, min_end_date, by = "1 day")
##     )
##     x <- x %>%
##       mutate(
##         time = as.Date(time),
##         year = lubridate::year(time),
##         month = lubridate::month(time)
##       )
##     x <- left_join(ts, x, by = "time")
##     missing_pct <- x %>%
##       group_by(year) %>%
##       ## group_by(year, month) %>%
##       summarize(missing_pct = sum(is.na(gdf)) / n())
##       ## summarize(missing_pct = sum(is.na(value)) / n())

##     if (any(missing_pct$missing_pct > (1 - min_data_availability))) {
##       exclude <- c(exclude, stn)
##     }

##     ## Update progress bar
##     setTxtProgressBar(pb, i)
##   }
##   close(pb)
##   return(exclude)
## }


rename_variables <- function(ds, ...) {
  rainfall_names <- c("single_level_tp", "rainfall", "tprate")
  temperature_names <- c("single_level_t2m", "tas", "t2m")
  ds <- ds %>%
    mutate(
      variable = case_when(
        variable %in% rainfall_names ~ "precip",
        variable %in% temperature_names ~ "temp",
        .default = variable
      )
    )
  ds
}


load_dataset <- function(root, partitioning, ...) {
  ds <-
    open_dataset(root, partitioning = partitioning) %>%
    collect() %>%
    rename_variables()
  ds
}


load_era5_dataset <- function(root, ...) {
  partitioning <- c("variable")
  ds <- load_dataset(root, partitioning)
  ds <- ds %>%
    mutate(
      value = case_when(
        variable %in% "precip" ~ value * 1000,
        .default = value
      ),
      product_type = "monthly_mean"
    ) %>%
    dplyr::select(-days_in_month)
  ds <- ds %>%
    unite(variable, variable, product_type, sep = "_")
  ds
}


load_haduk_dataset <- function(root, ...) {
  partitioning <- c("variable")
  ds <- load_dataset(root, partitioning)
  ds <- ds %>%
    mutate(
      value = case_when(
        variable %in% "precip" ~ value / days_in_month,
        variable %in% "temp" ~ value + 273.15,
        .default = value
      ),
      product_type = "monthly_mean"
    ) %>%
    dplyr::select(-days_in_month)
  ds <- ds %>%
    unite(variable, variable, product_type, sep = "_")
  ds
}


load_c3s_dataset <- function(root, ...) {
  partitioning <- c("variable")
  ds <- load_dataset(root, partitioning)
  ds <- ds %>%
    mutate(
      value = case_when(
        variable %in% "precip" ~ value * 86400 * 1000,
        .default = value
      )
    ) %>%
    dplyr::select(-days_in_month)
  ds <- ds %>%
    unite(variable, variable, product_type, sep = "_")
  ds
}


get_month_index <- function(season) {
  ## get_month_index("JJA")
  ## get_month_index("12")
  ## get_month_index("1")
  ## get_month_index("DJFM")
  ## get_month_index("SON")
  index <- suppressWarnings(as.numeric(season))
  if (!is.na(index)) {
    return(index)
  }
  season <- toupper(season)
  months <- "JFMAMJJASONDJFMAMJJASOND"
  match <- regexpr(toupper(season), months)
  if (match == -1) {
    stop("Invalid season!")
  }
  index <- seq(match, length.out = nchar(season))
  index[index > 12] <- index[index > 12] %% 12
  return(index)
}


## From the C3S data we want to take various statistics
summary_stats_fun <- function(ds, value_column, group_vars, probs, ...) {
  ## Quantiles
  funs <- list(mean = mean, minimum = min, maximum = max)
  if (length(probs) > 0) {
    p <- probs
    p_names <- map_chr(p, ~paste0("q", formatC(.x * 100, width = 2, flag = 0)))
    p_funs <-
      map(p, ~partial(quantile, probs = .x, na.rm = TRUE)) %>%
      set_names(nm = p_names)
    funs <- c(p_funs, funs)
  }
  ## Add summary statistics
  ds_sum <- ds %>%
    group_by(across(all_of(group_vars))) %>%
    summarize(across(all_of(value_column), funs, .names = "{.fn}")) %>%
    ungroup() %>%
    pivot_longer(-all_of(group_vars), names_to = "summary_statistic", values_to = "value") %>%
    mutate(across(all_of("value"), unname))
  ds_sum
}


get_season_ts <- function(x) {
  ts <- x %>%
    distinct(time) %>%
    mutate(year = lubridate::year(time), month = lubridate::month(time)) %>%
    filter(month %in% season_index) %>%
    mutate(start_time = case_when(month == 12 ~ time)) %>%
    mutate(start_time = zoo::na.locf(start_time, na.rm = FALSE)) %>%
    dplyr::select(-year, -month)
  ts
}


get_monthly_ts <- function(x, season_index) {
  ts <- x %>%
    distinct(id, time) %>%
    mutate(
      year = lubridate::year(time) %>% unname(),
      month = lubridate::month(time) %>% unname()
    ) %>%
    filter(month %in% season_index) %>%
    mutate(start_time = case_when(month == season_index[1] ~ time)) %>%
    mutate(start_time = zoo::na.locf(start_time, na.rm = FALSE)) %>%
    dplyr::select(-year, -month)

  ## Restrict to seasons with the correct number of months
  ts <- ts %>%
    group_by(start_time) %>%
    filter(n() == length(season_index))
  ts
}


## get_daily_ts <- function(x, season_index) {
##   ## NOTE streamflow data is daily data
##   ts <- x %>%
##     distinct(time) %>%
##     arrange(time) %>%
##     mutate(year = lubridate::year(time), month = lubridate::month(time)) %>%
##     filter(month %in% season_index) %>%
##     group_by(year, month) %>%
##     mutate(month_start = min(time)) %>%
##     ungroup() %>%
##     mutate(start_time = case_when(month == season_index[1] ~ month_start)) %>%
##     fill(start_time, .direction = "down") %>%
##     dplyr::select(-year, -month, -month_start)
##   ts
## }


## aggregate_q_data <- function(x, season_index) {
##   ts <- get_daily_ts(x, season_index)
##   ts <- expand_grid(ts, distinct(x, variable))
##   ts <- ts %>% filter(!is.na(start_time))
##   ts <- ts %>% left_join(x, by = c("time", "variable"))
##   probs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
##   ts <- summary_stats_fun(
##     ts, "value", c("start_time", "variable"),
##     probs = probs
##   )
##   ts
## }


## aggregate_era5_data <- function(x, season_index) {
##   ts <- get_monthly_ts(x, season_index)
##   ts <- expand_grid(ts, distinct(x, variable))
##   ts <- ts %>% left_join(x, by = c("time", "variable"))
##   ts <- summary_stats_fun(
##     ts, "value", c("start_time", "variable"), probs = NULL
##   )
##   ts
## }


aggregate_c3s_data <- function(x, season_index) {

  ts <- get_monthly_ts(x, season_index)

  ## Expand ts to include all possible lead times
  lead_times <- x$lead_time %>% unique() %>% sort()
  ts <- expand_grid(ts, tibble(latest_lead_time = lead_times)) %>%
    arrange(start_time, latest_lead_time)

  ## Given the lead time of the latest month, get the lead time of all
  ## prior months and use this to compute the initialization time of the
  ## forecast for each month
  ts <- ts %>%
    group_by(id, start_time, latest_lead_time) %>%
    mutate(lead_time = latest_lead_time - ((n() - 1):0)) %>%
    ungroup()

  ## Expand ts to include all source, member, variable combinations
  ## NOTE check this works when multiple source_id present
  ## NOTE this is needed to avoid a many-to-many join, which occurs
  ## because there are multiple members with the same values of `time`,
  ## `start_time` and `init_time`
  ## ts <- expand_grid(ts, distinct(x, source_id, member, variable))
  ts <- expand_grid(ts, distinct(x, system, member, variable))

  ## ## OPTION 1 - include forecasts from before the initialization date
  ## ts <- ts %>%
  ##   mutate(lead_time = pmax(lead_time, 0)) %>%
  ##   mutate(init_time = time - months(lead_time))
  ## ts <- ts %>% left_join(c3s_ds)
  ## ts <- ts %>%
  ##   group_by(start_time, latest_lead_time, source_id, member, variable) %>%
  ##   summarize(n=n()) # FIXME (n is calculated as a check)

  ## OPTION 2 - discard forecasts from before the initialization date
  ## NOTE at this stage we could also discard seasons that do not
  ## contain a forecast for each month
  ## ts <- ts %>% left_join(x, by = c("time", "lead_time", "source_id", "member", "variable"))
  ts <- ts %>% left_join(x, by = c("id", "time", "lead_time", "system", "member", "variable"))

  ## Only keep forecasts with a positive lead time
  ts <- ts %>% filter(lead_time >= 0)
  ts <- ts %>%
    ## group_by(start_time, latest_lead_time, source_id, member, variable) %>%
    group_by(id, start_time, latest_lead_time, system, member, variable) %>%
    summarize(
      n=n(),
      value = mean(value),
      start_time = min(start_time),
      init_time = min(init_time)
    ) %>%
    mutate(
      year = lubridate::year(start_time),
      ## month = lubridate::month(start_time),
      .after = start_time
    )

  ts <- ts %>% dplyr::select(id, system, member, start_time, year, n, init_time, latest_lead_time, variable, value)
  ts
}


## compute_era5_antecedent_climate <- function(ds, ...) {
##   ds <- ds %>%
##     mutate(
##       precip_monthly_mean_3 = zoo::rollmean(precip_monthly_mean, 3, na.pad = TRUE, align = "right"),
##       temp_monthly_mean_3 = zoo::rollmean(temp_monthly_mean, 3, na.pad = TRUE, align = "right"),
##       precip_monthly_mean_2 = zoo::rollmean(precip_monthly_mean, 2, na.pad = TRUE, align = "right"),
##       temp_monthly_mean_2 = zoo::rollmean(temp_monthly_mean, 2, na.pad = TRUE, align = "right")
##     )
##   antecedent_climate <- ds %>%
##     mutate(
##       ant_precip_monthly_mean_1 = lag(precip_monthly_mean, n = 1),
##       ant_precip_monthly_mean_2 = lag(precip_monthly_mean_2, n = 1),
##       ant_precip_monthly_mean_3 = lag(precip_monthly_mean_3, n = 1),
##       ant_temp_monthly_mean_1 = lag(temp_monthly_mean, n = 1),
##       ant_temp_monthly_mean_2 = lag(temp_monthly_mean_2, n = 1),
##       ant_temp_monthly_mean_3 = lag(temp_monthly_mean_3, n = 1)
##     ) %>%
##     dplyr::select(time, year, month, starts_with("ant_"))
##   antecedent_climate
## }


## compute_antecedent_streamflow <- function(ds, ...) {
##   ds <- ds %>%
##     group_by(variable) %>%
##     mutate(
##       value_2 = zoo::rollmean(value, 2, na.pad = TRUE, align = "right"),
##       value_3 = zoo::rollmean(value, 3, na.pad = TRUE, align = "right"),
##       ant_1 = lag(value, n = 1),
##       ant_2 = lag(value_2, n = 1),
##       ant_3 = lag(value_3, n = 1)
##     )
##   ## ds <- ds %>%
##   ##   mutate(
##   ##     Qd_mean_2 = zoo::rollmean(Qd_mean, 2, na.pad = TRUE, align = "right"),
##   ##     Qd_mean_3 = zoo::rollmean(Qd_mean, 3, na.pad = TRUE, align = "right")
##   ##   )
##   ## antecedent_streamflow <- ds %>%
##   ##   mutate(
##   ##     ant_Qd_mean_1 = lag(Qd_mean, n = 1),
##   ##     ant_Qd_mean_2 = lag(Qd_mean_2, n = 1),
##   ##     ant_Qd_mean_3 = lag(Qd_mean_3, n = 1)
##   ##   ) %>%
##   ##   dplyr::select(time, year, month, starts_with("ant_"))
##   ## antecedent_streamflow
##   ds
## }

compute_antecedent_values <- function(ds, ...) {
  ds <- ds %>%
    group_by(variable) %>%
    mutate(
      value_2 = zoo::rollmean(value, 2, na.pad = TRUE, align = "right"),
      value_3 = zoo::rollmean(value, 3, na.pad = TRUE, align = "right"),
      ant_1 = lag(value, n = 1),
      ant_2 = lag(value_2, n = 1),
      ant_3 = lag(value_3, n = 1)
    )
  ds <- ds %>% dplyr::select(-contains("value")) %>% pivot_longer(starts_with("ant_")) %>% unite(variable, variable, name) %>% pivot_wider(names_from=variable, values_from=value)
  ds
}

compute_antecedent_forecast_values <- function(ds, ...) {
  print("Hello, world")
}

## Compute indices
baseflow_sep <- function(strflow) {
  # Baseflow filter
  # Recursive digital filter techniques
  # ...

  # Initialization
  f1 <- 0.95
  f2 <- (1 + f1) / 2
  surfq <- strflow
  baseq <- matrix(rep(strflow, 3), ncol = 3, byrow = TRUE)
  surfq[1] <- strflow[1] * 0.5
  baseq[1, 1] <- strflow[1] - surfq[1]
  baseq[1, 2] <- baseq[1, 1]
  baseq[1, 3] <- baseq[1, 1]

  # First pass (forward)
  for (i in 2:length(strflow)) {
    surfq[i] <- f1 * surfq[i-1] + f2 * (strflow[i] - strflow[i-1])
    surfq[i] <- max(0, surfq[i])

    baseq[i, 1] <- strflow[i] - surfq[i]
    baseq[i, 1] <- max(0, baseq[i, 1])
    baseq[i, 1] <- min(strflow[i], baseq[i, 1])
  }

  # Second pass (backward)
  baseq[length(strflow)-1, 2] <- baseq[length(strflow)-1, 1]
  for (i in (length(strflow)-2):1) {
    surfq[i] <- f1 * surfq[i+1] + f2 * (baseq[i, 1] - baseq[i+1, 1])
    surfq[i] <- max(0, surfq[i])

    baseq[i, 2] <- baseq[i, 1] - surfq[i]
    baseq[i, 2] <- max(0, baseq[i, 2])
    baseq[i, 2] <- min(baseq[i, 1], baseq[i, 2])
  }

  # Third pass (forward)
  baseq[length(strflow)-1, 3] <- baseq[length(strflow)-1, 1]
  for (i in 2:length(strflow)) {
    surfq[i] <- f1 * surfq[i-1] + f2 * (baseq[i, 2] - baseq[i-1, 2])
    surfq[i] <- max(0, surfq[i])

    baseq[i, 3] <- baseq[i, 2] - surfq[i]
    baseq[i, 3] <- max(0, baseq[i, 3])
    baseq[i, 3] <- min(baseq[i, 2], baseq[i, 3])
  }

  return(baseq)
}

compute_bfi <- function(q) {
  qb <- baseflow_sep(q)[,3]
  bfi <- sum(qb, na.rm = TRUE) / sum(q, na.rm = TRUE)
  return(bfi)
}

compute_slope_fdc <- function(q, bins) {
  q[q == 0] <- 1e-6
  percentiles <- quantile(log(q), probs = bins, na.rm = TRUE)
  diff_percentiles <- diff(percentiles) %>% unname
  ## diff_percentiles <- log(percentiles[2]) - log(percentiles[1])
  diff_bins <- diff(bins)
  slp <- (diff_percentiles / (diff_bins * 100)) * 100
  return(slp)
}

compute_indices <- function(data) {
  Q5 <- data$Qd %>% quantile(probs = 0.05, na.rm = TRUE)
  Q95 <- data$Qd %>% quantile(probs = 0.95, na.rm = TRUE)
  Q50 <- data$Qd %>% quantile(probs = 0.5, na.rm = TRUE)
  q_mean <- data$Qd %>% mean(na.rm = TRUE)
  baseflow_index <- compute_bfi(data$Qd)
  ## slope_fdc <- compute_slope_fdc(data$Qd, c(0.33, 0.67))
  slope_fdc <- compute_slope_fdc(data$Qd, c(0.33, 0.67))
  ## return(tibble(q05 = q05, q95 = q95, q50 = q50, qmean = qmean, bfi = bfi, slope_fdc = slope_fdc))
  return(tibble(Q5_comp = Q5, Q95_comp = Q95, q_mean_comp = q_mean, baseflow_index_comp = baseflow_index, slope_fdc_comp = slope_fdc))
}

compute_lagged_indices <- function(x, min_years = 10) {
  x <- x %>% mutate(year = lubridate::year(time))
  years <- x$year %>% unique %>% sort
  n_years <- length(years)
  if (n_years <= min_years) {
    stop("Insufficient data!")
  }
  compute_years <- years[(min_years + 1):n_years]

  indices_list <- list()
  for (i in 1:length(compute_years)) {
    data <- x %>% filter(year < compute_years[i])
    indices <- compute_indices(data)
    indices <- indices %>% mutate(year = compute_years[i])
    indices_list[[length(indices_list) + 1]] <- indices
  }
  indices <- do.call("rbind", indices_list)
  return(indices)
}
## NOT USED:
##
## get_month_index <- function(season_str) {
##   months <- "JFMAMJJASOND"
##   months <- paste0(months, months)
##   if (!str_detect(months, season)) {
##     stop("Invalid season")
##   }
##   months_index <- str_locate(months, season)
##   months_index <- seq(months_index[1], months_index[2])
##   months_index[months_index > 12] <- months_index[months_index > 12] - 12
##   months_index
## }

## get_season_input_data <- function(lead_time,
##                                   dep_var,
##                                   subset,
##                                   formulas,
##                                   years,
##                                   station = NULL, ...) {

##   ## Load input dataset
##   ## ##################
##   input_dir <- file.path("results/input/combined")
##   ds <- open_dataset(
##     input_dir,
##     partitioning = c("id", "season", "lead_time")
##   ) %>%
##     collect()

##   ## Filter lead time
##   ds <- ds %>% filter(lead_time %in% !!lead_time)

##   ## Filter stations, if provided
##   if (!is.null(station)) {
##     ds <- ds %>% filter(id %in% !!station)
##   }

##   ## Rename dependent variable `Q`
##   ds <- ds %>% rename(Q = !!sym(dep_var))

##   ## One-hot encoding
##   ## ################
##   ds <- ds %>%
##     mutate(
##       season_DJF = ifelse(season == "DJF", 1, 0),
##       season_MAM = ifelse(season == "MAM", 1, 0),
##       season_JJA = ifelse(season == "JJA", 1, 0),
##       season_SON = ifelse(season == "SON", 1, 0)
##     )

##   ## Check for presence of NA values in dataset
##   ## ##########################################
##   get_vars <- function(formula) {
##     ## Function to retrieve formula variables
##     formula <- paste0(dep_var, "~", paste0(formula, collapse="+"))
##     formula <- as.formula(formula)
##     vars <- all.vars(formula)
##     vars
##   }
##   vars <- lapply(formulas, FUN = function(x) get_vars(x))
##   vars <- do.call("c", vars) %>% unique()
##   ## Some variables are computed dynamically - we only
##   ## check those that have to be present from the outset
##   vars <- vars[vars %in% names(ds)]
##   ## Remove rows with missing values for these variables
##   na_ix <- complete.cases(ds[, vars])
##   ds <- ds[na_ix, ]

##   ## Select training and testing dataset
##   ## ###################################
##   all_train_data <- ds %>%
##     filter(subset %in% !!subset) %>%
##     filter(season_year %in% years) %>%
##     filter(!is.na(Q)) %>%
##     arrange(season_year)

##   all_test_data <- ds %>%
##     filter(subset %in% !!subset) %>%
##     filter(season %in% !!season) %>%
##     filter(season_year %in% years)

##   ## Return list
##   list(train = all_train_data, test = all_test_data)
## }


## summarise_discharge_data <- function(x, season) {

##   ## ## Previous POT analysis
##   ## ## Neri et al [https://doi.org/10.1002/joc.5915]:
##   ## ## "To avoid double counting the same event, we only consider
##   ## ## one event in a window of +/- 5 days + logarithm of the
##   ## ## drainage area"
##   ## catchment_area = metadata[["catchment_area"]]
##   ## window_size = ((5 + log(catchment_area * 0.386102)) * 2) %>% round()
##   ## ## Solari et al [https://doi.org/10.1002/2016WR019426]:
##   ## ## "As the moving window travels through the series, each
##   ## ## time that the data maximum in the window is located at
##   ## ## its center, the maximum is regarded as a peak"
##   ## df =
##   ##   df %>% # Neri et al
##   ##   mutate(pot = roll_max(gdf, n=window_size, align="center", fill=NA)) %>%
##   ##   mutate(is_peak = Vectorize(isTRUE)(pot == gdf)) %>%
##   ##   mutate(peak = ifelse(is_peak, pot, NA))
##   ## ## Unsure whether we need to make this unique or not?
##   ## peaks = df$pot[df$is_peak] ##%>% unique()
##   ## n_years = length(df$year %>% unique())
##   ## peaks_sorted = sort(peaks, decreasing = TRUE)
##   ## threshold_1 = peaks_sorted[n_years]
##   ## threshold_2 = peaks_sorted[n_years * 2]
##   ## threshold_3 = peaks_sorted[n_years * 3]
##   ## threshold_4 = peaks_sorted[n_years * 4]
##   ## df =
##   ##   df %>%
##   ##   mutate(pot_1 = ifelse(Vectorize(isTRUE)(peak >= threshold_1), 1, 0)) %>%
##   ##   mutate(pot_2 = ifelse(Vectorize(isTRUE)(peak >= threshold_2), 1, 0)) %>%
##   ##   mutate(pot_3 = ifelse(Vectorize(isTRUE)(peak >= threshold_3), 1, 0)) %>%
##   ##   mutate(pot_4 = ifelse(Vectorize(isTRUE)(peak >= threshold_4), 1, 0))

##   ## Add climate season label to rows
##   months <- get_month_index(season)
##   x <- x %>% filter(month %in% months) %>% mutate(clim_season = season)
##   if (which.min(months) != 1) {
##     ## This adjusts season_year if the season covers two years (e.g.DJF)
##     next_year_months <- seq(min(months), months[length(months)])
##     x <- x %>% mutate(season_year = ifelse(month %in% next_year_months, year-1, year))
##   } else {
##     x <- x %>% mutate(season_year = year)
##   }

##   ## Quantiles per season
##   thresholds <- x %>%
##     group_by(ID, clim_season) %>%
##     summarize(
##       Q_99_threshold = quantile(Q, 0.99, na.rm = TRUE),
##       Q_95_threshold = quantile(Q, 0.95, na.rm = TRUE)
##     )

##   x <- x %>% left_join(thresholds, by=c("ID", "clim_season"))

##   ## Summarize to get flood counts
##   x <- x %>%
##     group_by(ID, clim_season, season_year) %>%
##     summarize(
##       missing_pct = (sum(is.na(Q)) / n()) * 100,
##       Q_max = max(Q, na.rm = TRUE),
##       Q_mean = mean(Q, na.rm = TRUE),
##       Q_01 = quantile(Q, probs = 0.01, na.rm = TRUE, names = FALSE),
##       Q_05 = quantile(Q, probs = 0.05, na.rm = TRUE, names = FALSE),
##       Q_50 = quantile(Q, probs = 0.50, na.rm = TRUE, names = FALSE),
##       Q_90 = quantile(Q, probs = 0.90, na.rm = TRUE, names = FALSE),
##       Q_95 = quantile(Q, probs = 0.95, na.rm = TRUE, names = FALSE),
##       Q_99 = quantile(Q, probs = 0.99, na.rm = TRUE, names = FALSE),
##       POT_1 = sum(Q > Q_99_threshold, na.rm = TRUE),
##       POT_2 = sum(Q > Q_95_threshold, na.rm = TRUE)
##     ) %>%
##     ungroup() %>%
##     mutate(across(Q_max:POT_2, ~ifelse(is.finite(.), ., NA)))
##   x
## }

## ## download_nrfa_data <- function(stn_id) {
## ##   ## Gauged daily flow [m3 s-1]
## ##   gdf <- get_ts(stn_id, "gdf") %>% as_tibble(rownames="time")
## ##   ## Create complete time series, in case the
## ##   ## raw time series has missing values.
## ##   start_date <- gdf$time[1]
## ##   end_date <- gdf$time[nrow(gdf)]
## ##   complete_ts <- seq.POSIXt(
## ##     as.POSIXct(start_date, tz="GMT", format="%Y-%m-%d"),
## ##     as.POSIXct(end_date, tz="GMT", format="%Y-%m-%d"),
## ##     by="1 day"
## ##   ) %>% as.Date() %>% as.character()
## ##   gdf <-
## ##     tibble(time = complete_ts) %>%
## ##     left_join(gdf, by="time")
## ##   availability = sum(!is.na(gdf$gdf)) / nrow(gdf) * 100
## ##   x <- gdf %>%
## ##     rename(Q = gdf) %>%
## ##     mutate(ID=stn_id, .after=time) %>%
## ##     mutate(time = as.Date(time))
## ##   x <- x %>%
## ##     mutate(year = format(time, "%Y") %>% as.integer) %>%
## ##     mutate(month = format(time, "%m") %>% as.integer)
## ##   x
## ## }

## make_qrf_forecast <- function(formula,
##                               train_data,
##                               test_data,
##                               quantiles = c(0.05, 0.25, 0.5, 0.75, 0.95),
##                               ...) {

##   ## https://stackoverflow.com/a/75455038
##   ## Use `ranger` to give probabilistic random forest prediction
##   model <- ranger(
##     formula = as.formula(formula),
##     data = train_data,
##     num.trees = 500,
##     min.node.size = 5,
##     ## min.bucket = 1,
##     max.depth = NULL,
##     quantreg = TRUE,
##     num.threads = parallel::detectCores() - 1
##   )

##   ## Make prediction
##   pred <- predict(
##     model,
##     data = test_data,
##     type = "quantiles",
##     quantiles = quantiles
##   )

##   nms <- (quantiles * 100) %>%
##     formatC(width = 2, flag = 0) %>%
##     paste0("Q", .)

##   pred <- pred$predictions %>%
##     as_tibble() %>%
##     set_names(nms) %>%
##     mutate(id = test_data$id) %>%
##     mutate(season_year = test_data$season_year)

##   df <- test_data %>%
##     dplyr::select(id, season_year, Q) %>%
##     rename(Q_obs = Q)

##   df <- df %>% left_join(pred, by = c("id", "season_year"))

##   ## Generate ensemble forecast
##   ens_fcst <- predict(
##     model,
##     data=test_data,
##     type = "quantiles",
##     quantiles = seq(0.01, 0.99, by = 0.01)
##   )
##   crps_ens_fcst <- EnsCrps(ens_fcst$predictions, test_data[["Q"]])
##   df[["crps_ens_fcst"]] <- crps_ens_fcst

##   ## Loop over training data to get climatology for each station
##   station_ids <- test_data$id %>% unique()
##   season <- test_data$season %>% unique()
##   if (length(season) != 1) {
##     stop("test data does not contain unique season")
##   }
##   df_list <- list()
##   for (k in 1:length(station_ids)) {
##     df_stn <- df %>%
##       filter(id %in% station_ids[k])
##     ## Get climatology for specific season
##     train_data_stn <- train_data %>%
##       filter(id %in% station_ids[k]) %>%
##       filter(season %in% !!season)
##     test_data_stn <- test_data %>%
##       filter(id %in% station_ids[k])
##     obs <- test_data_stn[["Q"]]
##     ens_climatology <- train_data_stn[["Q"]]
##     ens_climatology_mat <- t(matrix(
##       rep(ens_climatology, length(obs)), ncol = length(obs)
##     ))
##     crps_ens_climatology <- EnsCrps(ens_climatology_mat, obs)
##     df_stn[["crps_ens_climatology"]] <- crps_ens_climatology
##     df_list[[k]] <- df_stn
##   }
##   df <- do.call("rbind", df_list)

##   ## Add climatological prediction and model name
##   climatology <- train_data %>%
##     group_by(id) %>%
##     summarize(Q50_climatology = mean(Q))
##   df <- df %>% left_join(climatology, by = "id")
##   df
## }

## make_xgb_forecast <- function(formula,
##                               train_data,
##                               test_data,
##                               quantiles = c(0.05, 0.25, 0.5, 0.75, 0.95),
##                               ...) {

##   ## First we need to convert train_data to dense matrix
##   vars <- all.vars(as.formula(formula))
##   dep_var <- vars[1]
##   indep_vars <- vars[!vars %in% dep_var]
##   label <- train_data[[dep_var]] %>% as.numeric()
##   train_data_matrix <- train_data %>% dplyr::select(indep_vars) %>% as.matrix()
##   test_data_matrix <- test_data %>% dplyr::select(indep_vars) %>% as.matrix()

##   ## Train XGBoost model
##   ## TODO better understand these parameters
##   ## TODO hyperparameter optimization - use hyperopt?
##   ## TODO provide validation dataset
##   model <- xgboost(
##     data = train_data_matrix,
##     label = label,
##     eta = 0.3,
##     subsample = 0.5,
##     max_depth = 6,
##     ## eta = 1,
##     nrounds = 500,
##     ## nthread = 2,
##     objective = "reg:squarederror"
##   )
##   pred <- predict(model, newdata = test_data_matrix)
##   df <- test_data %>%
##     dplyr::select(id, season_year, Q) %>%
##     rename(Q_obs = Q) %>%
##     mutate(Q50 = pred)
##   df
## }

## make_forecast <- function(formula,
##                           train_data,
##                           test_data,
##                           test_period,
##                           test_season,
##                           model = "QRF",
##                           method = "forward",
##                           k = NA,
##                           multisite = FALSE, ...) {

##   if (method == "kfold") {
##     ## divide test period into k splits
##     flds <- createFolds(test_period, k, list = TRUE, returnTrain = FALSE)
##     test_period <- lapply(flds, FUN=function(index) test_period[index])
##   }
##   test_train_split <- function(x, test_period, index) {
##     if (method == "forward") {
##       test_period <- test_period[index]
##       train_period <- x[x < test_period]
##     } else if (method == "kfold") {
##       test_period <- test_period[[index]]
##       train_period <- x[!x %in% test_period]
##     } else if (method == "loocv") {
##       test_period <- test_period[index]
##       train_period <- x[!x %in% test_period]
##     }
##     list(train = train_period, test = test_period)
##   }

##   all_years <- train_data$season_year %>%
##     unique() %>%
##     sort()
##   prediction_list <- list()
##   for (i in 1:length(test_period)) {
##     split <- test_train_split(all_years, test_period, i)
##     train_split <- split$train
##     test_split <- split$test
##     train_data_i <- train_data %>%
##       filter(season_year %in% train_split)
##     test_data_i <- test_data %>%
##       filter(season_year %in% test_split) %>%
##       filter(season %in% test_season)

##     ## Handle case where no training/test data available
##     if (nrow(train_data_i) == 0 | nrow(test_data_i) == 0)
##       next

##     ## Add Q_mean from the training set
##     qmean <- train_data_i %>%
##       group_by(id) %>%
##       summarize(qmean = mean(Q_mean))
##     train_data_i <- train_data_i %>%
##       left_join(qmean, by = "id")
##     test_data_i <- test_data_i %>%
##       left_join(qmean, by = "id")

##     ## Make prediction for current time period
##     if (model == "QRF") {
##       prediction <- make_qrf_forecast(
##         formula, train_data_i, test_data_i,
##         quantiles = c(0.05, 0.25, 0.5, 0.75, 0.95)
##       )
##       prediction_list[[length(prediction_list) + 1]] <- prediction
##     } else if (model == "XGB") {
##       prediction <- make_xgb_forecast(
##         formula, train_data_i, test_data_i
##       )
##     }
##   }
##   out <- do.call("rbind", prediction_list)
##   out
## }

## ## TODO
## run_kfold_forecast <- function(...)

## run_loocv_forecast <- function(formula,
##                                train_data,
##                                test_data,
##                                test_period,
##                                test_season,
##                                model = "QRF",
##                                multisite = FALSE, ...) {

##   prediction_list <- list()
##   for (i in 1:length(test_period)) {
##     test_year <- test_period[i]
##     ## Remove test year from training data
##     train_data_i <- train_data %>% filter(!season_year %in% test_year)
##     test_data_i <- test_data %>%
##       filter(season_year %in% test_year) %>%
##       filter(season %in% test_season)

##     ## Handle case where no training/test data available
##     if (nrow(train_data_i) == 0 | nrow(test_data_i) == 0)
##       next


##   }
## }
