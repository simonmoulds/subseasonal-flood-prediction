#!/usr/bin/env Rscript

if (exists("snakemake")) {
  season <- snakemake@wildcards[["season"]]
  input_filename <- snakemake@input[[1]]
  output_filename <- snakemake@output[[1]]
  ens_output_filename <- snakemake@output[[2]]
  snakemake@source("utils.R")
} else {
  season <- "DJF"
  input_filename <- "results/preprocessing/EFAS/10003.parquet"
  output_filename <- ""
  source("workflow/scripts/utils.R")
}

season_index <- get_month_index(season)

x <- open_dataset(input_filename) %>% collect()

## Restrict to months in the current season
x <- x %>%
  mutate(month = lubridate::month(time) %>% unname()) %>%
  filter(month %in% season_index) %>%
  mutate(latest_lead_time = interval(init_time, time) %/% months(1)) %>%
  filter(latest_lead_time <= 5)

## Only include full seasons
x <- x %>%
  group_by(init_time, member) %>%
  mutate(n = length(unique(month))) %>%
  filter(n == length(season_index)) %>% ungroup()

## `start_time` is the start of the season
## `latest_lead_time` is the latest lead time, in months, used in the forecast
x <- x %>%
  group_by(id, init_time, member) %>%
  mutate(start_time = min(time)) %>%
  mutate(start_time = lubridate::floor_date(start_time, unit = "months")) %>%
  mutate(latest_lead_time = as.integer(max(latest_lead_time)))

x_seas <- x %>%
  summary_stats_fun("value", c("id", "member", "start_time", "init_time", "latest_lead_time", "variable"), probs = c()) %>%
  unite(variable, variable, summary_statistic) %>%
  mutate(
    year = lubridate::year(start_time) %>% unname(),
    season = as.character(season), .after = start_time
  )

write_parquet(x_seas, output_filename)

## Now compute ensemble statistics
probs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
x_seas_ens <- x_seas %>%
  summary_stats_fun("value", c("id", "start_time", "year", "season", "init_time", "latest_lead_time", "variable"), probs = probs) %>%
  unite(variable, variable, summary_statistic)

write_parquet(x_seas_ens, ens_output_filename)
