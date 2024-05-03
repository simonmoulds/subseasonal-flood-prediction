
## library(pacman)
## ## Load packages with pacman
## p_load(tidyverse, arrow, lubridate, zoo, yaml)
## p_load(sf, terra, exactextractr)
library(tidyverse)
library(arrow)
library(lubridate)
library(zoo)
library(yaml)
library(sf)
library(terra)
library(exactextractr)

## Increase cache size
terra::gdalCache(4096)


## get_catchment_poly <- function(path, ids, ...) {
##   stations_str <- paste0(ids, collapse=", ")
##   query <- paste0(
##     "SELECT * FROM catchment_poly WHERE ID IN (", stations_str, ")"
##   )
##   catchment_poly <- st_read(path, query = query)
##   catchment_poly
## }


parse_c3s_filename <- function(f) {
  f <- basename(f) # Remove file path
  f <- sub("(.*)\\..*$", "\\1", f) # Remove file extension
  parts <- str_split(f, "_")[[1]]
  institution <- parts[1]
  source_id <- parts[2]
  experiment <- parts[3]
  init_time <- gsub("S", "", parts[4])
  init_time <- as.Date(init_time, format = "%Y%m%d")
  variable <- parts[5]
  product_type <- gsub("^.*(monthly_mean|maximum|minimum).*$", "\\1", f)
  list(institution = institution,
       source_id = source_id,
       experiment = experiment,
       init_time = init_time,
       variable = variable,
       product_type = product_type)
}


construct_pattern <- function(start_year, end_year, system, variable) {
  years_regex <- paste0(seq(start_year, end_year), collapse = "|")
  ptn <- paste0(
    "(.*)_", system,
    "_hindcast_S(", years_regex, ")[0-9]{4}_",
    variable, "_monthly_(minimum|maximum|mean)", ".nc"
  )
  return(ptn)
}


extract_c3s_daily_catchment_data <- function(f,
                                             variable,
                                             polygons,
                                             polygons_id_column, ...) {

  ## TESTING
  ## TODO test whether applying remapcon is beneficial - could use to crop region
  ## cdo remapcon,gridfile.txt infile outfile
  ## https://www.sciencedirect.com/science/article/pii/S0022169423006352#s0010
  ## f <- "data/C3S_hindcast_daily_SEAS5//ecmf_SEAS5-v20171101_hindcast_S20100101_tp_daily.nc"
  ## variable <- "tp"
  ## polygons <- st_read("data/Catchment_Boundaries/CAMELS_GB_catchment_boundaries.shp")
  ## polygons_id_column <- "ID"

  metadata <- parse_c3s_filename(f)
  system <- paste0(metadata$institution, "_", metadata$source_id)

  r <- terra::rast(f)
  ## Convert longitudes from 0-360 to -180 to +180, if needed
  if (xmax(r) + xres(r) > 360) {
    r <- rotate(r)
  }
  crs(r) <- "epsg:4326"

  ## Add time to layer names
  nms <- names(r)
  member <- sub("^.*_number=([0-9]+)_([0-9]+)", "\\1", nms) %>% as.numeric
  tm <- time(r) %>% format("%Y%m%d")
  nms <- paste0(nms, "_", tm)
  names(r) <- nms

  ## Extract data
  x <- suppressWarnings(
    exact_extract(
      r, polygons, fun = "mean",
      append_cols = polygons_id_column, ...
    )
  )
  ## Add index data then pivot to long format
  x <- x |>
    as_tibble() |>
    rename(id = polygons_id_column) |>
    mutate(variable = variable)

  x <- x %>%
    pivot_longer(
      -all_of(c("id", "variable")),
      names_to = "name", values_to = "value"
    )

  ## Get member/lead_time from name
  x <- x %>%
    mutate(name = gsub("^.*_number=", "", name)) %>%
    separate(
      name, c("member", "lead_time", "time"),
      sep = "_", remove = TRUE, convert = FALSE
    ) %>%
    mutate(member = as.integer(member)) %>%
    dplyr::select(-lead_time)

  ## Recalculate lead time (instead of relying on layer names)
  x <- x %>%
    mutate(time = as.Date(time, format = "%Y%m%d")) %>%
    mutate(lead_time = interval(metadata$init_time, time) %/% days(1), .after = time)

  ## Add other metadata
  x <- x %>% mutate(init_time = metadata$init_time, source_id = system)

  ## Give columns a logical order
  x <- x %>%
    dplyr::select(
             id, source_id, member, init_time, lead_time, time,
             variable, value
           )
  x
}

extract_catchment_data <- function(f,
                                   variable,
                                   polygons,
                                   polygons_id_column,
                                   ...) {

  r <- terra::rast(f)
  tm <- time(r) %>% format("%Y%m%d")
  names(r) <- tm
  x <- suppressWarnings(
    exact_extract(
      r, polygons, fun = "mean",
      append_cols = polygons_id_column,
      full_colnames = TRUE, ...
    )
  )
  x <- x |>
    as_tibble() |>
    rename(id = polygons_id_column) |>
    mutate(variable = variable)

  x <- x |>
    pivot_longer(-all_of(c("id", "variable")), names_to="time", values_to="value")

  x <- x |>
    mutate(time = gsub("mean.", "", time)) |>
    mutate(time = as.Date(time, format = "%Y%m%d"))
  x
}

extract_haduk_daily_catchment_data <- function(f,
                                               variable,
                                               polygons,
                                               polygons_id_column,
                                               ...) {

  ## TESTING
  ## f <- "data/Haduk-Grid/v1.2.0.ceda/5km/rainfall/day/v20230328/rainfall_hadukgrid_uk_5km_day_20220101-20220131.nc"
  ## variable <- "rainfall"
  ## polygons <- st_read("data/Catchment_Boundaries/CAMELS_GB_catchment_boundaries.shp")
  ## polygons_id_column <- "ID_STRING"

  ## Extract data for polygons
  x <- extract_catchment_data(f, variable, polygons, polygons_id_column, ...)

  ## Give columns a logical order
  x <- x %>% dplyr::select(id, time, variable, value)
  x
}


## parse_c3s_filename <- function(f) {
##   f <- basename(f) # Remove file path
##   f <- sub("(.*)\\..*$", "\\1", f) # Remove file extension
##   parts <- str_split(f, "_")[[1]]
##   institution <- parts[1]
##   source_id <- parts[2]
##   experiment <- parts[3]
##   init_time <- gsub("S", "", parts[4])
##   init_time <- as.Date(init_time, format = "%Y%m%d")
##   variable <- parts[5]
##   product_type <- gsub("^.*(monthly_mean|maximum|minimum).*$", "\\1", f)
##   list(institution = institution,
##        source_id = source_id,
##        experiment = experiment,
##        init_time = init_time,
##        variable = variable,
##        product_type = product_type)
## }


## extract_c3s_catchment_data <- function(f,
##                                        variable,
##                                        polygons,
##                                        polygons_id_column,
##                                        ...) {
##   metadata <- parse_c3s_filename(f)
##   system <- paste0(metadata$institution, "_", metadata$source_id)

##   r <- terra::rast(f)
##   ## Convert longitudes from 0-360 to -180 to +180, if needed
##   if (xmax(r) + xres(r) > 360) {
##     r <- rotate(r)
##   }
##   crs(r) <- "epsg:4326"

##   ## Add time to layer names
##   nms <- names(r)
##   member <- sub("^.*_number=([0-9]+)_([0-9]+)", "\\1", nms) %>% as.numeric
##   tm <- time(r) %>% format("%Y%m01")
##   nms <- paste0(nms, "_", tm)
##   names(r) <- nms

##   ## Extract data
##   x <- suppressWarnings(
##     exact_extract(
##       r, polygons, fun = "mean",
##       append_cols = polygons_id_column, ...
##     )
##   )
##   ## Add index data then pivot to long format
##   x <- x %>%
##     rename(id = polygons_id_column) %>%
##     mutate(variable = metadata$variable)
##   x <- x %>%
##     pivot_longer(
##       -all_of(c("id", "variable")),
##       names_to = "name", values_to = "value"
##     )
##   ## Get member/lead_time from name
##   x <- x %>%
##     mutate(name = gsub("^.*_number=", "", name)) %>%
##     separate(
##       name, c("member", "lead_time", "time"),
##       sep = "_", remove = TRUE, convert = FALSE
##     ) %>%
##     mutate(member = as.integer(member)) %>%
##     dplyr::select(-lead_time)
##   ## Recalculate lead time (instead of relying on layer names)
##   x <- x %>%
##     mutate(time = as.Date(time, format = "%Y%m%d")) %>%
##     mutate(lead_time = interval(metadata$init_time, time) %/% months(1), .after = time)
##   ## Add other metadata
##   x <- x %>%
##     mutate(
##       year = year(time),
##       month = month(time),
##       days_in_month = lubridate::days_in_month(time) %>% unname(),
##       init_time = metadata$init_time,
##       init_month = month(metadata$init_time),
##       product_type = metadata$product_type,
##       ## source_id = metadata$source_id
##       source_id = system
##     )
##   ## Give columns a logical order
##   x <- x %>%
##     dplyr::select(
##              id, source_id, member, init_time, init_month, lead_time, time, year, month,
##              days_in_month, variable, product_type, value
##            )
##   x
## }

## ## extract_era5_catchment_data <- function(root,
## ##                                         variable,
## ##                                         polygons,
## ##                                         polygons_id_column,
## ##                                         ...) {
## ##   ## Get filename
## ##   fn <- sprintf("era5_reanalysis_%s.nc", variable)
## ##   r <- terra::rast(file.path(root, fn))
## ##   ## Convert longitudes from 0-360 to -180 to +180, if needed
## ##   if (xmax(r) + xres(r) > 360) {
## ##     r <- rotate(r)
## ##   }
## ##   crs(r) <- "epsg:4326"
## ##   ## Extract data for polygons
## ##   x <- extract_catchment_data(
## ##     r, variable, polygons,
## ##     polygons_id_column, progress = FALSE
## ##   )
## ##   ## Add time information
## ##   x <- x %>%
## ##     mutate(time = as.Date(time)) %>%
## ##     mutate(
## ##       year = year(time),
## ##       month = month(time),
## ##       days_in_month = lubridate::days_in_month(time)
## ##     ) %>%
## ##     dplyr::select(time, year, month, days_in_month, id, variable, value)
## ##   x
## ## }

## extract_era5_catchment_data <- function(f,
##                                         variable,
##                                         polygons,
##                                         polygons_id_column,
##                                         ...) {
##   r <- terra::rast(f)
##   if (xmax(r) + xres(r) > 360) {
##     r <- rotate(r)
##   }
##   crs(r) <- "epsg:4326"

##   ## Extract data for polygons
##   tm <- time(r) %>%
##     as.yearmon() %>%
##     format("%Y%m")
##   tm <- paste0(tm, "01")
##   names(r) <- tm
##   x <- suppressWarnings(
##     exact_extract(
##       r, polygons, fun = "mean",
##       append_cols = polygons_id_column,
##       full_colnames = TRUE
##     )
##   )
##   x <- x %>% rename(id = polygons_id_column)
##   x <- x %>% pivot_longer(-id, names_to="time", values_to="value")
##   x <- x %>%
##     mutate(time = gsub("mean.", "", time)) %>%
##     mutate(time = as.Date(time, format = "%Y%m%d"))

##   ## Add time information
##   x <- x %>%
##     mutate(time = as.Date(time)) %>%
##     mutate(
##       year = year(time),
##       month = month(time),
##       days_in_month = lubridate::days_in_month(time) %>% unname()
##     ) %>%
##     mutate(variable = variable, product_type = "monthly_mean") %>%
##     dplyr::select(id, time, year, month, days_in_month, variable, product_type, value)
##   x
## }


## get_haduk_path <- function(root, product, version, res, time_res, variable, ...) {
##   p <- file.path(root, product, version, res, variable, time_res)
##   vns <- list.files(p, pattern="^v[0-9]+") %>% sort(decreasing = TRUE)
##   latest_vn <- vns[1]
##   p <- file.path(p, latest_vn)
##   p
## }


