library(readr)
library(tidyverse)
library(dplyr)

dat_fish   <- readRDS("data/CPFV_PA_sp_decade_roms_n_ensemble_shiny.rds")
dat_inverts <- readRDS("data/Inverts_PA_sp_decade_roms_n_ensemble_shiny.rds")

assign_decade <- function(year) {
  dplyr::case_when(
    year >= 1995 & year <= 2000 ~ "1995_2000",
    year >= 2001 & year <= 2010 ~ "2001_2010",
    year >= 2011 & year <= 2020 ~ "2011_2020",
    year >= 2021 & year <= 2030 ~ "2021_2030",
    year >= 2031 & year <= 2040 ~ "2031_2040",
    year >= 2041 & year <= 2050 ~ "2041_2050",
    year >= 2051 & year <= 2060 ~ "2051_2060",
    year >= 2061 & year <= 2070 ~ "2061_2070",
    year >= 2071 & year <= 2080 ~ "2071_2080",
    year >= 2081 & year <= 2090 ~ "2081_2090",
    year >= 2091 & year <= 2100 ~ "2091_2100",
    TRUE ~ NA_character_
  )
}

decade_order <- c("1995_2000","2001_2010","2011_2020","2021_2030","2031_2040",
                  "2041_2050","2051_2060","2061_2070","2071_2080","2081_2090",
                  "2091_2100")

prepare_dat <- function(raw) {
  raw <- raw %>%
    mutate(
      lat = as.numeric(sub("-.*", "", cell_coord_id)),
      lon = -as.numeric(sub(".*-", "", cell_coord_id)),
      species = as.character(species),
      decade  = as.character(decade)
    ) %>%
    group_by(cell_coord_id, decade, species, lat, lon) %>%
    summarise(pa_decade_mean_pa = mean(pa_decade_mean_pa, na.rm = TRUE),
              .groups = "drop")
  raw
}

dat_fish    <- prepare_dat(dat_fish)
dat_inverts <- prepare_dat(dat_inverts)

compute_cog <- function(data) {
  lat_range <- range(data$lat, na.rm = TRUE)
  lon_range <- range(data$lon, na.rm = TRUE)
  
  data %>%
    filter(!is.na(pa_decade_mean_pa), pa_decade_mean_pa > 0) %>%
    group_by(species, decade) %>%
    summarise(
      # Sum of moments (position × weight) over sum of weights
      sum_moment_lat = sum(lat * pa_decade_mean_pa, na.rm = TRUE),
      sum_moment_lon = sum(lon * pa_decade_mean_pa, na.rm = TRUE),
      sum_weight     = sum(pa_decade_mean_pa, na.rm = TRUE),
      
      cog_lat = sum_moment_lat / sum_weight,
      cog_lon = sum_moment_lon / sum_weight,
      
      n_cells      = n(),           # diagnostic: how many cells contributed
      total_weight = sum_weight,    # diagnostic: total suitability mass
      .groups = "drop"
    ) %>%
    mutate(
      cog_lat = pmin(pmax(cog_lat, lat_range[1]), lat_range[2]),
      cog_lon = pmin(pmax(cog_lon, lon_range[1]), lon_range[2]),
      decade  = factor(decade, levels = decade_order)
    ) %>%
    arrange(species, decade)
}

cog_fish    <- compute_cog(dat_fish)
cog_inverts <- compute_cog(dat_inverts)
