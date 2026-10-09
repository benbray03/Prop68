library(readr)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(maps)

dat_fish    <- readRDS("data/CPFV_PA_sp_decade_roms_n_ensemble_shiny.rds")
dat_inverts <- readRDS("data/Inverts_PA_sp_decade_roms_n_ensemble_shiny.rds")

decade_order <- c("1995_2000","2001_2010","2011_2020","2021_2030","2031_2040",
                  "2041_2050","2051_2060","2061_2070","2071_2080","2081_2090",
                  "2091_2100")

prepare_dat <- function(raw) {
  raw %>%
    mutate(
      lat = as.numeric(sub("-.*", "", cell_coord_id)),
      lon = -as.numeric(sub(".*-", "", cell_coord_id)),
      species = as.character(species),
      decade  = as.character(decade)
    ) %>%
    group_by(cell_coord_id, decade, species, lat, lon) %>%
    summarise(pa_decade_mean_pa = mean(pa_decade_mean_pa, na.rm = TRUE),
              .groups = "drop")
}

dat_fish    <- prepare_dat(dat_fish)
dat_inverts <- prepare_dat(dat_inverts)

# ---------------------------------------------------------------------
# STEP 1 — find the westernmost point of the coast present in the grid
# ---------------------------------------------------------------------
# For each latitude row, the westernmost grid cell present is used as
# that row's "coastal reference" longitude. We also keep the single
# most-western point across the whole grid (e.g. near Pt. Conception) —
# that's the point every row's reference gets lined up against, so the
# shifted coordinates stay in real, interpretable degrees.
build_coast_ref <- function(data) {
  ref <- data %>%
    distinct(lat, lon) %>%
    group_by(lat) %>%
    summarise(west_lon = min(lon, na.rm = TRUE), .groups = "drop") %>%
    arrange(lat)
  
  list(
    ref = ref,
    fun = approxfun(ref$lat, ref$west_lon, rule = 2),  # interpolates between grid lats
    global_west_lon = min(ref$west_lon, na.rm = TRUE)
  )
}

# ---------------------------------------------------------------------
# STEP 2 — shift each row's longitudes so its coastal point lines up
# with the single westernmost coastal point
# ---------------------------------------------------------------------
shift_to_coast <- function(data, coast_ref) {
  data %>%
    mutate(
      row_west_lon = coast_ref$fun(lat),
      shifted_lon  = lon - row_west_lon + coast_ref$global_west_lon
    )
}

# ---------------------------------------------------------------------
# STEP 3 — CHECK: visualize every grid square, before vs after shifting
# ---------------------------------------------------------------------
plot_shift_check <- function(data, coast_ref, title_label) {
  grid_check <- data %>%
    distinct(lat, lon) %>%
    shift_to_coast(coast_ref)
  
  grid_check_long <- bind_rows(
    grid_check %>% transmute(lat, lon, panel = "Original"),
    grid_check %>% transmute(lat, lon = shifted_lon, panel = "Shifted")
  )
  
  ggplot(grid_check_long, aes(x = lon, y = lat, color = lat)) +
    geom_point(size = 0.7, alpha = 0.6) +
    facet_wrap(~panel, scales = "free_x") +
    scale_color_viridis_c(guide = "none") +
    labs(title = paste("Grid cells before vs after coastline straightening —", title_label),
         x = "Longitude", y = "Latitude") +
    theme_minimal(base_size = 13)
}

# ---------------------------------------------------------------------
# STEP 4 — compute center of gravity in the shifted (straightened) frame
# STEP 5 — shift the COG back using that specific cog_lat's transform
# ---------------------------------------------------------------------
compute_cog <- function(data, coast_ref) {
  lat_range <- range(data$lat, na.rm = TRUE)
  lon_range <- range(data$lon, na.rm = TRUE)
  
  data %>%
    filter(!is.na(pa_decade_mean_pa), pa_decade_mean_pa > 0) %>%
    shift_to_coast(coast_ref) %>%
    group_by(species, decade) %>%
    summarise(
      sum_weight      = sum(pa_decade_mean_pa, na.rm = TRUE),
      cog_lat         = sum(lat * pa_decade_mean_pa, na.rm = TRUE) / sum_weight,
      cog_shifted_lon = sum(shifted_lon * pa_decade_mean_pa, na.rm = TRUE) / sum_weight,
      n_cells         = n(),
      total_weight    = sum_weight,
      .groups = "drop"
    ) %>%
    mutate(
      cog_lat = pmin(pmax(cog_lat, lat_range[1]), lat_range[2]),
      # shift back OUT of straightened space, using the transform that
      # applies at this specific cog_lat (not a fixed row)
      cog_lon = cog_shifted_lon - coast_ref$global_west_lon + coast_ref$fun(cog_lat),
      cog_lon = pmin(pmax(cog_lon, lon_range[1]), lon_range[2]),
      decade  = factor(decade, levels = decade_order)
    ) %>%
    arrange(species, decade)
}

# ---------------------------------------------------------------------
# Run it
# ---------------------------------------------------------------------
coast_ref_fish    <- build_coast_ref(dat_fish)
coast_ref_inverts <- build_coast_ref(dat_inverts)

# STEP 3 check — look at this BEFORE trusting the COG below
print(plot_shift_check(dat_inverts, coast_ref_inverts, "Inverts grid"))

cog_fish    <- compute_cog(dat_fish, coast_ref_fish)
cog_inverts <- compute_cog(dat_inverts, coast_ref_inverts)

# need to change inverts to fish vice versa based on the species
cog_species_test <- cog_fish %>%
  filter(species == "California_halibut")

cog_species_test_long <- cog_species_test %>%
  select(decade, cog_lat, cog_lon) %>%
  pivot_longer(cols = c(cog_lat, cog_lon),
               names_to = "metric",
               values_to = "value") %>%
  mutate(metric = recode(metric, cog_lat = "Latitude", cog_lon = "Longitude"))

ggplot(cog_species_test_long, aes(x = decade, y = value, group = metric)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(color = "steelblue", size = 2) +
  facet_wrap(~ metric, scales = "free_y", ncol = 1) +
  labs(title = "Center of Gravity Over Time", x = "Decade", y = NULL) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold"))

# Map of COG
ca_coast <- map_data("state", region = "california")

ggplot() +
  geom_polygon(data = ca_coast, aes(x = long, y = lat, group = group),
               fill = "grey90", color = "grey60") +
  geom_path(data = cog_species_test, aes(x = cog_lon, y = cog_lat),
            color = "steelblue", linewidth = 0.6, alpha = 0.6) +
  geom_point(data = cog_species_test, aes(x = cog_lon, y = cog_lat, color = decade),
             size = 3) +
  coord_quickmap(xlim = range(cog_species_test$cog_lon) + c(-1, 1),
                 ylim = range(cog_species_test$cog_lat) + c(-1, 1)) +
  scale_color_viridis_d(name = "Decade") +
  labs(title = "Center of Gravity by Decade", x = "Longitude", y = "Latitude") +
  theme_minimal(base_size = 13)
