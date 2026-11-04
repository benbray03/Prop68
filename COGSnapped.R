library(readr)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(maps)
library(sf)
library(patchwork)   # install.packages("patchwork") if needed

# ---------------------------------------------------------------------
# 0. Data
# ---------------------------------------------------------------------
dat_fish    <- readRDS("data/CPFV_PA_sp_decade_roms_n_ensemble_shiny.rds")
dat_inverts <- readRDS("data/Inverts_PA_sp_decade_roms_n_ensemble_shiny.rds")

decade_order <- c("1995_2000","2001_2010","2011_2020","2021_2030","2031_2040",
                  "2041-2050","2051_2060","2061_2070","2071_2080","2081_2090",
                  "2091_2100")

CRS_PROJ <- 3310   # California Albers (metres)

prepare_dat <- function(raw) {
  raw %>%
    mutate(
      lat     = as.numeric(sub("-.*", "", cell_coord_id)),
      lon     = -as.numeric(sub(".*-", "", cell_coord_id)),
      species = as.character(species),
      decade  = as.character(decade)
    ) %>%
    group_by(cell_coord_id, decade, species, lat, lon) %>%
    summarise(pa_decade_mean_pa = mean(pa_decade_mean_pa, na.rm = TRUE),
              .groups = "drop")
}

dat_fish    <- prepare_dat(dat_fish)
dat_inverts <- prepare_dat(dat_inverts)

build_coast_ref <- function(data) {
  data %>%
    distinct(lat, lon) %>%
    group_by(lat) %>%
    summarise(coast_lon = max(lon, na.rm = TRUE), .groups = "drop") %>%
    arrange(lat)
}

coast_ref_fish    <- build_coast_ref(dat_fish)
coast_ref_inverts <- build_coast_ref(dat_inverts)

# CHECK: coastal reference over the grid
plot_coast_ref_check <- function(data, coast_ref, title_label) {
  ggplot() +
    geom_point(data = distinct(data, lat, lon), aes(x = lon, y = lat),
               size = 0.6, alpha = 0.4, colour = "grey60") +
    geom_path(data = coast_ref, aes(x = coast_lon, y = lat),
              colour = "firebrick", linewidth = 0.8) +
    coord_fixed() +
    labs(title = paste("Coastal reference over grid —", title_label),
         x = "Longitude", y = "Latitude") +
    theme_minimal(base_size = 13)
}

print(plot_coast_ref_check(dat_inverts, coast_ref_inverts, "Inverts grid"))

# ---------------------------------------------------------------------
# 2. COG with land check + snap to nearest coastal reference
# ---------------------------------------------------------------------
compute_cog_snapped <- function(data, coast_ref, crs_proj = CRS_PROJ) {
  
  # Plain weighted lat/lon COG
  cog <- data %>%
    filter(!is.na(pa_decade_mean_pa), pa_decade_mean_pa > 0) %>%
    group_by(species, decade) %>%
    summarise(
      sum_weight  = sum(pa_decade_mean_pa),
      cog_lat_raw = sum(lat * pa_decade_mean_pa) / sum_weight,
      cog_lon_raw = sum(lon * pa_decade_mean_pa) / sum_weight,
      n_cells     = n(),
      .groups = "drop"
    )
  
  # Coastal reference points in metres
  coast_ref <- coast_ref %>% arrange(lat)
  coast_xy  <- st_as_sf(coast_ref, coords = c("coast_lon", "lat"), crs = 4326) %>%
    st_transform(crs_proj) %>%
    st_coordinates()
  
  # Defaults: keep the raw COG
  cog$cog_lat      <- cog$cog_lat_raw
  cog$cog_lon      <- cog$cog_lon_raw
  cog$on_land      <- FALSE
  cog$snap_dist_km <- 0
  
  for (i in seq_len(nrow(cog))) {
    
    # Coastline longitude at this COG's latitude (interpolated between rows)
    coast_lon_here <- approx(coast_ref$lat, coast_ref$coast_lon,
                             xout = cog$cog_lat_raw[i], rule = 2)$y
    
    if (cog$cog_lon_raw[i] > coast_lon_here) {
      # East of the shoreward grid edge -> on land -> snap
      p <- st_as_sf(data.frame(lon = cog$cog_lon_raw[i], lat = cog$cog_lat_raw[i]),
                    coords = c("lon", "lat"), crs = 4326) %>%
        st_transform(crs_proj) %>%
        st_coordinates()
      
      d <- sqrt((coast_xy[, 1] - p[1])^2 + (coast_xy[, 2] - p[2])^2)
      j <- which.min(d)
      
      cog$cog_lat[i]      <- coast_ref$lat[j]
      cog$cog_lon[i]      <- coast_ref$coast_lon[j]
      cog$on_land[i]      <- TRUE
      cog$snap_dist_km[i] <- d[j] / 1000
    } else {
      # Already in the water -> leave as is
      next
    }
  }
  
  cog %>%
    mutate(decade = factor(decade, levels = decade_order)) %>%
    arrange(species, decade)
}

cog_fish    <- compute_cog_snapped(dat_fish,    coast_ref_fish)
cog_inverts <- compute_cog_snapped(dat_inverts, coast_ref_inverts)

# Diagnostic: how many COGs were moved, and how far
bind_rows(cog_fish    %>% mutate(group = "Fish"),
          cog_inverts %>% mutate(group = "Inverts")) %>%
  group_by(group) %>%
  summarise(n_cog       = n(),
            n_snapped   = sum(on_land),
            max_snap_km = max(snap_dist_km),
            .groups = "drop") %>%
  print()

# ---------------------------------------------------------------------
# 3. Single-species check
#    (switch cog_fish <-> cog_inverts depending on the species)
# ---------------------------------------------------------------------
cog_species_test <- cog_fish %>%
  filter(species == "California_halibut")

cog_species_test_long <- cog_species_test %>%
  select(decade, cog_lat, cog_lon) %>%
  pivot_longer(cols = c(cog_lat, cog_lon),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric, cog_lat = "Latitude", cog_lon = "Longitude"))

ggplot(cog_species_test_long, aes(x = decade, y = value, group = metric)) +
  geom_line(colour = "steelblue", linewidth = 1) +
  geom_point(colour = "steelblue", size = 2) +
  facet_wrap(~ metric, scales = "free_y", ncol = 1) +
  labs(title = "Center of Gravity Over Time", x = "Decade", y = NULL) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text  = element_text(face = "bold"))

ca_coast <- map_data("state", region = "california")

ggplot() +
  geom_polygon(data = ca_coast, aes(x = long, y = lat, group = group),
               fill = "grey90", colour = "grey60") +
  geom_path(data = cog_species_test, aes(x = cog_lon, y = cog_lat),
            colour = "steelblue", linewidth = 0.6, alpha = 0.6) +
  geom_point(data = cog_species_test, aes(x = cog_lon, y = cog_lat, colour = decade),
             size = 3) +
  geom_point(data = filter(cog_species_test, on_land),
             aes(x = cog_lon, y = cog_lat),
             shape = 21, size = 5, colour = "red", stroke = 1) +
  coord_quickmap(xlim = range(cog_species_test$cog_lon) + c(-1, 1),
                 ylim = range(cog_species_test$cog_lat) + c(-1, 1)) +
  scale_colour_viridis_d(name = "Decade") +
  labs(title = "Center of Gravity by Decade",
       subtitle = "Red ring = moved off land",
       x = "Longitude", y = "Latitude") +
  theme_minimal(base_size = 13)

# =====================================================================
# 4. Loop over every species
# =====================================================================
out_dir <- "figures/cog_by_species"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cog_all <- bind_rows(
  cog_fish    %>% mutate(group = "Fish"),
  cog_inverts %>% mutate(group = "Inverts")
) %>%
  arrange(group, species, decade)

metric_cols <- c("Latitude" = "#1b6ca8", "Longitude" = "#e07a1f")

# z-score that doesn't blow up when a metric never changes
safe_z <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

# Map window at least 2 deg wide/tall so small shifts still show context
pad_range <- function(x, min_span = 2) {
  r <- range(x, na.rm = TRUE)
  mid <- mean(r); half <- max(diff(r) / 2 + 0.5, min_span / 2)
  c(mid - half, mid + half)
}

plot_species_trends <- function(df, sp, grp) {
  sp_label <- gsub("_", " ", sp)
  
  long <- df %>%
    select(decade, Latitude = cog_lat, Longitude = cog_lon) %>%
    pivot_longer(-decade, names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(metric, levels = names(metric_cols))) %>%
    group_by(metric) %>%
    mutate(z = safe_z(value)) %>%
    ungroup()
  
  # Top left: standardised lat & lon overlaid
  p_overlay <- ggplot(long, aes(x = decade, y = z, colour = metric, group = metric)) +
    geom_hline(yintercept = 0, colour = "grey75", linetype = "dashed") +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_colour_manual(values = metric_cols, name = NULL) +
    labs(title = paste0(sp_label, "  (", grp, ")"),
         subtitle = "Standardised trends (z-score within each metric)",
         x = NULL, y = "Standardised value") +
    theme_minimal(base_size = 13) +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title  = element_text(face = "bold.italic"))
  
  # Bottom: raw values
  p_raw <- ggplot(long, aes(x = decade, y = value, colour = metric, group = metric)) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.8) +
    facet_wrap(~ metric, scales = "free_y", nrow = 1) +
    scale_colour_manual(values = metric_cols, guide = "none") +
    labs(subtitle = "Raw values", x = "Decade", y = NULL) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text  = element_text(face = "bold"))
  
  # Top right: map of snapped vs raw COG track
  map_df <- bind_rows(
    df %>% transmute(decade, on_land, lon = cog_lon,     lat = cog_lat,     method = "Snapped"),
    df %>% transmute(decade, on_land, lon = cog_lon_raw, lat = cog_lat_raw, method = "Raw lat/lon")
  ) %>%
    mutate(method = factor(method, levels = c("Snapped", "Raw lat/lon")))
  
  p_map <- ggplot() +
    geom_polygon(data = ca_coast, aes(x = long, y = lat, group = group),
                 fill = "grey90", colour = "grey60") +
    geom_path(data = map_df, aes(x = lon, y = lat, linetype = method),
              colour = "grey30", linewidth = 0.6, alpha = 0.7) +
    geom_point(data = map_df, aes(x = lon, y = lat, colour = decade, shape = method),
               size = 3) +
    geom_point(data = filter(map_df, method == "Snapped", on_land),
               aes(x = lon, y = lat), shape = 21, size = 5,
               colour = "red", stroke = 1) +
    coord_quickmap(xlim = pad_range(map_df$lon), ylim = pad_range(map_df$lat)) +
    scale_colour_viridis_d(name = "Decade", drop = FALSE) +
    scale_shape_manual(values = c("Snapped" = 16, "Raw lat/lon" = 17), name = "COG") +
    scale_linetype_manual(values = c("Snapped" = "solid", "Raw lat/lon" = "dashed"),
                          name = "COG") +
    labs(subtitle = "COG by decade (red ring = moved off land)",
         x = "Longitude", y = "Latitude") +
    theme_minimal(base_size = 12) +
    theme(legend.text = element_text(size = 8),
          legend.key.size = unit(0.4, "cm"))
  
  ((p_overlay | p_map) + plot_layout(widths = c(1.6, 1))) / p_raw +
    plot_layout(heights = c(1.3, 1))
}

species_keys <- cog_all %>% distinct(group, species)

plots <- map2(species_keys$group, species_keys$species, function(grp, sp) {
  df <- cog_all %>% filter(group == grp, species == sp) %>% arrange(decade)
  plot_species_trends(df, sp, grp)
})
names(plots) <- paste(species_keys$group, species_keys$species, sep = "_")

# One multi-page PDF ...
pdf(file.path(out_dir, "all_species_cog_trends.pdf"), width = 14, height = 9)
walk(plots, print)
dev.off()

# ... and one PNG per species
iwalk(plots, ~ ggsave(file.path(out_dir, paste0(.y, ".png")), .x,
                      width = 14, height = 9, dpi = 200, bg = "white"))

# ---------------------------------------------------------------------
# 5. Summary: net change first -> last decade, plus snapping diagnostics
# ---------------------------------------------------------------------
cog_change <- cog_all %>%
  group_by(group, species) %>%
  summarise(
    d_lat       = last(cog_lat) - first(cog_lat),
    d_lon       = last(cog_lon) - first(cog_lon),
    n_snapped   = sum(on_land),
    max_snap_km = max(snap_dist_km),
    .groups = "drop"
  )

write_csv(cog_change, file.path(out_dir, "cog_net_change_by_species.csv"))
