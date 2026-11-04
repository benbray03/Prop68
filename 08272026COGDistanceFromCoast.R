library(readr)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(maps)
library(sf)

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
# Coastal reference points: for each latitude row, the SHOREWARD
# (eastern / less-negative-longitude) edge of the grid — this is the
# coast itself, as opposed to the offshore boundary.
# ---------------------------------------------------------------------
build_coast_ref <- function(data) {
  data %>%
    distinct(lat, lon) %>%
    group_by(lat) %>%
    summarise(coast_lon = max(lon, na.rm = TRUE), .groups = "drop") %>%
    arrange(lat)
}

CRS_PROJ <- 3310

build_coast_line <- function(coast_ref, full_grid) {
  ref <- coast_ref %>% arrange(lat)   # south -> north
  
  pts_sf <- st_as_sf(ref, coords = c("coast_lon", "lat"), crs = 4326) %>%
    st_transform(CRS_PROJ)
  xy <- st_coordinates(pts_sf)
  
  seg_len <- c(0, sqrt(diff(xy[, 1])^2 + diff(xy[, 2])^2))
  cum_len <- cumsum(seg_len)
  
  n <- nrow(xy)
  tangent <- matrix(NA_real_, n, 2)
  for (i in seq_len(n)) {
    lo <- max(i - 1, 1); hi <- min(i + 1, n)
    d <- xy[hi, ] - xy[lo, ]
    tangent[i, ] <- d / sqrt(sum(d^2))
  }
  normal <- cbind(-tangent[, 2], tangent[, 1])  # rotate 90 deg; sign fixed below
  
  # fix each vertex's normal to point toward a known water cell at that lat
  for (i in seq_len(n)) {
    lat_i <- ref$lat[i]
    water_lon <- min(full_grid$lon[full_grid$lat == lat_i], na.rm = TRUE)
    if (is.finite(water_lon)) {
      wp <- st_as_sf(data.frame(lon = water_lon, lat = lat_i),
                     coords = c("lon", "lat"), crs = 4326) %>%
        st_transform(CRS_PROJ) %>% st_coordinates()
      v <- wp - xy[i, , drop = FALSE]
      if (sum(v * normal[i, ]) < 0) normal[i, ] <- -normal[i, ]
    }
  }
  
  list(lat = ref$lat, xy = xy, cum_len = cum_len,
       tangent = tangent, normal = normal, crs = CRS_PROJ)
}


project_alongshore <- function(data, coast_line) {
  pts_sf <- st_as_sf(data, coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>%
    st_transform(coast_line$crs)
  xy <- st_coordinates(pts_sf)
  
  idx <- match(data$lat, coast_line$lat)
  v <- xy - coast_line$xy[idx, , drop = FALSE]
  
  data %>%
    mutate(
      alongshore_dist = coast_line$cum_len[idx] +
        rowSums(v * coast_line$tangent[idx, , drop = FALSE]),
      crossshore_dist = rowSums(v * coast_line$normal[idx, , drop = FALSE])
    )
}


# CHECK 
plot_coast_line_check <- function(data, coast_line, title_label) {
  coast_ll <- st_as_sf(data.frame(x = coast_line$xy[, 1], y = coast_line$xy[, 2]),
                       coords = c("x", "y"), crs = coast_line$crs) %>%
    st_transform(4326) %>% st_coordinates() %>% as.data.frame()
  names(coast_ll) <- c("lon", "lat")
  
  ggplot() +
    geom_point(data = distinct(data, lat, lon), aes(x = lon, y = lat),
               size = 0.6, alpha = 0.4, color = "grey60") +
    geom_path(data = coast_ll, aes(x = lon, y = lat),
              color = "firebrick", linewidth = 0.8) +
    coord_fixed() +
    labs(title = paste("Coastline curve over grid —", title_label),
         x = "Longitude", y = "Latitude") +
    theme_minimal(base_size = 13)
}

compute_cog_coastline <- function(data, coast_line) {
  data %>%
    filter(!is.na(pa_decade_mean_pa), pa_decade_mean_pa > 0) %>%
    project_alongshore(coast_line) %>%
    group_by(species, decade) %>%
    summarise(
      sum_weight      = sum(pa_decade_mean_pa, na.rm = TRUE),
      cog_alongshore  = sum(alongshore_dist * pa_decade_mean_pa, na.rm = TRUE) / sum_weight,
      cog_crossshore  = sum(crossshore_dist * pa_decade_mean_pa, na.rm = TRUE) / sum_weight,
      n_cells         = n(),
      total_weight    = sum_weight,
      .groups = "drop"
    ) %>%
    mutate(decade = factor(decade, levels = decade_order)) %>%
    arrange(species, decade)
}


unproject_alongshore <- function(cog_df, coast_line) {
  cum_len <- coast_line$cum_len; xy <- coast_line$xy; normal <- coast_line$normal
  
  out_xy <- matrix(NA_real_, nrow(cog_df), 2)
  for (i in seq_len(nrow(cog_df))) {
    d <- pmin(pmax(cog_df$cog_alongshore[i], min(cum_len)), max(cum_len))
    j <- findInterval(d, cum_len, all.inside = TRUE)
    frac <- (d - cum_len[j]) / (cum_len[j + 1] - cum_len[j])
    
    coast_pt <- xy[j, ] + frac * (xy[j + 1, ] - xy[j, ])
    norm_dir <- normal[j, ] + frac * (normal[j + 1, ] - normal[j, ])
    norm_dir <- norm_dir / sqrt(sum(norm_dir^2))
    
    out_xy[i, ] <- coast_pt + cog_df$cog_crossshore[i] * norm_dir
  }
  
  out_ll <- st_as_sf(data.frame(x = out_xy[, 1], y = out_xy[, 2]),
                     coords = c("x", "y"), crs = coast_line$crs) %>%
    st_transform(4326) %>% st_coordinates()
  
  cog_df %>% mutate(cog_lon = out_ll[, 1], cog_lat = out_ll[, 2])
}

coast_ref_fish    <- build_coast_ref(dat_fish)
coast_ref_inverts <- build_coast_ref(dat_inverts)

coast_line_fish    <- build_coast_line(coast_ref_fish,    dat_fish    %>% distinct(lat, lon))
coast_line_inverts <- build_coast_line(coast_ref_inverts, dat_inverts %>% distinct(lat, lon))

# CHECK 
print(plot_coast_line_check(dat_inverts, coast_line_inverts, "Inverts grid"))

cog_fish    <- compute_cog_coastline(dat_fish, coast_line_fish)
cog_inverts <- compute_cog_coastline(dat_inverts, coast_line_inverts)

# need to change inverts to fish vice versa based on the species
cog_species_test <- cog_fish %>%
  filter(species == "California_halibut") %>%
  unproject_alongshore(coast_line_fish)

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


# =====================================================================
# Loop over every species: compare trends in lat / lon COG with
# distance from the coast (cross-shore COG).
# Run this AFTER cog_fish, cog_inverts, coast_line_fish and
# coast_line_inverts have been created.
# =====================================================================
library(patchwork)   # install.packages("patchwork") if needed

out_dir <- "figures/cog_by_species"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------
# 1. Un-project every species at once, each with its own coastline
#    (cross-shore is in metres in EPSG:3310; positive = offshore,
#     because the normals were flipped to point toward water)
# ---------------------------------------------------------------------
cog_all <- bind_rows(
  unproject_alongshore(cog_fish,    coast_line_fish)    %>% mutate(group = "Fish"),
  unproject_alongshore(cog_inverts, coast_line_inverts) %>% mutate(group = "Inverts")
) %>%
  mutate(
    crossshore_km = cog_crossshore / 1000,
    alongshore_km = cog_alongshore / 1000
  ) %>%
  arrange(group, species, decade)

# ---------------------------------------------------------------------
# 1b. Plain lat/lon COG (weighted mean of cell lat & lon, no coastline
#     projection) for comparison on the map
# ---------------------------------------------------------------------
compute_cog_latlon <- function(data) {
  data %>%
    filter(!is.na(pa_decade_mean_pa), pa_decade_mean_pa > 0) %>%
    group_by(species, decade) %>%
    summarise(
      cog_lat_ll = weighted.mean(lat, pa_decade_mean_pa, na.rm = TRUE),
      cog_lon_ll = weighted.mean(lon, pa_decade_mean_pa, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(decade = factor(decade, levels = decade_order))
}

cog_latlon_all <- bind_rows(
  compute_cog_latlon(dat_fish)    %>% mutate(group = "Fish"),
  compute_cog_latlon(dat_inverts) %>% mutate(group = "Inverts")
)

cog_all <- cog_all %>%
  left_join(cog_latlon_all, by = c("group", "species", "decade"))

ca_coast <- map_data("state", region = "california")

metric_cols <- c("Latitude"                 = "#1b6ca8",
                 "Longitude"                = "#e07a1f",
                 "Distance from coast (km)" = "#2a9d55")

# z-score that doesn't blow up when a metric never changes
safe_z <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

# ---------------------------------------------------------------------
# 2. Plot function for one species
#    Top:    all three metrics standardised and overlaid (compare shape)
#    Bottom: each metric in its real units
# ---------------------------------------------------------------------
plot_species_trends <- function(df, sp, grp) {
  sp_label <- gsub("_", " ", sp)
  
  long <- df %>%
    select(decade,
           Latitude                   = cog_lat,
           Longitude                  = cog_lon,
           `Distance from coast (km)` = crossshore_km) %>%
    pivot_longer(-decade, names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(metric, levels = names(metric_cols))) %>%
    group_by(metric) %>%
    mutate(z = safe_z(value)) %>%
    ungroup()
  
  p_overlay <- ggplot(long, aes(x = decade, y = z, colour = metric, group = metric)) +
    geom_hline(yintercept = 0, colour = "grey75", linetype = "dashed") +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_colour_manual(values = metric_cols, name = NULL) +
    labs(title = paste0(sp_label, "  (", grp, ")"),
         subtitle = "Standardised trends (z-score within each metric) — compare direction and timing",
         x = NULL, y = "Standardised value") +
    theme_minimal(base_size = 13) +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold.italic"))
  
  p_raw <- ggplot(long, aes(x = decade, y = value, colour = metric, group = metric)) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.8) +
    facet_wrap(~ metric, scales = "free_y", nrow = 1) +
    scale_colour_manual(values = metric_cols, guide = "none") +
    labs(subtitle = "Raw values", x = "Decade", y = NULL) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(face = "bold"))
  
  # Map of the COG track, coloured by decade
  # (window is at least 2 deg wide/tall so small shifts still show context)
  pad_range <- function(x, min_span = 2) {
    r <- range(x, na.rm = TRUE)
    mid <- mean(r); half <- max(diff(r) / 2 + 0.5, min_span / 2)
    c(mid - half, mid + half)
  }
  
  map_df <- bind_rows(
    df %>% transmute(decade, lon = cog_lon,    lat = cog_lat,    method = "Coastline-based"),
    df %>% transmute(decade, lon = cog_lon_ll, lat = cog_lat_ll, method = "Plain lat/lon")
  ) %>%
    mutate(method = factor(method, levels = c("Coastline-based", "Plain lat/lon")))
  
  p_map <- ggplot() +
    geom_polygon(data = ca_coast, aes(x = long, y = lat, group = group),
                 fill = "grey90", colour = "grey60") +
    geom_path(data = map_df, aes(x = lon, y = lat, linetype = method),
              colour = "grey30", linewidth = 0.6, alpha = 0.7) +
    geom_point(data = map_df, aes(x = lon, y = lat, colour = decade, shape = method),
               size = 3) +
    coord_quickmap(xlim = pad_range(map_df$lon), ylim = pad_range(map_df$lat)) +
    scale_colour_viridis_d(name = "Decade", drop = FALSE) +
    scale_shape_manual(values = c("Coastline-based" = 16, "Plain lat/lon" = 17),
                       name = "COG method") +
    scale_linetype_manual(values = c("Coastline-based" = "solid", "Plain lat/lon" = "dashed"),
                          name = "COG method") +
    labs(subtitle = "COG by decade", x = "Longitude", y = "Latitude") +
    theme_minimal(base_size = 12) +
    theme(legend.text = element_text(size = 8),
          legend.key.size = unit(0.4, "cm"))
  
  ((p_overlay | p_map) + plot_layout(widths = c(1.6, 1))) / p_raw +
    plot_layout(heights = c(1.3, 1))
}

# ---------------------------------------------------------------------
# 3. Loop over all species
# ---------------------------------------------------------------------
species_keys <- cog_all %>% distinct(group, species)

plots <- map2(species_keys$group, species_keys$species, function(grp, sp) {
  df <- cog_all %>% filter(group == grp, species == sp) %>% arrange(decade)
  plot_species_trends(df, sp, grp)
})
names(plots) <- paste(species_keys$group, species_keys$species, sep = "_")

# One multi-page PDF (easy to flip through) ...
pdf(file.path(out_dir, "all_species_cog_trends.pdf"), width = 14, height = 9)
walk(plots, print)
dev.off()

# ... and one PNG per species
iwalk(plots, ~ ggsave(file.path(out_dir, paste0(.y, ".png")), .x,
                      width = 14, height = 9, dpi = 200, bg = "white"))

# ---------------------------------------------------------------------
# 4. Optional: a quick summary table of net change, first -> last decade
# ---------------------------------------------------------------------
cog_change <- cog_all %>%
  group_by(group, species) %>%
  summarise(
    d_lat           = last(cog_lat)       - first(cog_lat),
    d_lon           = last(cog_lon)       - first(cog_lon),
    d_crossshore_km = last(crossshore_km) - first(crossshore_km),
    d_alongshore_km = last(alongshore_km) - first(alongshore_km),
    .groups = "drop"
  )

write_csv(cog_change, file.path(out_dir, "cog_net_change_by_species.csv"))