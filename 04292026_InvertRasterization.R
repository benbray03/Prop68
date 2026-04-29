library(terra)
library(readr)
library(tidyverse)
library(raster)

crs_str <- "+proj=longlat +datum=WGS84"

data <- read_rds("data/Inverts_PA_sp_decade_roms_n_ensemble_shiny.rds") %>%
  mutate(
    lat     = as.numeric(sub("-.*", "", cell_coord_id)),
    lon     = -as.numeric(sub(".*-", "", cell_coord_id)),
    species = as.character(species),
    decade  = as.character(decade)
  )

# Infer resolution from full dataset
res_x <- min(diff(sort(unique(data$lon))))
res_y <- min(diff(sort(unique(data$lat))))

# Average across roms_model and model_predictors (or filter first if preferred)
data_avg <- data %>%
  group_by(species, decade, lat, lon) %>%
  summarise(pa = mean(pa_decade_mean_pa, na.rm = TRUE), .groups = "drop")

# Get all combinations
combos <- data_avg %>% distinct(species, decade)
nrow(combos)  # should be 77

# Build raster template once (assumes same spatial extent across all combos)
r_template <- rast(
  xmin = min(data_avg$lon) - res_x / 2,
  xmax = max(data_avg$lon) + res_x / 2,
  ymin = min(data_avg$lat) - res_y / 2,
  ymax = max(data_avg$lat) + res_y / 2,
  resolution = c(res_x, res_y),
  crs = crs_str
)

# --- Loop ---
raster_list <- list()

for (i in seq_len(nrow(combos))) {
  sp  <- combos$species[i]
  dec <- combos$decade[i]
  
  sub <- data_avg %>% filter(species == sp, decade == dec)
  
  pts <- vect(sub, geom = c("lon", "lat"), crs = crs_str)
  r   <- rasterize(pts, r_template, field = "pa", fun = "mean")
  
  key <- paste0(sp, "_", dec)
  names(r) <- key
  raster_list[[key]] <- r
  
  message(i, "/", nrow(combos), " done: ", key)
}

raster_list_r <- lapply(raster_list, function(x) {
  raster::raster(as.matrix(x, wide = TRUE),
                 xmn = terra::xmin(x), xmx = terra::xmax(x),
                 ymn = terra::ymin(x), ymx = terra::ymax(x),
                 crs = raster::crs("+proj=longlat +datum=WGS84"))
})
saveRDS(raster_list_r, "data/invert_raster_list.rds")

