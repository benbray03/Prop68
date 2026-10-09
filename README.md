# Prop 68: Projected Species Distributions and MPA Coverage

Interactive Shiny app and supporting R scripts for projecting how environmental suitability for California fishery species shifts under climate change (1995–2100), and how much of each species' suitable habitat falls inside Marine Protected Areas (MPAs).

Suitability is modeled with ROMS-downscaled ocean projections from three Earth System Models (GFDL, HadGEM, IPSL) plus a historical run. Maps show the ensemble mean, summarized by decade.

**Live app:** https://casellelab.shinyapps.io/prop80/
**Contact:** Benjamin Bray — [email]

---

## Species

| Group | Species (data name) |
|---|---|
| **Fish** (CPFV) | `Black_rockfish`, `California_barracuda`, `California_halibut`, `California_sheephead`, `Copper_rockfish`, `Gopher_rockfish`, `Ocean_whitefish` |
| **Invertebrates** | `CA_Lobster`, `market_squid`, `pink_shrimp`, `red_urchin`, `Ridgeback_prawn`, `Sea_cucumber` |

Display names for each species are set in `fish_labels` and `invert_labels` in `app.R`.

---

## Pipeline

```
Model output (.rds)  →  Rasterization scripts  →  raster lists (.rds)  →  app.R
                                                                           ↑
                     all_sp_proportion_in_MPA_area_CA.rds  ────────────────┘
```

1. **Model output** (`data/*_PA_sp_decade_roms_n_ensemble_shiny.rds`): point-level decadal suitability for each species, ROMS model, and grid cell. *[Note where these files come from: the script, collaborator, or lab drive that produces them.]*
2. **Rasterization** (`04142026_Rasterization.R` for fish, `04292026_InvertRasterization.R` for inverts): averages across ROMS models, then converts each species × decade combination to a raster. The output is a named list saved as `.rds`.
3. **`app.R`**: loads the raster lists, the point data, and the MPA coverage data, and serves the map and plots.

---

## Repository structure

```
.
├── app.R                              # Shiny app (deployed to shinyapps.io)
├── 04142026_Rasterization.R           # fish point data → data/raster_list.rds
├── 04292026_InvertRasterization.R     # invert point data → data/invert_raster_list.rds
├── 07132026_COG.R                     # center-of-gravity analysis (exploratory)
├── MPA_Analysis.R                     # MPA coverage summaries (exploratory)
├── 03152026_PrePostComparison.R       # pre vs post heat wave predictor comparison
├── 04022026_DataExploration.R         # scratch data check
├── plot1_densities.png … plot5_delta_summary.png   # outputs of PrePostComparison
├── data/
│   ├── CPFV_PA_sp_decade_roms_n_ensemble_shiny.rds     # fish suitability (input)
│   ├── Inverts_PA_sp_decade_roms_n_ensemble_shiny.rds  # invert suitability (input)
│   ├── all_sp_proportion_in_MPA_area_CA.rds            # proportion of suitable area in MPAs, by species/year
│   ├── raster_list.rds                                 # generated: fish rasters
│   ├── invert_raster_list.rds                          # generated: invert rasters
│   └── test_preds_shiny.csv                            # predictors pre/post heat wave
├── shp/
│   └── California_Marine_Protected_Areas_[ds582].*     # MPA boundaries (keep all sidecar files together)
└── rsconnect/                         # shinyapps.io deployment config (casellelab account)
```

---

## Setup

Install the packages from R:

```r
install.packages(c(
  "shiny", "leaflet", "dplyr", "tidyr", "RColorBrewer", "ggplot2",
  "plotly", "jsonlite", "sf",                     # app.R
  "terra", "raster", "readr", "tidyverse",        # rasterization scripts
  "patchwork"                                      # PrePostComparison
))
```

---

## How to run

**Run the app locally.** Set the working directory to the repo root, then:

```r
shiny::runApp()
```

**Rebuild the raster lists.** You only need to do this when the input suitability `.rds` files change:

```r
source("04142026_Rasterization.R")       # fish   → data/raster_list.rds
source("04292026_InvertRasterization.R") # inverts → data/invert_raster_list.rds
```

**Redeploy to shinyapps.io.** You need access to the `casellelab` account:

```r
rsconnect::deployApp(appName = "prop80")
```

---

## Data format notes

**Suitability files** (`*_PA_sp_decade_roms_n_ensemble_shiny.rds`) need these columns:

| Column | Meaning |
|---|---|
| `cell_coord_id` | Grid cell ID in the form `"<lat>-<lon>"`, with longitude stored as positive degrees west. Scripts parse lat/lon from it and negate the longitude. |
| `species` | Species name, matching the names in the Species table above |
| `decade` | Decade label, e.g. `"2041_2050"` |
| `roms_model` | `gfdl`, `hadl`, `ipsl`, `ensemble`, or `histnew` |
| `pa_decade_mean_pa` | Mean decadal environmental suitability (0–1) |

**Raster lists** are named `<species>_<decade>`, e.g. `California_halibut_2041_2050`. The app looks rasters up by this exact key.

**Decades** are `1995_2000`, then 10-year bins from `2001_2010` through `2091_2100`. The bins are defined in `assign_decade()` in `app.R`.

**MPA coverage** (`all_sp_proportion_in_MPA_area_CA.rds`) needs the columns `species`, `year_factor`, and `prop_mpa`. The app bins these by decade and plots mean ± standard error.

---

## App features

- Toggle between the Fish and Invert species groups, then select a species
- Decade slider with play/pause animation
- Suitability map with a fixed 0–1 color scale, plus MPA boundary overlay
- Time series of decadal suitability for each ROMS model and the ensemble
- Proportion of suitable area within MPAs over time

---

## Known issues

- **`app.R` has leftover analysis code after `shinyApp(ui, server)`** (lines ~479–518, a center-of-gravity map and correlation test). It references `cog_fish` / `cog_inverts`, which `app.R` never defines. Shiny expects the last expression in `app.R` to be the app object, so this code can break local runs and redeploys. Move it into `07132026_COG.R` and make `shinyApp(ui, server)` the last line.
- **`MPA_Analysis.R` is inconsistent with the app.** It uses `"2041-2050"` (hyphen) where everything else uses `"2041_2050"`, and it overwrites `decade` with `floor(year / 10) * 10`. It also calls `case_when()` without loading `dplyr`. The app has its own working version, so treat this script as scratch.
- **`.Rhistory` and `.DS_Store` are committed.** Add them to `.gitignore`.

*[Delete this section once fixed.]*

---

## Exploratory scripts

These scripts aren't required to run the app:

- **`07132026_COG.R`**: computes a suitability-weighted center of gravity (lat/lon) for each species and decade, to track range shifts.
- **`03152026_PrePostComparison.R`**: compares model predictors (groundfish biomass, oxygen, salinity, temperature) before and after the marine heat wave, using Wilcoxon tests, and writes `plot1`–`plot5`.
- **`04022026_DataExploration.R`**: quick look at the fish data.
