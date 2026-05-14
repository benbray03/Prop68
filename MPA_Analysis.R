library(readr)

MPA <- readRDS("data/all_sp_proportion_in_MPA_area_CA.rds")

assign_decade <- function(year) {
  case_when(
    year >= 1995 & year <= 2000 ~ "1995_2000",
    year >= 2001 & year <= 2010 ~ "2001_2010",
    year >= 2011 & year <= 2020 ~ "2011_2020",
    year >= 2021 & year <= 2030 ~ "2021_2030",
    year >= 2031 & year <= 2040 ~ "2031_2040",
    year >= 2041 & year <= 2050 ~ "2041-2050",
    year >= 2051 & year <= 2060 ~ "2051_2060",
    year >= 2061 & year <= 2070 ~ "2061_2070",
    year >= 2071 & year <= 2080 ~ "2071_2080",
    year >= 2081 & year <= 2090 ~ "2081_2090",
    year >= 2091 & year <= 2100 ~ "2091_2100",
    TRUE ~ NA_character_
  )
}

decade_order <- c("1995_2000","2001_2010","2011_2020","2021_2030","2031_2040",
                  "2041-2050","2051_2060","2061_2070","2071_2080","2081_2090",
                  "2091_2100")


MPA <- MPA %>%
  mutate(
    species = as.character(species),
    year    = as.numeric(as.character(year_factor)),
    decade  = as.character(floor(year / 10) * 10)
  ) %>%
  filter(!is.na(prop_mpa))

mpa_decadal <- MPA %>%
  group_by(species, decade) %>%
  summarise(
    mean_prop = mean(prop_mpa, na.rm = TRUE),
    se_prop   = sd(prop_mpa, na.rm = TRUE) / sqrt(n()),
    .groups   = "drop"
  )

mpa_by_sp <- mpa_decadal %>%
  group_by(species) %>%
  group_split() %>%
  setNames(sort(unique(mpa_decadal$species)))
