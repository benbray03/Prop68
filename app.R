library(shiny)
library(leaflet)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggplot2)
library(plotly)
library(jsonlite)

# ── Load data ──────────────────────────────────────────────────────────────────
raster_list        <- readRDS("data/raster_list.rds")
raster_list_inverts <- readRDS("data/invert_raster_list.rds")

dat_fish   <- readRDS("data/CPFV_PA_sp_decade_roms_n_ensemble_shiny.rds")
dat_inverts <- readRDS("data/Inverts_PA_sp_decade_roms_n_ensemble_shiny.rds")

# ── Helper: prepare a dat object for either group ─────────────────────────────
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

# ── Raw data for time series ───────────────────────────────────────────────────
raw_fish    <- readRDS("data/CPFV_PA_sp_decade_roms_n_ensemble_shiny.rds") %>%
  mutate(species = as.character(species), decade = as.character(decade))
raw_inverts <- readRDS("data/Inverts_PA_sp_decade_roms_n_ensemble_shiny.rds") %>%
  mutate(species = as.character(species), decade = as.character(decade))

make_raw_by_species <- function(d) {
  d %>% group_by(species) %>% group_split() %>%
    setNames(sort(unique(d$species)))
}
make_ens_by_species <- function(d) {
  d %>% group_by(species) %>% group_split() %>%
    setNames(sort(unique(d$species)))
}

raw_fish_by_sp    <- make_raw_by_species(raw_fish)
raw_inverts_by_sp <- make_raw_by_species(raw_inverts)
ens_fish_by_sp    <- make_ens_by_species(dat_fish)
ens_inverts_by_sp <- make_ens_by_species(dat_inverts)

# ── Species name lookups ───────────────────────────────────────────────────────
fish_labels <- c(
  "Black_rockfish"       = "Black Rockfish",
  "California_barracuda" = "California Barracuda",
  "California_halibut"   = "California Halibut",
  "California_sheephead" = "California Sheephead",
  "Copper_rockfish"      = "Copper Rockfish",
  "Gopher_rockfish"      = "Gopher Rockfish",
  "Ocean_whitefish"      = "Ocean Whitefish"
)

invert_labels <- c(
  "red_urchin" = "Red Urchin",
  'CA_Lobster' = "CA Spiny Lobster",
  'market_squid' = "Market Squid",
  'pink_shrimp' = "Pink Shrimp",
  'Ridgeback_prawn' = "Ridgeback Prawn",
  'Sea_cucumber' = "Sea Cucumber"
)

# ── Pre-split helpers ──────────────────────────────────────────────────────────
make_split <- function(d) {
  keys <- sort(unique(paste(d$species, d$decade, sep = "||")))
  d %>%
    arrange(species, decade) %>%
    group_by(species, decade) %>%
    group_split() %>%
    setNames(keys)
}

fish_split    <- make_split(dat_fish)
inverts_split <- make_split(dat_inverts)

fish_species_list    <- sort(unique(dat_fish$species))
inverts_species_list <- sort(unique(dat_inverts$species))

decade_list <- sort(unique(dat_fish$decade))

# ── Colour palette ─────────────────────────────────────────────────────────────
pal_fun <- colorNumeric(
  palette  = rev(brewer.pal(11, "RdYlBu")),
  domain   = c(0, 1),
  na.color = "transparent"
)

# ── UI ─────────────────────────────────────────────────────────────────────────
ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      body { font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif;
             background: #f4f6f9; margin: 0; }
      .app-header { background: #1a3a5c; color: white; padding: 14px 24px; }
      .app-header h2 { margin: 0; font-size: 1.25rem; font-weight: 600; }
      .app-header .subtitle { font-size: 0.85rem; opacity: 0.75; }
      .sidebar-panel { background: white; border-radius: 8px;
                       box-shadow: 0 1px 4px rgba(0,0,0,.12);
                       padding: 18px; margin: 16px 8px 16px 16px; }
      .map-panel { margin: 16px 16px 16px 8px; border-radius: 8px;
                   overflow: hidden; box-shadow: 0 1px 4px rgba(0,0,0,.15); }
      .species-btn { display: block; width: 100%; text-align: left;
                     padding: 8px 12px; margin-bottom: 6px; border: none;
                     border-radius: 6px; cursor: pointer;
                     background: #eef2f7; color: #1a3a5c;
                     font-size: 0.88rem; font-weight: 500;
                     transition: background .15s, color .15s; }
      .species-btn:hover  { background: #c8d8ec; }
      .species-btn.active { background: #1a3a5c; color: white; }
      .section-label { font-size: 0.75rem; font-weight: 700; letter-spacing: .06em;
                       text-transform: uppercase; color: #7a90a8;
                       margin: 16px 0 8px; }
      .stat-box { background: #eef2f7; border-radius: 6px; padding: 10px 14px;
                  margin-top: 14px; font-size: 0.82rem; color: #1a3a5c; }
      .stat-box b { font-size: 1rem; }
      .irs-bar { background: #4a8abf !important; border-color: #4a8abf !important; }
      .irs-line { background: #c0cfe0 !important; }
      .irs-single { background: #1a3a5c !important; }
      #decade_select { width: 100%; border-radius: 6px;
                       border: 2px solid #c0cfe0; padding: 7px 10px;
                       color: #1a3a5c; font-size: 0.88rem; font-weight: 600;
                       background: white; }
      /* Group tab switcher */
      .group-tabs { display: flex; gap: 0; margin-bottom: 14px;
                    border-radius: 8px; overflow: hidden;
                    border: 2px solid #1a3a5c; }
      .group-tab  { flex: 1; padding: 7px 0; text-align: center;
                    font-size: 0.85rem; font-weight: 700; cursor: pointer;
                    background: white; color: #1a3a5c;
                    border: none; transition: background .15s, color .15s; }
      .group-tab:hover  { background: #c8d8ec; }
      .group-tab.active { background: #1a3a5c; color: white; }
    "))
  ),
  
  div(class = "app-header",
      h2("Species Distribution — Mean Decadal Environmental Suitability"),
      div(class = "subtitle", "ROMS ensemble model · ensemble Mean Decadal Environmental Suitability"),
  ),
  
  fluidRow(
    column(3,
           div(class = "sidebar-panel",
               
               # ── Group tab switcher ──
               div(class = "section-label", "Decade"),
               sliderInput("decade_select", label = NULL,
                           min   = 1,
                           max   = length(decade_list),
                           value = 1,
                           step  = 1,
                           width = "100%",
                           ticks = FALSE,
                           animate = FALSE),
               div(class = "group-tabs",
                   tags$button(id = "tab_fish",   class = "group-tab active", "Fish"),
                   tags$button(id = "tab_inverts", class = "group-tab",       "Inverts")
               ),
               div(class = "section-label", "Species"),
               uiOutput("species_buttons"),
               div(class = "stat-box", uiOutput("stats_out")),
               div(class = "section-label", "Mean Decadal Environmental Suitability"),
               plotlyOutput("ts_plot", height = "220px"),
               uiOutput("ts_legend")
           )
    ),
    column(9,
           div(class = "map-panel",
               leafletOutput("map", height = "82vh")
           )
    )
  ),
  
  tags$script(HTML("
    // Species buttons
    $(document).on('click', '.species-btn', function() {
      $('.species-btn').removeClass('active');
      $(this).addClass('active');
      Shiny.setInputValue('selected_species', $(this).data('val'), {priority: 'event'});
    });

    // Group tabs
    $(document).on('click', '.group-tab', function() {
      $('.group-tab').removeClass('active');
      $(this).addClass('active');
      var grp = $(this).attr('id') === 'tab_fish' ? 'fish' : 'inverts';
      Shiny.setInputValue('selected_group', grp, {priority: 'event'});
    });
  ")),
  
  tags$script(HTML(paste0("
    function updateDecadeSlider() {
      var idx = $('#decade_select').val() - 1;
      var decades = ", jsonlite::toJSON(decade_list), ";
      $('.irs-single').text(decades[idx]);
      $('.irs-min').text(decades[0]);
      $('.irs-max').text(decades[decades.length - 1]);
    }
    $(document).on('input change', '#decade_select', updateDecadeSlider);
    $(document).ready(function() { setTimeout(updateDecadeSlider, 300); });
  ")))
)

# ── Server ─────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {
  
  rv <- reactiveValues(
    group   = "fish",
    species = fish_species_list[1],
    decade  = decade_list[1]
  )
  
  observeEvent(input$selected_group, {
    rv$group <- input$selected_group
    # Switch to first species of the new group
    sp_list <- if (input$selected_group == "fish") fish_species_list else inverts_species_list
    rv$species <- sp_list[1]
  })
  
  observeEvent(input$selected_species, { rv$species <- input$selected_species })
  observeEvent(input$decade_select,    { rv$decade  <- decade_list[input$decade_select] })
  
  # Reactive helpers that depend on current group
  current_species_list <- reactive({
    if (rv$group == "fish") fish_species_list else inverts_species_list
  })
  current_labels <- reactive({
    if (rv$group == "fish") fish_labels else invert_labels
  })
  current_split <- reactive({
    if (rv$group == "fish") fish_split else inverts_split
  })
  current_raw_by_sp <- reactive({
    if (rv$group == "fish") raw_fish_by_sp else raw_inverts_by_sp
  })
  current_ens_by_sp <- reactive({
    if (rv$group == "fish") ens_fish_by_sp else ens_inverts_by_sp
  })
  current_raster_list <- reactive({
    if (rv$group == "fish") raster_list else raster_list_inverts
  })
  current_dat <- reactive({
    if (rv$group == "fish") dat_fish else dat_inverts
  })
  
  # Species buttons
  output$species_buttons <- renderUI({
    lapply(current_species_list(), function(sp) {
      cls   <- if (sp == rv$species) "species-btn active" else "species-btn"
      label <- current_labels()[[sp]]
      if (is.null(label) || is.na(label)) label <- gsub("_", " ", sp)
      tags$button(class = cls, `data-val` = sp, label)
    })
  })
  
  # Filtered subset
  subset_dat <- reactive({
    key <- paste(rv$species, rv$decade, sep = "||")
    current_split()[[key]]
  })
  
  # Stats box
  output$stats_out <- renderUI({
    d <- subset_dat()
    if (is.null(d) || nrow(d) == 0) return(p("No data for this selection."))
    label <- current_labels()[[rv$species]]
    if (is.null(label) || is.na(label)) label <- gsub("_", " ", rv$species)
    HTML(sprintf(
      "<b>%s</b><br>Decade: %s<br><br>
       Mean Decadal Environmental Suitability: <b>%.3f</b><br>",
      label, rv$decade,
      mean(d$pa_decade_mean_pa, na.rm = TRUE),
      min(d$pa_decade_mean_pa,  na.rm = TRUE),
      max(d$pa_decade_mean_pa,  na.rm = TRUE)
    ))
  })
  
  # Base map
  output$map <- renderLeaflet({
    leaflet() %>%
      addProviderTiles(providers$CartoDB.Positron) %>%
      setView(lng = -122, lat = 36, zoom = 6) %>%
      addLegend(
        position  = "bottomright",
        pal       = pal_fun,
        values    = seq(0, 1, by = 0.1),
        title     = "Mean Decadal<br>Environmental<br>Suitability",
        opacity   = 0.85,
        labFormat = labelFormat(transform = function(x) sort(x, decreasing = TRUE))
      )
  })
  
  observe({
    leafletProxy("map") %>% clearGroup("pa_layer")
    key <- paste0(rv$species, "_", rv$decade)
    r   <- current_raster_list()[[key]]
    if (is.null(r)) return()
    leafletProxy("map") %>%
      addRasterImage(r, colors = pal_fun, opacity = 0.85, group = "pa_layer")
  })
  
  # Time series plot
  output$ts_plot <- renderPlotly({
    decade_levels <- sort(unique(current_dat()$decade))
    
    ens <- current_ens_by_sp()[[rv$species]] %>%
      group_by(decade) %>%
      summarise(mean_pa = mean(pa_decade_mean_pa, na.rm = TRUE), .groups = "drop") %>%
      mutate(decade = factor(decade, levels = decade_levels))
    
    mod <- current_raw_by_sp()[[rv$species]] %>%
      group_by(decade, roms_model) %>%
      summarise(mean_pa = mean(pa_decade_mean_pa, na.rm = TRUE), .groups = "drop") %>%
      mutate(decade = factor(decade, levels = decade_levels))
    
    p <- ggplot() +
      geom_point(data = mod, aes(x = decade, y = mean_pa, color = roms_model),
                 size = 2.5, position = position_dodge(width = 0.4), show.legend = FALSE) +
      scale_color_manual(values = c("gfdl"     = "#2a9d8f",
                                    "hadl"     = "#e05c2a",
                                    "ipsl"     = "#e9c46a",
                                    "ensemble" = "#1a3a5c",
                                    "histnew"  = "#7a90a8"),
                         name = "Model") +
      scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
      labs(x = NULL, y = "Mean Decadal Environmental<br>Suitability") +
      theme_minimal(base_size = 11) +
      theme(
        axis.text.x        = element_text(angle = 45, hjust = 1, size = 8),
        axis.ticks.x       = element_blank(),
        axis.text.y        = element_text(size = 8),
        legend.position    = "none",
        panel.grid.major.x = element_blank(),
        plot.background    = element_rect(fill = "transparent", color = NA),
        panel.background   = element_rect(fill = "transparent", color = NA)
      )
    
    ggplotly(p, tooltip = c("y", "colour")) %>%
      layout(showlegend = FALSE) %>%
      config(displayModeBar = FALSE)
  })
  
  output$ts_legend <- renderUI({
    model_colors <- c(
      "gfdl"    = "#2a9d8f",
      "ipsl"    = "#e9c46a",
      "hadl"    = "#e05c2a",
      "ensemble"= "#1a3a5c",
      "histnew" = "#7a90a8"
    )
    items <- lapply(names(model_colors), function(m) {
      tags$div(style = "display: inline-flex; align-items: center; margin-right: 10px;",
               tags$span(style = paste0(
                 "display:inline-block; width:10px; height:10px; border-radius:50%;",
                 "background:", model_colors[[m]], "; margin-right:4px;"
               )),
               tags$span(style = "font-size:0.75rem; color:#1a3a5c;", m)
      )
    })
    tags$div(
      style = "display: flex; flex-wrap: wrap; gap: 4px; padding: 6px 0;",
      tags$span(style = "font-size:0.75rem; font-weight:700; color:#7a90a8;
                         margin-right:6px; align-self:center;", "MODEL"),
      items
    )
  })
}

shinyApp(ui, server)
