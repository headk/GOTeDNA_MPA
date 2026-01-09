ui <- fluidPage(
  shinyjs::useShinyjs(),

  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
  ),

  # ---- Data Request panel ABOVE the map ----
  div(
    id = "data_request_wrap",
    style = "padding: 10px 12px; background: rgba(255,255,255,0.92);
             border-radius: 6px; box-shadow: 0 2px 10px rgba(0,0,0,0.15);
             margin-bottom: 10px;",
    h4("Data Request"),
    fluidRow(
      column(
        width = 4,
        selectInput("req_location", "Location", choices = character(0), selected = NULL, selectize = FALSE),
        selectInput("req_protocol", "ProtocolID", choices = character(0), selected = NULL, selectize = FALSE)
      ),
      column(width = 8, uiOutput("protocol_details"))
    )
  ),

  div(
    id = "map_wrap",
    leafletOutput("map", height = "95vh"),
    absolutePanel(
      id="floating_panel",
      fixed = FALSE,
      draggable = TRUE,
      top = 10, left = 70, width = 360,

      h4("Selected area"),
      uiOutput("cell_summary"),
      hr(),

      h4("Filter"),
      selectInput("sel_year", "Year", choices = c("All","2022","2023","2024"), selected = "All"),

      h4("Group"),
      div(
        class = "filter-btn-grid",
        actionButton(
          "total_fish",
          label = tagList(
            tags$img(src = "img/fish.png", height = "28px"),
            tags$span("Fishes")
          ),
          class = "btn btn-default filter-btn"
        ),
        actionButton(
          "total_mammals",
          label = tagList(
            tags$img(src = "img/whale2.png", height = "28px"),
            tags$span("Mammals")
          ),
          class = "btn btn-default filter-btn"
        ),
        actionButton(
          "total_arthropods", 
          label = tagList(
            tags$img(src = "img/lobster.png", height = "28px"),
            tags$span("Arthropods")
          ),
          class = "btn btn-default filter-btn"
        ),
        actionButton(
          "total_reptiles",
          label = tagList(
            tags$img(src = "img/turtle.png", height = "28px"),
            tags$span("Reptiles")
          ),
          class = "btn btn-default filter-btn"
        ),
        actionButton(
          "total_birds",
          label = tagList(
            tags$img(src = "img/bird.png", height = "28px"),
            tags$span("Birds")
          ),
          class = "btn btn-default filter-btn"
        ), 
        actionButton(
          "total_molluscs",
          label = tagList(
            tags$img(src = "img/mollusc.png", height = "28px"),
            tags$span("Molluscs")
          ),
          class = "btn btn-default filter-btn"
        ), 
        actionButton(
          "total_plants",
          label = tagList(
            tags$img(src = "img/plant.png", height = "28px"),
            tags$span("Plants")
          ),
          class = "btn btn-default filter-btn"
        ),         
        actionButton("SARA", "SARA", class = "btn btn-default filter-btn"),
        actionButton("AIS", "AIS", class = "btn btn-default filter-btn")
      ),

      hr(),
      h4("Species list"),
      uiOutput("species_panel")
    )
  ),

  tabsetPanel(
    id = "tabs_below",
    tabPanel("SARA Details: Schedule 1", DT::DTOutput("sara_details")),
    tabPanel("AIS Details", DT::DTOutput("ais_details")),
    tabPanel("Detections", DT::DTOutput("detections_tbl"))
  )
)

