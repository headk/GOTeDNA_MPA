ui <- fluidPage(
  shinyjs::useShinyjs(),

  tags$head(
    tags$style(HTML("
      body { padding-top: 62px; }
      .navbar { margin-bottom: 10px; }
      #map { height: calc(95vh - 62px) !important; }

      #map_wrap { position: relative; }
      #floating_panel {
        background: rgba(255,255,255,0.92);
        border-radius: 6px;
        box-shadow: 0 2px 10px rgba(0,0,0,0.25);
        padding: 10px 12px;
        max-height: calc(95vh - 30px);
        overflow-y: auto;
        z-index: 999;
      }

      .filter-btn-grid{
        display: grid;
        grid-template-columns: repeat(3, 1fr);
        gap: 8px;
      }
      .filter-btn{
        width: 100% !important;
        height: 70px !important;
        padding: 6px 10px !important;
        font-size: 15px !important;
        white-space: nowrap;
        overflow: hidden;
        text-overflow: ellipsis;
      }

      #floating_toggle {
        width: 100%;
        text-align: left;
        display: flex;
        align-items: center;
        justify-content: space-between;
        font-size: 16px;
        font-weight: 600;
        padding: 10px 12px;
      }
      #floating_toggle .caret-icon { transition: transform 150ms ease; }
      #floating_toggle[aria-expanded='true'] .caret-icon { transform: rotate(180deg); }

      .scroll-section { scroll-margin-top: 80px; }
          /* Give the legend its “white box” back */
.legend-richness {
  background: rgba(255,255,255,0.92) !important;
  border-radius: 6px !important;
  box-shadow: 0 2px 10px rgba(0,0,0,0.25) !important;
  border: 1px solid rgba(0,0,0,0.15) !important;
}
.legend-richness {
  padding: 14px 16px !important;
  width: 160px !important;

  background: rgba(255,255,255,0.92) !important;
  border-radius: 6px !important;
  box-shadow: 0 2px 10px rgba(0,0,0,0.25) !important;
  border: 1px solid rgba(0,0,0,0.15) !important;
}
.legend-richness i {
  display: inline-block !important;
}
.legend-richness {
  padding: 14px 16px !important;
  width: 160px !important;

  font-size: 16px !important;
  line-height: 1.4 !important;
  font-weight: 400 !important;

  background: rgba(255,255,255,0.92) !important;
  border-radius: 6px !important;
  box-shadow: 0 2px 10px rgba(0,0,0,0.25) !important;
  border: 1px solid rgba(0,0,0,0.15) !important;
}
/* === protocol cards (GOTeDNA palette) === */
.protocol-grid {
  display: grid;
  grid-template-columns: repeat(auto-fill, minmax(140px, 1fr));
  gap: 12px;
}

.protocol-field-title {
  font-size: 13px;
  font-weight: 600;
  margin-bottom: 4px;
  color: #0b2a2a;
}

.protocol-card {
  background: #EAF7FB;            /* pastel baby blue */
  border: 1px solid rgba(49, 180, 219, 0.75); /* deepskyblue border */
  border-radius: 8px;
  padding: 10px 12px;
  min-height: 38px;
  font-size: 14px;
  font-weight: 600;
  color: #0b2a2a;                 /* readable dark text */
  box-shadow: 0 1px 6px rgba(44, 98, 203, 0.10);
  transition: background-color 250ms ease, border-color 250ms ease, box-shadow 250ms ease, color 250ms ease;
}

/* changed card = darker */
.protocol-card.changed {
  background: #9ADBE8;            /* royalblue (DARKEST = changed) */
  border-color: #2f9ae6;          /* dodgerblue border */
  color: #0b2a2a;
  box-shadow: 0 0 0 3px rgba(32, 224, 98, 0.35); /* spring green glow */
}

/* optional: for missing values */
.protocol-card.na {
  background: #EAF7FB;            /* pastel baby blue */
  border-color: #31b4db;          /* deepskyblue */
  color: #0b2a2a;
  opacity: 0.9;
}

/* === filter button grid + equal-size buttons === */
.filter-btn-grid{
  display: grid;
  grid-template-columns: repeat(3, 1fr);  /* 3 buttons per row (9 = 3 rows) */
  gap: 8px;
}

.filter-btn{
  width: 100% !important;   /* all same width within grid cell */
  height: 70px !important;  /* all same height */
  padding: 6px 10px !important;
  font-size: 15px !important;
  white-space: nowrap;      /* prevent wrapping */
  overflow: hidden;
  text-overflow: ellipsis;  /* ... if label too long */
}

#map_wrap { position: relative; }

#floating_panel {
background: rgba(255,255,255,0.92);
border-radius: 6px;
box-shadow: 0 2px 10px rgba(0,0,0,0.25);
padding: 10px 12px;
max-height: calc(95vh - 30px);
overflow-y: auto;
z-index: 999;
}

/* Header button */
  #floating_toggle {
  width: 100%;
text-align: left;
display: flex;
align-items: center;
justify-content: space-between;
font-size: 16px;
font-weight: 600;
padding: 10px 12px;
}
#floating_toggle .caret-icon {
transition: transform 150ms ease;
}
#floating_toggle[aria-expanded='true'] .caret-icon {
transform: rotate(180deg);
}

/* Optional: tighter body padding */
  #floating_panel .panel-body { padding-top: 10px; }

  /* Example: keep your existing button grid styles */
  .filter-btn-grid{
    display: grid;
    grid-template-columns: repeat(3, 1fr);
    gap: 8px;
  }
.filter-btn{
  width: 100% !important;
  height: 70px !important;
  padding: 6px 10px !important;
  font-size: 15px !important;
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
}
    ")),
    tags$script(HTML("
      $(document).on('click', 'a.nav-scroll', function(e){
        e.preventDefault();
        var target = $(this).data('target');
        var el = document.getElementById(target);
        if(el){
          el.scrollIntoView({behavior:'smooth', block:'start'});
        }
      });

      Shiny.addCustomMessageHandler('openFloating', function(message) {
        var $body = $('#floating_body');
        if ($body.length && !$body.hasClass('in')) {
          $body.collapse('show');
          $('#floating_toggle').attr('aria-expanded', 'true');
        }
      });
    "))
  ),

  # ---- REAL NAVBAR ----
  tags$nav(
    class = "navbar navbar-default navbar-fixed-top",
    tags$div(
      class = "container-fluid",
      tags$div(
        class = "navbar-header",
        tags$a(class = "navbar-brand nav-scroll", href = "#", `data-target` = "sec_map", "GOTeDNA-MPA")
      ),
      tags$ul(
        class = "nav navbar-nav",
        tags$li(tags$a(class="nav-scroll", href="#", `data-target`="sec_map",    "Map")),
        tags$li(tags$a(class="nav-scroll", href="#", `data-target`="sec_sara",   "SAR/AIS Detection")),
        tags$li(tags$a(class="nav-scroll", href="#", `data-target`="sec_method", "Method Comparison")),
        tags$li(tags$a(class="nav-scroll", href="#", `data-target`="sec_div",    "Diversity Metrics")),
        tags$li(tags$a(class="nav-scroll", href="#", `data-target`="sec_pie",    "Taxonomic Pie Chart"))
      )
    )
  ),

  # ---- MAP SECTION ----
  div(
    id = "sec_map", class = "scroll-section",
    div(
      id = "map_wrap",
      leafletOutput("map"),

      absolutePanel(
        id = "floating_panel",
        fixed = FALSE, draggable = TRUE,
        top = 10, left = 70, width = 360,

        tags$button(
          id = "floating_toggle",
          type = "button",
          class = "btn btn-default",
          `data-toggle` = "collapse",
          `data-target` = "#floating_body",
          `aria-expanded` = "false",
          `aria-controls` = "floating_body",
          tagList(
            tags$span("Select a site"),
            tags$span(class = "caret-icon", HTML("&#9662;"))
          )
        ),

        div(
          id = "floating_body",
          class = "collapse",
          div(
            class = "panel-body",
            h4("Filter"),
            selectInput("sel_year", "Year",
                        choices = c("All", "2022", "2023", "2024"),
                        selected = "All"),
            h4("Group"),
            div(
              class = "filter-btn-grid",
              actionButton("total_fish",       "Fishes",     class = "btn btn-default filter-btn"),
              actionButton("total_mammals",    "Mammals",    class = "btn btn-default filter-btn"),
              actionButton("total_reptiles",   "Reptiles",   class = "btn btn-default filter-btn"),
              actionButton("total_birds",      "Birds",      class = "btn btn-default filter-btn"),
              actionButton("total_molluscs",   "Molluscs",   class = "btn btn-default filter-btn"),
              actionButton("total_arthropods", "Arthropods", class = "btn btn-default filter-btn"),
              actionButton("total_plants",     "Plants",     class = "btn btn-default filter-btn"),
              actionButton("SARA",             "SARA",       class = "btn btn-default filter-btn"),
              actionButton("AIS",              "AIS",        class = "btn btn-default filter-btn")
            ),
            hr(),
            h4("Species list"),
            uiOutput("species_panel")
          )
        )
      )
    )
  ),

  # ---- SARA/AIS SECTION ----
  div(
    id = "sec_sara", class = "scroll-section",
    tabsetPanel(
      tabPanel("SARA Details: Schedule 1", DT::DTOutput("sara_details")),
      tabPanel("AIS Details",             DT::DTOutput("ais_details")),
      tabPanel("Detections",              DT::DTOutput("detections_tbl"))
    )
  ),

  # ---- METHOD COMPARISON SECTION ----
  div(
    id = "sec_method", class = "scroll-section",
    h3("Method Comparison"),

    div(
      id = "data_request_wrap",
      style = "padding: 10px 12px; background: rgba(255,255,255,0.92);
               border-radius: 6px; box-shadow: 0 2px 10px rgba(0,0,0,0.15);
               margin-bottom: 10px;",
      h4("Data Request"),
      fluidRow(
        column(
          width = 4,
          selectInput("req_location", "Location",
                      choices = character(0), selected = NULL, selectize = FALSE),
          selectInput("req_protocol", "ProtocolID",
                      choices = character(0), selected = NULL, selectize = FALSE)
        ),
        column(
          width = 8,
          uiOutput("protocol_details")
        )
      )
    )
  ),

  div(id="sec_div", class="scroll-section", h3("Diversity Metrics"), tags$div("...")),
  div(id="sec_pie", class="scroll-section", h3("Taxonomic Pie Chart"), tags$div("..."))
)
