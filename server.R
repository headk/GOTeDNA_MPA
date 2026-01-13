server <- function(input, output, session){
  
  # ---- Data Request: use occ_all directly ----
  meta_all <- reactive({
    req(occ_all)
    
    occ_all %>%
      dplyr::mutate(
        Location   = as.character(Location),
        ProtocolID = as.character(ProtocolID)
      )
  })
  
  # Populate Location dropdown once meta_all is available
  observeEvent(meta_all(), {
    df <- meta_all()
    
    locs <- df %>%
      dplyr::filter(!is.na(Location), Location != "") %>%
      dplyr::distinct(Location) %>%
      dplyr::arrange(Location) %>%
      dplyr::pull(Location)
    
    updateSelectInput(
      session, "req_location",
      choices  = locs,
      selected = if (length(locs)) locs[[1]] else NULL
    )
  }, ignoreInit = FALSE)
  
  # Populate ProtocolID dropdown (filtered by selected Location)
  observeEvent(list(meta_all(), input$req_location), {
    df <- meta_all()
    
    if (!is.null(input$req_location) && nzchar(input$req_location)) {
      df <- df %>% dplyr::filter(Location == input$req_location)
    }
    
    prots <- df %>%
      dplyr::filter(!is.na(ProtocolID), ProtocolID != "") %>%
      dplyr::distinct(ProtocolID) %>%
      dplyr::arrange(ProtocolID) %>%
      dplyr::pull(ProtocolID)
    
    updateSelectInput(
      session, "req_protocol",
      choices  = prots,
      selected = if (length(prots)) prots[[1]] else NULL
    )
  }, ignoreInit = FALSE)
  
  # helper: pick a single display value (unique or collapse)
  pick_display <- function(x) {
    x <- as.character(x)
    x <- x[!is.na(x) & nzchar(trimws(x))]
    if (length(x) == 0) return(NA_character_)
    ux <- unique(x)
    if (length(ux) == 1) ux else paste(ux, collapse = " | ")
  }
  
  # rows matching current Location + ProtocolID
  selected_protocol_rows <- reactive({
    req(input$req_location, input$req_protocol)
    
    meta_all() %>%
      dplyr::filter(
        Location   == input$req_location,
        ProtocolID == input$req_protocol
      )
  })
  
  # --- base protocol per Location (the "starting choice") ---
  base_protocol <- reactiveVal(NULL)
  
  # whenever Location changes, set base protocol to the first available protocol for that location
  observeEvent(list(meta_all(), input$req_location), {
    df <- meta_all()
    req(input$req_location)
    
    prots <- df %>%
      dplyr::filter(Location == input$req_location) %>%
      dplyr::filter(!is.na(ProtocolID), ProtocolID != "") %>%
      dplyr::distinct(ProtocolID) %>%
      dplyr::arrange(ProtocolID) %>%
      dplyr::pull(ProtocolID)
    
    base_protocol(if (length(prots)) as.character(prots[[1]]) else NULL)
  }, ignoreInit = FALSE)
  
  # live details card
  output$protocol_details <- renderUI({
    df <- selected_protocol_rows()
    
    if (nrow(df) == 0) {
      return(tags$div(style="margin-top:10px;", em("No rows found for this Location + ProtocolID.")))
    }
    
    field_map <- c(
      "Year"           = "Year",
      "Volume"         = "Volume",
      "Depth"          = "Depth",
      "Filter"         = "Filter",
      "Preservative"   = "Preservative",
      "Primers"        = "Primers",
      "Sequencer"      = "Sequencer",
      "Bioinformatics" = "Bioinformatics"
    )
    
    # --- current values for the selected protocol ---
    vals <- lapply(names(field_map), function(lbl) {
      col <- field_map[[lbl]]
      if (!col %in% names(df)) return(NA_character_)
      pick_display(df[[col]])
    })
    names(vals) <- names(field_map)
    
    # --- base values (compare against base_protocol for this Location) ---
    base_vals <- setNames(as.list(rep(NA_character_, length(field_map))), names(field_map))
    
    bp <- base_protocol()
    if (!is.null(bp) && nzchar(bp)) {
      df_base <- meta_all() %>%
        dplyr::filter(
          Location   == input$req_location,
          ProtocolID == bp
        )
      
      if (nrow(df_base) > 0) {
        base_vals <- lapply(names(field_map), function(lbl) {
          col <- field_map[[lbl]]
          if (!col %in% names(df_base)) return(NA_character_)
          pick_display(df_base[[col]])
        })
        names(base_vals) <- names(field_map)
      }
    }
    
    # --- compute which cards changed vs base ---
    changed_flags <- lapply(names(field_map), function(lbl) {
      cur <- vals[[lbl]]
      bas <- base_vals[[lbl]]
      
      cur <- ifelse(is.na(cur), "", as.character(cur))
      bas <- ifelse(is.na(bas), "", as.character(bas))
      
      nzchar(trimws(cur)) && nzchar(trimws(bas)) && !identical(trimws(cur), trimws(bas))
    })
    names(changed_flags) <- names(field_map)
    
    # --- UI ---
    tags$div(
      style = "margin-top:12px;",
      
      tags$div(
        style="font-weight:600; font-size:16px; margin-bottom:10px;",
        "Methods"
      ),
      
      tags$div(
        class = "protocol-grid",
        
        lapply(names(vals), function(lbl) {
          new_v <- ifelse(is.na(vals[[lbl]]), "", as.character(vals[[lbl]]))
          is_na <- !nzchar(trimws(new_v))
          changed <- isTRUE(changed_flags[[lbl]])
          
          tags$div(
            tags$div(class = "protocol-field-title", lbl),
            tags$div(
              class = paste(
                "protocol-card",
                if (changed) "changed",
                if (is_na) "na"
              ),
              if (is_na) "—" else new_v
            )
          )
        })
      )
    )
  })

  
  # --- SARA species set ---
  sara_set <- reactive(unique(na.omit(SARA$Scientific.Name)))
  
  # --- AIS species set ---
  ais_set  <- reactive(unique(na.omit(AIS$Scientific.Name)))   # adjust column if needed
  
  # --- toggle states ---
  filter_sara_on <- reactiveVal(FALSE)
  filter_ais_on  <- reactiveVal(FALSE)
  
  observeEvent(input$SARA, {
    new_state <- !isTRUE(filter_sara_on())
    filter_sara_on(new_state)
    
    if (new_state) shinyjs::addClass("SARA", "btn-sara-on")
    else           shinyjs::removeClass("SARA", "btn-sara-on")
  })
  
  observeEvent(input$AIS, {
    new_state <- !isTRUE(filter_ais_on())
    filter_ais_on(new_state)
    
    if (new_state) shinyjs::addClass("AIS", "btn-ais-on")
    else           shinyjs::removeClass("AIS", "btn-ais-on")
  })
  
  apply_interest_filter <- function(spp_vec) {
    out <- spp_vec
    if (isTRUE(filter_sara_on())) out <- intersect(out, sara_set())
    if (isTRUE(filter_ais_on()))  out <- intersect(out, ais_set())
    out
  }
  
  active_filters_label <- reactive({
    labs <- c()
    if (isTRUE(filter_sara_on())) labs <- c(labs, "SARA")
    if (isTRUE(filter_ais_on()))  labs <- c(labs, "AIS")
    if (length(labs) == 0) "None" else paste(labs, collapse = " + ")
  })
  
  output$sara_details <- DT::renderDT({
    if (!isTRUE(filter_sara_on())) {
      return(DT::datatable(
        data.frame(Message = "Click “Species at Risk (SARA) Schedule 1” to view SARA details."),
        rownames = FALSE,
        options = list(dom = "t")
      ))
    }
    
    det <- selected_detections()
    if (is.null(det) || nrow(det) == 0) {
      return(DT::datatable(
        data.frame(Message = "No detections in the current selection."),
        rownames = FALSE,
        options = list(dom = "t")
      ))
    }
    
    # 1) remove geometry + keep only SARA spp
    det <- det %>%
      sf::st_drop_geometry() %>%
      dplyr::mutate(scientificName = as.character(scientificName)) %>%
      dplyr::filter(scientificName %in% sara_set())
    
    if (nrow(det) == 0) {
      return(DT::datatable(
        data.frame(Message = "No SARA Schedule 1 detections for this selection."),
        rownames = FALSE,
        options = list(dom = "t")
      ))
    }
    
    # 2) join Rating + Common name (adjust Common.Name to your real column)
    det2 <- det %>%
      dplyr::left_join(
        SARA %>% dplyr::select(
          Scientific.Name,
          Common.Name,
          Rating
        ),
        by = c("scientificName" = "Scientific.Name")
      )
    
    out <- det2 %>%
      dplyr::group_by(scientificName, Rating) %>%
      dplyr::summarise(
        Common.Name = paste(sort(unique(na.omit(Common.Name))), collapse = " |OR| "),
        n_detections = dplyr::n(),
        n_samples    = dplyr::n_distinct(materialSampleID),
        samples      = paste(sort(unique(materialSampleID)), collapse = ", "),
        files        = paste(sort(unique(na.omit(source_file))), collapse = ", "),
        years        = paste(sort(unique(na.omit(year))), collapse = ", "),
        markers      = paste(sort(unique(na.omit(marker))), collapse = ", "),
        .groups = "drop"
      ) %>%
      dplyr::arrange(Rating, scientificName) %>%
      dplyr::relocate(Common.Name, .after = scientificName)
    
    DT::datatable(
      out,
      rownames = FALSE,
      options = list(pageLength = 10, scrollX = TRUE)
    )
  })
  
  output$ais_details <- DT::renderDT({
    if (!isTRUE(filter_ais_on())) {
      return(DT::datatable(
        data.frame(Message = "Click “Aquatic Invasive Species (AIS)” to view AIS details."),
        rownames = FALSE,
        options = list(dom = "t")
      ))
    }
    
    det <- selected_detections()
    if (is.null(det) || nrow(det) == 0) {
      return(DT::datatable(
        data.frame(Message = "No detections in the current selection."),
        rownames = FALSE,
        options = list(dom = "t")
      ))
    }
    
    det <- det %>%
      sf::st_drop_geometry() %>%
      dplyr::mutate(scientificName = as.character(scientificName)) %>%
      dplyr::filter(scientificName %in% ais_set())
    
    if (nrow(det) == 0) {
      return(DT::datatable(
        data.frame(Message = "No AIS detections for this selection."),
        rownames = FALSE,
        options = list(dom = "t")
      ))
    }
    
    # If AIS has extra columns you want to display, join them here
    # (edit these column names to match your AIS sheet)
    det2 <- det %>%
      dplyr::left_join(
        AIS %>% dplyr::select(Scientific.Name, dplyr::everything()),
        by = c("scientificName" = "Scientific.Name")
      )
    
    out <- det2 %>%
      dplyr::group_by(scientificName) %>%
      dplyr::summarise(
        n_detections = dplyr::n(),
        n_samples    = dplyr::n_distinct(materialSampleID),
        samples      = paste(sort(unique(materialSampleID)), collapse = ", "),
        files        = paste(sort(unique(na.omit(source_file))), collapse = ", "),
        years        = paste(sort(unique(na.omit(year))), collapse = ", "),
        markers      = paste(sort(unique(na.omit(marker))), collapse = ", "),
        .groups = "drop"
      ) %>%
      dplyr::arrange(scientificName)
    
    DT::datatable(out, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE))
  })
  
  
  `%||%` <- function(x, y) if (is.null(x)) y else x
  
  # ---- helper: convert leaflet.draw feature -> sf polygon (EPSG:4326) ----
  feature_to_sf <- function(feature) {
    req(feature$geometry$type)
    type <- feature$geometry$type
    coords <- feature$geometry$coordinates
    
    if (type == "Polygon") {
      ring <- coords[[1]]
      mat  <- do.call(rbind, lapply(ring, function(x) c(x[[1]], x[[2]])))
      poly <- sf::st_polygon(list(mat))
      return(sf::st_sfc(poly, crs = 4326) |> sf::st_sf())
    }
    
    if (type == "MultiPolygon") {
      polys <- lapply(coords, function(poly_i) {
        ring <- poly_i[[1]]
        mat  <- do.call(rbind, lapply(ring, function(x) c(x[[1]], x[[2]])))
        sf::st_polygon(list(mat))
      })
      return(sf::st_sfc(sf::st_multipolygon(polys), crs = 4326) |> sf::st_sf())
    }
    
    stop("Drawn feature type not supported: ", type)
  }
  
  # ---- store the most-recent drawn polygon ----
  drawn_poly <- reactiveVal(NULL)
  
  observeEvent(input$map_draw_new_feature, {
    drawn_poly(feature_to_sf(input$map_draw_new_feature))
  })
  
  observeEvent(input$map_draw_deleted_features, {
    drawn_poly(NULL)
  })
  
  # ---- Year selection (as character or "All") ----
  sel_year_chr <- reactive({
    yr <- input$sel_year %||% "All"
    as.character(yr)
  })
  
  # reactive filtered points
  filtered_sampling_pts <- reactive({
    yr <- sel_year_chr()  # "All" or "2022"/"2023"/"2024"
    
    if (yr == "All") return(sampling_pts)
    
    sampling_pts %>%
      filter(as.character(year) == yr)
  })
  
  # redraw points when year changes
  observeEvent(sel_year_chr(), {
    pts <- filtered_sampling_pts()
    
    leafletProxy("map") %>%
      clearGroup("Sampling points") %>%
      addCircleMarkers(
        data        = pts,
        group       = "Sampling points",
        radius      = 2,
        stroke      = TRUE,
        weight      = 1,
        opacity     = 1,
        fillOpacity = 0.8,
        label       = ~paste0(
          "Marker: ", marker,
          ifelse(is.na(year), "", paste0(" | Year: ", year)),
          ifelse(is.na(source_file), "", paste0(" | File: ", source_file)),
          ifelse(is.na(materialSampleID), "", paste0(" | Sample: ", materialSampleID))
        )
      )
  }, ignoreInit = TRUE)
  
  # ---- 1) Render the leaflet map ONCE ----
  output$map <- renderLeaflet({
    leaflet() %>%
      addProviderTiles(providers$CartoDB.Positron, group = "CartoDB Positron") %>%
      addProviderTiles(providers$Esri.OceanBasemap, group = "Esri Ocean Basemap") %>%
      addProviderTiles(providers$Esri.WorldImagery, group = "Esri World Imagery") %>%
      
      # placeholders: start with 2024 layers (or your preference)
      addPolygons(
        data        = grid_12S_by_year[["2024"]],
        group       = "Species detected by 12S",
        layerId     = ~cell_id,
        fillColor   = ~pal_rich(n_species),
        fillOpacity = ~ifelse(has_sampling, 0.8, 0.04),
        color       = NA,
        label       = ~ifelse(has_sampling,
                              paste("12S richness:", n_species),
                              "No sampling in this cell")
      ) %>%
      addPolygons(
        data        = grid_COI_by_year[["2024"]],
        group       = "Species detected by COI",
        layerId     = ~cell_id,
        fillColor   = ~pal_rich(n_species),
        fillOpacity = ~ifelse(has_sampling, 0.8, 0.04),
        color       = NA,
        label       = ~ifelse(has_sampling,
                              paste("COI richness:", n_species),
                              "No sampling in this cell")
      ) %>%
      addPolygons(
        data        = grid_ALL_static,
        group       = "Species detected by all markers",
        layerId     = ~cell_id,
        fillColor   = ~pal_rich(n_species_total),
        fillOpacity = ~ifelse(has_sampling, 0.8, 0.04),
        color       = NA,
        label       = ~ifelse(has_sampling,
                              paste("Total richness:", n_species_total),
                              "No sampling in this cell")
      ) %>%
      
      addCircleMarkers(
        data        = sampling_pts,
        group       = "Sampling points",
        radius      = 2,
        stroke      = TRUE,
        weight      = 1,
        opacity     = 1,
        fillOpacity = 0.8,
        label       = ~paste0(
          "Marker: ", marker,
          ifelse(is.na(year), "", paste0(" | Year: ", year)),
          ifelse(is.na(source_file), "", paste0(" | File: ", source_file)),
          ifelse(is.na(materialSampleID), "", paste0(" | Sample: ", materialSampleID))
        )
      ) %>%
      
      addPolygons(
        data        = all_polys_zones,
        group       = "MPA/AOI zone boundaries",
        fillOpacity = 0,
        color       = "black",
        weight      = 1,
        opacity     = 0.8,
        options     = pathOptions(clickable = FALSE)
      ) %>%
      
      addPolygons(
        data        = all_polys_click,
        group       = "Total species detected per MPA/AOI",
        layerId     = ~paste(site_type, site_name, sep="||"),
        fillOpacity = 0.01,
        color       = "black",
        weight      = 2,
        opacity     = 1,
        popup       = ~site_name
      ) %>%
      
      addLegend(
        position = "bottomright",
        pal      = pal_rich,
        values   = grid_ALL_static$n_species_total,
        title    = "Species richness from eDNA",
        opacity  = 1,
        className = "legend-richness"              #helps to adjust size of legend
      ) %>%
      
      addDrawToolbar(
        targetGroup = "drawn",
        polygonOptions = drawPolygonOptions(showArea = TRUE),
        rectangleOptions = drawRectangleOptions(),
        polylineOptions = FALSE,
        markerOptions = FALSE,
        circleOptions        = FALSE,
        circleMarkerOptions  = FALSE,
        editOptions = editToolbarOptions()
      ) %>%
      
      addLayersControl(
        baseGroups = c("CartoDB Positron", "Esri Ocean Basemap", "Esri World Imagery"),
        overlayGroups = c(
          "Total species detected per MPA/AOI",
          "Species detected by all markers",
          "Species detected by 12S",
          "Species detected by COI",
          "MPA/AOI zone boundaries",
          "Sampling points"
        ),
        options = layersControlOptions(collapsed = FALSE)
      )
  })
  
  # ---- 2) When the year changes, swap the richness layers ----
  observeEvent(sel_year_chr(), {
    yr <- sel_year_chr()
    yr_map <- if (yr == "All") "2024" else yr
    
    proxy <- leafletProxy("map")
    
    # --- helper to safely add a grid layer ---
    add_grid_layer <- function(data_sf, group_name, value_col, label_prefix) {
      if (is.null(data_sf) || nrow(data_sf) == 0) {
        proxy %>% clearGroup(group_name)
        return(invisible(NULL))
      }
      
      # Ensure has_sampling exists and is logical
      if (!"has_sampling" %in% names(data_sf)) data_sf$has_sampling <- FALSE
      data_sf$has_sampling <- as.logical(data_sf$has_sampling)
      
      # Ensure value_col exists (prevents crashes if a year-layer is missing a column)
      if (!value_col %in% names(data_sf)) data_sf[[value_col]] <- NA_real_
      
      # PRECOMPUTE what leaflet needs (no .data inside leaflet formulas)
      data_sf$.val   <- data_sf[[value_col]]
      data_sf$.label <- ifelse(data_sf$has_sampling,
                               paste0(label_prefix, ": ", data_sf$.val),
                               "No sampling in this cell")
      
      proxy %>%
        clearGroup(group_name) %>%
        addPolygons(
          data        = data_sf,
          group       = group_name,
          layerId     = ~cell_id,
          
          fill        = ~has_sampling,
          fillColor   = ~pal_rich(.val),
          fillOpacity = 0.8,
          
          color       = NA,
          label       = ~.label,
          highlightOptions = highlightOptions(weight = 2, bringToFront = TRUE)
        )
    }
    
    # 12S
    grid12 <- grid_12S_by_year[[yr_map]]
    add_grid_layer(grid12, "Species detected by 12S", "n_species", "12S richness")
    
    # COI
    gridCOI <- grid_COI_by_year[[yr_map]]
    add_grid_layer(gridCOI, "Species detected by COI", "n_species", "COI richness")
    
    # All markers (by year)
    if (exists("grid_ALL_by_year") && is.list(grid_ALL_by_year)) {
      gridALL <- grid_ALL_by_year[[yr_map]]
      add_grid_layer(gridALL, "Species detected by all markers", "n_species_total", "Total richness")
    } else {
      proxy %>% clearGroup("Species detected by all markers")
    }
    
  }, ignoreInit = TRUE)
  
  
  # ---- clicked cell helper ----
  clicked_cell <- reactive({
    click <- input$map_shape_click
    if (is.null(click) || is.null(click$id)) return(NA_integer_)
    suppressWarnings(as.integer(click$id))
  })
  
  selected_detections <- reactive({
    yr <- sel_year_chr()
    
    # start with the full detection points table
    pts <- species_sf_all
    if (yr != "All") pts <- pts %>% dplyr::filter(as.character(year) == yr)
    
    # A) drawn polygon selection takes priority
    poly <- drawn_poly()
    if (!is.null(poly)) {
      inside <- pts[sf::st_within(pts, poly, sparse = FALSE), , drop = FALSE]
      return(inside)
    }
    
    # B) click selection (grid or polygon)
    click <- input$map_shape_click
    if (is.null(click) || is.null(click$id)) return(NULL)
    
    # polygon click case: use your polygon->species list table if you have detection rows by polygon;
    # otherwise fallback to spatial within polygon geometry (recommended)
    if (grepl("\\|\\|", click$id)) {
      parts <- strsplit(click$id, "\\|\\|")[[1]]
      p_type <- parts[1]; p_name <- parts[2]
      
      poly_sel <- all_polys_click %>%
        dplyr::filter(site_type == p_type, site_name == p_name)
      
      if (nrow(poly_sel) == 0) return(NULL)
      
      inside <- pts[sf::st_within(pts, sf::st_geometry(poly_sel), sparse = FALSE), , drop = FALSE]
      return(inside)
    }
    
    # grid cell click: if you have a grid sf with cell_id polygons for the year, use it
    cid <- suppressWarnings(as.integer(click$id))
    if (is.na(cid)) return(NULL)
    
    # pick which grid geometry to use
    yr_map <- if (yr == "All")
      grid_all_year <- grid_ALL_static  # or a year-specific grid if you have one
    
    cell_poly <- grid_all_year %>% dplyr::filter(cell_id == cid)
    if (nrow(cell_poly) == 0) return(NULL)
    
    inside <- pts[sf::st_within(pts, sf::st_geometry(cell_poly), sparse = FALSE), , drop = FALSE]
    inside
  })
  
  
  # ---- species panel (drawn polygon OR click) ----
  output$species_panel <- renderUI({
    
    yr <- sel_year_chr()
    
    # A) drawn polygon selection
    poly <- drawn_poly()
    if (!is.null(poly)) {
      
      # choose which points table drives the list:
      # Use species_sf_all (sf points with scientificName, year, marker/target_gene)
      pts <- species_sf_all
      
      if (yr != "All") {
        pts <- pts %>% dplyr::filter(as.character(year) == yr)
      }
      
      inside <- pts[sf::st_within(pts, poly, sparse = FALSE), , drop = FALSE]
      
      spp <- inside %>%
        dplyr::distinct(scientificName) %>%
        dplyr::arrange(scientificName) %>%
        dplyr::pull(scientificName)
      
      spp <- apply_interest_filter(spp)
      
      if (length(spp) == 0) {
        return(em(if (isTRUE(filter_sara_on()))
          "No SARA Schedule 1 species fall inside the drawn polygon (with the current Year filter)."
          else
            "No species points fall inside the drawn polygon (with the current Year filter)."))
      }
      
      return(tagList(
        strong(paste0(
          "Drawn polygon selection (Year: ", yr, ")",
          if (isTRUE(filter_sara_on())) " — SARA only" else ""
        )),
        tags$br(),
        paste0("Points inside: ", nrow(inside)),
        tags$br(),
        paste0("Unique species: ", length(spp)),
        tags$hr(),
        tags$div(
          style = "max-height: 350px; overflow-y: auto; padding-left: 10px;",
          tags$ul(lapply(spp, tags$li))
        ),
        tags$hr(),
        em("Tip: delete the drawn shape to return to grid/MPA click lists.")
      ))
    }
    
    # B) click-based logic (your existing approach)
    click <- input$map_shape_click
    if (is.null(click) || is.null(click$id)) {
      return(em("Click a grid cell or an MPA/AOI outline, or draw a polygon."))
    }
    
    groups_on <- input$map_groups %||% character(0)
    show_total_grid <- "Species detected by all markers" %in% groups_on
    show_12S_grid   <- "Species detected by 12S" %in% groups_on
    show_COI_grid   <- "Species detected by COI" %in% groups_on
    show_poly_total <- "Total species detected per MPA/AOI" %in% groups_on
    
    make_list <- function(vec) {
      if (length(vec) == 0) return(em("No data available"))
      tags$div(
        style = "max-height: 180px; overflow-y: auto; padding-left: 10px;",
        tags$ul(lapply(vec, tags$li))
      )
    }
    
    # polygon click
    if (grepl("\\|\\|", click$id)) {
      
      if (!show_poly_total) {
        return(em("Turn ON “Total species detected per MPA/AOI” to view polygon species lists."))
      }
      
      parts <- strsplit(click$id, "\\|\\|")[[1]]
      p_type <- parts[1]
      p_name <- parts[2]
      
      row <- poly_species_total %>%
        dplyr::filter(site_type == p_type, site_name == p_name)
      
      spp <- if (nrow(row) == 0) character(0) else row$species[[1]]
      spp <- apply_interest_filter(spp)
      if (length(spp) == 0) return(em("No species records in this MPA/AOI polygon."))
      
      return(tagList(
        strong(paste0(p_type, ": ", p_name)),
        tags$br(),
        paste("Total unique species:", length(spp)),
        tags$hr(),
        tags$div(
          style = "max-height: 350px; overflow-y: auto; padding-left: 10px;",
          tags$ul(lapply(sort(spp), tags$li))
        )
      ))
    }
    
    # grid click
    cid <- suppressWarnings(as.integer(click$id))
    if (is.na(cid)) return(em("Click a grid cell or an MPA/AOI outline."))
    
    # IMPORTANT: year-aware species tables
    # If you have cell_species_2022_12S / cell_species_2023_12S, store them in lists like the grids:
    cell_12S_by_year <- list(
      "2022" = cell_species_2022_12S,
      "2023" = cell_species_2023_12S,
      "2024" = cell_species_2024_12S
    )
    cell_COI_by_year <- list(
      "2024" = cell_species_2024_COI
    )
    
    yr_map <- if (yr == "All") "2024" else yr
    
    spp_total <- cell_species_all %>%
      dplyr::filter(cell_id == cid) %>%
      { if (yr == "All") . else dplyr::filter(., as.character(year) == yr) } %>%
      dplyr::distinct(scientificName) %>%
      dplyr::arrange(scientificName) %>%
      dplyr::pull(scientificName)
    
    spp_total <- apply_interest_filter(spp_total)
    
    cell12_tbl <- cell_12S_by_year[[yr_map]]
    spp_12S <- if (is.null(cell12_tbl)) character(0) else {
      cell12_tbl %>%
        dplyr::filter(cell_id == cid) %>%
        dplyr::distinct(scientificName) %>%
        dplyr::arrange(scientificName) %>%
        dplyr::pull(scientificName)
    }
    
    spp_12S   <- apply_interest_filter(spp_12S)
    
    cellco_tbl <- cell_COI_by_year[[yr_map]]
    spp_COI <- if (is.null(cellco_tbl)) character(0) else {
      cellco_tbl %>%
        dplyr::filter(cell_id == cid) %>%
        dplyr::distinct(scientificName) %>%
        dplyr::arrange(scientificName) %>%
        dplyr::pull(scientificName)
    }
    
    spp_COI   <- apply_interest_filter(spp_COI)
    
    if (!show_total_grid && !show_12S_grid && !show_COI_grid) {
      return(em("Turn ON a richness layer (all markers / 12S / COI) to view the species list for grid cells."))
    }
    
    panels <- tagList()
    
    if (show_total_grid) {
      panels <- tagAppendChildren(
        panels,
        tags$details(open = TRUE,
                     tags$summary(strong(paste0("Total unique species (", length(spp_total), ")"))),
                     make_list(spp_total)
        )
      )
    }
    
    if (show_12S_grid) {
      panels <- tagAppendChildren(
        panels,
        tags$details(open = !show_total_grid,
                     tags$summary(strong(paste0("12S species (", length(spp_12S), ")"))),
                     make_list(spp_12S)
        )
      )
    }
    
    if (show_COI_grid) {
      panels <- tagAppendChildren(
        panels,
        tags$details(open = (!show_total_grid && !show_12S_grid),
                     tags$summary(strong(paste0("COI species (", length(spp_COI), ")"))),
                     make_list(spp_COI)
        )
      )
    }
    
    panels
  })
  
  output$detections_tbl <- DT::renderDT({
    if (!isTRUE(filter_sara_on()) && !isTRUE(filter_ais_on())) {
      return(DT::datatable(
        data.frame(Message = "Click “Species at Risk (SARA)” or “Aquatic Invasive Species (AIS)” to view filtered detections."),
        rownames = FALSE,
        options = list(dom = "t")
      ))
    }
    
    det <- selected_detections()
    if (is.null(det) || nrow(det) == 0) {
      return(DT::datatable(
        data.frame(Message = "No detections in the current selection."),
        rownames = FALSE,
        options = list(dom = "t")
      ))
    }
    
    det0 <- det %>%
      sf::st_drop_geometry() %>%
      dplyr::mutate(scientificName = as.character(scientificName))
    
    keep <- rep(TRUE, nrow(det0))
    if (isTRUE(filter_sara_on())) keep <- keep & det0$scientificName %in% sara_set()
    if (isTRUE(filter_ais_on()))  keep <- keep & det0$scientificName %in% ais_set()
    
    det_f <- det0[keep, , drop = FALSE]
    
    if (nrow(det_f) == 0) {
      return(DT::datatable(
        data.frame(Message = paste0("No detections match the current filter (", active_filters_label(), ").")),
        rownames = FALSE,
        options = list(dom = "t")
      ))
    }
    
    det_show <- det_f %>%
      dplyr::select(
        scientificName, materialSampleID, eventDate, year, month,
        marker, target_gene, organismQuantity, organismQuantityType,
        sampleSizeValue, sampleSizeUnit, source_file, occurrenceID
      ) %>%
      dplyr::arrange(scientificName, year, materialSampleID)
    
    DT::datatable(det_show, rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE))
  })
  
  
  # highlight selected cell (outline on top)
  observeEvent(input$map_shape_click, {
    cid <- clicked_cell()
    if (is.na(cid)) return(NULL)  # <- stops crash when click$id is NULL
    
    sel <- richness_grid_all_poly %>% filter(cell_id == cid)
    
    leafletProxy("map") %>%
      clearGroup("selected_cell") %>%
      addPolygons(
        data = sel,
        group = "selected_cell",
        fillOpacity = 0,
        color = "white",
        weight = 3
      )
  }, ignoreInit = TRUE)
  
}
