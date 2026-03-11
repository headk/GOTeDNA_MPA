library(sf)
library(dplyr)
library(leaflet)
library(leaflet.extras)
library(arcgislayers)
library(tidyr)
library(ggplot2)
library(scico)
library(wesanderson)
library(htmltools)
library(htmlwidgets)
library(DT)
library(shinyjs)
library(stringr)
library(openxlsx)
library(purrr)
library(readr)
library(robis)
library(bslib)
library(shiny)
library(phyloseq)
library(taxplore)  #Here is the link for taxplore tutorial once data is linked from OBIS: https://markschl.github.io/taxplore/articles/tutorial.html#shiny-apps
library(plotly)
library(worrms)
library(vegan)
library(forcats)
library(shinycssloaders)

#Canadian MPA Network Shapefiles and planning regions extracted from the Canadian Database of Protected and Conserved Areas. This data was filtered for the three MPAs targeted in the Maritimes region.

#Function to find St. Anns Bank, Musquash, Gully (MPA) and Eastern Shore Islands network, Fundian Channel-BrownsBank (AOI) polygons + transform to WGS84 for leaflet
read_poly_wgs84 <- function(...) {
  st_read(file.path("data", "polygons", ...), quiet = TRUE) %>%
    st_transform(4326)
}

mpa_targets <- read_poly_wgs84("mpa_targets.shp")
esi_poly    <- read_poly_wgs84("EasternShoreIslands_networksite.shp")
fcbb_poly   <- read_poly_wgs84("FCBB_Proposed_MPA_Boundary_zones.shp")

#Combine all polygons into one sf object
#Add a common name/type to each polygon set
mpa_polys <- mpa_targets %>%
  st_geometry() %>%                     # keep just geometry
  st_as_sf() %>%                        # convert back to sf
  mutate(
    site_name = mpa_targets$NAME_E,
    site_type = "MPA"
  )

esi_polys <- esi_poly %>%
  st_geometry() %>%
  st_as_sf() %>%
  mutate(
    site_name = "Eastern Shore Islands AOI",
    site_type = "AOI"
  )

fcbb_polys <- fcbb_poly %>%
  st_geometry() %>%
  st_as_sf() %>%
  mutate(
    site_name = "Fundian Channel–Browns Bank AOI",
    site_type = "AOI"
  )

#Bind everything together
all_polys <- bind_rows(
  mpa_polys,
  esi_polys,
  fcbb_polys
)

#Ensure it's still sf
all_polys <- st_as_sf(all_polys)

# ---- Clean once ----
all_polys_clean <- all_polys %>%
  st_make_valid() %>%
  st_collection_extract("POLYGON") %>%
  (\(x) x[!st_is_empty(x), ])() %>%
  st_as_sf()

# ---- A) zone outlines (keep every piece) ----
all_polys_zones <- all_polys_clean

# ---- B) clickable dissolved outlines (one feature per site) ----
all_polys_click <- all_polys_clean %>%
  group_by(site_type, site_name) %>%
  summarise(geometry = st_union(x), .groups = "drop") %>%
  st_as_sf()


#-----------------------------------------------------------------------------
read_data <- function(
    dataset_ids       = NULL,
    scientificname    = NULL,
    worms_id          = NULL,
    areaid            = "34", #Canada: North Atlantic
    #geometry          = "POLYGON ((-67.72 40.614, -56.821 40.614, -56.821 47.279, -67.72 47.279, -67.72 40.614))",
    join_by           = c("auto", "occurrenceID", "id"),
    require_absences  = TRUE
) {

  join_by <- match.arg(join_by)

  # ---- columns you want back ----
  occurrence_cols <- c(
    "recordedBy","bibliographicCitation","materialSampleID",
    "organismQuantity","organismQuantityType",
    "samp_size", "samp_size_unit", "decimalLatitude", "decimalLongitude",
    "minimumDepthInMeters","maximumDepthInMeters","month","year",
    "scientificNameID","kingdom","phylum","class","order","family","genus",
    "dataset_id", "bathymetry", "associatedSequences", "bibliographicCitation"
  )

  dna_cols <- c(
    "id","dna_sequence","target_gene","pcr_primer_forward", "pcr_primer_reverse",
    "samp_name", "env_broad_scale","env_local_scale","env_medium","samp_mat_process",
    "size_frac","samp_size","samp_size_unit","otu_db","seq_kit", "otu_seq_comp_appr",
    "pcr_primer_name_forward","pcr_primer_name_reverse", "pcr_primer_reference",
    "occurrenceID"
  )

  mof_cols <- c(
    "id","seq_id","samp_category","checkls_ver","assay_name","assay_type",
    "targetTaxonomicAssay","geo_loc_name","technical_rep_id","project_contact",
    "seq_run_id","lib_id","project_id","pcr_0_1","samp_store_sol","samp_store_temp",
    "platform","instrument","tax_assign_cat","LClabel","occurrenceID",
    "nucl_acid_ext","nucl_acid_ext_kit","filter_material"
  )

  required_ext_cols <- c(
    "samp_name",
    "target_gene",
    "pcr_primer_name_forward",
    "pcr_primer_name_reverse"
  )

  added_cols <- c("category", "flags")

  mandatory_obis <- c(
    "occurrenceID","eventDate","decimalLongitude","decimalLatitude",
    "scientificName","occurrenceStatus","basisOfRecord"
  )

  cols_included_from_OBIS <- unique(c(
    occurrence_cols, dna_cols, mof_cols, added_cols, mandatory_obis
  ))

  # ---- helper: enforce columns (no external function needed) ----
  enforce_cols <- function(occ_all, cols) {
    # add missing columns as NA
    missing <- setdiff(cols, names(occ_all))
    if (length(missing) > 0) {
      for (m in missing) occ_all[[m]] <- NA
    }
    # keep only requested columns in a consistent order
    occ_all <- occ_all[, intersect(cols, names(occ_all)), drop = FALSE]
    occ_all
  }

  # ---- helper: build extension-wide table (MOF wide + DNA) ----
  join_extensions <- function(rec_ext, cols_included_from_OBIS, join_by) {
    if (is.null(rec_ext) || nrow(rec_ext) == 0L) return(NULL)

    # base cores
    core_occ_ext <- dplyr::distinct(rec_ext, occurrenceID, .keep_all = TRUE)
    core_id_ext  <- dplyr::distinct(rec_ext, id, .keep_all = TRUE)

    # DNADerivedData extension
    dna_only <- robis::unnest_extension(rec_ext, "DNADerivedData")
    shared_dna_cols <- intersect(cols_included_from_OBIS, names(dna_only))
    dna_only <- dplyr::select(dna_only, dplyr::any_of(shared_dna_cols))

    # MeasurementOrFact extension -> wide
    mof_only <- robis::unnest_extension(rec_ext, "MeasurementOrFact") %>%
      dplyr::group_by(occurrenceID, measurementType) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()

    wide_mof <- tidyr::pivot_wider(
      mof_only,
      id_cols    = c(occurrenceID, id),
      names_from = measurementType,
      values_from = measurementValue
    )

    shared_mof_cols <- intersect(cols_included_from_OBIS, names(wide_mof))
    wide_mof <- dplyr::select(wide_mof, dplyr::any_of(shared_mof_cols))

    mof_and_dna <- dplyr::left_join(wide_mof, dna_only, by = "id")

    # choose join key
    join_choice <- join_by
    if (join_choice == "auto") {
      can_occ <- "occurrenceID" %in% names(core_occ_ext) &&
        "occurrenceID" %in% names(mof_and_dna) &&
        any(!is.na(core_occ_ext$occurrenceID)) &&
        any(!is.na(mof_and_dna$occurrenceID))
      can_id  <- "id" %in% names(core_id_ext) &&
        "id" %in% names(mof_and_dna) &&
        any(!is.na(core_id_ext$id)) &&
        any(!is.na(mof_and_dna$id))
      if (can_occ) join_choice <- "occurrenceID"
      else if (can_id) join_choice <- "id"
      else stop("Neither occurrenceID nor id can be used to join core and extensions.")
    }

    out_ext <- if (join_choice == "occurrenceID") {
      dplyr::left_join(core_occ_ext, mof_and_dna, by = "occurrenceID")
    } else {
      dplyr::left_join(core_id_ext, mof_and_dna, by = "id")
    }

    out_ext <- dplyr::select(out_ext, dplyr::any_of(cols_included_from_OBIS))
    out_ext
  }

  # ---- 0) discover dataset_ids if NULL ----
  if (is.null(dataset_ids)) {
    message("Discovering datasets with DNADerivedData and MeasurementOrFact extensions ...")
    ds_tbl <- robis::dataset(
      scientificname = scientificname,
      areaid         = areaid,
      #geometry       = geometry,
      taxonid        = worms_id,
      hasextensions  = c("DNADerivedData", "MeasurementOrFact")
    ) %>%
      dplyr::filter(statistics$absence != 0)

    if (nrow(ds_tbl) == 0L) {
      warning("No datasets found that match scientificname/areaid AND have DNADerivedData + absences.")
      return(tibble::tibble())
    }

    if ("id" %in% names(ds_tbl)) dataset_ids <- unique(ds_tbl$id)
    else if ("datasetid" %in% names(ds_tbl)) dataset_ids <- unique(ds_tbl$datasetid)
    else stop("Could not find dataset id column in dataset() output.")

    message("Found ", length(dataset_ids), " dataset(s).")
  }

  dataset_ids <- as.character(dataset_ids)

  exclude_list <- c("NO_COORD","ZERO_COORD","LON_OUT_OF_RANGE","LAT_OUT_OF_RANGE","NO_MATCH")

  # ---- 1) loop datasets ----
  obis_list <- purrr::map(dataset_ids, function(ds) {

    Sys.sleep(3)
    message("Pulling OBIS dataset: ", ds)

    # defensive: check extensions exist
    ds_meta <- robis::dataset(datasetid = ds)
    exts <- tolower(unlist(ds_meta$extensions))
    if (!"dnaderiveddata" %in% exts) {
      warning("Dataset ", ds, " has no DNADerivedData extension; skipping.")
      return(NULL)
    }
    if (!"measurementorfact" %in% exts) {
      warning("Dataset ", ds, " has no MeasurementOrFact extension; skipping.")
      return(NULL)
    }

    # 1a) FULL CORE with absences (NO extensions)
    core_all <- tryCatch(
      robis::occurrence(
        datasetid      = ds,
        scientificname = scientificname,
        taxonid        = worms_id,
        areaid         = areaid,
        #geometry       = geometry,
        absence        = "include",
        dropped        = "include",
        exclude        = exclude_list
      ),
      error = function(e) {
        warning("Failed to fetch CORE (absences) for dataset ", ds, ": ", conditionMessage(e))
        NULL
      }
    )

    if (is.null(core_all) || nrow(core_all) == 0L) {
      warning("No core occurrence records returned for dataset ", ds, ".")
      return(NULL)
    }

    core_all <- dplyr::distinct(core_all, occurrenceID, .keep_all = TRUE)

    if (!"occurrenceStatus" %in% names(core_all)) {
      warning("Dataset ", ds, " has no occurrenceStatus column; skipping.")
      return(NULL)
    }

    status_vals <- unique(na.omit(core_all$occurrenceStatus))
    if (require_absences && !all(c("present","absent") %in% status_vals)) {
      warning("Dataset ", ds, " does not contain both 'present' and 'absent'; skipping.")
      return(NULL)
    }

    core_all <- dplyr::select(core_all, dplyr::any_of(cols_included_from_OBIS))
    core_all <- enforce_cols(core_all, cols_included_from_OBIS)

    # 1b) EXTENSIONS (present-only; do NOT request absence="include")
    rec_ext <- tryCatch(
      robis::occurrence(
        datasetid      = ds,
        scientificname = scientificname,
        taxonid        = worms_id,
        areaid         = areaid,
        #geometry       = geometry,
        extensions     = c("DNADerivedData", "MeasurementOrFact"),
        hasextensions  = c("DNADerivedData", "MeasurementOrFact"),
        dropped        = "include",
        exclude        = exclude_list
      ),
      error = function(e) {
        warning("Failed to fetch EXTENSIONS for dataset ", ds, ": ", conditionMessage(e))
        NULL
      }
    )

    # If extensions fetch fails, return core (absences preserved)
    if (is.null(rec_ext) || nrow(rec_ext) == 0L) {
      core_all$samp_name <- as.character(dplyr::coalesce(core_all$samp_name, core_all$materialSampleID))
      return(core_all)
    }

    ext_joined <- join_extensions(rec_ext, cols_included_from_OBIS, join_by)

    ext_joined <- ext_joined %>%
      dplyr::distinct(.data$occurrenceID, .keep_all = TRUE)

    if (is.null(ext_joined) || nrow(ext_joined) == 0L) {
      core_all$samp_name <- as.character(dplyr::coalesce(core_all$samp_name, core_all$materialSampleID))
      return(core_all)
    }

    missing_required <- setdiff(required_ext_cols, names(ext_joined))
    if (length(missing_required) > 0) {
      warning(
        "Dataset ", ds,
        " missing required column(s): ",
        paste(missing_required, collapse = ", "),
        " ; skipping."
      )
      return(NULL)
    }

    # 1c) merge extensions onto full core (absences keep NA extension fields)
    out <- dplyr::left_join(core_all, ext_joined, by = "occurrenceID", suffix = c("", ".ext"))

    # coalesce duplicated columns (prefer extension values)
    dup_cols <- intersect(names(core_all), names(ext_joined))
    dup_cols <- setdiff(dup_cols, "occurrenceID")
    for (nm in dup_cols) {
      ext_nm <- paste0(nm, ".ext")
      if (ext_nm %in% names(out)) {
        out[[nm]] <- dplyr::coalesce(out[[ext_nm]], out[[nm]])
        out[[ext_nm]] <- NULL
      }
    }

    out <- dplyr::select(out, dplyr::any_of(cols_included_from_OBIS))
    out <- enforce_cols(out, cols_included_from_OBIS)
    out$samp_name <- as.character(dplyr::coalesce(out$samp_name, out$materialSampleID))
    out
  })

  obis_list <- purrr::compact(obis_list)

  if (length(obis_list) == 0L) {
    warning("No OBIS datasets returned any records for these filters.")
    return(tibble::tibble())
  }

  GOTeDNA_occ_all <- dplyr::bind_rows(obis_list)
  rownames(GOTeDNA_occ_all) <- NULL
  GOTeDNA_occ_all
}

STORED_DATA <- read_data()
saveRDS(STORED_DATA, "./data/test_read_data_file.rds")
STORED_DATA <- readRDS("./data/test_read_data_file.rds")

################DO NOT CHANGE ABOVE CODE



#Connect code below to stored data object

# ---- standardize types early ----
occ_all <- STORED_DATA %>%
  dplyr::mutate(
    year             = as.character(year),
    samp_name         = as.character(samp_name),
    occurrenceStatus  = tolower(as.character(occurrenceStatus)),
    decimalLatitude   = suppressWarnings(as.numeric(decimalLatitude)),
    decimalLongitude  = suppressWarnings(as.numeric(decimalLongitude)),
    target_gene       = dplyr::case_when(
      stringr::str_detect(tolower(target_gene), "12s") ~ "12S",
      stringr::str_detect(tolower(target_gene), "coi") ~ "COI",
      stringr::str_detect(tolower(target_gene), "16s") ~ "16S",
      stringr::str_detect(tolower(target_gene), "18s") ~ "18S",
      TRUE ~ as.character(target_gene)
    )
  )

# ---- Build available gene-year keys dynamically ----
KEY_TBL <- occ_all %>%
  dplyr::filter(!is.na(year), year != "", !is.na(target_gene), target_gene != "") %>%
  dplyr::distinct(target_gene, year) %>%
  dplyr::mutate(
    target_gene = as.character(target_gene),
    year        = as.character(year),
    key         = paste(target_gene, year, sep = "_")
  ) %>%
  dplyr::arrange(target_gene, year)

KEYS <- KEY_TBL$key

# ---- Split the *data* by key ----
DATA_BY_KEY <- occ_all %>%
  filter(!is.na(year), year != "", !is.na(target_gene), target_gene != "") %>%
  mutate(
    target_gene = as.character(target_gene),
    year        = as.character(year),
    key         = paste(target_gene, year, sep = "_")
  ) %>%
  group_by(key) %>%
  group_split(.keep = TRUE)

names(DATA_BY_KEY) <- KEYS

# ---- Convert each group to SF points (PRESENT only) ----
points_sf_from_occ_all <- function(occ_all) {
  occ_all %>%
    filter(!is.na(decimalLatitude), !is.na(decimalLongitude),
           tolower(as.character(occurrenceStatus)) == "present") %>%
    st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
}

SPECIES_SF_BY_KEY <- purrr::imap(DATA_BY_KEY, ~{
  parts <- strsplit(.y, "_")[[1]]
  gene <- parts[1]
  yr   <- parts[2]

  points_sf_from_occ_all(.x) %>%
    mutate(target_gene = gene, year = yr)
})

species_sf_all <- dplyr::bind_rows(SPECIES_SF_BY_KEY)


# --- Sampling points layer for leaflet (keeps year + source_file if present) ---

sampling_pts <- species_sf_all %>%
  mutate(
    year        = if ("year" %in% names(.)) as.character(year) else NA_character_,
    #source_file = if ("source_file" %in% names(.)) source_file else NA
  ) %>%
  distinct(target_gene, year, samp_name, geometry, .keep_all = TRUE)


#############################################################
#Turn species data into sf points  for the spatial join


# ---- ONE join: species points -> clickable MPA/AOI polygons ----
species_in_polys_all <- sf::st_join(
  species_sf_all,
  all_polys_click,          # dissolved clickable polygons
  join = sf::st_within,
  left = FALSE
) %>%
  dplyr::mutate(
    year = as.character(year),
    target_gene = as.character(target_gene),
    scientificName = as.character(scientificName)
  ) %>%
  dplyr::filter(!is.na(scientificName), scientificName != "")

# year-aware species list per polygon
poly_species_year <- species_in_polys_all %>%
  sf::st_drop_geometry() %>%
  dplyr::distinct(site_name, site_type, year, scientificName) %>%
  dplyr::group_by(site_name, site_type, year) %>%
  dplyr::summarise(species = list(sort(unique(scientificName))), .groups = "drop")

# all-years species list per polygon
poly_species_all <- species_in_polys_all %>%
  sf::st_drop_geometry() %>%
  dplyr::distinct(site_name, site_type, scientificName) %>%
  dplyr::group_by(site_name, site_type) %>%
  dplyr::summarise(species = list(sort(unique(scientificName))), .groups = "drop")

total_species_year <- species_in_polys_all %>%
  sf::st_drop_geometry() %>%
  dplyr::group_by(site_name, site_type, year) %>%
  dplyr::summarise(n_species_total = dplyr::n_distinct(scientificName), .groups = "drop")

total_species_all <- species_in_polys_all %>%
  sf::st_drop_geometry() %>%
  dplyr::group_by(site_name, site_type) %>%
  dplyr::summarise(n_species_total = dplyr::n_distinct(scientificName), .groups = "drop")

species_by_class_year <- species_in_polys_all %>%
  sf::st_drop_geometry() %>%
  dplyr::group_by(site_name, site_type, year, class) %>%
  dplyr::summarise(n_species = dplyr::n_distinct(scientificName), .groups = "drop")

species_by_class_wide_year <- species_by_class_year %>%
  tidyr::pivot_wider(names_from = class, values_from = n_species, values_fill = 0)

species_by_class_all <- species_in_polys_all %>%
  sf::st_drop_geometry() %>%
  dplyr::group_by(site_name, site_type, class) %>%
  dplyr::summarise(n_species = dplyr::n_distinct(scientificName), .groups = "drop")

species_by_class_wide_all <- species_by_class_all %>%
  tidyr::pivot_wider(names_from = class, values_from = n_species, values_fill = 0)






# optional convenience table for "All years combined"

species_in_polys_all %>%
  st_drop_geometry() %>%
  count(site_name, site_type, year, target_gene)

poly_species_all_from_year <- poly_species_year %>%
  group_by(site_name, site_type) %>%
  summarise(species = list(sort(unique(unlist(species)))), .groups = "drop")


#Summary report per polygon




##Species Richness Polygons

# 1) Make a grid over polygons
# 1) Work in a projected CRS (Canada Lambert is a good default)
##Species Richness Polygons (Option A: projected CRS to avoid s2 errors)

##Species Richness Polygons (Option A: build grid in EPSG:4326 so it stays upright in leaflet)

# 0) Work in lon/lat (leaflet native)
crs_ll <- 4326

# 1) Clean + dissolve polygons in 4326
poly_union_ll <- all_polys_click %>%
  sf::st_make_valid() %>%
  sf::st_transform(crs_ll) %>%
  sf::st_union() %>%
  sf::st_as_sf()

# 2) Choose an "approx 2000 m" grid size expressed in degrees at your latitude
cell_m <- 2000

# centroid latitude (used to approximate meters->degrees)
cent <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(poly_union_ll)))
lat0 <- mean(cent[,2], na.rm = TRUE)

deg_per_m_lat <- 1 / 111320
deg_per_m_lon <- 1 / (111320 * cos(lat0 * pi/180))

cellsize_deg <- c(cell_m * deg_per_m_lon, cell_m * deg_per_m_lat)  # c(lon_deg, lat_deg)

# 3) Build grid in 4326 (upright in leaflet)
grid_ll <- sf::st_make_grid(
  poly_union_ll,
  cellsize = cellsize_deg,
  square   = TRUE
) %>%
  sf::st_as_sf() %>%
  dplyr::mutate(cell_id = dplyr::row_number())

# 4) Clip grid ONCE (still 4326)
grid_clip_ll <- sf::st_intersection(grid_ll, poly_union_ll) %>%
  sf::st_make_valid() %>%
  sf::st_collection_extract("POLYGON") %>%
  (\(x) x[!sf::st_is_empty(x), ])() %>%
  sf::st_as_sf()

# (Leaflet uses 4326 anyway)
grid_clip <- grid_clip_ll

# 5) Points stay in 4326 too
SPECIES_SF_BY_KEY_ll <- SPECIES_SF_BY_KEY %>%
  purrr::map(~sf::st_transform(.x, crs_ll))

#Make an "all points" sf in 4326 (PRESENT only)
species_sf_all_ll <- dplyr::bind_rows(SPECIES_SF_BY_KEY_ll)

# 6) Richness in 4326 (no CRS mismatch)
make_richness_layer_fast <- function(grid_sf, pts_sf) {
  stopifnot(inherits(grid_sf, "sf"), inherits(pts_sf, "sf"))
  if (!"scientificName" %in% names(pts_sf)) {
    stop("pts_sf must contain a 'scientificName' column.")
  }

  idx <- sf::st_intersects(grid_sf, pts_sf)

  n_sp <- vapply(idx, function(i) {
    if (length(i) == 0) return(0L)
    length(unique(as.character(pts_sf$scientificName[i])))
  }, integer(1))

  grid_sf %>%
    dplyr::mutate(
      n_species    = n_sp,
      has_sampling = n_sp > 0
    )
}

RICHNESS_BY_KEY <- purrr::map(
  SPECIES_SF_BY_KEY_ll,
  ~make_richness_layer_fast(grid_clip_ll, .x)
)

# optional: all-years combined richness
RICHNESS_ALL <- make_richness_layer_fast(
  grid_clip_ll,
  dplyr::bind_rows(SPECIES_SF_BY_KEY_ll)
)

# species_sf_all: sf POINTS with at least cell_id, target_gene, scientificName (and year)
# grid_clip: sf POLYGONS with cell_id + geometry  (your analysis grid)

build_gene_all_grid <- function(gene, grid_sf, pts_sf, tax_col = "scientificName") {

  pts_g <- pts_sf %>%
    dplyr::filter(as.character(target_gene) == gene) %>%
    dplyr::mutate(
      cell_id = as.integer(cell_id),
      taxon   = as.character(.data[[tax_col]])
    ) %>%
    dplyr::filter(!is.na(cell_id), !is.na(taxon), taxon != "")

  # unique taxa per cell across ALL years
  rich_tbl <- pts_g %>%
    sf::st_drop_geometry() %>%
    dplyr::distinct(cell_id, taxon) %>%
    dplyr::count(cell_id, name = "n_species")

  out <- grid_sf %>%
    dplyr::mutate(cell_id = as.integer(cell_id)) %>%
    dplyr::left_join(rich_tbl, by = "cell_id") %>%
    dplyr::mutate(
      n_species    = tidyr::replace_na(n_species, 0L),
      has_sampling = n_species > 0
    )

  out
}

# Use the lon/lat versions you created for the grid + points
# grid_clip_ll and species_sf_all_ll are both EPSG:4326

species_sf_all_cell <- sf::st_join(
  species_sf_all_ll,
  grid_clip_ll %>% dplyr::select(cell_id),
  join = sf::st_within,
  left = FALSE
)

# Now this will work because pts_sf has cell_id
RICHNESS_GENE_ALL <- list(
  "12S" = build_gene_all_grid("12S", grid_clip_ll, species_sf_all_cell, tax_col = "scientificName"),
  "COI" = build_gene_all_grid("COI", grid_clip_ll, species_sf_all_cell, tax_col = "scientificName"),
  "16S" = build_gene_all_grid("16S", grid_clip_ll, species_sf_all_cell, tax_col = "scientificName"),
  "18S" = build_gene_all_grid("18S", grid_clip_ll, species_sf_all_cell, tax_col = "scientificName")
)

# optional: keep dissolved polygon union for leaflet
poly_union <- poly_union_ll


## ===== Unified richness palettes + leaflet map (Option A: 4326-only) =====

# Assumes these already exist from your Option A block:
# - grid_clip_ll (sf, EPSG:4326) with column cell_id
# - SPECIES_SF_BY_KEY_ll (named list of sf points, EPSG:4326) with scientificName, year, target_gene
# - RICHNESS_BY_KEY (named list of sf grids, EPSG:4326) with n_species
# - RICHNESS_ALL (sf grid, EPSG:4326) with n_species (from make_richness_layer_fast)
# - poly_union_ll or all_polys_click for outlines

# 0) Ensure grid has cell_id
if (!"cell_id" %in% names(grid_clip_ll)) {
  grid_clip_ll <- grid_clip_ll %>% dplyr::mutate(cell_id = dplyr::row_number())
}

# 1) Make an "all points" sf in 4326 (present only; should already be present-only, but keep safe)
species_sf_all_ll <- dplyr::bind_rows(SPECIES_SF_BY_KEY_ll) %>%
  dplyr::mutate(
    year             = as.character(year),
    target_gene      = as.character(target_gene),
    scientificName   = as.character(scientificName),
    occurrenceStatus = tolower(as.character(occurrenceStatus))
  ) %>%
  dplyr::filter(
    is.na(occurrenceStatus) | occurrenceStatus == "present",
    !is.na(scientificName), scientificName != ""
  )

# 2) Rename ALL grid column consistently for mapping (“total across whatever you used to compute it”)
# Your RICHNESS_ALL currently has n_species from make_richness_layer_fast().
RICHNESS_ALL <- RICHNESS_ALL %>%
  dplyr::rename(n_species_total = n_species)

# (Optional) If you want per-year ALL-markers layers in 4326:
years_all <- sort(unique(na.omit(species_sf_all_ll$year)))

RICHNESS_ALL_BY_YEAR <- setNames(
  lapply(years_all, function(yr) {
    pts_y <- species_sf_all_ll %>% dplyr::filter(year == yr)
    make_richness_layer_fast(grid_clip_ll, pts_y) %>%
      dplyr::rename(n_species_total = n_species)
  }),
  years_all
)

# 3) Shared richness palette domain across ALL layers
max_rich <- max(
  RICHNESS_ALL$n_species_total,
  unlist(purrr::map(RICHNESS_BY_KEY, ~ .x$n_species)),
  unlist(purrr::map(RICHNESS_ALL_BY_YEAR, ~ .x$n_species_total)),
  na.rm = TRUE
)

rich_domain <- c(0, max_rich)

wes_cont <- function(name, n = 100) {
  grDevices::colorRampPalette(
    wesanderson::wes_palette(name, type = "continuous")
  )(n)
}

pal_vec <- wes_cont("Zissou1", 100)

pal_rich <- leaflet::colorNumeric(
  palette  = pal_vec,
  domain   = rich_domain,
  na.color = "transparent"
)

# 4) Cell -> species list tables (ALL + by key) using *grid_clip_ll* (not projected!)
idx_all <- sf::st_intersects(grid_clip_ll, species_sf_all_ll)

CELL_SPECIES_ALL <- tibble::tibble(cell_id = grid_clip_ll$cell_id) %>%
  dplyr::mutate(
    spp = purrr::map(idx_all, \(i) sort(unique(species_sf_all_ll$scientificName[i])))
  )

CELL_SPECIES_BY_KEY <- purrr::imap(SPECIES_SF_BY_KEY_ll, ~{
  idx <- sf::st_intersects(grid_clip_ll, .x)
  tibble::tibble(cell_id = grid_clip_ll$cell_id) %>%
    dplyr::mutate(
      spp = purrr::map(idx, \(i) sort(unique(.x$scientificName[i])))
    )
})

# Long table: one row per (cell_id, scientificName, target_gene, year)
cell_species_all <- purrr::imap_dfr(SPECIES_SF_BY_KEY_ll, ~{
  key <- .y
  parts <- strsplit(key, "_")[[1]]
  gene <- parts[1]
  yr   <- parts[2]

  idx <- sf::st_intersects(grid_clip_ll, .x)

  tibble::tibble(cell_id = grid_clip_ll$cell_id) %>%
    dplyr::mutate(scientificName = purrr::map(idx, \(i) unique(.x$scientificName[i]))) %>%
    tidyr::unnest(scientificName) %>%
    dplyr::mutate(target_gene = gene, year = yr) %>%
    dplyr::filter(!is.na(scientificName), scientificName != "")
})

cell_species_total <- cell_species_all %>%
  dplyr::group_by(cell_id) %>%
  dplyr::summarise(n_species_total = dplyr::n_distinct(scientificName), .groups = "drop")

# 5) Convenience: split RICHNESS_BY_KEY into gene buckets
grid_12S_by_year <- RICHNESS_BY_KEY[grep("^12S_", names(RICHNESS_BY_KEY))]
grid_COI_by_year <- RICHNESS_BY_KEY[grep("^COI_", names(RICHNESS_BY_KEY))]
grid_16S_by_year <- RICHNESS_BY_KEY[grep("^16S_", names(RICHNESS_BY_KEY))]
grid_18s_by_year <- RICHNESS_BY_KEY[grep("^18S_", names(RICHNESS_BY_KEY))]

# 6) Helper: pick which richness layer to display
# gene: "All" / "12S" / "COI" / "16S" / "18S"
# year: "All" or specific year string
get_richness_layer <- function(gene = "All", year = "All") {
  gene <- as.character(gene)
  year <- as.character(year)

  if (gene == "All") {
    if (year == "All") return(RICHNESS_ALL)
    if (year %in% names(RICHNESS_ALL_BY_YEAR)) return(RICHNESS_ALL_BY_YEAR[[year]])
    return(RICHNESS_ALL)
  }

  # gene-specific layers are only defined for specific years (no "All-years" gene grid)
  if (year == "All") return(NULL)

  key <- paste(gene, year, sep = "_")
  if (key %in% names(RICHNESS_BY_KEY)) return(RICHNESS_BY_KEY[[key]])

  NULL
}

# 7) Defaults
default_gene <- "All"   # or "12S", "COI"
default_year <- "All"   # or "2024", etc.

init_layer <- get_richness_layer(default_gene, default_year)

# Optionally: a safe init if gene-specific All-years returns NULL
if (is.null(init_layer)) init_layer <- RICHNESS_ALL



#Define polygon layers

# Richness layers keyed by gene+year
selected_rich_key <- function(gene, year) {
  gene <- as.character(gene); year <- as.character(year)

  if (gene == "All") return("All")        # handle via RICHNESS_ALL / RICHNESS_ALL_BY_YEAR
  if (year == "All") return(NULL)         # no gene-specific layer for All-years
  paste(gene, year, sep = "_")            # e.g., "12S_2024"
}


#Organize per year layers into names lists

# --- Richness grids by year ---
grid_12S_by_year <- RICHNESS_BY_KEY[grep("^12S_", names(RICHNESS_BY_KEY))]
grid_COI_by_year <- RICHNESS_BY_KEY[grep("^COI_", names(RICHNESS_BY_KEY))]
grid_16S_by_year <- RICHNESS_BY_KEY[grep("^16S_", names(RICHNESS_BY_KEY))]
grid_18S_by_year <- RICHNESS_BY_KEY[grep("^18S_", names(RICHNESS_BY_KEY))]

grid_ALL_static <- RICHNESS_ALL


## ===== Shared richness colour scale across ALL layers =====

# 1) Collect maxima safely (handles NULLs / missing cols)
max_from_key_layers <- function(lst, col) {
  vals <- unlist(purrr::map(lst, ~{
    if (is.null(.x) || !col %in% names(.x)) return(NA_real_)
    suppressWarnings(as.numeric(.x[[col]]))
  }))
  suppressWarnings(max(vals, na.rm = TRUE))
}

max_all_markers <- max(
  suppressWarnings(max(as.numeric(RICHNESS_ALL$n_species_total), na.rm = TRUE)),
  max_from_key_layers(RICHNESS_ALL_BY_YEAR, "n_species_total"),
  na.rm = TRUE
)

max_gene_layers <- max_from_key_layers(RICHNESS_BY_KEY, "n_species")

max_rich <- max(max_all_markers, max_gene_layers, na.rm = TRUE)

# 2) Build ONE palette with ONE domain
rich_domain <- c(0, max_rich)


##NOTE: Right now the SARA Schedule 1 filter is matching based on the data inputted into the app. Once linking the code to OBIS data, edit so that it matches based on WoRMS AphiaID


#Code for the SARA Schedule 1 list (Cleanup and WoRMS linkage)

#Call in file from OBIS_Prep folder
#SARA <- read.xlsx(file.path("data", "sara_ais", "SARA_Clean_Schedule1_specieslist.xlsx"))

SARA <- read.xlsx(file.path("data", "sara_ais", "SARA_Clean_Schedule1_specieslist_noPacific.xlsx"))
AIS  <- read.xlsx(file.path("data", "sara_ais", "Target_AIS_List_Claudio.xlsx"))

# helper: safe lookup for one name
worms_lookup_one <- function(x) {
  if (is.na(x) || !nzchar(x)) return(tibble(AphiaID = NA_integer_, worms_status = NA_character_, worms_valid_name = NA_character_))

  res <- tryCatch(
    worrms::wm_records_name(name = x, fuzzy = TRUE, marine_only = TRUE),
    error = function(e) NULL
  )

  if (is.null(res) || nrow(res) == 0) {
    tibble(AphiaID = NA_integer_, worms_status = NA_character_, worms_valid_name = NA_character_)
  } else {
    tibble(
      AphiaID         = res$AphiaID[[1]],
      worms_status    = res$status[[1]],      # e.g., "accepted", "unaccepted"
      worms_valid_name= res$valid_name[[1]]   # accepted name WoRMS points to
    )
  }
}

##Clean up columns

worms_tbl_AIS <- AIS %>%
  distinct(Scientific.Name) %>%
  mutate(worms = map(Scientific.Name, worms_lookup_one)) %>%
  tidyr::unnest(worms)

AIS <- AIS %>%
  left_join(worms_tbl_AIS, by = "Scientific.Name") %>%
  mutate(
    worms_match = !is.na(AphiaID)
  )


sample_tag <- function(occ_all) {
  occ_all %>%
    dplyr::mutate(
      yr = dplyr::if_else(is.na(year), "", as.character(year)),
      mk = dplyr::if_else(is.na(target_gene), "", as.character(target_gene)),
      #sf = dplyr::if_else(is.na(source_file), "", as.character(source_file)),
      tag = paste0(samp_name,
                   dplyr::if_else(yr != "" | mk != "" | sf != "",
                                  paste0(" (",
                                         paste(c(yr, mk, sf)[c(yr, mk, sf) != ""], collapse = ", "),
                                         ")"),
                                  ""))
    ) %>%
    dplyr::pull(tag) %>%
    unique() %>%
    sort()
}


##Fix grid heatmap by year for all target_genes

grid_ALL_by_year <- RICHNESS_ALL_BY_YEAR   # already in 4326
grid_ALL_static  <- RICHNESS_ALL           # all years combined, already in 4326


# ----------------------------
# Depth colour mapping (3-zone ramps)
# ----------------------------

cool <- c(
  "#ccebc5", "#a8ddb5", "#7fcdbb",
  "#41b6c4", "#1d91c0", "#225ea8", "#0c2c84"
)

depth_breaks <- c(0, 25, 50, 75, 100, 500, 1000, Inf)
depth_labels <- c(
  "0–25 m", "26–50 m", "51–75 m", "76–100 m",
  "101–500 m", "501–1000 m", "1000+ m"
)

orange3 <- c("#fed976", "#feb24c", "#fd8d3c")                   # mixed depth overlay

orange_pal    <- grDevices::colorRampPalette(orange3)


get_depth_num <- function(occ_all) {
  depth_col <- dplyr::case_when(
    "minimumDepthInMeters" %in% names(occ_all) ~ "minimumDepthInMeters",
    "Depth" %in% names(occ_all) ~ "Depth",
    TRUE ~ NA_character_
  )
  if (is.na(depth_col)) return(occ_all %>% mutate(depth_m = NA_real_))

  occ_all %>% mutate(depth_m = suppressWarnings(as.numeric(.data[[depth_col]])))
}


# map depth (m) -> binned
depth_to_col <- function(d) {
  out <- rep(NA_character_, length(d))
  ok  <- !is.na(d)

  # findInterval returns 1..7 for these breaks
  idx <- findInterval(d[ok], vec = depth_breaks, rightmost.closed = TRUE)

  # idx could be 0 if d < 0; clamp just in case
  idx <- pmax(1, pmin(length(cool), idx))

  out[ok] <- cool[idx]
  out
}


# mixed depth -> orange shade (more orange with bigger within-cell range)
mixed_to_orange <- function(depth_range, cap_range = 200) {
  out <- rep(NA_character_, length(depth_range))
  i <- !is.na(depth_range)
  if (!any(i)) return(out)

  n <- 80
  rcap <- pmin(depth_range[i], cap_range)
  idx <- floor(rcap / cap_range * (n - 1)) + 1
  idx <- pmax(1, pmin(n, idx))
  out[i] <- orange_pal(n)[idx]
  out
}

DATA_BY_YEAR_GENE <- DATA_BY_KEY


# --- sample points by dataset (ONE per materialSampleID) ---
SAMPLE_PTS_BY_KEY <- purrr::imap(DATA_BY_YEAR_GENE, ~{
  parts <- strsplit(.y, "_")[[1]]
  gene <- parts[1]
  year <- parts[2]

  .x %>%
    get_depth_num() %>%
    filter(!is.na(decimalLatitude), !is.na(decimalLongitude)) %>%
    distinct(samp_name, decimalLatitude, decimalLongitude, depth_m, .keep_all = TRUE) %>%
    st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) %>%
    mutate(target_gene = gene, year = year)
})

sample_pts_all <- bind_rows(SAMPLE_PTS_BY_KEY) %>%
  mutate(year = as.character(year))

#Compute per-cell depth stats(min,max,median + mixed depth flag)
depth_on_grid_stats <- function(grid_sf, sample_pts_sf, depth_col = "depth_m",
                                fun_med = stats::median,
                                mixed_delta_m = 1.0) {  # <- threshold: min/max must differ by >= 1 m to count as "mixed"
  stopifnot("cell_id" %in% names(grid_sf))
  stopifnot(inherits(grid_sf, "sf"))
  stopifnot(inherits(sample_pts_sf, "sf"))

  if (sf::st_crs(sample_pts_sf) != sf::st_crs(grid_sf)) {
    sample_pts_sf <- sf::st_transform(sample_pts_sf, sf::st_crs(grid_sf))
  }
  if (sf::st_crs(all_polys) != sf::st_crs(grid_sf)) {
    all_polys2 <- sf::st_transform(all_polys, sf::st_crs(grid_sf))
  } else {
    all_polys2 <- all_polys
  }

  haspt <- lengths(sf::st_intersects(grid_sf, sample_pts_sf)) > 0
  grid_with <- sf::st_join(grid_sf, sample_pts_sf, join = sf::st_intersects, left = TRUE)

  stats_tbl <- grid_with %>%
    sf::st_drop_geometry() %>%
    dplyr::group_by(cell_id) %>%
    dplyr::summarise(
      depth_min = if (all(is.na(.data[[depth_col]]))) NA_real_ else min(.data[[depth_col]], na.rm = TRUE),
      depth_max = if (all(is.na(.data[[depth_col]]))) NA_real_ else max(.data[[depth_col]], na.rm = TRUE),
      depth_med = if (all(is.na(.data[[depth_col]]))) NA_real_ else fun_med(.data[[depth_col]], na.rm = TRUE),
      .groups   = "drop"
    )

  out <- grid_sf %>%
    dplyr::left_join(stats_tbl, by = "cell_id") %>%
    dplyr::mutate(
      has_sampling = haspt,
      depth_range  = depth_max - depth_min,
      mixed_depth  = dplyr::if_else(has_sampling & !is.na(depth_range) & depth_range >= mixed_delta_m, TRUE, FALSE)
    ) %>%
    sf::st_as_sf()

  sf::st_intersection(out, all_polys2) %>%
    sf::st_collection_extract("POLYGON") %>%
    (\(x) x[!sf::st_is_empty(x), ])()
}

#Compute depth grids by year only
# years available in sample points

# use the clipped grid that definitely has cell_id
YEARS_DEPTH <- sort(unique(na.omit(as.character(sample_pts_all$year))))

grid_for_depth <- grid_clip

depth_layers <- c(
  list("All" = depth_on_grid_stats(grid_for_depth, sample_pts_all)),
  setNames(
    purrr::map(YEARS_DEPTH, \(yr) depth_on_grid_stats(grid_for_depth, dplyr::filter(sample_pts_all, year == yr))),
    YEARS_DEPTH
  )
)

depth_layers <- lapply(depth_layers, function(occ_all) {
  occ_all %>%
    dplyr::mutate(
      mixed_depth = as.logical(mixed_depth),
      depth_fill  = depth_to_col(depth_med),
      mixed_fill  = mixed_to_orange(depth_range),
      final_fill  = dplyr::if_else(mixed_depth, mixed_fill, depth_fill, missing = depth_fill),
      final_fill  = dplyr::if_else(!has_sampling | is.na(depth_med), "#7bccc4", final_fill),
      final_alpha = dplyr::if_else(!has_sampling | is.na(depth_med), 0.05, 1)
    )
})


depth_legend_cols <- cool
depth_legend_labs <- depth_labels

#For later
#choices_gene <- sort(unique(KEY_TBL$target_gene))
#choices_year <- sort(unique(KEY_TBL$year))
#choices_key  <- KEY_TBL$key

#layer_sf <- RICHNESS_BY_KEY[[input$key]]

species_sf_min <- species_sf_all %>%
  dplyr::select(
    scientificName, target_gene, year,
    kingdom, phylum, class, order,
    geometry
  )

species_sf_by_year <- split(species_sf_min, as.character(species_sf_min$year))
species_sf_by_year[["All"]] <- species_sf_min

standardize_month_col <- function(df) {
  df %>%
    dplyr::mutate(
      month = dplyr::case_when(
        !is.na(suppressWarnings(as.integer(as.character(month)))) &
          suppressWarnings(as.integer(as.character(month))) %in% 1:12 ~
          month.abb[suppressWarnings(as.integer(as.character(month)))],
        as.character(month) %in% month.name ~ substr(as.character(month), 1, 3),
        as.character(month) %in% month.abb ~ as.character(month),
        TRUE ~ as.character(month)
      )
    )
}

occ_all        <- standardize_month_col(occ_all)
species_sf_all <- standardize_month_col(species_sf_all)

#Bundle the outputs in one list:
APP_DATA <- list(
  occ_all = occ_all,
  KEY_TBL = KEY_TBL,
  sampling_pts = sampling_pts,
  species_sf_all = species_sf_all,
  species_sf_min = species_sf_min,
  species_sf_by_year = species_sf_by_year,
  grid_clip = grid_clip,
  RICHNESS_BY_KEY = RICHNESS_BY_KEY,
  RICHNESS_ALL = RICHNESS_ALL,
  RICHNESS_ALL_BY_YEAR = RICHNESS_ALL_BY_YEAR,
  depth_layers = depth_layers,
  all_polys_click = all_polys_click,
  all_polys_zones = all_polys_zones,
  pal_rich = pal_rich
)

############################################################################################################
#Server + UI attempt using above code


## ---- 2) UI ----

ui <- fluidPage(
  shinyjs::useShinyjs(),
    theme = bslib::bs_theme(version = 3),

    tags$head(

      # jQuery UI (needed for draggable/resizable)
      tags$link(
        rel  = "stylesheet",
        href = "https://code.jquery.com/ui/1.13.2/themes/base/jquery-ui.css"
      ),
      tags$script(src = "https://code.jquery.com/ui/1.13.2/jquery-ui.min.js"),

      tags$style(HTML("
      body { padding-top: 62px; }
      .navbar { margin-bottom: 10px; }

      .navbar{
      background-color: #2241a7 !important;
      border-color: #2241a7 !important;
      z-index: 20000 !important;
      position: fixed !important;
      width: 100%;
      }

/* navbar text */
.navbar .navbar-brand,
.navbar .navbar-nav > li > a{
  color: #ffffff !important;
}

/* hover text */
.navbar .navbar-nav > li > a:hover,
.navbar .navbar-brand:hover{
  color: #ffffff !important;
  background-color: rgba(255,255,255,0.08) !important;
}

/* active tab */
.navbar .navbar-nav > .active > a,
.navbar .navbar-nav > .active > a:focus,
.navbar .navbar-nav > .active > a:hover{
  color: #ffffff !important;
  background-color: rgba(0,0,0,0.15) !important;
}

      /* Map container MUST be the positioning context + CLIP panel to map */
      #map_wrap{
        position: relative;
        height: calc(95vh - 62px);
        overflow: hidden; /* prevents panel overlapping other sections */
      }
      #map{ height: 100% !important; }

      /* Make the overlays panel wide enough */
      .leaflet-control-layers-overlays{ min-width: 220px; }

      /* Leaflet layer toggles: left aligned */
      .leaflet-control-layers-overlays{
        min-width: 220px;
        text-align: left !important;
      }
      .leaflet-control-layers-overlays label{
        display: flex !important;
        align-items: center !important;
        justify-content: flex-start !important;
        gap: 10px;
        width: 100%;
        margin: 6px 0;
      }
      .leaflet-control-layers-overlays label span{
        flex: 1 1 auto;
        text-align: left !important;
        white-space: normal;
        overflow-wrap: anywhere;
        line-height: 1.2;
      }

      /* Switch itself */
      .leaflet-control-layers-overlays input[type='checkbox']{
        flex: 0 0 auto;
        margin: 0 !important;

        appearance:none;
        -webkit-appearance:none;
        width:28px;
        height:16px;
        border-radius:999px;
        background:#cfd6dd;
        position:relative;
        cursor:pointer;
        outline:none;
        transition:background 0.15s ease-in-out;
      }
      .leaflet-control-layers-overlays input[type='checkbox']::after{
        content:'';
        position:absolute;
        top:2px; left:2px;
        width:12px; height:12px;
        border-radius:50%;
        background:white;
        box-shadow:0 1px 3px rgba(0,0,0,0.25);
        transition:transform 0.15s ease-in-out;
      }
      .leaflet-control-layers-overlays input[type='checkbox']:checked{ background:#2c7be5; }
      .leaflet-control-layers-overlays input[type='checkbox']:checked::after{
        transform:translateX(12px);
      }

     /* ---- Floating panel ---- */
#floating_panel{
  position: absolute;
  background: rgba(255,255,255,0.92);
  z-index: 9999;

  display: flex;
  flex-direction: column;

  /* keep panel inside map bounds */
  max-width:  calc(100% - 24px);
  max-height: calc(100% - 24px);

  /* panel itself should NOT scroll */
  overflow: hidden;

  /* optional nice defaults */
  border-radius: 6px;
  box-shadow: 0 2px 10px rgba(0,0,0,0.25);
}

/* Put Leaflet controls (esp. draw toolbar) above the floating panel */
#map_wrap .leaflet-top,
#map_wrap .leaflet-bottom{
  z-index: 12000 !important;
}

/* Specifically ensure the draw toolbar is above */
#map_wrap .leaflet-draw,
#map_wrap .leaflet-draw-toolbar,
#map_wrap .leaflet-control-draw{
  z-index: 13000 !important;
  position: relative; /* makes z-index apply reliably */
}

/* Let bootstrap control collapse visibility; we only improve layout */
#floating_body.collapse.in,
#floating_body.collapsing{
  display: flex !important;
  flex: 1 1 auto;
  flex-direction: column;

  /* REQUIRED for nested scrolling in flex layouts */
  min-height: 0;
}

/* Actual scrolling container */
#floating_body .panel-body{
  flex: 1 1 auto;

  /* REQUIRED for nested scrolling in flex layouts */
  min-height: 0;

  /* scrollbar comes back here */
  overflow-y: auto;

  /* guarantees a scroll region even if flex sizing is finicky */
  max-height: calc(100% - 52px);
}

/* Remove collapse animation delay (optional) */
#floating_body.collapse,
#floating_body.collapsing{
  -webkit-transition: none !important;
  transition: none !important;
}
#floating_panel{ transition: none !important; }


      /* resizable handles */
      #floating_panel .ui-resizable-handle{ z-index: 10000; }
      #floating_panel .ui-resizable-se{
        width: 14px; height: 14px;
        right: 2px; bottom: 2px;
        cursor: se-resize;
        background: rgba(0,0,0,0.20);
        border-radius: 3px;
      }
      #floating_panel .ui-resizable-e{
        width: 10px; right: 0px;
        top: 0; bottom: 0;
        cursor: e-resize;
      }
      #floating_panel .ui-resizable-s{
        height: 10px; bottom: 0px;
        left: 0; right: 0;
        cursor: s-resize;
      }

      /* Button grids */
      .filter-btn-grid-4{
        display: grid;
        grid-template-columns: repeat(4, 1fr);
        gap: 8px;
      }
      .filter-btn-grid-2{
        display: grid;
        grid-template-columns: repeat(2, 1fr);
        gap: 8px;
      }
      @media (max-width: 500px){
        .filter-btn-grid-4{ grid-template-columns: repeat(2, 1fr); }
        .filter-btn-grid-2{ grid-template-columns: repeat(2, 1fr); }
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
      .filter-btn-short{
        height: 40px !important;
        padding: 6px 10px !important;
        font-size: 15px !important;
      }

/* =========================
   GROUP FILTER BUTTONS
   ========================= */

/* Remove Bootstrap orange focus border */
.filter-btn:focus,
.filter-btn:active,
.filter-btn.active,
.filter-btn:focus:active{
  outline: none !important;
  box-shadow: none !important;
}

/* OFF button border */
.filter-btn{
  border-color: #d0d0d0 !important;
}

/* ON button border = darker version of fill */
.filter-btn.btn-group-on{
  background-color: #2241a7 !important;
  border-color: #2241a7 !important;  /* darker than fill */
  color: #0b2a2a !important;
}

/* OFF state (default look) */
.filter-btn{
  background-color: #f7f7f7 !important;
  border-color: #d0d0d0 !important;
  color: #0b2a2a !important;
}

/* ON state (group buttons) */
.filter-btn.btn-group-on{
  background-color: #2241a7 !important;
  border-color: #2241a7 !important;
  color: #FFFFFF !important;   /* <-- white text */
}

/* OFF hover */
.filter-btn:hover{
  background-color: #eeeeee !important;
}

/* ON hover (different colour) */
.filter-btn.btn-group-on:hover{
  background-color: #2241a7 !important;
  border-color: #2241a7 !important;
  color: #FFFFFF !important;
}

/* SARA / AIS buttons */
.filter-btn.btn-sara-on,
.filter-btn.btn-ais-on{
  background-color: #2241a7 !important;
  border-color: #2241a7 !important;
  color: #FFFFFF !important;   /* white text */
}

/* SARA / AIS hover when ON */
.filter-btn.btn-sara-on:hover,
.filter-btn.btn-ais-on:hover{
  background-color: #2241a7 !important;
  border-color: #2241a7 !important;
}

/* Download button (grey button inside blue panel) */
.btn-download-got,
.btn-download-got:focus,
.btn-download-got:active{
  background-color: #f2f2f2 !important;
  border-color: #d0d0d0 !important;
  color: #000000 !important;
  font-size: 16px !important;
  padding: 10px 22px !important;
  border-radius: 6px !important;
  box-shadow: none !important;
}

/* hover state */
.btn-download-got:hover{
  background-color: #e6e6e6 !important;
  border-color: #d0d0d0 !important;
  color: #000000 !important;
}

      /* Floating toggle header */
      #floating_toggle{
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

      /* Legends */
      .leaflet-control.legend-base{
        font-family: Arial !important;
        font-size: 15px !important;
        font-weight: 400 !important;
        line-height: 1.45 !important;
        color: #222 !important;
      }
      .leaflet-control.legend-base i{
        width: 18px !important;
        height: 18px !important;
        display: inline-block !important;
        margin-right: 8px !important;
        vertical-align: middle !important;
        float: none !important;
      }
      .leaflet-control.legend-richness-box,
      .leaflet-control.legend-depth-box{
        background: rgba(255,255,255,0.92) !important;
        border-radius: 6px !important;
        box-shadow: 0 2px 10px rgba(0,0,0,0.25) !important;
        border: 1px solid rgba(0,0,0,0.15) !important;
        padding: 14px 16px !important;
        width: 190px !important;
      }
      .legend-hidden{ display:none !important; }

/* ---- Data selection layout ---- */

.data-select-grid{
  display: grid;
  grid-template-columns: 250px 300px 300px 110px;
  column-gap: 50px;
  row-gap: 0;
  align-items: start;
}

.data-select-item{
  min-width: 0;
}

.data-select-item .form-group{
  margin-bottom: 0 !important;
}

.data-select-grid .selectize-control{
  margin-bottom: 0 !important;
}

.primer-btn-row{
  margin-top: 6px;
  display: flex;
  gap: 6px;
}

/* right-side action buttons */
.confirm-slot{
  align-self: start;
  margin-top: 24px;
}

.confirm-btn-row{
  display: flex;
  align-items: stretch;
  gap: 8px;
  width: 100%;
}

.confirm-btn-row .btn{
  height: 38px !important;   /* same height for both */
  display: inline-flex;
  align-items: center;
  justify-content: center;
  margin: 0 !important;
}

.confirm-btn-row .btn-primary{
  flex: 0 0 110px;           /* Confirm width */
}

.confirm-btn-row .btn-download-got{
  flex: 0 0 160px;           /* Download width */
  padding: 6px 12px !important;
  font-size: 14px !important;
}

@media (max-width: 1200px){
  .data-select-grid{
    grid-template-columns: 1fr 1fr;
    row-gap: 12px;
  }

  .confirm-slot{
    margin-top: 0;
  }

  .confirm-btn-row{
    width: 100%;
    flex-wrap: wrap;
  }
}

/* optional: reduce spacing under labels a bit */
.data-select-grid label{
  margin-bottom: 4px !important;
}

@media (max-width: 1200px){
  .data-select-grid{
    grid-template-columns: 1fr 1fr;
    row-gap: 12px;
  }

  .confirm-slot{
    margin-top: 0;
  }

  .confirm-slot .btn{
    width: auto;
  }
}

table.dataTable.nowrap td,
table.dataTable.nowrap th {
  white-space: nowrap;
}

/* circular monthly plot panel */
#monthly_plot_control{
  position: absolute;
  z-index: 11000;
  width: 350px;
  background: rgba(255,255,255,0.92);
  border-radius: 6px;
  box-shadow: 0 2px 10px rgba(0,0,0,0.25);
  border: 1px solid rgba(0,0,0,0.15);
  padding: 8px 8px 2px 8px;
  display: none;
}

#monthly_plot_title{
  font-size: 15px;
  font-weight: 600;
  margin-bottom: 4px;
  line-height: 1.2;
}

#monthly_plot_subtitle{
  font-size: 12px;
  color: #444;
  margin-bottom: 2px;
  line-height: 1.2;
}

/* remove default leaflet control margins when we manually position items */
#monthly_plot_control{
  margin: 0 !important;
}

.leaflet-control.legend-richness-box,
.leaflet-control.legend-depth-box{
  margin: 0 !important;
}

/* remove leaflet's default control spacing */
#monthly_plot_control,
.leaflet-control.legend-richness-box,
.leaflet-control.legend-depth-box{
  margin: 0 !important;
}

/* when we manually reposition legends, absolute positioning behaves more predictably */
.leaflet-right .leaflet-control.legend-richness-box,
.leaflet-right .leaflet-control.legend-depth-box{
  position: absolute !important;
}
    ")),

      tags$script(HTML("
$(function(){

  // Smooth scroll navbar
  $(document).on('click', 'a.nav-scroll', function(e){
    e.preventDefault();
    var target = $(this).data('target');
    var el = document.getElementById(target);
    if(el){
      el.scrollIntoView({behavior:'smooth', block:'start'});
    }
  });

    // Ensure Bootstrap 3 collapse is initialized on our target (prevents silent no-op)
  function ensureCollapseInit(){
    var $body = $('#floating_body');
    if(!$body.length) return;

    if (typeof $body.collapse === 'function') {
      // Initialize plugin without toggling
      $body.collapse({ toggle: false });
    }
  }

  function clampFloatingToMap(){
    var $p = $('#floating_panel');
    var $wrap = $('#map_wrap');
    if(!$p.length || !$wrap.length) return;

    var PAD = 14;

    var wrapW = $wrap.innerWidth();
    var wrapH = $wrap.innerHeight();

    $p.css({
      'max-width':  (wrapW - PAD*2) + 'px',
      'max-height': (wrapH - PAD*2) + 'px'
    });

    var pos = $p.position();
    var pW  = $p.outerWidth();
    var pH  = $p.outerHeight();

    var left = pos.left;
    var top  = pos.top;

    var minL = PAD;
    var minT = PAD;
    var maxL = wrapW - pW - PAD;
    var maxT = wrapH - pH - PAD;

    if(maxL < minL) maxL = minL;
    if(maxT < minT) maxT = minT;

    left = Math.min(Math.max(left, minL), maxL);
    top  = Math.min(Math.max(top,  minT), maxT);

    $p.css({ left: left + 'px', top: top + 'px' });
  }

  function setFloatingCollapsedUI(isCollapsed){
    var $p = $('#floating_panel');
    var $toggle = $('#floating_toggle');
    if(!$p.length || !$toggle.length) return;

    if(isCollapsed){
      if($('#floating_body').hasClass('in')){
        $p.data('open_h', $p.outerHeight());
      }

      $p.addClass('is-collapsed');

      var headerH = $toggle.outerHeight(true) + 2;
      if($p.hasClass('ui-resizable')){
        $p.resizable('option', 'minHeight', headerH);
      }
      $p.css({ height: headerH + 'px' });

    } else {
    $p.removeClass('is-collapsed');

    if($p.hasClass('ui-resizable')){
      $p.resizable('option', 'minHeight', 180);
    }

    // If we've never been opened before, pick a sane default height
    var oh = $p.data('open_h');
    if(!oh){
      var $wrap = $('#map_wrap');
      var wrapH = $wrap.length ? $wrap.innerHeight() : 700;
      oh = wrapH * 0.99;
      $p.data('open_h', oh);
    }
    $p.css({ height: oh + 'px' });
  }

    setTimeout(clampFloatingToMap, 0);
  }

  function enableFloatingResize(){
    var $p = $('#floating_panel');
    if(!$p.length) return;

    // IMPORTANT: we control draggable here. So set absolutePanel(draggable=FALSE).
    if($p.hasClass('ui-draggable')){
      $p.draggable('option', 'containment', '#map_wrap');
      $p.draggable('option', 'handle', '#floating_toggle');
    } else {
      $p.draggable({ containment: '#map_wrap', handle: '#floating_toggle' });
    }

    if(!$p.hasClass('ui-resizable')){
      $p.resizable({
        handles: 'e,s,se',
        minWidth: 320,
        minHeight: 200,
        containment: '#map_wrap'
      });
    } else {
      $p.resizable('option', 'containment', '#map_wrap');
    }

    $p.off('dragstop.clamp').on('dragstop.clamp', clampFloatingToMap);
    $p.off('resizestop.clamp').on('resizestop.clamp', clampFloatingToMap);

    $(window).off('scroll.clampFloating').on('scroll.clampFloating', clampFloatingToMap);
    $(window).off('resize.clampFloating').on('resize.clampFloating', clampFloatingToMap);

    setTimeout(clampFloatingToMap, 0);
    setTimeout(clampFloatingToMap, 250);
  }

  function initFloatingPanel(){
    ensureCollapseInit();
    enableFloatingResize();
    setFloatingCollapsedUI(true); // start collapsed
  }

  $(document).on('shiny:connected', initFloatingPanel);
  initFloatingPanel();

  // Bootstrap 3 collapse events
  $(document).on('hide.bs.collapse', '#floating_body', function(){
    setFloatingCollapsedUI(true);
  });

  $(document).on('hidden.bs.collapse', '#floating_body', function(){
    enableFloatingResize();
    setFloatingCollapsedUI(true);
    clampFloatingToMap();
  });

  $(document).on('shown.bs.collapse', '#floating_body', function(){
    enableFloatingResize();
    setFloatingCollapsedUI(false);
    clampFloatingToMap();
  });

  // Custom message to open the panel
  Shiny.addCustomMessageHandler('openFloating', function(message){
    var $p    = $('#floating_panel');
    var $body = $('#floating_body');
    var $tog  = $('#floating_toggle');
    if(!$p.length || !$body.length || !$tog.length) return;

     ensureCollapseInit();  // <<< add this too (safe)

    if($body.hasClass('in')){
      $tog.attr('aria-expanded','true');
      setFloatingCollapsedUI(false);
      enableFloatingResize();
      clampFloatingToMap();
      return;
    }

    $body.one('shown.bs.collapse.openFloating', function(){
      $tog.attr('aria-expanded','true');
      requestAnimationFrame(function(){
        enableFloatingResize();
        setFloatingCollapsedUI(false);
        clampFloatingToMap();
      });
    });

    if(typeof $body.collapse === 'function'){
      $body.collapse('show');
    } else {
      $body.addClass('in').css('display','block');
      $body.trigger('shown.bs.collapse');
    }
  });
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
        tags$li(tags$a(class="nav-scroll", href="#", `data-target`="sec_sara",   "Detection Details")),
        #tags$li(tags$a(class="nav-scroll", href="#", `data-target`="sec_method", "Method Comparison")),
        tags$li(tags$a(class="nav-scroll", href="#", `data-target`="sec_datsel", "Data Selection and Download")),
        tags$li(tags$a(class="nav-scroll", href="#", `data-target`="sec_div",    "Diversity Metrics")),
        tags$li(tags$a(class="nav-scroll", href="#", `data-target`="sec_pie",    "Taxonomic Pie Chart")),
        #tags$li(tags$a(class="nav-scroll", href="#", `data-target`="sec_dwnld",    "Download Data File")),
        #tags$li(tags$a(class="nav-scroll", href="#", `data-target`="sec_refdat",    "Reference Data Authorship")),
      )
    )
  ),

  # ---- MAP SECTION ----
  div(
    id = "sec_map", class = "scroll-section",
    div(
      id = "map_wrap",
      leafletOutput("map"),

      div(
        id = "monthly_plot_control",
        class = "leaflet-control",
        div(id = "monthly_plot_title", "Monthly samples collected"),
        div(id = "monthly_plot_subtitle", textOutput("monthly_plot_subtitle", inline = TRUE)),
        plotOutput("monthly_circular_plot", height = "320px", width = "340px")
      ),

      absolutePanel(
        id = "floating_panel",
        fixed = FALSE, draggable = FALSE,
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
            tags$span("Select a Site"),
            tags$span(class = "caret-icon", HTML("&#9662;"))
          )
        ),

        div(
          id = "floating_body",
          class = "collapse",
          div(
            class = "panel-body",
            h4("Filter"),
            selectInput(
                "sel_year",
                "Year",
                choices = c("All"), #use this code if you don't want to be able to select more than one year at a time
                selected = "All"),

                #choices  = "All",  #use this code to be able to select more than one year at a time
                #selected = "All",
                #multiple = TRUE
                #),

            h4("Group"),

            # --- Row 1: 4 across ---
            div(
              class = "filter-btn-grid-4",
              actionButton("total_fish",     "Fishes",        class = "btn btn-default filter-btn"),
              actionButton("total_sharks",   "Sharks & Rays", class = "btn btn-default filter-btn"),
              actionButton("total_mammals",  "Mammals",       class = "btn btn-default filter-btn"),
              actionButton("total_reptiles", "Turtles",      class = "btn btn-default filter-btn")
            ),

            tags$div(style="height:8px;"),  # optional spacing between rows

            # --- Row 2: 4 across ---
            div(
              class = "filter-btn-grid-4",
              actionButton("total_birds",      "Birds",       class = "btn btn-default filter-btn"),
              actionButton("total_molluscs",   "Molluscs",    class = "btn btn-default filter-btn"),
              actionButton("total_arthropods", "Arthropods",  class = "btn btn-default filter-btn"),
              actionButton("total_plants",     "Plants",      class = "btn btn-default filter-btn")
            ),

            tags$div(style="height:8px;"),

            # --- Row 3: 2 across ---
            div(
              class = "filter-btn-grid-2",
              actionButton("SARA", "SARA", class = "btn btn-default filter-btn filter-btn-short"),
              actionButton("AIS",  "AIS",  class = "btn btn-default filter-btn filter-btn-short")
            ),

            hr(),
            h4("Species List"),
              uiOutput("species_panel")
          )
        )
      )
    )
  ),

  # ---- DETECTION TABLE SECTION ----
  div(
    id = "sec_sara", class = "scroll-section",
    tabsetPanel(
      tabPanel("Detection Details",              DT::DTOutput("detections_tbl")),
      tabPanel("Species At Risk Act (SARA): Schedule 1 Details", DT::DTOutput("sara_details")),
      tabPanel("Aquatic Invasive Species (AIS) Details",             DT::DTOutput("ais_details"))
    )
  ),

    # ---- DATA SELECTION SECTION ----
  div(
    id = "sec_datsel", class = "scroll-section",
    h3("Data Selection and Download"),

    div(
      class = "data-select-grid",

      div(
        class = "data-select-item",
        selectInput(
          "tax_rank",
          "Taxonomy selection",
          choices = c(
            "Kingdom" = "kingdom",
            "Phylum"  = "phylum",
            "Class"   = "class",
            "Order"   = "order",
            "Family"  = "family",
            "Genus"   = "genus",
            "Species" = "scientificName"
          ),
          selected = "scientificName"
        )
      ),

      div(
        class = "data-select-item",
        selectizeInput(
          "div_target_gene",
          "Target gene",
          choices = NULL,
          selected = NULL,
          multiple = TRUE,
          options = list(
            plugins = list("remove_button"),
            placeholder = "Select target gene(s)"
          )
        )
      ),

      div(
        class = "data-select-item",
        selectizeInput(
          "div_primer",
          "Primer",
          choices = NULL,
          selected = NULL,
          multiple = TRUE,
          options = list(
            plugins = list("remove_button"),
            placeholder = "Select primer(s)"
          )
        ),
        div(
          class = "primer-btn-row",
          actionButton("div_primer_all", "Select all", class = "btn btn-default btn-sm"),
          actionButton("div_primer_none", "Deselect all", class = "btn btn-default btn-sm")
        )
      ),
      div(
        class = "data-select-item confirm-slot",
        div(
          class = "confirm-btn-row",
          actionButton(
            "div_apply",
            "Confirm",
            class = "btn btn-primary"
          ),
          downloadButton(
            "downloadData",
            "Download Data",
            class = "btn-download-got"
          )
        )
      )
    )
  ),

    #   #tags$p(tags$strong("Alpha diversity (boxplot): "), "I want to include the following alpha metrics as a dropdown: Observed, Chao1, Shannon, Simpson, Fisher, InvSimpson, ACE"),   #include if we want written discriptions as part of this section
    #   #tags$p(tags$strong("Beta diversity (PCoA): "), "I want to include the following beta metrics as a dropdown: Bray, Jaccard, Euclidean, Aitchison"),

  # ---- DIVERSITY METRICS SECTION ----
  div(
    id = "sec_div", class = "scroll-section",
    h3("Diversity Metrics"),

    div(
      class = "data-select-grid",

      div(
        class = "data-select-item",
        selectInput(
          "alpha_metric",
          "Alpha Diversity",
          choices = c(
            "Observed (Richness)" = "observed",
            "Shannon"             = "shannon",
            "Simpson"             = "simpson",
            "InvSimpson"          = "invsimpson",
            "ACE"                 = "ace"
          ),
          selected = "observed"
        )
      ),

      div(class = "data-select-item"),
      div(class = "data-select-item"),
      div(class = "data-select-item")
    )
  ),

      #tags$p(tags$strong("Alpha diversity (boxplot): "), "I want to include the following alpha metrics as a dropdown: Observed, Chao1, Shannon, Simpson, Fisher, InvSimpson, ACE"),   #include if we want written discriptions as part of this section
      #tags$p(tags$strong("Beta diversity (PCoA): "), "I want to include the following beta metrics as a dropdown: Bray, Jaccard, Euclidean, Aitchison"),

    # ---- Alpha plot (full width row) ----
  fluidRow(
    column(
      width = 8,
      offset = 2,
      shinycssloaders::withSpinner(
        plotly::plotlyOutput("alpha_boxplot", height = "700px"),
        type = 4
      )
    )
  ),

    # ---- Beta plot (full width row) ----
  div(
    class = "data-select-grid",

    div(
      class = "data-select-item",
      selectInput(
        "beta_metric",
        "Beta Diversity",
        choices = c(
          "Bray-Curtis" = "bray",
          "Jaccard"     = "jaccard",
          "Euclidean"   = "euclidean",
          "Robust Aitchison"  = "robust.aitchison"
        ),
        selected = "bray"
      )
    ),

    div(class = "data-select-item"),
    div(class = "data-select-item"),
    div(class = "data-select-item")
  ),

  # ---- Beta plot ----
  fluidRow(
    column(
      width = 9,
      offset = 2,
      shinycssloaders::withSpinner(
        plotly::plotlyOutput("beta_pcoa", height = "700px"),
        type = 4
      )
    )
  ),

  # ---- Taxonomic Pie Chart ----
  div(
    id = "sec_pie", class = "scroll-section",
    h3("Taxonomic Pie Chart"),

    shinycssloaders::withSpinner(
      taxplore::KronaChartOutput("tax_krona", height = "700px"),
      type = 4
    )
  )
)

#   # ---- Download Data Option ----
#   div(
#     id = "sec_dwnld", class = "scroll-section",
#     style = "
#     min-height: 820px;
#     width: 100%;
#     max-width: 1000px;
#     margin: 0 auto;
#     padding: 40px 20px;
#     display: flex;
#     align-items: center;
#     justify-content: center;
#   ",
#
#     div(
#       style = "
#       width: 100%;
#       max-width: 700px;
#       min-height: 260px;
#       text-align: center;
#       padding: 40px 30px;
#       background: #2241a7;
#       color: white;
#       border-radius: 8px;
#       box-shadow: 0 2px 8px rgba(0,0,0,0.08);
#       display: flex;
#       flex-direction: column;
#       justify-content: center;
#     ",
#       h3("Download Data File"),
#       p("Download a CSV filtered by the current floating-panel and diversity-metrics selections."),
#       br(),
#       downloadButton("downloadData", "Download Filtered Data", class = "btn-download-got")
#     )
#   )
# )
#
#   # ---- Reference Data Authorship ----
#
#   div(
#     id="sec_refdat", class="scroll-section",
#     h3("Reference Data Authorship"),
#
#   )
# )


## ---- 3) Server ----

server <- function(input, output, session){

  `%||%` <- function(x, y) if (is.null(x)) y else x

  species_cache <- reactiveValues(data = list())

  make_species_cache_key <- function(selected_key, year, groups, sara_on, ais_on) {
    paste(
      selected_key,
      year,
      paste(sort(groups), collapse = "|"),
      paste0("sara=", sara_on),
      paste0("ais=", ais_on),
      sep = "~~"
    )
  }

  get_species_by_layer_cached <- function(selected_key, det_sf, year, groups, sara_on, ais_on) {

    cache_key <- make_species_cache_key(
      selected_key = selected_key,
      year         = year,
      groups       = groups,
      sara_on      = sara_on,
      ais_on       = ais_on
    )

    if (!is.null(species_cache$data[[cache_key]])) {
      return(species_cache$data[[cache_key]])
    }

    det_f <- apply_species_filters(det_sf) %>%
      sf::st_drop_geometry() %>%
      dplyr::mutate(
        scientificName = as.character(scientificName),
        target_gene    = as.character(target_gene)
      )

    out <- list(
      All = det_f %>%
        dplyr::pull(scientificName) %>%
        unique() %>%
        stats::na.omit() %>%
        sort()
    )

    for (g in c("12S", "COI", "16S", "18S")) {
      out[[g]] <- det_f %>%
        dplyr::filter(target_gene == g) %>%
        dplyr::pull(scientificName) %>%
        unique() %>%
        stats::na.omit() %>%
        sort()
    }

    species_cache$data[[cache_key]] <- out
    out
  }

  inside_cache <- reactiveValues(data = list())

  make_inside_cache_key <- function(selected_key, year) {
    paste(selected_key, year, sep = "~~")
  }

  get_inside_cached <- function(selected_key, pts, geom, year) {
    cache_key <- make_inside_cache_key(selected_key, year)

    if (!is.null(inside_cache$data[[cache_key]])) {
      return(inside_cache$data[[cache_key]])
    }

    inside <- pts[within_any(pts, geom), , drop = FALSE]
    inside_cache$data[[cache_key]] <- inside
    inside
  }

  get_forward_primer_col <- function(df) {
    cand <- c("pcr_primer_name_forward", "pcr_primer_forward")
    hit <- intersect(cand, names(df))
    if (length(hit) == 0) return(NULL)
    hit[1]
  }

  get_reverse_primer_col <- function(df) {
    cand <- c("pcr_primer_name_reverse", "pcr_primer_reverse")
    hit <- intersect(cand, names(df))
    if (length(hit) == 0) return(NULL)
    hit[1]
  }

  add_primer_combo <- function(df) {
    fwd_col <- get_forward_primer_col(df)
    rev_col <- get_reverse_primer_col(df)

    if (is.null(fwd_col) && is.null(rev_col)) {
      df$primer_combo <- NA_character_
      return(df)
    }

    fwd <- if (!is.null(fwd_col)) as.character(df[[fwd_col]]) else rep(NA_character_, nrow(df))
    rev <- if (!is.null(rev_col)) as.character(df[[rev_col]]) else rep(NA_character_, nrow(df))

    fwd <- trimws(fwd)
    rev <- trimws(rev)

    fwd[fwd == ""] <- NA_character_
    rev[rev == ""] <- NA_character_

    df$primer_combo <- dplyr::case_when(
      !is.na(fwd) & !is.na(rev) ~ paste(fwd, rev, sep = " | "),
      !is.na(fwd) &  is.na(rev) ~ fwd,
      is.na(fwd) & !is.na(rev) ~ rev,
      TRUE ~ NA_character_
    )

    df
  }

  apply_diversity_dropdown_filters <- function(df, filters) {
    df <- add_primer_combo(df)

    if (length(filters$target_gene) > 0) {
      df <- df %>%
        dplyr::filter(as.character(target_gene) %in% filters$target_gene)
    }

    if (length(filters$primers) > 0) {
      df <- df %>%
        dplyr::filter(primer_combo %in% filters$primers)
    }

    df
  }

  within_any <- function(x_sf, geom) {
    if (is.null(geom) || nrow(x_sf) == 0) {
      return(rep(FALSE, nrow(x_sf)))
    }

    if (inherits(geom, "sfc")) {
      geom <- sf::st_sf(geometry = geom)
    }

    if (sf::st_crs(x_sf) != sf::st_crs(geom)) {
      geom <- sf::st_transform(geom, sf::st_crs(x_sf))
    }

    lengths(sf::st_intersects(x_sf, geom)) > 0
  }

  # store ALL drawn polygons (one row per polygon), with a stable id
  empty_drawn_sf <- sf::st_sf(
    draw_id    = character(0),
    draw_label = character(0),
    geometry   = sf::st_sfc(crs = 4326)
  )

  drawn_polys <- reactiveVal(empty_drawn_sf)
  selected_draw_id <- reactiveVal(NULL)

  detections_filtered <- reactive({
    det <- selected_detections()
    if (is.null(det) || nrow(det) == 0) return(NULL)

    det %>%
      apply_species_filters()
  })

  # ---- helper: convert leaflet.draw feature -> sf polygon (EPSG:4326) ----
  feature_to_sf <- function(feature) {
    req(feature$geometry$type)
    type <- feature$geometry$type
    coords <- feature$geometry$coordinates

    close_ring <- function(mat) {
      if (!all(mat[1, ] == mat[nrow(mat), ])) {
        mat <- rbind(mat, mat[1, , drop = FALSE])
      }
      mat
    }

    if (type == "Polygon") {
      rings <- lapply(coords, function(ring) {
        mat <- do.call(rbind, lapply(ring, function(x) c(x[[1]], x[[2]])))
        close_ring(mat)
      })

      poly <- sf::st_polygon(rings)
      out  <- sf::st_sf(geometry = sf::st_sfc(poly, crs = 4326))
      return(sf::st_make_valid(out))
    }

    if (type == "MultiPolygon") {
      mp <- lapply(coords, function(poly_i) {
        lapply(poly_i, function(ring) {
          mat <- do.call(rbind, lapply(ring, function(x) c(x[[1]], x[[2]])))
          close_ring(mat)
        })
      })

      geom <- sf::st_multipolygon(mp)
      out  <- sf::st_sf(geometry = sf::st_sfc(geom, crs = 4326))
      return(sf::st_make_valid(out))
    }

    stop("Drawn feature type not supported: ", type)
  }

  # ---- selection geometry (drawn polygon OR clicked polygon OR clicked grid cell) ----
  selection_geom <- reactive({
    sel_id <- selected_draw_id()
    if (!is.null(sel_id)) {
      polys <- drawn_polys()
      hit <- polys %>% dplyr::filter(draw_id == sel_id)
      if (nrow(hit) > 0) return(sf::st_geometry(hit))
    }

    click <- input$map_shape_click
    if (is.null(click) || is.null(click$id)) return(NULL)

    if (grepl("\\|\\|", click$id)) {
      parts <- strsplit(click$id, "\\|\\|")[[1]]
      p_type <- parts[1]
      p_name <- parts[2]
      poly_sel <- all_polys_click %>% dplyr::filter(site_type == p_type, site_name == p_name)
      if (nrow(poly_sel) == 0) return(NULL)
      return(sf::st_geometry(poly_sel))
    }

    cid <- suppressWarnings(as.integer(click$id))
    if (is.na(cid)) return(NULL)

    cell_poly <- grid_clip %>% dplyr::filter(cell_id == cid)
    if (nrow(cell_poly) == 0) return(NULL)

    sf::st_geometry(cell_poly)
  })

  observe({
    click <- input$map_shape_click
    proxy <- leafletProxy("map")

    proxy %>% clearGroup("Selection outline")

    if (is.null(click) || is.null(click$id)) return()

    id <- as.character(click$id)

    # only draw white outline for clicked grid cells
    cid <- suppressWarnings(as.integer(id))
    if (is.na(cid)) return()

    sel_sf <- grid_clip %>%
      dplyr::filter(cell_id == cid)

    if (nrow(sel_sf) == 0) return()

    proxy %>%
      addPolygons(
        data        = sel_sf,
        group       = "Selection outline",
        layerId     = ~paste0("selected_cell_", cell_id),
        fill        = FALSE,
        color       = "white",
        weight      = 3,
        opacity     = 1,
        options     = pathOptions(
          pane = "pane_selected_top",
          interactive = FALSE
        )
      )
  })

  # ---- Year selection (as character or "All") ----
  sel_year_chr <- reactive({
    yr <- input$sel_year %||% "All"
    as.character(yr)
  })

  species_gene_summary <- function(det_sf) {
    det_sf %>%
      sf::st_drop_geometry() %>%
      dplyr::mutate(
        scientificName = as.character(scientificName),
        target_gene    = as.character(target_gene),
        samp_name      = as.character(samp_name)
      ) %>%
      dplyr::filter(!is.na(scientificName), scientificName != "") %>%
      dplyr::group_by(scientificName) %>%
      dplyr::summarise(
        genes       = paste(sort(unique(na.omit(target_gene))), collapse = ", "),
        n_detections = dplyr::n(),
        n_samples    = dplyr::n_distinct(na.omit(samp_name)),
        .groups = "drop"
      ) %>%
      dplyr::arrange(scientificName)
  }

  # --- Is the Sampling points layer currently visible? ---
  sampling_points_layer_on <- reactive({
    groups_on <- input$map_groups %||% character(0)
    "Sampling points" %in% groups_on
  })

  # ---- meta table (ensures join keys are character) ----
  meta_all <- reactive({
    req(occ_all)
    occ_all %>%
      dplyr::mutate(
        occurrenceID = as.character(occurrenceID)
      )
  })

  # Which richness layers are currently ON (including "All")
  active_richness_layers <- reactive({
    groups_on <- input$map_groups %||% character(0)
    intersect(groups_on, c("All","12S","COI","16S","18S"))
  })

  make_list <- function(vec, max_h = 220) {
    if (length(vec) == 0) return(em("No species detected for this layer."))
    tags$div(
      style = paste0("max-height:", max_h, "px; overflow-y:auto; padding-left: 10px;"),
      tags$ul(lapply(vec, tags$li))
    )
  }

  # Return a named list: each name is a layer ("12S", "COI", ... or "All"),
  # each value is the species vector for that layer
  species_by_active_layers <- function(det_sf, layers_on, apply_filters_fn) {
    det_sf <- apply_filters_fn(det_sf)

    out <- list()

    # "All" means no gene filter
    if ("All" %in% layers_on) {
      spp_all <- det_sf %>%
        sf::st_drop_geometry() %>%
        dplyr::pull(scientificName) %>%
        as.character() %>%
        unique() %>%
        stats::na.omit() %>%
        sort()

      out[["All"]] <- spp_all
    }

    # gene-specific lists
    genes <- setdiff(layers_on, "All")
    for (g in genes) {
      spp_g <- det_sf %>%
        dplyr::filter(as.character(target_gene) == g) %>%
        sf::st_drop_geometry() %>%
        dplyr::pull(scientificName) %>%
        as.character() %>%
        unique() %>%
        stats::na.omit() %>%
        sort()

      out[[g]] <- spp_g
    }

    out
  }

  #remove cache of deleted polygons
  observeEvent(input$map_draw_deleted_features, {
    del <- input$map_draw_deleted_features
    if (is.null(del) || is.null(del$features) || length(del$features) == 0) return()

    del_ids <- vapply(del$features, function(f) {
      if (!is.null(f$properties) && !is.null(f$properties$`_leaflet_id`)) {
        as.character(f$properties$`_leaflet_id`)
      } else {
        NA_character_
      }
    }, character(1))

    del_ids <- del_ids[!is.na(del_ids) & nzchar(del_ids)]
    if (length(del_ids) == 0) return()

    old_keys <- names(species_cache$data)
    if (length(old_keys) > 0) {
      keep_keys <- old_keys[!vapply(old_keys, function(k) {
        any(vapply(del_ids, function(id) grepl(paste0("draw:", id), k, fixed = TRUE), logical(1)))
      }, logical(1))]

      species_cache$data <- species_cache$data[keep_keys]
    }

    old_inside_keys <- names(inside_cache$data)
    if (length(old_inside_keys) > 0) {
      keep_inside_keys <- old_inside_keys[!vapply(old_inside_keys, function(k) {
        any(vapply(del_ids, function(id) grepl(paste0("draw:", id), k, fixed = TRUE), logical(1)))
      }, logical(1))]

      inside_cache$data <- inside_cache$data[keep_inside_keys]
    }

  })

  # ---- year dropdown: only show years present in current selection geom ----
  observeEvent(selection_geom(), {
    pts <- species_sf_min
    shiny::req(pts)
    shiny::req(nrow(pts) > 0)

    g <- selection_geom()

    sel_id <- selected_draw_id()
    click  <- input$map_shape_click

    selected_key <- NULL
    if (!is.null(sel_id)) {
      selected_key <- paste0("draw:", sel_id)
    } else if (!is.null(click$id)) {
      selected_key <- paste0("click:", click$id)
    }

    inside <- if (!is.null(g) && !is.null(selected_key)) {
      get_inside_cached(
        selected_key = selected_key,
        pts          = pts,
        geom         = g,
        year         = "AllYearsForDropdown"
      )
    } else {
      pts
    }

    yrs <- inside %>%
      sf::st_drop_geometry() %>%
      dplyr::pull(year) %>%
      as.character() %>%
      unique() %>%
      stats::na.omit() %>%
      sort()

    cur <- as.character(input$sel_year %||% "All")
    new_choices  <- c("All", yrs)
    new_selected <- if (cur %in% new_choices) cur else "All"

    updateSelectInput(session, "sel_year", choices = new_choices, selected = new_selected)
  }, ignoreInit = FALSE)

  observeEvent(input$map_draw_new_feature, {
    feat <- input$map_draw_new_feature
    req(feat)

    poly_sf <- tryCatch(
      feature_to_sf(feat),
      error = function(e) {
        showNotification(
          paste("Could not read drawn polygon:", e$message),
          type = "error"
        )
        return(NULL)
      }
    )

    req(!is.null(poly_sf), nrow(poly_sf) > 0)

    fid <- NULL
    if (!is.null(feat$properties) && !is.null(feat$properties$`_leaflet_id`)) {
      fid <- as.character(feat$properties$`_leaflet_id`)
    }
    if (is.null(fid) || !nzchar(fid)) {
      fid <- paste0("draw_", as.integer(Sys.time()), "_", sample.int(1e6, 1))
    }

    cur <- drawn_polys()

    # keep label numbering stable by draw order
    next_num <- if (nrow(cur) == 0) 1L else {
      old_nums <- suppressWarnings(as.integer(sub("^Polygon\\s+", "", cur$draw_label)))
      max(old_nums, na.rm = TRUE) + 1L
    }
    if (!is.finite(next_num)) next_num <- 1L

    poly_sf$draw_id    <- fid
    poly_sf$draw_label <- paste0("Polygon ", next_num)

    if (nrow(cur) > 0) {
      cur <- cur[!(cur$draw_id %in% fid), , drop = FALSE]
    }

    drawn_polys(dplyr::bind_rows(cur, poly_sf))
    selected_draw_id(fid)

    session$sendCustomMessage("openFloating", list(id = fid))
  }, ignoreInit = TRUE)

  observeEvent(input$map_draw_deleted_features, {
    del <- input$map_draw_deleted_features
    if (is.null(del) || is.null(del$features) || length(del$features) == 0) return()

    del_ids <- vapply(del$features, function(f){
      if (!is.null(f$properties) && !is.null(f$properties$`_leaflet_id`)) {
        as.character(f$properties$`_leaflet_id`)
      } else {
        NA_character_
      }
    }, character(1))

    del_ids <- del_ids[!is.na(del_ids) & nzchar(del_ids)]
    if (length(del_ids) == 0) return()

    cur <- drawn_polys()
    cur <- cur[!(cur$draw_id %in% del_ids), , drop = FALSE]
    drawn_polys(cur)

    sid <- selected_draw_id()
    if (!is.null(sid) && sid %in% del_ids) {
      if (nrow(cur) > 0) {
        selected_draw_id(cur$draw_id[nrow(cur)])
      } else {
        selected_draw_id(NULL)
      }
    }
  })

  div_gene_initialized <- reactiveVal(FALSE)
  div_primer_initialized <- reactiveVal(FALSE)

  # Detections that fall inside any MPA/AOI polygon (drops outside points)
  detections_in_mpa <- reactive({
    yr <- sel_year_chr()

    pts <- species_sf_all
    if (yr != "All") pts <- pts %>% dplyr::filter(as.character(year) == yr)

    # spatial join: keep only detections that fall within an MPA/AOI polygon
    joined <- sf::st_join(
      pts,
      all_polys_click %>% dplyr::select(site_name, site_type),
      join = sf::st_within,
      left = FALSE
    )

    # If polygons overlap, a detection could match multiple polygons.
    # Pick the first match per detection row (safest quick fix).
    joined %>%
      arrange(occurrenceID) %>%
      dplyr::group_by(occurrenceID, samp_name, scientificName, year, target_gene) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()
  })

  # redraw points when year changes
  observeEvent(
    list(sel_year_chr(), sampling_points_layer_on()),
    {
      proxy <- leafletProxy("map")

      if (!isTRUE(sampling_points_layer_on())) {
        proxy %>% clearGroup("Sampling points")
        return(NULL)
      }

      yr <- sel_year_chr()
      pts <- sampling_pts

      # Year filter ONLY
      if (yr != "All") {
        pts <- pts %>% dplyr::filter(as.character(year) == yr)
      }

      proxy %>%
        clearGroup("Sampling points") %>%
        addCircleMarkers(
          data        = pts,
          group       = "Sampling points",
          radius      = 2,
          stroke      = TRUE,
          weight      = 1,
          opacity     = 1,
          fillOpacity = 0.8,
          options     = pathOptions(pane = "pane_points"),
          label       = ~paste0(
            "Marker: ", target_gene,
            ifelse(is.na(year), "", paste0(" | Year: ", year)),
            ifelse(is.na(samp_name), "", paste0(" | Sample: ", samp_name))
          )
        )
    },
    ignoreInit = TRUE
  )

  #Monthly sampling
  monthly_sample_counts <- reactive({
    det <- selected_detections()
    req(det)

    det0 <- det %>%
      sf::st_drop_geometry() %>%
      dplyr::mutate(
        samp_name = as.character(samp_name),
        month_chr = as.character(month)
      )

    shiny::validate(
      shiny::need(nrow(det0) > 0, "No detections in current selection.")
    )

    month_levels <- month.abb

    counts <- det0 %>%
      dplyr::filter(
        !is.na(samp_name), samp_name != "",
        !is.na(month_chr), month_chr %in% month_levels
      ) %>%
      dplyr::distinct(samp_name, month_chr) %>%
      dplyr::count(month_chr, name = "n_samples") %>%
      dplyr::mutate(month_chr = factor(month_chr, levels = month_levels))

    data.frame(
      month_chr = factor(month_levels, levels = month_levels)
    ) %>%
      dplyr::left_join(counts, by = "month_chr") %>%
      dplyr::mutate(n_samples = dplyr::coalesce(n_samples, 0L))
  })

  output$monthly_plot_subtitle <- renderText({
    click <- input$map_shape_click
    sel_id <- selected_draw_id()
    yr <- sel_year_chr()

    if (!is.null(sel_id)) {
      polys <- drawn_polys()
      hit <- polys %>% dplyr::filter(draw_id == sel_id)
      lab <- if (nrow(hit) > 0) hit$draw_label[1] else "Drawn polygon"
      return(paste0(lab, " | Year: ", yr))
    }

    if (!is.null(click$id)) {
      if (grepl("\\|\\|", click$id)) {
        parts <- strsplit(click$id, "\\|\\|")[[1]]
        return(paste0(parts[2], " | Year: ", yr))
      }

      cid <- suppressWarnings(as.integer(click$id))
      if (!is.na(cid)) {
        return(paste0("Grid cell ", cid, " | Year: ", yr))
      }
    }

    paste0("No selection | Year: ", yr)
  })

  output$monthly_circular_plot <- renderPlot({
    dat <- monthly_sample_counts()

    shiny::validate(
      shiny::need(nrow(dat) > 0, "No monthly data available."),
      shiny::need(sum(dat$n_samples, na.rm = TRUE) > 0, "No samples available for this selection.")
    )

    dat <- dat %>%
      dplyr::mutate(
        month_chr = factor(month_chr, levels = month.abb),
        fill_group = ifelse(n_samples == 0, "zero", "nonzero")
      )

    ymax <- max(dat$n_samples, na.rm = TRUE)

    # overall scale for the outer ring area
    outer_max <- max(1, ceiling(ymax * 1.10))

    # creates the empty hole in the middle
    inner_offset <- outer_max * 0.42

    # where month labels sit (further outside the circle)
    label_radius <- outer_max * 1.28

    # extra space beyond labels so they do not get clipped
    top_pad <- outer_max * 0.06

    # positions of circular guide rings
    ring_vals <- c(0, 0.25, 0.50, 0.75, 1.00) * outer_max + inner_offset

    ggplot2::ggplot(
      dat,
      ggplot2::aes(
        x = month_chr,
        y = n_samples + inner_offset,
        fill = fill_group
      )
    ) +

      # circular guide rings
      ggplot2::geom_hline(
        yintercept = ring_vals,
        colour = "grey80",
        linewidth = 0.5
      ) +

      # bars
      ggplot2::geom_col(
        width = 0.88,
        colour = NA
      ) +

      # value labels near the bar ends
      ggplot2::geom_text(
        ggplot2::aes(
          y = n_samples + inner_offset + outer_max * 0.10,
          label = n_samples
        ),
        size = 3.5,
        color = "black"
      ) +

      # month labels farther outside the rings
      ggplot2::geom_text(
        data = dat,
        ggplot2::aes(
          x = month_chr,
          y = inner_offset + label_radius,
          label = month_chr
        ),
        inherit.aes = FALSE,
        size = 5,
        color = "black"
      ) +

      ggplot2::coord_polar(start = -pi / 12) +

      ggplot2::scale_y_continuous(
        limits = c(0, inner_offset + label_radius + top_pad),
        expand = c(0, 0)
      ) +

      ggplot2::scale_fill_manual(
        values = c(
          zero = "white",     # grey for zero months
          nonzero = "#2241a7"   # blue for sampled months
        ),
        guide = "none"
      ) +

      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        axis.title = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        plot.background = ggplot2::element_rect(fill = NA, colour = NA),
        panel.background = ggplot2::element_rect(fill = NA, colour = NA),
        plot.margin = ggplot2::margin(2, 6, 2, 6)
      )
  }, res = 110)

  # your existing outputs here:
  output$cell_summary   <- renderUI({ tags$div("...") })

  # helper: pick a single display value (unique or collapse)
  pick_display <- function(x) {
    x <- as.character(x)
    x <- x[!is.na(x) & nzchar(trimws(x))]
    if (length(x) == 0) return(NA_character_)
    ux <- unique(x)
    if (length(ux) == 1) ux else paste(ux, collapse = " | ")
  }

  selected_detections <- reactive({
    yr <- sel_year_chr()

  pts <- species_sf_all
  if (yr != "All") {
  pts <- pts %>% dplyr::filter(as.character(year) == yr)
  }

    g <- selection_geom()
    if (is.null(g) || nrow(pts) == 0) {
      return(NULL)
    }

    keep <- tryCatch(
      within_any(pts, g),
      error = function(e) {
        showNotification(
          paste("Polygon selection failed:", e$message),
          type = "error"
        )
        return(rep(FALSE, nrow(pts)))
      }
    )

    pts[keep, , drop = FALSE]
  })

  selected_detections_min <- reactive({
    yr <- sel_year_chr()

    pts <- species_sf_by_year[[yr]]
    if (is.null(pts)) pts <- species_sf_min[0, ]

    g <- selection_geom()
    if (is.null(g) || nrow(pts) == 0) {
      return(NULL)
    }

    keep <- tryCatch(
      within_any(pts, g),
      error = function(e) {
        showNotification(
          paste("Polygon selection failed:", e$message),
          type = "error"
        )
        return(rep(FALSE, nrow(pts)))
      }
    )

    pts[keep, , drop = FALSE]
  })

  # ---- group filter helper (multi-select; union across selected groups) ----
  apply_group_filter <- function(occ_all, groups) {

    group_map <- list(
      "Fishes"        = list(col = "class",   vals = c("Teleostei")),
      "Sharks & Rays" = list(col = "class",   vals = c("Elasmobranchii")),
      "Mammals"       = list(col = "class",   vals = c("Mammalia")),
      "Turtles"      = list(col = "order",   vals = c("Testudines")),
      "Birds"         = list(col = "class",   vals = c("Aves")),
      "Molluscs"      = list(col = "phylum",  vals = c("Mollusca")),
      "Arthropods"    = list(col = "phylum",  vals = c("Arthropoda")),
      "Plants"        = list(col = "kingdom", vals = c("Plantae"))
    )

    groups <- as.character(groups %||% character(0))
    groups <- intersect(groups, names(group_map))
    if (length(groups) == 0) return(occ_all)

    cols_needed <- unique(vapply(group_map[groups], `[[`, character(1), "col"))
    missing_cols <- setdiff(cols_needed, names(occ_all))

    if (length(missing_cols) > 0) {
      # choose ONE behavior:
      # return(occ_all[0, , drop = FALSE])  # loud fail
      return(occ_all)                       # silent safe
    }

    keep <- rep(FALSE, nrow(occ_all))
    for (g in groups) {
      spec <- group_map[[g]]
      xcol <- tolower(as.character(occ_all[[spec$col]]))
      keep <- keep | (!is.na(xcol) & xcol %in% tolower(spec$vals))
    }

    occ_all[keep, , drop = FALSE]
  }

  # ---- active groups (multi-select) ----
  active_groups <- reactiveVal(character(0))

  toggle_group <- function(g) {
    cur <- active_groups()
    if (g %in% cur) active_groups(setdiff(cur, g)) else active_groups(c(cur, g))
  }

  observeEvent(input$total_fish,       { toggle_group("Fishes") }, ignoreInit = TRUE)
  observeEvent(input$total_sharks,     { toggle_group("Sharks & Rays") }, ignoreInit = TRUE)
  observeEvent(input$total_mammals,    { toggle_group("Mammals") }, ignoreInit = TRUE)
  observeEvent(input$total_reptiles,   { toggle_group("Turtles") }, ignoreInit = TRUE)
  observeEvent(input$total_birds,      { toggle_group("Birds") }, ignoreInit = TRUE)
  observeEvent(input$total_molluscs,   { toggle_group("Molluscs") }, ignoreInit = TRUE)
  observeEvent(input$total_arthropods, { toggle_group("Arthropods") }, ignoreInit = TRUE)
  observeEvent(input$total_plants,     { toggle_group("Plants") }, ignoreInit = TRUE)

  sync_group_button_classes <- function() {
    cur <- active_groups()

    btn_map <- c(
      total_fish       = "Fishes",
      total_sharks       = "Sharks & Rays",
      total_mammals    = "Mammals",
      total_reptiles   = "Turtles",
      total_birds      = "Birds",
      total_molluscs   = "Molluscs",
      total_arthropods = "Arthropods",
      total_plants     = "Plants"
    )

    for (btn_id in names(btn_map)) {
      g <- btn_map[[btn_id]]
      if (g %in% cur) shinyjs::addClass(btn_id, "btn-group-on") else shinyjs::removeClass(btn_id, "btn-group-on")
    }
  }

  # IMPORTANT: keep button visuals synced with state
  observeEvent(active_groups(), {
    sync_group_button_classes()
  }, ignoreInit = FALSE)


  # helper: current group label for display/debug if needed
  active_groups_label <- reactive({
    cur <- active_groups()
    if (length(cur) == 0) "All" else paste(cur, collapse = " + ")
  })

  # ---- SARA/AIS sets + toggles (assumes SARA & AIS exist globally) ----
  sara_set <- reactive(unique(na.omit(SARA$Scientific.Name)))
  ais_set  <- reactive(unique(na.omit(AIS$Scientific.Name)))  # adjust column if needed

  filter_sara_on <- reactiveVal(FALSE)
  filter_ais_on  <- reactiveVal(FALSE)

  apply_species_filters <- function(occ_all) {

    # 1) group buttons
    occ_all <- apply_group_filter(occ_all, active_groups())

    # 2) SARA/AIS union logic
    if ("scientificName" %in% names(occ_all)) {
      spp_keep <- character(0)

      if (isTRUE(filter_sara_on())) spp_keep <- union(spp_keep, sara_set())
      if (isTRUE(filter_ais_on()))  spp_keep <- union(spp_keep, ais_set())

      if (length(spp_keep) > 0) {
        occ_all <- occ_all %>% dplyr::filter(scientificName %in% spp_keep)
      }
    }

    occ_all
  }

  diversity_dropdown_data <- reactive({
    pts <- species_sf_all

    yr <- sel_year_chr()
    if (yr != "All") {
      pts <- pts %>% dplyr::filter(as.character(year) == yr)
    }

    pts <- apply_species_filters(pts)

    pts %>%
      sf::st_drop_geometry() %>%
      add_primer_combo() %>%
      dplyr::mutate(
        target_gene  = as.character(target_gene),
        primer_combo = as.character(primer_combo)
      )
  })

  # ---- source data for diversity dropdowns ----
  observeEvent(diversity_dropdown_data(), {
    dd <- diversity_dropdown_data()

    gene_choices <- dd %>%
      dplyr::filter(!is.na(target_gene), target_gene != "") %>%
      dplyr::pull(target_gene) %>%
      unique() %>%
      sort()

    cur_gene <- input$div_target_gene %||% character(0)
    sel_gene <- intersect(cur_gene, gene_choices)

    # first load = select all genes
    if (!div_gene_initialized()) {
      sel_gene <- gene_choices
      div_gene_initialized(TRUE)
    }

    freezeReactiveValue(input, "div_target_gene")
    updateSelectizeInput(
      session  = session,
      inputId  = "div_target_gene",
      choices  = gene_choices,
      selected = sel_gene,
      server   = TRUE
    )
  }, ignoreInit = FALSE)

  primer_choices_reactive <- reactive({
    dd <- diversity_dropdown_data()

    genes_selected <- input$div_target_gene %||% character(0)

    if (length(genes_selected) > 0) {
      dd <- dd %>%
        dplyr::filter(target_gene %in% genes_selected)
    } else {
      dd <- dd[0, , drop = FALSE]
    }

    dd %>%
      dplyr::filter(!is.na(primer_combo), primer_combo != "") %>%
      dplyr::pull(primer_combo) %>%
      unique() %>%
      sort()
  })

  # Primer choices
  observeEvent(
    list(primer_choices_reactive(), input$div_target_gene),
    {
      primer_choices <- primer_choices_reactive()

      cur_primer <- input$div_primer %||% character(0)
      sel_primer <- intersect(cur_primer, primer_choices)

      # first load = select all available primers
      if (!div_primer_initialized()) {
        sel_primer <- primer_choices
        div_primer_initialized(TRUE)
      }

      freezeReactiveValue(input, "div_primer")
      updateSelectizeInput(
        session  = session,
        inputId  = "div_primer",
        choices  = primer_choices,
        selected = sel_primer,
        server   = TRUE
      )
    },
    ignoreInit = FALSE
  )

  observeEvent(input$div_primer_all, {
    primer_choices <- isolate(primer_choices_reactive())

    updateSelectizeInput(
      session  = session,
      inputId  = "div_primer",
      choices  = primer_choices,
      selected = primer_choices,
      server   = TRUE
    )
  })

  observeEvent(input$div_primer_none, {
    primer_choices <- isolate(primer_choices_reactive())

    updateSelectizeInput(
      session  = session,
      inputId  = "div_primer",
      choices  = primer_choices,
      selected = character(0),
      server   = TRUE
    )
  })

  # ---- confirmed diversity controls ----
  div_filters <- eventReactive(input$div_apply, {
    list(
      tax_rank    = input$tax_rank %||% "scientificName",
      target_gene = input$div_target_gene %||% character(0),
      primers     = input$div_primer %||% character(0)
    )
  }, ignoreInit = FALSE)

  add_polygon_selection <- function(pts_sf) {
    if (is.null(pts_sf) || nrow(pts_sf) == 0) {
      pts_sf$polygon_selection <- character(0)
      return(pts_sf)
    }

    # start with empty labels
    poly_labels <- rep(NA_character_, nrow(pts_sf))

    # ---- existing MPA/AOI polygons ----
    if (!is.null(all_polys_click) && nrow(all_polys_click) > 0) {
      mpa_hits <- sf::st_intersects(pts_sf, all_polys_click)

      mpa_labels <- vapply(seq_along(mpa_hits), function(i) {
        hit <- mpa_hits[[i]]
        if (length(hit) == 0) return(NA_character_)

        labs <- paste0(
          all_polys_click$site_type[hit],
          ": ",
          all_polys_click$site_name[hit]
        )

        paste(unique(labs), collapse = " | ")
      }, character(1))

      poly_labels <- mpa_labels
    }

    # ---- user-drawn polygons ----
    polys_drawn <- drawn_polys()
    if (!is.null(polys_drawn) && nrow(polys_drawn) > 0) {
      drawn_hits <- sf::st_intersects(pts_sf, polys_drawn)

      drawn_labels <- vapply(seq_along(drawn_hits), function(i) {
        hit <- drawn_hits[[i]]
        if (length(hit) == 0) return(NA_character_)

        labs <- polys_drawn$draw_label[hit]
        paste(unique(labs), collapse = " | ")
      }, character(1))

      poly_labels <- ifelse(
        !is.na(poly_labels) & !is.na(drawn_labels),
        paste(poly_labels, drawn_labels, sep = " | "),
        dplyr::coalesce(poly_labels, drawn_labels)
      )
    }

    pts_sf$polygon_selection <- poly_labels
    pts_sf
  }


  download_data_reactive <- reactive({
    filters <- div_filters()
    yr <- sel_year_chr()

    dat <- species_sf_all %>%
      add_polygon_selection() %>%
      dplyr::mutate(
        decimalLongitude = sf::st_coordinates(.)[, 1],
        decimalLatitude  = sf::st_coordinates(.)[, 2]
      ) %>%
      sf::st_drop_geometry() %>%
      add_primer_combo()

    if (yr != "All") {
      dat <- dat %>%
        dplyr::filter(as.character(year) == yr)
    }

    dat <- apply_species_filters(dat)
    dat <- apply_diversity_dropdown_filters(dat, filters)

    rank_col <- filters$tax_rank %||% "scientificName"
    if (rank_col %in% names(dat)) {
      dat <- dat %>%
        dplyr::mutate(selected_taxonomy = as.character(.data[[rank_col]]))
    } else {
      dat <- dat %>%
        dplyr::mutate(selected_taxonomy = NA_character_)
    }

    dat
  })

  species_list_occ_all <- reactive({
    occ_all <- selected_detections()   # or occ_all_filtered(), etc.
    req(occ_all)

    occ_all <- apply_species_filters(occ_all)

    occ_all %>% dplyr::distinct(scientificName, .keep_all = TRUE) %>% dplyr::arrange(dplyr::coalesce(worms_valid_name, scientificName))
  })

  observeEvent(input$SARA, {
    new_state <- !isTRUE(filter_sara_on())
    filter_sara_on(new_state)
    if (new_state) shinyjs::addClass("SARA", "btn-sara-on") else shinyjs::removeClass("SARA", "btn-sara-on")
  })

  observeEvent(input$AIS, {
    new_state <- !isTRUE(filter_ais_on())
    filter_ais_on(new_state)
    if (new_state) shinyjs::addClass("AIS", "btn-ais-on") else shinyjs::removeClass("AIS", "btn-ais-on")
  })

    # UI

  apply_interest_filter <- function(spp_vec) {
    if (!isTRUE(filter_sara_on()) && !isTRUE(filter_ais_on())) {
      return(spp_vec)
    }

    keep <- character(0)

    if (isTRUE(filter_sara_on())) keep <- union(keep, sara_set())
    if (isTRUE(filter_ais_on()))  keep <- union(keep, ais_set())

    intersect(spp_vec, keep)
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
        data.frame(Message = "Select the “SARA” button to view species at risk details."),
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
      dplyr::filter(scientificName %in% sara_set())

    if (nrow(det) == 0) {
      return(DT::datatable(
        data.frame(Message = "No SARA Schedule 1 detections for this selection."),
        rownames = FALSE,
        options = list(dom = "t")
      ))
    }

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
        Common.Name  = paste(sort(unique(na.omit(Common.Name))), collapse = " |OR| "),
        n_detections = dplyr::n(),
        n_samples    = dplyr::n_distinct(samp_name),
        samples      = paste(sort(unique(samp_name)), collapse = ", "),
        years        = paste(sort(unique(na.omit(year))), collapse = ", "),
        months       = paste(sort(unique(na.omit(month))), collapse = ", "),
        markers      = paste(sort(unique(na.omit(target_gene))), collapse = ", "),
        primers      = paste(
          sort(unique(na.omit(
            paste(
              trimws(pcr_primer_name_forward),
              trimws(pcr_primer_name_reverse),
              sep = " | "
            )
          ))),
          collapse = ", "
        ),
        .groups = "drop"
      ) %>%
      dplyr::select(
        scientificName,
        Common.Name,
        Rating,
        n_detections,
        n_samples,
        samples,
        years,
        months,
        markers,
        primers
      ) %>%
      dplyr::arrange(Rating, scientificName)

    DT::datatable(
      out,
      rownames = FALSE,
      colnames = c(
        "Species",
        "Common Name",
        "SARA Rating",
        "Number of Detections",
        "Number of Samples",
        "Samples",
        "Years",
        "Months",
        "Target Genes",
        "Primers"
      ),
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        autoWidth = TRUE
      ),
      class = "nowrap"
    )
  })

  output$ais_details <- DT::renderDT({
    if (!isTRUE(filter_ais_on())) {
      return(DT::datatable(
        data.frame(Message = "Select the “AIS” button to view invasive species details."),
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

    det2 <- det %>%
      dplyr::left_join(
        AIS %>% dplyr::select(Scientific.Name, Common.Name, Type),
        by = c("scientificName" = "Scientific.Name")
      )

    out <- det2 %>%
      dplyr::group_by(scientificName) %>%
      dplyr::summarise(
        Common.Name   = paste(sort(unique(na.omit(Common.Name))), collapse = " |OR| "),
        Type          = paste(sort(unique(na.omit(Type))), collapse = " |OR| "),
        n_detections  = dplyr::n(),
        n_samples     = dplyr::n_distinct(samp_name),
        samples       = paste(sort(unique(samp_name)), collapse = ", "),
        years         = paste(sort(unique(na.omit(year))), collapse = ", "),
        months        = paste(sort(unique(na.omit(month))), collapse = ", "),
        markers       = paste(sort(unique(na.omit(target_gene))), collapse = ", "),
        primers       = paste(
          sort(unique(na.omit(
            paste(
              trimws(pcr_primer_name_forward),
              trimws(pcr_primer_name_reverse),
              sep = " | "
            )
          ))),
          collapse = ", "
        ),
        .groups = "drop"
      ) %>%
      dplyr::select(
        scientificName,
        Common.Name,
        Type,
        n_detections,
        n_samples,
        samples,
        years,
        months,
        markers,
        primers
      ) %>%
      dplyr::arrange(scientificName)

    DT::datatable(
      out,
      rownames = FALSE,
      colnames = c(
        "Species",
        "Common Name",
        "Type",
        "Number of Detections",
        "Number of Samples",
        "Samples",
        "Years",
        "Months",
        "Target Genes",
        "Primers"
      ),
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        autoWidth = TRUE
      ),
      class = "nowrap"
    )
  })

  output$downloadData <- downloadHandler(
    filename = function() {
      paste0("GOTeDNA_filtered_data_", Sys.Date(), ".csv")
    },
    content = function(file) {
      out <- download_data_reactive()

      # optional: put a few important columns first, but keep ALL columns
      preferred_cols <- c(
        "occurrenceID", "scientificName", "selected_taxonomy",
        "kingdom", "phylum", "class", "order", "family", "genus",
        "samp_name", "year", "month", "eventDate",
        "target_gene", "primer_combo",
        "organismQuantity", "organismQuantityType",
        "decimalLatitude", "decimalLongitude", "polygon_selection",
        "minimumDepthInMeters", "maximumDepthInMeters",
        "samp_size", "samp_size_unit"
      )

      keep_cols <- c(
        intersect(preferred_cols, names(out)),
        setdiff(names(out), preferred_cols)
      )

      out <- out %>%
        dplyr::select(dplyr::all_of(keep_cols))

      utils::write.csv(out, file, row.names = FALSE, na = "")
    }
  )

  # ---- 1) Render the leaflet map ONCE ----
  output$map <- renderLeaflet({

    # ---- choose initial year + initial layers safely ----
    yrs <- sort(unique(na.omit(as.character(KEY_TBL$year))))

    init_12S <- if (default_year == "All") RICHNESS_GENE_ALL[["12S"]] else RICHNESS_BY_KEY[[paste0("12S_", default_year)]]
    init_COI <- if (default_year == "All") RICHNESS_GENE_ALL[["COI"]] else RICHNESS_BY_KEY[[paste0("COI_", default_year)]]
    init_16S <- if (default_year == "All") RICHNESS_GENE_ALL[["16S"]] else RICHNESS_BY_KEY[[paste0("16S_", default_year)]]
    init_18S <- if (default_year == "All") RICHNESS_GENE_ALL[["18S"]] else RICHNESS_BY_KEY[[paste0("18S_", default_year)]]

    init_ALL <- {
      if (default_year == "All") {
        RICHNESS_ALL
      } else if (default_year %in% names(RICHNESS_ALL_BY_YEAR)) {
        RICHNESS_ALL_BY_YEAR[[default_year]]
      } else {
        RICHNESS_ALL
      }
    }

    # ---- defensive: ensure has_sampling exists where you use it ----
    if (!is.null(init_12S) && !"has_sampling" %in% names(init_12S)) init_12S$has_sampling <- FALSE
    if (!is.null(init_COI) && !"has_sampling" %in% names(init_COI)) init_COI$has_sampling <- FALSE
    if (!is.null(init_16S) && !"has_sampling" %in% names(init_16S)) init_16S$has_sampling <- FALSE
    if (!is.null(init_18S) && !"has_sampling" %in% names(init_18S)) init_18S$has_sampling <- FALSE
    if (!is.null(init_ALL) && !"has_sampling" %in% names(init_ALL)) init_ALL$has_sampling <- FALSE

    m <- leaflet() %>%
      addMapPane("pane_polys",      zIndex = 400) %>%
      addMapPane("pane_zones",      zIndex = 410) %>%
      addMapPane("pane_poly_total", zIndex = 420) %>%
      addMapPane("pane_points",     zIndex = 430) %>%
      addMapPane("pane_drawn_top",    zIndex = 900) %>%   #new code
      addMapPane("pane_selected_top", zIndex = 950) %>%   #new code
      addProviderTiles(providers$CartoDB.Positron, group = "CartoDB Positron") %>%
      addProviderTiles(providers$Esri.OceanBasemap, group = "Esri Ocean Basemap") %>%
      addProviderTiles(providers$Esri.WorldImagery, group = "Esri World Imagery")

    # ---- add initial richness layers if they exist ----
    if (!is.null(init_12S) && nrow(init_12S) > 0) {
      m <- m %>% addPolygons(
        data        = init_12S,
        group       = "12S",
        layerId     = ~cell_id,
        fillColor   = ~pal_rich(n_species),
        fillOpacity = ~ifelse(has_sampling, 0.8, 0.04),
        color       = NA,
        label       = ~ifelse(has_sampling, paste("12S richness:", n_species), "No sampling in this cell"),
        options     = pathOptions(pane = "pane_polys")
      )
    }

    if (!is.null(init_COI) && nrow(init_COI) > 0) {
      m <- m %>% addPolygons(
        data        = init_COI,
        group       = "COI",
        layerId     = ~cell_id,
        fillColor   = ~pal_rich(n_species),
        fillOpacity = ~ifelse(has_sampling, 0.8, 0.04),
        color       = NA,
        label       = ~ifelse(has_sampling, paste("COI richness:", n_species), "No sampling in this cell"),
        options     = pathOptions(pane = "pane_polys")
      )
    }

    if (!is.null(init_16S) && nrow(init_16S) > 0) {
      m <- m %>% addPolygons(
        data        = init_16S,
        group       = "16S",
        layerId     = ~cell_id,
        fillColor   = ~pal_rich(n_species),
        fillOpacity = ~ifelse(has_sampling, 0.8, 0.04),
        color       = NA,
        label       = ~ifelse(has_sampling, paste("16S richness:", n_species), "No sampling in this cell"),
        options     = pathOptions(pane = "pane_polys")
      )
    }

    if (!is.null(init_18S) && nrow(init_18S) > 0) {
      m <- m %>% addPolygons(
        data        = init_18S,
        group       = "18S",
        layerId     = ~cell_id,
        fillColor   = ~pal_rich(n_species),
        fillOpacity = ~ifelse(has_sampling, 0.8, 0.04),
        color       = NA,
        label       = ~ifelse(has_sampling, paste("18S richness:", n_species), "No sampling in this cell"),
        options     = pathOptions(pane = "pane_polys")
      )
    }
    # Always add ALL
    m <- m %>% addPolygons(
      data        = init_ALL,
      group       = "All",
      layerId     = ~cell_id,
      fillColor   = ~pal_rich(n_species_total),
      fillOpacity = ~ifelse(has_sampling, 0.8, 0.04),
      color       = NA,
      label       = ~ifelse(has_sampling, paste("Total richness:", n_species_total), "No sampling in this cell"),
      options     = pathOptions(pane = "pane_polys"),
      highlightOptions = highlightOptions(weight = 2, bringToFront = TRUE)
    )

    # ---- rest of your map layers ----
    m <- m %>%
      addPolygons(
        data = depth_layers[["All"]],
        group = "Sampling Depth",
        layerId = ~cell_id,
        fillColor   = ~final_fill,
        fillOpacity = ~final_alpha,
        opacity     = 1,
        color       = NA,
        label       = ~paste0(
          "Depth (m, median): ", ifelse(is.na(depth_med), "NA", round(depth_med, 1)), "<br>",
          "Min: ", ifelse(is.na(depth_min), "NA", round(depth_min, 1)), " | ",
          "Max: ", ifelse(is.na(depth_max), "NA", round(depth_max, 1)), "<br>"
        ) %>% lapply(htmltools::HTML),
        highlightOptions = leaflet::highlightOptions(weight = 2, bringToFront = TRUE),
        options = pathOptions(pane = "pane_polys")
      ) %>%
      addCircleMarkers(
        data        = sampling_pts,
        group       = "Sampling points",
        radius      = 2,
        stroke      = TRUE,
        weight      = 1,
        opacity     = 1,
        fillOpacity = 0.8,
        options = pathOptions(pane = "pane_points"),
        label       = ~paste0(
          "Marker: ", target_gene,
          ifelse(is.na(year), "", paste0(" | Year: ", year)),
          #ifelse(is.na(source_file), "", paste0(" | File: ", source_file)),
          ifelse(is.na(samp_name), "", paste0(" | Sample: ", samp_name))
        )
      ) %>%
      addPolygons(
        data        = all_polys_zones,
        group       = "MPA/AOI zone boundaries",
        fillOpacity = 0,
        color       = "black",
        weight      = 1,
        opacity     = 0.8,
        options     = pathOptions(clickable = FALSE, pane = "pane_zones")
      ) %>%
      addPolygons(
        data        = all_polys_click,
        group       = "Total species detected per MPA/AOI",
        layerId     = ~paste(site_type, site_name, sep="||"),
        fillOpacity = 0.05,
        color       = "black",
        weight      = 2,
        opacity     = 1,
        popup       = ~site_name,
        options     = pathOptions(pane = "pane_poly_total"),
        highlightOptions = highlightOptions(weight = 3, bringToFront = TRUE)
      )

    # ---- legends + controls ----
    m <- m %>%
      leaflet::addLegend(
        position = "bottomright",
        pal      = pal_rich,
        values   = init_ALL$n_species_total %||% numeric(0),
        title    = "Species richness from eDNA",
        opacity  = 1,
        className = "legend-base legend-richness-box"
      ) %>%
      leaflet::addLegend(
        position  = "bottomright",
        colors    = c(depth_legend_cols, "#feb24c"),
        labels    = c(depth_legend_labs, "Mixed depth (orange overlay)"),
        title     = "Sampling depth",
        opacity   = 1,
        className = "legend-base legend-depth-box"
      ) %>%
      addCircleMarkers(
        lng = 0, lat = 0,
        radius = 1,
        opacity = 0,
        fillOpacity = 0,
        stroke = FALSE,
        group = "Monthly Sampling Plot"
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
          "All", "12S", "COI", "16S", "18S",
          "MPA/AOI zone boundaries",
          "Sampling points",
          "Sampling Depth",
          "Monthly Sampling Plot"
        ),
        options = layersControlOptions(collapsed = FALSE)
      # )
      # addLayersControl(
      #   baseGroups = c("CartoDB Positron", "Esri Ocean Basemap", "Esri World Imagery"),
      #   overlayGroups = c(
      #     "Total species detected per MPA/AOI",
      #     "All", "12S", "COI", "16S", "18S",
      #     "MPA/AOI zone boundaries",
      #     "Sampling points",
      #     "Sampling Depth"
      #   ),
      #   options = layersControlOptions(collapsed = FALSE)
      ) %>%
      htmlwidgets::onRender("
                      function(el, x){

                        const richness = new Set([
                          'All',
                          '12S',
                          'COI',
                          '16S',
                          '18S'
                        ]);

                        const context = new Set([
                          'MPA/AOI zone boundaries',
                          'Sampling points',
                          'Sampling Depth',
                          'Monthly Sampling Plot'
                        ]);

                        const depthName = 'Sampling Depth';

                        const monthlyPlotName = 'Monthly Sampling Plot';

                        function getOverlayRows(){
                          const ctrl = el.querySelector('.leaflet-control-layers');
                          if(!ctrl) return [];
                          const rows = ctrl.querySelectorAll('.leaflet-control-layers-overlays label');
                          return Array.from(rows);
                        }

                        function labelText(row){
                          return row.textContent.replace(/\\s+/g,' ').trim();
                        }

                        function inputOf(row){
                          return row.querySelector('input[type=checkbox]');
                        }

                        function clickOffByName(name){
                          const rows = getOverlayRows();
                          rows.forEach(r => {
                            if(labelText(r) === name){
                              const cb = inputOf(r);
                              if(cb && cb.checked) cb.click();
                            }
                          });
                        }

                        function clickOffSet(nameSet){
                          const rows = getOverlayRows();
                          rows.forEach(r => {
                            const nm = labelText(r);
                            if(nameSet.has(nm)){
                              const cb = inputOf(r);
                              if(cb && cb.checked) cb.click();
                            }
                          });
                        }

                        function anyChecked(nameSet){
                          const rows = getOverlayRows();
                          for(const r of rows){
                            const nm = labelText(r);
                            const cb = inputOf(r);
                            if(cb && cb.checked && nameSet.has(nm)) return true;
                          }
                          return false;
                        }

                        function isChecked(name){
                          const rows = getOverlayRows();
                          for(const r of rows){
                            if(labelText(r) === name){
                              const cb = inputOf(r);
                              return cb ? cb.checked : false;
                            }
                          }
                          return false;
                        }

                       function updateMonthlyPlotVisibility(){
  const plotBox = document.getElementById('monthly_plot_control');
  if(!plotBox) return;

  const monthlyOn = isChecked(monthlyPlotName);

  const richLegend  = el.querySelector('.legend-richness-box');
  const depthLegend = el.querySelector('.legend-depth-box');

  let activeLegend = null;
  if(richLegend && !richLegend.classList.contains('legend-hidden')){
    activeLegend = richLegend;
  } else if(depthLegend && !depthLegend.classList.contains('legend-hidden')){
    activeLegend = depthLegend;
  }

  const pad = 14;
  const gap = 14;
  const defaultLegendRight = 10;
  const defaultLegendBottom = 10;

  // Monthly plot OFF
  if(!monthlyOn){
    plotBox.style.display = 'none';
    plotBox.style.left = 'auto';
    plotBox.style.right = 'auto';
    plotBox.style.top = 'auto';
    plotBox.style.bottom = 'auto';

    if(richLegend){
      richLegend.style.left = 'auto';
      richLegend.style.right = defaultLegendRight + 'px';
      richLegend.style.top = 'auto';
      richLegend.style.bottom = defaultLegendBottom + 'px';
      richLegend.style.margin = '0px';
    }

    if(depthLegend){
      depthLegend.style.left = 'auto';
      depthLegend.style.right = defaultLegendRight + 'px';
      depthLegend.style.top = 'auto';
      depthLegend.style.bottom = defaultLegendBottom + 'px';
      depthLegend.style.margin = '0px';
    }

    return;
  }

  // Monthly plot ON
  plotBox.style.display = 'block';
  plotBox.style.left = 'auto';
  plotBox.style.right = pad + 'px';
  plotBox.style.top = 'auto';
  plotBox.style.bottom = pad + 'px';
  plotBox.style.margin = '0px';

  if(activeLegend){
    const plotW = plotBox.offsetWidth || 320;

    activeLegend.style.left = 'auto';
    activeLegend.style.right = (plotW + gap + pad) + 'px';
    activeLegend.style.top = 'auto';
    activeLegend.style.bottom = pad + 'px';
    activeLegend.style.margin = '0px';
  }
}

                        function insertHeadings(){
                          const ctrl = el.querySelector('.leaflet-control-layers');
                          if(!ctrl) return;
                          if(ctrl.querySelector('.layers-heading')) return;

                          const overlayBox = ctrl.querySelector('.leaflet-control-layers-overlays');
                          if(!overlayBox) return;

                          const rows = getOverlayRows();
                          if(rows.length === 0) return;

                          let firstRich = null, firstCtx = null;
                          rows.forEach(r => {
                            const nm = labelText(r);
                            if(!firstRich && richness.has(nm)) firstRich = r;
                            if(!firstCtx  && context.has(nm))  firstCtx  = r;
                          });

                          if(firstRich){
                            const h1 = document.createElement('div');
                            h1.className = 'layers-heading';
                            h1.textContent = 'Species richness by gene region (grid cells)';
                            overlayBox.insertBefore(h1, firstRich);
                          }

                          if(firstCtx){
                            const sep = document.createElement('div');
                            sep.className = 'layers-sep';
                            overlayBox.insertBefore(sep, firstCtx);

                            const h2 = document.createElement('div');
                            h2.className = 'layers-heading';
                            h2.textContent = 'Optional Layers';
                            overlayBox.insertBefore(h2, firstCtx);
                          }
                        }

                        // ---- LEGEND SWAP ---- //

                        function updateLegends(){
                        const richLegend  = el.querySelector('.legend-richness-box');
                        const depthLegend = el.querySelector('.legend-depth-box');

                        const depthOn = isChecked(depthName);
                        const anyRichOn = anyChecked(richness);

                        if(depthLegend){
                        depthLegend.classList.toggle('legend-hidden', !depthOn);
                        }

                        if(richLegend){
                        richLegend.classList.toggle('legend-hidden', depthOn || !anyRichOn);
                        }

                        updateMonthlyPlotVisibility();
                        }

                        function wireExclusivity(){
                          const rows = getOverlayRows();
                          rows.forEach(r => {
                            const nm = labelText(r);
                            const cb = inputOf(r);
                            if(!cb) return;

                            // prevent duplicate listeners if MutationObserver fires
                            if(cb.dataset.wired === '1') return;
                            cb.dataset.wired = '1';

                            cb.addEventListener('change', function(){

                              // Richness ON -> Depth OFF
                              if(this.checked && richness.has(nm)){
                                if(isChecked(depthName)) clickOffByName(depthName);
                              }

                              // Depth ON -> all Richness OFF
                              if(this.checked && nm === depthName){
                                clickOffSet(richness);
                              }

                              updateLegends();
                            });
                          });
                        }

                        // ---- FORCE DEPTH OFF AT STARTUP ----
                          function forceDepthOffStartup(){
                            clickOffByName(depthName);
                            updateLegends();
                          }

                        insertHeadings();
                        wireExclusivity();

                        // ensure initial state after the control fully exists
                        setTimeout(forceDepthOffStartup, 0);

                        const obs = new MutationObserver(() => {
                          insertHeadings();
                          wireExclusivity();
                          updateLegends();
                        });
                        obs.observe(el, {childList:true, subtree:true});
                      }
                      ")

    m
  })

  observe({
    polys <- drawn_polys()
    proxy <- leafletProxy("map")

    proxy %>% clearGroup("Drawn polygons")
    proxy %>% clearGroup("Selected drawn polygon")

    if (nrow(polys) == 0) return()

    # all drawn polygons
    proxy %>%
      addPolygons(
        data        = polys,
        group       = "Drawn polygons",
        layerId     = ~draw_id,
        fillColor   = "#2241a7",
        fillOpacity = 0.10,
        color       = "#2241a7",
        weight      = 2,
        opacity     = 1,
        label       = ~draw_label,
        options     = pathOptions(pane = "pane_drawn_top")
      )

    # selected polygon highlighted separately
    sid <- selected_draw_id()
    if (!is.null(sid) && sid %in% polys$draw_id) {
      sel <- polys %>% dplyr::filter(draw_id == sid)

      proxy %>%
        addPolygons(
          data        = sel,
          group       = "Selected drawn polygon",
          layerId     = ~draw_id,
          fillColor   = "#Fdd262",
          fillOpacity = 0.18,
          color       = "#Fdd262",
          weight      = 4,
          opacity     = 1,
          label       = ~draw_label,
          options     = pathOptions(pane = "pane_selected_top")
        )
    }
  })

  # ---- 2) When the year changes, swap the richness layers ----
  observeEvent(sel_year_chr(), {
    yr <- sel_year_chr()
    proxy <- leafletProxy("map")

    add_grid_layer <- function(data_sf, group_name, value_col, label_prefix) {
      if (is.null(data_sf) || nrow(data_sf) == 0) {
        proxy %>% clearGroup(group_name)
        return(invisible(NULL))
      }
      if (!"has_sampling" %in% names(data_sf)) data_sf$has_sampling <- FALSE
      data_sf$has_sampling <- as.logical(data_sf$has_sampling)
      if (!value_col %in% names(data_sf)) data_sf[[value_col]] <- NA_real_

      data_sf$.val   <- data_sf[[value_col]]
      data_sf$.label <- ifelse(
        data_sf$has_sampling,
        paste0(label_prefix, ": ", data_sf$.val),
        "No sampling in this cell"
      )

      proxy %>%
        clearGroup(group_name) %>%
        addPolygons(
          data        = data_sf,
          group       = group_name,
          layerId     = ~cell_id,
          fillColor   = ~pal_rich(.val),
          fillOpacity = ~ifelse(has_sampling, 0.8, 0.04),
          color       = NA,
          label       = ~.label,
          options     = pathOptions(pane = "pane_polys"),
          highlightOptions = highlightOptions(weight = 2, bringToFront = TRUE)
        )
    }

    grid12  <- if (yr == "All") RICHNESS_GENE_ALL[["12S"]] else RICHNESS_BY_KEY[[paste0("12S_", yr)]]
    gridCOI <- if (yr == "All") RICHNESS_GENE_ALL[["COI"]] else RICHNESS_BY_KEY[[paste0("COI_", yr)]]
    grid16  <- if (yr == "All") RICHNESS_GENE_ALL[["16S"]] else RICHNESS_BY_KEY[[paste0("16S_", yr)]]
    grid18  <- if (yr == "All") RICHNESS_GENE_ALL[["18S"]] else RICHNESS_BY_KEY[[paste0("18S_", yr)]]

    add_grid_layer(grid12,  "12S", "n_species", "12S richness")
    add_grid_layer(gridCOI, "COI", "n_species", "COI richness")
    add_grid_layer(grid16,  "16S", "n_species", "16S richness")
    add_grid_layer(grid18,  "18S", "n_species", "18S richness")

    gridALL <- if (yr == "All") RICHNESS_ALL else RICHNESS_ALL_BY_YEAR[[yr]] %||% RICHNESS_ALL
    add_grid_layer(gridALL, "All", "n_species_total", "Total richness")
  }, ignoreInit = TRUE)

  #Depth toggle
  selected_depth_layer <- function(year) {
    if (identical(year, "All")) return("All")
    as.character(year)
  }

  observe({
    req(input$sel_year)

    yr_key <- selected_depth_layer(input$sel_year)
    if (!yr_key %in% names(depth_layers)) return()

    occ_all <- depth_layers[[yr_key]]

    leafletProxy("map") %>%
      clearGroup("Sampling Depth") %>%
      addPolygons(
        data        = occ_all,
        group       = "Sampling Depth",
        layerId     = ~cell_id,
        fill        = TRUE,
        fillColor   = ~final_fill,
        fillOpacity = ~final_alpha,
        opacity     = 0,
        color       = NA,
        label       = ~paste0(
          "Depth (m, median): ", ifelse(is.na(depth_med), "NA", round(depth_med, 1)), "<br>",
          "Min: ", ifelse(is.na(depth_min), "NA", round(depth_min, 1)), " | ",
          "Max: ", ifelse(is.na(depth_max), "NA", round(depth_max, 1))
        ) %>% lapply(htmltools::HTML)
      )
  })

  # ---- clicked cell helper ----
  clicked_cell <- reactive({
    click <- input$map_shape_click
    if (is.null(click) || is.null(click$id)) return(NA_integer_)
    suppressWarnings(as.integer(click$id))
  })

  # ---- detections used by Diversity section only ----
  # Floating panel filters still apply first.
  # tax_rank and div_target_gene only affect diversity calculations.
  # diversity_detections <- reactive({
  #   det <- selected_detections()
  #   req(det)
  #
  #   det <- apply_species_filters(det)
  #
  #   gene_sel <- input$div_target_gene %||% character(0)
  #
  #   # If user deselects everything, keep all genes
  #   if (length(gene_sel) > 0) {
  #     det <- det %>%
  #       dplyr::filter(as.character(target_gene) %in% gene_sel)
  #   }
  #
  #   det
  # })

  # ---- diversity detections for all MPA/AOI polygons ----
  diversity_detections_mpa <- reactive({
    yr <- sel_year_chr()
    filters <- div_filters()

    pts <- species_sf_all
     if (yr != "All") {
       pts <- pts %>% dplyr::filter(as.character(year) == yr)
     }

    joined <- sf::st_join(
      pts,
      all_polys_click %>% dplyr::select(site_name, site_type),
      join = sf::st_within,
      left = FALSE
    )

    joined <- apply_species_filters(joined)
    joined <- apply_diversity_dropdown_filters(joined, filters)

    joined %>%
      dplyr::arrange(occurrenceID) %>%
      dplyr::group_by(
        occurrenceID, samp_name, scientificName, year, target_gene,
        site_name, site_type
      ) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()
  })

  #Diversity plots
  # --- build sample x taxon matrix from current selection ---
  comm_mat_mpa <- reactive({
    det <- diversity_detections_mpa()
    req(det)

    occ_all <- det %>% sf::st_drop_geometry()

    shiny::validate(
      shiny::need(nrow(occ_all) > 0, "No detections available for the current filters.")
    )

    rank_col <- div_filters()$tax_rank

    shiny::validate(
      shiny::need("organismQuantity" %in% names(occ_all),
                  "organismQuantity column is missing.")
    )

    occ_all2 <- occ_all %>%
      dplyr::mutate(
        samp_name = as.character(samp_name),
        taxon     = as.character(.data[[rank_col]]),
        value     = as.numeric(organismQuantity)
      ) %>%
      dplyr::filter(
        !is.na(samp_name), samp_name != "",
        !is.na(site_name), site_name != "",
        !is.na(taxon), taxon != "",
        !is.na(value)
      ) %>%
      dplyr::group_by(site_name, samp_name, taxon) %>%
      dplyr::summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
      dplyr::mutate(sample_group = paste(site_name, samp_name, sep = " || "))

    shiny::validate(
      shiny::need(nrow(occ_all2) > 0, "No abundance data available after filtering.")
    )

    mat_wide <- occ_all2 %>%
      dplyr::select(sample_group, taxon, value) %>%
      tidyr::pivot_wider(
        names_from  = taxon,
        values_from = value,
        values_fill = 0
      )

    mat <- mat_wide %>%
      dplyr::select(-sample_group) %>%
      as.data.frame()

    rownames(mat) <- mat_wide$sample_group
    mat
  })


  sample_meta_mpa <- reactive({
    det <- diversity_detections_mpa()
    req(det)

    det %>%
      sf::st_drop_geometry() %>%
      dplyr::mutate(
        samp_name = as.character(samp_name),
        site_name = as.character(site_name),
        site_type = as.character(site_type),
        year      = if ("year" %in% names(.)) as.character(year) else NA_character_,
        sample_group = paste(site_name, samp_name, sep = " || ")
      ) %>%
      dplyr::distinct(sample_group, samp_name, site_name, site_type, year)
  })

  # ---- base detections used for beta ordination ----
  # Includes:
  #   - all detections inside MPA/AOI polygons
  #   - plus detections inside any drawn polygons
  # Floating-panel filters and diversity target_gene filter still apply.
  diversity_detections_beta <- reactive({
    yr <- sel_year_chr()
    filters <- div_filters()

    pts <- species_sf_all
    if (yr != "All") {
       pts <- pts %>% dplyr::filter(as.character(year) == yr)
     }

    pts <- apply_species_filters(pts)
    pts <- apply_diversity_dropdown_filters(pts, filters)

    in_mpa <- sf::st_join(
      pts,
      all_polys_click %>% dplyr::select(site_name, site_type),
      join = sf::st_within,
      left = FALSE
    ) %>%
      dplyr::arrange(occurrenceID) %>%
      dplyr::group_by(
        occurrenceID, samp_name, scientificName, year, target_gene,
        site_name, site_type
      ) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()

    polys <- drawn_polys()

    if (is.null(polys) || nrow(polys) == 0) {
      return(in_mpa)
    }

    drawn_list <- lapply(seq_len(nrow(polys)), function(i) {
      g_i <- sf::st_geometry(polys[i, , drop = FALSE])
      lab <- polys$draw_label[i]

      inside_i <- pts[within_any(pts, g_i), , drop = FALSE]
      if (nrow(inside_i) == 0) return(NULL)

      inside_i %>%
        dplyr::mutate(
          site_name = lab,
          site_type = "User"
        )
    })

    in_drawn <- dplyr::bind_rows(drawn_list)

    dplyr::bind_rows(in_mpa, in_drawn)
  })

  # ---- unique sample x taxon matrix for ONE shared ordination ----
  beta_comm_mat <- reactive({
    det <- diversity_detections_beta()
    req(det)

    occ_all <- det %>%
      sf::st_drop_geometry()

    shiny::validate(
      shiny::need(nrow(occ_all) > 0, "No detections available for the current filters.")
    )

    rank_col <- div_filters()$tax_rank

    shiny::validate(
      shiny::need("organismQuantity" %in% names(occ_all),
                  "organismQuantity column is missing.")
    )

    # IMPORTANT:
    # Build matrix by UNIQUE sample only, not by group.
    # If a sample is in both an MPA and a drawn polygon,
    # it still gets one row in the ordination matrix.
    occ_all2 <- occ_all %>%
      dplyr::mutate(
        samp_name = as.character(samp_name),
        taxon     = as.character(.data[[rank_col]]),
        value     = as.numeric(organismQuantity)
      ) %>%
      dplyr::filter(
        !is.na(samp_name), samp_name != "",
        !is.na(taxon), taxon != "",
        !is.na(value)
      ) %>%
      dplyr::group_by(samp_name, taxon) %>%
      dplyr::summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

    shiny::validate(
      shiny::need(nrow(occ_all2) > 0, "No abundance data available after filtering.")
    )

    mat_wide <- occ_all2 %>%
      tidyr::pivot_wider(
        names_from  = taxon,
        values_from = value,
        values_fill = 0
      )

    mat <- mat_wide %>%
      dplyr::select(-samp_name) %>%
      as.data.frame()

    rownames(mat) <- mat_wide$samp_name
    mat
  })


  # ---- metadata for plotting: samples can belong to multiple groups ----
  beta_plot_meta <- reactive({
    det <- diversity_detections_beta()
    req(det)

    det %>%
      sf::st_drop_geometry() %>%
      dplyr::mutate(
        samp_name  = as.character(samp_name),
        site_name  = as.character(site_name),
        site_type  = as.character(site_type),
        year       = if ("year" %in% names(.)) as.character(year) else NA_character_,
        group_label = as.character(site_name)
      ) %>%
      dplyr::filter(!is.na(samp_name), samp_name != "",
                    !is.na(group_label), group_label != "") %>%
      dplyr::distinct(samp_name, group_label, site_type, year)
  })

  # --- sample metadata for grouping/hover (Location etc.) ---
  sample_meta <- reactive({
    det <- selected_detections()
    req(det)

    occ_all <- det %>% sf::st_drop_geometry()

    out <- occ_all %>%
      dplyr::mutate(
        samp_name = as.character(samp_name),
        site_name  = if ("site_name" %in% names(.)) as.character(site_name) else NA_character_,
        year      = if ("year" %in% names(.)) as.character(year) else NA_character_
      ) %>%
      dplyr::distinct(samp_name, site_name, year)

    out
  })

  #Alpha diversity
  alpha_metric_vec <- reactive({
    req(comm_mat_mpa())
    mat <- comm_mat_mpa()

    metric <- input$alpha_metric %||% "observed"

    vals <- switch(
      metric,
      observed   = vegan::specnumber(mat),
      shannon    = vegan::diversity(mat, index = "shannon"),
      simpson    = vegan::diversity(mat, index = "simpson"),
      invsimpson = vegan::diversity(mat, index = "invsimpson"),
      ace        = vegan::estimateR(mat)["S.ACE", ],
      vegan::specnumber(mat)
    )

    data.frame(
      sample_group = names(vals),
      alpha_val    = as.numeric(vals),
      stringsAsFactors = FALSE
    )
  })

  # --- Alpha boxplot data: selected area + optional drawn polygons ---
  alpha_boxplot_occ_all <- reactive({
    yr <- sel_year_chr()

    alpha <- alpha_metric_vec()
    meta  <- sample_meta_mpa()

    alpha_mpa <- alpha %>%
      dplyr::left_join(meta, by = "sample_group") %>%
      dplyr::filter(!is.na(site_name), site_name != "") %>%
      dplyr::mutate(group_label = as.character(site_name))

    # ---- add drawn polygons as extra groups ----
    polys <- drawn_polys()

    if (!is.null(polys) && nrow(polys) > 0) {

      pts <- species_sf_all
      if (yr != "All") {
      pts <- pts %>% dplyr::filter(as.character(year) == yr)
      }

      pts <- apply_species_filters(pts)
      pts <- apply_diversity_dropdown_filters(pts, div_filters())

      rank_col <- div_filters()$tax_rank
      metric   <- input$alpha_metric %||% "observed"

      alpha_poly_list <- lapply(seq_len(nrow(polys)), function(i) {
        g_i <- sf::st_geometry(polys[i, , drop = FALSE])

        inside <- pts[within_any(pts, g_i), , drop = FALSE] %>%
          sf::st_drop_geometry() %>%
          dplyr::mutate(
            samp_name = as.character(samp_name),
            taxon     = as.character(.data[[rank_col]]),
            value     = as.numeric(organismQuantity)
          ) %>%
          dplyr::filter(
            !is.na(samp_name), samp_name != "",
            !is.na(taxon), taxon != "",
            !is.na(value)
          ) %>%
          dplyr::group_by(samp_name, taxon) %>%
          dplyr::summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

        if (nrow(inside) == 0) return(NULL)

        mat_poly <- inside %>%
          tidyr::pivot_wider(
            names_from  = taxon,
            values_from = value,
            values_fill = 0
          ) %>%
          as.data.frame()

        rownames(mat_poly) <- mat_poly$samp_name
        mat_poly$samp_name <- NULL

        vals_poly <- switch(
          metric,
          observed   = vegan::specnumber(mat_poly),
          shannon    = vegan::diversity(mat_poly, index = "shannon"),
          simpson    = vegan::diversity(mat_poly, index = "simpson"),
          invsimpson = vegan::diversity(mat_poly, index = "invsimpson"),
          ace        = vegan::estimateR(mat_poly)["S.ACE", ],
          vegan::specnumber(mat_poly)
        )

        poly_label <- polys$draw_label[i]

        data.frame(
          sample_group = names(vals_poly),
          alpha_val    = as.numeric(vals_poly),
          samp_name    = names(vals_poly),
          site_name    = poly_label,
          site_type    = "User",
          year         = NA_character_,
          group_label  = poly_label,
          stringsAsFactors = FALSE
        )
      })

      alpha_poly <- dplyr::bind_rows(alpha_poly_list)

      if (nrow(alpha_poly) > 0) {
        alpha_mpa <- dplyr::bind_rows(alpha_mpa, alpha_poly)
      }
    }

    base_lvls <- unique(alpha_mpa$group_label[alpha_mpa$site_type != "User"])
    draw_lvls <- unique(alpha_mpa$group_label[alpha_mpa$site_type == "User"])

    alpha_mpa %>%
      dplyr::mutate(group_label = factor(group_label, levels = c(base_lvls, draw_lvls)))
  })

  # --- Alpha diversity boxplot (Plotly) ---
  color_vec <- c("#046c9a", "#5BBCD6", "#ABDDDE", "#446455", "#00A08A","#Fdd262")

  output$alpha_boxplot <- plotly::renderPlotly({
    alpha <- alpha_boxplot_occ_all()

    shiny::validate(
      shiny::need(nrow(alpha) > 0, "No samples available for the current selection/year.")
    )

    metric_label <- switch(
      input$alpha_metric %||% "observed",
      observed   = "Richness (count of taxa)",
      shannon    = "Shannon diversity",
      simpson    = "Simpson diversity",
      invsimpson = "Inverse Simpson diversity",
      ace        = "ACE estimated richness",
      "Alpha diversity"
    )

    rank_label <- tools::toTitleCase(gsub("_", " ", div_filters()$tax_rank))
    gene_label <- div_filters()$target_gene
    primer_label <- div_filters()$primers

    gene_suffix <- if (length(gene_label) > 0) {
      paste0(" | Genes: ", paste(gene_label, collapse = ", "))
    } else {
      ""
    }

    primer_suffix <- if (length(primer_label) > 0) {
      paste0(" | Primers: ", paste(primer_label, collapse = "; "))
    } else {
      ""
    }

    metric_label <- paste0(
      metric_label, " at ", rank_label, " level",
      gene_suffix, primer_suffix
    )


    plotly::plot_ly(
      data = alpha,
      x = ~group_label,
      y = ~alpha_val,
      type = "box",
      color = ~group_label,
      colors = color_vec,
      boxpoints = "all",
      jitter = 0.3,
      pointpos = 0,
      hovertemplate = paste(
        "<b>Sample:</b> %{customdata[0]}<br>",
        "<b>Group:</b> %{x}<br>",
        "<b>", metric_label, ":</b> %{y}<extra></extra>"
      ),
      customdata = ~cbind(samp_name)
    ) %>%
      plotly::layout(
        font = list(size = 18),
        xaxis = list(
          title = list(text = "Location (Polygon)", font = list(size = 20)),
          tickfont = list(size = 16),
          showgrid = FALSE,
          showline = TRUE,
          linecolor = "black",
          rangemode = "tozero"
        ),
        yaxis = list(
          title = list(text = metric_label, font = list(size = 20)),
          tickfont = list(size = 16),
          showgrid = FALSE,
          showline = TRUE,
          linecolor = "black",
          range = c(0, max(alpha$alpha_val, na.rm = TRUE) * 1.05),
          fixedrange = FALSE
        ),
        margin = list(l = 80, r = 30, t = 40, b = 120),
        showlegend = FALSE
      )
  })

  #Beta diversity
  output$beta_pcoa <- plotly::renderPlotly({
    mat  <- beta_comm_mat()
    meta <- beta_plot_meta()

    shiny::validate(
      shiny::need(nrow(mat) > 2, "At least 3 samples are required to compute a PCoA.")
    )


    beta_method <- input$beta_metric %||% "bray"
    mat_use <- as.matrix(mat)

    if (beta_method %in% c("bray", "jaccard")) {
      rs <- rowSums(mat_use, na.rm = TRUE)
      mat_rel <- mat_use
      mat_rel[rs > 0, ] <- mat_use[rs > 0, , drop = FALSE] / rs[rs > 0]

      if (beta_method == "jaccard") {
        mat_dist_input <- (mat_rel > 0) * 1
      } else {
        mat_dist_input <- mat_rel
      }

      d <- vegan::vegdist(mat_dist_input, method = beta_method)

    } else if (beta_method == "euclidean") {
      d <- stats::dist(mat_use, method = "euclidean")

    } else if (beta_method == "robust.aitchison") {
      shiny::validate(
        shiny::need(all(mat_use >= 0, na.rm = TRUE),
                    "Robust Aitchison distance requires non-negative values.")
      )

      d <- vegan::vegdist(mat_use, method = "robust.aitchison")
    }


    ord <- stats::cmdscale(d, k = 2, eig = TRUE)

    sample_ids <- rownames(mat_use)

    shiny::validate(
      shiny::need(!is.null(sample_ids), "Sample names are missing from the beta diversity matrix."),
      shiny::need(length(sample_ids) == nrow(ord$points), "Mismatch between sample names and ordination points.")
    )

    scores <- data.frame(
      samp_name = sample_ids,
      PC1 = ord$points[, 1],
      PC2 = ord$points[, 2],
      stringsAsFactors = FALSE
    )

    plot_df <- scores %>%
      dplyr::left_join(meta, by = "samp_name") %>%
      dplyr::mutate(
        group_label = as.character(group_label),
        site_type   = as.character(site_type)
      ) %>%
      dplyr::filter(!is.na(group_label), nzchar(group_label)) %>%
      dplyr::distinct(samp_name, group_label, site_type, .keep_all = TRUE)

    shiny::validate(
      shiny::need(nrow(plot_df) > 0, "No plotting metadata available for the current filters.")
    )

    eig <- ord$eig
    eig_pos <- eig[eig > 0]
    pct1 <- if (length(eig_pos) >= 1) round(100 * eig_pos[1] / sum(eig_pos), 1) else NA_real_
    pct2 <- if (length(eig_pos) >= 2) round(100 * eig_pos[2] / sum(eig_pos), 1) else NA_real_

    method_label <- switch(
      beta_method,
      bray             = "Bray-Curtis",
      jaccard          = "Jaccard",
      euclidean        = "Euclidean",
      `robust.aitchison` = "Robust Aitchison",
      beta_method
    )

    rank_label <- tools::toTitleCase(gsub("_", " ", div_filters()$tax_rank))
    gene_label <- div_filters()$target_gene
    primer_label <- div_filters()$primers

    gene_suffix <- if (length(gene_label) > 0) {
      paste0(" | Genes: ", paste(gene_label, collapse = ", "))
    } else {
      ""
    }

    primer_suffix <- if (length(primer_label) > 0) {
      paste0(" | Primers: ", paste(primer_label, collapse = "; "))
    } else {
      ""
    }

    plot_title <- paste0(
      "Beta diversity (PCoA; ", method_label, ") at ",
      rank_label, " level",
      gene_suffix, primer_suffix
    )

    # Build hover text explicitly so it always matches plot_df row count
    plot_df <- plot_df %>%
      dplyr::mutate(
        hover_text = paste0(
          "<b>Sample:</b> ", samp_name,
          "<br><b>Group:</b> ", group_label,
          "<br><b>Type:</b> ", site_type,
          "<br><b>PC1:</b> ", sprintf("%.3f", PC1),
          "<br><b>PC2:</b> ", sprintf("%.3f", PC2)
        )
      )

    plotly::plot_ly(
      data = plot_df,
      x = ~PC1,
      y = ~PC2,
      type = "scatter",
      mode = "markers",
      color = ~group_label,
      colors = color_vec,
      text = ~hover_text,
      hoverinfo = "text"
    ) %>%
      plotly::layout(
        title = list(
          text = plot_title,
          font = list(size = 20)
        ),
        font = list(size = 18),
        xaxis = list(
          title = list(
            text = paste0("PC1", if (!is.na(pct1)) paste0(" (", pct1, "%)") else ""),
            font = list(size = 20)
          ),
          tickfont = list(size = 16),
          showgrid = FALSE,
          showline = TRUE,
          linecolor = "black"
        ),
        yaxis = list(
          title = list(
            text = paste0("PC2", if (!is.na(pct2)) paste0(" (", pct2, "%)") else ""),
            font = list(size = 20)
          ),
          tickfont = list(size = 16),
          showgrid = FALSE,
          showline = TRUE,
          linecolor = "black"
        ),
        legend = list(
          font = list(size = 16)
        ),
        margin = list(l = 80, r = 30, t = 70, b = 80),
        showlegend = TRUE
      )
  })


  #Krona plot
  output$tax_krona <- taxplore::renderKronaChart({
    det <- selected_detections()
    req(det)

    # floating-panel filters
    det <- apply_species_filters(det)

    # data-selection filters: target gene + primer only
    det <- apply_diversity_dropdown_filters(
      det,
      list(
        target_gene = div_filters()$target_gene,
        primers     = div_filters()$primers
      )
    )

    shiny::validate(
      shiny::need(!is.null(det) && nrow(det) > 0,
                  "Select a cell/polygon (or draw a polygon) to display a Krona chart.")
    )

    det0 <- det %>%
      sf::st_drop_geometry()

    shiny::validate(
      shiny::need("organismQuantity" %in% names(det0),
                  "organismQuantity column is missing.")
    )

    tax_ranks <- c("kingdom", "phylum", "class", "order", "family", "genus")

    krona_occ_all <- det0 %>%
      dplyr::mutate(
        species = as.character(scientificName),
        qty     = as.numeric(organismQuantity)
      ) %>%
      dplyr::filter(
        !is.na(species), species != "",
        !is.na(qty), qty > 0
      ) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(c(tax_ranks, "species")))) %>%
      dplyr::summarise(
        magnitude = sum(qty, na.rm = TRUE),
        .groups   = "drop"
      )

    shiny::validate(
      shiny::need(
        nrow(krona_occ_all) > 0,
        "No detections with organismQuantity > 0 after the current target gene and primer filters."
      )
    )

    tax_occ_all <- krona_occ_all %>%
      dplyr::select(dplyr::all_of(c(tax_ranks, "species")))

    taxplore::plot_krona(
      tax_occ_all,
      magnitude   = krona_occ_all$magnitude,
      total_label = "Sum organismQuantity"
    )
  })

  #Helper to generate species vectors
  get_species_vec <- function(det_sf, gene = NULL) {
    req(det_sf)

    occ_all <- det_sf

    if (!is.null(gene)) {
      occ_all <- occ_all %>% dplyr::filter(as.character(target_gene) == as.character(gene))
    }

    occ_all <- apply_species_filters(occ_all)

    occ_all %>%
      sf::st_drop_geometry() %>%
      dplyr::pull(scientificName) %>%
      as.character() %>%
      unique() %>%
      stats::na.omit() %>%
      sort()
  }

  # ---- species panel (drawn polygon OR click) ----
  output$species_panel <- renderUI({

    active_groups()
    filter_sara_on()
    filter_ais_on()

    yr <- sel_year_chr()

    make_list_local <- function(vec, max_h = 220) {
      if (length(vec) == 0) return(em("No species detected for this layer."))
      tags$div(
        style = paste0("max-height: ", max_h, "px; overflow-y: auto; padding-left: 10px;"),
        tags$ul(lapply(vec, tags$li))
      )
    }

    render_by_layer_ui <- function(by_layer, header = NULL, max_h = 220) {
      if (is.null(by_layer) || length(by_layer) == 0) {
        return(em("Turn ON a richness layer (All / 12S / COI / 16S / 18S) to view species lists."))
      }

      ord <- c("All", "12S", "COI", "16S", "18S")
      nm  <- intersect(ord, names(by_layer))
      if (length(nm) == 0) nm <- names(by_layer)

      panels <- tagList()
      for (k in nm) {
        vec <- by_layer[[k]]
        panels <- tagAppendChildren(
          panels,
          tags$details(
            open = TRUE,
            tags$summary(strong(paste0(k, " (", length(vec), ")"))),
            make_list_local(vec, max_h = max_h),
            tags$hr()
          )
        )
      }

      if (!is.null(header)) tagList(header, tags$hr(), panels) else panels
    }

    layers_on <- active_richness_layers()
    if (length(layers_on) == 0) {
      return(em("Turn ON a richness layer (All / 12S / COI / 16S / 18S) to view species lists."))
    }

    sel_id <- selected_draw_id()
    click  <- input$map_shape_click

    # drawn polygon branch
    if (!is.null(sel_id)) {
      polys <- drawn_polys()
      hit <- polys %>% dplyr::filter(draw_id == sel_id)
      if (nrow(hit) == 0) return(em("No species detected for the current selection."))

      pts <- species_sf_by_year[[yr]]
      if (is.null(pts)) pts <- species_sf_min[0, ]

      inside <- get_inside_cached(
        selected_key = paste0("draw:", sel_id),
        pts          = pts,
        geom         = sf::st_geometry(hit),
        year         = yr
      )

      by_layer_all <- get_species_by_layer_cached(
        selected_key = paste0("draw:", sel_id),
        det_sf       = inside,
        year         = yr,
        groups       = active_groups(),
        sara_on      = isTRUE(filter_sara_on()),
        ais_on       = isTRUE(filter_ais_on())
      )

      by_layer <- by_layer_all[intersect(names(by_layer_all), layers_on)]
      poly_lab <- hit$draw_label[1] %||% sel_id

      return(render_by_layer_ui(by_layer, strong(paste0("Drawn polygon: ", poly_lab)), 180))
    }

    if (is.null(click) || is.null(click$id)) {
      return(em("Click a grid cell or an MPA/AOI outline, or draw a polygon."))
    }

    if (grepl("\\|\\|", click$id)) {
      groups_on <- input$map_groups %||% character(0)
      show_poly_total <- "Total species detected per MPA/AOI" %in% groups_on
      if (!show_poly_total) {
        return(em("Turn ON “Total species detected per MPA/AOI” to view polygon species lists."))
      }

      parts  <- strsplit(click$id, "\\|\\|")[[1]]
      p_type <- parts[1]
      p_name <- parts[2]

      poly_sel <- all_polys_click %>%
        dplyr::filter(site_type == p_type, site_name == p_name)

      if (nrow(poly_sel) == 0) return(em("Polygon not found."))

      pts <- species_sf_by_year[[yr]]
      if (is.null(pts)) pts <- species_sf_min[0, ]

      inside <- get_inside_cached(
        selected_key = paste0("click:", click$id),
        pts          = pts,
        geom         = sf::st_geometry(poly_sel),
        year         = yr
      )

      by_layer_all <- get_species_by_layer_cached(
        selected_key = paste0("click:", click$id),
        det_sf       = inside,
        year         = yr,
        groups       = active_groups(),
        sara_on      = isTRUE(filter_sara_on()),
        ais_on       = isTRUE(filter_ais_on())
      )

      by_layer <- by_layer_all[intersect(names(by_layer_all), layers_on)]
      return(render_by_layer_ui(by_layer, strong(paste0(p_type, ": ", p_name)), 220))
    }

    cid <- suppressWarnings(as.integer(click$id))
    if (is.na(cid)) {
      return(em("Click a grid cell or an MPA/AOI outline."))
    }

    pts <- species_sf_by_year[[yr]]
    if (is.null(pts)) pts <- species_sf_min[0, ]

    cell_poly <- grid_clip %>% dplyr::filter(cell_id == cid)

    det <- get_inside_cached(
      selected_key = paste0("click:", click$id),
      pts          = pts,
      geom         = sf::st_geometry(cell_poly),
      year         = yr
    )

    if (is.null(det) || nrow(det) == 0) {
      return(em("No species detected for the current selection."))
    }

    by_layer_all <- get_species_by_layer_cached(
      selected_key = paste0("click:", click$id),
      det_sf       = det,
      year         = yr,
      groups       = active_groups(),
      sara_on      = isTRUE(filter_sara_on()),
      ais_on       = isTRUE(filter_ais_on())
    )

    by_layer <- by_layer_all[intersect(names(by_layer_all), layers_on)]
    render_by_layer_ui(by_layer, strong(paste0("Grid cell: ", cid)), 220)
  })

  output$detections_tbl <- DT::renderDT({
    det <- selected_detections()

    if (is.null(det) || nrow(det) == 0) {
      return(DT::datatable(
        data.frame(Message = "Click an MPA/AOI polygon/grid cell, or draw a polygon, to view detection details"),
        rownames = FALSE,
        options = list(dom = "t")
      ))
    }

    det_f <- apply_species_filters(det)

    if (nrow(det_f) == 0) {
      return(DT::datatable(
        data.frame(Message = "No detections match the current Year and/or Group filters."),
        rownames = FALSE,
        options = list(dom = "t")
      ))
    }

    det_show <- det_f %>%
      dplyr::mutate(
        decimalLongitude = sf::st_coordinates(.)[, 1],
        decimalLatitude  = sf::st_coordinates(.)[, 2]
      ) %>%
      sf::st_drop_geometry()

    needed_cols <- c(
      "kingdom", "phylum", "class", "order", "family", "genus",
      "scientificName", "samp_name", "target_gene", "eventDate",
      "pcr_primer_name_forward", "pcr_primer_name_reverse",
      "pcr_primer_forward", "pcr_primer_reverse",
      "organismQuantity", "project_contact", "LClabel",
      "samp_size", "samp_size_unit",
      "minimumDepthInMeters", "maximumDepthInMeters",
      "bathymetry",
      "decimalLatitude", "decimalLongitude",
      "flags", "dataset_id"
    )

    missing_cols <- setdiff(needed_cols, names(det_show))
    for (nm in missing_cols) {
      det_show[[nm]] <- NA_character_
    }

    fwd_col <- get_forward_primer_col(det_show)
    rev_col <- get_reverse_primer_col(det_show)

    fwd_vals <- if (!is.null(fwd_col)) as.character(det_show[[fwd_col]]) else rep(NA_character_, nrow(det_show))
    rev_vals <- if (!is.null(rev_col)) as.character(det_show[[rev_col]]) else rep(NA_character_, nrow(det_show))

    fwd_vals <- trimws(fwd_vals)
    rev_vals <- trimws(rev_vals)
    fwd_vals[fwd_vals == ""] <- NA_character_
    rev_vals[rev_vals == ""] <- NA_character_

    det_show$Primers <- dplyr::case_when(
      !is.na(fwd_vals) & !is.na(rev_vals) ~ paste(fwd_vals, rev_vals, sep = " | "),
      !is.na(fwd_vals) ~ fwd_vals,
      !is.na(rev_vals) ~ rev_vals,
      TRUE ~ NA_character_
    )

    det_show <- det_show %>%
      dplyr::mutate(
        across(c(
          scientificName, samp_name, target_gene,
          bathymetry, flags, dataset_id, eventDate
        ), ~as.character(.)),
        Volume = dplyr::case_when(
          !is.na(samp_size) & !is.na(samp_size_unit) ~ paste(samp_size, samp_size_unit),
          !is.na(samp_size) ~ as.character(samp_size),
          !is.na(samp_size_unit) ~ as.character(samp_size_unit),
          TRUE ~ NA_character_
        ),
        `Depth Range (m)` = dplyr::case_when(
          !is.na(minimumDepthInMeters) & !is.na(maximumDepthInMeters) ~
            paste(minimumDepthInMeters, maximumDepthInMeters, sep = " | "),
          !is.na(minimumDepthInMeters) ~ as.character(minimumDepthInMeters),
          !is.na(maximumDepthInMeters) ~ as.character(maximumDepthInMeters),
          TRUE ~ NA_character_
        )
      ) %>%
      dplyr::arrange(scientificName, samp_name) %>%
      dplyr::select(
        kingdom,
        phylum,
        class,
        order,
        family,
        genus,
        scientificName,
        samp_name,
        eventDate,
        target_gene,
        Primers,
        organismQuantity,
        Volume,
        `Depth Range (m)`,
        bathymetry,
        decimalLatitude,
        decimalLongitude,
        flags,
        dataset_id,
        project_contact,
        LClabel
      )

    DT::datatable(
      det_show,
      rownames = FALSE,
      colnames = c(
        "Kingdom",
        "Phylum",
        "Class",
        "Order",
        "Family",
        "Genus",
        "Species",
        "Sample Name",
        "Date",
        "Target Gene",
        "Primers",
        "# of Sequence Reads",
        "Volume",
        "Minimum/Maximum Depth (m)",
        "Bathymetry",
        "Latitude",
        "Longitude",
        "OBIS Data Flags",
        "OBIS Dataset Identifier",
        "Project Contact",
        "Indigenous Acknowledgement & Contributions"
      ),
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        autoWidth = TRUE
      ),
      class = "nowrap"
    )
  })

  # highlight selected cell (outline on top)
  observeEvent(input$map_shape_click, {
    click <- input$map_shape_click
    req(click, click$id)

    id <- as.character(click$id)

    polys <- drawn_polys()

    # clicked one of the stored drawn polygons
    if (nrow(polys) > 0 && id %in% polys$draw_id) {
      selected_draw_id(id)
      session$sendCustomMessage("openFloating", list(id = id))
      return()
    }

    # otherwise switch away from drawn polygons
    selected_draw_id(NULL)

    is_mpa  <- grepl("\\|\\|", id)
    is_cell <- !is.na(suppressWarnings(as.integer(id)))

    if (is_mpa || is_cell) {
      session$sendCustomMessage("openFloating", list(id = id))
    }
  }, ignoreInit = TRUE)

}

shinyApp(ui, server)

###Updates from meeting

#edit organism quantity for read_data

###Wishlist items
#add polygons for other marine conservation regions in the Atlantic
#add NMDS plot for community structure - see code Nick provides - year, season, depth


#draft skeleton of paper with a sentence or two of why each was done
#Look into journals for publishing function manuscript (maybe ecological informatics)
