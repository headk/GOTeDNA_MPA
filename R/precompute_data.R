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
    dataset_ids       = "634a3b54-9e8e-4b20-80fa-23eaed15284d",
    scientificname    = NULL,
    worms_id          = NULL,
    areaid            = NULL,
    join_by           = c("auto", "occurrenceID", "id"),
    require_absences  = TRUE
) {

  join_by <- match.arg(join_by)

  # ---- columns you want back ----
  occurrence_cols <- c(
    "recordedBy","bibliographicCitation","materialSampleID",
    "organismQuantity","organismQuantityType",
    "sampleSizeValue","sampleSizeUnit","associatedSequences",
    "minimumDepthInMeters","maximumDepthInMeters","month","year",
    "scientificNameID","kingdom","phylum","class","order","family","genus"
  )

  dna_cols <- c(
    "id","dna_sequence","target_gene","pcr_primer_forward","samp_name",
    "env_broad_scale","env_local_scale","env_medium","samp_mat_process",
    "size_frac","samp_size","samp_size_unit","otu_db","seq_kit",
    "otu_seq_comp_appr","pcr_primer_name_forward","pcr_primer_name_reverse",
    "pcr_primer_reference","occurrenceID"
  )

  mof_cols <- c(
    "id","seq_id","samp_category","checkls_ver","assay_name","assay_type",
    "targetTaxonomicAssay","geo_loc_name","technical_rep_id","project_contact",
    "seq_run_id","lib_id","project_id","pcr_0_1","samp_store_sol","samp_store_temp",
    "platform","instrument","tax_assign_cat","LClabel","occurrenceID",
    "nucl_acid_ext","nucl_acid_ext_kit","filter_material"
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
  enforce_cols <- function(df, cols) {
    # add missing columns as NA
    missing <- setdiff(cols, names(df))
    if (length(missing) > 0) {
      for (m in missing) df[[m]] <- NA
    }
    # keep only requested columns in a consistent order
    df <- df[, intersect(cols, names(df)), drop = FALSE]
    df
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

  GOTeDNA_df <- dplyr::bind_rows(obis_list)
  rownames(GOTeDNA_df) <- NULL
  GOTeDNA_df
}

STORED_DATA <- read_data_new()
saveRDS(STORED_DATA, "./data/test_read_data_file.rds")
STORED_DATA <- readRDS("./data/test_read_data_file.rds")

################DO NOT CHANGE ABOVE CODE



#Connect code below to stored data object

occ_all <- STORED_DATA %>%
  mutate(
    year            = as.character(year),
    occurrenceStatus = tolower(as.character(occurrenceStatus)),
    decimalLatitude  = suppressWarnings(as.numeric(decimalLatitude)),
    decimalLongitude = suppressWarnings(as.numeric(decimalLongitude)),
    target_gene      = case_when(
      str_detect(tolower(target_gene), "12s") ~ "12S",
      str_detect(tolower(target_gene), "coi") ~ "COI",
      str_detect(tolower(target_gene), "16s") ~ "16S",
      str_detect(tolower(target_gene), "18s") ~ "18S",
      TRUE ~ as.character(target_gene)
    )
  )

# ---- Build available gene-year keys dynamically ----
KEY_TBL <- occ_all %>%
  filter(!is.na(year), year != "", !is.na(target_gene), target_gene != "") %>%
  distinct(target_gene, year) %>%
  mutate(
    target_gene = as.character(target_gene),
    year        = as.character(year),
    key         = paste(target_gene, year, sep = "_")
  ) %>%
  arrange(target_gene, year)

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
points_sf_from_df <- function(df) {
  df %>%
    filter(!is.na(decimalLatitude), !is.na(decimalLongitude),
           tolower(as.character(occurrenceStatus)) == "present") %>%
    st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
}

SPECIES_SF_BY_KEY <- purrr::imap(DATA_BY_KEY, ~{
  parts <- strsplit(.y, "_")[[1]]
  gene <- parts[1]
  yr   <- parts[2]

  points_sf_from_df(.x) %>%
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
SARA <- read.xlsx(file.path("data", "sara_ais", "SARA_Clean_Schedule1_specieslist.xlsx"))
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


sample_tag <- function(df) {
  df %>%
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


get_depth_num <- function(df) {
  depth_col <- dplyr::case_when(
    "minimumDepthInMeters" %in% names(df) ~ "minimumDepthInMeters",
    "Depth" %in% names(df) ~ "Depth",
    TRUE ~ NA_character_
  )
  if (is.na(depth_col)) return(df %>% mutate(depth_m = NA_real_))

  df %>% mutate(depth_m = suppressWarnings(as.numeric(.data[[depth_col]])))
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

  sf::st_intersection(out, all_polys) %>%
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

depth_layers <- lapply(depth_layers, function(df) {
  df %>%
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

#Bundle the outputs in one list:
APP_DATA <- list(
  occ_all = occ_all,
  KEY_TBL = KEY_TBL,
  sampling_pts = sampling_pts,
  species_sf_all = species_sf_all,
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

/* === Shared legend typography === */
.leaflet-control.legend-base {
  font-family: Arial !important;
  font-size: 15px !important;
  font-weight: 400 !important;
  line-height: 1.45 !important;
  color: #222 !important;
}

/* === Shared legend swatches === */
.leaflet-control.legend-base i{
  width: 18px !important;
  height: 18px !important;
  display: inline-block !important;
  margin-right: 8px !important;
  vertical-align: middle !important;
  float: none !important;
}

.leaflet-control.legend-richness-box{
  background: rgba(255,255,255,0.92) !important;
  border-radius: 6px !important;
  box-shadow: 0 2px 10px rgba(0,0,0,0.25) !important;
  border: 1px solid rgba(0,0,0,0.15) !important;
  padding: 14px 16px !important;
  width: 190px !important;
}


/* === protocol cards (GOTeDNA palette) === */
.protocol-grid {
  display: grid;
  grid-template-columns: repeat(auto-fill, minmax(140px, 1fr));
  gap: 12px;
}

.protocol-field-title {
  font-size: 13px;
  font-weight: 400 !important;
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
/* --- Leaflet legend container --- */
.leaflet-control.legend-richness{
  background: rgba(255,255,255,0.92) !important;
  border-radius: 6px !important;
  box-shadow: 0 2px 10px rgba(0,0,0,0.25) !important;
  border: 1px solid rgba(0,0,0,0.15) !important;
  padding: 14px 16px !important;
  width: 190px !important;
  font-size: 16px !important;
  line-height: 1.4 !important;
}

/* --- Swatches (this is the important part) --- */
.leaflet-control.legend-richness i{
  width: 18px !important;
  height: 18px !important;
  display: inline-block !important;
  margin-right: 8px !important;
  vertical-align: middle !important;
  float: none !important;            /* avoid float weirdness */
}

/* Each label is text after the swatch and before <br> */
.leaflet-control.legend-richness br{
  line-height: 22px !important;
}
/* Headings inside the layer control */
.leaflet-control-layers .layers-heading{
  font-weight: 700;
  margin: 6px 0 4px 0;
  font-size: 13px;
  opacity: 0.9;
}

/* A separator line */
.leaflet-control-layers .layers-sep{
  border-top: 1px solid rgba(0,0,0,0.15);
  margin: 6px 0;
}

.legend-hidden{ display:none !important; }

/* Depth legend box */
.leaflet-control.legend-depth-box{
  background: rgba(255,255,255,0.92) !important;
  border-radius: 6px !important;
  box-shadow: 0 2px 10px rgba(0,0,0,0.25) !important;
  border: 1px solid rgba(0,0,0,0.15) !important;
  padding: 14px 16px !important;
  width: 190px !important;
}

/* Depth legend color bins */
.leaflet-control.legend-depth-box i{
  width: 18px !important;
  height: 18px !important;
  display: inline-block !important;
  margin-right: 8px !important;
  vertical-align: middle !important;
  float: none !important;
}

.protocol-group-title{
  font-weight: 700 !important;
  font-size: 15px;
  margin: 14px 0 8px 0;
  color: #0b2a2a;
}

/* remove the extra top offset */
.protocol-details-wrap { margin-top: 0 !important; }

/* only the FIRST section title gets no top margin */
.protocol-details-wrap .protocol-group-title:first-of-type{
  margin-top: 0 !important;
}

/* normal spacing for every section title */
.protocol-group-title{
  font-weight: 700;
  margin: 14px 0 8px 0;   /* <-- this is the separation you lost */
}

/* but DON'T push the first title down */
.protocol-details-wrap .protocol-group-title:first-of-type{
  margin-top: 0;
}

.protocol-grid { margin-bottom: 10px; }

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
          selectInput("req_location", "Location",                                         #Remove later once linked to OBIS
                      choices = character(0), selected = NULL, selectize = FALSE),        #Remove later once linked to OBIS
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

    div(
      id = "sec_div", class = "scroll-section",
      h3("Diversity Metrics"),

      #tags$p(tags$strong("Alpha diversity (boxplot): "), "I want to include the following alpha metrics as a dropdown: Observed, Chao1, Shannon, Simpson, Fisher, InvSimpson, ACE"),   #include if we want written discriptions as part of this section
      #tags$p(tags$strong("Beta diversity (PCoA): "), "I want to include the following beta metrics as a dropdown: Bray, Jaccard, Euclidean, Aitchison"),

      fluidRow(
        column(
        width = 2,
        selectInput(
          "tax_rank",
          "Taxonomic Rank",
          choices = c(
            "Species"      = "scientificName",
            "Genus"        = "genus",
            "Family"       = "family",
            "Order"        = "order",
            "Class"        = "class",
            "Phylum"       = "phylum",
            "Kingdom"      = "kingdom"
          ),
          selected = "scientificName"
        )
      ),

    column(
      width = 2,
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
     )
    ),

    # ---- Alpha plot (full width row) ----
    fluidRow(
      column(
        width = 10,
        offset = 1,
        plotly::plotlyOutput("alpha_boxplot", height = "720px")
      )
    ),

    # ---- Beta plot (full width row) ----
    fluidRow(
      column(
        width = 10,
        offset = 1,
        plotly::plotlyOutput("beta_pcoa", height = "720px")
      )
    )
  ),

  div(
    id="sec_pie", class="scroll-section",
    h3("Taxonomic Pie Chart"),

    taxplore::KronaChartOutput("tax_krona", height = "670px"))
)


## ---- 3) Server ----

server <- function(input, output, session){
  
  `%||%` <- function(x, y) if (is.null(x)) y else x
  
  dat <- APP_DATA

  occ_all          <- dat$occ_all
  KEY_TBL          <- dat$KEY_TBL
  sampling_pts     <- dat$sampling_pts
  species_sf_all   <- dat$species_sf_all
  grid_clip        <- dat$grid_clip
  RICHNESS_BY_KEY  <- dat$RICHNESS_BY_KEY
  RICHNESS_ALL     <- dat$RICHNESS_ALL
  RICHNESS_ALL_BY_YEAR <- dat$RICHNESS_ALL_BY_YEAR
  depth_layers     <- dat$depth_layers
  all_polys_click  <- dat$all_polys_click
  all_polys_zones  <- dat$all_polys_zones
  pal_rich         <- dat$pal_rich

  # ---- Data Request: use occ_all directly ----

  meta_all <- reactive({
    req(occ_all)
    df <- occ_all

    df %>%
      dplyr::mutate(
        Location   = if ("Location" %in% names(.)) as.character(Location) else NA_character_,
        ProtocolID = if ("ProtocolID" %in% names(.)) as.character(ProtocolID) else NA_character_
      )
  })

  # ---- Dynamic Year dropdown from available data ----

  # ---- Year dropdown updates to only years present in current selection ----
  observeEvent(
    list(selection_geom(), species_sf_all),
    {
      pts <- species_sf_all
      shiny::req(pts)

      g <- selection_geom()

      # If there is a selection geometry, filter points inside it; else use all points
      if (!is.null(g)) {
        inside <- pts[sf::st_within(pts, g, sparse = FALSE), , drop = FALSE]
      } else {
        inside <- pts
      }

      yrs <- inside %>%
        sf::st_drop_geometry() %>%
        dplyr::pull(year) %>%
        as.character() %>%
        unique() %>%
        na.omit() %>%
        sort()

      # Keep current selection if still valid; otherwise reset to "All"
      cur <- input$sel_year %||% "All"
      cur <- as.character(cur)
      new_choices <- c("All", yrs)

      new_selected <- if (cur %in% new_choices) cur else "All"

      updateSelectInput(
        session,
        "sel_year",
        choices  = new_choices,
        selected = new_selected
      )
    },
    ignoreInit = FALSE
  )

  # Leaflet polygon click input is typically input$map_shape_click
  observeEvent(input$map_shape_click, {
    req(input$map_shape_click)

    # whenever a user selects an MPA/AOI (polygon click), open the floating panel
    session$sendCustomMessage("openFloating", list())
  })

  # your existing outputs here:
  output$cell_summary   <- renderUI({ tags$div("...") })


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

    method_groups <- list(
      "Field Methods" = c("samp_size","size_frac","filter_material","samp_mat_process",
                          "minimumDepthInMeters","maximumDepthInMeters"),
      "Storage Methods" = c("samp_store_temp","samp_store_sol"),
      "Lab Methods" = c("target_gene","pcr_primer_forward","pcr_primer_reverse","nucl_acid_ext_kit"),
      "Library Preparation" = c("platform","instrument","seq_kit"),
      "Bioinformatic Methods" = c("otu_db","tax_assign_cat","otu_seq_comp_appr")
    )

    pick_display <- function(x) {
      x <- as.character(x)
      x <- x[!is.na(x) & nzchar(trimws(x))]
      if (length(x) == 0) return(NA_character_)
      ux <- unique(x)
      if (length(ux) == 1) ux else paste(ux, collapse = " | ")
    }

    # Build base df (for comparison)
    bp <- base_protocol()
    df_base <- NULL
    if (!is.null(bp) && nzchar(bp)) {
      df_base <- meta_all() %>%
        dplyr::filter(Location == input$req_location, ProtocolID == bp)
      if (nrow(df_base) == 0) df_base <- NULL
    }

    # UI

    tags$div(
      class = "protocol-details-wrap",
      tagList(
        lapply(names(method_groups), function(group_name) {
          fields <- method_groups[[group_name]]

          # Keep only fields that actually exist in df
          fields <- fields[fields %in% names(df)]
          if (length(fields) == 0) return(NULL)  # hide empty groups

          tags$div(
            tags$h5(class = "protocol-group-title", group_name),

            tags$div(
              class="protocol-grid",

              lapply(fields, function(f) {
                cur <- pick_display(df[[f]])
                bas <- if (!is.null(df_base) && f %in% names(df_base)) pick_display(df_base[[f]]) else NA_character_

                cur2 <- ifelse(is.na(cur), "", trimws(as.character(cur)))
                bas2 <- ifelse(is.na(bas), "", trimws(as.character(bas)))

                changed <- nzchar(cur2) && nzchar(bas2) && !identical(cur2, bas2)
                is_na   <- !nzchar(cur2)

                tags$div(
                  tags$div(class="protocol-field-title", f),
                  tags$div(
                    class = paste("protocol-card", if (changed) "changed", if (is_na) "na"),
                    if (is_na) "—" else cur2
                  )
                )
              })
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
        n_samples    = dplyr::n_distinct(samp_name),
        samples      = paste(sort(unique(samp_name)), collapse = ", "),
        #files        = paste(sort(unique(na.omit(source_file))), collapse = ", "),
        years        = paste(sort(unique(na.omit(year))), collapse = ", "),
        markers      = paste(sort(unique(na.omit(target_gene))), collapse = ", "),
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
        n_samples    = dplyr::n_distinct(samp_name),
        samples      = paste(sort(unique(samp_name)), collapse = ", "),
        #files        = paste(sort(unique(na.omit(source_file))), collapse = ", "),
        years        = paste(sort(unique(na.omit(year))), collapse = ", "),
        markers      = paste(sort(unique(na.omit(target_gene))), collapse = ", "),
        .groups = "drop"
      ) %>%
      dplyr::arrange(scientificName)

    DT::datatable(out, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE))
  })


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

  # ---- selection geometry (drawn polygon OR clicked polygon OR clicked grid cell) ----
  selection_geom <- reactive({
    # A) drawn polygon takes priority
    poly <- drawn_poly()
    if (!is.null(poly) && nrow(poly) > 0) {
      return(sf::st_geometry(poly))
    }

    # B) click-based selection
    click <- input$map_shape_click
    if (is.null(click) || is.null(click$id)) return(NULL)

    # B1) MPA/AOI polygon click (you set layerId = paste(site_type, site_name, sep="||"))
    if (grepl("\\|\\|", click$id)) {
      parts <- strsplit(click$id, "\\|\\|")[[1]]
      p_type <- parts[1]
      p_name <- parts[2]
      poly_sel <- all_polys_click %>% dplyr::filter(site_type == p_type, site_name == p_name)
      if (nrow(poly_sel) == 0) return(NULL)
      return(sf::st_geometry(poly_sel))
    }

    # B2) grid cell click
    cid <- suppressWarnings(as.integer(click$id))
    if (is.na(cid)) return(NULL)
    cell_poly <- grid_clip %>% dplyr::filter(cell_id == cid)
    if (nrow(cell_poly) == 0) return(NULL)
    sf::st_geometry(cell_poly)
  })

  observeEvent(input$map_draw_new_feature, {
    drawn_poly(feature_to_sf(input$map_draw_new_feature))
  })

  observeEvent(input$map_draw_deleted_features, {
    drawn_poly(NULL)
  })

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
      dplyr::group_by(occurrenceID, samp_name, scientificName, year, target_gene) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()
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
        options = pathOptions(pane = "pane_points"),
        label       = ~paste0(
          "Marker: ", target_gene,
          ifelse(is.na(year), "", paste0(" | Year: ", year)),
          #ifelse(is.na(source_file), "", paste0(" | File: ", source_file)),
          ifelse(is.na(samp_name), "", paste0(" | Sample: ", samp_name))
        )
      )
  }, ignoreInit = TRUE)

  # ---- 1) Render the leaflet map ONCE ----
  output$map <- renderLeaflet({

    # ---- choose initial year + initial layers safely ----
    yrs <- sort(unique(na.omit(as.character(KEY_TBL$year))))

    init_12S <- {
      k <- paste0("12S_", default_year)
      if (k %in% names(RICHNESS_BY_KEY)) RICHNESS_BY_KEY[[k]] else NULL
    }

    init_COI <- {
      k <- paste0("COI_", default_year)
      if (k %in% names(RICHNESS_BY_KEY)) RICHNESS_BY_KEY[[k]] else NULL
    }

    init_16S <- {
      k <- paste0("16S_", default_year)
      if (k %in% names(RICHNESS_BY_KEY)) RICHNESS_BY_KEY[[k]] else NULL
    }

    init_18S <- {
      k <- paste0("18S_", default_year)
      if (k %in% names(RICHNESS_BY_KEY)) RICHNESS_BY_KEY[[k]] else NULL
    }

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
      options     = pathOptions(pane = "pane_polys")
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
        fillOpacity = 0,
        color       = "black",
        weight      = 2,
        opacity     = 1,
        popup       = ~site_name,
        options     = pathOptions(pane = "pane_poly_total")
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
          "Sampling Depth"
        ),
        options = layersControlOptions(collapsed = FALSE)
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
                          'Sampling Depth'
                        ]);

                        const depthName = 'Sampling Depth';

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
                            h1.textContent = 'Species richness by gene region';
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

                        // ---- LEGEND SWAP ----
                          function updateLegends(){
                            const richLegend  = el.querySelector('.legend-richness-box');
                            const depthLegend = el.querySelector('.legend-depth-box');

                            const depthOn = isChecked(depthName);
                            const anyRichOn = anyChecked(richness);

                            // show depth legend only when depth is on
                            if(depthLegend){
                              depthLegend.classList.toggle('legend-hidden', !depthOn);
                            }

                            // show richness legend when any richness layer is on (and depth is off)
                            if(richLegend){
                              richLegend.classList.toggle('legend-hidden', depthOn || !anyRichOn);
                            }
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

  # ---- 2) When the year changes, swap the richness layers ----
  observeEvent(sel_year_chr(), {
    yr <- sel_year_chr()  # "All" or "2023" etc

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

    # --- 12S (only if that key exists) ---
    grid12 <- NULL
    if (yr != "All") {
      k <- paste0("12S_", yr)
      if (k %in% names(RICHNESS_BY_KEY)) grid12 <- RICHNESS_BY_KEY[[k]]
    }
    add_grid_layer(grid12, "12S", "n_species", "12S richness")

    # --- COI ---
    gridCOI <- NULL
    if (yr != "All") {
      k <- paste0("COI_", yr)
      if (k %in% names(RICHNESS_BY_KEY)) gridCOI <- RICHNESS_BY_KEY[[k]]
    }
    add_grid_layer(gridCOI, "COI", "n_species", "COI richness")


    # --- 16S ---
    grid16 <- NULL
    if (yr != "All") {
      k <- paste0("16S_", yr)
      if (k %in% names(RICHNESS_BY_KEY)) grid16 <- RICHNESS_BY_KEY[[k]]
    }
    add_grid_layer(grid16, "16S", "n_species", "16S richness")

    # --- 18S ---
    grid18 <- NULL
    if (yr != "All") {
      k <- paste0("18S_", yr)
      if (k %in% names(RICHNESS_BY_KEY)) grid18 <- RICHNESS_BY_KEY[[k]]
    }
    add_grid_layer(grid18, "18S", "n_species", "18S richness")

    # --- ALL markers ---
    gridALL <- if (yr == "All") {
      RICHNESS_ALL
    } else if (yr %in% names(RICHNESS_ALL_BY_YEAR)) {
      RICHNESS_ALL_BY_YEAR[[yr]]
    } else {
      RICHNESS_ALL
    }
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

    df <- depth_layers[[yr_key]]

    leafletProxy("map") %>%
      clearGroup("Sampling Depth") %>%
      addPolygons(
        data        = df,
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
    grid_all_year <- grid_clip  # use the clipped grid you built once

    cell_poly <- grid_all_year %>% dplyr::filter(cell_id == cid)
    if (nrow(cell_poly) == 0) return(NULL)

    inside <- pts[sf::st_within(pts, sf::st_geometry(cell_poly), sparse = FALSE), , drop = FALSE]
    inside
  })

  #Diversity plots
  # --- build sample x taxon matrix from current selection ---
  comm_mat_mpa <- reactive({
    det <- detections_in_mpa()
    req(det)

    df <- det %>% sf::st_drop_geometry()

    shiny::validate(
      shiny::need(nrow(df) > 0, "No samples fall inside MPA/AOI boundaries for the current Year filter.")
    )

    val_col <- if ("organismQuantity" %in% names(df)) "organismQuantity" else NULL

    rank_col <- input$tax_rank %||% "scientificName"

    shiny::validate(
      shiny::need("organismQuantity" %in% names(df),
                  "organismQuantity column is missing.")
    )

    df2 <- df %>%
      dplyr::mutate(
        samp_name = as.character(samp_name),
        taxon     = as.character(.data[[rank_col]]),
        value     = as.numeric(organismQuantity)
      ) %>%
      dplyr::filter(
        !is.na(taxon), taxon != "",
        !is.na(value)
      ) %>%
      dplyr::group_by(samp_name, taxon) %>%
      dplyr::summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

    mat_wide <- df2 %>%
      tidyr::pivot_wider(
        names_from  = taxon,
        values_from = value,
        values_fill = 0
      )

    mat <- mat_wide %>% dplyr::select(-samp_name) %>% as.data.frame()
    rownames(mat) <- mat_wide$samp_name
    mat
  })

  sample_meta_mpa <- reactive({
    det <- detections_in_mpa()
    req(det)

    det %>%
      sf::st_drop_geometry() %>%
      dplyr::mutate(samp_name = as.character(samp_name)) %>%
      dplyr::distinct(samp_name, site_name, site_type, year)
  })

  # --- sample metadata for grouping/hover (Location etc.) ---
  sample_meta <- reactive({
    det <- selected_detections()
    req(det)

    df <- det %>% sf::st_drop_geometry()

    out <- df %>%
      dplyr::mutate(
        samp_name = as.character(samp_name),
        Location  = if ("Location" %in% names(.)) as.character(Location) else NA_character_,
        year      = if ("year" %in% names(.)) as.character(year) else NA_character_
      ) %>%
      dplyr::distinct(samp_name, Location, year)

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
      samp_name = names(vals),
      alpha_val = as.numeric(vals),
      stringsAsFactors = FALSE
    )
  })

  # --- Alpha boxplot data: MPA/AOI boxes + (optional) drawn polygon box ---

    # 1) Base: richness per sample by MPA/AOI (what you already do)
    alpha_boxplot_df <- reactive({
      yr <- sel_year_chr()

      mat   <- comm_mat_mpa()
      meta  <- sample_meta_mpa()
      alpha <- alpha_metric_vec()

      alpha_mpa <- alpha %>%
        dplyr::left_join(meta, by = "samp_name") %>%
        dplyr::filter(!is.na(site_name), site_name != "") %>%
        dplyr::mutate(group_label = as.character(site_name))


    # 2) add a "Drawn polygon" group (only if a polygon exists)
      poly <- drawn_poly()
      if (!is.null(poly)) {

        pts <- species_sf_all
        if (yr != "All") pts <- pts %>% dplyr::filter(as.character(year) == yr)

        rank_col <- input$tax_rank %||% "scientificName"

        # match comm_mat_mpa(): use organismQuantity if present
        val_col <- if ("organismQuantity" %in% names(pts)) "organismQuantity" else NULL

        inside <- pts[sf::st_within(pts, poly, sparse = FALSE), , drop = FALSE] %>%
          sf::st_drop_geometry() %>%
          dplyr::mutate(
            samp_name = as.character(samp_name),
            taxon     = as.character(.data[[rank_col]]),
            value     = as.numeric(organismQuantity)
          ) %>%
          dplyr::filter(!is.na(samp_name), samp_name != "",
                        !is.na(taxon), taxon != "",
                        !is.na(value)) %>%                     # keep numeric values only
          dplyr::group_by(samp_name, taxon) %>%
          dplyr::summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

        if (nrow(inside) > 0) {

          mat_poly <- inside %>%
            tidyr::pivot_wider(
              names_from  = taxon,
              values_from = value,
              values_fill = 0
            ) %>%
            as.data.frame()

          rownames(mat_poly) <- mat_poly$samp_name
          mat_poly$samp_name <- NULL

          metric <- input$alpha_metric %||% "observed"
          vals_poly <- switch(
            metric,
            observed   = vegan::specnumber(mat_poly),
            shannon    = vegan::diversity(mat_poly, index = "shannon"),
            simpson    = vegan::diversity(mat_poly, index = "simpson"),
            invsimpson = vegan::diversity(mat_poly, index = "invsimpson"),
            ace        = vegan::estimateR(mat_poly)["S.ACE", ],
            vegan::specnumber(mat_poly)
          )

          alpha_drawn <- data.frame(
            samp_name   = names(vals_poly),
            alpha_val   = as.numeric(vals_poly),
            site_name   = "Drawn polygon",
            site_type   = "User",
            group_label = "Drawn polygon",
            stringsAsFactors = FALSE
          )

          alpha_mpa <- dplyr::bind_rows(alpha_mpa, alpha_drawn)
        }
      }

    # 3) Make sure "Drawn polygon" appears at the end (only if present)
    lvls <- unique(alpha_mpa$group_label[alpha_mpa$group_label != "Drawn polygon"])
    if (any(alpha_mpa$group_label == "Drawn polygon")) lvls <- c(lvls, "Drawn polygon")

    alpha_mpa %>%
      dplyr::mutate(group_label = factor(group_label, levels = lvls))
  })

  # --- Alpha diversity boxplot (Plotly) ---
  color_vec <- c("#046c9a", "#5BBCD6", "#ABDDDE", "#446455", "#00A08A","#Fdd262")

  output$alpha_boxplot <- plotly::renderPlotly({
    alpha <- alpha_boxplot_df()

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

    rank_label <- tools::toTitleCase(gsub("_", " ", input$tax_rank))

    metric_label <- paste(metric_label, "at", rank_label, "level")

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
        xaxis = list(title = "Location (Polygon)", showgrid = FALSE, showline = TRUE, linecolor = "black", rangemode = "tozero"),
        yaxis = list(title = metric_label, showgrid = FALSE, showline = TRUE, linecolor = "black",
                     range = c(0, max(alpha$alpha_val, na.rm = TRUE) * 1.05), fixedrange = FALSE),
        margin = list(l = 60, r = 20, t = 50, b = 110),
        showlegend = FALSE   # remove legend if you don’t want repeated labels
      )
  })

  #Beta diversity
  output$beta_pcoa <- plotly::renderPlotly({
    mat <- comm_mat_mpa()
    meta <- sample_meta_mpa()

    shiny::validate(
      shiny::need(nrow(mat) > 2, "Select an area with at least 3 samples to compute a PCoA.")
    )

    # optional: transform to relative abundance to reduce library-size effects
    rs <- rowSums(mat)
    mat_rel <- mat
    mat_rel[rs > 0, ] <- mat[rs > 0, , drop = FALSE] / rs[rs > 0]

    # distance + ordination
    d <- vegan::vegdist(mat_rel, method = "bray")
    ord <- stats::cmdscale(d, k = 2, eig = TRUE)

    scores <- data.frame(
      samp_name = rownames(ord$points),
      PC1 = ord$points[, 1],
      PC2 = ord$points[, 2],
      stringsAsFactors = FALSE
    ) %>%
      dplyr::left_join(meta, by = "samp_name") %>%
      dplyr::mutate(Location = ifelse(is.na(Location) | Location == "", "Unknown", Location))

    # percent variance (cmdscale eigenvalues; can be negative with non-euclidean distances)
    eig <- ord$eig
    eig_pos <- eig[eig > 0]
    pct1 <- if (length(eig_pos) >= 1) round(100 * eig_pos[1] / sum(eig_pos), 1) else NA_real_
    pct2 <- if (length(eig_pos) >= 2) round(100 * eig_pos[2] / sum(eig_pos), 1) else NA_real_

    plotly::plot_ly(
      data = scores,
      x = ~PC1,
      y = ~PC2,
      type = "scatter",
      mode = "markers",
      color = ~Location,
      hovertemplate = paste(
        "<b>Sample:</b> %{customdata[0]}<br>",
        "<b>Location:</b> %{customdata[1]}<br>",
        "<b>PC1:</b> %{x:.3f}<br>",
        "<b>PC2:</b> %{y:.3f}<extra></extra>"
      ),
      customdata = ~cbind(samp_name, Location)
    ) %>%
      plotly::layout(
        title = "Beta diversity (PCoA; Bray–Curtis)",
        xaxis = list(title = paste0("PC1", if (!is.na(pct1)) paste0(" (", pct1, "%)") else "")),
        yaxis = list(title = paste0("PC2", if (!is.na(pct2)) paste0(" (", pct2, "%)") else "")),
        margin = list(l = 60, r = 20, t = 50, b = 60)
      )
  })



  #Krona plot
  output$tax_krona <- taxplore::renderKronaChart({
    det <- selected_detections()

    shiny::validate(
      shiny::need(!is.null(det) && nrow(det) > 0,
                  "Select a cell/polygon (or draw a polygon) to display a Krona chart.")
    )

    det0 <- det %>% sf::st_drop_geometry()

    shiny::validate(
      shiny::need("organismQuantity" %in% names(det0),
                  "organismQuantity column is missing.")
    )

    tax_ranks <- c("kingdom","phylum","class","order","family","genus")

    krona_df <- det0 %>%
      dplyr::mutate(
        species = as.character(scientificName),
        qty     = as.numeric(organismQuantity)
      ) %>%
      dplyr::filter(
        !is.na(species), species != "",
        !is.na(qty), qty > 0
      ) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(c(tax_ranks, "species")))) %>%
      dplyr::summarise(magnitude = sum(qty, na.rm = TRUE), .groups = "drop")

    shiny::validate(
      shiny::need(nrow(krona_df) > 0, "No detections with organismQuantity > 0 in the current selection.")
    )

    tax_df <- krona_df %>%
      dplyr::select(dplyr::all_of(c(tax_ranks, "species")))

    taxplore::plot_krona(
      tax_df,
      magnitude   = krona_df$magnitude,
      total_label = "Sum organismQuantity"
    )
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
    show_total_grid <- "All" %in% groups_on
    show_12S_grid   <- "12S" %in% groups_on
    show_COI_grid   <- "COI" %in% groups_on
    show_16S_grid   <- "16S" %in% groups_on
    show_18S_grid   <- "18S" %in% groups_on
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

      parts  <- strsplit(click$id, "\\|\\|")[[1]]
      p_type <- parts[1]
      p_name <- parts[2]

      # Year-aware lookup
      if (yr == "All") {
        row <- poly_species_all %>%
          dplyr::filter(site_type == p_type, site_name == p_name)
      } else {
        row <- poly_species_year %>%
          dplyr::filter(site_type == p_type, site_name == p_name, year == yr)
      }

      spp <- if (nrow(row) == 0) character(0) else row$species[[1]]
      spp <- apply_interest_filter(spp)

      if (length(spp) == 0) {
        return(em(
          if (yr == "All") "No species records in this MPA/AOI polygon."
          else paste0("No species records in this MPA/AOI polygon for Year = ", yr, ".")
        ))
      }

      return(tagList(
        strong(paste0(p_type, ": ", p_name)),
        tags$br(),
        paste0("Year: ", yr),
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
    # ---- dynamic species extraction for a cell + (optional) gene + year ----
    get_cell_spp <- function(cid, gene = NULL, year = NULL) {
      # gene=NULL + year=NULL = all genes/years
      if (is.null(gene) && is.null(year)) {
        row <- CELL_SPECIES_ALL %>% dplyr::filter(cell_id == cid)
        if (nrow(row) == 0) return(character(0))
        return(row$spp[[1]])
      }

      # gene+year key
      key <- paste(gene, year, sep = "_")
      if (!key %in% names(CELL_SPECIES_BY_KEY)) return(character(0))

      row <- CELL_SPECIES_BY_KEY[[key]] %>% dplyr::filter(cell_id == cid)
      if (nrow(row) == 0) return(character(0))
      row$spp[[1]]
    }

    # total species list respects Year filter
    spp_total <- if (yr == "All") {
      get_cell_spp(cid, gene = NULL, year = NULL)
    } else {
      # all genes but this year: easiest from detections sf (reliable)
      pts_y <- species_sf_all %>% dplyr::filter(as.character(year) == yr)
      inside <- pts_y[sf::st_within(pts_y, sf::st_geometry(grid_clip %>% dplyr::filter(cell_id == cid)), sparse = FALSE), , drop = FALSE]
      sort(unique(as.character(inside$scientificName)))
    }

    spp_total <- apply_interest_filter(spp_total)

    # gene-specific lists only if that gene-year exists
    spp_12S <- if (yr == "All") character(0) else get_cell_spp(cid, gene = "12S", year = yr)
    spp_COI <- if (yr == "All") character(0) else get_cell_spp(cid, gene = "COI", year = yr)
    spp_16S <- if (yr == "All") character(0) else get_cell_spp(cid, gene = "16S", year = yr)
    spp_18S <- if (yr == "All") character(0) else get_cell_spp(cid, gene = "18S", year = yr)

    spp_12S <- apply_interest_filter(spp_12S)
    spp_COI <- apply_interest_filter(spp_COI)
    spp_16S <- apply_interest_filter(spp_16S)
    spp_18S <- apply_interest_filter(spp_18S)

    if (!show_total_grid && !show_12S_grid && !show_COI_grid && !show_16S_grid && !show_18S_grid) {
      return(em("Turn ON a richness layer (all markers / 12S / COI / 16S / 18S) to view the species list for grid cells."))
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
        tags$details(open = (!show_total_grid && !show_12S_grid),
                     tags$summary(strong(paste0("12S species (", length(spp_12S), ")"))),
                     make_list(spp_12S)
        )
      )
    }

    if (show_COI_grid) {
      panels <- tagAppendChildren(
        panels,
        tags$details(open = (!show_total_grid && !show_COI_grid),
                     tags$summary(strong(paste0("COI species (", length(spp_COI), ")"))),
                     make_list(spp_COI)
        )
      )
    }

    if (show_16S_grid) {
      panels <- tagAppendChildren(
        panels,
        tags$details(open = (!show_total_grid && !show_16S_grid),
                     tags$summary(strong(paste0("16S species (", length(spp_16S), ")"))),
                     make_list(spp_16S)
        )
      )
    }

    if (show_18S_grid) {
      panels <- tagAppendChildren(
        panels,
        tags$details(open = (!show_total_grid && !show_18S_grid),
                     tags$summary(strong(paste0("18S species (", length(spp_18S), ")"))),
                     make_list(spp_18S)
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

    det <- det %>%
      sf::st_drop_geometry() %>%
      dplyr::mutate(
        scientificName = as.character(scientificName),
        samp_name      = if ("samp_name" %in% names(.)) as.character(samp_name) else NA_character_,
        #source_file    = if ("source_file" %in% names(.)) as.character(source_file) else NA_character_,
        year           = if ("year" %in% names(.)) as.character(year) else NA_character_,
        target_gene    = if ("target_gene" %in% names(.)) as.character(target_gene) else NA_character_
      )

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
      dplyr::select(dplyr::any_of(c(
        "scientificName","samp_name","eventDate","year","month",
        "target_gene","organismQuantity","organismQuantityType",
        "sampleSizeValue","sampleSizeUnit", "occurrenceID"
      )))%>%
      dplyr::arrange(scientificName, year, samp_name)

    DT::datatable(det_show, rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE))
  })


  # highlight selected cell (outline on top)
  observeEvent(input$map_shape_click, {
    cid <- clicked_cell()
    if (is.na(cid)) return(NULL)  # <- stops crash when click$id is NULL

    sel <- grid_clip %>% dplyr::filter(cell_id == cid)

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

shinyApp(ui, server)


###Updates from meeting

#add multiple polygon options to alpha diversity plot - START HERE
#fix beta diversity PCoA
#fix bugs
#add Luke's code for read_data and protocol cards +NMDS/barchart when he is done

#add NMDS plot for community structure - see code Nick provides - year, season, depth
#add button to download datasets based on selected MPA/AOI polygon and/or drawn polygons

#draft skeleton of paper with a sentence or two of why each was done
#draft email to meet with co-authors in two weeks (after march break)

#Look into journals for publishing function manuscript (maybe ecological informatics)
