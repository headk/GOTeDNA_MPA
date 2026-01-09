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

#Clean polygons without losing pieces
all_polys_click <- all_polys %>%
  st_make_valid() %>%
  st_collection_extract("POLYGON") %>%
  filter(!st_is_empty(.)) %>%
  group_by(site_type, site_name) %>%
  summarise(geometry = st_union(st_geometry(.)), .groups = "drop") %>%
  st_as_sf()


#-----------------------------------------------------------------------------
#TEMPORARY CODE UNTIL WE HAVE ACCESS TO THE DDD EXTENSION FROM OBIS!!!
read_occ_csv <- function(...) {
  read_csv(file.path("data", "temporary_occurrence", ...),
  show_col_types = FALSE
  )
}

SABGULFUN_24_COI <- read_occ_csv("OBIS_MCT_SABGULFUN2024_COI_occurrence.csv")
SABGULFUN_24_12S <- read_occ_csv("OBIS_MCT_SABGULFUN2024_12S_occurrence.csv")
MUSQ16_12S <- read_occ_csv("OBIS_MCT_Musq16_occurrence_FORMAP_notcomplete.csv")
MUSQ20_12S <- read_occ_csv("OBIS_MCT_Musq20_occurrence_FORMAP_notcomplete.csv")
MUSQ24_12S <- read_occ_csv("OBIS_MCT_Musq24_occurrence_FORMAP_notcomplete.csv")

SABGULFUN_24_COI <- SABGULFUN_24_COI %>%
  mutate(
    target_gene = "COI",
    scientificName = str_replace_all(scientificName, "_", " ")
  )

SABGULFUN_24_12S <- SABGULFUN_24_12S %>%
  mutate(
    target_gene = "12S"
  )

MUSQ16_12S <- MUSQ16_12S %>%
  mutate(
    target_gene = "12S",
    Location = "Musquash MPA",
    ProtocolID = "16",
    Year = "2024",
    Volume = "1",
    Depth = "0.5",
    Pump = "Smith Root Vacuum Pump",
    Filter = "1.2um Smith Root",
    Preservative = "Self-preserving filters",
    Primers = "MiFishU+248_F(12S)",
    Sequencer = "MiSeq",
    Bioinformatics = "QIIME2 with dada2"
  )

MUSQ20_12S <- MUSQ20_12S %>%
  mutate(
    target_gene = "12S",
    Location = "Musquash MPA",
    ProtocolID = "20",
    Year = "2022",
    Volume = "1",
    Depth = "0.5",
    Pump = "Peristaltic",
    Filter = "0.22um PVDF Sterivex",
    Preservative = "Longmire buffer",
    Primers = "MiFishU (12S)",
    Sequencer = "NovaSeq",
    Bioinformatics = "QIIME2 with dada2"
  )

MUSQ24_12S <- MUSQ24_12S %>%
  mutate(
    target_gene = "12S",
    Location = "Musquash MPA",
    ProtocolID = "24",
    Year = "2023",
    Volume = "1",
    Depth = "0.5",
    Pump = "Smith Root Vacuum Pump",
    Filter = "1.2um Smith Root",
    Preservative = "Self-preserving filters",
    Primers = "MiFishU (12S)",
    Sequencer = "NovaSeq",
    Bioinformatics = "QIIME2 with dada2"
  )

#Make sure columns are in the same integer/character format
fix_ids <- function(df) {
  df %>%
    mutate(
      materialSampleID = as.character(materialSampleID)
    )
}

data_2024_12S <- dplyr::bind_rows(
  fix_ids(SABGULFUN_24_12S) %>% dplyr::mutate(source_file = "SABGULFUN_24_12S"),
  fix_ids(MUSQ16_12S)       %>% dplyr::mutate(source_file = "MUSQ16_12S")
)

data_2024_COI <- dplyr::bind_rows(
  fix_ids(SABGULFUN_24_COI) %>% dplyr::mutate(source_file = "SABGULFUN_24_COI")
)

data_2023_12S <- dplyr::bind_rows(
  fix_ids(MUSQ24_12S)       %>% dplyr::mutate(source_file = "MUSQ24_12S")
)

data_2022_12S <- dplyr::bind_rows(
  fix_ids(MUSQ20_12S)       %>% dplyr::mutate(source_file = "MUSQ20_12S")
)

occ_all <- bind_rows(data_2024_12S, data_2024_COI, data_2023_12S, data_2022_12S)

species_sf_all <- occ_all %>%
  filter(!is.na(decimalLatitude), !is.na(decimalLongitude), occurrenceStatus == "present") %>%
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

# --- Sampling points layer for leaflet (keeps year + source_file if present) ---
sampling_pts <- species_sf_all %>%
  mutate(
    year        = if ("year" %in% names(.)) year else NA,
    source_file = if ("source_file" %in% names(.)) source_file else NA
  ) %>%
  distinct(target_gene, year, source_file, materialSampleID, geometry, .keep_all = TRUE)

sampling_pts <- sampling_pts %>%
  mutate(
    year = if ("year" %in% names(.)) as.character(year) else NA_character_
  )


#############################################################
#Turn species data into sf points  for the spatial join

species_sf_2024_COI <- data_2024_COI %>%
  filter(!is.na(decimalLatitude), !is.na(decimalLongitude),
         occurrenceStatus == "present") %>%
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) %>%
  mutate(marker = "COI")

species_sf_2024_12S <- data_2024_12S %>%
  filter(!is.na(decimalLatitude), !is.na(decimalLongitude),
         occurrenceStatus == "present") %>%
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) %>%
  mutate(marker = "12S")

species_sf_2023_12S <- data_2023_12S %>%
  filter(!is.na(decimalLatitude), !is.na(decimalLongitude),
         occurrenceStatus == "present") %>%
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) %>%
  mutate(marker = "12S")

species_sf_2022_12S <- data_2022_12S %>%
  filter(!is.na(decimalLatitude), !is.na(decimalLongitude),
         occurrenceStatus == "present") %>%
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) %>%
  mutate(marker = "12S")


sampling_pts <- bind_rows(species_sf_2024_COI, species_sf_2024_12S, species_sf_2023_12S, species_sf_2022_12S) %>%
  mutate(
    year        = if ("year" %in% names(.)) year else NA,
    source_file = if ("source_file" %in% names(.)) source_file else NA
  ) %>%
  distinct(marker, year, source_file, materialSampleID, geometry, .keep_all = TRUE)


#Assign each record to a polygon (spatial join)
species_in_polys_2024_COI <- st_join(
  species_sf_2024_COI,
  all_polys,
  join = st_within,
  left = FALSE   # drop species points not in any polygon
)

species_in_polys_2024_12S <- st_join(
  species_sf_2024_12S,
  all_polys,
  join = st_within,
  left = FALSE   # drop species points not in any polygon
)

species_in_polys_2023_12S <- st_join(
  species_sf_2023_12S,
  all_polys,
  join = st_within,
  left = FALSE   # drop species points not in any polygon
)

species_in_polys_2022_12S <- st_join(
  species_sf_2022_12S,
  all_polys,
  join = st_within,
  left = FALSE   # drop species points not in any polygon
)

poly_species_total <- bind_rows(
  species_in_polys_2024_COI %>% st_drop_geometry() %>% transmute(site_name, site_type, scientificName),
  species_in_polys_2024_12S %>% st_drop_geometry() %>% transmute(site_name, site_type, scientificName),
  species_in_polys_2023_12S %>% st_drop_geometry() %>% transmute(site_name, site_type, scientificName),
  species_in_polys_2022_12S %>% st_drop_geometry() %>% transmute(site_name, site_type, scientificName)
) %>%
  filter(!is.na(scientificName), scientificName != "") %>%
  distinct(site_name, site_type, scientificName) %>%
  group_by(site_name, site_type) %>%
  summarise(species = list(sort(unique(scientificName))), .groups = "drop")


species_in_polys_2024_COI %>%
  st_drop_geometry() %>%
  count(site_name, site_type)

species_in_polys_2024_12S %>%
  st_drop_geometry() %>%
  count(site_name, site_type)

species_in_polys_2023_12S %>%
  st_drop_geometry() %>%
  count(site_name, site_type)

species_in_polys_2022_12S %>%
  st_drop_geometry() %>%
  count(site_name, site_type)


#Summary report per polygon

#Total unique species per site
total_species_2024_COI <- species_in_polys_2024_COI %>%
  st_drop_geometry() %>%
  group_by(site_name, site_type) %>%
  summarise(
    n_species_total = n_distinct(scientificName),
    .groups = "drop"
  )

total_species_2024_12S <- species_in_polys_2024_12S %>%
  st_drop_geometry() %>%
  group_by(site_name, site_type) %>%
  summarise(
    n_species_total = n_distinct(scientificName),
    .groups = "drop"
  )

total_species_2023_12S <- species_in_polys_2023_12S %>%
  st_drop_geometry() %>%
  group_by(site_name, site_type) %>%
  summarise(
    n_species_total = n_distinct(scientificName),
    .groups = "drop"
  )

total_species_2022_12S <- species_in_polys_2022_12S %>%
  st_drop_geometry() %>%
  group_by(site_name, site_type) %>%
  summarise(
    n_species_total = n_distinct(scientificName),
    .groups = "drop"
  )


#Unique species by taxonomic class per site

species_by_class_2024_COI <- species_in_polys_2024_COI %>%
  st_drop_geometry() %>%
  group_by(site_name, site_type, class) %>%
  summarise(
    n_species = n_distinct(scientificName),
    .groups = "drop"
  ) %>%
  arrange(site_name, class)

species_by_class_2024_12S <- species_in_polys_2024_12S %>%
  st_drop_geometry() %>%
  group_by(site_name, site_type, class) %>%
  summarise(
    n_species = n_distinct(scientificName),
    .groups = "drop"
  ) %>%
  arrange(site_name, class)

species_by_class_2023_12S <- species_in_polys_2023_12S %>%
  st_drop_geometry() %>%
  group_by(site_name, site_type, class) %>%
  summarise(
    n_species = n_distinct(scientificName),
    .groups = "drop"
  ) %>%
  arrange(site_name, class)

species_by_class_2022_12S <- species_in_polys_2022_12S %>%
  st_drop_geometry() %>%
  group_by(site_name, site_type, class) %>%
  summarise(
    n_species = n_distinct(scientificName),
    .groups = "drop"
  ) %>%
  arrange(site_name, class)

species_by_class_wide_2024_COI <- species_by_class_2024_COI %>%
  pivot_wider(
    names_from  = class,
    values_from = n_species,
    values_fill = 0
  )

species_by_class_wide_2024_12S <- species_by_class_2024_12S %>%
  pivot_wider(
    names_from  = class,
    values_from = n_species,
    values_fill = 0
  )

species_by_class_wide_2023_12S <- species_by_class_2023_12S %>%
  pivot_wider(
    names_from  = class,
    values_from = n_species,
    values_fill = 0
  )

species_by_class_wide_2022_12S <- species_by_class_2022_12S %>%
  pivot_wider(
    names_from  = class,
    values_from = n_species,
    values_fill = 0
  )

#IUCN and AIS (when I add them from OBIS)
#species_by_iucn <- species_in_polys %>%
#  st_drop_geometry() %>%
#  group_by(site_name, site_type, iucn_status) %>%
#  summarise(
#    n_species = n_distinct(scientificName),
#    .groups = "drop"
#  )

#alien_summary <- species_in_polys %>%
#  st_drop_geometry() %>%
#  filter(alien_invasive == "yes") %>%
#  group_by(site_name, site_type) %>%
#  summarise(
#    n_alien_species = n_distinct(scientificName),
#    .groups = "drop"
#  )

#Combined summary table
summary_report_2024_COI <- total_species_2024_COI %>%
  #left_join(alien_summary,
  #         by = c("site_name", "site_type")) %>%
  #mutate(
  #  n_alien_species = ifelse(is.na(n_alien_species), 0L, n_alien_species)
  #) %>%
  left_join(species_by_class_wide_2024_COI,
            by = c("site_name", "site_type"))

summary_report_2024_12S <- total_species_2024_12S %>%
  #left_join(alien_summary,
  #         by = c("site_name", "site_type")) %>%
  #mutate(
  #  n_alien_species = ifelse(is.na(n_alien_species), 0L, n_alien_species)
  #) %>%
  left_join(species_by_class_wide_2024_12S,
            by = c("site_name", "site_type"))

summary_report_2023_12S <- total_species_2023_12S %>%
  #left_join(alien_summary,
  #         by = c("site_name", "site_type")) %>%
  #mutate(
  #  n_alien_species = ifelse(is.na(n_alien_species), 0L, n_alien_species)
  #) %>%
  left_join(species_by_class_wide_2023_12S,
            by = c("site_name", "site_type"))

summary_report_2022_12S <- total_species_2022_12S %>%
  #left_join(alien_summary,
  #         by = c("site_name", "site_type")) %>%
  #mutate(
  #  n_alien_species = ifelse(is.na(n_alien_species), 0L, n_alien_species)
  #) %>%
  left_join(species_by_class_wide_2022_12S,
            by = c("site_name", "site_type"))

#List unique species names per polygon

species_list_2024_COI <- species_in_polys_2024_COI %>%
  st_drop_geometry() %>%
  group_by(site_name, site_type) %>%
  summarise(
    species = list(unique(scientificName)),
    .groups = "drop"
  )

species_list_2024_12S <- species_in_polys_2024_12S %>%
  st_drop_geometry() %>%
  group_by(site_name, site_type) %>%
  summarise(
    species = list(unique(scientificName)),
    .groups = "drop"
  )

species_list_2023_12S <- species_in_polys_2023_12S %>%
  st_drop_geometry() %>%
  group_by(site_name, site_type) %>%
  summarise(
    species = list(unique(scientificName)),
    .groups = "drop"
  )

species_list_2022_12S <- species_in_polys_2022_12S %>%
  st_drop_geometry() %>%
  group_by(site_name, site_type) %>%
  summarise(
    species = list(unique(scientificName)),
    .groups = "drop"
  )

species_long_2024_COI <- species_in_polys_2024_COI %>%
  st_drop_geometry() %>%
  distinct(site_name, site_type, scientificName) %>%
  arrange(site_name, scientificName)

species_long_2024_12S <- species_in_polys_2024_12S %>%
  st_drop_geometry() %>%
  distinct(site_name, site_type, scientificName) %>%
  arrange(site_name, scientificName)

species_long_2023_12S <- species_in_polys_2023_12S %>%
  st_drop_geometry() %>%
  distinct(site_name, site_type, scientificName) %>%
  arrange(site_name, scientificName)

species_long_2022_12S <- species_in_polys_2022_12S %>%
  st_drop_geometry() %>%
  distinct(site_name, site_type, scientificName) %>%
  arrange(site_name, scientificName)


summary_report_with_species_2024_COI <- summary_report_2024_COI %>%
  left_join(
    species_list_2024_COI,
    by = c("site_name", "site_type")
  )

summary_report_with_species_2024_12S <- summary_report_2024_12S %>%
  left_join(
    species_list_2024_12S,
    by = c("site_name", "site_type")
  )

summary_report_with_species_2023_12S <- summary_report_2023_12S %>%
  left_join(
    species_list_2023_12S,
    by = c("site_name", "site_type")
  )

summary_report_with_species_2022_12S <- summary_report_2022_12S %>%
  left_join(
    species_list_2022_12S,
    by = c("site_name", "site_type")
  )



##Species Richness Polygons

#Turn occurrences into sf points

# Starting from your original occurrence dataframe

species_sf_2024_COI <- data_2024_COI %>%
  filter(
    !is.na(decimalLatitude),
    !is.na(decimalLongitude),
    occurrenceStatus == "present"       # only present detections
  ) %>%
  st_as_sf(
    coords = c("decimalLongitude", "decimalLatitude"),
    crs = 4326
  )

species_sf_2024_12S <- data_2024_12S %>%
  filter(
    !is.na(decimalLatitude),
    !is.na(decimalLongitude),
    occurrenceStatus == "present"       # only present detections
  ) %>%
  st_as_sf(
    coords = c("decimalLongitude", "decimalLatitude"),
    crs = 4326
  )

species_sf_2023_12S <- data_2023_12S %>%
  filter(
    !is.na(decimalLatitude),
    !is.na(decimalLongitude),
    occurrenceStatus == "present"       # only present detections
  ) %>%
  st_as_sf(
    coords = c("decimalLongitude", "decimalLatitude"),
    crs = 4326
  )

species_sf_2022_12S <- data_2022_12S %>%
  filter(
    !is.na(decimalLatitude),
    !is.na(decimalLongitude),
    occurrenceStatus == "present"       # only present detections
  ) %>%
  st_as_sf(
    coords = c("decimalLongitude", "decimalLatitude"),
    crs = 4326
  )


# Make a grid over polygons
# all_polys = your combined MPA + AOI polygons with site_name, site_type

grid <- st_make_grid(
  all_polys,
  cellsize = 0.018,           # tweak this: smaller = finer grid, more cells
  square   = TRUE
) %>%
  st_as_sf() %>%
  mutate(cell_id = row_number())

# ---- NEW: cell has >=1 sampling point? (by dataset/year) ----
# Use the same point sets you used to compute richness.

haspt_2024_COI <- lengths(st_intersects(grid, species_sf_2024_COI)) > 0
haspt_2024_12S <- lengths(st_intersects(grid, species_sf_2024_12S)) > 0
haspt_2023_12S <- lengths(st_intersects(grid, species_sf_2023_12S)) > 0
haspt_2022_12S <- lengths(st_intersects(grid, species_sf_2022_12S)) > 0

# For the "All markers" (total) grid:
haspt_ALL <- lengths(st_intersects(grid, species_sf_all)) > 0

#Compute species richness per grid cell
# join points to grid cells
grid_with_pts_2024_COI <- st_join(
  grid,
  species_sf_2024_COI,
  join = st_intersects,
  left = TRUE
)

grid_with_pts_2024_12S <- st_join(
  grid,
  species_sf_2024_12S,
  join = st_intersects,
  left = TRUE
)

grid_with_pts_2023_12S <- st_join(
  grid,
  species_sf_2023_12S,
  join = st_intersects,
  left = TRUE
)

grid_with_pts_2022_12S <- st_join(
  grid,
  species_sf_2022_12S,
  join = st_intersects,
  left = TRUE
)

# richness per cell = number of unique species in that cell
richness_grid_2024_COI <- grid_with_pts_2024_COI %>%
  st_drop_geometry() %>%
  group_by(cell_id) %>%
  summarise(
    n_species = n_distinct(scientificName[!is.na(scientificName)]),
    .groups   = "drop"
  ) %>%
  right_join(grid, by = "cell_id") %>%
  st_as_sf() %>%
  mutate(
    n_species     = ifelse(is.na(n_species), 0L, n_species),
    has_sampling  = haspt_2024_COI[cell_id]
  )

richness_grid_2024_12S <- grid_with_pts_2024_12S %>%
  st_drop_geometry() %>%
  group_by(cell_id) %>%
  summarise(
    n_species = n_distinct(scientificName[!is.na(scientificName)]),
    .groups   = "drop"
  ) %>%
  right_join(grid, by = "cell_id") %>%
  st_as_sf() %>%
  mutate(
    n_species     = ifelse(is.na(n_species), 0L, n_species),
    has_sampling  = haspt_2024_12S[cell_id]
  )

richness_grid_2023_12S <- grid_with_pts_2023_12S %>%
  st_drop_geometry() %>%
  group_by(cell_id) %>%
  summarise(
    n_species = n_distinct(scientificName[!is.na(scientificName)]),
    .groups   = "drop"
  ) %>%
  right_join(grid, by = "cell_id") %>%
  st_as_sf() %>%
  mutate(
    n_species     = ifelse(is.na(n_species), 0L, n_species),
    has_sampling  = haspt_2023_12S[cell_id]
  )

richness_grid_2022_12S <- grid_with_pts_2022_12S %>%
  st_drop_geometry() %>%
  group_by(cell_id) %>%
  summarise(
    n_species = n_distinct(scientificName[!is.na(scientificName)]),
    .groups   = "drop"
  ) %>%
  right_join(grid, by = "cell_id") %>%
  st_as_sf() %>%
  mutate(
    n_species     = ifelse(is.na(n_species), 0L, n_species),
    has_sampling  = haspt_2022_12S[cell_id]
  )

richness_grid_2024_COI <- richness_grid_2024_COI %>%
  mutate(n_species = ifelse(is.na(n_species), 0L, n_species))

richness_grid_2024_12S <- richness_grid_2024_12S %>%
  mutate(n_species = ifelse(is.na(n_species), 0L, n_species))

richness_grid_2023_12S <- richness_grid_2023_12S %>%
  mutate(n_species = ifelse(is.na(n_species), 0L, n_species))

richness_grid_2022_12S <- richness_grid_2022_12S %>%
  mutate(n_species = ifelse(is.na(n_species), 0L, n_species))


richness_grid_in_poly_2024_COI <- st_intersection(
  richness_grid_2024_COI,
  all_polys
)

richness_grid_in_poly_2024_12S <- st_intersection(
  richness_grid_2024_12S,
  all_polys
)

richness_grid_in_poly_2023_12S <- st_intersection(
  richness_grid_2023_12S,
  all_polys
)

richness_grid_in_poly_2022_12S <- st_intersection(
  richness_grid_2022_12S,
  all_polys
)

#Rough copy of interactive species richness plot below

richness_grid_in_poly_2024_COI_poly <- richness_grid_in_poly_2024_COI %>%
  st_collection_extract("POLYGON")

richness_grid_in_poly_2024_COI_poly <- richness_grid_in_poly_2024_COI_poly[
  !st_is_empty(richness_grid_in_poly_2024_COI_poly),
]


richness_grid_in_poly_2024_12S_poly <- richness_grid_in_poly_2024_12S %>%
  st_collection_extract("POLYGON")

richness_grid_in_poly_2024_12S_poly <- richness_grid_in_poly_2024_12S_poly[
  !st_is_empty(richness_grid_in_poly_2024_12S_poly),
]

richness_grid_in_poly_2023_12S_poly <- richness_grid_in_poly_2023_12S %>%
  st_collection_extract("POLYGON")

richness_grid_in_poly_2023_12S_poly <- richness_grid_in_poly_2023_12S_poly[
  !st_is_empty(richness_grid_in_poly_2023_12S_poly),
]

richness_grid_in_poly_2022_12S_poly <- richness_grid_in_poly_2022_12S %>%
  st_collection_extract("POLYGON")

richness_grid_in_poly_2022_12S_poly <- richness_grid_in_poly_2022_12S_poly[
  !st_is_empty(richness_grid_in_poly_2022_12S_poly),
]

# (optional) outlines, if needed
all_polys_poly <- all_polys %>%
  st_collection_extract("POLYGON")

all_polys_poly <- all_polys_poly[!st_is_empty(all_polys_poly), ]

st_geometry_type(richness_grid_in_poly_2024_COI_poly)
st_geometry_type(richness_grid_in_poly_2024_12S_poly)
st_geometry_type(richness_grid_in_poly_2023_12S_poly)
st_geometry_type(richness_grid_in_poly_2022_12S_poly)

#Combine 12s + coi points into one sf

species_sf_all <- bind_rows(
  species_sf_2024_COI %>% mutate(marker = "COI", year = "2024"),
  species_sf_2024_12S %>% mutate(marker = "12S", year = "2024"),
  species_sf_2023_12S %>% mutate(marker = "12S", year = "2023"),
  species_sf_2022_12S %>% mutate(marker = "12S", year = "2022")
)

#join combined points to the grid

grid_with_pts_all <- st_join(
  grid,
  species_sf_all,
  join = st_intersects,
  left = TRUE
)


#Total richness per cell (no double counting)

richness_grid_all <- grid_with_pts_all %>%
  st_drop_geometry() %>%
  group_by(cell_id) %>%
  summarise(
    n_species_total = n_distinct(scientificName[!is.na(scientificName)]),
    .groups = "drop"
  ) %>%
  right_join(grid, by = "cell_id") %>%
  st_as_sf() %>%
  mutate(
    n_species_total = ifelse(is.na(n_species_total), 0L, n_species_total),
    has_sampling    = haspt_ALL[cell_id]
  )

richness_grid_all_in_poly <- st_intersection(
  richness_grid_all,
  all_polys
)

#Add total richness as a leaflet layer
richness_grid_all_poly <- richness_grid_all_in_poly %>%
  st_collection_extract("POLYGON")
richness_grid_all_poly <- richness_grid_all_poly[!st_is_empty(richness_grid_all_poly), ]


## ===== Unified richness palettes + leaflet map (Wes Anderson, shared scale) =====

## 1. Make sure all richness layers & outlines are polygons (no geometry collections)

richness_grid_all_poly <- richness_grid_all_in_poly %>%
  st_collection_extract("POLYGON")
richness_grid_all_poly <- richness_grid_all_poly[!st_is_empty(richness_grid_all_poly), ]

richness_grid_in_poly_2024_COI_poly <- richness_grid_in_poly_2024_COI %>%
  st_collection_extract("POLYGON")
richness_grid_in_poly_2024_COI_poly <- richness_grid_in_poly_2024_COI_poly[!st_is_empty(richness_grid_in_poly_2024_COI_poly), ]

richness_grid_in_poly_2024_12S_poly <- richness_grid_in_poly_2024_12S %>%
  st_collection_extract("POLYGON")
richness_grid_in_poly_2024_12S_poly <- richness_grid_in_poly_2024_12S_poly[!st_is_empty(richness_grid_in_poly_2024_12S_poly), ]

richness_grid_in_poly_2023_12S_poly <- richness_grid_in_poly_2023_12S %>%
  st_collection_extract("POLYGON")
richness_grid_in_poly_2023_12S_poly <- richness_grid_in_poly_2023_12S_poly[!st_is_empty(richness_grid_in_poly_2023_12S_poly), ]

richness_grid_in_poly_2022_12S_poly <- richness_grid_in_poly_2022_12S %>%
  st_collection_extract("POLYGON")
richness_grid_in_poly_2022_12S_poly <- richness_grid_in_poly_2022_12S_poly[!st_is_empty(richness_grid_in_poly_2022_12S_poly), ]


all_polys_poly <- all_polys %>%
  st_collection_extract("POLYGON")
all_polys_poly <- all_polys_poly[!st_is_empty(all_polys_poly), ]


## 2. Shared richness domain across TOTAL, 12S, and COI

max_rich <- max(
  richness_grid_all_poly$n_species_total,
  richness_grid_in_poly_2024_COI_poly$n_species,
  richness_grid_in_poly_2024_12S_poly$n_species,
  richness_grid_in_poly_2023_12S_poly$n_species,
  richness_grid_in_poly_2022_12S_poly$n_species,
  na.rm = TRUE
)

rich_domain <- c(0, max_rich)


## 3. One Wes Anderson palette for everything (Zissou1 continuous)

wes_cont <- function(name, n = 100) {
  colorRampPalette(wes_palette(name, type = "continuous"))(n)
}

pal_vec <- wes_cont("Zissou1", 100)

pal_rich <- colorNumeric(
  palette  = pal_vec,
  domain   = rich_domain,      # SAME scale for all 3 layers
  na.color = "transparent"
)


###########################################

# A) zone outlines (keep every zone piece)
all_polys_zones <- all_polys %>%
  st_make_valid() %>%
  st_collection_extract("POLYGON") %>%
  filter(!st_is_empty(.)) %>%
  st_as_sf()

# B) clickable dissolved outlines (one feature per site_name/site_type)
all_polys_click <- all_polys_zones %>%
  group_by(site_type, site_name) %>%
  summarise(do_union = TRUE, .groups = "drop")


## ---- 1) Build cell -> species lookup tables (fast at click-time) ----
cell_species_all <- grid_with_pts_all %>%
  st_drop_geometry() %>%
  transmute(
    cell_id         = cell_id,
    target_gene     = target_gene,   # from species_sf_all
    year            = year,          # from species_sf_all
    scientificName  = scientificName
  ) %>%
  filter(!is.na(scientificName), scientificName != "") %>%
  distinct(cell_id, target_gene, year, scientificName)

cell_species_2024_COI <- grid_with_pts_2024_COI %>%
  st_drop_geometry() %>%
  transmute(
    cell_id = cell_id,
    scientificName = scientificName
  ) %>%
  filter(!is.na(scientificName), scientificName != "") %>%
  distinct(cell_id, scientificName)

cell_species_2024_12S <- grid_with_pts_2024_12S %>%
  st_drop_geometry() %>%
  transmute(
    cell_id = cell_id,
    scientificName = scientificName
  ) %>%
  filter(!is.na(scientificName), scientificName != "") %>%
  distinct(cell_id, scientificName)

cell_species_2023_12S <- grid_with_pts_2023_12S %>%
  st_drop_geometry() %>%
  transmute(
    cell_id = cell_id,
    scientificName = scientificName
  ) %>%
  filter(!is.na(scientificName), scientificName != "") %>%
  distinct(cell_id, scientificName)

cell_species_2022_12S <- grid_with_pts_2022_12S %>%
  st_drop_geometry() %>%
  transmute(
    cell_id = cell_id,
    scientificName = scientificName
  ) %>%
  filter(!is.na(scientificName), scientificName != "") %>%
  distinct(cell_id, scientificName)


# convenience: precompute total unique spp per cell (across markers)
cell_species_total <- cell_species_all %>%
  group_by(cell_id) %>%
  summarise(
    n_species_total = n_distinct(scientificName),
    .groups = "drop"
  )

#Define polygon layers

# Richness layers keyed by gene+year
rich_layers <- list(
  "All_All"   = richness_grid_all_poly,
  "COI_2024"  = richness_grid_in_poly_2024_COI_poly,
  "12S_2024"  = richness_grid_in_poly_2024_12S_poly,
  "12S_2023"  = richness_grid_in_poly_2023_12S_poly,
  "12S_2022"  = richness_grid_in_poly_2022_12S_poly
)

# helper: which layer should be visible for current filters?
selected_rich_key <- function(gene, year) {
  if (identical(gene, "All") && identical(year, "All")) return("TOTAL_ALL")
  if (identical(gene, "All") && !identical(year, "All")) return(NULL)  # no "all genes within a single year" layer
  if (!identical(gene, "All") && identical(year, "All")) return(NULL)  # no "all years within a gene" layer
  paste(gene, year, sep = "_")
}


#Organize per year layers into names lists

# --- Richness grids by year ---
grid_12S_by_year <- list(
  "2022" = richness_grid_in_poly_2022_12S_poly,
  "2023" = richness_grid_in_poly_2023_12S_poly,
  "2024" = richness_grid_in_poly_2024_12S_poly
)

# If you ONLY have COI for 2024 right now:
grid_COI_by_year <- list(
  "2024" = richness_grid_in_poly_2024_COI_poly
)

# If you have year-specific "all markers" grids, do the same:
# grid_ALL_by_year <- list("2022"=..., "2023"=..., "2024"=...)
# If you do NOT, keep one static object:

grid_ALL_static <- richness_grid_all_poly %>%
  sf::st_as_sf()

if (!"geometry" %in% names(richness_grid_all_poly) && "x" %in% names(richness_grid_all_poly)) {
  grid_ALL_static <- richness_grid_all_poly %>% dplyr::rename(geometry = x) %>% sf::st_as_sf()
}


# ---- generic helper to compute richness + has_sampling for a given grid and points ----
rich_on_grid <- function(grid_sf, pts_sf, value_col = "scientificName") {
  haspt <- lengths(st_intersects(grid_sf, pts_sf)) > 0
  
  grid_with <- st_join(grid_sf, pts_sf, join = st_intersects, left = TRUE)
  
  out <- grid_with |>
    st_drop_geometry() |>
    group_by(cell_id) |>
    summarise(
      n_species = n_distinct(.data[[value_col]][!is.na(.data[[value_col]])]),
      .groups = "drop"
    ) |>
    right_join(grid_sf, by = "cell_id") |>
    st_as_sf() |>
    mutate(
      n_species    = ifelse(is.na(n_species), 0L, n_species),
      has_sampling = haspt[cell_id]
    )
  
  st_intersection(out, all_polys) |>
    st_collection_extract("POLYGON") |>
    (\(x) x[!st_is_empty(x), ])()
}

##NOTE: Right now the SARA Schedule 1 filter is matching based on the data inputted into the app. Once linking the code to OBIS data, edit so that it matches based on WoRMS AphiaID
    
#Call for the SARA Schedule 1 list and AIS list
read_sara_ais_xlsx <- function(...) {
  readxl::read_xlsx(file.path("data", "sara_ais", ...))
}

SARA <- read_sara_ais_xlsx("SARA_Clean_Schedule1_specieslist.xlsx")
AIS  <- read_sara_ais_xlsx("Target_AIS_List_Claudio.xlsx")

##Fix grid heatmap by year for all markers

# If x is already an sfc column (most likely):
grid_ALL_static <- richness_grid_all_poly %>%
  rename(geometry = x) %>%
  st_as_sf()

# If CRS got lost, set it (use your real CRS; your app is 4326)
st_crs(grid_ALL_static) <- 4326

make_all_marker_grid_from_detections <- function(grid_template_sf, det_sf, year_val) {
  stopifnot(inherits(grid_template_sf, "sf"))
  stopifnot(inherits(det_sf, "sf"))
  stopifnot("cell_id" %in% names(grid_template_sf))
  
  # Keep geometry + IDs only (prevents .x/.y join suffix issues)
  grid_base <- grid_template_sf %>%
    dplyr::select(cell_id, geometry)
  
  det_year <- det_sf
  if (!identical(as.character(year_val), "All")) {
    det_year <- det_year %>% dplyr::filter(as.character(year) == as.character(year_val))
  }
  
  if (nrow(det_year) == 0) {
    return(grid_base %>% dplyr::mutate(n_species_total = 0L, has_sampling = FALSE))
  }
  
  # (Optional but recommended) project to a planar CRS before st_within
  # Pick something appropriate for your area; 3347 works well for Canada
  grid_planar <- sf::st_transform(grid_base, 3347)
  det_planar  <- sf::st_transform(det_year, 3347)
  
  det_with_cell <- sf::st_join(
    det_planar,
    grid_planar %>% dplyr::select(cell_id),
    join = sf::st_within,
    left = FALSE
  )
  
  counts <- det_with_cell %>%
    sf::st_drop_geometry() %>%
    dplyr::group_by(cell_id) %>%
    dplyr::summarise(
      n_species_total = dplyr::n_distinct(scientificName),
      has_sampling    = TRUE,
      .groups = "drop"
    )
  
  out <- grid_planar %>%
    dplyr::left_join(counts, by = "cell_id") %>%
    dplyr::mutate(
      has_sampling    = dplyr::coalesce(has_sampling, FALSE),
      n_species_total = dplyr::coalesce(n_species_total, 0L)
    )
  
  # return back in WGS84 for leaflet
  sf::st_transform(out, 4326)
}


grid_ALL_by_year <- list(
  "2022" = make_all_marker_grid_from_detections(grid_ALL_static, species_sf_all, "2022"),
  "2023" = make_all_marker_grid_from_detections(grid_ALL_static, species_sf_all, "2023"),
  "2024" = make_all_marker_grid_from_detections(grid_ALL_static, species_sf_all, "2024")
)

