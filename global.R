## global.R
source("R/packages.R")
source(file.path("R", "precompute_data.R"))

cache_file <- file.path("cache", "precomputed.rds")

# Auto-rebuild cache when inputs change
inputs <- c(
  list.files(file.path("data", "polygons"), full.names = TRUE, recursive = TRUE),
  list.files(file.path("data", "temporary_occurrence"), full.names = TRUE, recursive = TRUE)
)

cache_is_stale <- function(cache, files) {
  if (!file.exists(cache)) return(TRUE)
  cache_time <- file.info(cache)$mtime
  any(file.info(files)$mtime > cache_time, na.rm = TRUE)
}

if (cache_is_stale(cache_file, inputs)) {
  message("Cache stale (or missing). Recomputing...")
  pre <- precompute_all()
  saveRDS(pre, cache_file)
} else {
  message("Loading cached precomputed objects...")
  pre <- readRDS(cache_file)
}

# Export objects into global env (so server.R / app.R can use them directly)
list2env(pre, envir = .GlobalEnv)
