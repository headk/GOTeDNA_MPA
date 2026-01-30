# global.R
source("R/packages.R")
source(file.path("R", "precompute_data.R"))

cache_file <- file.path("cache", "precomputed.rds")

# Export objects into global env (so server.R / app.R can use them directly)
list2env(pre, envir = .GlobalEnv)
