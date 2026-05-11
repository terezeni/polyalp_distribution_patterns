# NC cost-distance
## from full_geo_lim_3035 (lithology shapefile)
## calculate conductance raster for 3 different categories of plants (substrate preference)
## using different conversion tables and plot
## TZ 23.1.25


## SETUP ----
library(tidyverse)
library(rio)

library(sf)
library(terra)
library(raster)
library(sp)
library(gdistance)

library(tidyterra)
library(patchwork)

# DIRECTORIES ----
project_dir <- here::here()
data_dir <- file.path(project_dir, "data_in")
gis_dir  <- file.path(project_dir, "GIS_data")


## DATA ----
# import GeoLiM Donnini et al 2019
original_geo_v <- terra::vect(paste0(gis_dir,"/full_geo_lim_3035_V2.gpkg"),
                              layer = "full_final_map_diss")

# reference raster (sensible extent and resolution 100m)
ref_raster <- terra::rast(xmin = 3900000, 
                          xmax = 4900000, 
                          ymin = 2300000, 
                          ymax = 2900000)
res(ref_raster) <- 100  # 100 but then aggregate 4 pixels

# import conversion table (3 sheet in conversion_tables.xlsx), conversion info in Supplementary Table 1
conversion_tables <- list("v0" = rio::import("data_in/conversion_tables.xlsx", which = 1), 
                          "vA" = rio::import("data_in/conversion_tables.xlsx", which = 2), 
                          "vB" = rio::import("data_in/conversion_tables.xlsx", which = 3))

## FUNCTIONS ----
# helper function for rasterizing and aggregating
process_raster <- function(geo, ref_raster, value_column, output_prefix) {
  raster <- rasterize(
    geo, 
    ref_raster,
    value_column,
    filename = paste0(output_prefix, ".tif"),
    overwrite = TRUE
  )
  # 200 m resolution
  aggregated_raster <- terra::aggregate(
    raster, 
    fact = 2,
    filename = paste0(output_prefix, "_200.tif"),
    fun = function(x) mean(x, na.rm = TRUE),
    overwrite = TRUE
  )
  return(list(raster = raster, aggregated_raster = aggregated_raster))
}

## PROCESSING ----
# Process each conversion table
for (name in names(conversion_tables)) {
  geo <- tidyterra::left_join(original_geo_v, 
                              conversion_tables[[name]])
  output_dir <- paste0("conductance_rast/", name,"/")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  
  
  # Rasterize and aggregate for each mapped value column
  process_raster(geo, ref_raster, "mapped_value_S", paste0(output_dir, "geo_S"))
  process_raster(geo, ref_raster, "mapped_value_K", paste0(output_dir, "geo_K"))
  process_raster(geo, ref_raster, "mapped_value_I", paste0(output_dir, "geo_I"))
}

## DATA ----

# import sites
points_sp <- sf::st_read(paste0(gis_dir,"/mountains_3035_46.gpkg"))

plot(points_sp)

# import refugia shapefiles
polygon_sp_K <- sf::st_read(paste0(gis_dir,"/K_refs_complete.shp"))  # on K 
polygon_sp_K <- polygon_sp_K[polygon_sp_K$id %in% c(2,3,4,5,8,9,10),] #!subset excluded not supported ref
polygon_sp_S <- sf::st_read(paste0(gis_dir,"/S_refs_complete.shp"))  # on S

polygon_sp_I <- rbind(polygon_sp_K, polygon_sp_S)  #maybe exclude some
polygon_sp_I$id <- paste0(polygon_sp_I$geology, polygon_sp_I$id) 

# Define versions and their corresponding folders
versions <- c("v0", "vA", "vB")
conductance_base_dirs <- paste0("conductance_rast/", versions)

# Helper function to set correct CRS and load raster
load_and_set_crs <- function(filepath, crs_code = "EPSG:3035") {
  raster <- terra::rast(filepath)
  crs(raster) <- crs_code
  return(raster)
}


# Outer loop for each version
for (version_idx in seq_along(versions)) {
  version <- versions[version_idx]
  base_dir <- conductance_base_dirs[version_idx]
  
  cat(sprintf("Processing %s...\n", version))
  
  # Load conductance rasters for the current version
  conductance_l <- list(
    I = load_and_set_crs(file.path(base_dir, "geo_I_200.tif")),
    K = load_and_set_crs(file.path(base_dir, "geo_K_200.tif")),
    S = load_and_set_crs(file.path(base_dir, "geo_S_200.tif"))
  )
  
  # Polygon list
  polygon_sp_l <- list(I = polygon_sp_I, K = polygon_sp_K, S = polygon_sp_S)
  
  # Prepare coordinates of points
  points_sp_points <- sf::st_cast(points_sp, "POINT")  # Convert MULTIPOINT to POINT
  from_coords_matrix <- sf::st_coordinates(points_sp_points)[, 1:2]
  
  # Function to calculate cost distance
  calculate_cost_distance <- function(conductance, polygon_sp, from_coords) {
    polygon_boundaries <- sf::st_boundary(polygon_sp)
    to_coords_matrix <- sf::st_coordinates(polygon_boundaries)[, 1:2]
    
    tr_matrix <- gdistance::transition(raster(conductance), 
                                       transitionFunction = mean, 
                                       directions = 8)
    cost_dist <- gdistance::costDistance(tr_matrix, 
                                         fromCoords = from_coords, 
                                         toCoords = to_coords_matrix)
    return(cost_dist)
  }
  
  # Calculate cost distances for all conductance types
  distance_l <- lapply(names(conductance_l), function(name) {
    cat(sprintf("Processing %s conductance...\n", name))
    calculate_cost_distance(conductance_l[[name]], polygon_sp_l[[name]], from_coords_matrix)
  })
  names(distance_l) <- names(conductance_l)
  
  # Save results for the current version
  save(distance_l, file = file.path(base_dir, "distance_l.RData"))
  distance_l_backup <- distance_l
  
  # Compute minimum distances
  minimum_distances <- lapply(distance_l, function(cost_dist) {
    apply(cost_dist, 1, FUN = min, na.rm = TRUE)
  })
  
  # Assign minimum distances to points_sp
  points_sp$cost_dist_I <- minimum_distances$I
  points_sp$cost_dist_K <- minimum_distances$K
  points_sp$cost_dist_S <- minimum_distances$S
  
  # Correct distances: if the point is in a refugium the distance is 0
  correct_distances <- function(points_sp, column_name, ref_polygon) {
    within_ref <- sf::st_within(points_sp, ref_polygon, sparse = FALSE) |> apply(1, any)
    corrected <- ifelse(within_ref, 0, points_sp[[column_name]])
    return(corrected)
  }
  
  points_sp$cost_dist_I_corr <- correct_distances(points_sp, "cost_dist_I", polygon_sp_I)
  points_sp$cost_dist_K_corr <- correct_distances(points_sp, "cost_dist_K", polygon_sp_K)
  points_sp$cost_dist_S_corr <- correct_distances(points_sp, "cost_dist_S", polygon_sp_S)
  
  # Calculate straight-line distances and correct
  calculate_straight_distance <- function(points_sp, polygon_sp, ref_name) {
    dist <- sf::st_distance(points_sp, polygon_sp) |> apply(1, min)
    within_ref <- sf::st_within(points_sp, polygon_sp, sparse = FALSE) |> apply(1, any)
    corrected <- ifelse(within_ref, 0, dist)
    return(corrected)
  }
  
  points_sp$dist_I_corr <- calculate_straight_distance(points_sp, polygon_sp_I, "I")
  points_sp$dist_K_corr <- calculate_straight_distance(points_sp, polygon_sp_K, "K")
  points_sp$dist_S_corr <- calculate_straight_distance(points_sp, polygon_sp_S, "S")
  
  # Save points_sp results for the current version
  sf::st_write(points_sp, file.path(base_dir, "points_sp_processed.gpkg"), delete_dsn = TRUE)
  write.csv(points_sp, file.path(base_dir, "points_sp_processed.csv"), row.names = T, sep = ";")
  cat(sprintf("Finished processing %s.\n", version))
}

# Reproducibility
sessionInfo()
