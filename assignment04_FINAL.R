################################################################################
# Assignment 04: Simple Precipitation Downscaling for SSP3-7.0
################################################################################

library(ncdf4)
library(rnaturalearth)
library(sf)
library(terra)

# Create output directory
out_dir <- "OUTPUTS_ASSIGNMENT04"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Use temperature 
nc_file <- "DATA/SSP370/all_SSP370_2015_2100.nc"

# Check file
nc_check <- nc_open(nc_file)
lon_check <- ncvar_get(nc_check, "lon")
lat_check <- ncvar_get(nc_check, "lat")
cat("  Longitude range:", min(lon_check), "to", max(lon_check), "(", length(lon_check), "points )\n")
cat("  Latitude range:", min(lat_check), "to", max(lat_check), "(", length(lat_check), "points )\n")
cat("  Time dimension:", nc_check$dim$time$len, "timesteps\n")
nc_close(nc_check)

# Years for warming thresholds
year_1.5C <- 2078
year_2.0C <- 2089


# Get annual mean for each year 
get_annual_mean <- function(year) {
  nc <- nc_open(nc_file)
  
  # Get dimensions
  lon <- ncvar_get(nc, "lon")
  lat <- ncvar_get(nc, "lat")
  time <- ncvar_get(nc, "time")
  
  # Parse time units 
  time_units <- ncatt_get(nc, "time", "units")$value
  time_origin <- as.Date(sub(".*since ", "", time_units))
  dates <- time_origin + as.numeric(time)
  
  # Find indices for this year
  idx <- which(format(dates, "%Y") == as.character(year))
  cat("  Year", year, "-", length(idx), "days (indices", min(idx), "to", max(idx), ")\n")
  
  if (length(idx) == 0) {
    nc_close(nc)
    stop("No data found for year ", year, ". Date range: ", min(dates), " to ", max(dates))
  }
  
  # Read only this year's data
  nlon <- length(lon)
  nlat <- length(lat)
  year_start <- min(idx)
  year_count <- length(idx)
  
  tas_year <- ncvar_get(nc, "tas", 
                        start = c(1, 1, year_start),
                        count = c(nlon, nlat, year_count))
  nc_close(nc)
  
  # Convert K to C and calculate annual mean
  tas_year <- tas_year - 273.15
  annual <- apply(tas_year, c(1,2), mean, na.rm = TRUE)
  
  r <- rast(nrows = nlat, ncols = nlon,
            xmin = min(lon), xmax = max(lon),
            ymin = min(lat), ymax = max(lat),
            crs = "EPSG:4326")
  
  # Transpose and flip for proper orientation
  values(r) <- as.vector(t(annual[, ncol(annual):1]))
  return(r)
}

r_1.5 <- get_annual_mean(year_1.5C)
r_2.0 <- get_annual_mean(year_2.0C)

cat("\nOriginal global data extent:")
cat("\n  Longitude:", round(ext(r_1.5)[1], 2), "to", round(ext(r_1.5)[2], 2), "°E")
cat("\n  Latitude:", round(ext(r_1.5)[3], 2), "to", round(ext(r_1.5)[4], 2), "°N")
cat("\n  Resolution:", paste(round(res(r_1.5), 4), collapse = " x "), "degrees (~", round(res(r_1.5)[1] * 111, 0), "km )")
cat("\n  Dimensions:", paste(dim(r_1.5)[1:2], collapse = " x "), "cells\n")

# Crop to Czech Republic extent 
cat("\nCropping to Czech Republic (12-19°E, 48-51°N)...\n")
cze_ext <- ext(12, 19, 48, 51)
r_1.5 <- crop(r_1.5, cze_ext)
r_2.0 <- crop(r_2.0, cze_ext)

cat("Cropped CZ extent:")
cat("\n  Longitude:", round(ext(r_1.5)[1], 2), "to", round(ext(r_1.5)[2], 2), "°E")
cat("\n  Latitude:", round(ext(r_1.5)[3], 2), "to", round(ext(r_1.5)[4], 2), "°N")
cat("\n  Dimensions:", paste(dim(r_1.5)[1:2], collapse = " x "), "cells")
cat(" (", round(dim(r_1.5)[1] * res(r_1.5)[1] * 111, 0), "x", round(dim(r_1.5)[2] * res(r_1.5)[2] * 111, 0), "km )\n")

cat("\nData range check:")
cat("\n  1.5°C warming:", round(global(r_1.5, "min", na.rm=TRUE)[,1], 2), "to", 
    round(global(r_1.5, "max", na.rm=TRUE)[,1], 2), "°C")
cat("\n  2.0°C warming:", round(global(r_2.0, "min", na.rm=TRUE)[,1], 2), "to", 
    round(global(r_2.0, "max", na.rm=TRUE)[,1], 2), "°C\n")

cat("\nDownscaling...\n")

# Bilinear
cat("  Bilinear 50->25km\n")
r_1.5_bi25 <- disagg(r_1.5, fact = 2, method = "bilinear")
r_2.0_bi25 <- disagg(r_2.0, fact = 2, method = "bilinear")

cat("  Bilinear 25->12.5km\n")
r_1.5_bi12 <- disagg(r_1.5_bi25, fact = 2, method = "bilinear")
r_2.0_bi12 <- disagg(r_2.0_bi25, fact = 2, method = "bilinear")

# Nearest neighbor 
cat("  Nearest neighbor 50->25km\n")
r_1.5_cb25 <- disagg(r_1.5, fact = 2, method = "near")
r_2.0_cb25 <- disagg(r_2.0, fact = 2, method = "near")

cat("  Nearest neighbor 25->12.5km\n")
r_1.5_cb12 <- disagg(r_1.5_cb25, fact = 2, method = "near")
r_2.0_cb12 <- disagg(r_2.0_cb25, fact = 2, method = "near")

cat("\nCalculating statistics...\n")

# Stats
get_stats <- function(r) {
  s <- global(r, c("min", "max", "mean", "sum"), na.rm = TRUE)
  c(s$min, s$max, s$mean, s$sum)
}

comparison <- data.frame(
  Warming = rep(c("1.5°C", "2.0°C"), each = 5),
  Year = rep(c(year_1.5C, year_2.0C), each = 5),
  Resolution = rep(c("50km", "25km", "25km", "12.5km", "12.5km"), 2),
  Method = rep(c("Original", "Bilinear", "NearestNeighbor", "Bilinear", "NearestNeighbor"), 2)
)

comparison$Min <- round(c(
  get_stats(r_1.5)[1], get_stats(r_1.5_bi25)[1], get_stats(r_1.5_cb25)[1],
  get_stats(r_1.5_bi12)[1], get_stats(r_1.5_cb12)[1],
  get_stats(r_2.0)[1], get_stats(r_2.0_bi25)[1], get_stats(r_2.0_cb25)[1],
  get_stats(r_2.0_bi12)[1], get_stats(r_2.0_cb12)[1]
), 2)

comparison$Max <- round(c(
  get_stats(r_1.5)[2], get_stats(r_1.5_bi25)[2], get_stats(r_1.5_cb25)[2],
  get_stats(r_1.5_bi12)[2], get_stats(r_1.5_cb12)[2],
  get_stats(r_2.0)[2], get_stats(r_2.0_bi25)[2], get_stats(r_2.0_cb25)[2],
  get_stats(r_2.0_bi12)[2], get_stats(r_2.0_cb12)[2]
), 2)

comparison$Mean <- round(c(
  get_stats(r_1.5)[3], get_stats(r_1.5_bi25)[3], get_stats(r_1.5_cb25)[3],
  get_stats(r_1.5_bi12)[3], get_stats(r_1.5_cb12)[3],
  get_stats(r_2.0)[3], get_stats(r_2.0_bi25)[3], get_stats(r_2.0_cb25)[3],
  get_stats(r_2.0_bi12)[3], get_stats(r_2.0_cb12)[3]
), 2)

comparison$Sum <- round(c(
  get_stats(r_1.5)[4], get_stats(r_1.5_bi25)[4], get_stats(r_1.5_cb25)[4],
  get_stats(r_1.5_bi12)[4], get_stats(r_1.5_cb12)[4],
  get_stats(r_2.0)[4], get_stats(r_2.0_bi25)[4], get_stats(r_2.0_cb25)[4],
  get_stats(r_2.0_bi12)[4], get_stats(r_2.0_cb12)[4]
), 0)

print(comparison)
write.csv(comparison, file.path(out_dir, "comparison_table.csv"), row.names = FALSE)

cat("\nCreating maps...\n")

# Czech Republic boundary
cze_boundary <- NULL
if (requireNamespace("rnaturalearth", quietly = TRUE)) {
  library(rnaturalearth)
  cze_boundary <- ne_countries(country = "czechia", scale = "medium", returnclass = "sf")
  cze_boundary <- st_transform(cze_boundary, "EPSG:4326")
} else if (requireNamespace("geodata", quietly = TRUE)) {
  library(geodata)
  cze_boundary <- gadm(country = "CZE", level = 0, path = tempdir())
}

# Plotting function with CZE boundary
plot_with_boundary <- function(r, title) {
  plot(r, main = title, col = terrain.colors(50), axes = TRUE)
  if (!is.null(cze_boundary)) {
    plot(st_geometry(cze_boundary), add = TRUE, border = "black", lwd = 2)
  }
}

# Maps - 1.5°C warming
png(file.path(out_dir, "maps_1.5C.png"), width = 1800, height = 1200, res = 150)
par(mfrow = c(2, 3), mar = c(3,3,3,6))
plot_with_boundary(r_1.5, sprintf("Original 50km\n1.5°C warming (%d)", year_1.5C))
plot_with_boundary(r_1.5_bi25, "Bilinear 25km")
plot_with_boundary(r_1.5_cb25, "Nearest Neighbor 25km")
plot_with_boundary(r_1.5_bi12, "Bilinear 12.5km")
plot_with_boundary(r_1.5_cb12, "Nearest Neighbor 12.5km")
plot(r_1.5_bi12 - r_1.5_cb12, main = "Diff (Bilinear-NN)\n°C", col = hcl.colors(50, "Blue-Red 3"), axes = TRUE)
if (!is.null(cze_boundary)) plot(st_geometry(cze_boundary), add = TRUE, border = "black", lwd = 2)
dev.off()

# Maps - 2.0°C warming
png(file.path(out_dir, "maps_2.0C.png"), width = 1800, height = 1200, res = 150)
par(mfrow = c(2, 3), mar = c(3,3,3,6))
plot_with_boundary(r_2.0, sprintf("Original 50km\n2.0°C warming (%d)", year_2.0C))
plot_with_boundary(r_2.0_bi25, "Bilinear 25km")
plot_with_boundary(r_2.0_cb25, "Nearest Neighbor 25km")
plot_with_boundary(r_2.0_bi12, "Bilinear 12.5km")
plot_with_boundary(r_2.0_cb12, "Nearest Neighbor 12.5km")
plot(r_2.0_bi12 - r_2.0_cb12, main = "Diff (Bilinear-NN)\n°C", col = hcl.colors(50, "Blue-Red 3"), axes = TRUE)
if (!is.null(cze_boundary)) plot(st_geometry(cze_boundary), add = TRUE, border = "black", lwd = 2)
dev.off()

cat("\n=== DONE ===\n")
cat("Files saved to:", out_dir, "\n")

