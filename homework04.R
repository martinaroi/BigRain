library(terra)
library(ncdf4)
library(dplyr)
library(lubridate)

# --- Historical baseline (1850–2014) ---
hist_file <- file.path("DATA", "ismip", "all_tas_cze_1850_2014.nc")
nc_hist <- nc_open(hist_file)
tas_hist <- ncvar_get(nc_hist, "tas") - 273.15
time_hist <- ncvar_get(nc_hist, "time")
time_units <- ncatt_get(nc_hist, "time", "units")$value
origin_date <- sub("days since ", "", time_units)
dates_hist <- as.Date(time_hist, origin = origin_date)
years_hist <- year(dates_hist)
nc_close(nc_hist)

# Baseline mean 1850–1900
baseline_temp <- mean(hist_annual$tasC[hist_annual$year >= 1995 & hist_annual$year <= 2014], na.rm = TRUE)

# --- SSP3-7.0 (2015–2100) ---
tas_file <- file.path("DATA", "SSP370", "all_SSP370_2015_2100.nc")
tas_ssp <- rast(tas_file, subds = "tas")

# Extract yearly means efficiently (CZ region already)
tvals <- terra::time(tas_ssp)
years_ssp <- as.integer(format(as.Date(tvals), "%Y"))

# Compute yearly mean over all pixels
yearly_means <- tapply(1:nlyr(tas_ssp), years_ssp, function(i) {
  mean(global(tas_ssp[[i]], mean, na.rm = TRUE)[[1]]) - 273.15
})

ssp_df <- data.frame(
  year = as.integer(names(yearly_means)),
  tas_mean = as.numeric(yearly_means)
)

# --- Threshold crossing ---
thresholds <- c(1.5, 2.0)

threshold_years <- sapply(thresholds, function(deltaT) {
  idx <- which(ssp_df$tas_mean >= (baseline_temp + deltaT))[1]
  if (is.na(idx)) return(NA)
  ssp_df$year[idx]
})

print(threshold_years)

# Check what range of values we actually have
cat("Baseline (1850–1900) mean:", baseline_temp, "\n")
cat("Min–max SSP370 mean (2015–2100):", min(ssp_df$tas_mean, na.rm = TRUE), "to", 
    max(ssp_df$tas_mean, na.rm = TRUE), "\n")

# Optional: see a few last years
tail(ssp_df, 10)


