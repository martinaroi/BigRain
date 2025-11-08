library(terra)
library(ncdf4)
library(dplyr)
library(lubridate)

# --- Historical baseline (1850–1900) ---
hist_file <- file.path("DATA", "ismip", "all_tas_cze_1850_2014.nc")
nc_hist <- nc_open(hist_file)

# Extract historical temperature (1x1x60265 array)
tas_hist <- ncvar_get(nc_hist, "tas")  
tas_hist_vec <- as.vector(tas_hist) - 273.15  # K -> °C

# Time handling
time_hist <- ncvar_get(nc_hist, "time")
time_units <- ncatt_get(nc_hist, "time", "units")$value
dates_hist <- as.Date(time_hist, origin = sub("days since ", "", time_units))
years_hist <- year(dates_hist)

# Compute annual mean per year
hist_annual <- data.frame(year = years_hist, tas = tas_hist_vec) %>%
  group_by(year) %>%
  summarise(tasC = mean(tas, na.rm = TRUE), .groups = "drop")

# Baseline 1850–1900
baseline_temp <- mean(hist_annual$tasC[hist_annual$year >= 1850 & hist_annual$year <= 1900])

# --- SSP3-7.0 future data ---
tas_file <- file.path("DATA", "SSP370", "all_SSP370_2015_2100.nc")
tas_ssp <- rast(tas_file, subds = "tas")  # SpatRaster (CZ only)

# Time info
tvals <- terra::time(tas_ssp)
years_ssp <- as.integer(format(as.Date(tvals), "%Y"))

# Annual mean per year (CZ only)
ssp_annual <- data.frame(year = years_ssp, 
                         tasC = global(tas_ssp, "mean", na.rm = TRUE)[,1] - 273.15)

# Combine historical + SSP
combined <- bind_rows(hist_annual, ssp_annual) %>%
  arrange(year) %>%
  mutate(anomaly = tasC - baseline_temp)

# --- Find first year exceeding thresholds ---
thresholds <- c(1.5, 2.0)
first_years <- sapply(thresholds, function(t) {
  idx <- which(combined$year >= 2015 & combined$anomaly >= t)[1]
  combined$year[idx]
})
names(first_years) <- paste0(thresholds, "C")

print(first_years)

