# Load libraries
library(terra)
library(ggplot2)
library(ncdf4)
library(here)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)

# --- Paths ---
base_dir <- here::here("DATA")
out_dir <- here::here("OUTPUTS_WarmingLevels")
if(!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)

# --- Load the correct file ---
pr_file <- file.path(base_dir, "Pr", "pre.nc")

# Load Czechia boundary
czech <- ne_countries(scale = "medium", country = "Czechia", returnclass = "sf")

# --- Open and inspect NetCDF ---
nc <- nc_open(pr_file)
print(nc)  # check variable names — we expect 'pre'
pr_array <- ncvar_get(nc, "pre")   # <-- use 'pre' instead of 'pr'
lon <- ncvar_get(nc, "lon")
lat <- ncvar_get(nc, "lat")
time <- ncvar_get(nc, "time")
nc_close(nc)

# --- Convert time to actual years ---
# The file says: units = "days since 1850-1-1 00:00:00"
origin <- as.Date("1850-01-01")
time_dates <- origin + time
years <- as.numeric(format(time_dates, "%Y"))

# --- Select a specific year (e.g., 2015) ---
time_index <- which(years == 2015)
# Average over all 2015 days
pr_2015 <- apply(pr_array[,,time_index], c(1,2), mean, na.rm = TRUE)

# --- Create raster manually ---
r <- rast(t(pr_2015), extent = c(min(lon), max(lon), min(lat), max(lat)), crs = "EPSG:4326")
r <- flip(r, direction = "vertical")

# --- Convert to data frame for ggplot ---
pr_df <- as.data.frame(r, xy = TRUE)
names(pr_df)[3] <- "precip_mm_day"

# Convert Czechia to terra format
czech_terra <- vect(czech)

# Crop and mask the raster to Czechia
r_czech <- crop(r, czech_terra)
r_czech <- mask(r_czech, czech_terra)

# --- Plot ---
p <- ggplot(pr_df) +
  geom_raster(aes(x = x, y = y, fill = precip_mm_day)) +
  scale_fill_viridis_c(name = "Precipitation (mm/day)", option = "C") +
  coord_equal() +
  theme_minimal() +
  labs(title = "Spatial distribution of precipitation (2015)",
       x = "Longitude", y = "Latitude")

plot(r_czech, main = "Precipitation over Czech Republic")
plot(czech_terra, add = TRUE, border = "black", lwd = 1.5)

# --- Save as PDF ---
pdf_file <- file.path(out_dir, "precipitation_2015.pdf")
ggsave(pdf_file, plot = p, width = 8, height = 6)

cat("PDF saved to:", pdf_file, "\n")
