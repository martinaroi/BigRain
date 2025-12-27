# install if needed
#install.packages(c("ncdf4","terra","sf","rnaturalearth","rnaturalearthdata",
 #             "ggplot2","viridis","magick","dplyr","lubridate","rmarkdown", "gganimate","tidyr"))

# load
library(ncdf4)
library(terra)
library(sf)
library(rnaturalearth)
library(viridis)
library(magick)
library(dplyr)
library(lubridate)
library(rmarkdown)
library(here)
library(ggplot2)
library(gganimate)
library(tidyr)

base_dir <- here::here("DATA")
out_dir  <- here::here("OUTPUTS_WarmingLevels")

pr_files <- list(
  ssp126 = file.path(base_dir, "Pr","cze", "ssp126_pr_cze_annual_2015_2100_mm_per_year.nc"),
  ssp370 = file.path(base_dir, "Pr", "cze", "ssp370_pr_cze_annual_2015_2100_mm_per_year.nc"),
  ssp585 = file.path(base_dir, "Pr", "cze", "ssp585_pr_cze_annual_2015_2100_mm_per_year.nc")
)

tas_files <- list(
  ssp126 = file.path(base_dir, "SSP126", "all_SSP126_2015_2100.nc"),
  ssp370 = file.path(base_dir, "SSP370", "all_SSP370_2015_2100.nc"),
  ssp585 = file.path(base_dir, "SSP585", "all_SSP585_2015_2100.nc")
)

# --- Load and summarize data for animation ---
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

precip_df <- lapply(names(pr_files), function(scn) {
  r <- terra::rast(pr_files[[scn]])
  tvals <- tryCatch(terra::time(r), error = function(e) NULL)
  yrs <- if (is.null(tvals)) seq(2015L, by = 1L, length.out = terra::nlyr(r)) else as.integer(format(as.Date(tvals), "%Y"))
  m <- terra::global(r, fun = "mean", na.rm = TRUE)[, 1]
  data.frame(year = yrs, precip = as.numeric(m), scenario = scn, stringsAsFactors = FALSE)
}) |>
  dplyr::bind_rows() |>
  dplyr::filter(!is.na(year)) |>
  dplyr::arrange(scenario, year)

# Focus on 2025–2100 
precip_df <- precip_df |>
  dplyr::filter(year >= 2025, year <= 2100) |>
  dplyr::group_by(scenario) |>
  dplyr::filter(dplyr::n() >= 2) |>
  dplyr::ungroup()

#  labels for scenarios
lab_map <- c(ssp126 = "SSP1-2.6", ssp370 = "SSP3-7.0", ssp585 = "SSP5-8.5")
precip_df$scenario <- lab_map[precip_df$scenario]
precip_df$scenario <- factor(precip_df$scenario, levels = c("SSP1-2.6","SSP3-7.0","SSP5-8.5"))

# Animated line plot
p <- ggplot(precip_df, aes(x = year, y = precip, color = scenario, group = scenario)) +
  geom_line(size = 1.1, alpha = 0.9) +
  geom_point(size = 2) +
  scale_color_viridis_d(option = "C") +
  labs(
    title = "Precipitation Evolution in Czech Republic (2025–2100)",
    x = "Year",
    y = "Precipitation (mm/year)",
    color = "Scenario"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom") +
  transition_reveal(along = year)

# Save animation as GIF (magick requires fps be a factor of 100; use 10)
gif_path <- file.path(out_dir, "precipitation_evolution_CzechRepublic.gif")
anim <- animate(p, nframes = length(unique(precip_df$year)), fps = 10, width = 800, height = 500)
anim_save(gif_path, animation = anim)
message("Saved GIF: ", gif_path)# Define target years

years <- list(
  current = 2025,
  `+1.5C` = c(ssp126=2035, ssp370=2028, ssp585=2025),
  `+2C`   = c(ssp126=2055, ssp370=2040, ssp585=2035),
  `+2.5C` = c(ssp126=2070, ssp370=2050, ssp585=2040),
  `+3C`   = c(ssp126=2090, ssp370=2060, ssp585=2050)
)

# Function to load only a single year's layer
load_precip <- function(file, year) {
  r <- rast(file)
  t <- time(r)
  idx <- which(format(t, "%Y") == as.character(year))
  subset(r, idx)
}

# Downscale function 
downscale_raster <- function(r, factor=2){
aggregate(r, fact=factor, fun=mean)
}

make_comparison_gifs_by_temp <- function(factor=2){
  for(temp_label in c("+1.5C","+2C","+2.5C","+3C")){
    maps <- list()
    for(scenario in names(pr_files)){
      year <- years[[temp_label]][scenario]
      cat("Processing", scenario, year, "\n")
      
      r <- load_precip(pr_files[[scenario]], year)
      r <- downscale_raster(r, factor=factor)
      
      df <- as.data.frame(r, xy=TRUE)
      colnames(df) <- c("lon","lat","pr")
      
      p <- ggplot(df) +
        geom_raster(aes(lon, lat, fill=pr)) +
        scale_fill_viridis_c(option="C", name="Precip (mm/year)") +
        labs(title=paste(scenario, temp_label, "(", year, ")")) +
        coord_equal() +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust=0.5, size=12),
          axis.text = element_text(size=8)
        )
      
      png_file <- paste0("map_", scenario, "_", temp_label, ".png")
      ggsave(png_file, p, width=6, height=4)
      maps[[scenario]] <- png_file
    }
    
    # Create one GIF per temperature increase
    imgs <- image_read(unlist(maps))
    anim <- image_animate(image_join(imgs), fps=1)
    out_gif <- paste0("precipitation_", gsub("\\+","plus", temp_label), ".gif")
    image_write(anim, out_gif)
    cat("GIF saved as", out_gif, "\n")
  }
}

# Run the function
make_comparison_gifs_by_temp(factor=2)

# SSP5-8.5 precipitation file
pr_file <- file.path(base_dir, "Pr","cze","ssp585_pr_cze_annual_2015_2100_mm_per_year.nc")

# Teplá catchment extent
lon_min <- 12.5
lon_max <- 13.5
lat_min <- 49.8
lat_max <- 50.2

# Load full SSP5-8.5 raster
r <- rast(pr_file)

# Crop to Teplá catchment
r_tepla <- crop(r, ext(lon_min, lon_max, lat_min, lat_max))

# Optional: downscale by factor of 2 to reduce grid size for plotting
r_tepla_small <- aggregate(r_tepla, fact = 2, fun = mean)

years_ssp585 <- c(
  `+1.5C` = 2025,
  `+2C`   = 2035,
  `+2.5C` = 2040,
  `+3C`   = 2050
)

year_to_layer <- function(year) {
  year - 2015 + 1
}

# Compute global min/max across all years for Teplá
all_values <- values(r_tepla_small)
all_values <- all_values[!is.na(all_values)]
pr_min <- min(all_values)
pr_max <- max(all_values)

plot_precip <- function(year, threshold_label) {
  idx <- year_to_layer(year)
  r_layer <- r_tepla_small[[idx]]
  df <- as.data.frame(r_layer, xy = TRUE)
  colnames(df) <- c("lon", "lat", "pr")
  df <- df[complete.cases(df), ]  # remove NAs
  
  p <- ggplot(df) +
    geom_raster(aes(lon, lat, fill = pr)) +
    scale_fill_viridis_c(option = "C", name = "Precip (mm/year)", limits = c(pr_min, pr_max)) +
    coord_equal() +
    theme_minimal() +
    labs(title = paste0("SSP5-8.5 Precipitation (Teplá) - ", threshold_label, " (", year, ")"))
  
  png_file <- file.path(out_dir, paste0("Tepla_", threshold_label, ".png"))
  ggsave(png_file, p, width = 6, height = 4)
  return(png_file)
}

png_files <- sapply(names(years_ssp585), function(th) {
  plot_precip(years_ssp585[th], th)
})

# Create animated GIF
imgs <- image_read(png_files)
anim <- image_animate(image_join(imgs), fps = 1)
gif_file <- file.path(out_dir, "Tepla_precip_SSP585.gif")
image_write(anim, gif_file)

cat("GIF saved to:", gif_file, "\n")

# Compute mean precipitation for each warming level in Teplá catchment
tepla_means <- data.frame(
  Warming_Level = names(years_ssp585),
  Year = unname(years_ssp585),
  Mean_Precip_mm = sapply(unname(years_ssp585), function(year) {
    idx <- year_to_layer(year)
    r_layer <- r_tepla_small[[idx]]
    mean(values(r_layer), na.rm = TRUE)
  })
)

# Save to CSV
write.csv(tepla_means, file.path(out_dir, "Tepla_precip_means_SSP585.csv"), row.names = FALSE)

# Print the table
print(tepla_means)

# === Table: Mean precipitation for each SSP and warming level ===
precip_summary <- precip_table %>%
  mutate(
    Mean_Precip_mm = round(Mean_Precip_mm, 1)
  ) %>%
  arrange(Warming_Level, Scenario)

print(precip_summary)

# Optionally, write to CSV for record
write.csv(precip_summary, file.path(out_dir, "CzechRepublic_Precip_Summary.csv"), row.names = FALSE)
