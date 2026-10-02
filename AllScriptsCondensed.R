# STATUS AND TRENDS MAPS *******************************************************
# LOADING THE CONFIGURATION AND DATA
# Load centralized paths and secrets safely
source("paths_config.R")

# load all of the packages
library(ebirdst)
library(fields)
library(gridExtra)
library(mccf1)
library(ranger)
library(scam)
library(sf)
library(terra)
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthhires)


#download YBCU data
# request key every 6 months from https://ebird.org/st/request
set_ebirdst_access_key(EBIRD_KEY, overwrite = TRUE)
path <- ebirdst_download_status(species = "yellow-billed cuckoo", force = TRUE)

# load relative abundance raster stack with 52 layers, one for each week, 13 through 27
abd <- load_raster(species = "yellow-billed cuckoo", product = "abundance", resolution = "3km")

# the date that each raster layer corresponds to is stored within the labels
weeks <- names(abd)
print(weeks)

# select a week, 13 through 27, change this and repeat code to get each map
abd <- abd[[20]]

# ZOOMING IN/MAKING THE EXTENT
#make a boundary polygon for the border states we're studying
# project to same coordinate reference system as the raster data
us_border_states <- ne_states(iso_a2 = "US", returnclass = "sf") %>% filter(
  region == "South") %>% filter(postal != "SC") %>% filter(postal != "NC") %>% filter(
    postal != "VA") %>% filter(postal != "DC") %>% filter(postal != "MD") %>% filter(
      postal != "DE") %>% filter(postal != "OK") %>% filter(postal != "KY") %>% filter(
        postal != "TN") %>% filter(postal != "WV") %>% st_transform(st_crs(abd))

# crop data to the border states we're studying
abd_us_border_states <- crop(abd, us_border_states)

#FIXING THE PROJECTION
# load the mapping parameters to fix distortion
fac_parameters <- load_fac_map_parameters(species = "yellow-billed cuckoo")
crs <- fac_parameters$custom_projection

# transform to the custom projection using nearest neighbor resampling
abd_projected <- project(abd_us_border_states, crs, method = "near")

#EVENLY DISTRIBUTING THE ABUNDANCE BINS WITH CELLS ON MAP
# quantiles of non-zero values
v <- values(abd_projected)
v <- v[!is.na(v) & v > 0]
bins <- quantile(v, seq(0, 1, by = 0.1))

# add a bin for 0
bins <- c(0, bins)

# status and trends palette
pal <- abundance_palette(length(bins) - 1)

# add a color for zero
pal <- c(NA, pal)

# map using the quantile bins
plot(abd_projected, breaks = bins, col = pal, axes = FALSE, legend = FALSE)

#MAKING A BASEMAP BEHIND THE DATA
# natural earth boundaries
countries <- ne_countries(returnclass = "sf") %>% 
  st_geometry() %>% 
  st_transform(crs)
states <- ne_states(iso_a2 = "US", returnclass = "sf") %>% 
  st_geometry() %>% 
  st_transform(crs)

# add basemap
plot(countries, col = "grey", border = "#888888", add = TRUE, legend = FALSE)
plot(states, col = "grey", border = "#888888", add = TRUE, legend = FALSE)

# add data
plot(abd_projected, 
     breaks = bins, col = pal, 
     axes = FALSE, legend = FALSE, add = TRUE)

# add title, change for each week
title(sub = "May 17th", line = -2, font.sub = 2, family = "sans") 

# add legend
# label the bottom, middle, and top (this isn't working??)
labels <- quantile(bins, c(0, 0.5, 1))
plot(abd_projected, zlim = c(0, 1), legend.only = TRUE, 
     col = pal, breaks = seq(0, 1, length.out = length(bins)), 
     legend.shrink = 0.9, legend.width = 2,  
     axis.args = list(at = seq(0, 1, length.out = length(labels)), 
                      labels = signif(labels, 3),
                      col.axis = "black", fg = NA,
                      cex.axis = 0.9, lwd.ticks = 0,
                      line = -0.5))


# EBD HABITAT SELECTION AND MAP ************************************************
# install all the packages
library(auk)
library(ggplot2)
library(gridExtra)
library(lubridate)
library(readr)
library(sf)
library(dplyr)
library(rnaturalearth)
library(exactextractr)
library(landscapemetrics)
library(stringr)
library(terra)
library(tidyr)
library(units)
library(viridis)
library(luna)
library(rnaturalearthhires)
library(geodata)

# PREPARE GIS DATA
# file to save spatial data
gpkg_file <- file.path(DROPBOX_DIR, "data/relative_abundance_environmental_data/gis-data.gpkg")
dir.create(dirname(gpkg_file), showWarnings = FALSE, recursive = TRUE)

# political boundaries
# land border with lakes removed
ne_land <- ne_download(scale = 50, category = "cultural",
                       type = "admin_0_countries_lakes",
                       returnclass = "sf") %>% 
  filter(CONTINENT %in% c("North America", "South America")) %>% 
  st_set_precision(1e6) %>% 
  st_union()

# country boundaries
ne_countries <- ne_download(scale = 50, category = "cultural",
                            type = "admin_0_countries_lakes",
                            returnclass = "sf") %>% 
  select(country = ADMIN, country_code = ISO_A2)

# state boundaries for united states
ne_states <- ne_download(scale = 50, category = "cultural",
                         type = "admin_1_states_provinces",
                         returnclass = "sf") %>%  
  filter(iso_a2 == "US") %>%  
  select(state = name, state_code = iso_3166_2)

# country lines
# downloaded globally then filtered to north america with st_intersect()
ne_country_lines <- ne_download(scale = 50, category = "cultural",
                                type = "admin_0_boundary_lines_land",
                                returnclass = "sf") %>%  
  st_geometry()
lines_on_land <- st_intersects(ne_country_lines, ne_land, sparse = FALSE) %>% 
  as.logical()
ne_country_lines <- ne_country_lines[lines_on_land]
# states, north america
ne_state_lines <- ne_download(scale = 50, category = "cultural",
                              type = "admin_1_states_provinces_lines",
                              returnclass = "sf") %>% 
  filter(ADM0_A3 %in% c("USA", "CAN")) %>% 
  mutate(iso_a2 = recode(ADM0_A3, USA = "US", CAN = "CAN")) %>%  
  select(country = ADM0_NAME, country_code = iso_a2)

# save all layers to a geopackage
unlink(gpkg_file)
write_sf(ne_land, gpkg_file, "ne_land")
write_sf(ne_countries, gpkg_file, "ne_countries")
write_sf(ne_states, gpkg_file, "ne_states")
write_sf(ne_country_lines, gpkg_file, "ne_country_lines")
write_sf(ne_state_lines, gpkg_file, "ne_state_lines")

# CLEAN UP THE DATA
# import the SED and EBD
f_sed <- file.path(DROPBOX_DIR, "data_raw/ebd_US_yebcuc_202303_202307_smp_relSep-2023/ebd_US_yebcuc_202303_202307_smp_relSep-2023_sampling.txt")
checklists <- read_sampling(f_sed)
f_ebd <- file.path(DROPBOX_DIR, "data_raw/ebd_US_yebcuc_202303_202307_smp_relSep-2023/ebd_US_yebcuc_202303_202307_smp_relSep-2023.txt")
observations <- read_ebd(f_ebd)

# filter the SED and EBD to only include complete checklists
checklists <- checklists %>% 
  filter(all_species_reported)
observations <- observations %>% 
  filter(all_species_reported)

# convert checklist locations to points geometries
checklists_sf <- checklists %>%  
  dplyr:::select(checklist_id, latitude, longitude) %>%  
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

# boundary of study region, buffered by 1 km
study_region_buffered <- read_sf(file.path(DROPBOX_DIR, "data/relative_abundance_environmental_data/gis-data.gpkg", layer = "ne_states")) %>% 
  filter(state_code == "US-GA" | state_code == "US-FL" | state_code == "US-AL" | state_code == "US-MS" | state_code == "US-LA" | state_code == "US-TX" | state_code == "US-AR") %>% 
  st_transform(crs = st_crs(checklists_sf)) %>% 
  st_buffer(dist = 1000)

# spatially subset the checklists to those in the study region
in_region <- checklists_sf[study_region_buffered, ]

# join to checklists and observations to remove checklists outside region
checklists <- semi_join(checklists, in_region, by = "checklist_id")
observations <- semi_join(observations, in_region, by = "checklist_id")

# remove observations without matching checklists
observations <- semi_join(observations, checklists, by = "checklist_id")

# combine checklist and observation data together to create zero-filled data
zf <- auk_zerofill(observations, checklists, collapse = TRUE)

# create function to convert time to hours since midnight
time_to_decimal <- function(x){
  x <- hms(x, quiet = TRUE)
  hour(x) + minute(x) / 60 + second(x) / 3600
}

# clean up time, stationary distance traveled, observation count, and create speed variable
zf <- zf %>% 
  mutate(
    # convert X to NA
    observation_count = as.integer(observation_count),
    # effort_distance_km to 0 for non-travelling counts
    effort_distance_km = if_else(protocol_type == "Stationary", 
                                 0, effort_distance_km),
    # convert duration to hours
    effort_hours = duration_minutes / 60,
    # speed km/h
    effort_speed_kmph = effort_distance_km / effort_hours,
    # convert time to decimal hours since midnight
    hours_of_day = time_to_decimal(time_observations_started),
    # split date into year and day of year
    year = year(observation_date),
    day_of_year = yday(observation_date)
  )

# filter effort data to add a more consistent structure
zf_filtered <- zf %>% 
  filter(protocol_type %in% c("Stationary", "Traveling"),
         effort_hours <= 6,
         effort_distance_km <= 10,
         effort_speed_kmph <= 100,
         number_observers <= 10) %>% filter(state == "Florida" | state == "Georgia" | state == "Alabama" | state == "Mississippi" | state == "Louisiana" | state == "Arkansas" | state == "Texas")

# test-train split for the model (0.2 vs 0.8)
zf_filtered$type <- if_else(runif(nrow(zf_filtered)) <= 0.8, "train", "test")
# confirm the proportion in each set is correct
table(zf_filtered$type) / nrow(zf_filtered)

# get rid of redundant and unnecessary variables
checklists <- zf_filtered %>% 
  dplyr:::select(checklist_id, observer_id,
                 observation_count, species_observed, 
                 state_code, locality_id, latitude, longitude,
                 protocol_type, all_species_reported,
                 observation_date, year, day_of_year,
                 hours_of_day, 
                 effort_hours, effort_distance_km, effort_speed_kmph,
                 number_observers)
write.csv(checklists, file.path(DROPBOX_DIR, "data/relative_abundance_zero_filtered_data/relative_abundance_zero_filtered_data", na = ""))

# install all the packages
library(auk)
library(ggplot2)
library(gridExtra)
library(lubridate)
library(readr)
library(sf)
library(dplyr)
library(rnaturalearth)
library(exactextractr)
library(landscapemetrics)
library(stringr)
library(terra)
library(tidyr)
library(units)
library(viridis)
library(luna)
library(rnaturalearthhires)
library(geodata)
library(raster)

# GET ENVIRONMENTAL VARIABLES FOR MAPS (let's use landcover first)
# set the date
product <- "MCD12Q1"
start <- "2022-03-01"
end <- "2022-07-31"

#define the area of interest
south <- geodata::gadm("United States", level = 1, path = ".")
i <- south$NAME_1 == "Florida" | south$NAME_1 == "Georgia" | south$NAME_1 == "Alabama" | south$NAME_1 == "Mississippi" | south$NAME_1 == "Louisiana" | south$NAME_1 == "Texas" | south$NAME_1 == "Arkansas"
aoi <- south[i, ]

# find what MODIS data is available
mf <- luna::getNASA(product, start_date = start, end_date = end, aoi = aoi, download = FALSE)

# set up a directory
datadir <- file.path(DROPBOX_DIR, "data_raw/modis_landcover_data_tiles_01_2024")
dir.create(datadir, showWarnings = FALSE)

# download the data
mf <- luna::getNASA(product, start, end, aoi = aoi, download = TRUE, path = datadir, username = NASA_USER, password = NASA_PASS)

# https://rspatial.org/modis/3-explore.html has further instructions for visualizing the data
# quality control
# define quality assessment bands to reject
from <- c(1, 3, 11, 14)
to <- c(2, 3, 11, 14)
reject <- c("01,10", "1", "1", "1")
qa_bits <- cbind(from, to, reject)

# stitch MODIS data together
hdf1 <- file.path(datadir, "MCD12Q1.A2022001.h08v05.061.2023243090807.hdf")
hdf1 <- rast(hdf1)
hdf2 <- file.path(datadir, "MCD12Q1.A2022001.h08v06.061.2023243092049.hdf")
hdf2 <- rast(hdf2)
hdf3 <- file.path(datadir, "MCD12Q1.A2022001.h09v05.061.2023243094009.hdf")
hdf3 <- rast(hdf3)
hdf4 <- file.path(datadir, "MCD12Q1.A2022001.h09v06.061.2023243075235.hdf")
hdf4 <- rast(hdf4)
hdf5 <- file.path(datadir, "MCD12Q1.A2022001.h10v05.061.2023243103913.hdf")
hdf5 <- rast(hdf5)
hdf6 <- file.path(datadir, "MCD12Q1.A2022001.h10v06.061.2023243073808.hdf")
hdf6 <- rast(hdf6)
hdf7 <- file.path(datadir, "MCD12Q1.A2022001.h11v05.061.2023243093838.hdf")
hdf7 <- rast(hdf7)
hdf8 <- file.path(datadir, "MCD12Q1.A2022001.h11v06.061.2023243063834.hdf")
hdf8 <- rast(hdf8)
r <- mosaic(hdf1, hdf2, hdf3, hdf4, hdf5, hdf6, hdf7, hdf8)
plot(r)

# generate the quality mask using band 12 sur_refl_state_500m
qc <- r[[12]]
plot(qc, main = "Quality")
quality_mask <- modis_mask(qc, 16, qa_bits)
plot(quality_mask, main = "Quality mask")

# apply to all bands
rmask <- mask(r, quality_mask)
plot(rmask)

# change vector boundary coordinate system
# rmask <- rast(rmask)
prj <- crs(rmask)
poly <- project(aoi, prj)

# crop the image using the transformed vector boundaries
# rcrop <- crop(rmask, poly)

# plot MODIS and add the boundary
plot(rmask, y = 2)
lines(poly, col="deepskyblue4")

# save the data!
writeRaster(rmask, filename = "processedMCD12q1.tif")


# load the landcover data
landcover <- rast(file.path(DROPBOX_DIR, "data_raw/modis_landcover_data_tiles_01_2024/processedMCD12q1.tif"))
crs(landcover)

# look up table of the landcover classes
lc_classes <- read_csv(file.path(DROPBOX_DIR, "data_raw/modis_landcover_data_tiles_01_2024/mcd12q1_umd_classes.csv"))

# ebird checklist locations
checklists <- read_csv(file.path(DROPBOX_DIR, "data/relative_abundance_zero_filtered_data/relative_abundance_zero_filtered_data"))

# generate circular neighborhoods for all checklists
checklists_sf <- checklists %>% distinct(locality_id, latitude, longitude) %>% st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
buffers <- st_buffer(checklists_sf, dist = set_units(1.5, km))

# crop and mask landcover layer according to checklist layer to the circular neighborhood
# then calculate percent landcover and edge density within each neighborhood
# this may take 30 minutes to run!
lsm <- NULL
for (i in seq_len(nrow(buffers))) {
  buffer_i <- st_transform(buffers[i, ], crs = crs(landcover))
  
  # crop and mask landcover raster so all values outside buffer are missing
  lsm[[i]] <- crop(landcover, buffer_i) %>%  
    mask(buffer_i) %>%  
    # calculate landscape metrics
    calculate_lsm(level = "class", metric = c("pland", "ed"))  %>% 
    # add variables to uniquely identify each point
    mutate(locality_id = buffer_i$locality_id)  %>% 
    dplyr:::select(locality_id, class, metric, value)
}
lsm <- bind_rows(lsm)

# transform the data so one row per location and environmental data appears as columns
lsm_wide <- lsm %>%  
  # fill missing classes with zeros
  complete(nesting(locality_id),
           class = lc_classes$class,
           metric = c("ed", "pland"),
           fill = list(value = 0)) %>%  
  # bring in more descriptive names
  inner_join(dplyr:::select(lc_classes, class, label), by = "class") %>%  
  # transform from long to wide format
  pivot_wider(values_from = value,
              names_from = c(class, label, metric),
              names_glue = "{metric}_c{str_pad(class, 2, pad = '0')}_{label}",
              names_sort = TRUE) %>%  
  arrange(locality_id)

# elevation raster
elevation <- rast(file.path(DROPBOX_DIR, "data_raw/earthenv_topography_data/elevation_1KMmd_GMTEDmd.tif"))

# mean and standard deviation within each circular neighborhood
elev_buffer <- exact_extract(elevation, buffers, fun = c("mean", "stdev"),
                             progress = FALSE) %>%  
  # add variables to uniquely identify each point
  mutate(locality_id = buffers$locality_id) %>%  
  dplyr:::select(locality_id, 
                 elevation_mean = mean,
                 elevation_sd = stdev)

# combine elevation and landcover
env_variables <- inner_join(elev_buffer, lsm_wide,
                            by = "locality_id")

# attach and expand to checklists
env_variables <- checklists %>%  
  dplyr:::select(checklist_id, locality_id) %>%  
  inner_join(env_variables, by = "locality_id") %>%  
  dplyr:::select(-locality_id)

# save to csv, dropping any rows with missing variables
write_csv(drop_na(env_variables), file.path(DROPBOX_DIR, "data/relative_abundance_environmental_data/environmental-variables_checklists_spring_us-south.csv"), na = "")

# find the center coordinates of the checklists data
checklists_sf %>%  
  st_union() %>%  
  st_centroid() %>%  
  st_coordinates() %>%  
  round(1)

# lambert's azimuthal equal area projection for the south
laea_crs <- st_crs("+proj=laea +lat_0=30.9 +lon_0=-89.6")

study_region <- read_sf(file.path(DROPBOX_DIR, "data/relative_abundance_environmental_data/gis-data.gpkg"), layer = "ne_states") %>% filter(state_code == "US-GA" | state_code == "US-FL" | state_code == "US-AL" | state_code == "US-MS" | state_code == "US-LA" | state_code == "US-TX" | state_code == "US-AR") %>% st_transform(crs = laea_crs)

# create a raster template covering the region with a 3km resolution
r <- rast(study_region, res = c(3000, 3000))

# fill the raster with 1s inside the study region
r <- rasterize(study_region, r, values = 1) %>%  
  setNames("study_region")

# save for later use
r <- writeRaster(r, file.path(DROPBOX_DIR, "data/relative_abundance_environmental_data/prediction-grid_us-south.tif"),
                 overwrite = TRUE,
                 gdal = "COMPRESS=DEFLATE")

# generate neighborhoods for the prediction grid cell centers
buffers_pg <- as.data.frame(r, cells = TRUE, xy = TRUE) %>%  
  dplyr:::select(cell_id = cell, x, y) %>%  
  st_as_sf(coords = c("x", "y"), crs = laea_crs, remove = FALSE) %>%  
  st_transform(crs = 4326) %>%  
  st_buffer(set_units(3, "km"))

# estimate landscape metrics for each cell in the prediction grid
lsm_pg <- NULL
for (i in seq_len(nrow(buffers_pg))) {
  buffer_i <- st_transform(buffers_pg[i, ], crs = crs(landcover))
  
  # crop and mask landcover raster so all values outside buffer are missing
  lsm_pg[[i]] <- crop(landcover, buffer_i) %>%  
    mask(buffer_i) %>%  
    # calcualte landscape metrics
    calculate_lsm(level = "class", metric = c("pland", "ed")) %>%  
    # add variable to uniquely identify each point
    mutate(cell_id = buffer_i$cell_id) %>%  
    dplyr:::select(cell_id, class, metric, value)
}
lsm_pg <- bind_rows(lsm_pg)

# transform to wide format
lsm_wide_pg <- lsm_pg %>%  
  # fill missing classes with zeros
  complete(cell_id,
           class = lc_classes$class,
           metric = c("ed", "pland"),
           fill = list(value = 0)) %>%  
  # bring in more descriptive names
  inner_join(dplyr:::select(lc_classes, class, label), by = "class") %>%  
  # transform from long to wide format
  pivot_wider(values_from = value,
              names_from = c(class, label, metric),
              names_glue = "{metric}_c{str_pad(class, 2, pad = '0')}_{label}",
              names_sort = TRUE, values_fill = 0) %>% arrange(cell_id)

# mean and standard deviation of elevation
elev_buffer_pg <- exact_extract(elevation, buffers_pg, 
                                fun = c("mean", "stdev"),
                                progress = FALSE) %>%  
  # add variables to uniquely identify each point
  mutate(cell_id = buffers_pg$cell_id) %>%  
  dplyr:::select(cell_id, elevation_mean = mean, elevation_sd = stdev)

# combine landcover and elevation
env_variables_pg <- inner_join(elev_buffer_pg, lsm_wide_pg, by = "cell_id")

# attach the xy coordinates of the cell centers
env_variables_pg <- buffers_pg %>%  
  st_drop_geometry() %>%  
  dplyr:::select(cell_id, x, y) %>%  
  inner_join(env_variables_pg, by = "cell_id")

# save as csv, dropping any rows with missing variables
write_csv(drop_na(env_variables_pg),
          file.path(DROPBOX_DIR, "data/relative_abundance_environmental_data/environmental-variables_prediction-grid_us-south.csv"), 
          na = "")

#encounter rate definition: probability an eBirder encounters a species on a standard eBird checklist
# load all the packages
library(dplyr)
library(ebirdst)
library(fields)
library(ggplot2)
library(gridExtra)
library(lubridate)
library(mccf1)
library(ranger)
library(readr)
library(scam)
library(sf)
library(terra)
library(tidyr)
library(stringr)
library(GGally)

# set random number seed for reproducibility
set.seed(1)

# environmental variables: landcover and elevation
checklists_env3 <- read_csv(file.path(DOCUMENTS_DIR, "data/relative_abundance_environmental_data/checklists_env2.csv"))

# prediction grid
env_variables_pg <- read_csv(file.path(DOCUMENTS_DIR, "data/relative_abundance_environmental_data/environmental-variables_prediction-grid_us-ga2.csv"))
pred_grid <- env_variables_pg
pred_mat<-matrix(nrow=nrow(checklists_env3), ncol=ncol(checklists_env3[, 22:51]))

# raster template for the grid
r <- rast(file.path(DOCUMENTS_DIR, "data/relative_abundance_environmental_data/prediction-grid_us-south.tif"))
# get the coordinate reference system of the prediction grid
crs <- st_crs(r)

# load gis data for making maps
study_region <- read_sf(file.path(DOCUMENTS_DIR, "data/relative_abundance_environmental_data/gis-data.gpkg"), "ne_states") %>%  
  dplyr::filter(state_code == "US-GA" | state_code == "US-FL" | state_code == "US-AL" | state_code == "US-MS" | state_code == "US-LA" | state_code == "US-TX" | state_code == "US-AR") %>% 
  st_transform(crs = crs) %>%  
  st_geometry()
ne_land <- read_sf(file.path(DOCUMENTS_DIR, "data/relative_abundance_environmental_data/gis-data.gpkg"), "ne_land") %>%  
  st_transform(crs = crs) %>%  
  st_geometry()
ne_country_lines <- read_sf(file.path(DOCUMENTS_DIR, "data/relative_abundance_environmental_data/gis-data.gpkg"), "ne_country_lines") %>%  
  st_transform(crs = crs) %>%  
  st_geometry()
ne_state_lines <- read_sf(file.path(DOCUMENTS_DIR, "data/relative_abundance_environmental_data/gis-data.gpkg"), "ne_state_lines") %>%  
  st_transform(crs = crs) %>%  
  st_geometry()

# remake the test and train types from before
checklists_env3$type <- if_else(runif(nrow(checklists_env3)) <= 0.8, "train", "test")

# sample one checklist per 3km x 3km x 1 week grid for each year
# sample detection/non-detection independently 
checklists_ss <- grid_sample_stratified(checklists_env3,
                                        obs_column = "species_observed",
                                        sample_by = "type", maximum_ss = 122804)

# ADD THE DISTANCE TO COAST
# load the necessary packages
library(dplyr)
library(ebirdst)
library(fields)
library(ggplot2)
library(gridExtra)
library(lubridate)
library(mccf1)
library(ranger)
library(readr)
library(scam)
library(sf)
library(terra)
library(tidyr)
library(stringr)

# make the points into an sf object
checklists_ss_sf <- st_as_sf(checklists_ss, coords = c("longitude", "latitude"), crs = 4326)

# load the raster data from PacIOOS
distance_raster <- terra::rast(file.path(DOCUMENTS_DIR, "data/relative_abundance_environmental_data/dist2coast.txt"))
terra::crs(distance_raster) <- "epsg:4326"

# crop the raster to fit to the checklist points
distance_raster <- terra::crop(distance_raster, checklists_ss_sf)

# extract the value at each point
distance_extracted <- terra::extract(distance_raster, checklists_ss_sf)

# add the values to the original checklist data as variable
checklists_ss$dist2coast <- distance_extracted$dist2coast

# convert checklists to spatial features
all_pts <- checklists_env3 %>%   
  filter(type == "train") %>%  
  st_as_sf(coords = c("longitude","latitude"), crs = 4326) %>% 
  st_transform(crs = crs) %>%  
  dplyr:::select(species_observed)
ss_pts <- checklists_ss %>%   
  filter(type == "train") %>%  
  st_as_sf(coords = c("longitude","latitude"), crs = 4326) %>% 
  st_transform(crs = crs) %>%  
  dplyr:::select(species_observed)
both_pts <- list(before_ss = all_pts, after_ss = ss_pts)

# map
p <- par(mfrow = c(1, 2))
for (i in seq_along(both_pts)) {
  par(mar = c(0.25, 0.25, 0.25, 0.25))
  # set up plot area
  plot(st_geometry(both_pts[[i]]), col = NA)
  # contextual gis data
  plot(ne_land, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)
  #plot(study_region, col = "#cccccc", border = NA, add = TRUE)
  plot(ne_state_lines, col = "#ffffff", lwd = 0.75, add = TRUE)
  plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
  # ebird observations
  # not observed
  plot(st_geometry(both_pts[[i]]),
       pch = 19, cex = 0.1, col = alpha("#555555", 0.25),
       add = TRUE)
  # observed
  plot(filter(both_pts[[i]], species_observed) |> st_geometry(),
       pch = 19, cex = 0.3, col = alpha("#4daf4a", 0.5),
       add = TRUE)
  # legend
  legend("bottomright", bty = "n",
         col = c("#555555", "#4daf4a"),
         legend = c("Non-detection", "Detection"),
         pch = 19)
  box()
  par(new = TRUE, mar = c(0, 0, 3, 0))
  if (names(both_pts)[i] == "before_ss") {
    title("Yellow-Billed Cuckoo eBird Observations\nBefore subsampling")
  } else {
    title("After subsampling")
  }
}
par(p)

# random forest model to relate detection/non-detection to the environmental variables
checklists_train <- checklists_ss %>%  
  filter(type == "train") %>%  
  # select only the columns to be used in the model
  dplyr:::select(species_observed,
                 year, day_of_year, hours_of_day, dist2coast,
                 effort_hours, effort_distance_km, effort_speed_kmph,
                 number_observers, 
                 starts_with("pland_"),
                 starts_with("ed_"),
                 starts_with("elevation_"))

# start the balanced random forest approach
detection_freq <- mean(checklists_train$species_observed)
checklists_train <- na.omit(checklists_train)

# ranger requires a factor response to do classification
er_model <- ranger(formula =  as.factor(species_observed) ~ ., 
                   data = checklists_train,
                   importance = "impurity",
                   probability = TRUE,
                   replace = TRUE, 
                   sample.fraction = c(detection_freq, detection_freq))

# predicted encounter rate based on out of bag samples
er_pred <- er_model$predictions[, 2]
# observed detection, converted back from factor
det_obs <- as.integer(checklists_train$species_observed)
# construct a data frame to train the scam model
obs_pred <- data.frame(obs = det_obs, pred = er_pred)

# train calibration model
calibration_model <- scam(obs ~ s(pred, k = 6, bs = "mpi"), 
                          gamma = 2,
                          data = obs_pred)

# group the predicted encounter rate into bins of width 0.02
# then calculate the mean observed encounter rates in each bin
er_breaks <- seq(0, 1, by = 0.02)
mean_er <- obs_pred %>% 
  mutate(er_bin = cut(pred, breaks = er_breaks, include.lowest = TRUE)) %>% 
  group_by(er_bin) %>% 
  summarise(n_checklists = n(),
            pred = mean(pred), 
            obs = mean(obs),
            .groups = "drop")

# make predictions from the calibration model
calibration_curve <- data.frame(pred = er_breaks)
cal_pred <- predict(calibration_model, calibration_curve, type = "response")
calibration_curve$calibrated <- cal_pred

# compared binned mean encounter rates to calibration model
ggplot(calibration_curve) +
  aes(x = pred, y = calibrated) + theme(axis.text = element_text(size = 14, color = "black", family = "Arial")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_line(color = "blue") +
  geom_point(data = mean_er, 
             aes(x = pred, y = obs),
             size = 2, alpha = 1,
             show.legend = FALSE) +
  labs(x = "Estimated Encounter Rate",
       y = "Observed Encounter Rate") + theme(axis.title = element_text(size = 16, color = "black", family = "Arial")) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1))

# create a threshold encounter rate to make it binary (yes vs no)
# mcc and fscore calculation for various thresholds
mcc_f1 <- mccf1(
  # observed detection/non-detection
  response = obs_pred$obs,
  # predicted encounter rate from random forest
  predictor = obs_pred$pred)

# identify best threshold
mcc_f1_summary <- summary(mcc_f1)
threshold <- mcc_f1_summary$best_threshold[1]

# assessment of the quality of the calibrated random forest model with the train
# get the test set held out from training
checklists_test <- filter(checklists_ss, type == "test") |> 
  mutate(species_observed = as.integer(species_observed))

# make missing values = 0
checklists_test$dist2coast[is.na(checklists_test$dist2coast)] <- 0

# predict to test data using random forest model
pred_er <- predict(er_model, data = checklists_test, type = "response")
# extract probability of detection
pred_er <- pred_er$predictions[, 2]
# convert predictions to binary (presence/absence) using the threshold
pred_binary <- as.integer(pred_er > threshold)
# calibrate
pred_calibrated <- predict(calibration_model, 
                           newdata = data.frame(pred = pred_er), 
                           type = "response") |> 
  as.numeric()
# constrain probabilities to 0-1
pred_calibrated[pred_calibrated < 0] <- 0
pred_calibrated[pred_calibrated > 1] <- 1
# combine observations and estimates
obs_pred_test <- data.frame(id = seq_along(pred_calibrated),
                            # actual detection/non-detection
                            obs = as.integer(checklists_test$species_observed),
                            # binary detection/on-detection prediction
                            pred_binary = pred_binary,
                            # calibrated encounter rate
                            pred_calibrated = pred_calibrated)

# mean squared error (mse)
mse <- mean((obs_pred_test$obs - obs_pred_test$pred_calibrated)^2, na.rm = TRUE)

# precision-recall auc
em <- precrec::evalmod(scores = obs_pred_test$pred_binary, 
                       labels = obs_pred_test$obs)
pr_auc <- precrec::auc(em) |> 
  filter(curvetypes == "PRC") |> 
  pull(aucs)

# calculate metrics for binary prediction: sensitivity, specificity
pa_metrics <- obs_pred_test |> 
  dplyr::select(id, obs, pred_binary) |> 
  PresenceAbsence::presence.absence.accuracy(na.rm = TRUE, st.dev = FALSE)

# mcc and f1
mcc_f1 <- calculate_mcc_f1(obs_pred_test$obs, obs_pred_test$pred_binary)

# combine ppms together
ppms <- data.frame(
  mse = mse,
  sensitivity = pa_metrics$sensitivity,
  specificity = pa_metrics$specificity,
  pr_auc = pr_auc,
  mcc = mcc_f1$mcc,
  f1 = mcc_f1$f1
)
knitr::kable(pivot_longer(ppms, everything()), digits = 3)

# extract predictor importance from the random forest model object
pi <- er_model$variable.importance
pi <- data.frame(predictor = names(pi), importance = pi) |> 
  arrange(desc(importance)) 
head(pi, 15)
# plot importance of top 15 predictors: environmental variables
ggplot(pi[c(3,4,8,9,10,11,12,13,14,15),]) + theme_classic() +
  aes(x = reorder(predictor, importance), y = importance) + theme(axis.text = element_text(size = 20, family = "Arial", colour = "black")) +
  geom_col(fill = "black") +
  geom_hline(yintercept = 0, lwd = 2, colour = "black") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(labels = c("elevation_mean" = "Elevation Mean",
                              "dist2coast" = "Distance To Coast (km)",
                              "pland_c06_closed_shrubland" = "Closed Shrubland (Land Cover)",
                              "elevation_sd" = "Elevation (SD)",
                              "pland_c13_urban" = "Urban (Land Cover)",
                              "ed_c06_closed_shrubland" = "Closed Shrubland (Edge)",
                              "pland_c09_savanna" = "Savanna (Land Cover)",
                              "pland_c10_grassland" = "Grassland (Land Cover)",
                              "pland_c08_woody_savanna" = "Woody Savanna (Land Cover)",
                              "pland_c04_deciduous_broadleaf" = "Deciduous Broadleaf (Land Cover)")) +
  coord_flip() +
  labs(x = NULL, 
       y = "Predictor Importance (Gini Index)") +
  theme(panel.grid = element_blank(),
        panel.grid.major.x = element_line(colour = "black"), axis.title.x=element_text(size=25, family = "Arial"), axis.text.y=element_text(size=25, family = "Arial"))

# plot importance of top 15 predictors: effort variables
ggplot(pi[c(1,2,5,6,7),]) + theme_classic() +
  aes(x = reorder(predictor, importance), y = importance) + theme(axis.text = element_text(size = 20, family = "Arial", colour = "black")) +
  geom_col(fill = "black") +
  geom_hline(yintercept = 0, lwd = 2, colour = "black") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(labels = c("day_of_year" = "Day Of Year",
                              "effort_hours" = "Effort Hours",
                              "effort_speed_kmph" = "Effort Speed(km/h)",
                              "effort_distance_km" = "Effort Distance (km)",
                              "hours_of_day" = "Hours Of Day")) +
  coord_flip() +
  labs(x = NULL, 
       y = "Predictor Importance (Gini Index)") +
  theme(panel.grid = element_blank(),
        panel.grid.major.x = element_line(colour = "black"), axis.title.x=element_text(size=25, family = "Arial"), axis.text.y=element_text(size=25, family = "Arial"))


# function to calculate partial dependence for a single predictor
calculate_pd <- function(predictor, er_model, count_model, cal_model, data, 
                         x_res = 25, n = 1000) {
  # create prediction grid using quantiles
  x_grid <- quantile(data[[predictor]],
                     probs = seq(from = 0, to = 1, length = x_res),
                     na.rm = TRUE)
  # remove duplicates
  x_grid <- x_grid[!duplicated(signif(x_grid, 8))]
  x_grid <- unname(unique(x_grid))
  grid <- data.frame(predictor = predictor, x = x_grid)
  names(grid) <- c("predictor", predictor)
  
  # subsample training data
  n <- min(n, nrow(data))
  data <- data[sample(seq.int(nrow(data)), size = n, replace = FALSE), ]
  
  # drop focal predictor from data
  data <- data[names(data) != predictor]
  grid <- merge(grid, data, all = TRUE)
  
  # predict
  p_er <- predict(er_model, data = grid)
  grid$predicted_er <- p_er$predictions[, 2]
  p_count <- predict(count_model, data = grid)
  
  # summarize
  pd <- grid[, c("predictor", predictor)]
  names(pd) <- c("predictor", "x")
  pd$encounter_rate <- p_er$predictions[, 2]
  pd$count <- p_count$predictions
  
  er_cal <- stats::predict(cal_model, newdata = tibble(pred = pd$encounter_rate), 
                           type = "response")
  er_cal[er_cal < 0] <- 0
  er_cal[er_cal > 1] <- 1
  pd$er_cal <- er_cal
  
  pd$relative_abundance <- pd$er_cal * pd$count
  pd <- dplyr::group_by(pd, predictor, x) %>% 
    dplyr::summarise(relative_abundance = mean(relative_abundance, na.rm = TRUE),
                     .groups = "drop")
  
  return(pd)
}

# top 9 predictors other than date
top_pred <- pi %>% 
  top_n(n = 15, wt = importance) %>% 
  arrange(desc(importance))

# calculate partial dependence for each predictor
# map is used to iteratively apply calculate_pd to each predictor
pd <- NULL
for (predictor in top_pred$predictor) {
  pd <- calculate_pd(predictor, er_model = er_model, count_model = count_model, 
                     cal_model = calibration_model, data = checklists_train) %>% 
    bind_rows(pd, .)
}

library(scales)
# plot
ggplot(pd[76:100,]) +
  aes(x = x, y = relative_abundance) + scale_y_continuous(breaks = seq(0, 0.12, by = 0.03), limits = c(0, 0.12)) + 
  geom_line() +
  geom_point() +
  facet_wrap(~ factor(predictor), ncol = 2, scales = "free") +
  theme_minimal() +
  theme_minimal() + labs(title = "A. Distance To Coast (km)") +
  theme(panel.grid = element_blank(), strip.text = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks  = element_line(color = "black"), text = element_text(size = 45), plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(),
        axis.title.y = element_blank(),axis.text.y = element_text(size = 50, colour = "black", family = "Arial"), axis.text.x = element_text(size = 50, colour = "black", family = "Arial"))

ggplot(pd[160:176,]) +
  aes(x = x, y = relative_abundance) + scale_y_continuous(breaks = seq(0, 0.12, by = 0.03), limits = c(0, 0.12)) + 
  geom_line() +
  geom_point() +
  facet_wrap(~ factor(predictor), ncol = 2, scales = "free") +
  theme_minimal() +
  theme_minimal() + labs(title = "B. Closed Shrubland (Land Cover)") +
  theme(panel.grid = element_blank(), strip.text = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks  = element_line(color = "black"), text = element_text(size = 45), plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(),
        axis.title.y = element_blank(),axis.text.y = element_text(size = 50, colour = "black", family = "Arial"), axis.text.x = element_text(size = 50, colour = "black", family = "Arial"))
  
ggplot(pd[177:200,]) +
  aes(x = x, y = relative_abundance) + scale_y_continuous(breaks = seq(0, 0.12, by = 0.03), limits = c(0, 0.12)) + 
  geom_line() +
  geom_point() +
  facet_wrap(~ factor(predictor), ncol = 2, scales = "free") +
  theme_minimal() +
  theme_minimal() + labs(title = "C. Elevation (SD)") +
  theme(panel.grid = element_blank(), strip.text = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks  = element_line(color = "black"), text = element_text(size = 45), plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(),
        axis.title.y = element_blank(),axis.text.y = element_text(size = 50, colour = "black", family = "Arial"), axis.text.x = element_text(size = 50, colour = "black", family = "Arial"))

ggplot(pd[202:213,]) +
  aes(x = x, y = relative_abundance) + scale_y_continuous(breaks = seq(0, 0.12, by = 0.03), limits = c(0, 0.12)) + 
  geom_line() +
  geom_point() +
  facet_wrap(~ factor(predictor), ncol = 2, scales = "free") +
  theme_minimal() +
  theme_minimal() + labs(title = "D. Urban (Land Cover)") +
  theme(panel.grid = element_blank(), strip.text = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks  = element_line(color = "black"), text = element_text(size = 45), plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(),
        axis.title.y = element_blank(),axis.text.y = element_text(size = 50, colour = "black", family = "Arial"), axis.text.x = element_text(size = 50, colour = "black", family = "Arial"))

ggplot(pd[214:230,]) +
  aes(x = x, y = relative_abundance) + scale_y_continuous(breaks = seq(0, 0.12, by = 0.03), limits = c(0, 0.12)) + 
  geom_line() +
  geom_point() +
  facet_wrap(~ factor(predictor), ncol = 2, scales = "free") +
  theme_minimal() +
  theme_minimal() + labs(title = "E. Closed Shrubland (Edge)") +
  theme(panel.grid = element_blank(), strip.text = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks  = element_line(color = "black"), text = element_text(size = 45), plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(),
        axis.title.y = element_blank(),axis.text.y = element_text(size = 50, colour = "black", family = "Arial"), axis.text.x = element_text(size = 50, colour = "black", family = "Arial"))

ggplot(pd[231:250,]) +
  aes(x = x, y = relative_abundance) + scale_y_continuous(breaks = seq(0, 0.12, by = 0.03), limits = c(0, 0.12)) + 
  geom_line() +
  geom_point() +
  facet_wrap(~ factor(predictor), ncol = 2, scales = "free") +
  theme_minimal() +
  theme_minimal() + labs(title = "F. Savanna (Land Cover)") +
  theme(panel.grid = element_blank(), strip.text = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks  = element_line(color = "black"), text = element_text(size = 45), plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(),
        axis.title.y = element_blank(),axis.text.y = element_text(size = 50, colour = "black", family = "Arial"), axis.text.x = element_text(size = 50, colour = "black", family = "Arial"))

ggplot(pd[251:272,]) +
  aes(x = x, y = relative_abundance) + scale_y_continuous(breaks = seq(0, 0.12, by = 0.03), limits = c(0, 0.12)) + 
  geom_line() +
  geom_point() +
  facet_wrap(~ factor(predictor), ncol = 2, scales = "free") +
  theme_minimal() +
  theme_minimal() + labs(title = "G. Grassland (Land Cover)") +
  theme(panel.grid = element_blank(), strip.text = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks  = element_line(color = "black"), text = element_text(size = 45), plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(),
        axis.title.y = element_blank(),axis.text.y = element_text(size = 50, colour = "black", family = "Arial"), axis.text.x = element_text(size = 50, colour = "black", family = "Arial"))

ggplot(pd[273:292,]) +
  aes(x = x, y = relative_abundance) + scale_y_continuous(breaks = seq(0, 0.12, by = 0.03), limits = c(0, 0.12)) + 
  geom_line() +
  geom_point() +
  facet_wrap(~ factor(predictor), ncol = 2, scales = "free") +
  theme_minimal() +
  theme_minimal() + labs(title = "H. Woody Savanna (Land Cover)") +
  theme(panel.grid = element_blank(), strip.text = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks  = element_line(color = "black"), text = element_text(size = 45), plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(),
        axis.title.y = element_blank(),axis.text.y = element_text(size = 50, colour = "black", family = "Arial"), axis.text.x = element_text(size = 50, colour = "black", family = "Arial"))

ggplot(pd[293:312,]) +
  aes(x = x, y = relative_abundance) + scale_y_continuous(breaks = seq(0, 0.12, by = 0.03), limits = c(0, 0.12)) + 
  geom_line() +
  geom_point() +
  facet_wrap(~ factor(predictor), ncol = 2, scales = "free") +
  theme_minimal() +
  theme_minimal() + labs(title = "I. Deciduous Broadleaf (Land Cover)") +
  theme(panel.grid = element_blank(), strip.text = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks  = element_line(color = "black"), text = element_text(size = 45), plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(),
        axis.title.y = element_blank(),axis.text.y = element_text(size = 50, colour = "black", family = "Arial"), axis.text.x = element_text(size = 50, colour = "black", family = "Arial"))

ggplot(pd[1:25,]) +
  aes(x = x, y = relative_abundance) + scale_y_continuous(breaks = seq(0, 0.12, by = 0.03), limits = c(0, 0.12)) + 
  geom_line() +
  geom_point() +
  facet_wrap(~ factor(predictor), ncol = 2, scales = "free") +
  theme_minimal() +
  theme_minimal() + labs(title = "J. Day Of Year") +
  theme(panel.grid = element_blank(), strip.text = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks  = element_line(color = "black"), text = element_text(size = 45), plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(),
        axis.title.y = element_blank(),axis.text.y = element_text(size = 50, colour = "black", family = "Arial"), axis.text.x = element_text(size = 50, colour = "black", family = "Arial"))

ggplot(pd[26:50,]) +
  aes(x = x, y = relative_abundance) + scale_y_continuous(breaks = seq(0, 0.12, by = 0.03), limits = c(0, 0.12)) + 
  geom_line() +
  geom_point() +
  facet_wrap(~ factor(predictor), ncol = 2, scales = "free") +
  theme_minimal() +
  theme_minimal() + labs(title = "K. Effort Hours") +
  theme(panel.grid = element_blank(), strip.text = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks  = element_line(color = "black"), text = element_text(size = 45), plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(),
        axis.title.y = element_blank(),axis.text.y = element_text(size = 50, colour = "black", family = "Arial"), axis.text.x = element_text(size = 50, colour = "black", family = "Arial"))

ggplot(pd[51:75,]) +
  aes(x = x, y = relative_abundance) + scale_y_continuous(breaks = seq(0, 0.12, by = 0.03), limits = c(0, 0.12)) + 
  geom_line() +
  geom_point() +
  facet_wrap(~ factor(predictor), ncol = 2, scales = "free") +
  theme_minimal() +
  theme_minimal() + labs(title = "L. Elevation Mean") +
  theme(panel.grid = element_blank(), strip.text = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks  = element_line(color = "black"), text = element_text(size = 45), plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(),
        axis.title.y = element_blank(),axis.text.y = element_text(size = 50, colour = "black", family = "Arial"), axis.text.x = element_text(size = 50, colour = "black", family = "Arial"))

ggplot(pd[101:117,]) +
  aes(x = x, y = relative_abundance) + scale_y_continuous(breaks = seq(0, 0.12, by = 0.03), limits = c(0, 0.12)) + 
  geom_line() +
  geom_point() +
  facet_wrap(~ factor(predictor), ncol = 2, scales = "free") +
  theme_minimal() +
  theme_minimal() + labs(title = "M. Effort Speed (km/h)") +
  theme(panel.grid = element_blank(), strip.text = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks  = element_line(color = "black"), text = element_text(size = 45), plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(),
        axis.title.y = element_blank(),axis.text.y = element_text(size = 50, colour = "black", family = "Arial"), axis.text.x = element_text(size = 50, colour = "black", family = "Arial"))

ggplot(pd[118:134,]) +
  aes(x = x, y = relative_abundance) + scale_y_continuous(breaks = seq(0, 0.12, by = 0.03), limits = c(0, 0.12)) + 
  geom_line() +
  geom_point() +
  facet_wrap(~ factor(predictor), ncol = 2, scales = "free") +
  theme_minimal() +
  theme_minimal() + labs(title = "N. Effort Distance (km)") +
  theme(panel.grid = element_blank(), strip.text = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks  = element_line(color = "black"), text = element_text(size = 45), plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(),
        axis.title.y = element_blank(),axis.text.y = element_text(size = 50, colour = "black", family = "Arial"), axis.text.x = element_text(size = 50, colour = "black", family = "Arial"))

ggplot(pd[135:159,]) +
  aes(x = x, y = relative_abundance) + scale_y_continuous(breaks = seq(0, 0.12, by = 0.03), limits = c(0, 0.12)) + scale_x_continuous(breaks = seq(0, 24, by = 6)) + 
  geom_line() +
  geom_point() +
  facet_wrap(~ factor(predictor), ncol = 2) +
  theme_minimal() +
  theme_minimal() + labs(title = "O. Hours Of Day") +
  theme(panel.grid = element_blank(), strip.text = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks  = element_line(color = "black"), text = element_text(size = 45), plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(),
        axis.title.y = element_blank(),axis.text.y = element_text(size = 50, colour = "black", family = "Arial"), axis.text.x = element_text(size = 50, colour = "black", family = "Arial"))

# find peak time of day from partial dependence
pd_time <- calculate_pd("hours_of_day",
                        er_model = er_model, 
                        calibration_model = calibration_model,
                        data = checklists_train) |> 
  dplyr::select(hours_of_day = x, encounter_rate)

# histogram
g_hist <- ggplot(checklists_train) +
  aes(x = hours_of_day) +
  geom_histogram(binwidth = 1, center = 0.5, color = "grey30",
                 fill = "grey50") +
  scale_x_continuous(breaks = seq(0, 24, by = 3)) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Hours since midnight",
       y = "# checklists",
       title = "Distribution of observation start times")

# partial dependence plot
g_pd <- ggplot(pd_time) +
  aes(x = hours_of_day, y = encounter_rate) +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 24, by = 3)) +
  labs(x = "Hours since midnight",
       y = "Encounter rate",
       title = "Observation start time partial dependence")

# combine
grid.arrange(g_hist, g_pd)

# calculate detection frequency for the balance random forest
detection_freq <- mean(checklists_train$species_observed)

checklists_train <- checklists_ss %>%  
  filter(type == "train") %>%  
  # select only the columns to be used in the model
  dplyr:::select(species_observed, observation_count,
                 year, day_of_year, hours_of_day,
                 effort_hours, effort_distance_km, effort_speed_kmph,
                 number_observers, 
                 starts_with("pland_"),
                 starts_with("ed_"),
                 starts_with("elevation_"))

# train a random forest model for encounter rate
train_er <- dplyr::select(checklists_train, -observation_count)
er_model <- ranger(formula =  as.factor(species_observed) ~ ., 
                   data = train_er,
                   importance = "impurity",
                   probability = TRUE,
                   replace = TRUE,
                   sample.fraction = c(detection_freq, detection_freq))

# select the mcc-f1 optimizing occurrence threshold
obs_pred <- data.frame(obs = as.integer(train_er$species_observed), 
                       pred = er_model$predictions[, 2])
mcc_f1 <- mccf1(response = obs_pred$obs, predictor = obs_pred$pred)
mcc_f1_summary <- summary(mcc_f1)
threshold <- mcc_f1_summary$best_threshold[1]

# calibration model
calibration_model <- scam(obs ~ s(pred, k = 6, bs = "mpi"), 
                          gamma = 2,
                          data = obs_pred)

train_count <- checklists_train
train_count$pred_er <- er_model$predictions[, 2]
# subset to only observed or predicted detections
train_count <- train_count |> 
  filter(!is.na(observation_count),
         observation_count > 0 | pred_er > threshold) |> 
  select(-species_observed, -pred_er)

predicted_er <- predict(er_model, data = train_count, type = "response")
predicted_er <- predicted_er$predictions[, 2]
train_count$predicted_er <- predicted_er

count_model <- ranger(formula = observation_count ~ .,
                      data = train_count,
                      importance = "impurity",
                      replace = TRUE)


pred_grid_eff <- pred_grid |> 
  mutate(observation_date = ymd("2023-05-01"),
         year = year(observation_date),
         day_of_year = yday(observation_date),
         # determined as optimal time for detection in previous chapter
         hours_of_day = 6.6,
         effort_distance_km = 2,
         effort_hours = 1,
         effort_speed_kmph = 2,
         number_observers = 1)

# estimate encounter rate
pred_er <- predict(er_model, data = pred_grid_eff, type = "response")
pred_er <- pred_er$predictions[, 2]
# define range-boundary
pred_binary <- as.integer(pred_er > threshold)
# apply calibration
pred_er_cal <- predict(calibration_model, 
                       data.frame(pred = pred_er), 
                       type = "response") |> 
  as.numeric()
# constrain to 0-1
pred_er_cal[pred_er_cal < 0] <- 0
pred_er_cal[pred_er_cal > 1] <- 1

# add predicted encounter rate required for count estimates
pred_grid_eff$predicted_er <- pred_er
# estimate count
pred_count <- predict(count_model, data = pred_grid_eff, type = "response")
pred_count <- pred_count$predictions

# combine predictions with coordinates from prediction grid
predictions <- data.frame(cell_id = pred_grid_eff$cell_id,
                          x = pred_grid_eff$x,
                          y = pred_grid_eff$y,
                          in_range = pred_binary, 
                          encounter_rate = pred_er_cal,
                          count = pred_count)

# add relative abundance estimate
predictions$abundance <- predictions$encounter_rate * predictions$count

# rasterize
layers <- c("in_range", "encounter_rate", "count", "abundance")
r_pred <- predictions |> 
  # convert to spatial features
  st_as_sf(coords = c("x", "y"), crs = crs) |> 
  dplyr::select(all_of(layers)) |> 
  # rasterize
  terra::rasterize(r, field = layers)

par(mar = c(0.25, 0.25, 1.25, 0.25))
# set up plot area
plot(study_region,
     main = "Yellow-Billed Cuckoo Range",
     col = NA, border = NA)
plot(ne_land, col = "#cfcfcf", border = "#888888", lwd = 0.5, add = TRUE)

# convert binary prediction to categorical
r_range <- as.factor(r_pred[["in_range"]])
plot(r_range, col = c("#e6e6e6", "forestgreen"),
     maxpixels = ncell(r_range),
     legend = FALSE, axes = FALSE, bty = "n",
     add = TRUE)

# borders
plot(ne_state_lines, col = "#ffffff", lwd = 0.75, add = TRUE)
plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
plot(study_region, border = "#000000", col = NA, lwd = 1, add = TRUE)
box()

dev.off()

par(mar = c(4, 0.25, 0.25, 0.25))
# set up plot area
plot(study_region, col = NA, border = NA)
plot(ne_land, col = "#cfcfcf", border = "#888888", lwd = 0.5, add = TRUE)

# define quantile breaks
brks <- global(r_pred[["encounter_rate"]], fun = quantile, 
               probs = seq(0, 1, 0.1), na.rm = TRUE) |> 
  as.numeric() |> 
  unique()
# label the bottom, middle, and top value
lbls <- round(c(0, median(brks), max(brks)), 2)
# ebird status and trends color palette
pal <- ebirdst_palettes(length(brks) - 1)
plot(r_pred[["encounter_rate"]], 
     col = pal, breaks = brks, 
     maxpixels = ncell(r_pred),
     legend = FALSE, axes = FALSE, bty = "n",
     add = TRUE)

# borders
plot(ne_state_lines, col = "#ffffff", lwd = 0.75, add = TRUE)
plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
plot(study_region, border = "#000000", col = NA, lwd = 1, add = TRUE)
box()

# legend
par(new = TRUE, mar = c(0, 0, 0, 0))
title <- "Yellow-Billed Cuckoo Encounter Rate"
image.plot(zlim = c(0, 1), legend.only = TRUE, 
           col = pal, breaks = seq(0, 1, length.out = length(brks)),
           smallplot = c(0.25, 0.75, 0.03, 0.06),
           horizontal = TRUE,
           axis.args = list(at = c(0, 0.5, 1), labels = lbls,
                            fg = "black", col.axis = "black",
                            cex.axis = 0.75, lwd.ticks = 0.5,
                            padj = -1.5),
           legend.args = list(text = title,
                              side = 3, col = "black",
                              cex = 1, line = 0))

# download all the packages
library(dplyr)
library(ebirdst)
library(fields)
library(ggplot2)
library(gridExtra)
library(lubridate)
library(mccf1)
library(ranger)
library(readr)
library(scam)
library(sf)
library(terra)
library(tidyr)
library(geosphere)

# DATA PREPARATION
# set random number seed for reproducibility
set.seed(1)

# environmental variables: landcover and elevation
checklists <- read_csv(file.path(DOCUMENTS_DIR, "data/relative_abundance_environmental_data/checklists_env2.csv"))

# prediction grid
env_variables_pg <- read_csv(file.path(DOCUMENTS_DIR, "data/relative_abundance_environmental_data/environmental-variables_prediction-grid_us-ga2.csv"))
pred_grid <- env_variables_pg
pred_mat<-matrix(nrow=nrow(checklists), ncol=ncol(checklists[, 22:51]))

# raster template for the grid
r <- rast(file.path(DOCUMENTS_DIR, "data/relative_abundance_environmental_data/prediction-grid_us-south.tif"))
# get the coordinate reference system of the prediction grid
crs <- st_crs(r)

# load gis data for making maps
study_region <- read_sf(file.path(DOCUMENTS_DIR, "data/relative_abundance_environmental_data/gis-data.gpkg"), "ne_states") %>%  
  filter(state_code == "US-GA" | state_code == "US-FL" | state_code == "US-AL" | state_code == "US-MS" | state_code == "US-LA" | state_code == "US-TX" | state_code == "US-AR") %>% 
  st_transform(crs = crs) %>%  
  st_geometry()
ne_land <- read_sf(file.path(DOCUMENTS_DIR, "data/relative_abundance_environmental_data/gis-data.gpkg"), "ne_land") %>%  
  st_transform(crs = crs) %>%  
  st_geometry()
ne_country_lines <- read_sf(file.path(DOCUMENTS_DIR, "data/relative_abundance_environmental_data/gis-data.gpkg"), "ne_country_lines") %>%  
  st_transform(crs = crs) %>%  
  st_geometry()
ne_state_lines <- read_sf(file.path(DOCUMENTS_DIR, "data/relative_abundance_environmental_data/gis-data.gpkg"), "ne_state_lines") %>%  
  st_transform(crs = crs) %>%  
  st_geometry()

# remake the test and train types from before
checklists$type <- if_else(runif(nrow(checklists)) <= 0.8, "train", "test")

# sample one checklist per 3km x 3km x 1 week grid for each year
# sample detection/non-detection independently 
checklists_ss <- grid_sample_stratified(checklists,
                                        obs_column = "species_observed",
                                        sample_by = "type", maximum_ss = 122804)

# ADD THE DISTANCE TO COAST
# load the necessary packages
library(dplyr)
library(ebirdst)
library(fields)
library(ggplot2)
library(gridExtra)
library(lubridate)
library(mccf1)
library(ranger)
library(readr)
library(scam)
library(sf)
library(terra)
library(tidyr)
library(stringr)

# make the points into an sf object
checklists_ss_sf <- st_as_sf(checklists_ss, coords = c("longitude", "latitude"), crs = 4326)

# load the raster data from PacIOOS
distance_raster <- terra::rast(file.path(DOCUMENTS_DIR, "data/relative_abundance_environmental_data/dist2coast.txt"))
terra::crs(distance_raster) <- "epsg:4326"

# crop the raster to fit to the checklist points
distance_raster <- terra::crop(distance_raster, checklists_ss_sf)

# extract the value at each point
distance_extracted <- terra::extract(distance_raster, checklists_ss_sf)

# add the values to the original checklist data as variable
checklists_ss$dist2coast <- distance_extracted$dist2coast

checklists_train <- checklists_ss |> 
  filter(type == "train") |> 
  # select only the columns to be used in the model
  dplyr::select(species_observed, observation_count,
                year, day_of_year, hours_of_day,
                effort_hours, effort_distance_km, effort_speed_kmph,
                number_observers, 
                starts_with("pland_"),
                starts_with("ed_"),
                starts_with("elevation_"))

# ENCOUNTER RATE
# calculate detection frequency for the balance random forest
detection_freq <- mean(checklists_train$species_observed)

# train a random forest model for encounter rate
train_er <- dplyr::select(checklists_train, -observation_count)
er_model <- ranger(formula =  as.factor(species_observed) ~ ., 
                   data = train_er,
                   importance = "impurity",
                   probability = TRUE,
                   replace = TRUE,
                   sample.fraction = c(detection_freq, detection_freq))

# select the mcc-f1 optimizing occurrence threshold
obs_pred <- data.frame(obs = as.integer(train_er$species_observed), 
                       pred = er_model$predictions[, 2])
mcc_f1 <- mccf1(response = obs_pred$obs, predictor = obs_pred$pred)
mcc_f1_summary <- summary(mcc_f1)
threshold <- mcc_f1_summary$best_threshold[1]

# calibration model
calibration_model <- scam(obs ~ s(pred, k = 6, bs = "mpi"), 
                          gamma = 2,
                          data = obs_pred)

# COUNT
# attach the predicted encounter rate based on out of bag samples
train_count <- checklists_train
train_count$pred_er <- er_model$predictions[, 2]
# subset to only observed or predicted detections
train_count <- train_count |> 
  filter(!is.na(observation_count),
         observation_count > 0 | pred_er > threshold) |> 
  select(-species_observed, -pred_er)

# include estimated encounter rate as a predictor
predicted_er <- predict(er_model, data = train_count, type = "response")
predicted_er <- predicted_er$predictions[, 2]
train_count$predicted_er <- predicted_er

count_model <- ranger(formula = observation_count ~ .,
                      data = train_count,
                      importance = "impurity",
                      replace = TRUE)

# ASSESSMENT
# get the test set held out from training
checklists_test <- filter(checklists_ss, type == "test") |> 
  mutate(species_observed = as.integer(species_observed)) |> 
  filter(!is.na(observation_count))

# estimate encounter rate for test data
pred_er <- predict(er_model, data = checklists_test, type = "response")
# extract probability of detection
pred_er <- pred_er$predictions[, 2]
# convert to binary using the threshold
pred_binary <- as.integer(pred_er > threshold)
# calibrate
pred_calibrated <- predict(calibration_model, 
                           newdata = data.frame(pred = pred_er), 
                           type = "response") |> 
  as.numeric()
# constrain probabilities to 0-1
pred_calibrated[pred_calibrated < 0] <- 0
pred_calibrated[pred_calibrated > 1] <- 1

# add predicted encounter rate required for count estimates
checklists_test$predicted_er <- pred_er
# estimate count
pred_count <- predict(count_model, data = checklists_test, type = "response")
pred_count <- pred_count$predictions

# relative abundance is the product of encounter rate and count
pred_abundance <- pred_calibrated * pred_count

# combine observations and estimates
obs_pred_test <- data.frame(
  id = seq_along(pred_abundance),
  # actual detection/non-detection
  obs_detected = as.integer(checklists_test$species_observed),
  obs_count = checklists_test$observation_count,
  # model estimates
  pred_binary = pred_binary,
  pred_er = pred_calibrated,
  pred_count = pred_count,
  pred_abundance = pred_abundance
)

# subset to only those checklists where detect occurred
detections_test <- filter(obs_pred_test, obs_detected > 0)

# count metrics, based only on checklists where detect occurred
count_spearman <- cor(detections_test$pred_count, 
                      detections_test$obs_count,
                      method = "spearman")
log_count_pearson <- cor(log(detections_test$pred_count + 1),
                         log(detections_test$obs_count + 1),
                         method = "pearson")

# abundance metrics
abundance_spearman <- cor(detections_test$pred_abundance, 
                          detections_test$obs_count,
                          method = "spearman")
log_abundance_pearson <- cor(log(detections_test$pred_abundance + 1),
                             log(detections_test$obs_count + 1),
                             method = "pearson")

# combine ppms together
ppms <- data.frame(
  count_spearman = count_spearman,
  log_count_pearson = log_count_pearson,
  abundance_spearman = abundance_spearman,
  log_abundance_pearson = log_abundance_pearson
)
knitr::kable(pivot_longer(ppms, everything()), digits = 3)

# PREDICTION 
pred_grid_eff <- pred_grid |> 
  mutate(observation_date = ymd("2023-05-01"),
         year = year(observation_date),
         day_of_year = yday(observation_date),
         # determined as optimal time for detection in previous chapter
         hours_of_day = 6.6,
         effort_distance_km = 2,
         effort_hours = 1,
         effort_speed_kmph = 2,
         number_observers = 1)

# estimate encounter rate
pred_er <- predict(er_model, data = pred_grid_eff, type = "response")
pred_er <- pred_er$predictions[, 2]
# define range-boundary
pred_binary <- as.integer(pred_er > threshold)
# apply calibration
pred_er_cal <- predict(calibration_model, 
                       data.frame(pred = pred_er), 
                       type = "response") |> 
  as.numeric()
# constrain to 0-1
pred_er_cal[pred_er_cal < 0] <- 0
pred_er_cal[pred_er_cal > 1] <- 1

# add predicted encounter rate required for count estimates
pred_grid_eff$predicted_er <- pred_er
# estimate count
pred_count <- predict(count_model, data = pred_grid_eff, type = "response")
pred_count <- pred_count$predictions

# combine predictions with coordinates from prediction grid
predictions <- data.frame(cell_id = pred_grid_eff$cell_id,
                          x = pred_grid_eff$x,
                          y = pred_grid_eff$y,
                          in_range = pred_binary, 
                          encounter_rate = pred_er_cal,
                          count = pred_count)

# add relative abundance estimate
predictions$abundance <- predictions$encounter_rate * predictions$count

# rasterize
layers <- c("in_range", "encounter_rate", "count", "abundance")
r_pred <- predictions |> 
  # convert to spatial features
  st_as_sf(coords = c("x", "y"), crs = crs) |> 
  select(all_of(layers)) |> 
  # rasterize
  rasterize(r, field = layers)
print(r_pred)

# project background layers to match the raster's CRS
study_region_proj <- project(vect(study_region), crs(r_pred))
ne_land_proj      <- project(vect(ne_land), crs(r_pred))

# plot the map
par(mar = c(4, 0.25, 0.25, 0.25))
# set up plot area
plot(study_region_proj, axes = FALSE, bty = "n", col = NA, border = NA)
plot(ne_land_proj, col = "#cfcfcf", border = "#888888", lwd = 0.5, add = TRUE)

# define quantile breaks, excluding zeros
brks <- ifel(r_pred[["abundance"]] > 0, r_pred[["abundance"]], NA) |> 
  global(fun = quantile, 
         probs = seq(0, 1, 0.1), na.rm = TRUE) |> 
  as.numeric() |> 
  unique()
# label the bottom, middle, and top value
lbls <- round(c(min(brks), median(brks), max(brks)), 2)
# ebird status and trends color palette
pal <- rev(ebirdst_palettes(length(brks) - 1))
plot(r_pred[["abundance"]], 
     col = c("#e6e6e6", pal), breaks = c(0, brks), 
     maxpixels = ncell(r_pred),
     legend = FALSE, axes = FALSE, bty = "n",
     add = TRUE)

# borders
plot(ne_state_lines, col = "#ffffff", lwd = 0.75, add = TRUE)
plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
plot(study_region, border = "#000000", col = NA, lwd = 1, add = TRUE)

# legend
par(new = TRUE, mar = c(0, 0, 0, 0))
#title <- "Yellow-Billed Cuckoo Relative Abundance (May 2023)"
image.plot(zlim = c(0, 1), legend.only = TRUE, 
           col = pal, breaks = seq(0, 1, length.out = length(brks)),
           smallplot = c(0.25, 0.75, 0.3, 0.33),
           horizontal = TRUE,
           axis.args = list(at = c(0, 0.5, 1), labels = lbls,
                            fg = "black", col.axis = "black",
                            cex.axis = 0.75, lwd.ticks = 0.5,
                            padj = -1.5))


# NEED CODE FOR SATELLITE TRACKING HABITAT SELECTION ***************************


# MORAN I FOR STATUS AND TRENDS ************************************************
# load all the packages
library(spdep)
library(sf)
library(ebirdst)
library(dplyr)
library(terra)
library(raster)
library(stars)
library(tmap)
library(rnaturalearth)

#download YBCU data
set_ebirdst_access_key(EBIRD_KEY)
path <- ebirdst_download_status(species = "yellow-billed cuckoo", path = file.path(DOCUMENTS_DIR, "data/yebcuc2022"))

# load relative abundance raster stack with 52 layers, one for each week
abd <- load_raster(species = "yebcuc", path = file.path(DOCUMENTS_DIR, "data/yebcuc2022"), product = "abundance", period = "week", resolution = "3km")

# ZOOMING IN/MAKING THE EXTENT
#make a boundary polygon for the border states we're studying
# project to same coordinate reference system as the raster data
us_border_states <- ne_states(iso_a2 = "US", returnclass = "sf") %>% dplyr::filter(
  region == "South") %>% dplyr::filter(postal != "SC") %>% dplyr::filter(postal != "NC") %>% dplyr::filter(
    postal != "VA") %>% dplyr::filter(postal != "DC") %>% dplyr::filter(postal != "MD") %>% dplyr::filter(
      postal != "DE") %>% dplyr::filter(postal != "OK") %>% dplyr::filter(postal != "KY") %>% dplyr::filter(
        postal != "TN") %>% dplyr::filter(postal != "WV") %>% st_transform(st_crs(abd))

# crop data to the border states we're studying
abd_us_border_states <- crop(abd, us_border_states)

# create a sf object
abd_stars <- st_as_stars(abd_us_border_states)
abd_sf <- st_as_sf(abd_stars)

# obtain the neighbors of each area
nb <- poly2nb(abd_sf, queen = TRUE)

# create a list of the neighbors
nbw <- nb2listw(nb, style = "W", zero.policy = TRUE)

# set NA values to 0
abd_sf <- abd_sf %>% replace(is.na(.), values = 0)

# Global Moran's I
gmoran <- moran.test(abd_sf$`2022-07-05`, nbw, alternative = "greater")
gmoran

# Moran's I
gmoran[["estimates"]][["Moran I statistic"]]

# z-score
gmoran[["statistic"]]

# p-value
gmoran[["p.value"]]


# MORAN I FOR EBD **************************************************************
# Global Moran's I for the EBD/SED dataset
# load packages
library(sf)
library(dplyr)
library(spdep)

# Use the predicted relative abundance for month of May using EBD/SED
# Define the custom eBird/SED projection (Albers)
ebird_projection <- "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

#get predictions from Encounter Rate script
predictions_sf <- predictions %>%
  # Use the meter-based projection directly
  st_as_sf(coords = c("x", "y"), crs = ebird_projection) %>%
  
  # Buffer by 1500m
  st_buffer(dist = 1500, nQuadSegs = 1, endCapStyle = "SQUARE")

# define neighboring polygons
nb1 <- poly2nb(predictions_sf, queen = TRUE)

# assign weights to neighbors
lw1 <- nb2listw(nb1, style = "W", zero.policy = TRUE)

# compute the statistic
I <- moran.test(predictions_sf$abundance, lw1, alternative = "greater")
I


# KERNEL DENSITY PLOTS FOR SATELLITE TRACKING **********************************
# download the packages we need
library(sf)
library(eks)
library(colorspace)
library(ebirdst)
library(scales)

# download the data
rawdata <- st_read(file.path(DOCUMENTS_DIR, "data_raw/retrackingdatapart2/centroidsGulf9Jul25.shp"))
crs <- st_crs(rawdata)

# add breeding regions 
rawdata <- rawdata %>% mutate(region = case_when(
  id == 198831 ~ "East",
  id == 198847 ~ "Central",
  id == 198852 ~ "East",
  id == 200491 ~ "Central",
  id == 237990 ~ "Central",
  id == 237992 ~ "Central",
  id == 237994 ~ "Central",
  id == 237998 ~ "Central",
  id == 238000 ~ "Central",
  id == 238001 ~ "Central",
  id == 238002 ~ "Central",
  id == 238003 ~ "East",
  id == 238075 ~ "Central",
  id == 238078 ~ "East",
  id == 238079 ~ "Central",
  id == 238080 ~ "East",
  id == 238081 ~ "Central",
  id == 238084 ~ "East",
  id == 238087 ~ "Central",
  id == 238088 ~ "Central",
  id == 254071 ~ "Central",
  id == 254073 ~ "East",
  id == 238076 ~ "Did Not Complete",
  id == 254085 ~ "Did Not Complete",
  id == "179003-2" ~ "Central",
  TRUE ~ "Other"
))

# only one cluster per stopover
rawdata <- rawdata %>% dplyr::filter(!(clusterId == "237992-43" | clusterId == "237992-42" | clusterId == "238001-18" | clusterId == "238078-33" | clusterId == "238078-34" | clusterId == "238088-21" | clusterId == "238088-22" | clusterId == "238088-24" | clusterId == "237990-32" | clusterId == "237990-36" | clusterId == "237992-210" | clusterId == "238001-168" | clusterId == "238079-139" | clusterId == "238079-142" | clusterId == "238079-144" | clusterId == "238080-133" | clusterId == "238080-134" | clusterId == "238081-146" | clusterId == "238087-29" | clusterId == "254071-35" | clusterId == "254071-37" | clusterId == "254073-29"))

# load gis data for making maps
study_region <- read_sf(file.path(DOCUMENTS_DIR, "data/relative_abundance_environmental_data/gis-data.gpkg"), "ne_states") %>%  
  filter(state_code == "US-GA" | state_code == "US-FL" | state_code == "US-AL" | state_code == "US-MS" | state_code == "US-LA" | state_code == "US-TX" | state_code == "US-AR") %>% 
  st_transform(crs = crs) %>% 
  st_geometry()
ne_land <- read_sf(file.path(DOCUMENTS_DIR, "data/relative_abundance_environmental_data/gis-data.gpkg"), "ne_land") %>%  
  st_transform(crs = crs) %>%  
  st_geometry()
ne_country_lines <- read_sf(file.path(DOCUMENTS_DIR, "data/relative_abundance_environmental_data/gis-data.gpkg"), "ne_country_lines") %>%  
  st_transform(crs = crs) %>%  
  st_geometry()
ne_state_lines <- read_sf(file.path(DOCUMENTS_DIR, "data/relative_abundance_environmental_data/gis-data.gpkg"), "ne_state_lines") %>%  
  st_transform(crs = crs) %>%  
  st_geometry()

# calculate KDE
skde1 <- st_kde(rawdata)

pal <- ebirdst_palettes(10)

my_colors <- c("Central" = "black", "East" = "steelblue3", "Did Not Complete" = "lavender")
par(family = "Arial")
# set up plot area
plot(study_region, col = NA, border = NA, reset = FALSE)
plot(ne_land, col = "#cfcfcf", border = "#888888", lwd = 0.5, add = TRUE)
# borders
plot(ne_state_lines, col = "#ffffff", lwd = 0.75, add = TRUE)
plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
# contours
skde_sf <- st_get_contour(skde1, cont = c(50, 75, 95))
plot(st_geometry(skde_sf), 
     col = adjustcolor(c("#E2E4BD", "#E8982C", "#D23e69"), alpha.f = 0.5),
     border = "white",                           
     lwd = 1, 
     add = TRUE)
# data
plot(rawdata, col = my_colors[rawdata$region], pch = 18, cex = 1, add = TRUE)
# add legend
legend("bottomleft",
       legend = c("95%", "75%", "50%"),
       col = c("#E2E4BD", "#E8982C", "#D23e69"),
       pch = 15,
       bty = "n",
       title = "Density")
legend("bottom",
       legend = names(my_colors),
       col = my_colors,
       pch = 18,
       bty = "n",
       y.intersp = 0.7,
       cex = 1,
       title = "Breeding Regions")


# QUADRAT TEST FOR SATELLITE TRACKING ******************************************
# download the packages we need
library(spatstat)
library(dplyr)
library(terra)
library(sf)
library(rgdal)

# download the data
rawdata <- st_read(file.path(DOCUMENTS_DIR, "data_raw/retrackingdatapart2/centroidsGulf9Jul25.shp"))

# only one cluster per stopover
rawdata <- rawdata %>% dplyr::filter(!(clusterId == "237992-43" | clusterId == "237992-42" | clusterId == "238001-18" | clusterId == "238078-33" | clusterId == "238078-34" | clusterId == "238088-21" | clusterId == "238088-22" | clusterId == "238088-24" | clusterId == "237990-32" | clusterId == "237990-36" | clusterId == "237992-210" | clusterId == "238001-168" | clusterId == "238079-139" | clusterId == "238079-142" | clusterId == "238079-144" | clusterId == "238080-133" | clusterId == "238080-134" | clusterId == "238081-146" | clusterId == "238087-29" | clusterId == "254071-35" | clusterId == "254071-37" | clusterId == "254073-29"))

# make the study region
crs <- crs(rawdata)
study_region <- read_sf(file.path(DOCUMENTS_DIR, "data/relative_abundance_environmental_data/gis-data.gpkg"), "ne_states") %>%  
  filter(state_code == "US-GA" | state_code == "US-FL" | state_code == "US-AL" | state_code == "US-MS" | state_code == "US-LA" | state_code == "US-TX" | state_code == "US-AR") %>% 
  st_transform(crs = crs) %>%  
  st_geometry()

# convert to the correct format
study_region <- as.owin(study_region)


# convert to a ppp object
rawdata <- as.ppp(rawdata, study_region)

# divide the window
Q <- quadratcount(rawdata, nx = 8, ny = 7)
plot(Q)

# conduct the quadrat test
quadrat.test(Q, alternative = "clustered", method = "MonteCarlo")