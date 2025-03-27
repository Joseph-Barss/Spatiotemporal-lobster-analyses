
############## Applied analysis ###################
# Libraries used
library(dplyr)
library(readr)
library(tidyr)
library(forcats)
library(lubridate)
library(tibble)
library(ggplot2)
library(sf)
library(terra)
library(rnaturalearth)
library(rnaturalearthhires)
library(sdmTMB)
library(sdmTMBextra)
library(fmesher)
library(future)
library(future.apply)
library(blockCV)
library(visreg)

# Source functions
source("functions.R")

################ Data processing
# Make shapefile for study area
map_data <- rnaturalearth::ne_countries(scale = "large", returnclass = "sf", country = "canada")
box <- c(xmin = -68, ymin = 42, xmax = -63, ymax = 46)
map_fundy <- st_crop(map_data, box)
# Transform the CRS (needed for sdmTMB)
crs_WG84 <- st_crs(map_data)
crs_UTM20N <- 32620
map_data <- st_transform(map_data, crs_UTM20N)
map_fundy <- st_transform(map_fundy, crs_UTM20N)
shape <- make_shape()

# Load the depth data
bathymetry <- rast("ofi_DEM_UTMZ20.tif")

# Read in and process the trawl survey data
lobster_data <- read_count_data(shape, bathymetry)

#Make prediction grids, removing depths under 8m
grid <- make_prediction_grid(shape = shape, bathymetry = bathymetry, 
                             sweptarea = 0.01, depth_cutoff = 8)
grid_yrs <- replicate_df(grid, "year", unique(lobster_data$year))
grid_yrs$year_fac <- as_factor(grid_yrs$year)

# Make nonconvex hull mesh
barrier_mesh <- make_barrier_mesh(lobster_data, map_data)

############### Descriptive plots

# Map of depth
ggplot(map_fundy)+
  geom_sf()+
  geom_tile(data = grid, aes(xmetres, ymetres, fill = -depth))+
  scale_fill_continuous()+
  theme_light()+
  labs(fill = "Depth (m)", x = "Longitude", y = "Latitude")

# Number of tows per year, by gear type
ggplot(lobster_data, aes(year, fill = fct_rev(gear)))+
  geom_histogram(position = "identity", alpha = 0.4, bins = 29) +
  labs(x = "Year", y = "Tows", fill = "Gear type") +
  scale_fill_viridis_d(option = "D", direction = -1)+
  theme_light()

# Number of tows per month, by gear type
ggplot(lobster_data, aes(fct_relevel(month, "5"), fill = fct_rev(gear)))+
  geom_bar(alpha = 0.4) +
  scale_x_discrete(labels = c("May", "June", "July", "August", "September", "October"))+
  labs(x = "Month", y = "Tows", fill = "Gear type") +
  scale_fill_viridis_d(option = "D", direction = -1)+
  theme_light()

# Map of tows by gear type
ggplot(map_fundy) + geom_sf() +
  geom_point(data = lobster_data, aes(x = X * 1000, y = Y * 1000, col = fct_rev(gear)), alpha = 0.4, size = 0.01) +
  scale_colour_viridis_d(option = "D", direction = -1)+
  theme_light() +
  guides(colour = guide_legend(override.aes = list(size=2)))+
  labs(x = "Longitude", y = "Latitude", col = "Gear")

# Histogram of counts, by gear type
ggplot(lobster_data, aes(count, fill = fct_rev(gear)))+
  geom_histogram(position = "identity", alpha = 0.4, bins = 50) +
  scale_fill_viridis_d(option = "D", direction = -1)+
  scale_y_continuous(transform = "log1p", breaks = 10^(0:4))+
  labs(x = "Count of lobsters", y = "Frequency (log scale)", fill = "Gear type") +
  theme_light()

# Histogram of swept areas, by gear type
ggplot(lobster_data, aes(area_km2, fill = fct_rev(gear)))+
  geom_histogram(position = "identity", alpha = 0.4, bins = 50) +
  scale_fill_viridis_d(option = "D", direction = -1)+
  scale_y_continuous(transform = "log1p", breaks = 10^(0:4))+
  labs(x = "Swept area (km2)", y = "Frequency (log scale)", fill = "Gear type") +
  theme_light()

# Map of counts
ggplot(map_fundy) + geom_sf() +
  geom_point(data = lobster_data, aes(x = X * 1000, y = Y * 1000, col = count), alpha = 0.4, size = 0.01) +
  scale_colour_viridis_c(option = "B", direction = -1, transform = "log1p", breaks = 10^(0:3))+
  theme_light()+
  labs(x = "Longitude", y = "Latitude", col = "Lobster count")

# Plot of mesh
plot(barrier_mesh)

############### Applied Analysis

# Model with no random fields
norfmod <- sdmTMB(count ~ 0+year_fac+logdepth+gear+moult, data = lobster_data,
                  mesh = barrier_mesh, family = nbinom2(), offset = "log_area",
                  time = "year", spatial = "off", spatiotemporal = "off")
# Model with spatial random field
spatmod <- sdmTMB(count ~ 0+year_fac+logdepth+gear+moult, data = lobster_data,
                  mesh = barrier_mesh, family = nbinom2(), offset = "log_area",
                  time = "year", spatial = "on", spatiotemporal = "off")
# Model with IID spatiotemporal random fields
iidmod <- sdmTMB(count ~ 0+year_fac+logdepth+gear+moult, data = lobster_data,
                 mesh = barrier_mesh, family = nbinom2(), offset = "log_area",
                 time = "year", spatial = "on", spatiotemporal = "iid")
# Model with random walk spatiotemporal random fields
rwmod <- sdmTMB(count ~ 0+year_fac+logdepth+gear+moult, data = lobster_data,
                mesh = barrier_mesh, family = nbinom2(), offset = "log_area",
                time = "year", spatial = "on", spatiotemporal = "rw")
# Model with AR(1) spatiotemporal random fields
ar1mod <- sdmTMB(count ~ 0+year_fac+logdepth+gear+moult, data = lobster_data,
                 mesh = barrier_mesh, family = nbinom2(), offset = "log_area",
                 time = "year", spatial = "on", spatiotemporal = "ar1")

# Check marginal AICs
AIC(norfmod, spatmod, iidmod, rwmod, ar1mod)
# Check BICs
BIC(norfmod, spatmod, iidmod, rwmod, ar1mod)
# Check conditional AICs
cAIC(spatmod)
cAIC(iidmod) 
cAIC(rwmod)
cAIC(ar1mod)

# Cross-validation with blockCV
# Create spatial blocks
lobster_dat_sf <- st_as_sf(lobster_data, coords = c("xmetres", "ymetres"), crs = crs_UTM20N)
set.seed(940)
spatial_blocks <- cv_spatial(lobster_dat_sf, k = 5, size = 50000)
cv_plot(spatial_blocks, lobster_dat_sf)
# Run cross-validation in parallel
plan(multisession, workers = 5)
# No random fields
norf_cv <- sdmTMB_cv(count ~ 0+year_fac+logdepth+gear+moult, data = lobster_data,
                     mesh = barrier_mesh, family = nbinom2(), offset = "log_area",
                     time = "year", spatial = "off", spatiotemporal = "off", parallel = TRUE, 
                     fold_ids = spatial_blocks$folds_ids)
# Spatial random field
spat_cv <- sdmTMB_cv(count ~ 0+year_fac+logdepth+gear+moult, data = lobster_data,
                     mesh = barrier_mesh, family = nbinom2(), offset = "log_area",
                     time = "year", spatial = "on", spatiotemporal = "off", parallel = TRUE, 
                     fold_ids = spatial_blocks$folds_ids)
# IID random fields
iid_cv <- sdmTMB_cv(count ~ 0+year_fac+logdepth+gear+moult, data = lobster_data,
                    mesh = barrier_mesh, family = nbinom2(), offset = "log_area",
                    time = "year", spatial = "on", spatiotemporal = "iid", parallel = TRUE, 
                    fold_ids = spatial_blocks$folds_ids)
# Random walk random fields
rw_cv <- sdmTMB_cv(count ~ 0+year_fac+logdepth+gear+moult, data = lobster_data,
                   mesh = barrier_mesh, family = nbinom2(), offset = "log_area",
                   time = "year", spatial = "on", spatiotemporal = "rw", parallel = TRUE, 
                   fold_ids = spatial_blocks$folds_ids)
# AR(1) random fields
ar1_cv <- sdmTMB_cv(count ~ 0+year_fac+logdepth+gear+moult, data = lobster_data,
                    mesh = barrier_mesh, family = nbinom2(), offset = "log_area",
                    time = "year", spatial = "on", spatiotemporal = "ar1", parallel = TRUE, 
                    fold_ids = spatial_blocks$folds_ids)
plan(sequential)
# Check sum of log-likelihoods
norf_cv$sum_loglik
spat_cv$sum_loglik
iid_cv$sum_loglik
rw_cv$sum_loglik
ar1_cv$sum_loglik

# We choose the model with a spatial random field
summary(spatmod)
# Plot the effect of log(depth)
logdepthplot <- visreg(spatmod, xvar = "logdepth", data = lobster_data, gg = TRUE)
logdepthplot+xlab("log(Depth)")+ylab("Link scale")+theme_light()

# Residual analysis
set.seed(11)
lobster_data$spresids <- residuals(spatmod)

plot_resids_spatial(lobster_data, plot_type = "histogram")
plot_resids_spatial(lobster_data, plot_type = "qq-plot")
plot_resids_spatial(lobster_data, plot_type = "map")

# Make predictions on a grid
spatialpreds <- predict(spatmod, newdata = grid_yrs, type = "link", 
                        offset = grid_yrs$log_area, return_tmb_object = TRUE)
# Map the predictions (response scale) and spatial random field (link scale)
plot_map(spatialpreds$data, years = 2023, return = "response")
plot_map(spatialpreds$data, years = 2023, return = "spatial")
# Make the SD and CV for the predictions using simulations
sims <- predict(spatmod, newdata = grid_yrs, type = "link", offset = grid_yrs$log_area, nsim = 100)
spatialpreds$data$sd <- round(apply(exp(sims), 1, function(x) sd(x)), 2)
spatialpreds$data$cv <- round(apply(exp(sims), 1, function(x) sd(x) / mean(x)), 2)
# Map the SD and CV
plot_map(spatialpreds$data, years = 2023, return = "sd")
plot_map(spatialpreds$data, years = 2023, return = "cv")

# Create abundance index
spatialindex <- get_index(spatialpreds, area = 400, bias_correct = TRUE)
# Plot the index
plot_index(spatialindex)

# Plot for all years (Appendix A)
plot_map(spatialpreds$data, return = "response", axes = FALSE)
plot_map(spatialpreds$data, return = "spatial", axes = FALSE)
plot_map(spatialpreds$data, return = "sd", axes = FALSE)
plot_map(spatialpreds$data, return = "cv", axes = FALSE)

# An alternative model uses a smooth term on depth (Appendix B)
spatmodsmooth <- sdmTMB(count ~ 0+year_fac+s(depth)+gear+moult, data = lobster_data,
                        mesh = barrier_mesh, family = nbinom2(), offset = "log_area",
                        time = "year", spatial = "on", spatiotemporal = "off")
summary(spatmodsmooth)
smoothplot <- visreg(spatmodsmooth, xvar = "depth", data = lobster_data, gg = TRUE)
smoothplot+xlab("Depth (m)")+ylab("Link scale")+theme_light()
