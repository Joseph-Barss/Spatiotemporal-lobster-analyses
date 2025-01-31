
################ Simulation studies ################

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
library(knitr)
library(future)
library(foreach)
library(doFuture)
library(doRNG)
library(future.apply)
library(fda)

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

########### Simulation studies

# Fit the chosen model with the spatial random field
spatmod <- sdmTMB(count ~ 0+year_fac+logdepth+gear+moult, data = lobster_data,
                  mesh = barrier_mesh, family = nbinom2(), offset = "area_km2",
                  time = "year", spatial = "on", spatiotemporal = "off")
# Make predictions on a grid
spatialpreds <- predict(spatmod, newdata = grid_yrs, type = "link", 
                        offset = grid_yrs$area_km2, return_tmb_object = TRUE)
# Create abundance index
spatialindex <- get_index(spatialpreds, area = 400, bias_correct = TRUE)
# Extract parameters from the spatial model, making required transformations
pars <- get_pars(spatmod)
betas <- pars$b_j
ln_phi <- pars$ln_phi
ln_kappa <- pars$ln_kappa[1,]
ln_tau_O <- pars$ln_tau_O
phi <- exp(ln_phi)
range <- sqrt(8)/exp(ln_kappa)
sigma_O <- 1 / sqrt(4 * pi * exp(2 * ln_tau_O + 2 * ln_kappa))
omega_s_nodes <- pars$omega_s
# Use the index estimates as "true" abundance values (include scaled estimates)
sim_true_index <- spatialindex[, 1:2]
sim_true_index$est_scaled <- sim_true_index$est/100000000

# Make prediction grid
sim_grid <- make_prediction_grid(shape = shape, bathymetry = bathymetry,
                                 sweptarea = 0.01, depth_cutoff = 8)
sim_grid_yrs <- replicate_df(sim_grid, "year", unique(lobster_data$year))
sim_grid_yrs$year_fac <- as_factor(sim_grid_yrs$year)

# Simulate negative binomial distributed data sets 100 times
registerDoFuture()
registerDoRNG(10)
future::plan(multisession, workers = 10)
nsims <- 100
sim_nb_list <- foreach(i = 1:nsims, .options.future = list(seed = TRUE)) %dopar%
  get_nb_sims(lobster_data, barrier_mesh)
future::plan(sequential)

# For each simulated data set, fit different models and make indices
future::plan(multisession, workers = 10)
# Indices based on a negative binomial model
set.seed(7)
sim_nb_index_list <- future_lapply(sim_nb_list, get_nb_index_sims, 
                                   mesh = barrier_mesh, grid = sim_grid_yrs, 
                                   area = 400, future.seed = TRUE)

# Indices based on a Tweedie model
set.seed(7)
sim_tw_index_list <- future_lapply(sim_nb_list, get_tw_index_sims, 
                                   mesh = barrier_mesh, grid = sim_grid_yrs, 
                                   area = 4, future.seed = TRUE)

# Indices based on a delta-gamma model
set.seed(7)
sim_dg_index_list <- future_lapply(sim_nb_list, get_dg_index_sims, 
                                   mesh = barrier_mesh, grid = sim_grid_yrs, 
                                   area = 4, future.seed = TRUE)

# Indices based on a delta-lognormal model
set.seed(7)
sim_dln_index_list <- future_lapply(sim_nb_list, get_dln_index_sims, 
                                    mesh = barrier_mesh, grid = sim_grid_yrs, 
                                    area = 4, future.seed = TRUE)

future::plan(sequential)

# Simulating data at only scallop survey locations requires recreating the mesh
sim_scallop_template <- dplyr::filter(lobster_data, gear == "DREDGE")
scallop_mesh <- make_mesh(sim_scallop_template, c("X", "Y"), mesh = barrier_mesh$mesh)
scallop_barrier_mesh <- add_barrier_mesh(scallop_mesh, map_data, plot = TRUE, proj_scaling = 1000)
# Simulate negative binomial distributed data sets 100 times
registerDoFuture()
registerDoRNG(10)
future::plan(multisession, workers = 10)
nsims <- 100
scallop_list <- foreach(i = 1:nsims, .options.future = list(seed = TRUE)) %dopar%
  get_scallop_sims(sim_scallop_template, scallop_barrier_mesh)
future::plan(sequential)

# Simulating data at only lobster survey locations requires recreating the mesh
sim_lobster_template <- dplyr::filter(lobster_data, gear != "DREDGE")
lobster_mesh <- make_mesh(sim_lobster_template, c("X", "Y"), mesh = barrier_mesh$mesh)
lobster_barrier_mesh <- add_barrier_mesh(lobster_mesh, map_data, plot = TRUE, proj_scaling = 1000)
# Simulate negative binomial distributed data sets 100 times
registerDoFuture()
registerDoRNG(10)
future::plan(multisession, workers = 10)
nsims <- 100
lobster_list <- foreach(i = 1:nsims, .options.future = list(seed = TRUE)) %dopar%
  get_lobster_sims(sim_lobster_template, lobster_barrier_mesh)
future::plan(sequential)

# Fit models to simulated data sets and make indices
future::plan(multisession, workers = 10)
# Scallop survey locations only
set.seed(7)
sim_scallop_index_list <- future_lapply(scallop_list, get_scallop_index_sims, 
                                        mesh = scallop_barrier_mesh, 
                                        grid = sim_grid_yrs[sim_grid_yrs$year != 2020, ], 
                                        area = 400, future.seed = TRUE)
# Lobster survey locations only
set.seed(7)
sim_lobster_index_list <- future_lapply(lobster_list, get_nb_index_sims, 
                                        mesh = lobster_barrier_mesh, 
                                        grid = sim_grid_yrs, 
                                        area = 400, future.seed = TRUE)
future::plan(sequential)

# Include scaled estimate, lower, and upper columns for each index
# Also, scallop survey indices require explicitly numeric columns due to NA value
sim_nb_index_list <- lapply(sim_nb_index_list, function(x) dplyr::mutate(x, est_scaled = est/100000000,
                                                                         lwr_scaled = lwr/100000000,
                                                                         upr_scaled = upr/100000000))
sim_tw_index_list <- lapply(sim_tw_index_list, function(x) dplyr::mutate(x, est_scaled = est/100000000,
                                                                         lwr_scaled = lwr/100000000,
                                                                         upr_scaled = upr/100000000))
sim_dg_index_list <- lapply(sim_dg_index_list, function(x) dplyr::mutate(x, est_scaled = est/100000000,
                                                                         lwr_scaled = lwr/100000000,
                                                                         upr_scaled = upr/100000000))
sim_dln_index_list <- lapply(sim_dln_index_list, function(x) dplyr::mutate(x, est_scaled = est/100000000,
                                                                           lwr_scaled = lwr/100000000,
                                                                           upr_scaled = upr/100000000))
sim_lobster_index_list <- lapply(sim_lobster_index_list, function(x) dplyr::mutate(x, est_scaled = est/100000000,
                                                                                   lwr_scaled = lwr/100000000,
                                                                                   upr_scaled = upr/100000000))
sim_scallop_index_list <- lapply(sim_scallop_index_list, function(x) dplyr::mutate(x, est = as.numeric(est),
                                                                                   lwr = as.numeric(lwr),
                                                                                   upr = as.numeric(upr),
                                                                                   log_est = as.numeric(log_est),
                                                                                   se = as.numeric(se),
                                                                                   est_scaled = est/100000000,
                                                                                   lwr_scaled = lwr/100000000,
                                                                                   upr_scaled = upr/100000000))


# Plot the simulated indices against the true index for each scenario
plot_sim_indices(sim_true_index$est, sim_nb_index_list)
plot_sim_indices(sim_true_index$est, sim_tw_index_list)
plot_sim_indices(sim_true_index$est, sim_dg_index_list)
plot_sim_indices(sim_true_index$est, sim_dln_index_list)
plot_sim_indices(sim_true_index$est, sim_scallop_index_list)
plot_sim_indices(sim_true_index$est, sim_lobster_index_list)

# Get relative errors for each scenario
nb_rel_errors <- get_rel_errors(sim_true_index$est_scaled, sim_nb_index_list)
tw_rel_errors <- get_rel_errors(sim_true_index$est_scaled, sim_tw_index_list)
dg_rel_errors <- get_rel_errors(sim_true_index$est_scaled, sim_dg_index_list)
dln_rel_errors <- get_rel_errors(sim_true_index$est_scaled, sim_dln_index_list)
scallop_rel_errors <- get_rel_errors(sim_true_index$est_scaled, sim_scallop_index_list)
scallop_rel_errors[26, ] <- rep(FALSE, 100) # Needs this manual adjustment
lobster_rel_errors <- get_rel_errors(sim_true_index$est_scaled, sim_lobster_index_list)
# Calculate Mean Relative Errors for each year
nb_MREs <- round(rowMeans(nb_rel_errors), 3)
tw_MREs <- round(rowMeans(tw_rel_errors), 3)
dg_MREs <- round(rowMeans(dg_rel_errors), 3)
dln_MREs <- round(rowMeans(dln_rel_errors), 3)
scallop_MREs <- round(rowMeans(scallop_rel_errors), 3)
scallop_MREs[26] <- NA # Needs this manual adjustment
lobster_MREs <- round(rowMeans(lobster_rel_errors), 3)
MREs <- data.frame(year = 1995:2023, nb = nb_MREs, tw = tw_MREs,
                   dg = dg_MREs, dln = dln_MREs, 
                   scallop = scallop_MREs, lobster = lobster_MREs)

# Calculate Root Mean Square Errors for each year
nb_RMSEs <- get_RMSEs(sim_true_index$est_scaled, sim_nb_index_list)
tw_RMSEs <- get_RMSEs(sim_true_index$est_scaled, sim_tw_index_list)
dg_RMSEs <- get_RMSEs(sim_true_index$est_scaled, sim_dg_index_list)
dln_RMSEs <- get_RMSEs(sim_true_index$est_scaled, sim_dln_index_list)
scallop_RMSEs <- get_RMSEs(sim_true_index$est_scaled, sim_scallop_index_list)
lobster_RMSEs <- get_RMSEs(sim_true_index$est_scaled, sim_lobster_index_list)
RMSEs <- data.frame(year = 1995:2023, nb = nb_RMSEs, tw = tw_RMSEs,
                    dg = dg_RMSEs, dln = dln_RMSEs, 
                    scallop = scallop_RMSEs, lobster = lobster_RMSEs)

# Calculate 95% CI coverage rates for each year
nb_CRs <- get_CovRates(sim_true_index$est_scaled, sim_nb_index_list)
tw_CRs <- get_CovRates(sim_true_index$est_scaled, sim_tw_index_list)
dg_CRs <- get_CovRates(sim_true_index$est_scaled, sim_dg_index_list)
dln_CRs <- get_CovRates(sim_true_index$est_scaled, sim_dln_index_list)
scallop_CRs <- get_CovRates(sim_true_index$est_scaled, sim_scallop_index_list)
lobster_CRs <- get_CovRates(sim_true_index$est_scaled, sim_lobster_index_list)
CRs <- data.frame(year = 1995:2023, nb = nb_CRs, tw = tw_CRs,
                  dg = dg_CRs, dln = dln_CRs, 
                  scallop = scallop_CRs, lobster = lobster_CRs)

# Calculate 95% CI average widths for each year
nb_widths <- get_AveWidths(sim_nb_index_list)
tw_widths <- get_AveWidths(sim_tw_index_list)
dg_widths <- get_AveWidths(sim_dg_index_list)
dln_widths <- get_AveWidths(sim_dln_index_list)
scallop_widths <- get_AveWidths(sim_scallop_index_list)
lobster_widths <- get_AveWidths(sim_lobster_index_list)
AveWidths <- data.frame(year = 1995:2023, nb = nb_widths, tw = tw_widths,
                        dg = dg_widths, dln = dln_widths, 
                        scallop = scallop_widths, lobster = lobster_widths)

# Make perfomance metric tables (Appendix D)
# Simulation study 1
kable(MREs[,1:5], format = "latex", digits = 3, booktabs= TRUE,
      col.names = c("Year", "Negative binomial",
                    "Tweedie", "Delta-gamma",
                    "Delta-lognormal"))
kable(RMSEs[,1:5], format = "latex", digits = 3, booktabs= TRUE,
      col.names = c("Year", "Negative binomial",
                    "Tweedie", "Delta-gamma",
                    "Delta-lognormal"))
kable(CRs[,1:5], format = "latex", booktabs= TRUE,
      col.names = c("Year", "Negative binomial",
                    "Tweedie", "Delta-gamma",
                    "Delta-lognormal"))
kable(AveWidths[,1:5], format = "latex", digits = 3, booktabs= TRUE,
      col.names = c("Year", "Negative binomial",
                    "Tweedie", "Delta-gamma",
                    "Delta-lognormal"))
# Simulation study 2
kable(MREs[,c(1,2,6,7)], format = "latex", digits = 3, booktabs= TRUE,
      col.names = c("Year", "Both surveys",
                    "Scallop survey", "Lobster survey"))
kable(RMSEs[,c(1,2,6,7)], format = "latex", digits = 3, booktabs= TRUE,
      col.names = c("Year", "Both surveys",
                    "Scallop survey", "Lobster survey"))
kable(CRs[,c(1,2,6,7)], format = "latex", booktabs= TRUE,
      col.names = c("Year", "Both surveys",
                    "Scallop survey", "Lobster survey"))
kable(AveWidths[,c(1,2,6,7)], format = "latex", digits = 3, booktabs= TRUE,
      col.names = c("Year", "Both surveys",
                    "Scallop survey", "Lobster survey"))

# Plot performance metrics
# Begin by pivoting the tables of metrics to long format
MREs <- pivot_longer(MREs, cols = nb:lobster, 
                     names_to = "model", values_to = "MRE")
RMSEs <- pivot_longer(RMSEs, cols = nb:lobster, 
                      names_to = "model", values_to = "RMSE")
CRs <- pivot_longer(CRs, cols = nb:lobster, 
                    names_to = "model", values_to = "CR")
AWs <- pivot_longer(AveWidths, cols = nb:lobster, 
                    names_to = "model", values_to = "AW")
# Join them together
metrics <- MREs %>% full_join(RMSEs) %>% 
  full_join(CRs) %>% full_join(AWs) %>% 
  pivot_longer(cols = MRE:AW, names_to = "metric", values_to = "value") %>% 
  mutate(model = as_factor(model), metric = as_factor(metric))
# Make labels
labels <- c(
  MRE = "Mean relative error",
  RMSE = "Root mean square error",
  CR = "Coverage rate",
  AW = "Average width"
)
# Plots for simulation study 1
ggplot(filter(metrics, model %in% c("nb", "tw", "dg", "dln")), aes(year, value, colour = model))+
  geom_line()+
  facet_wrap(facets = vars(metric), scales = "free_y", 
             labeller = labeller(metric = labels))+
  scale_colour_discrete(name = "Model", labels = c("NB", "TW", "DG", "DLN"))+
  labs(x = "Year", y = "Metric value")+
  theme_light()
# Plots for simulation study 2
ggplot(filter(metrics, model %in% c("nb", "scallop", "lobster")), aes(year, value, colour = model))+
  geom_line()+
  facet_wrap(facets = vars(metric), scales = "free_y",
             labeller = labeller(metric = labels))+
  scale_colour_discrete(name = "Model", labels = c("Both", "Scallop", "Lobster"))+
  labs(x = "Year", y = "Metric value")+
  theme_light()

# Make functional boxplots (Appendix E)
# Boxplots for simulation study 1
par(mfrow = c(2,2))
fbplot(nb_rel_errors, ylim = c(-1, 4), xlab = "Year", ylab = "Relative Bias", 
       main = "Negative binomial", xaxt = "n")
axis(1, 1:29, 1995:2023)
abline(h = 0, lty = "dashed")
fbplot(tw_rel_errors, ylim = c(-1, 4), xlab = "Year", ylab = "Relative Bias", 
       main = "Tweedie", xaxt = "n")
axis(1, 1:29, 1995:2023)
abline(h = 0, lty = "dashed")
fbplot(dg_rel_errors, ylim = c(-1, 4), xlab = "Year", ylab = "Relative Bias",
       main = "Delta-gamma", xaxt = "n")
axis(1, 1:29, 1995:2023)
abline(h = 0, lty = "dashed")
fbplot(dln_rel_errors, ylim = c(-1, 4), xlab = "Year", ylab = "Relative Bias", 
       main = "Delta-lognormal", xaxt = "n")
axis(1, 1:29, 1995:2023)
abline(h = 0, lty = "dashed")
# Boxplots for simulation study 2
par(mfrow = c(3,1))
fbplot(nb_rel_errors, ylim = c(-1, 1.2), xlab = "Year", ylab = "Relative Bias", 
       main = "Both surveys", xaxt = "n")
axis(1, 1:29, 1995:2023)
abline(h = 0, lty = "dashed")
fbplot(scallop_rel_errors, ylim = c(-1, 1.2), xlab = "Year", ylab = "Relative Bias", 
       main = "Scallop survey", xaxt = "n")
axis(1, 1:29, 1995:2023)
abline(h = 0, lty = "dashed")
fbplot(lobster_rel_errors, ylim = c(-1, 1.2), xlab = "Year", ylab = "Relative Bias",
       main = "Lobster survey", xaxt = "n")
axis(1, 1:29, 1995:2023)
abline(h = 0, lty = "dashed")



