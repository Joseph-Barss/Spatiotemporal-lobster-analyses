
##################### Bayesian analysis ################

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
library(tmbstan)
library(bayesplot)
library(parallel)

# Source functions
source("functions.R")

############ Data processing 

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

########## Bayesian model

#### Fit model with priors
mod_for_tmbstan <- sdmTMB(count ~ 0+year_fac+logdepth+gear+moult, data = lobster_data,
                          mesh = barrier_mesh, family = nbinom2(), offset = "log_area",
                          time = "year", spatial = "on", spatiotemporal = "off",
                          bayesian = TRUE,
                          priors = sdmTMBpriors(
                            b = normal(c(rep(10, 29), -1, 0, 0, 0), c(rep(2, 29), 1, 3, 3, 3)),
                            phi = halfnormal(0, 10),
                            matern_s = pc_matern(range_gt = 50, sigma_lt = 35)
                          ))
# Fit the model in tmbstan (four chains in parallel)
options(mc.cores = parallel::detectCores()-1)
stanmod <- tmbstan::tmbstan(mod_for_tmbstan$tmb_obj, iter = 3000, warmup = 1000,
                            chains = 4, seed = 510
)
options(mc.cores = 1)
# Show model summary
summary(stanmod)
# Make trace plots for the model (Appendix C)
pars_plot1 <- c("b_j[1]","b_j[2]","b_j[3]","b_j[4]","b_j[5]",
                "b_j[6]","b_j[7]","b_j[8]","b_j[9]","b_j[10]",
                "b_j[11]","b_j[12]","b_j[13]","b_j[14]","b_j[15]",
                "b_j[16]","b_j[17]","b_j[18]")
pars_plot2 <- c("b_j[19]","b_j[20]",
                "b_j[21]","b_j[22]","b_j[23]","b_j[24]","b_j[25]",
                "b_j[26]","b_j[27]","b_j[28]","b_j[29]","b_j[30]",
                "b_j[31]","b_j[32]","b_j[33]",
                "ln_tau_O", "ln_kappa", "ln_phi")
bayesplot::mcmc_trace(stanmod, pars = pars_plot1)
bayesplot::mcmc_trace(stanmod, pars = pars_plot2)
# Extract posterior simulations
post <- rstan::extract(stanmod)
# Transformations for shape and random field parameters
ln_phi_b <- post$ln_phi
ln_kappa_b <- post$ln_kappa
ln_tau_O_b <- post$ln_tau_O
phi_b <- exp(ln_phi_b)
range_b <- sqrt(8)/exp(ln_kappa_b)
sigma_O_b <- 1 / sqrt(4 * pi * exp(2 * ln_tau_O_b + 2 * ln_kappa_b))
mean(phi_b)
mean(range_b)
mean(sigma_O_b)

# Make predictions and SDs on the grid using the Bayesian model
set.seed(300)
samps <- sdmTMBextra::extract_mcmc(stanmod)
pred_grid_bayes23 <- dplyr::filter(grid_yrs, year == 2023)
pred_bayes <- predict(mod_for_tmbstan, newdata = grid_yrs, 
                      offset = grid_yrs$log_area, mcmc_samples = samps)
pred_grid_bayes23$est <- apply(exp(pred_bayes[229237:237423, ]), 1, mean)
pred_grid_bayes23$sd <- apply(exp(pred_bayes[229237:237423, ]), 1, sd)


# Map the predictions and SDs
plot_map(pred_grid_bayes23, years = 2023, return = "bayes")
plot_map(pred_grid_bayes23, years = 2023, return = "bayesSD")

# Make Bayesian index
pred_bayes <- as_tibble(pred_bayes, rownames = "Year")
pred_abun <- pred_bayes %>%
  mutate(across("V1":"V8000", exp), .keep = "unused") %>% 
  group_by(Year) %>% 
  summarise(across("V1":"V8000", sum)) %>% 
  dplyr::select(-Year)
bayes_index <- data.frame(year = 1995:2023,
                          est = 400*apply(pred_abun, 1, mean),
                          lwr = 400*apply(pred_abun, 1, quantile, probs = 0.025),
                          upr = 400*apply(pred_abun, 1, quantile, probs = 0.975))
# Plot Bayesian index
ggplot(bayes_index, aes(year, est)) + geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4) +
  labs(x = "Year", y = "Index value")+
  scale_y_continuous(labels = c(0, 1, 2, 3, 4, 5))+
  theme_light()
