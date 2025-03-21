
############### Functions used in the analyses ################

######## Data processing functions

# Reads in and processes LFA shapefile
make_shape <- function(){
  lfas <- readRDS("LFAPolysSF.rds")%>% 
    filter(LFA %in% 34:38) %>% 
    mutate(area = as.numeric(st_area(geometry))/1000) 
  shape <- lfas %>% 
    filter(area > 200000) %>% 
    st_union() %>% 
    st_simplify(dTolerance = 500) %>% 
    st_transform(crs_UTM20N)
  islands <-  lfas %>% 
    filter(area <= 200000 & area > 2000) %>% 
    st_union() %>% 
    st_simplify(dTolerance = 500) %>% 
    st_transform(crs_UTM20N)
  shape <- st_difference(shape, islands)
  return(shape)
}

# Reads in and processes all trawl survey data and joins with depth data
read_count_data <- function(shape, bathymetry){
  ILTS <- read_csv("TotalLobsterCatch_ILTS.csv", col_types = "icciTdddci") %>% 
    dplyr::select(TRIP_ID, SET_NO, SET_DATE, YEAR, SET_LONG, SET_LAT, GEAR, sweptArea, NUM_CAUGHT) %>%
    rename(trip = TRIP_ID, tow = SET_NO, date = SET_DATE, year = YEAR, long = SET_LONG, lat = SET_LAT, 
           gear = GEAR, area_km2 = sweptArea, count = NUM_CAUGHT) %>% 
    mutate(month = as_factor(month(date)), .after = year)
  ISAS <- read_csv("scallop_dredge_lobster_abund_data_July2024.csv", col_types = "cccTifiddcid") %>% 
    dplyr::select(CRUISE, TOW_NO, TOW_DATE, year, month, SLONG, SLAT, GEAR, sample_area_m2, abun_raw) %>% 
    rename(trip = CRUISE, tow = TOW_NO, date = TOW_DATE, long = SLONG, lat = SLAT, gear = GEAR, count = abun_raw) %>% 
    mutate(area_km2 = sample_area_m2/1000000, .keep = "unused", .after = gear)
  dat <- bind_rows(ILTS, ISAS) %>% 
    arrange(year) %>% 
    mutate(year_fac = as_factor(year), .after = year) %>% 
    mutate(gear = fct(gear, levels = c("NEST", "280 BALLOON", "DREDGE")), 
           moult = fct_rev(fct_collapse(month, Post = c("8", "9", "10"), Pre = c("5", "6", "7"))),
           .before = area_km2) %>% 
    mutate(log_area = log(area_km2), .after = area_km2) %>% 
    mutate(count_per_km2 = count/area_km2)
  # Add UTM columns to data
  crs_UTM20N <- 32620
  dat <- add_utm_columns(dat, ll_names = c("long", "lat"), utm_crs = crs_UTM20N) %>% 
    mutate(xmetres = X*1000, ymetres = Y*1000)
  dat_sf <- st_as_sf(dat, coords = c("xmetres", "ymetres")) 
  st_crs(dat_sf) <- crs_UTM20N
  dat_sf <- st_intersection(dat_sf, shape)
  dat <- as.data.frame(dat_sf) %>% 
    mutate(xmetres = X*1000, ymetres = Y*1000) %>% 
    dplyr::select(-geometry)
  # Get depth and add to data
  sp_dat <- terra::vect(dat, geom = c("xmetres", "ymetres"))
  sample_depths <- terra::extract(bathymetry, sp_dat)
  dat <- dat %>% mutate(depth = -sample_depths[,2], logdepth = log(depth)) %>% 
    dplyr::select(trip:lat, X, Y, xmetres:logdepth, gear, moult, area_km2:count_per_km2)
  dat <- dat[!is.na(dat$depth) & dat$depth > 0,]
  return(dat)
}

# Makes a prediction grid for the study area
make_prediction_grid <- function(shape, bathymetry, density_m = 2000, sweptarea = 1, depth_cutoff){
  grid <- shape %>% 
    st_make_grid(cellsize = c(density_m, density_m), what = "centers") %>% # 4km^2 size cells default
    st_intersection(shape) %>% 
    st_coordinates() %>% 
    as.data.frame() %>% 
    rename(xmetres = X, ymetres = Y) %>% 
    mutate(X = xmetres/1000, Y = ymetres/1000, gear = "NEST", moult = "Post")
  sp_grid <- vect(grid, geom = c("xmetres", "ymetres"))
  grid_depths <- terra::extract(bathymetry, sp_grid)
  grid <- mutate(grid, "depth" = -grid_depths[, 2], .before = gear)
  grid <- mutate(grid, "area_km2" = sweptarea)
  grid <- mutate(grid, "log_area" = log(sweptarea))
  grid <- grid[!is.na(grid$depth) & grid$depth > depth_cutoff,] # don't want too shallow
  grid <- mutate(grid, "logdepth" = log(depth), .after = depth)
  return(grid)
}

# Makes a barrier mesh 
make_barrier_mesh <- function(dat, map){
  mesh_dat <- as.matrix(dat[dat$year == 2023, c("X", "Y")])
  bnd <- fm_nonconvex_hull(mesh_dat, convex = -0.075)
  nonconvex_mesh <- fm_mesh_2d(
    boundary = bnd,
    cutoff = 4,
    max.edge = c(20, 60),
    offset = c(20, 60)
  )
  mesh <- make_mesh(dat, c("X", "Y"), mesh = nonconvex_mesh)
  plot(mesh)
  barrier_mesh <- add_barrier_mesh(mesh, map, plot = TRUE, proj_scaling = 1000)
  return(barrier_mesh)
}

######## Simulation study functions

# Simulates negative binomial-distributed count data at all sampling locations
get_nb_sims <- function(template, mesh){
  sim_dat <- sdmTMB_simulate(~0+year_fac+logdepth+gear+moult, template, mesh, family = nbinom2(link = "log"),
                             offset = "log_area", time = "year", B = betas, range = range, sigma_O = sigma_O,
                             phi = phi,
                             spatiotemporal = "off",
                             fixed_re = list(omega_s = omega_s_nodes
                             ))
  sim_dat <- sim_dat[c(1:7, 37)]
  sim_dat$year_fac <- template$year_fac
  sim_dat$depth <- template$depth
  sim_dat$gear <- template$gear
  sim_dat$moult <- template$moult
  sim_dat <- rename(sim_dat, count = observed)
  sim_dat$area_km2 <- template$area_km2
  sim_dat$log_area <- template$log_area
  sim_dat$count_per_km2 <- sim_dat$count/sim_dat$area_km2
  return(sim_dat)
}

# Simulates negative binomial-distributed count data at scallop survey locations
get_scallop_sims <- function(template, mesh){
  sim_dat <- sdmTMB_simulate(~0+year_fac+logdepth+moult, template, mesh, family = nbinom2(link = "log"),
                             offset = "log_area", time = "year", 
                             B = betas[c(1:25, 27:30, 33)], # No gear, and 2020 missing 
                             range = range, sigma_O = sigma_O, 
                             phi = phi,
                             spatiotemporal = "off",
                             fixed_re = list(omega_s = omega_s_nodes
                             ))
  sim_dat <- sim_dat[c(1:7, 36)]
  sim_dat$year_fac <- template$year_fac
  sim_dat$depth <- template$depth
  sim_dat$gear <- template$gear
  sim_dat$moult <- template$moult
  sim_dat <- rename(sim_dat, count = observed)
  sim_dat$area_km2 <- template$area_km2
  sim_dat$count_per_km2 <- sim_dat$count/sim_dat$area_km2
  return(sim_dat)
}

# Simulates negative binomial-distributed count data at lobster survey locations
get_lobster_sims <- function(template, mesh){
  sim_dat <- sdmTMB_simulate(~0+year_fac+logdepth+gear+moult, template, mesh, family = nbinom2(link = "log"),
                             offset = "log_area", time = "year", B = betas[c(1:31, 33)], range = range, sigma_O = sigma_O,
                             phi = phi,
                             spatiotemporal = "off",
                             fixed_re = list(omega_s = omega_s_nodes
                             ))
  sim_dat <- sim_dat[c(1:7, 37)]
  sim_dat$year_fac <- template$year_fac
  sim_dat$depth <- template$depth
  sim_dat$gear <- template$gear
  sim_dat$moult <- template$moult
  sim_dat <- rename(sim_dat, count = observed)
  sim_dat$area_km2 <- template$area_km2
  sim_dat$count_per_km2 <- sim_dat$count/sim_dat$area_km2
  return(sim_dat)
}

# Simulates abundance index using negative binomial model
get_nb_index_sims <- function(dat, mesh, grid, area){
  sim_mod_nb <- sdmTMB(count ~ 0+year_fac+logdepth+gear+moult, data = dat,
                       mesh = mesh, family = nbinom2(), offset = "log_area", time = "year",
                       spatial = "on", spatiotemporal = "off", do_index = TRUE, 
                       predict_args = list(newdata = grid, offset = grid$area_km2), index_args = list(area = area))
  gc()
  sim_nb_index <- get_index(sim_mod_nb, bias_correct = TRUE)
  rm(sim_mod_nb)
  gc()
  return(sim_nb_index)
}

# Simulates abundance index using negative binomial model, with adjustment for missing year
get_scallop_index_sims <- function(dat, mesh, grid, area){
  sim_mod_scallop <- sdmTMB(count ~ 0+year_fac+logdepth+moult, data = dat,
                            mesh = mesh, family = nbinom2(), offset = "log_area", time = "year",
                            spatial = "on", spatiotemporal = "off", do_index = TRUE, 
                            predict_args = list(newdata = grid, offset = grid$area_km2), index_args = list(area = area))
  gc()
  scallop_index <- get_index(sim_mod_scallop, bias_correct = TRUE)
  rm(sim_mod_scallop)
  scallop_index <- rbind(scallop_index[1:25, ], c(2020, NA, NA, NA, NA, NA, "index"), scallop_index[26:28,])
  dplyr::mutate(scallop_index, est = as.numeric(est),
                lwr = as.numeric(lwr),
                upr = as.numeric(upr),
                log_est = as.numeric(log_est),
                se = as.numeric(se))
  rownames(scallop_index) <- 1:nrow(scallop_index)
  gc()
  return(scallop_index)
}

# Simulates abundance index using Tweedie model
get_tw_index_sims <- function(dat, mesh, grid, area){
  sim_mod_tw <- sdmTMB(count_per_km2 ~ 0+year_fac+logdepth+gear+moult, data = dat,
                       mesh = mesh, family = tweedie(), time = "year",
                       spatial = "on", spatiotemporal = "off", do_index = TRUE, 
                       predict_args = list(newdata = grid), index_args = list(area = area))
  gc()
  sim_tw_index <- get_index(sim_mod_tw, bias_correct = TRUE)
  rm(sim_mod_tw)
  gc()
  return(sim_tw_index)
}

# Simulates abundance index using delta-gamma model
get_dg_index_sims <- function(dat, mesh, grid, area){
  sim_mod_dg <- sdmTMB(count_per_km2 ~ 0+year_fac+logdepth+gear+moult, data = dat,
                       mesh = mesh, family = delta_gamma(), time = "year",
                       spatial = "on", spatiotemporal = "off", do_index = TRUE, 
                       predict_args = list(newdata = grid), index_args = list(area = area))
  gc()
  sim_dg_index <- get_index(sim_mod_dg, bias_correct = TRUE)
  rm(sim_mod_dg)
  gc()
  return(sim_dg_index)
}

# Simulates abundance index using delta-lognormal model
get_dln_index_sims <- function(dat, mesh, grid, area){
  sim_mod_dln <- sdmTMB(count_per_km2 ~ 0+year_fac+logdepth+gear+moult, data = dat,
                        mesh = mesh, family = delta_lognormal(), time = "year",
                        spatial = "on", spatiotemporal = "off", do_index = TRUE, 
                        predict_args = list(newdata = grid), index_args = list(area = area))
  gc()
  sim_dln_index <- get_index(sim_mod_dln, bias_correct = TRUE)
  rm(sim_mod_dln)
  gc()
  return(sim_dln_index)
}

# Extract the relative errors of the simulated indices vs. the true index
get_rel_errors <- function(true, index_list){
  est_list <- lapply(index_list, function(x) dplyr::select(x, est_scaled))
  rel_errors <- lapply(est_list, function(x) (x-true)/true)
  rel_errors <- matrix(unlist(rel_errors), ncol = 100, byrow = FALSE)
}

# Extract the RMSEs
get_RMSEs <- function(true, index_list){
  est_list <- lapply(index_list, function(x) dplyr::select(x, est_scaled))
  sq_errors <- lapply(est_list, function(x) (x-true)^2)
  sq_errors <- matrix(unlist(sq_errors), ncol = 100, byrow = FALSE)
  RMSEs <- rowMeans(sq_errors)
}

# Extract the 95% confidence interval coverage rates
get_CovRates <- function(true, index_list){
  lower_bounds <- lapply(index_list, function(x) dplyr::select(x, lwr_scaled))
  upper_bounds <- lapply(index_list, function(x) dplyr::select(x, upr_scaled))
  contains_true <- mapply(function(x, y) x <= true & y >= true, lower_bounds, upper_bounds, SIMPLIFY = FALSE)
  contains_true <- matrix(unlist(contains_true), ncol = 100, byrow = FALSE)
  coverage_rates <- 100*rowMeans(contains_true)
}

# Extract the 95% confidence interval average widths
get_AveWidths <- function(index_list){
  lower_bounds <- lapply(index_list, function(x) dplyr::select(x, lwr_scaled))
  upper_bounds <- lapply(index_list, function(x) dplyr::select(x, upr_scaled))
  widths <- mapply(function(x, y) y-x, lower_bounds, upper_bounds, SIMPLIFY = FALSE)
  widths <- matrix(unlist(widths), ncol = 100, byrow = FALSE)
  ave_widths <- rowMeans(widths)
}

########### Plotting functions

# Plot the residuals with a histogram, QQ-plot, or map
plot_resids_spatial <- function(data, plot_type = c("histogram", "qq-plot", "map")){
  if(plot_type == "histogram"){
    plot <- ggplot(data, aes(spresids))+
      geom_histogram()+
      labs(x = "Randomized quantile residuals")
  } else if(plot_type == "qq-plot"){
    plot <- ggplot(data, aes(sample = spresids))+
      geom_qq()+
      geom_qq_line()+
      labs(x = "Theoretical", y = "Sample")
  } else {
    plot <- ggplot(map_fundy) + geom_sf() + 
      geom_point(data = data, mapping = aes(x = xmetres, y = ymetres, col = spresids), alpha = 0.5, size = 0.5) + 
      scale_colour_gradient2() +
      labs(x = "Longitude", y = "Latitude", col = "Residual")
  }
  plot <- plot+theme_light()
  return(plot)
}

# Map the predicted abundance, standard deviation, coefficient of variation, or random field values
plot_map <- function(dat, years = "all", 
                     return = c("link", "response", "spatial", "spatiotemporal", 
                                "sd", "cv", "bayes", "bayesSD"),
                     axes = TRUE) {
  if(return == "link"){
    dat$response <- dat$est
    label = "Log counts"
    trans = "identity"
  } else if (return == "spatial") {
    dat$response <- dat$omega_s
    label = "Spatial RF"
    trans = "identity"
  } else if (return == "spatiotemporal"){
    dat$response <- dat$epsilon_st
    label = "Spatiotemporal RF"
    trans = "identity"
  } else if (return == "sd"){
    dat$response <- dat$sd
    label = "Standard Deviation"
    trans = "identity"
  } else if (return == "cv"){
    dat$response <- dat$cv
    label = "Coefficient of Variation"
    trans = "identity"
  } else if (return == "bayes") {
    dat$response <- dat$est
    label = "Response (counts)"
    trans = "log10"
  } else if (return == "bayesSD") {
    dat$response <- dat$sd
    label = "Standard Deviation"
    trans = "log10"
  } else {
    dat$response <- exp(dat$est)
    label = "Response (counts)"
    trans = "log10"
  }
  if(is.numeric(years)){
    dat <- filter(dat, year %in% years)
  }
  p <- ggplot(map_fundy) +
    geom_sf() +
    geom_tile(data = dat, aes(xmetres, ymetres, fill = response)) +
    facet_wrap(~year) +
    coord_sf()+
    scale_fill_viridis_c(option = "B", direction = -1, transform = trans)+
    labs(x = "Longitude", y = "Latitude", fill = label)+
    theme_light()
  if(axes == FALSE){
    p <- p+theme(axis.line=element_blank(),axis.text.x=element_blank(),
                 axis.text.y=element_blank(),axis.ticks=element_blank(),
                 axis.title.x=element_blank(),
                 axis.title.y=element_blank())
  }
  return(p)
}

# Plot the abundance index
plot_index <- function(dat){
  p <- ggplot(dat, aes(year, est)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4) +
    labs(x = "Year", y = "Index value")+
    scale_y_continuous(labels=function(x)x/100000000)+
    theme_light()
  return(p)
}

# Plot the simulated indices against the true index
plot_sim_indices <- function(true, index_list){
  est_list <- lapply(index_list, function(x) dplyr::select(x, est))
  dat <- data.frame(year = 1995:2023, true = true)
  for(i in seq(length(est_list))){
    est <- est_list[[i]]
    names(est) <- paste("est", as.character(i))
    dat <- dplyr::mutate(dat, est)
  }
  gc()
  dat <- pivot_longer(dat, !year, names_to = "series", values_to = "index")
  p <- ggplot()+
    geom_line(data = dat[dat[ , 2] != "true", ], mapping = aes(year, index, col = series, alpha = 0.2)) +
    geom_line(data = dat[dat[ , 2] == "true", ], mapping = aes(year, index), colour = "black", linewidth = 2)+
    scale_y_continuous(breaks = c(0, 1e+08, 2e+08, 3e+08, 4e+08, 5e+08),
                       labels = c(0, 1, 2, 3, 4, 5))+
    theme_light()+
    theme(legend.position="none")+
    labs(x = "Year", y = "Index value")
  return(p)
}