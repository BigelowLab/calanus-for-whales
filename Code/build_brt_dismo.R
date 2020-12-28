# Camille Ross
# 10 June, 2020
# Purpose: Run a boosted regression tree model for zooplankton species with variable environmental covariates

# -------- Load libraries --------
require(dplyr)
require(readr)
require(dismo)
require(gbm)
require(viridis)
require(Metrics)
require(mapdata)
require(maps)
require(ggplot2)
require(stats)

# -------- Source outside files --------
# Source data formatting function
source("./projects/calanus4whales/calanus_for_whales/format_model_data.R")
# Source covariate loading function
source("./projects/calanus4whales/calanus_for_whales/load_covars.R")
# Source data binding function
source("./projects/calanus4whales/calanus_for_whales/bind_years.R")

# -------- Main function --------
#'@param version <chr> version of model
#'@param fp_md <chr> file path to formatted data used in model
#'@param species <chr> species to model; choices are "cfin", "ctyp", or "pseudo"
#'@param fp_covars <chr> file path to environmental covariate data
#'@param env_covars <vector> vector of covariates to include in the model
#'@param years <vectors> years for which to run the model
#'@param fp_out <chr> file path save the data to 
#'@param format_data <logical> if true, data is formatted within function; only used if model_data is NULL
#'@param fp_zpd <chr> filepath to the zooplankton database if data is formatted within function
build_brt <- function(version, fp_md, species, fp_covars, env_covars, threshold, years, proj_year, fp_out,
                      format_data = FALSE, fp_zpd = NULL) {
  
  # -------- Create output directories --------
  dir.create(fp_out, showWarnings = FALSE) 
  dir.create(file.path(fp_out, species), showWarnings = FALSE)
  dir.create(file.path(fp_out, species, "BRTs"), showWarnings = FALSE)
  dir.create(file.path(fp_out, species, "BRTs", version), showWarnings = FALSE)
  # -------- Create directory for projections --------
  dir.create(file.path(fp_out, species, "BRTs", version, "Projections"), showWarnings = FALSE)
  # -------- Create directory for plots --------
  dir.create(file.path(fp_out, species, "BRTs", version, "Plots"), showWarnings = FALSE)
  # -------- Create directory for interactions --------
  dir.create(file.path(fp_out, species, "BRTs", version, "Interactions"), showWarnings = FALSE)
  # -------- Create directory for evaluations --------
  dir.create(file.path(fp_out, species, "BRTs", version, "Evals"), showWarnings = FALSE)
  
  # -------- Format model data --------
  if (format_data) {
    # Format zooplankton and environmental covariate data
    format_model_data(fp_data = fp_zpd, fp_covars = fp_covars, env_covars = "all", years = years, fp_out = fp_md)
  }
  
  # -------- Load model data --------
  if (length(years) == 1) {
    md <- readr::read_csv(file.path(fp_md, paste0(years[1], ".csv"))) %>% 
      dplyr::filter(dataset %in% biomod_dataset)
  } else {
    md <- bind_years(fp = file.path(fp_md), years = years) %>%
      dplyr::filter(dataset %in% biomod_dataset)
  }
  
  # -------- Compute anomaly --------
  if (species == "cfin") {
    md <- md %>% dplyr::group_by(dataset) %>%
      dplyr::mutate(mean = mean(log10(`cfin_CV_VI` + 1), na.rm = TRUE),
                    sd = sd(log10(`cfin_CV_VI` + 1), na.rm = TRUE),
                    anomaly = (log10(`cfin_CV_VI` + 1) - mean) / sd) %>%
      dplyr::ungroup()
  } else if (species == "ctyp") {
    md <- md %>% dplyr::group_by(dataset) %>%
      dplyr::mutate(mean = mean(log10(`ctyp_total` + 1), na.rm = TRUE),
                    sd = sd(log10(`ctyp_total` + 1), na.rm = TRUE),
                    anomaly = (log10(`ctyp_total` + 1) - mean) / sd) %>%
      dplyr::ungroup()
  } else if (species == "pcal") {
    md <- md %>% dplyr::group_by(dataset) %>%
      dplyr::mutate(mean = mean(log10(`pcal_total` + 1), na.rm = TRUE),
                    sd = sd(log10(`pcal_total` + 1), na.rm = TRUE),
                    anomaly = (log10(`pcal_total` + 1) - mean) / sd) %>%
      dplyr::ungroup()
  }
  
  # -------- Take the log of count data --------
  if (species == "cfin") {
    md$abund <- as.data.frame(md[paste0(species, "_CV_VI")])$cfin_CV_VI
  } else if (species == "ctyp") {
    md$abund <- as.data.frame(md[paste0(species, "_total")])$ctyp_total
  } else if (species == "pcal") {
    md$abund <- as.data.frame(md[paste0(species, "_total")])$pcal_total
  }
  
  # -------- Exclude NAs and select columns --------
  md <- md %>% dplyr::select(lat, lon, year, month, abund, wind, fetch, chl, int_chl, bots, bott, sss, sst, lag_sst, uv, bat, dist, slope) %>%
    as.data.frame() %>%
    na.exclude() %>%
    dplyr::mutate(season = if_else(month %in% c(1:3), 1,
                                   if_else(month %in% c(4:6), 2,
                                           if_else(month %in% c(7:9), 3, 4))))
  
  # -------- Take log of bathymetry and chlorophyll --------
  md$chl <- log(abs(md$chl))
  md$int_chl <- log(abs(md$int_chl))
  md$bat <- log(abs(md$bat))
  
  
  # -------- Zero out depths below 1000 --------
  md$wind[md$bat >= 1000] <- 0
  md$fetch[md$bat >= 1000] <- 0
  md$chl[md$bat >= 1000] <- 0
  md$bots[md$bat >= 1000] <- 0
  md$bott[md$bat >= 1000] <- 0
  md$sss[md$bat >= 1000] <- 0
  md$sst[md$bat >= 1000] <- 0
  md$uv[md$bat >= 1000] <- 0
  md$bat[md$bat >= 1000] <- 0
  
  # -------- Initialize brt arguments --------
  brt_args <- list(n.trees = 1000, interaction.depth = 5,
                   n.minobsinnode = 3,
                   shrinkage = 0.001, bag.fraction = 0.5, 
                   train.fraction = 1, cv.folds = 10)
  
  # -------- Load world map data --------
  worldmap <- ggplot2::map_data("world")
  
  # -------- Initialize projection year --------
  proj_year <- 2017
  
  # -------- Determine which covariates to use --------
  if ("all" %in% env_covars) {
    env_covars <- c("wind", "fetch", "chl",
                    "bots", "bott", "sss",
                    "sst", "uv", "bat")
  }
  
  env_covars <- c("chl", "sss", "sst")
  
  # -------- Loop over months --------
  for (i in 4:12) {
    # -------- Isolate month data --------
    month_md <- md %>% dplyr::filter(month == i) %>%
      mutate(abund = if_else(abund < threshold, 0, 1))
    
    # -------- Divide into training and testing data --------
    # indices <- sample.int(n = nrow(month_md), size = floor(x = 0.7*nrow(month_md)), replace = FALSE)
    # train <- month_md[indices,]
    # test <- month_md[!duplicated(rbind(train, month_md))[-(1:nrow(train))],]
    
    # -------- Build BRT with all covariates --------
    brt_sdm <- dismo::gbm.step(data = month_md, gbm.x = c(
                                                          "bots", "bott", "sss",
                                                           "bat", "int_chl",
                                                          "lag_sst"), gbm.y = 5,
                               family = "gaussian", tree.complexity = 5,
                               learning.rate = 0.001, bag.fraction = 0.5,
                               n.minobsinnode = 2, nTrain = 1)
    
    
    # -------- Plot results --------
    plot(brt_sdm)
    # -------- Load summary of model --------
    png(file.path(fp_out, species, "BRTs", version, "Plots", paste0("var_cont_", proj_year, "_", i, ".png")))
    summary(brt_sdm)
    dev.off()
    
    # -------- Test interactions between covariates --------
    find_int <- as.data.frame(dismo::gbm.interactions(brt_sdm)$interactions) %>% 
       write_csv(path = file.path(fp_out, species, "BRTs", version, "Interactions", paste0("var_int_", proj_year, "_", i, ".csv")))
    
    # -------- Determine which covariates to use --------
    # brt_simp <- gbm.simplify(brt_sdm, n.drops = 5)
    # 
    # brt_simp$pred.list
    
    # -------- Plot response curves --------
    #dismo::gbm.plot(brt_sdm, n.plots=4, plot.layout=c(2, 2))
    
    # -------- Plot fitted values in relation to predictors --------
    # gbm.plot.fits(brt_sdm)
    
    # -------- Predict onto testing data --------
    # tree_list <- seq(1000, 3000, by=100)
    # pred <- predict.gbm(brt_sdm, test, n.trees = tree_list, "response")
    
    # -------- Load environmental covariates for projection --------
    covars <- load_covars(fp_covars = fp_covars, year = proj_year, month = i,
                          env_covars = env_covars,
                          as_raster = TRUE, lag_sst = TRUE)
    
    # Load bathymetry layer 
    bat <- load_covars(fp_covars = fp_covars, year = proj_year, month = i,
                       env_covars = "bat",
                       as_raster = TRUE)
    
    # -------- Zero out below 1000m --------
    if ("wind" %in% env_covars) {
      covars$wind[bat >= 1000] <- 0
    }
    if ("fetch" %in% env_covars) {
      covars$fetch[bat >= 1000] <- 0
    }
    if ("chl" %in% env_covars) {
      covars$chl[bat >= 1000] <- 0
    }
    if ("bots" %in% env_covars) {
      covars$bots[bat >= 1000] <- 0
    }
    if ("bott" %in% env_covars) {
      covars$bott[bat >= 1000] <- 0
    }
    if ("sss" %in% env_covars) {
      covars$sss[bat >= 1000] <- 0
    }
    if ("sst" %in% env_covars) {
      covars$sst[bat >= 1000] <- 0
    }
    if ("uv" %in% env_covars) {
      covars$uv[bat >= 1000] <- 0
    }
    if ("bat" %in% env_covars) {
      covars$bat[bat >= 1000] <- 0
    }
    
    # -------- Take log of bathymetry and chlorophyll --------
    if ("chl" %in% env_covars) {
      covars$chl <- log(abs(covars$chl))
    }
    if ("bat" %in% env_covars) {
      covars$bat <- log(abs(covars$bat))
    }
    
    # -------- Project model onto covariates --------
    proj <- raster::predict(covars, brt_sdm, n.trees = brt_sdm$gbm.call$best.trees,
                            filename = file.path(fp_out, species, "BRTs", version, "Projections", paste0("proj_", proj_year, "_", i)), progress = "text",
                            overwrite = TRUE)
    
    # Zero out projections below 1000m
    proj[bat >= 1000] <- 0
    
    # Convert proj to data frame
    proj_df <- as.data.frame(proj, xy = TRUE)
    # Add column names
    names(proj_df) <- c("x", "y", "pred")
    
    # -------- Plot projection --------
    ggplot() + 
      #Add projection data
      geom_tile(data = proj_df, aes(x, y, fill = pred)) +
      #Add projection color gradient and label
      scale_fill_gradientn(colors = inferno(500), limits = c(min(md$abund), max(md$abund)), na.value = "white") +
      labs(x = "", 
           y = "") +
      # geom_point(data = month_md, aes(lon, lat, color = abund)) +
      # scale_color_gradientn(colors = inferno(500), limits = c(min(md$abund), max(md$abund)), na.value = "white") +
      # Add world map data
      geom_polygon(data = worldmap, aes(long, lat, group = group), fill = NA, colour = "gray43") +
      coord_quickmap(xlim = c(round(min(proj_df$x)), round(max(proj_df$x))), 
                     ylim = c(round(min(proj_df$y)), round(max(proj_df$y))),
                     expand = TRUE) +
      #Remove grid lines
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      #Save plot to hard drive
      ggsave(filename = file.path(fp_out, species, "BRTs", version, "Plots", paste0("proj_", proj_year, "_", i, ".png")), width = 7, height = 7)
    
    # -------- Extract predicted values at locations of original data --------
    month_md$pred <- as.data.frame(raster::extract(x = proj, y = month_md[c("lon", "lat")]))
    
    # -------- Convert NAs to zeros --------
    month_md$pred[is.na(month_md$pred)] <- 0
    month_md$abund[is.na(month_md$abund)] <- 0
    
    # -------- Get RMSE --------
    Metrics::rmse(actual = month_md$abund, 
                  predicted = month_md$pred)
    
    # -------- Get r^2 --------
    rsq <- function (x, y) cor(x, y) ^ 2
    write.table(rsq(month_md$abund, month_md$pred), file = file.path(fp_out, species, "BRTs", version, "Evals", paste0("rsq_", i)))
    
    # -------- Unlist variables -------- 
    month_md$abund <- unlist(month_md$abund)
    month_md$pred <- unlist(month_md$pred)
    
    # -------- Plot actual vs. predicted values
    ggplot(data = month_md, aes(x = abund, y = pred)) +
      geom_point() +
      ylim(c(min(md$abund), max(md$abund))) +
      ggsave(filename = file.path(fp_out, species, "BRTs", version, "Plots", paste0("actualvspred_", proj_year, "_", i, ".png")))
    
  }
}






