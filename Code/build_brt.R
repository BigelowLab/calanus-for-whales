# Camille Ross
# 10 June, 2020
# Purpose: Run a boosted regression tree model for zooplankton species with variable environmental covariates

# -------- Load libraries --------
require(dplyr)
require(readr)
require(dismo)
require(gbm)
require(pals)
require(viridis)
require(Metrics)
require(mapdata)
require(maps)
require(ggplot2)
require(stats)
require(devtools)
require(rgdal)

# -------- Source outside files --------
# Source data formatting function
source("./calanus-for-whales/Code/format_model_data.R")
# Source covariate loading function
source("./calanus-for-whales/Code/load_covars.R")
# Source data binding function
source("./calanus_data/Code/bind_years.R")

# -------- Main function --------
#'@param version <chr> version of model
#'@param fp_md <chr> file path to formatted data used in model
#'@param fp_covars <chr> file path to environmental covariate data
#'@param env_covars <vector> vector of covariates to include in the model
#'@param years <vectors> years for which to run the model
#'@param fp_out <chr> file path save the data to 
#'@param species <chr> species to model; choices are "cfin", "ctyp", or "pcal"
#'@param anomaly <logical> if true, model is run using calanus anomaly
#'@param format_data <logical> if true, data is formatted within function; only used if model_data is NULL
#'@param fp_zpd <chr> filepath to the zooplankton database if data is formatted within function
build_brt <- function(version, fp_md, datasets, fp_covars, env_covars, years, fp_out, 
                      species = "cfin", anomaly = FALSE,
                      format_data = FALSE, fp_zpd = NULL) {
  
  # -------- Create output directories --------
  dir.create(fp_out, showWarnings = FALSE) 
  dir.create(file.path(fp_out, species), showWarnings = FALSE)
  dir.create(file.path(fp_out, species, version), showWarnings = FALSE)
  dir.create(file.path(fp_out, species, version, "BRTs"), showWarnings = FALSE)
  # -------- Create directory for projections --------
  dir.create(file.path(fp_out, species, version, "BRTs", "Projections"), showWarnings = FALSE)
  # -------- Create directory for plots --------
  dir.create(file.path(fp_out, species, version, "BRTs", "Plots"), showWarnings = FALSE)
  # -------- Create directory for evaluations --------
  dir.create(file.path(fp_out, species, version, "BRTs", "Evals"), showWarnings = FALSE)
  
  # -------- Format model data --------
  if (format_data) {
    # Format zooplankton and environmental covariate data
    format_model_data(fp_data = fp_zpd, fp_covars = fp_covars, fp_out = fp_md,
                      env_covars = "all", years = years)
  }
  
  # -------- Load model data --------
  if (length(years) == 1) {
    if (anomaly) {
      md <- readr::read_csv(file.path(fp_md, paste0(years[1], ".csv")))
    } else {
      md <- readr::read_csv(file.path(fp_md, paste0(years[1], ".csv"))) %>% 
        dplyr::filter(dataset %in% datasets)
    }
  } else {
    if (anomaly) {
      md <- bind_years(fp = file.path(fp_md), years = years)
    } else {
      md <- bind_years(fp = file.path(fp_md), years = years) %>%
        dplyr::filter(dataset %in% datasets)
    }
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
  if (anomaly) {
    md$abund <- md$anomaly
  } else {
    if (species == "cfin") {
      md$abund <- as.data.frame(log10(md[paste0(species, "_CV_VI")] + 1))$cfin_CV_VI
    } else if (species == "ctyp") {
      md$abund <- as.data.frame(log10(md[paste0(species, "_total")] + 1))$ctyp_total
    } else if (species == "pcal") {
      md$abund <- as.data.frame(log10(md[paste0(species, "_total")] + 1))$pcal_total
    }
  }

  # -------- Exclude NAs and select columns --------
  md <- md %>% dplyr::select(lat, lon, year, month, abund, wind, fetch, chl, int_chl, bots, bott, sss, sst, lag_sst, uv, bat, dist, slope) %>%
    as.data.frame() %>%
    na.exclude() %>%
    dplyr::filter(!is.infinite(abund)) %>%
    dplyr::mutate(season = if_else(month %in% c(1:3), 1,
                                   if_else(month %in% c(4:6), 2,
                                           if_else(month %in% c(7:9), 3, 4))))
  
  # -------- Take log of bathymetry and chlorophyll --------
  md$chl <- log(abs(md$chl))
  md$int_chl <- log(abs(md$int_chl))
  md$bat <- log(abs(md$bat))

  # -------- Initialize brt arguments --------
  brt_args <- list(n.trees = 1000, interaction.depth = 5,
                   n.minobsinnode = 3,
                   shrinkage = 0.001, bag.fraction = 0.5, 
                   train.fraction = 1, cv.folds = 10)
  
  # -------- Load world map data --------
  worldmap <- ggplot2::map_data("world")
  
  # Load bathymetry layer
  bat <- load_covars(fp_covars = fp_covars, year = 2000, month = 1,
                     env_covars = "bat",
                     as_raster = TRUE)
  
  # -------- Loop over months --------
  for (i in years) {
    for (j in 1:12) {
      # -------- Isolate month data --------
      month_md <- md %>% dplyr::filter(month == j)
      
      if ((nrow(month_md) < 15) | (length(unique(month_md$abund)) == 1)) {
        print("skipped")
        next
      }
    
      # -------- Build GAM with all covariates --------
      brt_sdm <- gbm::gbm(formula = stats::reformulate(env_covars, "abund"),
                          distribution = "gaussian",
                          data = month_md,
                          n.trees = brt_args[['n.trees']],
                          interaction.depth = brt_args[['interaction.depth']],
                          n.minobsinnode = brt_args[['n.minobsinnode']],
                          shrinkage = brt_args[['shrinkage']],
                          bag.fraction = brt_args[['bag.fraction']],
                          train.fraction = brt_args[['train.fraction']],
                          cv.folds = brt_args[['cv.folds']])
      
      # -------- Plot results --------
      plot(brt_sdm)
      # -------- Load summary of model --------
      png(file.path(fp_out, species, version, "BRTs", "Plots", paste0("var_cont_", i, "_", j, ".png")))
      summary(brt_sdm)
      dev.off()
        
      # -------- Load environmental covariates for projection --------
      covars <- load_covars(fp_covars = fp_covars, year = i, month = j,
                            env_covars = env_covars,
                            as_raster = TRUE)
      
      if ("chl" %in% env_covars) {
        covars$chl <- log(covars$chl)
      }
      if ("int_chl" %in% env_covars) {
        covars$int_chl <- log(covars$int_chl)
      }
      if ("bat" %in% env_covars) {
        covars$bat <- log(covars$bat)
      }
      
      # -------- Project model onto covariates --------
      proj <- raster::predict(covars, brt_sdm,
                              filename = file.path(fp_out, species, version, "BRTs", "Projections", paste0("proj_", i, "_", j)), progress = "text",
                              overwrite = TRUE,
                              format = "GTiff")
      
      # Zero out projections below 1000m
      proj[bat >= 1000] <- 0
      
      if (anomaly) {
        color_scale <- ocean.balance(500)
        max_val <- 7
        min_val <- -7
      } else {
        color_scale <- inferno(500)
        max_val <- max(md$abund, na.rm = TRUE) + 4
        min_val <- 0
        # Zero out negative values
        proj[proj < 0] <- 0
      }
                        
      # Convert proj to data frame
      proj_df <- as.data.frame(proj, xy = TRUE)
      # Add column names
      names(proj_df) <- c("x", "y", "pred")
      
      
      
      # -------- Plot projection --------
      ggplot() + 
        #Add projection data
        geom_tile(data = proj_df, aes(x, y, fill = pred)) +
        #Add projection color gradient and label
        scale_fill_gradientn(colors = ocean.balance(500), limits = c(min(md$abund), max(md$abund)), na.value = "white") +
        labs(x = "", 
             y = "") +
        # geom_point(data = month_md, aes(lon, lat, color = abund)) +
        # scale_color_gradientn(colors = ocean.balance(500), limits = c(min(md$abund), max(md$abund)), na.value = "white") +
        # Add world map data
        geom_polygon(data = worldmap, aes(long, lat, group = group), fill = NA, colour = "gray43") +
        coord_quickmap(xlim = c(round(min(proj_df$x)), round(max(proj_df$x))), 
                       ylim = c(round(min(proj_df$y)), round(max(proj_df$y))),
                       expand = TRUE) +
        #Remove grid lines
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        #Save plot to hard drive
        ggsave(filename = file.path(fp_out, species, version, "BRTs", "Plots", paste0("proj_", i, "_", j, ".png")), width = 7, height = 7)
      
      # -------- Extract predicted values at locations of original data --------
      month_md$pred <- as.data.frame(brt_sdm$fit)
      
      # -------- Convert NAs to zeros --------
      month_md$pred[is.na(month_md$pred)] <- 0
      month_md$abund[is.na(month_md$abund)] <- 0
      
      # -------- Get RMSE --------
      write.table(Metrics::rmse(actual = month_md$abund, 
                                predicted = as.numeric(unlist(month_md$pred))), file = file.path(fp_out, species, version, "BRTs", "Evals", paste0("rmse_", i, "_", j, ".csv")))
      
      
      # -------- Get r^2 --------
      rsq <- function (x, y) cor(x, y) ^ 2
      write.table(rsq(month_md$abund, month_md$pred), file = file.path(fp_out, species, version, "BRTs", "Evals", paste0("rsq_", i, "_", j, ".csv")))
      
      
      # -------- Unlist variables -------- 
      month_md$abund <- unlist(month_md$abund)
      month_md$pred <- unlist(month_md$pred)
      
      readr::write_csv(month_md %>% dplyr::select(abund, pred), file.path(fp_out, species, version, "BRTs", "Projections", paste0("abund_vs_pred_", i, "_", j, ".csv")))
      
      
      # -------- Plot actual vs. predicted values
      ggplot(data = month_md, aes(x = abund, y = pred)) +
        geom_point() +
        ylim(c(min(md$abund), max(md$abund))) +
        ggsave(filename = file.path(fp_out, species, version, "BRTs", "Plots", paste0("actual_vs_pred_", i, "_", j, ".png")))
    }
  }
}






