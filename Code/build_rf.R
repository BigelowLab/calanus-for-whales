# Camille Ross
# 10 June, 2020
# Purpose: Run a random forest model for zooplankton species with variable environmental covariates

# -------- Load libraries --------
require(dplyr)
require(readr)
require(dismo)
require(randomForest)
library(viridis)
library(ggplot2)

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
#'@param species <chr> species to model; choices are "cfin", "ctyp", or "pseudo"
#'@param fp_covars <chr> file path to environmental covariate data
#'@param env_covars <vector> vector of covariates to include in the model
#'@param years <vectors> years for which to run the model
#'@param fp_out <chr> file path save the data to 
#'@param format_data <logical> if true, data is formatted within function; only used if model_data is NULL
#'@param fp_zpd <chr> filepath to the zooplankton database if data is formatted within function
build_gam <- function(version, fp_md, species, fp_covars, env_covars, years, proj_year, fp_out,
                      format_data = FALSE, fp_zpd = NULL) {
  
  # -------- Create output directory --------
  dir.create(fp_out, showWarnings = FALSE) 
  dir.create(file.path(fp_out, species, version), showWarnings = FALSE)
  dir.create(file.path(fp_out, species, version, "RFs"), showWarnings = FALSE)
  # -------- Create directory for projections --------
  dir.create(file.path(fp_out, species, version, "RFs", "Projections"), showWarnings = FALSE)
  # -------- Create directory for plots ---------
  dir.create(file.path(fp_out, species, version, "RFs", "Plots"), showWarnings = FALSE)
  
  # -------- Format model data --------
  if (format_data) {
    # Format zooplankton and environmental covariate data
    format_model_data(fp_data = fp_zpd, fp_covars = fp_covars, env_covars = env_covars, years = years, fp_out = fp_md)
  }
  
  # -------- Load model data --------
  if (length(years) == 1) {
    md <- readr::read_csv(file.path(fp_md, paste0(years[1], ".csv")))
  } else {
    md <- bind_years(fp = file.path(fp_md), years = years)
  }
  # Take the log of count data
  md$abund <- log10(md[paste0(species, "_total")] + 1)
  
  # -------- Initialize gam arguments --------
  # gam_args <- list(k = -1, bs = 'cs')
  
  proj_year <- 2009
  
  # -------- Loop over months --------
  for (i in 1:12) {
    # -------- Isolate month data --------
    month_md <- md %>% dplyr::filter(month == i)
    # -------- Exclude NAs from dataset --------
    month_md <- na.exclude(month_md)
    
    # -------- Build RF with all covariates --------
    rf_sdm <- randomForest(stats::reformulate(env_covars, "abund"),
                data = month_md)
    
    # -------- Plot results --------
    plot(rf_sdm, pages = 1)
    # Load summary of model
    summary(rf_sdm)
    
    # -------- Load environmental covariates for projection --------
    covars <- load_covars(fp_covars = fp_covars, year = proj_year, month = i,
                          env_covars = env_covars,
                          as_raster = TRUE)
    
    # -------- Project model onto covariates --------
    proj <- raster::predict(covars, rf_sdm, filename = file.path(fp_out, species, version, "RFs", "Projections", paste0("proj_", proj_year, "_", i)), progress = "text",
                            overwrite = TRUE,
                            format = "GTif")
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
      # Add world map data
      geom_polygon(data = worldmap, aes(long, lat, group = group), fill = NA, colour = "gray43") +
      coord_quickmap(xlim = c(round(min(proj_df$x)), round(max(proj_df$x))), 
                     ylim = c(round(min(proj_df$y)), round(max(proj_df$y))),
                     expand = TRUE) +
      #Remove grid lines
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      #Save plot to hard drive
      ggsave(filename = file.path(fp_out, species, version, "RFs", "Plots", paste0("proj_", proj_year, "_", i, ".png")), width = 7, height = 7)
    
    # -------- Extract predicted values at locations of original data --------
    month_md$pred <- as.data.frame(raster::extract(x = proj, y = month_md[c("lon", "lat")])) 
    
    # -------- Convert NAs to zeros --------
    month_md$pred[is.na(month_md$pred)] <- 0
    
    # -------- Get RMSE --------
    Metrics::rmse(actual = month_md$abund, 
                  predicted = month_md$pred)
    
    # -------- Get r^2 --------
    rsq <- function (x, y) cor(x, y) ^ 2
    rsq(month_md$abund, month_md$pred)
    
    # -------- Unlist variables --------
    month_md$abund <- unlist(month_md$abund)
    month_md$pred <- unlist(month_md$pred)
    
    # -------- Plot actual vs. predicted values --------
    ggplot(data = month_md, aes(x = abund, y = pred)) +
      geom_point() +
      ylim(c(min(md$abund), max(md$abund))) +
      ggsave(filename = file.path(fp_out, species, version, "RFs", "Plots", paste0("actualvspred_", proj_year, "_", i, ".png")))
    
  }
}






