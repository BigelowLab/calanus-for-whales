# Camille Ross
# 10 June, 2020
# Purpose: Build monthly ensembles of GAMs and BRTs

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
# Source data binding function
source("./calanus_data/Code/bind_years.R")

# -------- Main function --------
#'@param version <chr> version of model
#'@param fp_md <chr> file path to formatted data used in model
#'@param species <chr> species to model; choices are "cfin", "ctyp", or "pseudo"
#'@param fp_covars <chr> file path to environmental covariate data
#'@param years <vectors> years for which to run the model
#'@param fp_out <chr> file path save the data to 
build_ensemble <- function(version, fp_md, species, fp_covars, years, proj_year, fp_out) {
  
  # -------- Create output directories --------
  dir.create(fp_out, showWarnings = FALSE) 
  dir.create(file.path(fp_out, species, version), showWarnings = FALSE)
  dir.create(file.path(fp_out, species, version, "Ensembles"), showWarnings = FALSE)
  # -------- Create directory for projections --------
  dir.create(file.path(fp_out, species, version, "Ensembles", "Projections"), showWarnings = FALSE)
  # -------- Create directory for plots --------
  dir.create(file.path(fp_out, species, version, "Ensembles", "Plots"), showWarnings = FALSE)
  # -------- Create directory for evaluations --------
  dir.create(file.path(fp_out, species, version, "Ensembles", "Evals"), showWarnings = FALSE)
  
  # -------- Load model data --------
  if (length(years) == 1) {
    md <- readr::read_csv(file.path(fp_md, paste0(years[1], ".csv"))) %>% 
      dplyr::filter(dataset == "CPR" | dataset == "NOAA_CPR")
  } else {
    md <- bind_years(fp = file.path(fp_md), years = years) %>% 
      dplyr::filter(dataset == "CPR" | dataset == "NOAA_CPR")
  }
  
  # -------- Take the log of count data --------
  md$abund <- as.data.frame(log10(md[paste0(species, "_CV_VI")] + 1))$cfin_CV_VI
  
  # -------- Exclude NAs and select columns --------
  md <- md %>% dplyr::select(lat, lon, month, abund, wind, fetch, chl, int_chl, bots, bott, sss, sst, uv, bat, dist, slope) %>%
    as.data.frame() %>%
    na.exclude() %>%
    dplyr::mutate(season = if_else(month %in% c(1:3), 1,
                                   if_else(month %in% c(4:6), 2,
                                           if_else(month %in% c(7:9), 3, 4))))
  
  # -------- Take log of chlorophyll --------
  md$chl <- log(abs(md$chl))
  md$int_chl <- log(abs(md$int_chl))
  
  # -------- Load world map data --------
  worldmap <- ggplot2::map_data("world")
  
  # -------- Initialize projection year --------
  proj_year <- 2010
  
  # -------- Loop over months --------
  for (i in 1:12) {
    
    # -------- Test if projection exists --------
    if (paste0("proj_", proj_year, "_", i, ".grd") %in% list.files(file.path(fp_out, species, version, "GAMs", "Projections")) &
        paste0("proj_", proj_year, "_", i, ".grd") %in% list.files(file.path(fp_out, species, version, "BRTs", "Projections"))) {
      
      # -------- Isolate month data --------
      month_md <- md %>% dplyr::filter(month == i)
      
      # -------- Load projections --------
      brt_proj <- raster::raster(file.path(fp_out, species, version, "BRTs", "Projections", paste0("proj_", proj_year, "_", i, ".grd")))
      gam_proj <- raster::raster(file.path(fp_out, species, version, "GAMs", "Projections", paste0("proj_", proj_year, "_", i, ".grd")))
      
      ensemble <- mean(brt_proj, gam_proj)
      
      raster::writeRaster(x = ensemble, filename = file.path(fp_out, species, version, "Ensembles", "Projections", paste0("proj_", proj_year, "_", i, ".grd")),
                          overwrite = TRUE)
      
      # Convert proj to data frame
      proj_df <- as.data.frame(ensemble, xy = TRUE)
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
        ggsave(filename = file.path(fp_out, species, version, "Ensembles", "Plots", paste0("proj_", proj_year, "_", i, ".png")), width = 7, height = 7)
      
      # -------- Extract predicted values at locations of original data --------
      month_md$pred <- as.data.frame(raster::extract(x = ensemble, y = month_md[c("lon", "lat")]))
      
      # -------- Convert NAs to zeros --------
      month_md$pred[is.na(month_md$pred)] <- 0
      month_md$abund[is.na(month_md$abund)] <- 0
      
      # -------- Get RMSE --------
      write.table(Metrics::rmse(actual = month_md$abund, 
                                predicted = month_md$pred), file = file.path(fp_out, species, version, "Ensembles", "Evals", paste0("rmse_", i, ".csv")))
      
      
      # -------- Get r^2 --------
      rsq <- function (x, y) cor(x, y) ^ 2
      write.table(rsq(month_md$abund, month_md$pred), file = file.path(fp_out, species, version, "Ensembles", "Evals", paste0("rsq_", i, ".csv")))
      
      
      # -------- Unlist variables -------- 
      month_md$abund <- unlist(month_md$abund)
      month_md$pred <- unlist(month_md$pred)
      
      # -------- Plot actual vs. predicted values
      ggplot(data = month_md, aes(x = abund, y = pred)) +
        geom_point() +
        ylim(c(min(md$abund), max(md$abund))) +
        ggsave(filename = file.path(fp_out, species, version, "Ensembles", "Plots", paste0("actualvspred_", proj_year, "_", i, ".png")))
    }
  }
}







  