# Camille Ross
# 10 June, 2020
# Purpose: Build monthly climatologies of GAMs and BRTs

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
require(pals)

# -------- Main function --------
#'@param version <chr> version of model
#'@param fp_out <chr> file path save the data to 
#'@param years <vectors> years for which to run the model
#'@param species <chr> species to model; choices are "cfin", "ctyp", or "pcal"
build_climatology <- function(version, fp_out, years, species = "cfin", anomaly = FALSE) {
  
  # -------- Create output directories --------
  dir.create(fp_out, showWarnings = FALSE) 
  dir.create(file.path(fp_out, species, version, "Climatologies"), showWarnings = FALSE)
  # -------- Create directory for projections --------
  dir.create(file.path(fp_out, species, version, "Climatologies", "Projections"), showWarnings = FALSE)
  # -------- Create directory for plots --------
  dir.create(file.path(fp_out, species, version, "Climatologies", "Plots"), showWarnings = FALSE)
  
  # -------- Load world map data --------
  worldmap <- ggplot2::map_data("world")
  
  # -------- Load bathymetry layer ---------
  bat <- load_covars(fp_covars = fp_covars, year = 2000, month = 1,
                     env_covars = "bat",
                     as_raster = TRUE)
  
  # -------- Loop over months --------
  for (i in 1:12) {
    print(i)
    
    # -------- Test if projection exists --------
    for (year in years) {
      if (paste0("proj_", year, "_", i, ".tif") %in% list.files(file.path(fp_out, species, version, "GAMs", "Projections")) &
          paste0("proj_", year, "_", i, ".tif") %in% list.files(file.path(fp_out, species, version, "BRTs", "Projections"))) {
        if (year == 2000 & i == 1) {
          gam_proj <- raster::raster(file.path(fp_out, species, version, "GAMs", "Projections", paste0("proj_", year, "_", i, ".tif")))
          brt_proj <- raster::raster(file.path(fp_out, species, version, "BRTs", "Projections", paste0("proj_", year, "_", i, ".tif")))
        } else {
          gam_temp <- raster::raster(file.path(fp_out, species, version, "GAMs", "Projections", paste0("proj_", year, "_", i, ".tif")))
          gam_proj <- raster::stack(gam_proj, gam_temp)
          
          brt_temp <- raster::raster(file.path(fp_out, species, version, "BRTs", "Projections", paste0("proj_", year, "_", i, ".tif")))
          brt_proj <- raster::stack(brt_proj, brt_temp)
        }
      }
    }
    
    gam_clim <- calc(gam_proj, fun = mean, na.rm = TRUE)
    brt_clim <- calc(brt_proj, fun = mean, na.rm = TRUE)
      
    raster::writeRaster(x = gam_clim, filename = file.path(fp_out, species, version, "Climatologies", "Projections", paste0("gam_proj_", i, ".tif")),
                        overwrite = TRUE)
    raster::writeRaster(x = brt_clim, filename = file.path(fp_out, species, version, "Climatologies", "Projections", paste0("brt_proj_", i, ".tif")),
                        overwrite = TRUE)
    
    # Zero out projections below 1000m
    if (species == "cfin") {
      gam_clim[bat >= 1000] <- 0
      brt_clim[bat >= 1000] <- 0
    }
    
    # Convert proj to data frame
    gam_df <- as.data.frame(gam_clim, xy = TRUE)
    brt_df <- as.data.frame(brt_clim, xy = TRUE)
    # Add column names
    names(gam_df) <- c("x", "y", "pred")
    names(brt_df) <- c("x", "y", "pred")
    
    if (anomaly) {
      color_scale <- ocean.balance(500)
      max_val <- 2
      min_val <- -2
    } else {
      color_scale <- inferno(500)
      max_val <- 9
      min_val <- -2
      # Zero out negative values
      gam_clim[gam_clim < 0] <- 0
      brt_clim[brt_clim < 0] <- 0
    }
    
    # -------- Plot projection --------
    # ---- GAMs ----
    ggplot() + 
      #Add projection data
      geom_tile(data = gam_df, aes(x, y, fill = pred)) +
      #Add projection color gradient and label
      scale_fill_gradientn(colors = color_scale, limits = c(min_val, max_val), na.value = "white") +
      labs(x = "", 
           y = "") +
      # Add world map data
      geom_polygon(data = worldmap, aes(long, lat, group = group), fill = NA, colour = "gray43") +
      coord_quickmap(xlim = c(round(min(gam_df$x)), round(max(gam_df$x))), 
                     ylim = c(round(min(gam_df$y)), round(max(gam_df$y))),
                     expand = TRUE) +
      #Remove grid lines
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      #Save plot to hard drive
      ggsave(filename = file.path(fp_out, species, version, "Climatologies", "Plots", paste0("gam_proj_", i, ".png")), width = 7, height = 7)
    
    # ---- BRTs ----
    ggplot() + 
      #Add projection data
      geom_tile(data = brt_df, aes(x, y, fill = pred)) +
      #Add projection color gradient and label
      scale_fill_gradientn(colors = color_scale, limits = c(min_val, max_val), na.value = "white") +
      labs(x = "", 
           y = "") +
      # Add world map data
      geom_polygon(data = worldmap, aes(long, lat, group = group), fill = NA, colour = "gray43") +
      coord_quickmap(xlim = c(round(min(brt_df$x)), round(max(brt_df$x))), 
                     ylim = c(round(min(brt_df$y)), round(max(brt_df$y))),
                     expand = TRUE) +
      #Remove grid lines
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      #Save plot to hard drive
      ggsave(filename = file.path(fp_out, species, version, "Climatologies", "Plots", paste0("brt_proj_", i, ".png")), width = 7, height = 7)
  }
}


