require(viridis)
require(Metrics)
require(mapdata)
require(maps)
require(ggplot2)
require(stats)
require(pals)
require(dplyr)
require(raster)

# -------- Main function --------
#'@param version <chr> version of model
#'@param fp_out <chr> file path save the data to 
#'@param years <vectors> years for which to run the model
#'@param species <chr> species to model; choices are "cfin", "ctyp", or "pcal"
build_biomod_climatology <- function(version, fp_out, years, species = "cfin") {

  # -------- Create output directories --------
  dir.create(fp_out, showWarnings = FALSE) 
  dir.create(file.path(fp_out, species, version, "Biomod_Climatologies"), showWarnings = FALSE)
  # -------- Create directory for projections --------
  dir.create(file.path(fp_out, species, version, "Biomod_Climatologies", "Projections"), showWarnings = FALSE)
  # -------- Create directory for plots --------
  dir.create(file.path(fp_out, species, version, "Biomod_Climatologies", "Plots"), showWarnings = FALSE)
  
  # -------- Load world map data --------
  worldmap <- ggplot2::map_data("world")
  
  months <- c("January", "February", "March",
              "April", "May", "June", 
              "July", "August", "September",
              "October", "November", "December")
  
  # -------- Loop over months --------
  for (i in 1:12) {
    
    # -------- Test if projection exists --------
    for (year in years) {
      if (paste0("proj_Ensemble_", year, "_", i, "_", species, version, "_ensemble.grd") %in% list.files(file.path(fp_out, species, version, "Biomod", "Projections", paste0(species, version), paste0("proj_Ensemble_", year, "_", i))) &
          paste0("proj_GAM_", year, "_", i, "_", species, version, ".grd") %in% list.files(file.path(fp_out, species, version, "Biomod", "Projections", paste0(species, version), paste0("proj_GAM_", year, "_", i)))) {
        if (year == 2000 & i == 1) {
          # Load ensemble as raster
          # Divide by 1000 to convert probabilities to percentages
          ensemble_proj <- raster(file.path(fp_out, species, version, "Biomod", "Projections", paste0(species, version), paste0("proj_Ensemble_", year, "_", i), paste0("proj_Ensemble_", year, "_", i, "_", species, version, "_ensemble.grd"))) %>%
            `/`(1000)
          
          gam_proj <- raster(file.path(fp_out, species, version, "Biomod", "Projections", paste0(species, version), paste0("proj_GAM_", year, "_", i), paste0("proj_GAM_", year, "_", i, "_", species, version, ".grd"))) %>%
            `/`(1000)
          
          brt_proj <- raster(file.path(fp_out, species, version, "Biomod", "Projections", paste0(species, version), paste0("proj_GBM_", year, "_", i), paste0("proj_GBM_", year, "_", i, "_", species, version, ".grd"))) %>%
            `/`(1000)
          
          rf_proj <- raster(file.path(fp_out, species, version, "Biomod", "Projections", paste0(species, version), paste0("proj_RF_", year, "_", i), paste0("proj_RF_", year, "_", i, "_", species, version, ".grd"))) %>%
            `/`(1000)
          
          crs(ensemble_proj) <- '+init=epsg:4121 +proj=longlat +ellps=GRS80 +datum=GGRS87 +no_defs +towgs84=-199.87,74.79,246.62'
          crs(gam_proj) <- '+init=epsg:4121 +proj=longlat +ellps=GRS80 +datum=GGRS87 +no_defs +towgs84=-199.87,74.79,246.62'
          crs(brt_proj) <- '+init=epsg:4121 +proj=longlat +ellps=GRS80 +datum=GGRS87 +no_defs +towgs84=-199.87,74.79,246.62'
          crs(rf_proj) <- '+init=epsg:4121 +proj=longlat +ellps=GRS80 +datum=GGRS87 +no_defs +towgs84=-199.87,74.79,246.62'
          
        } else {
          # Load ensemble as raster
          # Divide by 1000 to convert probabilities to percentages
          ensemble_temp <- raster(file.path(fp_out, species, version, "Biomod", "Projections", paste0(species, version), paste0("proj_Ensemble_", year, "_", i), paste0("proj_Ensemble_", year, "_", i, "_", species, version, "_ensemble.grd"))) %>%
            `/`(1000)
          
          gam_temp <- raster(file.path(fp_out, species, version, "Biomod", "Projections", paste0(species, version), paste0("proj_GAM_", year, "_", i), paste0("proj_GAM_", year, "_", i, "_", species, version, ".grd"))) %>%
            `/`(1000)
          
          brt_temp <- raster(file.path(fp_out, species, version, "Biomod", "Projections", paste0(species, version), paste0("proj_GBM_", year, "_", i), paste0("proj_GBM_", year, "_", i, "_", species, version, ".grd"))) %>%
            `/`(1000)
          
          rf_temp <- raster(file.path(fp_out, species, version, "Biomod", "Projections", paste0(species, version), paste0("proj_RF_", year, "_", i), paste0("proj_RF_", year, "_", i, "_", species, version, ".grd"))) %>%
            `/`(1000)
          
          crs(ensemble_temp) <- '+init=epsg:4121 +proj=longlat +ellps=GRS80 +datum=GGRS87 +no_defs +towgs84=-199.87,74.79,246.62'
          crs(gam_temp) <- '+init=epsg:4121 +proj=longlat +ellps=GRS80 +datum=GGRS87 +no_defs +towgs84=-199.87,74.79,246.62'
          crs(brt_temp) <- '+init=epsg:4121 +proj=longlat +ellps=GRS80 +datum=GGRS87 +no_defs +towgs84=-199.87,74.79,246.62'
          crs(rf_temp) <- '+init=epsg:4121 +proj=longlat +ellps=GRS80 +datum=GGRS87 +no_defs +towgs84=-199.87,74.79,246.62'
          
          # Stack rasters
          ensemble_proj <- raster::stack(ensemble_proj, ensemble_temp)
          gam_proj <- raster::stack(gam_proj, gam_temp)
          brt_proj <- raster::stack(brt_proj, brt_temp)
          rf_proj <- raster::stack(rf_proj, rf_temp)
        }
      }
    }
    
    ensemble_clim <- calc(ensemble_proj, fun = mean, na.rm = TRUE)
    gam_clim <- calc(gam_proj, fun = mean, na.rm = TRUE)
    brt_clim <- calc(brt_proj, fun = mean, na.rm = TRUE)
    rf_clim <- calc(rf_proj, fun = mean, na.rm = TRUE)
    
    raster::writeRaster(x = ensemble_clim, filename = file.path(fp_out, species, version, "Biomod_Climatologies", "Projections", paste0("ensemble_proj_", i, ".grd")),
                        overwrite = TRUE)
    raster::writeRaster(x = gam_clim, filename = file.path(fp_out, species, version, "Biomod_Climatologies", "Projections", paste0("gam_proj_", i, ".grd")),
                        overwrite = TRUE)
    raster::writeRaster(x = brt_clim, filename = file.path(fp_out, species, version, "Biomod_Climatologies", "Projections", paste0("brt_proj_", i, ".grd")),
                        overwrite = TRUE)
    raster::writeRaster(x = rf_clim, filename = file.path(fp_out, species, version, "Biomod_Climatologies", "Projections", paste0("rf_proj_", i, ".grd")),
                        overwrite = TRUE)
    
    # Convert proj to data frame
    # Save projection raster as data frame with xy coords and no NAs
    ensemble_df <- as.data.frame(ensemble_clim, xy = TRUE, na.rm = TRUE)
    gam_df <- as.data.frame(gam_clim, xy = TRUE, na.rm = TRUE)
    brt_df <- as.data.frame(brt_clim, xy = TRUE, na.rm = TRUE)
    rf_df <- as.data.frame(rf_clim, xy = TRUE, na.rm = TRUE)
    
    # Assign column names
    names(ensemble_df) <- c('x', 'y', "pred")
    names(gam_df) <- c('x', 'y', "pred")
    names(brt_df) <- c('x', 'y', "pred")
    names(rf_df) <- c('x', 'y', "pred")

    # -------- Plot projection --------
    # ---- Ensemble ----
    ggplot() + 
      # Add projection data
      geom_tile(data = ensemble_df, aes(x, y, fill = pred)) +
      # Add projection color gradient and label
      scale_fill_gradientn(colors = inferno(500), limits = c(0, 0.5), na.value = "white") +
      labs(x = "", 
           y = "",
           title = months[i]) +
      # Add world map data
      geom_polygon(data = worldmap, aes(long, lat, group = group), fill = NA, colour = "gray43") +
      coord_quickmap(xlim = c(round(min(gam_df$x)), round(max(gam_df$x))), 
                     ylim = c(round(min(gam_df$y)), round(max(gam_df$y))),
                     expand = TRUE) +
      # Remove grid lines
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            plot.title = element_text(size=50)) +
      # Save plot to hard drive
      ggsave(filename = file.path(fp_out, species, version, "Biomod_Climatologies", "Plots", paste0("ensemble_proj_", i, ".png")), width = 7, height = 7)
    
    # ---- GAM ----
    ggplot() + 
      # Add projection data
      geom_tile(data = gam_df, aes(x, y, fill = pred)) +
      # Add projection color gradient and label
      scale_fill_gradientn(colors = inferno(500), limits = c(0, 0.5), na.value = "white") +
      labs(x = "", 
           y = "",
           title = months[i]) +
      # Add world map data
      geom_polygon(data = worldmap, aes(long, lat, group = group), fill = NA, colour = "gray43") +
      coord_quickmap(xlim = c(round(min(gam_df$x)), round(max(gam_df$x))), 
                     ylim = c(round(min(gam_df$y)), round(max(gam_df$y))),
                     expand = TRUE) +
      # Remove grid lines
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            plot.title = element_text(size=50)) +
      # Save plot to hard drive
      ggsave(filename = file.path(fp_out, species, version, "Biomod_Climatologies", "Plots", paste0("gam_proj_", i, ".png")), width = 7, height = 7)
    
    # ---- BRT ----
    ggplot() + 
      # Add projection data
      geom_tile(data = brt_df, aes(x, y, fill = pred)) +
      # Add projection color gradient and label
      scale_fill_gradientn(colors = inferno(500), limits = c(0, 0.5), na.value = "white") +
      labs(x = "", 
           y = "",
           title = months[i]) +
      # Add world map data
      geom_polygon(data = worldmap, aes(long, lat, group = group), fill = NA, colour = "gray43") +
      coord_quickmap(xlim = c(round(min(gam_df$x)), round(max(gam_df$x))), 
                     ylim = c(round(min(gam_df$y)), round(max(gam_df$y))),
                     expand = TRUE) +
      # Remove grid lines
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            plot.title = element_text(size=50)) +
      # Save plot to hard drive
      ggsave(filename = file.path(fp_out, species, version, "Biomod_Climatologies", "Plots", paste0("brt_proj_", i, ".png")), width = 7, height = 7)
    
    # ---- RF ----
    ggplot() + 
      # Add projection data
      geom_tile(data = rf_df, aes(x, y, fill = pred)) +
      # Add projection color gradient and label
      scale_fill_gradientn(colors = inferno(500), limits = c(0, 0.5), na.value = "white") +
      labs(x = "", 
           y = "",
           title = months[i]) +
      # Add world map data
      geom_polygon(data = worldmap, aes(long, lat, group = group), fill = NA, colour = "gray43") +
      coord_quickmap(xlim = c(round(min(gam_df$x)), round(max(gam_df$x))), 
                     ylim = c(round(min(gam_df$y)), round(max(gam_df$y))),
                     expand = TRUE) +
      # Remove grid lines
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            plot.title = element_text(size=50)) +
      # Save plot to hard drive
      ggsave(filename = file.path(fp_out, species, version, "Biomod_Climatologies", "Plots", paste0("rf_proj_", i, ".png")), width = 7, height = 7)
    
  }
}
