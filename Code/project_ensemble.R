library(biomod2)

project_ensemble <- function(modelOut, biomodEM, j, fp_out, species, version) {
  
  # -------- Project ensemble model --------
  biomodProjFuture <- BIOMOD_Projection(modeling.output = modelOut,
                                        new.env = covars,
                                        proj.name = 'all',
                                        selected.models = 'all',
                                        binary.meth = 'ROC',
                                        compress = 'xz',
                                        build.clamping.mask = FALSE,
                                        output.format = '.grd')
  
  # Build ensemble forecast
  myBiomodEF <- BIOMOD_EnsembleForecasting(EM.output = biomodEM,
                                           projection.output = biomodProjFuture,
                                           proj.name = paste0("Ensemble_", j))
  
  # Load ensemble forecast as raster
  # Divide by 1000 to convert probabilities to percentages
  ensemble_proj_raster <- raster(file.path(paste0(species, version), paste0("proj_Ensemble_", j), paste0("proj_Ensemble_", j, "_", species, version, "_ensemble.grd"))) %>%
    `/`(1000)
  
  crs(ensemble_proj_raster) <- '+init=epsg:4121 +proj=longlat +ellps=GRS80 +datum=GGRS87 +no_defs +towgs84=-199.87,74.79,246.62'
  
  # Save projection raster as data frame with xy coords and no NAs
  ensemble_proj_df <- as.data.frame(ensemble_proj_raster, xy = TRUE, na.rm = TRUE)
  
  # Assign column names
  names(ensemble_proj_df) <- c('x', 'y', 'pred')
  
  # -------- Plot projection --------
  ggplot() + 
    # Add projection data
    geom_tile(data = ensemble_proj_df, aes(x, y, fill = pred)) +
    # Add projection color gradient and label
    scale_fill_gradientn(colors = inferno(500), limits = c(0,1), na.value = "white") +
    labs(x = "", 
         y = "",
         fill = "Probability") +
    # Add world map data
    geom_polygon(data = worldmap, aes(long, lat, group = group), fill = NA, colour = "gray43") +
    coord_quickmap(xlim = c(round(min(ensemble_proj_df$x)), round(max(ensemble_proj_df$x))), 
                   ylim = c(round(min(ensemble_proj_df$y)), round(max(ensemble_proj_df$y))),
                   expand = TRUE) +
    # Remove grid lines
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    # Save plot to hard drive
    ggsave(filename = file.path(fp_out, species, version, "Biomod", "Plots", paste0("ensemble_proj_", j, ".png")), width = 7, height = 7)
  
  
  # # -------- Extract predicted values --------
  # month_md$ensemble_pred <- raster::extract(ensemble_proj_raster, month_md$abund)
  # 
  # # -------- Convert NAs to zeros --------
  # month_md$ensemble_pred[is.na(month_md$ensemble_pred)] <- 0
  # month_md$abund[is.na(month_md$abund)] <- 0
  # 
  # # -------- Unlist variables --------
  # month_md$abund <- unlist(month_md$abund)
  # month_md$ensemble_pred <- unlist(month_md$ensemble_pred)
  # 
  # readr::write_csv(month_md %>% dplyr::select(abund, ensemble_pred), file.path(fp_out, species, version, "Biomod", "Projections", paste0("ensemble_abund_vs_pred_", i, "_", j, ".csv")))
  # 
  # # -------- Plot actual vs. predicted values --------
  # ggplot(data = month_md, aes(x = log10(abund + 1), y = ensemble_pred)) +
  #   geom_point() +
  #   ylim(c(0, 1)) +
  #   ggsave(filename = file.path(fp_out, species, version, "Biomod", "Plots", paste0("ensemble_actualvspred_", j, ".png")))
  # 
}