
library(biomod2)

project_biomod <- function(model, modelOut, j, species, version, fp_out) {
  
  months <- c("January", "February", "March",
              "April", "May", "June",
              "July", "August", "September",
              "October", "November", "December")
  # ------- Project GAMS -------
  # Create vector all runs of model algorithm for projection
  select_models <- c()
  for (k in 1:10) {
    select_models[k] <- paste0(species, version, "_AllData_RUN", k, "_", model)
  }
  
  proj <- BIOMOD_Projection(modeling.output = modelOut,
                            new.env = covars,
                            proj.name = paste0(model, "_", j),
                            selected.models = select_models,
                            binary.meth = 'ROC',
                            compress = 'xz',
                            build.clamping.mask = TRUE,
                            output.format = '.grd')
  
  
  # Load ensemble forecast as raster
  # Divide by 1000 to convert probabilities to percentages
  proj_raster <- raster(file.path(paste0(species, version), paste0("proj_", model, "_", j), paste0("proj_", model, "_", j, "_", species, version, '.grd'))) %>%
    `/`(1000)
  
  crs(proj_raster) <- '+init=epsg:4121 +proj=longlat +ellps=GRS80 +datum=GGRS87 +no_defs +towgs84=-199.87,74.79,246.62'
  
  # Save projection raster as data frame with xy coords and no NAs
  proj_df <- as.data.frame(proj_raster, xy = TRUE, na.rm = TRUE)
  
  # Assign column names
  names(proj_df) <- c('pred', 'x', 'y')
  
  # -------- Plot projection --------
  ggplot() +
    # Add projection data
    geom_tile(data = proj_df, aes(x, y, fill = pred)) +
    # Add projection color gradient and label
    scale_fill_gradientn(colors = inferno(500), limits = c(0,1), na.value = "white") +
    labs(x = "",
         y = "",
         main = months[j]) +
    # Add world map data
    geom_polygon(data = worldmap, aes(long, lat, group = group), fill = NA, colour = "gray43") +
    coord_quickmap(xlim = c(round(min(proj_df$x)), round(max(proj_df$x))),
                   ylim = c(round(min(proj_df$y)), round(max(proj_df$y))),
                   expand = TRUE) +
    # Remove grid lines
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    # Save plot to hard drive
    ggsave(filename = file.path(fp_out, species, version, "Biomod", "Plots", paste0(model, "_proj_", j, ".png")), width = 7, height = 7)
  
}