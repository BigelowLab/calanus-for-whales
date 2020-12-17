# Camille Ross
# 10 November 2020
# Purpose: Create table of zooplankton abundance per month per year

# Camille Ross
# 10 August, 2020
# Purpose: Divide model results into three regions and plot against observed abundance data.

# -------- Load libraries --------

require(lubridate)
require(readr)
require(dplyr)
require(raster)
require(ggplot2)
require(gridExtra)

# Source data binding function
source("./calanus_data/Code/bind_years.R")

# -------- Main function --------
#'@param version <chr> version of model
#'@param fp_out <chr> file path save the data to 
#'@param threshold <numeric> right whale feeding threshold
#'@param biomod_dataset <chr> dataset used to create biomod2 models
#'@param species <chr> species to model; choices are "cfin", "ctyp", or "pseudo"
plot_regions_biomod <- function(version, fp_out, threshold, biomod_dataset, species = "cfin") {
  
  # -------- Create output directories --------
  dir.create(fp_out, showWarnings = FALSE) 
  dir.create(file.path(fp_out, species, version), showWarnings = FALSE)
  dir.create(file.path(fp_out, species, version, "Climatologies"), showWarnings = FALSE)
  # -------- Create directory for evaluations --------
  dir.create(file.path(fp_out, species, version, "Climatologies", "Plots"), showWarnings = FALSE)
  
  # -------- CLIMATOLOGICAL AVERAGE --------
  
  # -------- Load model data --------
  if (length(years) == 1) {
    if (anomaly) {
      md <- readr::read_csv(file.path(fp_md, paste0(years[1], ".csv")))
    } else {
      md <- readr::read_csv(file.path(fp_md, paste0(years[1], ".csv"))) %>% 
        dplyr::filter(dataset %in% c(biomod_dataset, "CPR"))
    }
  } else {
    if (anomaly) {
      md <- bind_years(fp = file.path(fp_md), years = years)
    } else {
      md <- bind_years(fp = file.path(fp_md), years = years) %>%
        dplyr::filter(dataset %in% c(biomod_dataset, "CPR"))
    }
  }
  
  # -------- Compute anomaly --------
  if (species == "cfin") {
    md <- md %>% dplyr::group_by(dataset) %>%
      dplyr::mutate(mean = mean(log10(`cfin_CV_VI` + 1), na.rm = TRUE),
                    sd = sd(log10(`cfin_CV_VI` + 1), na.rm = TRUE),
                    var = var(log10(`cfin_CV_VI` + 1), na.rm = TRUE),
                    anomaly = (log10(`cfin_CV_VI` + 1) - mean) / sd) %>%
      dplyr::ungroup()
  } else if (species == "ctyp") {
    md <- md %>% dplyr::group_by(dataset) %>%
      dplyr::mutate(mean = mean(log10(`ctyp_total` + 1), na.rm = TRUE),
                    sd = sd(log10(`ctyp_total` + 1), na.rm = TRUE),
                    var = var(log10(`ctyp_total` + 1), na.rm = TRUE),
                    anomaly = (log10(`ctyp_total` + 1) - mean) / sd) %>%
      dplyr::ungroup()
  } else if (species == "pseudo") {
    md <- md %>% dplyr::group_by(dataset) %>%
      dplyr::mutate(mean = mean(log10(`pseudo_total` + 1), na.rm = TRUE),
                    sd = sd(log10(`pseudo_total` + 1), na.rm = TRUE),
                    var = var(log10(`pseudo_total` + 1), na.rm = TRUE),
                    anomaly = (log10(`pseudo_total` + 1) - mean) / sd) %>%
      dplyr::ungroup()
  }
  
  # -------- Take the log of count data --------
  if (anomaly) {
    md$abund <- md$anomaly
  } else {
    if (species == "cfin") {
      md$abund <- as.data.frame(md[paste0(species, "_CV_VI")])$cfin_CV_VI
    } else if (species == "ctyp") {
      md$abund <- as.data.frame(md[paste0(species, "_total")])$ctyp_total
    } else if (species == "pseudo") {
      md$abund <- as.data.frame(md[paste0(species, "_total")])$pseudo_total
    }
  }
  
  # -------- Exclude NAs and select columns --------
  md <- md %>% dplyr::select(lat, lon, year, month, dataset, abund, wind, fetch, chl, int_chl, bots, bott, sss, sst, lag_sst, uv, bat, dist, slope) %>%
    as.data.frame() %>%
    na.exclude() %>%
    dplyr::filter(!is.infinite(abund)) %>%
    dplyr::mutate(season = if_else(month %in% c(1:3), 1,
                                   if_else(month %in% c(4:6), 2,
                                           if_else(month %in% c(7:9), 3, 4))),
                  region = if_else(lat <= 41.5 & lon < -70, "MAB", 
                                   if_else(lat >= 39 & lat <= 42 & lon >= -70 & lon <= -68, "GBK", "GOM"))) %>%
    dplyr::mutate(pa = if_else(abund >= threshold, 1, 0),
                  abund = log10(abund + 1)) %>%
    dplyr::group_by(region, month, dataset) %>%
    dplyr::summarize(mean = mean(abund, na.rm=TRUE),
                     stdev = sd(abund, na.rm = TRUE),
                     var = var(abund, na.rm = TRUE),
                     var_pa = var(pa, na.rm = TRUE),
                     pa = mean(pa, na.rm = TRUE))
  
  
  for (year in 2000:2017) {
    for (i in 1:12) {
  
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
          
          ensemble_proj <- as.data.frame(ensemble_proj, xy = TRUE) %>%
            dplyr::mutate(month = i, year = year)
          gam_proj <- as.data.frame(gam_proj, xy = TRUE) %>%
            dplyr::mutate(month = i, year = year)
          brt_proj <- as.data.frame(brt_proj, xy = TRUE) %>%
            dplyr::mutate(month = i, year = year)
          rf_proj <- as.data.frame(rf_proj, xy = TRUE) %>%
            dplyr::mutate(month = i, year = year)
        
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
          
          ensemble_temp <- as.data.frame(ensemble_temp, xy = TRUE) %>%
            dplyr::mutate(month = i, year = year)
          gam_temp <- as.data.frame(gam_temp, xy = TRUE) %>%
            dplyr::mutate(month = i, year = year)
          brt_temp <- as.data.frame(brt_temp, xy = TRUE) %>%
            dplyr::mutate(month = i, year = year)
          rf_temp <- as.data.frame(rf_temp, xy = TRUE) %>%
            dplyr::mutate(month = i, year = year)
          
          # Stack rasters
          ensemble_proj <- rbind(ensemble_proj, ensemble_temp)
          gam_proj <- rbind(gam_proj, gam_temp)
          brt_proj <- rbind(brt_proj, brt_temp)
          rf_proj <- rbind(rf_proj, rf_temp)
        }
      }
    }
  }
  
  names(ensemble_proj) <- c("x", "y", "proj", "month", "year")
  names(gam_proj) <- c("x", "y", "proj", "month", "year")
  names(brt_proj) <- c("x", "y", "proj", "month", "year")
  names(rf_proj) <- c("x", "y", "proj", "month", "year")

  # -------- MONTHLY VARIABILITY --------
  # -------- Compute regions --------
  ensemble_proj_monthly <- ensemble_proj %>%
    dplyr::mutate(region = if_else(y <= 41 & x < -70, "MAB", 
                                   if_else(y >= 40 & x <= 42 & y >= -70 & x <= -68, "GBK", "GOM"))) %>%
    dplyr::group_by(region, month) %>%
    dplyr::summarize(mean = mean(proj, na.rm = TRUE),
                     stdev = sd(proj, na.rm = TRUE),
                     var = var(proj, na.rm = TRUE))
  
  gam_proj_monthly <- gam_proj %>%
    dplyr::mutate(region = if_else(y <= 41 & x < -70, "MAB", 
                                   if_else(y >= 40 & x <= 42 & y >= -70 & x <= -68, "GBK", "GOM"))) %>%
    dplyr::group_by(region, month) %>%
    dplyr::summarize(mean = mean(proj, na.rm = TRUE),
                     stdev = sd(proj, na.rm = TRUE),
                     var = var(proj, na.rm = TRUE))
  
  brt_proj_monthly <- brt_proj %>%
    dplyr::mutate(region = if_else(y <= 41 & x < -70, "MAB", 
                                   if_else(y >= 40 & x <= 42 & y >= -70 & x <= -68, "GBK", "GOM"))) %>%
    dplyr::group_by(region, month) %>%
    dplyr::summarize(mean = mean(proj, na.rm = TRUE),
                     stdev = sd(proj, na.rm = TRUE),
                     var = var(proj, na.rm = TRUE))
  
  rf_proj_monthly <- rf_proj %>%
    dplyr::mutate(region = if_else(y <= 41 & x < -70, "MAB", 
                                   if_else(y >= 40 & x <= 42 & y >= -70 & x <= -68, "GBK", "GOM"))) %>%
    dplyr::group_by(region, month) %>%
    dplyr::summarize(mean = mean(proj, na.rm = TRUE),
                     stdev = sd(proj, na.rm = TRUE),
                     var = var(proj, na.rm = TRUE))
  
  # -------- Initialize legend colors --------
  colors <- c("Actual" = "red", "Predicted" = "blue")
  
  # -------- Plot MAB --------
  if ("MAB" %in% unique(md$region)) {
    
    # ---- Plot Ensembles ----
    abund <- ggplot(data = md %>% dplyr::filter(region == "MAB", dataset %in% biomod_dataset), mapping = aes(x = month, y = mean, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Climatological Abundance",
           color = "Legend") +
      ggtitle("Ensemble Mid Atlantic Bight") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    abund_pa <- ggplot(data = md %>% dplyr::filter(region == "MAB", dataset %in% biomod_dataset), mapping = aes(x = month, y = pa, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = pa - var_pa, ymax = pa + var_pa), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Exceeding of Feeding Threshold",
           color = "Legend") +
      ggtitle("Ensemble Mid Atlantic Bight") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white"))
    
    abund_cpr <- ggplot(data = md %>% dplyr::filter(region == "MAB", dataset == "CPR"), mapping = aes(x = month, y = mean, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Climatological Abundance",
           color = "Legend") +
      ggtitle("Ensemble Mid Atlantic Bight") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white"))
    
    abund_cpr_pa <- ggplot(data = md %>% dplyr::filter(region == "MAB", dataset == "CPR"), mapping = aes(x = month, y = pa, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = pa - var_pa, ymax = pa + var_pa), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Exceeding of Feeding Threshold",
           color = "Legend") +
      ggtitle("Ensemble Mid Atlantic Bight") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white"))
    
    pred <- ggplot(data = ensemble_proj_monthly %>% dplyr::filter(region == "MAB"), mapping = aes(x = month, y = mean, color = "Predicted")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "blue", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "Month",
           y = "Predicted Probability of Patch",
           color = "Legend") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'MAB_ensemble_abund_pred_monthly.png'))
    grid.arrange(abund, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'MAB_ensemble_pa_pred_monthly.png'))
    grid.arrange(abund_pa, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'MAB_ensemble_abund_pred_monthly_cpr.png'))
    grid.arrange(abund_cpr, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'MAB_ensemble_pa_pred_monthly_cpr.png'))
    grid.arrange(abund_cpr_pa, pred)
    dev.off()
    
    # ---- Plot GAMs ----
    abund <- ggplot(data = md %>% dplyr::filter(region == "MAB", dataset %in% biomod_dataset), mapping = aes(x = month, y = mean, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Climatological Abundance",
           color = "Legend") +
      ggtitle("GAM Mid Atlantic Bight") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    abund_pa <- ggplot(data = md %>% dplyr::filter(region == "MAB", dataset %in% biomod_dataset), mapping = aes(x = month, y = pa, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = pa - var_pa, ymax = pa + var_pa), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Exceeding of Feeding Threshold",
           color = "Legend") +
      ggtitle("GAM Mid Atlantic Bight") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    abund_cpr <- ggplot(data = md %>% dplyr::filter(region == "MAB", dataset == "CPR"), mapping = aes(x = month, y = mean, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Climatological Abundance",
           color = "Legend") +
      ggtitle("GAM Mid Atlantic Bight") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white"))
    
    abund_cpr_pa <- ggplot(data = md %>% dplyr::filter(region == "MAB", dataset == "CPR"), mapping = aes(x = month, y = pa, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = pa - var_pa, ymax = pa + var_pa), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Exceeding of Feeding Threshold",
           color = "Legend") +
      ggtitle("GAM Mid Atlantic Bight") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white"))
    
    pred <- ggplot(data = gam_proj_monthly %>% dplyr::filter(region == "MAB"), mapping = aes(x = month, y = mean, color = "Predicted")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "blue", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "Month",
           y = "Predicted Probability of Patch",
           color = "Legend") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'MAB_gam_abund_pred_monthly.png'))
    grid.arrange(abund, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'MAB_gam_pa_pred_monthly.png'))
    grid.arrange(abund_pa, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'MAB_gam_abund_pred_monthly_cpr.png'))
    grid.arrange(abund_cpr, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'MAB_gam_pa_pred_monthly_cpr.png'))
    grid.arrange(abund_cpr_pa, pred)
    dev.off()
    
    # ---- Plot BRTs ----
    abund <- ggplot(data = md %>% dplyr::filter(region == "MAB", dataset %in% biomod_dataset), mapping = aes(x = month, y = mean, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Climatological Abundance",
           color = "Legend") +
      ggtitle("BRT Mid Atlantic Bight") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    abund_pa <- ggplot(data = md %>% dplyr::filter(region == "MAB", dataset %in% biomod_dataset), mapping = aes(x = month, y = pa, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = pa - var_pa, ymax = pa + var_pa), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Exceeding of Feeding Threshold",
           color = "Legend") +
      ggtitle("BRT Mid Atlantic Bight") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white"))
    
    abund_cpr <- ggplot(data = md %>% dplyr::filter(region == "MAB", dataset == "CPR"), mapping = aes(x = month, y = mean, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Climatological Abundance",
           color = "Legend") +
      ggtitle("BRT Mid Atlantic Bight") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    abund_cpr_pa <- ggplot(data = md %>% dplyr::filter(region == "MAB", dataset == "CPR"), mapping = aes(x = month, y = pa, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = pa - var_pa, ymax = pa + var_pa), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Exceeding of Feeding Threshold",
           color = "Legend") +
      ggtitle("BRT Mid Atlantic Bight") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    pred <- ggplot(data = brt_proj_monthly %>% dplyr::filter(region == "MAB"), mapping = aes(x = month, y = mean, color = "Predicted")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "blue", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "Month",
           y = "Predicted Probability of Patch",
           color = "Legend") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'MAB_brt_abund_pred_monthly.png'))
    grid.arrange(abund, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'MAB_brt_pa_pred_monthly.png'))
    grid.arrange(abund_pa, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'MAB_brt_abund_pred_monthly_cpr.png'))
    grid.arrange(abund_cpr, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'MAB_brt_pa_pred_monthly_cpr.png'))
    grid.arrange(abund_cpr_pa, pred)
    dev.off()
    
    # ---- Plot RF ----
    abund <- ggplot(data = md %>% dplyr::filter(region == "MAB", dataset %in% biomod_dataset), mapping = aes(x = month, y = mean, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Climatological Abundance",
           color = "Legend") +
      ggtitle("RF Mid Atlantic Bight") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    abund_pa <- ggplot(data = md %>% dplyr::filter(region == "MAB", dataset %in% biomod_dataset), mapping = aes(x = month, y = pa, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = pa - var_pa, ymax = pa + var_pa), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Exceeding of Feeding Threshold",
           color = "Legend") +
      ggtitle("RF Mid Atlantic Bight") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    abund_cpr <- ggplot(data = md %>% dplyr::filter(region == "MAB", dataset == "CPR"), mapping = aes(x = month, y = mean, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Climatological Abundance",
           color = "Legend") +
      ggtitle("RF Mid Atlantic Bight") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    abund_cpr_pa <- ggplot(data = md %>% dplyr::filter(region == "MAB", dataset == "CPR"), mapping = aes(x = month, y = pa, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = pa - var_pa, ymax = pa + var_pa), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Exceeding of Feeding Threshold",
           color = "Legend") +
      ggtitle("RF Mid Atlantic Bight") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    pred <- ggplot(data = rf_proj_monthly %>% dplyr::filter(region == "MAB"), mapping = aes(x = month, y = mean, color = "Predicted")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "blue", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "Month",
           y = "Predicted Probability of Patch",
           color = "Legend") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'MAB_rf_abund_pred_monthly.png'))
    grid.arrange(abund, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'MAB_rf_pa_pred_monthly.png'))
    grid.arrange(abund_pa, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'MAB_rf_abund_pred_monthly_cpr.png'))
    grid.arrange(abund_cpr, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'MAB_rf_pa_pred_monthly_cpr.png'))
    grid.arrange(abund_cpr_pa, pred)
    dev.off()
    
  }
  
  # -------- Plot GBK --------
  if ("GBK" %in% unique(md$region)) {
    # ---- Plot Ensembles ----
    abund <- ggplot(data = md %>% dplyr::filter(region == "GBK", dataset %in% biomod_dataset), mapping = aes(x = month, y = mean, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Climatological Abundance",
           color = "Legend") +
      ggtitle("Ensemble George's Bank") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    abund_pa <- ggplot(data = md %>% dplyr::filter(region == "GBK", dataset %in% biomod_dataset), mapping = aes(x = month, y = pa, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = pa - var_pa, ymax = pa + var_pa), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Exceeding of Feeding Threshold",
           color = "Legend") +
      ggtitle("Ensemble George's Bank") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white"))
    
    pred <- ggplot(data = ensemble_proj_monthly %>% dplyr::filter(region == "GBK"), mapping = aes(x = month, y = mean, color = "Predicted")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "blue", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "Month",
           y = "Predicted Probability of Patch",
           color = "Legend") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'GBK_ensemble_abund_pred_monthly.png'))
    grid.arrange(abund, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'GBK_ensemble_pa_pred_monthly.png'))
    grid.arrange(abund_pa, pred)
    dev.off()
    
    # ---- Plot GAMs ----
    abund <- ggplot(data = md %>% dplyr::filter(region == "GBK", dataset %in% biomod_dataset), mapping = aes(x = month, y = mean, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Climatological Abundance",
           color = "Legend") +
      ggtitle("GAM George's Bank") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    abund_pa <- ggplot(data = md %>% dplyr::filter(region == "GBK", dataset %in% biomod_dataset), mapping = aes(x = month, y = pa, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = pa - var_pa, ymax = pa + var_pa), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Exceeding of Feeding Threshold",
           color = "Legend") +
      ggtitle("GAM George's Bank") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    pred <- ggplot(data = gam_proj_monthly %>% dplyr::filter(region == "GBK"), mapping = aes(x = month, y = mean, color = "Predicted")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "blue", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "Month",
           y = "Predicted Probability of Patch",
           color = "Legend") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'GBK_gam_abund_pred_monthly.png'))
    grid.arrange(abund, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'GBK_gam_pa_pred_monthly.png'))
    grid.arrange(abund_pa, pred)
    dev.off()
    
    # ---- Plot BRTs ----
    abund <- ggplot(data = md %>% dplyr::filter(region == "GBK", dataset %in% biomod_dataset), mapping = aes(x = month, y = mean, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Climatological Abundance",
           color = "Legend") +
      ggtitle("BRT George's Bank") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    abund_pa <- ggplot(data = md %>% dplyr::filter(region == "GBK", dataset %in% biomod_dataset), mapping = aes(x = month, y = pa, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = pa - var_pa, ymax = pa + var_pa), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Exceeding of Feeding Threshold",
           color = "Legend") +
      ggtitle("BRT George's Bank") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white"))
    
    pred <- ggplot(data = brt_proj_monthly %>% dplyr::filter(region == "GBK"), mapping = aes(x = month, y = mean, color = "Predicted")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "blue", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "Month",
           y = "Predicted Probability of Patch",
           color = "Legend") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'GBK_brt_abund_pred_monthly.png'))
    grid.arrange(abund, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'GBK_brt_pa_pred_monthly.png'))
    grid.arrange(abund_pa, pred)
    dev.off()
    
    # ---- Plot RF ----
    abund <- ggplot(data = md %>% dplyr::filter(region == "GBK", dataset %in% biomod_dataset), mapping = aes(x = month, y = mean, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Climatological Abundance",
           color = "Legend") +
      ggtitle("RF George's Bank") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    abund_pa <- ggplot(data = md %>% dplyr::filter(region == "GBK", dataset %in% biomod_dataset), mapping = aes(x = month, y = pa, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = pa - var_pa, ymax = pa + var_pa), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Exceeding of Feeding Threshold",
           color = "Legend") +
      ggtitle("RF George's Bank") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    pred <- ggplot(data = rf_proj_monthly %>% dplyr::filter(region == "GBK"), mapping = aes(x = month, y = mean, color = "Predicted")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "blue", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "Month",
           y = "Predicted Probability of Patch",
           color = "Legend") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'GBK_rf_abund_pred_monthly.png'))
    grid.arrange(abund, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'GBK_rf_pa_pred_monthly.png'))
    grid.arrange(abund_pa, pred)
    dev.off()
    
  }
  
  # -------- Plot GOM --------
  if ("GOM" %in% unique(md$region)) {
    # ---- Plot Ensembles ----
    abund <- ggplot(data = md %>% dplyr::filter(region == "GOM", dataset %in% biomod_dataset), mapping = aes(x = month, y = mean, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Climatological Abundance",
           color = "Legend") +
      ggtitle("Ensemble Gulf of Maine") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white"))
    
    abund_pa <- ggplot(data = md %>% dplyr::filter(region == "GOM", dataset %in% biomod_dataset), mapping = aes(x = month, y = pa, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = pa - var_pa, ymax = pa + var_pa), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Exceeding of Feeding Threshold",
           color = "Legend") +
      ggtitle("Ensemble Gulf of Maine") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    abund_cpr <- ggplot(data = md %>% dplyr::filter(region == "GOM", dataset == "CPR"), mapping = aes(x = month, y = mean, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Climatological Abundance",
           color = "Legend") +
      ggtitle("Ensemble Gulf of Maine") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    abund_cpr_pa <- ggplot(data = md %>% dplyr::filter(region == "GOM", dataset == "CPR"), mapping = aes(x = month, y = pa, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = pa - var_pa, ymax = pa + var_pa), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Exceeding of Feeding Threshold",
           color = "Legend") +
      ggtitle("Ensemble Gulf of Maine") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    pred <- ggplot(data = ensemble_proj_monthly %>% dplyr::filter(region == "GOM"), mapping = aes(x = month, y = mean, color = "Predicted")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "blue", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "Month",
           y = "Predicted Probability of Patch",
           color = "Legend") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'GOM_ensemble_abund_pred_monthly.png'))
    grid.arrange(abund, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'GOM_ensemble_pa_pred_monthly.png'))
    grid.arrange(abund_pa, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'GOM_ensemble_abund_pred_monthly_cpr.png'))
    grid.arrange(abund_cpr, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'GOM_ensemble_pa_pred_monthly_cpr.png'))
    grid.arrange(abund_cpr_pa, pred)
    dev.off()
    
    # ---- Plot GAMs ----
    abund <- ggplot(data = md %>% dplyr::filter(region == "GOM", dataset %in% biomod_dataset), mapping = aes(x = month, y = mean, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Climatological Abundance",
           color = "Legend") +
      ggtitle("GAM Gulf of Maine") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    abund_pa <- ggplot(data = md %>% dplyr::filter(region == "GOM", dataset %in% biomod_dataset), mapping = aes(x = month, y = pa, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = pa - var_pa, ymax = pa + var_pa), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Exceeding of Feeding Threshold",
           color = "Legend") +
      ggtitle("GAM Gulf of Maine") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    abund_cpr <- ggplot(data = md %>% dplyr::filter(region == "GOM", dataset == "CPR"), mapping = aes(x = month, y = mean, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Climatological Abundance",
           color = "Legend") +
      ggtitle("GAM Gulf of Maine") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    abund_cpr_pa <- ggplot(data = md %>% dplyr::filter(region == "GOM", dataset == "CPR"), mapping = aes(x = month, y = pa, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = pa - var_pa, ymax = pa + var_pa), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Exceeding of Feeding Threshold",
           color = "Legend") +
      ggtitle("GAM Gulf of Maine") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white"))
    
    pred <- ggplot(data = gam_proj_monthly %>% dplyr::filter(region == "GOM"), mapping = aes(x = month, y = mean, color = "Predicted")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "blue", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "Month",
           y = "Predicted Probability of Patch",
           color = "Legend") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'GOM_gam_abund_pred_monthly.png'))
    grid.arrange(abund, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'GOM_gam_pa_pred_monthly.png'))
    grid.arrange(abund_pa, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'GOM_gam_abund_pred_monthly_cpr.png'))
    grid.arrange(abund_cpr, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'GOM_gam_pa_pred_monthly_cpr.png'))
    grid.arrange(abund_cpr_pa, pred)
    dev.off()
    
    # ---- Plot BRTs ----
    abund <- ggplot(data = md %>% dplyr::filter(region == "GOM", dataset %in% biomod_dataset), mapping = aes(x = month, y = mean, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Climatological Abundance",
           color = "Legend") +
      ggtitle("BRT Gulf of Maine") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    abund_pa <- ggplot(data = md %>% dplyr::filter(region == "GOM", dataset %in% biomod_dataset), mapping = aes(x = month, y = pa, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = pa - var_pa, ymax = pa + var_pa), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Exceeding of Feeding Threshold",
           color = "Legend") +
      ggtitle("BRT Gulf of Maine") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    abund_cpr <- ggplot(data = md %>% dplyr::filter(region == "GOM", dataset == "CPR"), mapping = aes(x = month, y = mean, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Climatological Abundance",
           color = "Legend") +
      ggtitle("BRT Gulf of Maine") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    abund_cpr_pa <- ggplot(data = md %>% dplyr::filter(region == "GOM", dataset == "CPR"), mapping = aes(x = month, y = pa, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = pa - var_pa, ymax = pa + var_pa), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Exceeding of Feeding Threshold",
           color = "Legend") +
      ggtitle("BRT Gulf of Maine") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    pred <- ggplot(data = brt_proj_monthly %>% dplyr::filter(region == "GOM"), mapping = aes(x = month, y = mean, color = "Predicted")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "blue", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "Month",
           y = "Predicted Probability of Patch",
           color = "Legend") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'GOM_brt_abund_pred_monthly.png'))
    grid.arrange(abund, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'GOM_brt_pa_pred_monthly.png'))
    grid.arrange(abund_pa, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'GOM_brt_abund_pred_monthly_cpr.png'))
    grid.arrange(abund_cpr, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'GOM_brt_pa_pred_monthly_cpr.png'))
    grid.arrange(abund_cpr_pa, pred)
    dev.off()
    
    # ---- Plot RF ----
    abund <- ggplot(data = md %>% dplyr::filter(region == "GOM", dataset %in% biomod_dataset), mapping = aes(x = month, y = mean, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Climatological Abundance",
           color = "Legend") +
      ggtitle("RF Gulf of Maine") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    abund_pa <- ggplot(data = md %>% dplyr::filter(region == "GOM", dataset %in% biomod_dataset), mapping = aes(x = month, y = pa, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = pa - var_pa, ymax = pa + var_pa), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Exceeding of Feeding Threshold",
           color = "Legend") +
      ggtitle("RF Gulf of Maine") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    abund_cpr <- ggplot(data = md %>% dplyr::filter(region == "GOM", dataset == "CPR"), mapping = aes(x = month, y = mean, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Climatological Abundance",
           color = "Legend") +
      ggtitle("RF Gulf of Maine") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white"))
    
    abund_cpr_pa <- ggplot(data = md %>% dplyr::filter(region == "GOM", dataset == "CPR"), mapping = aes(x = month, y = pa, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = pa - var_pa, ymax = pa + var_pa), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Exceeding of Feeding Threshold",
           color = "Legend") +
      ggtitle("RF Gulf of Maine") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white"))
    
    pred <- ggplot(data = rf_proj_monthly %>% dplyr::filter(region == "GOM"), mapping = aes(x = month, y = mean, color = "Predicted")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "blue", color = NA) +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "Month",
           y = "Predicted Probability of Patch",
           color = "Legend") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'GOM_rf_abund_pred_monthly.png'))
    grid.arrange(abund, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'GOM_rf_pa_pred_monthly.png'))
    grid.arrange(abund_pa, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'GOM_rf_abund_pred_monthly_cpr.png'))
    grid.arrange(abund_cpr, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'GOM_rf_pa_pred_monthly_cpr.png'))
    grid.arrange(abund_cpr_pa, pred)
    dev.off()
  }
  
  synth_dat <- full_join(md %>% dplyr::filter(dataset %in% biomod_dataset), ensemble_proj_monthly, by = c("region", "month"))
  
  ggplot(data = synth_dat, mapping = aes(x = mean.x, y = mean.y, color = region)) +
    geom_point() +
    geom_errorbar(aes(ymin = mean.y - var.y, ymax = mean.y + var.y), width=.2) +
    geom_errorbarh(aes(xmin = mean.x - var.x, xmax = mean.x + var.x), height = 0.025) +
    labs(y = "Predicted Probability of Patch",
         x = "Climatological abundance",
         color = "Region")  +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key = element_rect(color = "transparent", fill = "white")) +
    ggsave(file.path(fp_out, species, version, "Biomod", "Plots", 'ensemble_abund_vs_pred_monthly.png'))
  
  ggplot(data = synth_dat, mapping = aes(x = pa, y = mean.y, color = region)) +
    geom_point() +
    geom_errorbar(aes(ymin = mean.y - var.y, ymax = mean.y + var.y), width=.01) +
    geom_errorbarh(aes(xmin = pa - var_pa, xmax = pa + var_pa), height = 0.025) +
    labs(y = "Predicted Probability of Patch",
         x = "Climatological Actual Probability of Patch",
         color = "Region")  +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key = element_rect(color = "transparent", fill = "white")) +
    ggsave(file.path(fp_out, species, version, "Biomod", "Plots", 'ensemble_pa_vs_pred_monthly.png'))
  
  synth_dat <- full_join(md %>% dplyr::filter(dataset == "CPR"), ensemble_proj_monthly, by = c("region", "month"))
  
  ggplot(data = synth_dat, mapping = aes(x = mean.x, y = mean.y, color = region)) +
    geom_point() +
    geom_errorbar(aes(ymin = mean.y - var.y, ymax = mean.y + var.y), width=.2) +
    geom_errorbarh(aes(xmin = mean.x - var.x, xmax = mean.x + var.x), height = 0.025) +
    labs(y = "Predicted Probability of Patch",
         x = "Climatological abundance",
         color = "Region")  +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key = element_rect(color = "transparent", fill = "white")) +
    ggsave(file.path(fp_out, species, version, "Biomod", "Plots", 'ensemble_abund_vs_pred_monthly_cpr.png'))
  
  ggplot(data = synth_dat, mapping = aes(x = pa, y = mean.y, color = region)) +
    geom_point() +
    geom_errorbar(aes(ymin = mean.y - var.y, ymax = mean.y + var.y), width=.01) +
    geom_errorbarh(aes(xmin = pa - var_pa, xmax = pa + var_pa), height = 0.025) +
    labs(y = "Predicted Probability of Patch",
         x = "Climatological Actual Probability of Patch",
         color = "Region")  +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key = element_rect(color = "transparent", fill = "white")) +
    ggsave(file.path(fp_out, species, version, "Biomod", "Plots", 'ensemble_pa_vs_pred_monthly_cpr.png'))
  
  synth_dat <- full_join(md %>% dplyr::filter(dataset %in% biomod_dataset), gam_proj_monthly, by = c("region", "month"))
  
  ggplot(data = synth_dat, mapping = aes(x = mean.x, y = mean.y, color = region)) +
    geom_point() +
    geom_errorbar(aes(ymin = mean.y - var.y, ymax = mean.y + var.y), width=.2) +
    geom_errorbarh(aes(xmin = mean.x - var.x, xmax = mean.x + var.x), height = 0.025) +
    labs(y = "Predicted Probability of Patch",
         x = "Climatological abundance",
         color = "Region")  +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key = element_rect(color = "transparent", fill = "white")) +
    ggsave(file.path(fp_out, species, version, "Biomod", "Plots", 'gam_abund_vs_pred_monthly.png'))
  
  ggplot(data = synth_dat, mapping = aes(x = pa, y = mean.y, color = region)) +
    geom_point() +
    geom_errorbar(aes(ymin = mean.y - var.y, ymax = mean.y + var.y), width=.01) +
    geom_errorbarh(aes(xmin = pa - var_pa, xmax = pa + var_pa), height = 0.025) +
    labs(y = "Predicted Probability of Patch",
         x = "Climatological Actual Probability of Patch",
         color = "Region")  +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key = element_rect(color = "transparent", fill = "white")) +
    ggsave(file.path(fp_out, species, version, "Biomod", "Plots", 'gam_pa_vs_pred_monthly.png'))
  
  synth_dat <- full_join(md %>% dplyr::filter(dataset == "CPR"), gam_proj_monthly, by = c("region", "month"))
  
  ggplot(data = synth_dat, mapping = aes(x = mean.x, y = mean.y, color = region)) +
    geom_point() +
    geom_errorbar(aes(ymin = mean.y - var.y, ymax = mean.y + var.y), width=.2) +
    geom_errorbarh(aes(xmin = mean.x - var.x, xmax = mean.x + var.x), height = 0.025) +
    labs(y = "Predicted Probability of Patch",
         x = "Climatological abundance",
         color = "Region")  +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key = element_rect(color = "transparent", fill = "white")) +
    ggsave(file.path(fp_out, species, version, "Biomod", "Plots", 'gam_abund_vs_pred_monthly_cpr.png'))
  
  ggplot(data = synth_dat, mapping = aes(x = pa, y = mean.y, color = region)) +
    geom_point() +
    geom_errorbar(aes(ymin = mean.y - var.y, ymax = mean.y + var.y), width=.01) +
    geom_errorbarh(aes(xmin = pa - var_pa, xmax = pa + var_pa), height = 0.025) +
    labs(y = "Predicted Probability of Patch",
         x = "Climatological Actual Probability of Patch",
         color = "Region")  +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key = element_rect(color = "transparent", fill = "white")) +
    ggsave(file.path(fp_out, species, version, "Biomod", "Plots", 'gam_pa_vs_pred_monthly_cpr.png'))
  
  synth_dat <- full_join(md %>% dplyr::filter(dataset %in% biomod_dataset), brt_proj_monthly, by = c("region", "month"))
  
  ggplot(data = synth_dat, mapping = aes(x = mean.x, y = mean.y, color = region)) +
    geom_point() +
    geom_errorbar(aes(ymin = mean.y - var.y, ymax = mean.y + var.y), width=.2) +
    geom_errorbarh(aes(xmin = mean.x - var.x, xmax = mean.x + var.x), height = 0.025) +
    labs(y = "Predicted Probability of Patch",
         x = "Climatological abundance",
         color = "Region")  +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key = element_rect(color = "transparent", fill = "white")) +
    ggsave(file.path(fp_out, species, version, "Biomod", "Plots", 'brt_abund_vs_pred_monthly.png'))
  
  ggplot(data = synth_dat, mapping = aes(x = pa, y = mean.y, color = region)) +
    geom_point() +
    geom_errorbar(aes(ymin = mean.y - var.y, ymax = mean.y + var.y), width=.01) +
    geom_errorbarh(aes(xmin = pa - var_pa, xmax = pa + var_pa), height = 0.025) +
    labs(y = "Predicted Probability of Patch",
         x = "Climatological Actual Probability of Patch",
         color = "Region")  +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key = element_rect(color = "transparent", fill = "white")) +
    ggsave(file.path(fp_out, species, version, "Biomod", "Plots", 'brt_pa_vs_pred_monthly.png'))
  
  synth_dat <- full_join(md %>% dplyr::filter(dataset == "CPR"), brt_proj_monthly, by = c("region", "month"))
  
  ggplot(data = synth_dat, mapping = aes(x = mean.x, y = mean.y, color = region)) +
    geom_point() +
    geom_errorbar(aes(ymin = mean.y - var.y, ymax = mean.y + var.y), width=.2) +
    geom_errorbarh(aes(xmin = mean.x - var.x, xmax = mean.x + var.x), height = 0.025) +
    labs(y = "Predicted Probability of Patch",
         x = "Climatological abundance",
         color = "Region")  +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key = element_rect(color = "transparent", fill = "white")) +
    ggsave(file.path(fp_out, species, version, "Biomod", "Plots", 'brt_abund_vs_pred_monthly_cpr.png'))
  
  ggplot(data = synth_dat, mapping = aes(x = pa, y = mean.y, color = region)) +
    geom_point() +
    geom_errorbar(aes(ymin = mean.y - var.y, ymax = mean.y + var.y), width=.01) +
    geom_errorbarh(aes(xmin = pa - var_pa, xmax = pa + var_pa), height = 0.025) +
    labs(y = "Predicted Probability of Patch",
         x = "Climatological Actual Probability of Patch",
         color = "Region")  +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key = element_rect(color = "transparent", fill = "white")) +
    ggsave(file.path(fp_out, species, version, "Biomod", "Plots", 'brt_pa_vs_pred_monthly_cpr.png'))
  
  synth_dat <- full_join(md %>% dplyr::filter(dataset %in% biomod_dataset), rf_proj_monthly, by = c("region", "month"))
  
  ggplot(data = synth_dat, mapping = aes(x = mean.x, y = mean.y, color = region)) +
    geom_point() +
    geom_errorbar(aes(ymin = mean.y - var.y, ymax = mean.y + var.y), width=.2) +
    geom_errorbarh(aes(xmin = mean.x - var.x, xmax = mean.x + var.x), height = 0.025) +
    labs(y = "Predicted Probability of Patch",
         x = "Climatological abundance",
         color = "Region")  +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key = element_rect(color = "transparent", fill = "white")) +
    ggsave(file.path(fp_out, species, version, "Biomod", "Plots", 'rf_abund_vs_pred_monthly.png'))
  
  ggplot(data = synth_dat, mapping = aes(x = pa, y = mean.y, color = region)) +
    geom_point() +
    geom_errorbar(aes(ymin = mean.y - var.y, ymax = mean.y + var.y), width=.01) +
    geom_errorbarh(aes(xmin = pa - var_pa, xmax = pa + var_pa), height = 0.025) +
    labs(y = "Predicted Probability of Patch",
         x = "Climatological Actual Probability of Patch",
         color = "Region")  +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key = element_rect(color = "transparent", fill = "white")) +
    ggsave(file.path(fp_out, species, version, "Biomod", "Plots", 'rf_pa_vs_pred_monthly.png'))
  
  synth_dat <- full_join(md %>% dplyr::filter(dataset == "CPR"), rf_proj_monthly, by = c("region", "month"))
  
  ggplot(data = synth_dat, mapping = aes(x = mean.x, y = mean.y, color = region)) +
    geom_point() +
    geom_errorbar(aes(ymin = mean.y - var.y, ymax = mean.y + var.y), width=.2) +
    geom_errorbarh(aes(xmin = mean.x - var.x, xmax = mean.x + var.x), height = 0.025) +
    labs(y = "Predicted Probability of Patch",
         x = "Climatological abundance",
         color = "Region")  +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key = element_rect(color = "transparent", fill = "white")) +
    ggsave(file.path(fp_out, species, version, "Biomod", "Plots", 'rf_abund_vs_pred_monthly_cpr.png'))
  
  ggplot(data = synth_dat, mapping = aes(x = pa, y = mean.y, color = region)) +
    geom_point() +
    geom_errorbar(aes(ymin = mean.y - var.y, ymax = mean.y + var.y), width=.01) +
    geom_errorbarh(aes(xmin = pa - var_pa, xmax = pa + var_pa), height = 0.025) +
    labs(y = "Predicted Probability of Patch",
         x = "Climatological Actual Probability of Patch",
         color = "Region")  +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key = element_rect(color = "transparent", fill = "white")) +
    ggsave(file.path(fp_out, species, version, "Biomod", "Plots", 'rf_pa_vs_pred_monthly_cpr.png'))
  
  # -------- INTERANNUAL VARIABILITY --------
  
  # -------- Load model data --------
  if (length(years) == 1) {
    if (anomaly) {
      md <- readr::read_csv(file.path(fp_md, paste0(years[1], ".csv")))
    } else {
      md <- readr::read_csv(file.path(fp_md, paste0(years[1], ".csv"))) %>% 
        dplyr::filter(dataset %in% c(biomod_dataset, "CPR"))
    }
  } else {
    if (anomaly) {
      md <- bind_years(fp = file.path(fp_md), years = years)
    } else {
      md <- bind_years(fp = file.path(fp_md), years = years) %>%
        dplyr::filter(dataset %in% c(biomod_dataset, "CPR"))
    }
  }
  
  # -------- Compute anomaly --------
  if (species == "cfin") {
    md <- md %>% dplyr::group_by(dataset) %>%
      dplyr::mutate(mean = mean(log10(`cfin_CV_VI` + 1), na.rm = TRUE),
                    sd = sd(log10(`cfin_CV_VI` + 1), na.rm = TRUE),
                    var = var(log10(`cfin_CV_VI` + 1), na.rm = TRUE),
                    anomaly = (log10(`cfin_CV_VI` + 1) - mean) / sd) %>%
      dplyr::ungroup()
  } else if (species == "ctyp") {
    md <- md %>% dplyr::group_by(dataset) %>%
      dplyr::mutate(mean = mean(log10(`ctyp_total` + 1), na.rm = TRUE),
                    sd = sd(log10(`ctyp_total` + 1), na.rm = TRUE),
                    var = var(log10(`ctyp_total` + 1), na.rm = TRUE),
                    anomaly = (log10(`ctyp_total` + 1) - mean) / sd) %>%
      dplyr::ungroup()
  } else if (species == "pseudo") {
    md <- md %>% dplyr::group_by(dataset) %>%
      dplyr::mutate(mean = mean(log10(`pseudo_total` + 1), na.rm = TRUE),
                    sd = sd(log10(`pseudo_total` + 1), na.rm = TRUE),
                    var = var(log10(`pseudo_total` + 1), na.rm = TRUE),
                    anomaly = (log10(`pseudo_total` + 1) - mean) / sd) %>%
      dplyr::ungroup()
  }
  
  # -------- Take the log of count data --------
  if (anomaly) {
    md$abund <- md$anomaly
  } else {
    if (species == "cfin") {
      md$abund <- as.data.frame(md[paste0(species, "_CV_VI")])$cfin_CV_VI
    } else if (species == "ctyp") {
      md$abund <- as.data.frame(md[paste0(species, "_total")])$ctyp_total
    } else if (species == "pseudo") {
      md$abund <- as.data.frame(md[paste0(species, "_total")])$pseudo_total
    }
  }
  
  # -------- Exclude NAs and select columns --------
  md <- md %>% dplyr::select(lat, lon, year, month, dataset, abund, wind, fetch, chl, int_chl, bots, bott, sss, sst, lag_sst, uv, bat, dist, slope) %>%
    as.data.frame() %>%
    na.exclude() %>%
    dplyr::filter(!is.infinite(abund)) %>%
    dplyr::mutate(season = if_else(month %in% c(1:3), 1,
                                   if_else(month %in% c(4:6), 2,
                                           if_else(month %in% c(7:9), 3, 4))),
                  region = if_else(lat <= 41.5 & lon < -70, "MAB", 
                                   if_else(lat >= 39 & lat <= 42 & lon >= -70 & lon <= -68, "GBK", "GOM"))) %>%
    dplyr::mutate(pa = if_else(abund >= threshold, 1, 0),
                  abund = log10(abund + 1)) %>%
    dplyr::group_by(region, year, dataset) %>%
    dplyr::summarize(mean = mean(abund, na.rm=TRUE),
                     stdev = sd(abund, na.rm = TRUE),
                     var = var(abund, na.rm = TRUE),
                     var_pa = var(pa, na.rm = TRUE),
                     pa = mean(pa, na.rm = TRUE))
  
  
  # -------- Compute regions --------
  ensemble_proj_yearly <- ensemble_proj %>%
    dplyr::mutate(region = if_else(y <= 41 & x < -70, "MAB", 
                                   if_else(y >= 40 & x <= 42 & y >= -70 & x <= -68, "GBK", "GOM"))) %>%
    dplyr::group_by(region, year) %>%
    dplyr::summarize(mean = mean(proj, na.rm = TRUE),
                     stdev = sd(proj, na.rm = TRUE),
                     var = var(proj, na.rm = TRUE))
  gam_proj_yearly <- gam_proj %>%
    dplyr::mutate(region = if_else(y <= 41 & x < -70, "MAB", 
                                   if_else(y >= 40 & x <= 42 & y >= -70 & x <= -68, "GBK", "GOM"))) %>%
    dplyr::group_by(region, year) %>%
    dplyr::summarize(mean = mean(proj, na.rm = TRUE),
                     stdev = sd(proj, na.rm = TRUE),
                     var = var(proj, na.rm = TRUE))
  
  brt_proj_yearly <- brt_proj %>%
    dplyr::mutate(region = if_else(y <= 41 & x < -70, "MAB", 
                                   if_else(y >= 40 & x <= 42 & y >= -70 & x <= -68, "GBK", "GOM"))) %>%
    dplyr::group_by(region, year) %>%
    dplyr::summarize(mean = mean(proj, na.rm = TRUE),
                     stdev = sd(proj, na.rm = TRUE),
                     var = var(proj, na.rm = TRUE))
  
  rf_proj_yearly <- rf_proj %>%
    dplyr::mutate(region = if_else(y <= 41 & x < -70, "MAB", 
                                   if_else(y >= 40 & x <= 42 & y >= -70 & x <= -68, "GBK", "GOM"))) %>%
    dplyr::group_by(region, year) %>%
    dplyr::summarize(mean = mean(proj, na.rm = TRUE),
                     stdev = sd(proj, na.rm = TRUE),
                     var = var(proj, na.rm = TRUE))
  
  # -------- Plot MAB --------
  if ("MAB" %in% unique(md$region)) {
    
    # ---- Plot Ensembles ----
    abund <- ggplot(data = md %>% dplyr::filter(region == "MAB", dataset %in% biomod_dataset), mapping = aes(x = year, y = mean, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - stdev, ymax = mean + stdev), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Annual Abundance",
           color = "Legend") +
      ggtitle("Ensemble Mid Atlantic Bight") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    abund_pa <- ggplot(data = md %>% dplyr::filter(region == "MAB", dataset %in% biomod_dataset), mapping = aes(x = year, y = pa, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = pa - var_pa, ymax = pa + var_pa), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Exceeding of Feeding Threshold",
           color = "Legend") +
      ggtitle("Ensemble Mid Atlantic Bight") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    abund_cpr <- ggplot(data = md %>% dplyr::filter(region == "MAB", dataset == "CPR"), mapping = aes(x = year, y = mean, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - stdev, ymax = mean + stdev), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Annual Abundance",
           color = "Legend") +
      ggtitle("Ensemble Mid Atlantic Bight") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    abund_cpr_pa <- ggplot(data = md %>% dplyr::filter(region == "MAB", dataset == "CPR"), mapping = aes(x = year, y = pa, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = pa - var_pa, ymax = pa + var_pa), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Exceeding of Feeding Threshold",
           color = "Legend") +
      ggtitle("Ensemble Mid Atlantic Bight") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    pred <- ggplot(data = ensemble_proj_yearly %>% dplyr::filter(region == "MAB"), mapping = aes(x = year, y = mean, color = "Predicted")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "blue", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "Year",
           y = "Predicted Probability of Patch",
           color = "Legend") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'MAB_ensemble_abund_pred_yearly.png'))
    grid.arrange(abund, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'MAB_ensemble_pa_pred_yearly.png'))
    grid.arrange(abund_pa, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'MAB_ensemble_abund_pred_yearly_cpr.png'))
    grid.arrange(abund_cpr, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'MAB_ensemble_pa_pred_yearly_cpr.png'))
    grid.arrange(abund_cpr_pa, pred)
    dev.off()
    
    # ---- Plot GAMs ----
    abund <- ggplot(data = md %>% dplyr::filter(region == "MAB", dataset %in% biomod_dataset), mapping = aes(x = year, y = mean, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Annual Abundance",
           color = "Legend") +
      ggtitle("GAM Mid Atlantic Bight") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    abund_pa <- ggplot(data = md %>% dplyr::filter(region == "MAB", dataset %in% biomod_dataset), mapping = aes(x = year, y = pa, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = pa - var_pa, ymax = pa + var_pa), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Exceeding of Feeding Threshold",
           color = "Legend") +
      ggtitle("GAM Mid Atlantic Bight") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    abund_cpr <- ggplot(data = md %>% dplyr::filter(region == "MAB", dataset == "CPR"), mapping = aes(x = year, y = mean, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Annual Abundance",
           color = "Legend") +
      ggtitle("GAM Mid Atlantic Bight") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    abund_cpr_pa <- ggplot(data = md %>% dplyr::filter(region == "MAB", dataset == "CPR"), mapping = aes(x = year, y = pa, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = pa - var_pa, ymax = pa + var_pa), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Exceeding of Feeding Threshold",
           color = "Legend") +
      ggtitle("GAM Mid Atlantic Bight") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    pred <- ggplot(data = gam_proj_yearly %>% dplyr::filter(region == "MAB"), mapping = aes(x = year, y = mean, color = "Predicted")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "blue", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "Year",
           y = "Predicted Probability of Patch",
           color = "Legend") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'MAB_gam_abund_pred_yearly.png'))
    grid.arrange(abund, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'MAB_gam_pa_pred_yearly.png'))
    grid.arrange(abund_pa, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'MAB_gam_abund_pred_yearly_cpr.png'))
    grid.arrange(abund_cpr, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'MAB_gam_pa_pred_yearly_cpr.png'))
    grid.arrange(abund_cpr_pa, pred)
    dev.off()
    
    # ---- Plot BRTs ----
    abund <- ggplot(data = md %>% dplyr::filter(region == "MAB", dataset %in% biomod_dataset), mapping = aes(x = year, y = mean, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Annual Abundance",
           color = "Legend") +
      ggtitle("BRT Mid Atlantic Bight") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    abund_pa <- ggplot(data = md %>% dplyr::filter(region == "MAB", dataset %in% biomod_dataset), mapping = aes(x = year, y = pa, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = pa - var_pa, ymax = pa + var_pa), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Exceeding of Feeding Threshold",
           color = "Legend") +
      ggtitle("BRT Mid Atlantic Bight") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    abund_cpr <- ggplot(data = md %>% dplyr::filter(region == "MAB", dataset == "CPR"), mapping = aes(x = year, y = mean, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Annual Abundance",
           color = "Legend") +
      ggtitle("BRT Mid Atlantic Bight") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    abund_cpr_pa <- ggplot(data = md %>% dplyr::filter(region == "MAB", dataset == "CPR"), mapping = aes(x = year, y = pa, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = pa - var_pa, ymax = pa + var_pa), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Exceeding of Feeding Threshold",
           color = "Legend") +
      ggtitle("BRT Mid Atlantic Bight") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white"))
    
    pred <- ggplot(data = brt_proj_yearly %>% dplyr::filter(region == "MAB"), mapping = aes(x = year, y = mean, color = "Predicted")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "blue", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "Year",
           y = "Predicted Probability of Patch",
           color = "Legend") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'MAB_brt_abund_pred_yearly.png'))
    grid.arrange(abund, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'MAB_brt_pa_pred_yearly.png'))
    grid.arrange(abund_pa, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'MAB_brt_abund_pred_yearly_cpr.png'))
    grid.arrange(abund_cpr, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'MAB_brt_pa_pred_yearly_cpr.png'))
    grid.arrange(abund_cpr_pa, pred)
    dev.off()
    
    # ---- Plot RF ----
    abund <- ggplot(data = md %>% dplyr::filter(region == "MAB", dataset %in% biomod_dataset), mapping = aes(x = year, y = mean, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Annual Abundance",
           color = "Legend") +
      ggtitle("RF Mid Atlantic Bight") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    abund_pa <- ggplot(data = md %>% dplyr::filter(region == "MAB", dataset %in% biomod_dataset), mapping = aes(x = year, y = pa, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = pa - var_pa, ymax = pa + var_pa), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Exceeding of Feeding Threshold",
           color = "Legend") +
      ggtitle("RF Mid Atlantic Bight") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    abund_cpr <- ggplot(data = md %>% dplyr::filter(region == "MAB", dataset == "CPR"), mapping = aes(x = year, y = mean, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Annual Abundance",
           color = "Legend") +
      ggtitle("RF Mid Atlantic Bight") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    abund_cpr_pa <- ggplot(data = md %>% dplyr::filter(region == "MAB", dataset == "CPR"), mapping = aes(x = year, y = pa, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = pa - var_pa, ymax = pa + var_pa), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Exceeding of Feeding Threshold",
           color = "Legend") +
      ggtitle("RF Mid Atlantic Bight") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    pred <- ggplot(data = rf_proj_yearly %>% dplyr::filter(region == "MAB"), mapping = aes(x = year, y = mean, color = "Predicted")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "blue", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "Year",
           y = "Predicted Probability of Patch",
           color = "Legend") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'MAB_rf_abund_pred_yearly.png'))
    grid.arrange(abund, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'MAB_rf_pa_pred_yearly.png'))
    grid.arrange(abund_pa, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'MAB_rf_abund_pred_yearly_cpr.png'))
    grid.arrange(abund_cpr, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'MAB_rf_pa_pred_yearly_cpr.png'))
    grid.arrange(abund_cpr_pa, pred)
    dev.off()
    
  }
  
  # -------- Plot GBK --------
  if ("GBK" %in% unique(md$region)) {
    # ---- Plot Ensembles ----
    abund <- ggplot(data = md %>% dplyr::filter(region == "GBK", dataset %in% biomod_dataset), mapping = aes(x = year, y = mean, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Annual Abundance",
           color = "Legend") +
      ggtitle("Ensemble George's Bank") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    abund_pa <- ggplot(data = md %>% dplyr::filter(region == "GBK", dataset %in% biomod_dataset), mapping = aes(x = year, y = pa, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = pa - var_pa, ymax = pa + var_pa), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Exceeding of Feeding Threshold",
           color = "Legend") +
      ggtitle("Ensemble George's Bank") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    pred <- ggplot(data = ensemble_proj_yearly %>% dplyr::filter(region == "GBK"), mapping = aes(x = year, y = mean, color = "Predicted")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "blue", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "Year",
           y = "Predicted Probability of Patch",
           color = "Legend") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'GBK_ensemble_abund_pred_yearly.png'))
    grid.arrange(abund, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'GBK_ensemble_pa_pred_yearly.png'))
    grid.arrange(abund_pa, pred)
    dev.off()
    
    # ---- Plot GAMs ----
    abund <- ggplot(data = md %>% dplyr::filter(region == "GBK", dataset %in% biomod_dataset), mapping = aes(x = year, y = mean, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Annual Abundance",
           color = "Legend") +
      ggtitle("GAM George's Bank") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white"))
    
    abund_pa <- ggplot(data = md %>% dplyr::filter(region == "GBK", dataset %in% biomod_dataset), mapping = aes(x = year, y = pa, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = pa - var_pa, ymax = pa + var_pa), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Exceeding of Feeding Threshold",
           color = "Legend") +
      ggtitle("GAM George's Bank") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white"))
    
    pred <- ggplot(data = gam_proj_yearly %>% dplyr::filter(region == "GBK"), mapping = aes(x = year, y = mean, color = "Predicted")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "blue", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "Year",
           y = "Predicted Probability of Patch",
           color = "Legend") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'GBK_gam_abund_pred_yearly.png'))
    grid.arrange(abund, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'GBK_gam_pa_pred_yearly.png'))
    grid.arrange(abund_pa, pred)
    dev.off()
    
    # ---- Plot BRTs ----
    abund <- ggplot(data = md %>% dplyr::filter(region == "GBK", dataset %in% biomod_dataset), mapping = aes(x = year, y = mean, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Annual Abundance",
           color = "Legend") +
      ggtitle("BRT George's Bank") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white"))
    
    abund_pa <- ggplot(data = md %>% dplyr::filter(region == "GBK", dataset %in% biomod_dataset), mapping = aes(x = year, y = pa, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = pa - var_pa, ymax = pa + var_pa), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Exceeding of Feeding Threshold",
           color = "Legend") +
      ggtitle("BRT George's Bank") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white"))
    
    pred <- ggplot(data = brt_proj_yearly %>% dplyr::filter(region == "GBK"), mapping = aes(x = year, y = mean, color = "Predicted")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "blue", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "Year",
           y = "Predicted Probability of Patch",
           color = "Legend") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'GBK_brt_abund_pred_yearly.png'))
    grid.arrange(abund, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'GBK_brt_pa_pred_yearly.png'))
    grid.arrange(abund_pa, pred)
    dev.off()
    
    # ---- Plot RF ----
    abund <- ggplot(data = md %>% dplyr::filter(region == "GBK", dataset %in% biomod_dataset), mapping = aes(x = year, y = mean, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Annual Abundance",
           color = "Legend") +
      ggtitle("RF George's Bank") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    abund_pa <- ggplot(data = md %>% dplyr::filter(region == "GBK", dataset %in% biomod_dataset), mapping = aes(x = year, y = pa, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = pa - var_pa, ymax = pa + var_pa), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Exceeding of Feeding Threshold",
           color = "Legend") +
      ggtitle("RF George's Bank") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    pred <- ggplot(data = rf_proj_yearly %>% dplyr::filter(region == "GBK"), mapping = aes(x = year, y = mean, color = "Predicted")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "blue", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "Year",
           y = "Predicted Probability of Patch",
           color = "Legend") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'GBK_rf_abund_pred_yearly.png'))
    grid.arrange(abund, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'GBK_rf_pa_pred_yearly.png'))
    grid.arrange(abund_pa, pred)
    dev.off()
    
  }
  
  # -------- Plot GOM --------
  if ("GOM" %in% unique(md$region)) {
    # ---- Plot Ensembles ----
    abund <- ggplot(data = md %>% dplyr::filter(region == "GOM", dataset %in% biomod_dataset), mapping = aes(x = year, y = mean, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Annual Abundance",
           color = "Legend") +
      ggtitle("Ensemble Gulf of Maine") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    abund_pa <- ggplot(data = md %>% dplyr::filter(region == "GOM", dataset %in% biomod_dataset), mapping = aes(x = year, y = pa, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = pa - var_pa, ymax = pa + var_pa), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Exceeding of Feeding Threshold",
           color = "Legend") +
      ggtitle("Ensemble Gulf of Maine") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    abund_cpr <- ggplot(data = md %>% dplyr::filter(region == "GOM", dataset == "CPR"), mapping = aes(x = year, y = mean, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Annual Abundance",
           color = "Legend") +
      ggtitle("Ensemble Gulf of Maine") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    abund_cpr_pa <- ggplot(data = md %>% dplyr::filter(region == "GOM", dataset == "CPR"), mapping = aes(x = year, y = pa, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = pa - var_pa, ymax = pa + var_pa), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Exceeding of Feeding Threshold",
           color = "Legend") +
      ggtitle("Ensemble Gulf of Maine") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    pred <- ggplot(data = ensemble_proj_yearly %>% dplyr::filter(region == "GOM"), mapping = aes(x = year, y = mean, color = "Predicted")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "blue", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "Year",
           y = "Predicted Probability of Patch",
           color = "Legend") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'GOM_ensemble_abund_pred_yearly.png'))
    grid.arrange(abund, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'GOM_ensemble_pa_pred_yearly.png'))
    grid.arrange(abund_pa, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'GOM_ensemble_abund_pred_yearly_cpr.png'))
    grid.arrange(abund_cpr, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'GOM_ensemble_pa_pred_yearly_cpr.png'))
    grid.arrange(abund_cpr_pa, pred)
    dev.off()
    
    # ---- Plot GAMs ----
    abund <- ggplot(data = md %>% dplyr::filter(region == "GOM", dataset %in% biomod_dataset), mapping = aes(x = year, y = mean, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Annual Abundance",
           color = "Legend") +
      ggtitle("GAM Geulf of Maine") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    abund_pa <- ggplot(data = md %>% dplyr::filter(region == "GOM", dataset %in% biomod_dataset), mapping = aes(x = year, y = pa, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = pa - var_pa, ymax = pa + var_pa), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Exceeding of Feeding Threshold",
           color = "Legend") +
      ggtitle("GAM Geulf of Maine") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white"))
    
    abund_cpr <- ggplot(data = md %>% dplyr::filter(region == "GOM", dataset == "CPR"), mapping = aes(x = year, y = mean, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Annual Abundance",
           color = "Legend") +
      ggtitle("GAM Geulf of Maine") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white"))
    
    abund_cpr_pa <- ggplot(data = md %>% dplyr::filter(region == "GOM", dataset == "CPR"), mapping = aes(x = year, y = pa, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = pa - var_pa, ymax = pa + var_pa), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Exceeding of Feeding Threshold",
           color = "Legend") +
      ggtitle("GAM Geulf of Maine") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white"))
    
    pred <- ggplot(data = gam_proj_yearly %>% dplyr::filter(region == "GOM"), mapping = aes(x = year, y = mean, color = "Predicted")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "blue", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "Year",
           y = "Predicted Probability of Patch",
           color = "Legend") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'GOM_gam_abund_pred_yearly.png'))
    grid.arrange(abund, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'GOM_gam_pa_pred_yearly.png'))
    grid.arrange(abund_pa, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'GOM_gam_abund_pred_yearly_cpr.png'))
    grid.arrange(abund_cpr, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'GOM_gam_pa_pred_yearly_cpr.png'))
    grid.arrange(abund_cpr_pa, pred)
    dev.off()
    
    # ---- Plot BRTs ----
    abund <- ggplot(data = md %>% dplyr::filter(region == "GOM", dataset %in% biomod_dataset), mapping = aes(x = year, y = mean, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Annual Abundance",
           color = "Legend") +
      ggtitle("BRT Gulf of Maine") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    abund_pa <- ggplot(data = md %>% dplyr::filter(region == "GOM", dataset %in% biomod_dataset), mapping = aes(x = year, y = pa, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = pa - var_pa, ymax = pa + var_pa), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Exceeding of Feeding Threshold",
           color = "Legend") +
      ggtitle("BRT Gulf of Maine") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white"))
    
    abund_cpr <- ggplot(data = md %>% dplyr::filter(region == "GOM", dataset == "CPR"), mapping = aes(x = year, y = mean, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Annual Abundance",
           color = "Legend") +
      ggtitle("BRT Gulf of Maine") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white"))
    
    abund_cpr_pa <- ggplot(data = md %>% dplyr::filter(region == "GOM", dataset == "CPR"), mapping = aes(x = year, y = pa, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = pa - var_pa, ymax = pa + var_pa), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Exceeding of Feeding Threshold",
           color = "Legend") +
      ggtitle("BRT Gulf of Maine") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white"))
    
    pred <- ggplot(data = brt_proj_yearly %>% dplyr::filter(region == "GOM"), mapping = aes(x = year, y = mean, color = "Predicted")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "blue", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "Year",
           y = "Predicted Probability of Patch",
           color = "Legend") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'GOM_brt_abund_pred_yearly.png'))
    grid.arrange(abund, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'GOM_brt_pa_pred_yearly.png'))
    grid.arrange(abund_pa, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'GOM_brt_abund_pred_yearly_cpr.png'))
    grid.arrange(abund_cpr, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'GOM_brt_pa_pred_yearly_cpr.png'))
    grid.arrange(abund_cpr_pa, pred)
    dev.off()
    
    # ---- Plot RF ----
    abund <- ggplot(data = md %>% dplyr::filter(region == "GOM", dataset %in% biomod_dataset), mapping = aes(x = year, y = mean, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Annual Abundance",
           color = "Legend") +
      ggtitle("RF Gulf of Maine") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    abund_pa <- ggplot(data = md %>% dplyr::filter(region == "GOM", dataset %in% biomod_dataset), mapping = aes(x = year, y = pa, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = pa - var_pa, ymax = pa + var_pa), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Exceeding of Feeding Threshold",
           color = "Legend") +
      ggtitle("RF Gulf of Maine") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    abund_cpr <- ggplot(data = md %>% dplyr::filter(region == "GOM", dataset == "CPR"), mapping = aes(x = year, y = mean, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Annual Abundance",
           color = "Legend") +
      ggtitle("RF Gulf of Maine") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    abund_cpr_pa <- ggplot(data = md %>% dplyr::filter(region == "GOM", dataset == "CPR"), mapping = aes(x = year, y = pa, color = "Actual")) +
      geom_path() +
      geom_ribbon(aes(ymin = pa - var_pa, ymax = pa + var_pa), alpha=0.2, fill = "red", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "",
           y = "Exceeding of Feeding Threshold",
           color = "Legend") +
      ggtitle("RF Gulf of Maine") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    pred <- ggplot(data = rf_proj_yearly %>% dplyr::filter(region == "GOM"), mapping = aes(x = year, y = mean, color = "Predicted")) +
      geom_path() +
      geom_ribbon(aes(ymin = mean - var, ymax = mean + var), alpha=0.2, fill = "blue", color = NA) +
      scale_x_continuous(breaks = c(2000, 2005, 2010, 2015)) +
      scale_color_manual(values = colors) +
      labs(x = "Year",
           y = "Predicted Probability of Patch",
           color = "Legend") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) 
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'GOM_rf_abund_pred_yearly.png'))
    grid.arrange(abund, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'GOM_rf_pa_pred_yearly.png'))
    grid.arrange(abund_pa, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'GOM_rf_abund_pred_yearly_cpr.png'))
    grid.arrange(abund_cpr, pred)
    dev.off()
    
    png(file.path(fp_out, species, version, "Biomod", "Plots", 'GOM_rf_pa_pred_yearly_cpr.png'))
    grid.arrange(abund_cpr_pa, pred)
    dev.off()
  }
  
  synth_dat <- full_join(md %>% dplyr::filter(dataset %in% biomod_dataset), ensemble_proj_yearly, by = c("region", "year"))
  
  ggplot(data = synth_dat, mapping = aes(x = mean.x, y = mean.y, color = region)) +
    geom_point() +
    geom_errorbar(aes(ymin = mean.y - var.y, ymax = mean.y + var.y), width=.2) +
    geom_errorbarh(aes(xmin = mean.x - var.x, xmax = mean.x + var.x), height = 0.025) +
    labs(y = "Predicted Probability of Patch",
         x = "Inter-annual abundance",
         color = "Region")  +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key = element_rect(color = "transparent", fill = "white")) +
    ggsave(file.path(fp_out, species, version, "Biomod", "Plots", 'ensemble_abund_vs_pred_yearly.png'))
  
  ggplot(data = synth_dat, mapping = aes(x = pa, y = mean.y, color = region)) +
    geom_point() +
    geom_errorbar(aes(ymin = mean.y - var.y, ymax = mean.y + var.y), width=.01) +
    geom_errorbarh(aes(xmin = pa - var_pa, xmax = pa + var_pa), height = 0.025) +
    labs(y = "Predicted Probability of Patch",
         x = "Inter-annual Actual Probability of Patch",
         color = "Region")  +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key = element_rect(color = "transparent", fill = "white")) +
    ggsave(file.path(fp_out, species, version, "Biomod", "Plots", 'ensemble_pa_vs_pred_yearly.png'))
  
  synth_dat <- full_join(md %>% dplyr::filter(dataset == "CPR"), ensemble_proj_yearly, by = c("region", "year"))
  
  ggplot(data = synth_dat, mapping = aes(x = mean.x, y = mean.y, color = region)) +
    geom_point() +
    geom_errorbar(aes(ymin = mean.y - var.y, ymax = mean.y + var.y), width=.2) +
    geom_errorbarh(aes(xmin = mean.x - var.x, xmax = mean.x + var.x), height = 0.025) +
    labs(y = "Predicted Probability of Patch",
         x = "Inter-annual abundance",
         color = "Region")  +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key = element_rect(color = "transparent", fill = "white")) +
    ggsave(file.path(fp_out, species, version, "Biomod", "Plots", 'ensemble_abund_vs_pred_yearly_cpr.png'))
  
  ggplot(data = synth_dat, mapping = aes(x = pa, y = mean.y, color = region)) +
    geom_point() +
    geom_errorbar(aes(ymin = mean.y - var.y, ymax = mean.y + var.y), width=.01) +
    geom_errorbarh(aes(xmin = pa - var_pa, xmax = pa + var_pa), height = 0.025) +
    labs(y = "Predicted Probability of Patch",
         x = "Inter-annual Actual Probability of Patch",
         color = "Region")  +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key = element_rect(color = "transparent", fill = "white")) +
    ggsave(file.path(fp_out, species, version, "Biomod", "Plots", 'ensemble_pa_vs_pred_yearly_cpr.png'))
  
  synth_dat <- full_join(md %>% dplyr::filter(dataset %in% biomod_dataset), gam_proj_yearly, by = c("region", "year"))
  
  ggplot(data = synth_dat, mapping = aes(x = mean.x, y = mean.y, color = region)) +
    geom_point() +
    geom_errorbar(aes(ymin = mean.y - var.y, ymax = mean.y + var.y), width=.2) +
    geom_errorbarh(aes(xmin = mean.x - var.x, xmax = mean.x + var.x), height = 0.025) +
    labs(y = "Predicted Probability of Patch",
         x = "Inter-annual abundance",
         color = "Region")  +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key = element_rect(color = "transparent", fill = "white")) +
    ggsave(file.path(fp_out, species, version, "Biomod", "Plots", 'gam_abund_vs_pred_yearly.png'))
  
  ggplot(data = synth_dat, mapping = aes(x = pa, y = mean.y, color = region)) +
    geom_point() +
    geom_errorbar(aes(ymin = mean.y - var.y, ymax = mean.y + var.y), width=.01) +
    geom_errorbarh(aes(xmin = pa - var_pa, xmax = pa + var_pa), height = 0.025) +
    labs(y = "Predicted Probability of Patch",
         x = "Inter-annual Actual Probability of Patch",
         color = "Region")  +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key = element_rect(color = "transparent", fill = "white")) +
    ggsave(file.path(fp_out, species, version, "Biomod", "Plots", 'gam_pa_vs_pred_yearly.png'))
  
  synth_dat <- full_join(md %>% dplyr::filter(dataset == "CPR"), gam_proj_yearly, by = c("region", "year"))
  
  ggplot(data = synth_dat, mapping = aes(x = mean.x, y = mean.y, color = region)) +
    geom_point() +
    geom_errorbar(aes(ymin = mean.y - var.y, ymax = mean.y + var.y), width=.2) +
    geom_errorbarh(aes(xmin = mean.x - var.x, xmax = mean.x + var.x), height = 0.025) +
    labs(y = "Predicted Probability of Patch",
         x = "Inter-annual abundance",
         color = "Region")  +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key = element_rect(color = "transparent", fill = "white")) +
    ggsave(file.path(fp_out, species, version, "Biomod", "Plots", 'gam_abund_vs_pred_yearly_cpr.png'))
  
  ggplot(data = synth_dat, mapping = aes(x = pa, y = mean.y, color = region)) +
    geom_point() +
    geom_errorbar(aes(ymin = mean.y - var.y, ymax = mean.y + var.y), width=.01) +
    geom_errorbarh(aes(xmin = pa - var_pa, xmax = pa + var_pa), height = 0.025) +
    labs(y = "Predicted Probability of Patch",
         x = "Inter-annual Actual Probability of Patch",
         color = "Region")  +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key = element_rect(color = "transparent", fill = "white")) +
    ggsave(file.path(fp_out, species, version, "Biomod", "Plots", 'gam_pa_vs_pred_yearly_cpr.png'))
  
  synth_dat <- full_join(md %>% dplyr::filter(dataset %in% biomod_dataset), brt_proj_yearly, by = c("region", "year"))
  
  ggplot(data = synth_dat, mapping = aes(x = mean.x, y = mean.y, color = region)) +
    geom_point() +
    geom_errorbar(aes(ymin = mean.y - var.y, ymax = mean.y + var.y), width=.2) +
    geom_errorbarh(aes(xmin = mean.x - var.x, xmax = mean.x + var.x), height = 0.025) +
    labs(y = "Predicted Probability of Patch",
         x = "Inter-annual abundance",
         color = "Region")  +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key = element_rect(color = "transparent", fill = "white")) +
    ggsave(file.path(fp_out, species, version, "Biomod", "Plots", 'brt_abund_vs_pred_yearly.png'))
  
  ggplot(data = synth_dat, mapping = aes(x = pa, y = mean.y, color = region)) +
    geom_point() +
    geom_errorbar(aes(ymin = mean.y - var.y, ymax = mean.y + var.y), width=.01) +
    geom_errorbarh(aes(xmin = pa - var_pa, xmax = pa + var_pa), height = 0.025) +
    labs(y = "Predicted Probability of Patch",
         x = "Inter-annual Actual Probability of Patch",
         color = "Region")  +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key = element_rect(color = "transparent", fill = "white")) +
    ggsave(file.path(fp_out, species, version, "Biomod", "Plots", 'brt_pa_vs_pred_yearly.png'))
  
  synth_dat <- full_join(md %>% dplyr::filter(dataset == "CPR"), brt_proj_yearly, by = c("region", "year"))
  
  ggplot(data = synth_dat, mapping = aes(x = mean.x, y = mean.y, color = region)) +
    geom_point() +
    geom_errorbar(aes(ymin = mean.y - var.y, ymax = mean.y + var.y), width=.2) +
    geom_errorbarh(aes(xmin = mean.x - var.x, xmax = mean.x + var.x), height = 0.025) +
    labs(y = "Predicted Probability of Patch",
         x = "Inter-annual abundance",
         color = "Region")  +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key = element_rect(color = "transparent", fill = "white")) +
    ggsave(file.path(fp_out, species, version, "Biomod", "Plots", 'brt_abund_vs_pred_yearly_cpr.png'))
  
  ggplot(data = synth_dat, mapping = aes(x = pa, y = mean.y, color = region)) +
    geom_point() +
    geom_errorbar(aes(ymin = mean.y - var.y, ymax = mean.y + var.y), width=.01) +
    geom_errorbarh(aes(xmin = pa - var_pa, xmax = pa + var_pa), height = 0.025) +
    labs(y = "Predicted Probability of Patch",
         x = "Inter-annual Actual Probability of Patch",
         color = "Region")  +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key = element_rect(color = "transparent", fill = "white")) +
    ggsave(file.path(fp_out, species, version, "Biomod", "Plots", 'brt_pa_vs_pred_yearly_cpr.png'))
  
  synth_dat <- full_join(md %>% dplyr::filter(dataset %in% biomod_dataset), rf_proj_yearly, by = c("region", "year"))
  
  ggplot(data = synth_dat, mapping = aes(x = mean.x, y = mean.y, color = region)) +
    geom_point() +
    geom_errorbar(aes(ymin = mean.y - var.y, ymax = mean.y + var.y), width=.2) +
    geom_errorbarh(aes(xmin = mean.x - var.x, xmax = mean.x + var.x), height = 0.025) +
    labs(y = "Predicted Probability of Patch",
         x = "Inter-annual abundance",
         color = "Region") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key = element_rect(color = "transparent", fill = "white")) +
    ggsave(file.path(fp_out, species, version, "Biomod", "Plots", 'rf_abund_vs_pred_yearly.png'))
  
  ggplot(data = synth_dat, mapping = aes(x = pa, y = mean.y, color = region)) +
    geom_point() +
    geom_errorbar(aes(ymin = mean.y - var.y, ymax = mean.y + var.y), width=.01) +
    geom_errorbarh(aes(xmin = pa - var_pa, xmax = pa + var_pa), height = 0.025) +
    labs(y = "Predicted Probability of Patch",
         x = "Inter-annual Actual Probability of Patch",
         color = "Region") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key = element_rect(color = "transparent", fill = "white")) +
    ggsave(file.path(fp_out, species, version, "Biomod", "Plots", 'rf_pa_vs_pred_yearly.png'))
  
  synth_dat <- full_join(md %>% dplyr::filter(dataset == "CPR"), rf_proj_yearly, by = c("region", "year"))
  
  ggplot(data = synth_dat, mapping = aes(x = mean.x, y = mean.y, color = region)) +
    geom_point() +
    geom_errorbar(aes(ymin = mean.y - var.y, ymax = mean.y + var.y), width=.2) +
    geom_errorbarh(aes(xmin = mean.x - var.x, xmax = mean.x + var.x), height = 0.025) +
    labs(y = "Predicted Probability of Patch",
         x = "Inter-annual abundance",
         color = "Region") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key = element_rect(color = "transparent", fill = "white")) +
    ggsave(file.path(fp_out, species, version, "Biomod", "Plots", 'rf_abund_vs_pred_yearly_cpr.png'))
  
  ggplot(data = synth_dat, mapping = aes(x = pa, y = mean.y, color = region)) +
    geom_point() +
    geom_errorbar(aes(ymin = mean.y - var.y, ymax = mean.y + var.y), width=.01) +
    geom_errorbarh(aes(xmin = pa - var_pa, xmax = pa + var_pa), height = 0.025) +
    labs(y = "Predicted Probability of Patch",
         x = "Inter-annual Actual Probability of Patch",
         color = "Region") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key = element_rect(color = "transparent", fill = "white")) +
    ggsave(file.path(fp_out, species, version, "Biomod", "Plots", 'rf_pa_vs_pred_yearly_cpr.png'))
  
}
  
