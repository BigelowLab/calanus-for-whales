# Camille Ross
# 10 August, 2020
# Purpose: Divide model results into three regions and plot against observed abundance data.

# -------- Load libraries --------

require(lubridate)
require(readr)
require(dplyr)

# -------- Main function --------
#'@param version <chr> version of model
#'@param fp_out <chr> file path save the data to 
#'@param species <chr> species to model; choices are "cfin", "ctyp", or "pseudo"
plot_regions <- function(version, fp_out, datasets, species = "cfin") {
  
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
  } else if (species == "pseudo") {
    md <- md %>% dplyr::group_by(dataset) %>%
      dplyr::mutate(mean = mean(log10(`pseudo_total` + 1), na.rm = TRUE),
                    sd = sd(log10(`pseudo_total` + 1), na.rm = TRUE),
                    anomaly = (log10(`pseudo_total` + 1) - mean) / sd) %>%
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
    } else if (species == "pseudo") {
      md$abund <- as.data.frame(log10(md[paste0(species, "_total")] + 1))$pseudo_total
    }
  }
  
  # -------- Exclude NAs and select columns --------
  md <- md %>% dplyr::select(lat, lon, year, month, abund, wind, fetch, chl, int_chl, bots, bott, sss, sst, lag_sst, uv, bat, dist, slope) %>%
    as.data.frame() %>%
    na.exclude() %>%
    dplyr::filter(!is.infinite(abund)) %>%
    dplyr::mutate(season = if_else(month %in% c(1:3), 1,
                                   if_else(month %in% c(4:6), 2,
                                           if_else(month %in% c(7:9), 3, 4))),
                  region = if_else(lat <= 41.5 & lon < -70, "MAB", 
                                   if_else(lat >= 39 & lat <= 42 & lon >= -70 & lon <= -68, "GBK", "GOM"))) %>%
    dplyr::group_by(region, month) %>%
    dplyr::summarize(mean = mean(abund, na.rm=TRUE),
                     stdev = sd(abund, na.rm = TRUE))
  
  for (i in 1:12) {
    # -------- Test if projection exists --------
    if (paste0("gam_proj_", i, ".tif") %in% list.files(file.path(fp_out, species, version, "Climatologies", "Projections")) &
        paste0("brt_proj_", i, ".tif") %in% list.files(file.path(fp_out, species, version, "Climatologies", "Projections"))) {
      if (i == 1) {
        gam_proj_df <- raster::raster(file.path(fp_out, species, version, "Climatologies", "Projections", paste0("gam_proj_", i, ".tif"))) %>%
          as.data.frame(xy = TRUE) %>%
          dplyr::mutate(month = i)
        
        names(gam_proj_df) <- c("x", "y", "proj", "month")
        
        brt_proj_df <- raster::raster(file.path(fp_out, species, version, "Climatologies", "Projections", paste0("brt_proj_", i, ".tif"))) %>%
          as.data.frame(xy = TRUE) %>%
          dplyr::mutate(month = i)
        
        names(brt_proj_df) <- c("x", "y", "proj", "month")
        
      } else {
        gam_proj <- raster::raster(file.path(fp_out, species, version, "Climatologies", "Projections", paste0("gam_proj_", i, ".tif"))) %>%
          as.data.frame(xy = TRUE) %>%
          dplyr::mutate(month = i)
        
        names(gam_proj) <- c("x", "y", "proj", "month")
        
        gam_proj_df <- rbind(gam_proj_df, gam_proj)
        
        brt_proj <- raster::raster(file.path(fp_out, species, version, "Climatologies", "Projections", paste0("brt_proj_", i, ".tif"))) %>%
          as.data.frame(xy = TRUE) %>%
          dplyr::mutate(month = i)
        
        names(brt_proj) <- c("x", "y", "proj", "month")
        
        brt_proj_df <- rbind(brt_proj_df, brt_proj)
      }
    }
  }
  
  # -------- Compute regions --------
  gam_proj_df <- gam_proj_df %>%
    dplyr::mutate(region = if_else(y <= 41 & x < -70, "MAB", 
                                   if_else(y >= 40 & x <= 42 & y >= -70 & x <= -68, "GBK", "GOM"))) %>%
    dplyr::group_by(region, month) %>%
    dplyr::summarize(mean = mean(proj, na.rm = TRUE),
                     stdev = sd(proj, na.rm = TRUE))
  
  brt_proj_df <- brt_proj_df %>%
    dplyr::mutate(region = if_else(y <= 41 & x < -70, "MAB", 
                                   if_else(y >= 40 & x <= 42 & y >= -70 & x <= -68, "GBK", "GOM"))) %>%
    dplyr::group_by(region, month) %>%
    dplyr::summarize(mean = mean(proj, na.rm = TRUE),
                     stdev = sd(proj, na.rm = TRUE))
  
  # -------- Initialize legend colors --------
  colors <- c("Actual" = "red", "Predicted" = "blue")

  # -------- Plot MAB --------
  if ("MAB" %in% unique(md$region)) {
    # ---- Plot GAMs ----
    ggplot(data = md %>% dplyr::filter(region == "MAB"), mapping = aes(x = month, y = mean, color = "Actual")) +
      geom_smooth(fill = "red") +
      geom_smooth(data = gam_proj_df %>% dplyr::filter(region == "MAB"), mapping = aes(x = month, y = mean, color = "Predicted"), fill = "blue") +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "Month",
           y = "Climatological Abundance",
           color = "Legend") +
      ggtitle("Mid Atlantic Bight") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) +
      ggsave(filename = "MAB_gam_abund_pred_climatological.png", path = file.path(fp_out, species, version, "Climatologies", "Plots"))
    
    # ---- Plot BRTs ----
    ggplot(data = md %>% dplyr::filter(region == "MAB"), mapping = aes(x = month, y = mean, color = "Actual")) +
      geom_smooth(fill = "red") +
      geom_smooth(data = brt_proj_df %>% dplyr::filter(region == "MAB"), mapping = aes(x = month, y = mean, color = "Predicted"), fill = "blue") +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) +
      labs(x = "Month",
           y = "Climatological Abundance",
           color = "Legend") +
      ggtitle("Mid Atlantic Bight") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) +
      ggsave(filename = "MAB_brt_abund_pred_climatological.png", path = file.path(fp_out, species, version, "Climatologies", "Plots"))
    
  }
  
  # -------- Plot GBK --------
  if ("GBK" %in% unique(md$region)) {
    # ---- Plot GAMs ----
    ggplot(data = md %>% dplyr::filter(region == "GBK"), mapping = aes(x = month, y = mean, color = "Actual")) +
      geom_smooth(fill = "red") +
      geom_smooth(data = gam_proj_df %>% dplyr::filter(region == "GBK"), mapping = aes(x = month, y = mean, color = "Predicted"), fill = "blue") +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) + 
      labs(x = "Month",
           y = "Climatological Abundance",
           color = "Legend") +
      ggtitle("George's Bank") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) +
      ggsave(filename = "GBK_gam_abund_pred_climatological.png", path = file.path(fp_out, species, version, "Climatologies", "Plots"))
    
    # ---- Plot BRTs ----
    ggplot(data = md %>% dplyr::filter(region == "GBK"), mapping = aes(x = month, y = mean, color = "Actual")) +
      geom_smooth(fill = "red") +
      geom_smooth(data = brt_proj_df %>% dplyr::filter(region == "GBK"), mapping = aes(x = month, y = mean, color = "Predicted"), fill = "blue") +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) + 
      labs(x = "Month",
           y = "Climatological Abundance",
           color = "Legend") +
      ggtitle("George's Bank") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) +
      ggsave(filename = "GBK_brt_abund_pred_climatological.png", path = file.path(fp_out, species, version, "Climatologies", "Plots"))
    
  }
  
  # -------- Plot GOM --------
  if ("GOM" %in% unique(md$region)) {
    # ---- Plot GAMs ----
    ggplot(data = md %>% dplyr::filter(region == "GOM"), mapping = aes(x = month, y = mean, color = "Actual")) +
      geom_smooth(fill = "red") +
      geom_smooth(data = gam_proj_df %>% dplyr::filter(region == "GOM"), mapping = aes(x = month, y = mean, color = "Predicted"), fill = "blue") +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) + 
      labs(x = "Month",
           y = "Climatological Abundance",
           color = "Legend") +
      ggtitle("Gulf of Maine") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) +
      ggsave(filename = "GOM_gam_abund_pred_climatological.png", path = file.path(fp_out, species, version, "Climatologies", "Plots"))
    
    # ---- Plot BRTs ----
    ggplot(data = md %>% dplyr::filter(region == "GOM"), mapping = aes(x = month, y = mean, color = "Actual")) +
      geom_smooth(fill = "red") +
      geom_smooth(data = brt_proj_df %>% dplyr::filter(region == "GOM"), mapping = aes(x = month, y = mean, color = "Predicted"), fill = "blue") +
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      scale_color_manual(values = colors) + 
      labs(x = "Month",
           y = "Climatological Abundance",
           color = "Legend") +
      ggtitle("Gulf of Maine") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) +
      ggsave(filename = "GOM_brt_abund_pred_climatological.png", path = file.path(fp_out, species, version, "Climatologies", "Plots"))
    
  }
  
  # -------- ANNUAL AVERAGE --------
  
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
  } else if (species == "pseudo") {
    md <- md %>% dplyr::group_by(dataset) %>%
      dplyr::mutate(mean = mean(log10(`pseudo_total` + 1), na.rm = TRUE),
                    sd = sd(log10(`pseudo_total` + 1), na.rm = TRUE),
                    anomaly = (log10(`pseudo_total` + 1) - mean) / sd) %>%
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
    } else if (species == "pseudo") {
      md$abund <- as.data.frame(log10(md[paste0(species, "_total")] + 1))$pseudo_total
    }
  }
  
  # -------- Exclude NAs and select columns --------
  md <- md %>% dplyr::select(lat, lon, year, month, abund, wind, fetch, chl, int_chl, bots, bott, sss, sst, lag_sst, uv, bat, dist, slope) %>%
    as.data.frame() %>%
    na.exclude() %>%
    dplyr::filter(!is.infinite(abund)) %>%
    dplyr::mutate(season = if_else(month %in% c(1:3), 1,
                                   if_else(month %in% c(4:6), 2,
                                           if_else(month %in% c(7:9), 3, 4))),
                  region = if_else(lat <= 41.5 & lon < -70, "MAB", 
                                   if_else(lat >= 39 & lat <= 42 & lon >= -70 & lon <= -68, "GBK", "GOM"))) %>%
    dplyr::group_by(region, year) %>%
    dplyr::summarize(mean = mean(abund, na.rm=TRUE),
                     stdev = sd(abund, na.rm = TRUE))
  
  for (i in 1:12) {
    for (j in 2000:2017) {
      # -------- Test if projection exists --------
      if (paste0("proj_", j, "_", i, ".tif") %in% list.files(file.path(fp_out, species, version, "GAMs", "Projections")) &
          paste0("proj_", j, "_", i, ".tif") %in% list.files(file.path(fp_out, species, version, "BRTs", "Projections"))) {
        if (i == 1 & j == 2000) {
          gam_proj_df <- raster::raster(file.path(fp_out, species, version, "GAMs", "Projections", paste0("proj_", j, "_", i, ".tif"))) %>%
            as.data.frame(xy = TRUE) %>%
            dplyr::mutate(month = i,
                          year = j)
          
          names(gam_proj_df) <- c("x", "y", "proj", "month", "year")
          
          brt_proj_df <- raster::raster(file.path(fp_out, species, version, "BRTs", "Projections", paste0("proj_", j, "_", i, ".tif"))) %>%
            as.data.frame(xy = TRUE) %>%
            dplyr::mutate(month = i,
                          year = j)
          
          names(brt_proj_df) <- c("x", "y", "proj", "month", "year")
          
        } else {
          gam_proj <- raster::raster(file.path(fp_out, species, version, "GAMs", "Projections", paste0("proj_", j, "_", i, ".tif"))) %>%
            as.data.frame(xy = TRUE) %>%
            dplyr::mutate(month = i,
                          year = j)
          
          names(gam_proj) <- c("x", "y", "proj", "month", "year")
          
          gam_proj_df <- rbind(gam_proj_df, gam_proj)
          
          brt_proj <- raster::raster(file.path(fp_out, species, version, "BRTs", "Projections", paste0("proj_", j, "_", i, ".tif"))) %>%
            as.data.frame(xy = TRUE) %>%
            dplyr::mutate(month = i,
                          year = j)
          
          names(brt_proj) <- c("x", "y", "proj", "month", "year")
          
          brt_proj_df <- rbind(brt_proj_df, brt_proj)
        }
      }
    }
  }
  
  # -------- Compute regions --------
  gam_proj_df <- gam_proj_df %>%
    dplyr::mutate(region = if_else(y <= 41 & x < -70, "MAB", 
                                   if_else(y >= 40 & x <= 42 & y >= -70 & x <= -68, "GBK", "GOM"))) %>%
    dplyr::group_by(region, year) %>%
    dplyr::summarize(mean = mean(proj, na.rm = TRUE),
                     stdev = sd(proj, na.rm = TRUE))
  
  brt_proj_df <- brt_proj_df %>%
    dplyr::mutate(region = if_else(y <= 41 & x < -70, "MAB", 
                                   if_else(y >= 40 & x <= 42 & y >= -70 & x <= -68, "GBK", "GOM"))) %>%
    dplyr::group_by(region, year) %>%
    dplyr::summarize(mean = mean(proj, na.rm = TRUE),
                     stdev = sd(proj, na.rm = TRUE))
  
  # -------- Plot MAB --------
  if ("MAB" %in% unique(md$region)) {
    # ---- GAMs ----
    ggplot(data = md %>% dplyr::filter(region == "MAB"), mapping = aes(x = year, y = mean, color = "Actual")) +
      geom_smooth(fill = "red") +
      geom_smooth(data = gam_proj_df %>% dplyr::filter(region == "MAB"), mapping = aes(x = year, y = mean, color = "Predicted"), fill = "blue") +
      scale_x_continuous(breaks = c(2000, 2004, 2008, 2012, 2016)) +
      scale_color_manual(values = colors) +
      labs(x = "Year",
           y = "Annual Average Abundance",
           color = "")+
      ggtitle("Mid Atlantic Bight") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) +
      ggsave(filename = "MAB_gam_abund_pred_annual.png", path = file.path(fp_out, species, version, "GAMs", "Plots"))
    
    # ---- BRTs ----
    ggplot(data = md %>% dplyr::filter(region == "MAB"), mapping = aes(x = year, y = mean, color = "Actual")) +
      geom_smooth(fill = "red") +
      geom_smooth(data = brt_proj_df %>% dplyr::filter(region == "MAB"), mapping = aes(x = year, y = mean, color = "Predicted"), fill = "blue") +
      scale_x_continuous(breaks = c(2000, 2004, 2008, 2012, 2016)) +
      scale_color_manual(values = colors) +
      labs(x = "Year",
           y = "Annual Average Abundance",
           color = "") +
      ggtitle("Mid Atlantic Bight") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) +
      ggsave(filename = "MAB_brt_abund_pred_annual.png", path = file.path(fp_out, species, version, "BRTs", "Plots"))
    
  }
  
  # -------- Plot GBK --------
  if ("GBK" %in% unique(md$region)) {
    # ---- GAMs ----
    ggplot(data = md %>% dplyr::filter(region == "GBK"), mapping = aes(x = year, y = mean, color = "Actual")) +
      geom_smooth(fill = "red") +
      geom_smooth(data = gam_proj_df %>% dplyr::filter(region == "GBK"), mapping = aes(x = year, y = mean, color = "Predicted"), fill = "blue") +
      scale_x_continuous(breaks = c(2000, 2004, 2008, 2012, 2016)) +
      scale_color_manual(values = colors) + 
      labs(x = "Year",
           y = "Annual Average Abundance",
           color = "") +
      ggtitle("George's Bank") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) +
      ggsave(filename = "GBK_gam_abund_pred_annual.png", path = file.path(fp_out, species, version, "GAMs", "Plots"))
    
    # ---- BRTs ----
    ggplot(data = md %>% dplyr::filter(region == "GBK"), mapping = aes(x = year, y = mean, color = "Actual")) +
      geom_smooth(fill = "red") +
      geom_smooth(data = brt_proj_df %>% dplyr::filter(region == "GBK"), mapping = aes(x = year, y = mean, color = "Predicted"), fill = "blue") +
      scale_x_continuous(breaks = c(2000, 2004, 2008, 2012, 2016)) +
      scale_color_manual(values = colors) + 
      labs(x = "Year",
           y = "Annual Average Abundance",
           color = "")+
      ggtitle("George's Bank") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) +
      ggsave(filename = "GBK_brt_abund_pred_annual.png", path = file.path(fp_out, species, version, "BRTs", "Plots"))
    
  }
  
  # -------- Plot GOM --------
  if ("GOM" %in% unique(md$region)) {
    # ---- GAMs ----
    ggplot(data = md %>% dplyr::filter(region == "GOM"), mapping = aes(x = year, y = mean, color = "Actual")) +
      geom_smooth(fill = "red") +
      geom_smooth(data = gam_proj_df %>% dplyr::filter(region == "GOM"), mapping = aes(x = year, y = mean, color = "Predicted"), fill = "blue") +
      scale_x_continuous(breaks = c(2000, 2004, 2008, 2012, 2016)) +
      scale_color_manual(values = colors) + 
      labs(x = "Year",
           y = "Annual Average Abundance",
           color = "") +
      ggtitle("Gulf of Maine") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) +
      ggsave(filename = "GOM_gam_abund_pred_annual.png", path = file.path(fp_out, species, version, "GAMs", "Plots"))
    
    # ---- BRTs ----
    ggplot(data = md %>% dplyr::filter(region == "GOM"), mapping = aes(x = year, y = mean, color = "Actual")) +
      geom_smooth(fill = "red") +
      geom_smooth(data = brt_proj_df %>% dplyr::filter(region == "GOM"), mapping = aes(x = year, y = mean, color = "Predicted"), fill = "blue") +
      scale_x_continuous(breaks = c(2000, 2004, 2008, 2012, 2016)) +
      scale_color_manual(values = colors) + 
      labs(x = "Year",
           y = "Annual Average Abundance",
           color = "") +
      ggtitle("Gulf of Maine") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) +
      ggsave(filename = "GOM_brt_abund_pred_annual.png", path = file.path(fp_out, species, version, "BRTs", "Plots"))
    
  }
  
  # -------- MONTHLY TIME SERIES --------
  
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
  } else if (species == "pseudo") {
    md <- md %>% dplyr::group_by(dataset) %>%
      dplyr::mutate(mean = mean(log10(`pseudo_total` + 1), na.rm = TRUE),
                    sd = sd(log10(`pseudo_total` + 1), na.rm = TRUE),
                    anomaly = (log10(`pseudo_total` + 1) - mean) / sd) %>%
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
    } else if (species == "pseudo") {
      md$abund <- as.data.frame(log10(md[paste0(species, "_total")] + 1))$pseudo_total
    }
  }
  
  # -------- Exclude NAs and select columns --------
  md <- md %>% dplyr::select(lat, lon, year, month, abund, wind, fetch, chl, int_chl, bots, bott, sss, sst, lag_sst, uv, bat, dist, slope) %>%
    as.data.frame() %>%
    na.exclude() %>%
    dplyr::filter(!is.infinite(abund)) %>%
    dplyr::mutate(season = if_else(month %in% c(1:3), 1,
                                   if_else(month %in% c(4:6), 2,
                                           if_else(month %in% c(7:9), 3, 4))),
                  region = if_else(lat <= 41.5 & lon < -70, "MAB", 
                                   if_else(lat >= 39 & lat <= 42 & lon >= -70 & lon <= -68, "GBK", "GOM"))) %>%
    dplyr::group_by(region, year, month) %>%
    dplyr::summarize(mean = mean(abund, na.rm=TRUE),
                     stdev = sd(abund, na.rm = TRUE)) %>%
    dplyr::mutate(date = zoo::as.Date(zoo::as.yearmon(paste(month, year), "%m %Y")))
  
  # -------- Load projections -------
  for (i in 1:12) {
    for (j in years) {
      # -------- Test if projection exists --------
      if (paste0("proj_", j, "_", i, ".tif") %in% list.files(file.path(fp_out, species, version, "GAMs", "Projections")) &
          paste0("proj_", j, "_", i, ".tif") %in% list.files(file.path(fp_out, species, version, "BRTs", "Projections"))) {
        if (i == 1 & j == 2000) {
          gam_proj_df <- raster::raster(file.path(fp_out, species, version, "GAMs", "Projections", paste0("proj_", j, "_", i, ".tif"))) %>%
            as.data.frame(xy = TRUE) %>%
            dplyr::mutate(month = i,
                          year = j)
          
          names(gam_proj_df) <- c("x", "y", "proj", "month", "year")
          
          brt_proj_df <- raster::raster(file.path(fp_out, species, version, "BRTs", "Projections", paste0("proj_", j, "_", i, ".tif"))) %>%
            as.data.frame(xy = TRUE) %>%
            dplyr::mutate(month = i,
                          year = j)
          
          names(brt_proj_df) <- c("x", "y", "proj", "month", "year")
          
        } else {
          gam_proj <- raster::raster(file.path(fp_out, species, version, "GAMs", "Projections", paste0("proj_", j, "_", i, ".tif"))) %>%
            as.data.frame(xy = TRUE) %>%
            dplyr::mutate(month = i,
                          year = j)
          
          names(gam_proj) <- c("x", "y", "proj", "month", "year")
          
          gam_proj_df <- rbind(gam_proj_df, gam_proj)
          
          brt_proj <- raster::raster(file.path(fp_out, species, version, "BRTs", "Projections", paste0("proj_", j, "_", i, ".tif"))) %>%
            as.data.frame(xy = TRUE) %>%
            dplyr::mutate(month = i,
                          year = j)
          
          names(brt_proj) <- c("x", "y", "proj", "month", "year")
          
          brt_proj_df <- rbind(brt_proj_df, brt_proj)
        }
      }
    }
  }
  
  # -------- Compute regions --------
  gam_proj_df <- gam_proj_df %>%
    dplyr::mutate(region = if_else(y <= 41 & x < -70, "MAB", 
                                   if_else(y >= 40 & x <= 42 & y >= -70 & x <= -68, "GBK", "GOM"))) %>%
    dplyr::group_by(region, year, month) %>%
    # Compute average and standard deviation
    dplyr::summarize(mean = mean(proj, na.rm = TRUE),
                     stdev = sd(proj, na.rm = TRUE)) %>%
    # Add date
    dplyr::mutate(date = zoo::as.Date(zoo::as.yearmon(paste(month, year), "%m %Y")))
  
  brt_proj_df <- brt_proj_df %>%
    dplyr::mutate(region = if_else(y <= 41 & x < -70, "MAB", 
                                   if_else(y >= 40 & x <= 42 & y >= -70 & x <= -68, "GBK", "GOM"))) %>%
    dplyr::group_by(region, year, month) %>%
    # Compute average and standard deviation
    dplyr::summarize(mean = mean(proj, na.rm = TRUE),
                     stdev = sd(proj, na.rm = TRUE)) %>%
    # Add date
    dplyr::mutate(date = zoo::as.Date(zoo::as.yearmon(paste(month, year), "%m %Y")))
  
  # -------- Plot MAB --------
  if ("MAB" %in% unique(md$region)) {
    # ---- GAMs ----
    ggplot(data = md %>% dplyr::filter(region == "MAB"), mapping = aes(x = date, y = mean, color = "Actual")) +
      geom_smooth(fill = "red") +
      geom_smooth(data = gam_proj_df %>% dplyr::filter(region == "MAB"), mapping = aes(x = date, y = mean, color = "Predicted"), fill = "blue") +
      scale_color_manual(values = colors) +
      labs(x = "Date",
           y = "Average Abundance",
           color = "")+
      ggtitle("Mid Atlantic Bight") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) +
      ggsave(filename = "MAB_gam_abund_pred_monthly.png", path = file.path(fp_out, species, version, "GAMs", "Plots"))
    
    # ---- BRTs ----
    ggplot(data = md %>% dplyr::filter(region == "MAB"), mapping = aes(x = date, y = mean, color = "Actual")) +
      geom_smooth(fill = "red") +
      geom_smooth(data = brt_proj_df %>% dplyr::filter(region == "MAB"), mapping = aes(x = date, y = mean, color = "Predicted"), fill = "blue") +
      scale_color_manual(values = colors) +
      labs(x = "Date",
           y = "Average Abundance",
           color = "") +
      ggtitle("Mid Atlantic Bight") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) +
      ggsave(filename = "MAB_brt_abund_pred_monthly.png", path = file.path(fp_out, species, version, "BRTs", "Plots"))
  
  }
  
  # -------- Plot GBK --------
  if ("GBK" %in% unique(md$region)) {
    # ---- GAMs ----
    ggplot(data = md %>% dplyr::filter(region == "GBK"), mapping = aes(x = date, y = mean, color = "Actual")) +
      geom_smooth(fill = "red") +
      geom_smooth(data = gam_proj_df %>% dplyr::filter(region == "GBK"), mapping = aes(x = date, y = mean, color = "Predicted"), fill = "blue") +
      scale_color_manual(values = colors) + 
      labs(x = "Date",
           y = "Average Abundance",
           color = "") +
      ggtitle("George's Bank") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) +
      ggsave(filename = "GBK_gam_abund_pred_monthly.png", path = file.path(fp_out, species, version, "GAMs", "Plots"))
    
    # ---- BRTs ----
    ggplot(data = md %>% dplyr::filter(region == "GBK"), mapping = aes(x = date, y = mean, color = "Actual")) +
      geom_smooth(fill = "red") +
      geom_smooth(data = brt_proj_df %>% dplyr::filter(region == "GBK"), mapping = aes(x = date, y = mean, color = "Predicted"), fill = "blue") +
      scale_color_manual(values = colors) + 
      labs(x = "Date",
           y = "Average Abundance",
           color = "")+
      ggtitle("George's Bank") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) +
      ggsave(filename = "GBK_brt_abund_pred_monthly.png", path = file.path(fp_out, species, version, "BRTs", "Plots"))
    
  }
  
  # -------- Plot GOM --------
  if ("GOM" %in% unique(md$region)) {
    # ---- GAMs ----
    ggplot(data = md %>% dplyr::filter(region == "GOM"), mapping = aes(x = date, y = mean, color = "Actual")) +
      geom_smooth(fill = "red") +
      geom_smooth(data = gam_proj_df %>% dplyr::filter(region == "GOM"), mapping = aes(x = date, y = mean, color = "Predicted"), fill = "blue") +
      scale_color_manual(values = colors) + 
      labs(x = "Date",
           y = "Average Abundance",
           color = "") +
      ggtitle("Gulf of Maine") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) +
      ggsave(filename = "GOM_gam_abund_pred_monthly.png", path = file.path(fp_out, species, version, "GAMs", "Plots"))
    
    # ---- BRTs ----
    ggplot(data = md %>% dplyr::filter(region == "GOM"), mapping = aes(x = date, y = mean, color = "Actual")) +
      geom_smooth(fill = "red") +
      geom_smooth(data = brt_proj_df %>% dplyr::filter(region == "GOM"), mapping = aes(x = date, y = mean, color = "Predicted"), fill = "blue") +
      scale_color_manual(values = colors) + 
      labs(x = "Date",
           y = "Average Abundance",
           color = "") +
      ggtitle("Gulf of Maine") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(color = "transparent", fill = "white")) +
      ggsave(filename = "GOM_brt_abund_pred_monthly.png", path = file.path(fp_out, species, version, "BRTs", "Plots"))
    
  }
  
}
  
