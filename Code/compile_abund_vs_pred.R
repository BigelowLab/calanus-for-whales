# Camille Ross
# 20 July, 2020
# Purpose: Compile and plot actual vs. predicted values for multiple years of model runs

# -------- Load libraries --------
require(plyr)
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
#'@param species <chr> species to model; choices are "cfin", "ctyp", or "pseudo"
compile_abund_vs_pred <- function(version, fp_out, threshold, years, species = "cfin") {
  
  # -------- Create output directories --------
  dir.create(fp_out, showWarnings = FALSE) 
  dir.create(file.path(fp_out, species, version), showWarnings = FALSE)
  dir.create(file.path(fp_out, species, version, "Climatologies"), showWarnings = FALSE)
  # -------- Create directory for evaluations --------
  dir.create(file.path(fp_out, species, version, "Climatologies", "Evals"), showWarnings = FALSE)
  
  # -------- Test if projection exists --------
  for (year in years) {
    for (month in 1:12) {
      if (paste0("abund_vs_pred_", year, "_", month, ".csv") %in% list.files(file.path(fp_out, species, version, "GAMs", "Projections")) &
          paste0("abund_vs_pred_", year, "_", month, ".csv") %in% list.files(file.path(fp_out, species, version, "BRTs", "Projections"))) {
        if (year == 2000 & month == 1) {
          gam_abunds <- readr::read_csv(file.path(fp_out, species, version, "GAMs", "Projections", "abund_vs_pred_2000_1.csv"))
          gam_full_data <- data.frame("year" = year, "month" = month, "abund" = gam_abunds$abund, "pred" = gam_abunds$pred)
          
          brt_abunds <- readr::read_csv(file.path(fp_out, species, version, "BRTs", "Projections", "abund_vs_pred_2000_1.csv"))
          brt_full_data <- data.frame("year" = year, "month" = month, "abund" = brt_abunds$abund, "pred" = brt_abunds$pred)
          
          if (paste0("ensemble_abund_vs_pred_", year, "_", month, ".csv") %in% list.files(file.path(fp_out, species, version, "Biomod", "Projections"))) {
            biomod_abunds <- readr::read_csv(file.path(fp_out, species, version, "Biomod", "Projections", "ensemble_abund_vs_pred_2000_1.csv"))
            biomod_full_data <- data.frame("year" = year, "month" = month, "abund" = biomod_abunds$abund, "pred" = biomod_abunds$ensemble_pred)
          }
        } else {
          gam_abunds <- readr::read_csv(file.path(fp_out, species, version, "GAMs", "Projections", paste0("abund_vs_pred_", year, "_", month, ".csv")))
          
          gam_temp <- data.frame("year" = year, "month" = month, "abund" = gam_abunds$abund, "pred" = gam_abunds$pred)
          gam_full_data <- rbind(gam_full_data, gam_temp)
          
          brt_abunds <- readr::read_csv(file.path(fp_out, species, version, "BRTs", "Projections", paste0("abund_vs_pred_", year, "_", month, ".csv")))
          
          brt_temp <- data.frame("year" = year, "month" = month, "abund" = brt_abunds$abund, "pred" = brt_abunds$pred)
          brt_full_data <- rbind(brt_full_data, brt_temp)
          
          if (paste0("ensemble_abund_vs_pred_", year, "_", month, ".csv") %in% list.files(file.path(fp_out, species, version, "Biomod", "Projections"))) {
            biomod_abunds <- readr::read_csv(file.path(fp_out, species, version, "Biomod", "Projections", paste0("ensemble_abund_vs_pred_", year, "_", month, ".csv")))
            
            biomod_temp <- data.frame("year" = year, "month" = month, "abund" = biomod_abunds$abund, "pred" = biomod_abunds$ensemble_pred)
            biomod_full_data <- rbind(biomod_full_data, biomod_temp)
          }
        }
      }
    }
  }
  
  # -------- Write data to CSV --------
  gam_full_data <- gam_full_data %>% 
    readr::write_csv(file.path(fp_out, species, version, "Climatologies", "Projections", "gam_compiled_abund_vs_pred.csv")) %>%
    count(abund, pred)
  
  brt_full_data <- brt_full_data %>% 
    readr::write_csv(file.path(fp_out, species, version, "Climatologies", "Projections", "brt_compiled_abund_vs_pred.csv")) %>%
    count(abund, pred)
  
  biomod_full_data <- biomod_full_data %>% 
    readr::write_csv(file.path(fp_out, species, version, "Biomod_Climatologies", "Projections", "biomod_compiled_abund_vs_pred.csv")) %>%
    dplyr::count(abund, pred)
  
  # -------- Plot actual vs. predicted --------
  for (i in 1:12) {
    # ---- GAMs ----
    ggplot(gam_full_data %>% dplyr::filter(month == i), aes(x = abund, y = pred, color = as.factor(year)), alpha = 0.01) +
      geom_point() +
      scale_color_viridis(name = "Year", discrete=TRUE) +
      ylab("Predicted value") +
      xlab("Actual Value") +
      ylim(c(0,5)) +
      xlim(c(0,5)) +
      guides(color=guide_legend(ncol=2)) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black")) +
      ggsave(filename = file.path(fp_out, species, version, "Climatologies", "Plots", paste0("gam_abund_vs_pred_", i, ".png")))
  
    # ---- BRTs ----
    ggplot(brt_full_data %>% dplyr::filter(month == i), aes(x = abund, y = pred, color = as.factor(year)), alpha = 0.01) +
      geom_point() +
      scale_color_viridis(name = "Year", discrete=TRUE) +
      ylab("Predicted value") +
      xlab("Actual Value") +
      ylim(c(0,5)) +
      xlim(c(0,5)) +
      guides(color=guide_legend(ncol=2)) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black")) +
      ggsave(filename = file.path(fp_out, species, version, "Climatologies", "Plots", paste0("brt_abund_vs_pred_", i, ".png")))
    
    # ---- Biomod ensemble ----
    ggplot(biomod_full_data %>% dplyr::filter(month == i), aes(x = log10(abund + 1), y = pred, color = as.factor(count)), alpha = 0.01) +
      geom_point() +
      scale_color_viridis(name = "Year", discrete=TRUE) +
      ylab("Predicted value") +
      xlab("Actual Value") +
      ylim(c(0,1)) +
      xlim(c(0,5)) +
      geom_vline(xintercept = log10(threshold + 1)) +
      guides(color=guide_legend(ncol=2)) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black")) +
      ggsave(filename = file.path(fp_out, species, version, "Biomod_Climatologies", "Plots", paste0("biomod_abund_vs_pred_", i, ".png")))
  }
}


