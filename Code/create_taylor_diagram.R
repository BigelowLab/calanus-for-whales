# Camille Ross
# 24 Nov, 2020
# Purpose: Compile and plot actual vs. predicted values for multiple years of model runs

# -------- Load libraries --------
require(plyr)
require(dplyr)
require(readr)
require(gbm)
require(viridis)
require(Metrics)
require(mapdata)
require(maps)
require(ggplot2)
require(stats)
require(pals)
require(plotrix)


# -------- Main function --------
#'@param version <chr> version of model
#'@param fp_out <chr> file path save the data to 
#'@param species <chr> species to model; choices are "cfin", "ctyp", or "pcal"
create_taylor_diagram <- function(version, fp_out, threshold, years, species = "cfin") {
  
  # -------- Create output directories --------
  dir.create(fp_out, showWarnings = FALSE) 
  dir.create(file.path(fp_out, species, version), showWarnings = FALSE)
  dir.create(file.path(fp_out, species, version, "Climatologies"), showWarnings = FALSE)
  # -------- Create directory for evaluations --------
  dir.create(file.path(fp_out, species, version, "Climatologies", "Evals"), showWarnings = FALSE)
  
  # Create color blind friendly pallete
  cbp <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
           "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  # Initialize months
  months <- c("January", "February", "March",
              "April", "May", "June",
              "July", "August", "September",
              "October", "November", "December")
  
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
            
            biomod_gam_abunds <- readr::read_csv(file.path(fp_out, species, version, "Biomod", "Projections", "gam_abund_vs_pred_2000_1.csv"))
            biomod_gam_full_data <- data.frame("year" = year, "month" = month, "abund" = biomod_gam_abunds$abund, "pred" = biomod_gam_abunds$gam_pred)

            biomod_brt_abunds <- readr::read_csv(file.path(fp_out, species, version, "Biomod", "Projections", "brt_abund_vs_pred_2000_1.csv"))
            biomod_brt_full_data <- data.frame("year" = year, "month" = month, "abund" = biomod_brt_abunds$abund, "pred" = biomod_brt_abunds$brt_pred)

            biomod_rf_abunds <- readr::read_csv(file.path(fp_out, species, version, "Biomod", "Projections", "rf_abund_vs_pred_2000_1.csv"))
            biomod_rf_full_data <- data.frame("year" = year, "month" = month, "abund" = biomod_rf_abunds$abund, "pred" = biomod_rf_abunds$rf_pred)

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
            
            biomod_gam_abunds <- readr::read_csv(file.path(fp_out, species, version, "Biomod", "Projections", paste0("gam_abund_vs_pred_", year, "_", month, ".csv")))

            biomod_gam_temp <- data.frame("year" = year, "month" = month, "abund" = biomod_gam_abunds$abund, "pred" = biomod_gam_abunds$gam_pred)
            biomod_gam_full_data <- rbind(biomod_gam_full_data, biomod_gam_temp)

            biomod_brt_abunds <- readr::read_csv(file.path(fp_out, species, version, "Biomod", "Projections", paste0("brt_abund_vs_pred_", year, "_", month, ".csv")))

            biomod_brt_temp <- data.frame("year" = year, "month" = month, "abund" = biomod_brt_abunds$abund, "pred" = biomod_brt_abunds$brt_pred)
            biomod_brt_full_data <- rbind(biomod_brt_full_data, biomod_brt_temp)

            biomod_rf_abunds <- readr::read_csv(file.path(fp_out, species, version, "Biomod", "Projections", paste0("rf_abund_vs_pred_", year, "_", month, ".csv")))

            biomod_rf_temp <- data.frame("year" = year, "month" = month, "abund" = biomod_rf_abunds$abund, "pred" = biomod_rf_abunds$rf_pred)
            biomod_rf_full_data <- rbind(biomod_rf_full_data, biomod_rf_temp)
          }
        }
      }
    }
  }
    
  for (i in 1:12) {
  
    # -------- Plot Taylor diagram --------
    # ---- GAMs ----
    abund_gam <- (gam_full_data %>% 
                    dplyr::filter(month == i) %>% 
                    dplyr::select(abund))$abund
    pred_gam <- (gam_full_data %>% 
                   dplyr::filter(month == i) %>% 
                   dplyr::select(pred))$pred
    
    taylor.diagram(abund_gam, pred_gam, pos.cor = TRUE, col = cbp[1], main = months[i])
    
    # ---- BRTs ----
    abund_brt <- (brt_full_data %>% 
                    dplyr::filter(month == i) %>% 
                    dplyr::select(abund))$abund
    pred_brt <- (brt_full_data %>% 
                   dplyr::filter(month == i) %>% 
                   dplyr::select(pred))$pred
      
    taylor.diagram(abund_brt, pred_brt, add = TRUE, col = cbp[2])
    
    # ---- Biomod ensemble ----
    abund_ensemble <- (biomod_full_data %>% 
                         dplyr::filter(month == i) %>% 
                         dplyr::select(abund))$abund
    pred_ensemble <- (biomod_full_data %>% 
                        dplyr::filter(month == i) %>% 
                        dplyr::select(pred))$pred
    
    taylor.diagram(abund_ensemble, pred_ensemble, add = TRUE, col = cbp[3])
    
    # ---- Biomod gam ----
    abund_gam <- (biomod_gam_full_data %>% 
                         dplyr::filter(month == i) %>% 
                         dplyr::select(abund))$abund
    pred_gam <- (biomod_gam_full_data %>% 
                        dplyr::filter(month == i) %>% 
                        dplyr::select(pred))$pred
    
    taylor.diagram(abund_gam, pred_gam, add = TRUE, col = cbp[4])
    
    # ---- Biomod brt ----
    abund_brt <- (biomod_brt_full_data %>% 
                         dplyr::filter(month == i) %>% 
                         dplyr::select(abund))$abund
    pred_brt <- (biomod_brt_full_data %>% 
                        dplyr::filter(month == i) %>% 
                        dplyr::select(pred))$pred
    
    taylor.diagram(abund_brt, pred_brt, add = TRUE, col = cbp[5])
    
    # ---- Biomod rf ----
    abund_rf <- (biomod_rf_full_data %>% 
                         dplyr::filter(month == i) %>% 
                         dplyr::select(abund))$abund
    pred_rf <- (biomod_rf_full_data %>% 
                        dplyr::filter(month == i) %>% 
                        dplyr::select(pred))$pred
    
    png(file.path(fp_out, species, version, "Climatologies", "Plots", paste0("taylor_diagram_", i, ".png")))
    taylor.diagram(abund_rf, pred_rf, add = TRUE, col = cbp[6])
    legend(1.5, 1.5, cex = 1.2, pt.cex = 1.2, legend = c("GAM", "BRT", "Ensemble",
                                                         "Biomod GAM", "Biomod BRT", "Biomod RF"), pch=19,
           col = c(cbp[1], cbp[2], cbp[3],
                   cbp[4], cbp[5], cbp[6]))
    dev.off()
    
  }
}
