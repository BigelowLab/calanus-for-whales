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
compile_abund_vs_pred <- function(version, fp_out, species = "cfin") {
  
  version = "v0.1.9"
  # -------- Create output directories --------
  dir.create(fp_out, showWarnings = FALSE) 
  dir.create(file.path(fp_out, species, "Climatologies"), showWarnings = FALSE)
  dir.create(file.path(fp_out, species, "Climatologies", version), showWarnings = FALSE)
  # -------- Create directory for evaluations --------
  dir.create(file.path(fp_out, species, "Climatologies", version, "Evals"), showWarnings = FALSE)
  
  
  # -------- Test if projection exists --------
  for (year in 2000:2017) {
    for (month in 1:12) {
      if (paste0("abund_vs_pred_", year, "_", month, ".csv") %in% list.files(file.path(fp_out, species, "GAMs", version, "Projections")) &
          paste0("abundvspred_", year, "_", month, ".csv") %in% list.files(file.path(fp_out, species, "BRTs", version, "Projections"))) {
        if (year == 2000 & month == 1) {
          gam_abunds <- readr::read_csv(file.path(fp_out, species, "GAMs", version, "Projections", "abund_vs_pred_2000_1.csv"))
          gam_full_data <- data.frame("year" = year, "month" = month, "abund" = gam_abunds$abund, "pred" = gam_abunds$pred)
          
          brt_abunds <- readr::read_csv(file.path(fp_out, species, "BRTs", version, "Projections", "abundvspred_2000_1.csv"))
          brt_full_data <- data.frame("year" = year, "month" = month, "abund" = brt_abunds$abund, "pred" = brt_abunds$pred)
          
        } else {
          gam_abunds <- readr::read_csv(file.path(fp_out, species, "GAMs", version, "Projections", paste0("abund_vs_pred_", year, "_", month, ".csv")))
          
          gam_temp <- data.frame("year" = year, "month" = month, "abund" = gam_abunds$abund, "pred" = gam_abunds$pred)
          gam_full_data <- rbind(gam_full_data, gam_temp)
          
          brt_abunds <- readr::read_csv(file.path(fp_out, species, "BRTs", version, "Projections", paste0("abundvspred_", year, "_", month, ".csv")))
          
          brt_temp <- data.frame("year" = year, "month" = month, "abund" = brt_abunds$abund, "pred" = brt_abunds$pred)
          brt_full_data <- rbind(brt_full_data, brt_temp)
        }
      }
    }
  }
  
  # -------- Write data to CSV --------
  gam_full_data <- gam_full_data %>% 
    readr::write_csv(file.path(fp_out, species, "Climatologies", version, "Projections", "gam_compiled_abund_vs_pred.csv")) %>%
    duplicated()
    dplyr::mutate(count = dplyr::count_(gam_full_data, vars = c('abund','pred'))$n)
  
  brt_full_data <- brt_full_data %>% 
    readr::write_csv(file.path(fp_out, species, "Climatologies", version, "Projections", "brt_compiled_abund_vs_pred.csv")) %>%
    dplyr::group_by(abund, pred)
  
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
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.position = "none") +
      ggsave(filename = file.path(fp_out, species, "Climatologies", version, "Plots", paste0("gam_abund_vs_pred_", i, ".png")))
  
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
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.position = "none") +
      ggsave(filename = file.path(fp_out, species, "Climatologies", version, "Plots", paste0("gam_abund_vs_pred_", i, ".png")))
  }
}


