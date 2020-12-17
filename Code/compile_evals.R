# Camille Ross
# 17 July, 2020
# Purpose: Compile evaluations from multiple years of model runs

# -------- Load libraries --------
require(dplyr)
require(readr)
require(viridis)
require(Metrics)
require(mapdata)
require(maps)
require(ggplot2)
require(stats)
require(pals)
require(tidyverse)

# -------- Main function --------
#'@param version <chr> version of model
#'@param fp_out <chr> file path save the data to 
#'@param years <vectors> years for which to run the model
#'@param species <chr> species to model; choices are "cfin", "ctyp", or "pseudo"
compile_evals <- function(version, fp_out, years, species = "cfin", anomaly = FALSE) {
  
  # -------- Create output directories --------
  dir.create(fp_out, showWarnings = FALSE) 
  dir.create(file.path(fp_out, species, version, "Climatologies"), showWarnings = FALSE)
  # ---- Create directory for evaluations ----
  dir.create(file.path(fp_out, species, version, "Climatologies", "Evals"), showWarnings = FALSE)
  # -------- Create output directories for biomod --------
  dir.create(file.path(fp_out, species, version, "Biomod_Climatologies"), showWarnings = FALSE)
  # ---- Create directory for evaluations ----
  dir.create(file.path(fp_out, species, version, "Biomod_Climatologies", "Evals"), showWarnings = FALSE)
  
  # -------- Compile evals --------
  for (year in years) {
    for (month in 1:12) {
      # ---- Test if projection exists ----
      if (paste0("proj_", year, "_", month, ".tif") %in% list.files(file.path(fp_out, species, version, "GAMs", "Projections")) &
          paste0("proj_", year, "_", month, ".tif") %in% list.files(file.path(fp_out, species, version, "BRTs", "Projections"))) {
        if (year == 2000 & month == 1) {
          aic <- as.numeric(substring(readr::read_csv(file.path(fp_out, species, version, "GAMs", "Evals", paste0("aic_", year, "_", month, ".csv"))), 4))
          gam_rmse <- as.numeric(substring(readr::read_csv(file.path(fp_out, species, version,"GAMs", "Evals", paste0("rmse_", year, "_", month, ".csv"))), 4))
          gam_rsq <- as.numeric(substring(readr::read_csv(file.path(fp_out, species, version, "GAMs", "Evals", paste0("rsq_", year, "_", month, ".csv"))), 4))[1]
          gam_evals <- data.frame("year" = year, "month" = month, "aic" = aic, "rmse" = gam_rmse, "rsq" = gam_rsq)
          
          brt_rmse <- as.numeric(substring(readr::read_csv(file.path(fp_out, species, version, "BRTs", "Evals", paste0("rmse_", year, "_", month, ".csv"))), 4))
          brt_rsq <- as.numeric(substring(readr::read_csv(file.path(fp_out, species, version, "BRTs", "Evals", paste0("rsq_", year, "_", month, ".csv"))), 4))[1]
          brt_evals <- data.frame("year" = year, "month" = month, "rmse" = brt_rmse, "rsq" = brt_rsq)

          biomod_gam_auc <- as.numeric(readr::read_csv(file.path(fp_out, species, version, "Biomod", "Evals", paste0("evals_", year, "_", month, ".csv"))) %>%
            dplyr::filter(stringr::str_detect(Model.name, "GAM"), Eval.metric == "ROC") %>%
            dplyr::summarize(Testing.data = mean(Testing.data, na.rm = TRUE)))
          biomod_gam_tss <- as.numeric(readr::read_csv(file.path(fp_out, species, version, "Biomod", "Evals", paste0("evals_", year, "_", month, ".csv"))) %>%
            dplyr::filter(stringr::str_detect(Model.name, "GAM"), Eval.metric == "TSS") %>%
            dplyr::summarize(Testing.data = mean(Testing.data, na.rm = TRUE)))
          biomod_gam_evals <- data.frame("year" = year, "month" = month, "auc" = biomod_gam_auc, "tss" = biomod_gam_tss)
          
          biomod_brt_auc <- as.numeric(readr::read_csv(file.path(fp_out, species, version, "Biomod", "Evals", paste0("evals_", year, "_", month, ".csv"))) %>%
            dplyr::filter(stringr::str_detect(Model.name, "GBM"), Eval.metric == "ROC") %>%
            dplyr::summarize(Testing.data = mean(Testing.data, na.rm = TRUE)))
          biomod_brt_tss <- as.numeric(readr::read_csv(file.path(fp_out, species, version, "Biomod", "Evals", paste0("evals_", year, "_", month, ".csv"))) %>%
            dplyr::filter(stringr::str_detect(Model.name, "GBM"), Eval.metric == "TSS") %>%
            dplyr::summarize(Testing.data = mean(Testing.data, na.rm = TRUE)))
          biomod_brt_evals <- data.frame("year" = year, "month" = month, "auc" = biomod_brt_auc, "tss" = biomod_brt_tss)
          
          biomod_rf_auc <- as.numeric(readr::read_csv(file.path(fp_out, species, version, "Biomod", "Evals", paste0("evals_", year, "_", month, ".csv"))) %>%
            dplyr::filter(stringr::str_detect(Model.name, "RF"), Eval.metric == "ROC") %>%
            dplyr::summarize(Testing.data = mean(Testing.data, na.rm = TRUE)))
          biomod_rf_tss <- as.numeric(readr::read_csv(file.path(fp_out, species, version, "Biomod", "Evals", paste0("evals_", year, "_", month, ".csv"))) %>%
            dplyr::filter(stringr::str_detect(Model.name, "RF"), Eval.metric == "TSS") %>%
            dplyr::summarize(Testing.data = mean(Testing.data, na.rm = TRUE)))
          biomod_rf_evals <- data.frame("year" = year, "month" = month, "auc" = biomod_rf_auc, "tss" = biomod_rf_tss)
          
          ensemble_auc <- as.numeric(readr::read_csv(file.path(fp_out, species, version, "Biomod", "Evals", paste0("ensemble_evals_", year, "_", month, ".csv"))) %>%
            dplyr::filter(stringr::str_detect(Model.name, "EMmeanByROC"), Eval.metric == "ROC") %>%
            dplyr::summarize(Testing.data = mean(Testing.data, na.rm = TRUE)))
          ensemble_tss <- as.numeric(readr::read_csv(file.path(fp_out, species, version, "Biomod", "Evals", paste0("ensemble_evals_", year, "_", month, ".csv"))) %>%
            dplyr::filter(stringr::str_detect(Model.name, "EMmeanByROC"), Eval.metric == "TSS") %>%
            dplyr::summarize(Testing.data = mean(Testing.data, na.rm = TRUE)))
          ensemble_evals <- data.frame("year" = year, "month" = month, "auc" = ensemble_auc, "tss" = ensemble_tss)
          
        } else {
          aic <- as.numeric(substring(readr::read_csv(file.path(fp_out, species, version, "GAMs", "Evals", paste0("aic_", year, "_", month, ".csv"))), 4))
          gam_rmse <- as.numeric(substring(readr::read_csv(file.path(fp_out, species, version, "GAMs", "Evals", paste0("rmse_", year, "_", month, ".csv"))), 4))
          gam_rsq <- as.numeric(substring(readr::read_csv(file.path(fp_out, species, version, "GAMs", "Evals", paste0("rsq_", year, "_", month, ".csv"))), 4))[1]
          
          gam_temp <- data.frame("year" = year, "month" = month, "aic" = aic, "rmse" = gam_rmse, "rsq" = gam_rsq)
          gam_evals <- rbind(gam_evals, gam_temp)
          
          brt_rmse <- as.numeric(substring(readr::read_csv(file.path(fp_out, species, version, "BRTs", "Evals", paste0("rmse_", year, "_", month, ".csv"))), 4))
          brt_rsq <- as.numeric(substring(readr::read_csv(file.path(fp_out, species, version, "BRTs", "Evals", paste0("rsq_", year, "_", month, ".csv"))), 4))[1]
          
          brt_temp <- data.frame("year" = year, "month" = month, "rmse" = brt_rmse, "rsq" = brt_rsq)
          brt_evals <- rbind(brt_evals, brt_temp)
          
          if (paste0("evals_", year, "_", month, ".csv") %in% list.files(file.path(fp_out, species, version, "Biomod", "Evals"))) {
            biomod_gam_auc <- as.numeric(readr::read_csv(file.path(fp_out, species, version, "Biomod", "Evals", paste0("evals_", year, "_", month, ".csv"))) %>%
              dplyr::filter(stringr::str_detect(Model.name, "GAM"), Eval.metric == "ROC") %>%
              dplyr::summarize(Testing.data = mean(Testing.data)))
            biomod_gam_tss <- as.numeric(readr::read_csv(file.path(fp_out, species, version, "Biomod", "Evals", paste0("evals_", year, "_", month, ".csv"))) %>%
              dplyr::filter(stringr::str_detect(Model.name, "GAM"), Eval.metric == "TSS") %>%
              dplyr::summarize(Testing.data = mean(Testing.data)))
              
            biomod_gam_temp <- data.frame("year" = year, "month" = month, "auc" = biomod_gam_auc, "tss" = biomod_gam_tss)
            biomod_gam_evals <- rbind(biomod_gam_evals, biomod_gam_temp)
            
            biomod_brt_auc <- as.numeric(readr::read_csv(file.path(fp_out, species, version, "Biomod", "Evals", paste0("evals_", year, "_", month, ".csv"))) %>%
                                           dplyr::filter(stringr::str_detect(Model.name, "GBM"), Eval.metric == "ROC") %>%
                                           dplyr::summarize(Testing.data = mean(Testing.data)))
            biomod_brt_tss <- as.numeric(readr::read_csv(file.path(fp_out, species, version, "Biomod", "Evals", paste0("evals_", year, "_", month, ".csv"))) %>%
                                           dplyr::filter(stringr::str_detect(Model.name, "GBM"), Eval.metric == "TSS") %>%
                                           dplyr::summarize(Testing.data = mean(Testing.data)))
            
            biomod_brt_temp <- data.frame("year" = year, "month" = month, "auc" = biomod_brt_auc, "tss" = biomod_brt_tss)
            biomod_brt_evals <- rbind(biomod_brt_evals, biomod_brt_temp)
            
            biomod_rf_auc <- as.numeric(readr::read_csv(file.path(fp_out, species, version, "Biomod", "Evals", paste0("evals_", year, "_", month, ".csv"))) %>%
                                           dplyr::filter(stringr::str_detect(Model.name, "RF"), Eval.metric == "ROC") %>%
                                           dplyr::summarize(Testing.data = mean(Testing.data)))
            biomod_rf_tss <- as.numeric(readr::read_csv(file.path(fp_out, species, version, "Biomod", "Evals", paste0("evals_", year, "_", month, ".csv"))) %>%
                                           dplyr::filter(stringr::str_detect(Model.name, "RF"), Eval.metric == "TSS") %>%
                                           dplyr::summarize(Testing.data = mean(Testing.data)))
            
            biomod_rf_temp <- data.frame("year" = year, "month" = month, "auc" = biomod_rf_auc, "tss" = biomod_rf_tss)
            biomod_rf_evals <- rbind(biomod_rf_evals, biomod_rf_temp)
          }
          
          if (paste0("ensemble_evals_", year, "_", month, ".csv") %in% list.files(file.path(fp_out, species, version, "Biomod", "Evals"))) {
            ensemble_auc <- as.numeric(readr::read_csv(file.path(fp_out, species, version, "Biomod", "Evals", paste0("ensemble_evals_", year, "_", month, ".csv"))) %>%
                                       dplyr::filter(stringr::str_detect(Model.name, "EMmeanByROC"), Eval.metric == "ROC") %>%
                                       dplyr::summarize(Testing.data = mean(Testing.data, na.rm = TRUE)))
            ensemble_tss <- as.numeric(readr::read_csv(file.path(fp_out, species, version, "Biomod", "Evals", paste0("ensemble_evals_", year, "_", month, ".csv"))) %>%
                                       dplyr::filter(stringr::str_detect(Model.name, "EMmeanByROC"), Eval.metric == "TSS") %>%
                                       dplyr::summarize(Testing.data = mean(Testing.data, na.rm = TRUE)))
            
            ensemble_temp <- data.frame("year" = year, "month" = month, "auc" = ensemble_auc, "tss" = ensemble_tss)
            ensemble_evals <- rbind(ensemble_evals, ensemble_temp)
          }
        }
      }
    }
  }
  
  # -------- Compute mean eval values --------
  gam_evals <- gam_evals %>% 
    dplyr::group_by(month) %>% 
    dplyr::mutate(mean_aic = mean(aic), mean_rmse = mean(rmse), mean_rsq = mean(rsq)) %>%
    # ---- Write to CSV ----
    readr::write_csv(file.path(fp_out, species, version, "Climatologies", "Evals", "gam_compiled_evals.csv"))
  
  brt_evals <- brt_evals %>% 
    dplyr::group_by(month) %>% 
    dplyr::mutate(mean_rmse = mean(rmse), mean_rsq = mean(rsq)) %>%
    # ---- Write to CSV ----
    readr::write_csv(file.path(fp_out, species, version, "Climatologies", "Evals", "brt_compiled_evals.csv"))
  
  biomod_gam_evals <- biomod_gam_evals %>% 
    dplyr::group_by(month) %>% 
    dplyr::mutate(mean_auc = mean(auc, na.rm = TRUE), mean_tss = mean(tss, na.rm = TRUE)) %>%
    # ---- Write to CSV ----
    readr::write_csv(file.path(fp_out, species, version, "Biomod_Climatologies", "Evals", "gam_compiled_evals.csv"))
  
  biomod_brt_evals <- biomod_brt_evals %>% 
    dplyr::group_by(month) %>% 
    dplyr::mutate(mean_auc = mean(auc, na.rm = TRUE), mean_tss = mean(tss, na.rm = TRUE)) %>%
    # ---- Write to CSV ----
  readr::write_csv(file.path(fp_out, species, version, "Biomod_Climatologies", "Evals", "brt_compiled_evals.csv"))
  
  biomod_rf_evals <- biomod_rf_evals %>% 
    dplyr::group_by(month) %>% 
    dplyr::mutate(mean_auc = mean(auc, na.rm = TRUE), mean_tss = mean(tss, na.rm = TRUE)) %>%
    # ---- Write to CSV ----
  readr::write_csv(file.path(fp_out, species, version, "Biomod_Climatologies", "Evals", "rf_compiled_evals.csv"))
  
  ensemble_evals <- ensemble_evals %>% 
    dplyr::group_by(month) %>% 
    dplyr::mutate(mean_auc = mean(auc, na.rm = TRUE), mean_tss = mean(tss, na.rm = TRUE)) %>%
    # ---- Write to CSV ----
    readr::write_csv(file.path(fp_out, species, version, "Biomod_Climatologies", "Evals", "ensemble_compiled_evals.csv"))
  
  
  # -------- Plot AIC --------
  ggplot(gam_evals, aes(x = month, y = aic, color = as.factor(year))) +
    geom_path() +
    scale_color_viridis(name = "Year", discrete=TRUE) +
    ylab("AIC") +
    xlab("Month") +
    guides(color = guide_legend(ncol = 2)) +
    scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position = "none") +
    ggsave(filename = file.path(fp_out, species, version, "Climatologies", "Evals", "gam_aic.png"))
  
  # -------- Plot RMSE --------
  # ---- GAMs ----
  ggplot(gam_evals, aes(x = month, y = rmse, color = as.factor(year))) +
    geom_path() +
    scale_color_viridis(name = "Year", discrete=TRUE) +
    ylab("RMSE") +
    xlab("Month") +
    guides(color = guide_legend(ncol = 2)) +
    scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position = "none") +
    ggsave(filename = file.path(fp_out, species, version, "Climatologies", "Evals", "gam_rmse.png"))
  
  # ---- BRTs ----
  ggplot(brt_evals, aes(x = month, y = rmse, color = as.factor(year))) +
    geom_path() +
    scale_color_viridis(name = "Year", discrete=TRUE) +
    ylab("RMSE") +
    xlab("Month") +
    guides(color = guide_legend(ncol = 2)) +
    scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position = "none") +
    ggsave(filename = file.path(fp_out, species, version, "Climatologies", "Evals", "brt_rmse.png"))
  
  # -------- Plot RSQ --------
  # ---- GAMs ----
  ggplot(gam_evals, aes(x = month, y = rsq, color = as.factor(year))) +
    geom_path() + 
    scale_color_viridis(name = "Year", discrete=TRUE) +
    ylab("RSQ") +
    xlab("Month") +
    guides(color = guide_legend(ncol = 2)) +
    scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position = "none") +
    ggsave(filename = file.path(fp_out, species, version, "Climatologies", "Evals", "gam_rsq.png"))
  
  # ---- BRTs ----
  ggplot(brt_evals, aes(x = month, y = rsq, color = as.factor(year))) +
    geom_path() + 
    scale_color_viridis(name = "Year", discrete=TRUE) +
    ylab("RSQ") +
    xlab("Month") +
    guides(color = guide_legend(ncol = 2)) +
    scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position = "none") +
    ggsave(filename = file.path(fp_out, species, version, "Climatologies", "Evals", "brt_rsq.png"))
  
  # -------- Plot AUC --------
  # ---- GAMs ----
  ggplot(biomod_gam_evals, aes(x = month, y = auc, color = as.factor(year))) +
    geom_path() +
    scale_color_viridis(name = "Year", discrete=TRUE) +
    ylim(c(0.6,1)) +
    ylab("AUC") +
    xlab("Month") +
    guides(color = guide_legend(ncol = 2)) +
    scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position = "none") +
    ggsave(filename = file.path(fp_out, species, version, "Biomod_Climatologies", "Evals", "gam_AUC.png"))
  
  # ---- BRTs ----
  ggplot(biomod_brt_evals, aes(x = month, y = auc, color = as.factor(year))) +
    geom_path() +
    scale_color_viridis(name = "Year", discrete=TRUE) +
    ylim(c(0.6,1)) +
    ylab("AUC") +
    xlab("Month") +
    guides(color = guide_legend(ncol = 2)) +
    scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position = "none") +
    ggsave(filename = file.path(fp_out, species, version, "Biomod_Climatologies", "Evals", "brt_AUC.png"))
  
  # ---- RF ----
  ggplot(biomod_rf_evals, aes(x = month, y = auc, color = as.factor(year))) +
    geom_path() +
    scale_color_viridis(name = "Year", discrete=TRUE) +
    ylim(c(0.6,1)) +
    ylab("AUC") +
    xlab("Month") +
    guides(color = guide_legend(ncol = 2)) +
    scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position = "none") +
    ggsave(filename = file.path(fp_out, species, version, "Biomod_Climatologies", "Evals", "rf_AUC.png"))
  
  
  # ---- Ensemble ----
  ggplot(ensemble_evals, aes(x = month, y = auc, color = as.factor(year))) +
    geom_path() +
    scale_color_viridis(name = "Year", discrete=TRUE) +
    ylim(c(0.6,1)) +
    ylab("AUC") +
    xlab("Month") +
    guides(color = guide_legend(ncol = 2)) +
    scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position = "none") +
    ggsave(filename = file.path(fp_out, species, version, "Biomod_Climatologies", "Evals", "ensemble_AUC.png"))
  
  # -------- Plot TSS --------
  # ---- GAMs ----
  ggplot(biomod_gam_evals, aes(x = month, y = tss, color = as.factor(year))) +
    geom_path() +
    scale_color_viridis(name = "Year", discrete=TRUE) +
    ylim(c(0,1)) +
    ylab("TSS") +
    xlab("Month") +
    guides(color = guide_legend(ncol = 2)) +
    scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position = "none") +
    ggsave(filename = file.path(fp_out, species, version, "Biomod_Climatologies", "Evals", "gam_TSS.png"))
  
  # ---- BRTs ----
  ggplot(biomod_brt_evals, aes(x = month, y = tss, color = as.factor(year))) +
    geom_path() +
    scale_color_viridis(name = "Year", discrete=TRUE) +
    ylim(c(0,1)) +
    ylab("TSS") +
    xlab("Month") +
    guides(color = guide_legend(ncol = 2)) +
    scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position = "none") +
    ggsave(filename = file.path(fp_out, species, version, "Biomod_Climatologies", "Evals", "brt_TSS.png"))
  
  # ---- GAMs ----
  ggplot(biomod_rf_evals, aes(x = month, y = tss, color = as.factor(year))) +
    geom_path() +
    scale_color_viridis(name = "Year", discrete=TRUE) +
    ylim(c(0,1)) +
    ylab("TSS") +
    xlab("Month") +
    guides(color = guide_legend(ncol = 2)) +
    scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position = "none") +
    ggsave(filename = file.path(fp_out, species, version, "Biomod_Climatologies", "Evals", "rf_TSS.png"))
  
  # ---- Ensemble ----
  ggplot(ensemble_evals, aes(x = month, y = tss, color = as.factor(year))) +
    geom_path() +
    scale_color_viridis(name = "Year", discrete=TRUE) +
    ylim(c(0,1)) +
    ylab("TSS") +
    xlab("Month") +
    guides(color = guide_legend(ncol = 2)) +
    scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position = "none") +
    ggsave(filename = file.path(fp_out, species, version, "Biomod_Climatologies", "Evals", "ensemble_TSS.png"))
  
}


