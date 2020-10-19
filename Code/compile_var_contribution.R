# Camille Ross
# 19 October, 2020
# Purpose: Compile variable contribution from multiple years of model runs

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
#'@param env_covars <list> list of environmental covariates
#'@param species <chr> species to model; choices are "cfin", "ctyp", or "pcal"
#'@param anomaly <logical> whether or not to convert data to an anomaly
compile_var_contribution <- function(version, fp_out, years, env_covars, species = "cfin", anomaly = FALSE) {
  
  # -------- Create output directories --------
  dir.create(fp_out, showWarnings = FALSE) 
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
          # Read in intial variable importance dataframe
          df_var <- read_csv(file.path(fp_out, species, version, "Biomod", "Evals", '1_var_importance.csv')) %>%
            dplyr::select(contains('RUN'))
          
          # Save backup
          var_cont <- df_var
          
          # dplyr::select BRT models (GBM)
          var_cont_brt <- df_var %>% dplyr::select(contains("GBM"))
          # Reformat
          var_cont_brt <- var_cont_brt %>% gather() %>% mutate(Var = rep(env_covars, 10),
                                                               Month = 1)
          # Rename columns
          names(var_cont_brt) <- c("Model", "Perc.Cont", "Var", "Month") 
          # Compute normalized percent contribution
          var_cont_brt <- var_cont_brt %>% mutate(Model = "BRT", 
                                                  Perc.Cont.Temp = rep(rowMeans(df_var %>% dplyr::select(contains("GBM")), na.rm = TRUE), 10)) %>%
            mutate(Perc.Cont.Norm = Perc.Cont/sum(unique(Perc.Cont.Temp)))                  
          
          # dplyr::select GAMs
          var_cont_gam <- df_var %>% dplyr::select(contains("GAM"))
          # Reformat
          var_cont_gam <- var_cont_gam %>% gather() %>% mutate(Var = rep(env_covars, 10),
                                                               Month = 1)
          # Rename columns
          names(var_cont_gam) <- c("Model", "Perc.Cont", "Var", "Month") 
          # Compute normalized percent contribution
          var_cont_gam <- var_cont_gam %>% mutate(Model = "GAM", 
                                                  Perc.Cont.Temp = rep(rowMeans(df_var %>% dplyr::select(contains("GAM")), na.rm = TRUE), 10)) %>%
            mutate(Perc.Cont.Norm = Perc.Cont/sum(unique(Perc.Cont.Temp))) 
          
          # dplyr::select RFs
          var_cont_rf <- df_var %>% dplyr::select(contains("RF"))
          # Reformat
          var_cont_rf <- var_cont_rf %>% gather() %>% mutate(Var = rep(c("sst", "btm_temp", "bat", "cal", "chlor"),10),
                                                               Month = 1)
          # Rename columns
          names(var_cont_rf) <- c("Model", "Perc.Cont", "Var", "Month") 
          # Compute normalized percent contribution
          var_cont_rf <- var_cont_rf %>% mutate(Model = "RF", 
                                                  Perc.Cont.Temp = rep(rowMeans(df_var %>% dplyr::select(contains("RF")), na.rm = TRUE), 10)) %>%
            mutate(Perc.Cont.Norm = Perc.Cont/sum(unique(Perc.Cont.Temp))) 
          
          # Bind individual model dataframes together
          var_cont <- rbind(var_cont_brt, var_cont_gam, var_cont_rf)
          
        } else {
          
          # Load variable importance for month
          df_var <- read_csv(file.path(fp_out, species, version, "Biomod", "Evals", paste0(i,'_var_importance.csv'))) %>%
            dplyr::select(contains('RUN'))
          
          # Select BRTs (GBMs)
          var_cont_brt <- df_var %>% dplyr::select(contains("GBM"))
          # Reformat
          var_cont_brt <- var_cont_brt %>% gather() %>% mutate(Var = rep(c("sst", "btm_temp", "bat", "cal", "chlor"),10),
                                                               Month = i)
          # Rename columns
          names(var_cont_brt) <- c("Model", "Perc.Cont", "Var", "Month") 
          # Compute normalized percent contribution
          var_cont_brt <- var_cont_brt %>% mutate(Model = "BRT", 
                                                  Perc.Cont.Temp = rep(rowMeans(df_var %>% dplyr::select(contains("GBM")), na.rm = TRUE), 10)) %>%
            mutate(Perc.Cont.Norm = Perc.Cont/sum(unique(Perc.Cont.Temp)))
          
          # Select GAMs
          var_cont_gam <- df_var %>% dplyr::select(contains("GAM"))
          # Reformat
          var_cont_gam <- var_cont_gam %>% gather() %>% mutate(Var = rep(c("sst", "btm_temp", "bat", "cal", "chlor"),10),
                                                               Month = i)
          # Rename columns
          names(var_cont_gam) <- c("Model", "Perc.Cont", "Var", "Month") 
          # Compute normalized percent contribution
          var_cont_gam <- var_cont_gam %>% mutate(Model = "GAM", 
                                                  Perc.Cont.Temp = rep(rowMeans(df_var %>% dplyr::select(contains("GAM")), na.rm = TRUE), 10)) %>%
            mutate(Perc.Cont.Norm = Perc.Cont/sum(unique(Perc.Cont.Temp)))
          
          # Select RFs
          var_cont_rf <- df_var %>% dplyr::select(contains("RF"))
          # Reformat
          var_cont_rf <- var_cont_rf %>% gather() %>% mutate(Var = rep(c("sst", "btm_temp", "bat", "cal", "chlor"),10),
                                                               Month = i)
          # Rename columns
          names(var_cont_rf) <- c("Model", "Perc.Cont", "Var", "Month") 
          # Compute normalized percent contribution
          var_cont_rf <- var_cont_rf %>% mutate(Model = "RF", 
                                                  Perc.Cont.Temp = rep(rowMeans(df_var %>% dplyr::select(contains("RF")), na.rm = TRUE), 10)) %>%
            mutate(Perc.Cont.Norm = Perc.Cont/sum(unique(Perc.Cont.Temp)))  
          
          # Bind dataframes togther
          var_cont_temp <- rbind(var_cont_brt, var_cont_gam, var_cont_rf)
          
          # Bind to initial variable importance dataframe
          var_cont <- rbind(var_cont, var_cont_temp)

        }
      }
    }
  }
  
  
  # Plot BRT variable contributions
  ggplot(data = (var_cont %>% filter(Model == "BRT")), aes(x = Month, fill = factor(Var, env_covars))) +
    geom_col(aes(y = Perc.Cont.Norm)) +
    scale_x_continuous(breaks = seq(1,12,2), labels = c("Jan", "Mar", "May", "Jul", "Sep", "Nov")) + 
    scale_fill_manual(labels = env_covars, values = cbp[c(3,5,6,8,4)]) +
    labs(y = "Contribution", fill = "Covariate") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    ggsave(filename = file.path(fp_out, species, version, "Biomod_Climatologies", "Evals", "gam_var_cont.png"))
  
  # Plot GAM variable contributions
  ggplot(data = (var_cont %>% filter(Model == "GAM")), aes(x = Month, fill = factor(Var, env_covars))) +
    geom_col(aes(y = Perc.Cont.Norm)) +
    scale_x_continuous(breaks = seq(1,12,2), labels = c("Jan", "Mar", "May", "Jul", "Sep", "Nov")) + 
    scale_fill_manual(labels = env_covars, values = cbp[c(3,5,6,8,4)]) +
    labs(y = "Contribution", fill = "Covariate") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    ggsave(filename = file.path(fp_out, species, version, "Biomod_Climatologies", "Evals", "brt_var_cont.png"))
  
  # Plot RF variable contributions
  ggplot(data = (var_cont %>% filter(Model == "RF")), aes(x = Month, fill = factor(Var, env_covars))) +
    geom_col(aes(y = Perc.Cont.Norm)) +
    scale_x_continuous(breaks = seq(1,12,2), labels = c("Jan", "Mar", "May", "Jul", "Sep", "Nov")) + 
    scale_fill_manual(labels = env_covars, values = cbp[c(3,5,6,8,4)]) +
    labs(x = "Month", y = "Contribution", fill = "Covariate") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    ggsave(filename = file.path(fp_out, species, version, "Biomod_Climatologies", "Evals", "rf_var_cont.png"))
  
}



