# Camille Ross
# 18 August, 2020
# Purpose: Code to run zooplankton models and summarize output.

library(yaml)

DIR <- "/Users/camille/Calanus_Project"
setwd(DIR)

# ---- Source necessary scripts ----
source("./calanus-for-whales/Code/build_gam.R")
source("./calanus-for-whales/Code/build_brt.R")
source("./calanus-for-whales/Code/build_biomod.R")
source("./calanus-for-whales/Code/build_climatology.R")
source("./calanus-for-whales/Code/build_biomod_climatology.R")
source("./calanus-for-whales/Code/compile_evals.R")
source("./calanus-for-whales/Code/plot_regions.R")
source("./calanus-for-whales/Code/plot_regions_biomod.R")
source("./calanus-for-whales/Code/compile_abund_vs_pred.R")
source("./calanus-for-whales/Code/compile_var_contribution.R")

# Function to run all models and analysis and build report
#'@param version <chr> name of .yaml file with model configuration information
run_models <- function(version) {

  # ---- Read in config file ----
  config <- read_yaml(file.path("./calanus-for-whales/Versions", version, paste0(version, ".yaml")))
  
  # ---- Build GAM ----
  build_gam(version = config$version, fp_md = config$fp_md, datasets = config$datasets, fp_covars = config$fp_covars, env_covars = config$env_covars, 
            years = config$years, fp_out = config$fp_out, species = config$species, anomaly = config$anomaly, 
            format_data = config$format_data, fp_zpd = config$fp_zpd)
  
  # ---- Build BRT ----
  build_brt(version = config$version, fp_md = config$fp_md, datasets = config$datasets, fp_covars = config$fp_covars, env_covars = config$env_covars, 
            years = config$years, fp_out = config$fp_out, species = config$species, anomaly = config$anomaly, 
            format_data = config$format_data, fp_zpd = config$fp_zpd)
  
  # ---- Build Biomod2 models ----
  build_biomod(version = config$version, fp_md = config$fp_md, biomod_dataset = config$biomod_dataset, fp_covars = config$fp_covars, env_covars = config$env_covars, 
               years = config$years, fp_out = config$fp_out, species = config$species, threshold = config$threshold,
               format_data = config$format_data, fp_zpd = config$fp_zpd)
  
  # ---- Build climatology ----
  build_climatology(version = config$version, fp_out = config$fp_out, years = config$years, species = config$species, 
                    anomaly = config$anomaly)
  
  # ---- Build biomod2 ensemble climatology ----
  build_biomod_climatology(version = config$version, fp_out = config$fp_out, years = config$years, species = config$species)
  
  # ---- Compile evaluations ----
  compile_evals(version = config$version, fp_out = config$fp_out, years = config$years, species = config$species, 
                anomaly = config$anomaly)
  
  # ---- Compile biomod variable contributions ----
  compile_var_contribution(version = config$version, fp_out = config$fp_out, years = config$years, 
                           env_covars = config$env_covars, species = config$species, anomaly = config$anomaly)
  
  # ---- Create region-specific actual vs. predicted plots ----
  plot_regions(version = config$version, fp_out = config$fp_out, datasets = config$datasets, species = config$species)
  
  # ---- Create region-specific actual vs. predicted plots for biomod ----
  plot_regions_biomod(version = config$version, fp_out = config$fp_out, threshold = config$threshold, biomod_dataset = config$biomod_dataset, species = config$species)
    
  # ---- Compile actual vs. predicted values ----
  compile_abund_vs_pred(version = config$version, fp_out = config$fp_out, threshold = config$threshold, years = config$years, species = config$species)
  
  # ---- Render model summary ----
  rmarkdown::render(file.path(DIR, "calanus-for-whales/Code/build_summary.Rmd"), 
                    output_file = file.path(DIR, "calanus-for-whales/Versions", config$version, paste0(config$version, ".html")),
                    params = list(set_title = config$version,
                                  fp_out = config$fp_out,
                                  species = config$species,
                                  version = config$version,
                                  env_covars = config$env_covars,
                                  threshold = config$threshold,
                                  years = config$years,
                                  datasets = config$datasets,
                                  biomod_dataset = config$biomod_dataset))
  
  # ---- Render model summary for Jason ----
  rmarkdown::render(file.path(DIR, "calanus-for-whales/Code/report_for_jason.Rmd"), 
                    output_file = file.path(DIR, "calanus-for-whales/Versions", config$version, paste0(config$version, "_for_jason.html")),
                    params = list(set_title = config$version,
                                  fp_out = config$fp_out,
                                  species = config$species,
                                  version = config$version,
                                  env_covars = config$env_covars,
                                  threshold = config$threshold,
                                  years = config$years,
                                  datasets = config$datasets,
                                  biomod_dataset = config$biomod_dataset))
  
}
