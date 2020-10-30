# Camille Ross
# 20 August, 2020
# Purpose: write configuration files

# ---- Load libraries ----
library(yaml)

# ---- Function to write configuration file ----
#'@param version <chr> model version in "vx.x.x" format
#'@param species <chr> species for which the model is run; one of "cfin", "ctyp", or "pcal"
#'@param env_covars <list> list of environmental covariates; options are wind, fetch, sst, sss,
# bott, bots, bat, slope, dist, chl, lag_sst, int_chl
#'@param anomaly <logical> whether or not model converts data to anomalies
#'@param threshold <num> right whale feeding thershold for biomod models
#'@param years <vector> years to run model for between 2000 and 2017
#'@param fp_md <filepath> filepath to model data
#'@param fp_covars <filepath> filepath to environmental covariates
#'@param fp_out <filepath> filepath to store model output
#'@param datasets <chr> datasets to use for mgcv gams and gbm brts
#'@param biomod_dataset <chr> datasets to use for biomod models
#'@param format_data <logical> whether or not to format sightings and covar data
#'@param fp_zpd <filepath> filepath to zooplankton data if format_data is TRUE
write_config <- function(version, species, env_covars, anomaly,
                         threshold, years, fp_md, fp_covars,
                         fp_out, datasets, biomod_dataset, format_data, 
                         fp_zpd) {

  # ---- Initialize yaml parameters ----
  x <- list(version = version, species = species, env_covars = env_covars, anomaly = anomaly, threshold = threshold, years = years, 
           fp_md = fp_md, datasets = datasets, biomod_dataset = biomod_dataset, fp_covars = fp_covars,
           fp_out = fp_out, format_data = format_data, fp_zpd = fp_zpd)
  
  # ---- Create directory ----
  dir.create(file.path("./calanus-for-whales/Versions/", version))
  
  # ---- Write yaml ----
  write_yaml(x, file.path("./calanus-for-whales/Versions", version, paste0(version, ".yaml")))

}

