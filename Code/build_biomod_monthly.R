# Camille Ross
# 10 June, 2020
# Purpose: Run a boosted regression tree model for zooplankton species with variable environmental covariates

# -------- Load libraries --------
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
require(biomod2)
require(ff)
library(rgdal)

# -------- Source outside files --------
# Source data binding function
source("./calanus-for-whales/Code/load_md.R")
# Source data binding function
source("./calanus-for-whales/Code/load_biomod_data.R")
source("./calanus-for-whales/Code/compile_biomod_evals.R")
source("./calanus-for-whales/Code/build_biomod_ensemble.R")
source("./calanus-for-whales/Code/project_ensemble.R")
source("./calanus-for-whales/Code/project_biomod.R")

# -------- Main function --------
#'@param version <chr> version of model
#'@param fp_md <chr> file path to formatted data used in model
#'@param species <chr> species to model; choices are "cfin", "ctyp", or "pcal"
#'@param fp_covars <chr> file path to environmental covariate data
#'@param env_covars <vector> vector of covariates to include in the model
#'@param years <vectors> years for which to run the model
#'@param fp_out <chr> file path save the data to 
#'@param format_data <logical> if true, data is formatted within function; only used if model_data is NULL
#'@param fp_zpd <chr> filepath to the zooplankton database if data is formatted within function
build_biomod <- function(version, fp_md, biomod_dataset, fp_covars, env_covars, 
                         years, fp_out, species, threshold, 
                         format_data, fp_zpd) {
  
  md <- load_md(version, fp_md, biomod_dataset, fp_covars, env_covars, 
                years, fp_out, species,
                format_data, fp_zpd)
  
  # -------- Load world map data --------
  worldmap <- ggplot2::map_data("world")
  
  modelOptions <- BIOMOD_ModelingOptions(GAM = list(k = 4))
  
  for (j in 1:12) {
    
    biomodData <- load_biomod_data(md, j, threshold, env_covars, species, version) 
    
    if (!is.null(biomodData)) {
    
      plot(biomodData)
      
      # ---------------------- BUILD MODELS ----------------------
      modelOut <- BIOMOD_Modeling(data = biomodData,
                                  models = c("GAM", "GBM", "RF"),
                                  models.options = modelOptions,
                                  NbRunEval = 10,
                                  DataSplit = 70,
                                  Prevalence = 0.5,
                                  VarImport = 5,
                                  models.eval.meth = c('ROC', 'TSS', 'KAPPA'),
                                  SaveObj = TRUE,
                                  rescal.all.models = FALSE,
                                  do.full.models = FALSE,
                                  modeling.id = paste0(species, version))
      
      
      compile_biomod_evals(modelOut, j, fp_out, species, version)
      
      #---------------------- BUILD ENSEMBLE MODEL ----------------------
      biomodEM <- build_biomod_ensemble(modelOut, j, fp_out, species, version) 
      
      
      # -------- Load environmental covariates for projection --------
      
      covars <- stack(file.path(fp_covars, "Climatology", paste0("climatology_all_covars_", j))) %>%
        subset(env_covars)
      
      # Project the ensemble and save associated plot
      project_ensemble(modelOut, biomodEM, j, fp_out, species, version)
      
      # Project individual models
      models <- c("GAM", "GBM", "RF")
      
      for (model in models) {
        project_biomod(model, modelOut, j, species, version, fp_out)
      }
    }
  }
}





