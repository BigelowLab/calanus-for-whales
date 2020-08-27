# Camille Ross
# 16 June, 2020
# Purpose: Run models for zooplankton species with variable environmental covariates

# -------- Source outside files --------
# Source function to build model
source("./calanus_for_whales/build_model.R")

# Run monthly BRT models using dismo
build_brt_dismo(version = "v0.0.0", fp_md = "./projectdata/calanus4whales/Databases/zooplankton_covar_data", species = "cfin", 
            fp_covars = "./projectdata/calanus4whales/Env_Covars", env_covars = c("chl", "bott"), 
            years = 2000:2009, proj_year = 2009, fp_out = "./projects/calanus4whales/Models")

# Run monthly BRT models
build_brt(version = "v0.0.0", fp_md = "./projectdata/calanus4whales/Databases/zooplankton_covar_data", species = "cfin", 
          fp_covars = "./projectdata/calanus4whales/Env_Covars", env_covars = c("chl", "bott"), 
          years = 2000:2009, proj_year = 2009, fp_out = "./projects/calanus4whales/Models")

# Run monthly GAMs
build_gam(version = "v0.0.0", fp_md = "./projectdata/calanus4whales/Databases/zooplankton_covar_data", species = "cfin", 
          fp_covars = "./projectdata/calanus4whales/Env_Covars", env_covars = c("chl", "bott"), 
          years = 2000:2009, proj_year = 2009, fp_out = "./projects/calanus4whales/Models")

# Run monthly RF models
build_rf(version = "v0.0.0", fp_md = "./projectdata/calanus4whales/Databases/zooplankton_covar_data", species = "cfin", 
        fp_covars = "./projectdata/calanus4whales/Env_Covars", env_covars = c("chl", "bott"), 
        years = 2000:2009, proj_year = 2009, fp_out = "./projects/calanus4whales/Models")






