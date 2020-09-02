# Camille Ross
# 20 August, 2020
# Purpose: write configuration files

# ---- Load libraries ----
library(yaml)

# ---- Initialize parameters ----

version <- "v0.2.4"
species <- c("cfin", "ctyp", "pseudo")[1]
env_covars <- c("wind", "lag_sst", "int_chl", "sss")
anomaly <- FALSE
years <- 2000:2017

# ---- Initialize relative filepaths ----
fp_md <- "./calanus_data/Data/Databases/zooplankton_covar_data"
fp_covars <- "./Env_Covars"
# Must be full filepath to use rmarkdown::render(); see line 51
fp_out <- "~/Desktop/Calanus_Project/projects/calanus4whales/Models"

# ---- Initialize data formatting options ----
format_data <- FALSE
fp_zpd <- NULL

# ---- Initialize yaml parameters ----
x <- list(version = version, species = species, env_covars = env_covars, anomaly = anomaly, years = years, 
         fp_md = fp_md, fp_covars = fp_covars,
         fp_out = fp_out, format_data = format_data, fp_zpd = fp_zpd)

# ---- Write yaml ----
write_yaml(x, file.path("./calanus-for-whales/Versions", paste0(version, ".yaml")))



