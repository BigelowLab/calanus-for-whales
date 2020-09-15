# Camille Ross
# 20 August, 2020
# Purpose: write configuration files

# ---- Load libraries ----
library(yaml)

# ---- Initialize parameters ----

version <- "v0.2.6"
species <- c("cfin", "ctyp", "pseudo")[1]
env_covars <- c("sst", "chl", "bat", "sss", "wind",
                "lag_sst", "int_chl", "slope",
                "bott", "bots", "dist")
anomaly <- FALSE
threshold <- 40000
years <- 2000:2017

# ---- Initialize relative filepaths ----
fp_md <- "./calanus_data/Data/Databases/zooplankton_covar_data"
fp_covars <- "./Env_Covars"
# Must be full filepath to use rmarkdown::render(); see line 51
fp_out <- "~/Desktop/Calanus_Project/projects/calanus4whales/Models"

# ---- Initialize dataset ----
datasets <- c("CPR", "NOAA_CPR")
biomod_dataset <- "ECOMON"

# ---- Initialize data formatting options ----
format_data <- FALSE
fp_zpd <- NULL

# ---- Initialize yaml parameters ----
x <- list(version = version, species = species, env_covars = env_covars, anomaly = anomaly, threshold = threshold, years = years, 
         fp_md = fp_md, datasets = datasets, biomod_dataset = biomod_dataset, fp_covars = fp_covars,
         fp_out = fp_out, format_data = format_data, fp_zpd = fp_zpd)

# ---- Create directory ----
dir.create(file.path("./calanus-for-whales/Versions/", version))

# ---- Write yaml ----
write_yaml(x, file.path("./calanus-for-whales/Versions", version, paste0(version, ".yaml")))



