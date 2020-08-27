# Camille Ross
# 18 August, 2020
# Purpose: Example configuration file to run zooplankton models and summarize output.

# ---- Load libraries ----
library(yaml)

# ---- Source necessary scripts ----
source("./calanus_for_whales/Code/build_gam.R")
source("./calanus_for_whales/Code/build_brt.R")
source("./calanus_for_whales/Code/build_biomod.R")
source("./calanus_for_whales/Code/build_climatology.R")
source("./calanus_for_whales/Code/build_biomod_climatology.R")
source("./calanus_for_whales/Code/compile_evals.R")
source("./calanus_for_whales/Code/plot_regions.R")

# ---- Read in config file ----
config <- read_yaml("./calanus_for_whales/Versions/v0.2.3.yaml")

# ---- Build GAM ----
build_gam(version = config$version, fp_md = config$fp_md, fp_covars = config$fp_covars, env_covars = config$env_covars, 
          years = config$years, fp_out = config$fp_out, species = config$species, anomaly = config$anomaly, 
          format_data = config$format_data, fp_zpd = config$fp_zpd)

# ---- Build BRT ----
build_brt(version = config$version, fp_md = config$fp_md, fp_covars = config$fp_covars, env_covars = config$env_covars, 
          years = config$years, fp_out = config$fp_out, species = config$species, anomaly = config$anomaly, 
          format_data = config$format_data, fp_zpd = config$fp_zpd)

# ---- Build Biomod2 models ----
build_biomod(version = config$version, fp_md = config$fp_md, fp_covars = config$fp_covars, env_covars = config$env_covars, 
             years = config$years, fp_out = config$fp_out, species = config$species, 
             format_data = config$format_data, fp_zpd = config$fp_zpd)

# ---- Build climatology ----
build_climatology(version = config$version, fp_out = config$fp_out, years = config$years, species = config$species, 
                  anomaly = config$anomaly)

# ---- Build biomod2 ensemble climatology ----
build_biomod_climatology(version = config$version, fp_out = config$fp_out, years = config$years, species = config$species)

# ---- Compile evaluations ----
compile_evals(version = config$version, fp_out = config$fp_out, years = config$years, species = config$species, 
              anomaly = config$anomaly)

# ---- Create region-specific actual vs. predicted plots ----
plot_regions(version = config$version, fp_out = config$fp_out, species = config$species)

# ---- Render model summary ----
rmarkdown::render("~/Desktop/Calanus_Project/projects/calanus4whales/calanus_for_whales/Code/build_summary.Rmd", 
                  output_file = file.path("~/Desktop/Calanus_Project/projects/calanus4whales/calanus_for_whales/Versions", species, version, paste0(version, ".html")),
                  params = list(set_title = config$version,
                                fp_out = config$fp_out,
                                species = config$species,
                                version = config$version,
                                env_covars = config$env_covars))

