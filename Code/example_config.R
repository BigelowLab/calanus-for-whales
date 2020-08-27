# Camille Ross
# 18 August, 2020
# Purpose: Example configuration file to run zooplankton models and summarize output.

# ---- Source necessary scripts ----
source("./calanus_for_whales/Code/build_gam.R")
source("./calanus_for_whales/Code/build_brt.R")
source("./calanus_for_whales/Code/build_biomod.R")
source("./calanus_for_whales/Code/build_climatology.R")
source("./calanus_for_whales/Code/build_biomod_climatology.R")
source("./calanus_for_whales/Code/compile_evals.R")
source("./calanus_for_whales/Code/plot_regions.R")

# ---- Initialize parameters ----
version <- "v0.2.3"
species <- c("cfin", "ctyp", "pseudo")[1]
env_covars <- c("wind", "lag_sst", "int_chl", "sss")
anomaly <- FALSE
years <- 2000:2017

# ---- Initialize relative filepaths ----
fp_md <- "./calanus_data/Data/Databases/zooplankton_covar_data"
fp_covars <- "./Env_Covars"
# Must be full filepath to use rmarkdown::render(); see line 51
fp_out <- "~/Desktop/Calanus_Project/projects/calanus4whales/calanus_for_whales/Versions"

# ---- Initialize data formatting options ----
format_data <- FALSE
fp_zpd <- NULL

# ---- Build GAM ----
build_gam(version = version, fp_md = fp_md, fp_covars = fp_covars, env_covars = env_covars, 
          years = years, fp_out = fp_out, species = species, anomaly = anomaly, 
          format_data = format_data, fp_zpd = fp_zpd)

# ---- Build BRT ----
build_brt(version = version, fp_md = fp_md, fp_covars = fp_covars, env_covars = env_covars, 
          years = years, fp_out = fp_out, species = species, anomaly = anomaly, 
          format_data = format_data, fp_zpd = fp_zpd)

# ---- Build Biomod2 models ----
build_biomod(version = version, fp_md = fp_md, fp_covars = fp_covars, env_covars = env_covars, 
             years = years, fp_out = fp_out, species = species, 
             format_data = format_data, fp_zpd = fp_zpd)

# ---- Build climatology ----
build_climatology(version = version, fp_out = fp_out, years = years, species = species, 
                  anomaly = anomaly)

# ---- Build biomod2 ensemble climatology ----
build_biomod_climatology(version = version, fp_out = fp_out, years = years, species = species)

# ---- Compile evaluations ----
compile_evals(version = version, fp_out = fp_out, years = years, species = species, 
              anomaly = anomaly)

# ---- Create region-specific actual vs. predicted plots ----
plot_regions(version = version, fp_out = fp_out, species = species)

# ---- Render model summary ----
rmarkdown::render("~/Desktop/Calanus_Project/projects/calanus4whales/calanus_for_whales/Code/build_summary.Rmd", 
                  output_file = file.path("~/Desktop/Calanus_Project/projects/calanus4whales/calanus_for_whales/Versions", species, version, paste0(version, ".html")),
                  params = list(set_title = version,
                                fp_out = fp_out,
                                species = species,
                                version = version,
                                env_covars = env_covars))

