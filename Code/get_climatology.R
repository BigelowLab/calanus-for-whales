
# Source covariate loading function
source("./calanus-for-whales/Code/load_covars.R")

# -------- Main function --------
#'@param version <chr> version of model
#'@param fp_md <chr> file path to formatted data used in model
#'@param fp_covars <chr> file path to environmental covariate data
#'@param env_covars <vector> vector of covariates to include in the model
#'@param years <vector> years for which to run the model
#'@param fp_out <chr> file path save the data to 
#'@param species <chr> species to model; choices are "cfin", "ctyp", or "pseudo"
#'@param anomaly <logical> if true, model is run using calanus anomaly
#'@param format_data <logical> if true, data is formatted within function; only used if model_data is NULL
#'@param fp_zpd <chr> file path to the zooplankton database if data is formatted within function
get_climatology <- function(fp_covars, env_covars, month) {
  
  clim <- load_covars(fp_covars = fp_covars, year = 2000, month = month,
                      env_covars = env_covars,
                      as_raster = TRUE)
  
  for (i in 2001:2017) {
    temp <- load_covars(fp_covars = fp_covars, year = i, month = month,
                        env_covars = env_covars,
                        as_raster = TRUE)
    
    clim <- raster::stack(clim, temp)
    
  }
  
  clim <- raster::calc(clim, fun = mean, na.rm = TRUE)
  
  return(clim)
  
}