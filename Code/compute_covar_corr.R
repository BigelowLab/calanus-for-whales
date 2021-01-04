
source("./calanus-for-whales/Code/load_covars.R")

for (year in years) {
  for (month in 1:12) {
    if (year == 2000 & month == 1) {
      covars <- load_covars(fp_covars, year, month, env_covars = "all", as_raster = TRUE)
      
      wind <- covars$wind
      fetch <- covars$fetch
      chl <- covars$chl
      int_chl <- covars$int_chl
      sst <- covars$sst
      lag_sst <- covars$lag_sst
      uv <- covars$uv
      bat <- covars$bat
      dist <- covars$dist
      slope <- covars$slope
      wind <- covars$wind
      bots <- covars$bots
      bott <- covars$bott
      sss <- covars$sss
      
    } else {
      covars <- load_covars(fp_covars, year, month, env_covars = "all", as_raster = TRUE)
      
      wind <- stack(wind, covars$wind)
      fetch <- stack(fetch, covars$fetch)
      chl <- stack(chl, covars$chl)
      int_chl <- stack(int_chl, covars$int_chl)
      sst <- stack(sst, covars$sst)
      lag_sst <- stack(lag_sst, covars$lag_sst)
      uv <- stack(uv, covars$uv)
      bat <- stack(bat, covars$bat)
      dist <- stack(dist, covars$dist)
      slope <- stack(slope, covars$slope)
      wind <- stack(wind, covars$wind)
      bots <- stack(bots, covars$bots)
      bott <- stack(bott, covars$bott)
      sss <- stack(sss, covars$sss)

    }
  }
}


covars_clim <- calc(covars, fun = mean, na.rm = TRUE)


