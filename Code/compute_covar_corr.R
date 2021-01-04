
source("./calanus-for-whales/Code/load_covars.R")

for (year in years) {
  for (month in 1:12) {
    if (year == 2000 & month == 1) {
      covars <- load_covars(fp_covars, year, month, env_covars = "all", as_raster = TRUE)
      
      wind <- covars$wind
      fetch <- covars$fetch
      chl <- covars$chl
      int_chl_var <- covars$int_chl
      sst <- covars$sst
      lag_sst <- covars$lag_sst
      uv <- covars$uv
      bat <- covars$bat
      dist <- covars$dist
      slope <- covars$slope
      bots <- covars$bots
      bott <- covars$bott
      sss <- covars$sss
      
    } else {
      covars <- load_covars(fp_covars, year, month, env_covars = "all", as_raster = TRUE)
      
      wind <- stack(wind, covars$wind)
      fetch <- stack(fetch, covars$fetch)
      chl <- stack(chl, covars$chl)
      int_chl_var <- stack(int_chl_var, covars$int_chl)
      sst <- stack(sst, covars$sst)
      lag_sst <- stack(lag_sst, covars$lag_sst)
      uv <- stack(uv, covars$uv)
      bat <- stack(bat, covars$bat)
      dist <- stack(dist, covars$dist)
      slope <- stack(slope, covars$slope)
      bots <- stack(bots, covars$bots)
      bott <- stack(bott, covars$bott)
      sss <- stack(sss, covars$sss)

    }
  }
}


wind_clim <- calc(wind, fun = mean, na.rm = TRUE)
fetch_clim <- calc(fetch, fun = mean, na.rm = TRUE)
chl_clim <- calc(chl, fun = mean, na.rm = TRUE)
int_chl_clim <- calc(int_chl_var, fun = mean, na.rm = TRUE)
sst_clim <- calc(sst, fun = mean, na.rm = TRUE)
lag_sst_clim <- calc(lag_sst, fun = mean, na.rm = TRUE)
uv_clim <- calc(uv, fun = mean, na.rm = TRUE)
bat_clim <- calc(bat, fun = mean, na.rm = TRUE)
dist_clim <- calc(dist, fun = mean, na.rm = TRUE)
slope_clim <- calc(slope, fun = mean, na.rm = TRUE)
bots_clim <- calc(bots, fun = mean, na.rm = TRUE)
bott_clim <- calc(bott, fun = mean, na.rm = TRUE)
sss_clim <- calc(sss, fun = mean, na.rm = TRUE)

covars_clim



