
source("./calanus-for-whales/Code/load_covars.R")

for (year in years) {
  print(year)
  for (month in 1:12) {
    print(month)
    if (year == 2000 & month == 1) {
      covars <- load_covars(fp_covars, year, month, env_covars = "all", as_raster = TRUE)
      
      wind <- covars$wind
      fetch <- covars$fetch
      chl <- covars$chl
      int_chl_var <- covars$int_chl
      sst <- covars$sst
      lag_sst <- covars$lag_sst
      sst_grad <- covars$sst_grad
      uv <- covars$uv
      uv_grad <- covars$uv_grad
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
      sst_grad <- stack(sst_grad, covars$sst_grad)
      uv <- stack(uv, covars$uv)
      uv_grad <- stack(uv_grad, covars$uv_grad)
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
sst_grad_clim <- calc(sst_grad, fun = mean, na.rm = TRUE)
uv_clim <- calc(uv, fun = mean, na.rm = TRUE)
uv_grad_clim <- calc(uv_grad, fun = mean, na.rm = TRUE)
bat_clim <- calc(bat, fun = mean, na.rm = TRUE)
dist_clim <- calc(dist, fun = mean, na.rm = TRUE)
slope_clim <- calc(slope, fun = mean, na.rm = TRUE)
bots_clim <- calc(bots, fun = mean, na.rm = TRUE)
bott_clim <- calc(bott, fun = mean, na.rm = TRUE)
sss_clim <- calc(sss, fun = mean, na.rm = TRUE)

covars_clim <- stack(wind_clim, fetch_clim, chl_clim,
                     int_chl_clim, sst_clim, lag_sst_clim, sst_grad_clim,
                     uv_clim, uv_grad_clim, bat_clim, dist_clim,
                     slope_clim, bots_clim, bott_clim,
                     sss_clim)

names(covars_clim) <- c("wind", "fetch", "chl",
                        "int_chl", "sst", "lag_sst", "sst_grad",
                        "uv", "uv_grad", "bat", "dist",
                        "slope", "bots", "bott",
                        "sss")

stats <- layerStats(covars_clim, 'pearson', na.rm = TRUE)
corr_matrix <- stats$'pearson correlation coefficient'


# -------- Load monthly environmental covariates --------
covars <- stack(file.path(fp_covars, "Climatology", paste0("climatology_all_covars_", 1)))

stats <- layerStats(covars, 'pearson', na.rm = TRUE)
corr_matrix <- stats$'pearson correlation coefficient'


