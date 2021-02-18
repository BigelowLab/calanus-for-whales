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
# Source data formatting function
source("./calanus-for-whales/Code/format_model_data.R")
# Source covariate loading function
source("./calanus-for-whales/Code/load_covars.R")
# Source covariate loading function
source("./calanus-for-whales/Code/get_climatology.R")
# Source data binding function
source("./calanus_data/Code/bind_years.R")

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
  
  # -------- Create output directories --------
  dir.create(fp_out, showWarnings = FALSE) 
  dir.create(file.path(fp_out, species), showWarnings = FALSE)
  dir.create(file.path(fp_out, species, version), showWarnings = FALSE)
  dir.create(file.path(fp_out, species, version, "Biomod"), showWarnings = FALSE)
  # -------- Create directory for projections --------
  dir.create(file.path(fp_out, species, version, "Biomod", "Projections"), showWarnings = FALSE)
  # -------- Create directory for plots --------
  dir.create(file.path(fp_out, species, version, "Biomod", "Plots"), showWarnings = FALSE)
  # -------- Create directory for interactions --------
  dir.create(file.path(fp_out, species, version, "Biomod", "Evals"), showWarnings = FALSE)
  
  # -------- Format model data --------
  if (format_data) {
    # Format zooplankton and environmental covariate data
    format_model_data(fp_data = fp_zpd, fp_covars = fp_covars, env_covars = "all", years = years, fp_out = fp_md)
  }
  
  # -------- Load model data --------
  if (length(years) == 1) {
    md <- readr::read_csv(file.path(fp_md, paste0(years[1], ".csv"))) %>% 
      dplyr::filter(dataset %in% biomod_dataset)
  } else {
    md <- bind_years(fp = file.path(fp_md), years = years) %>%
      dplyr::filter(dataset %in% biomod_dataset)
  }
  
  # -------- Compute anomaly --------
  if (species == "cfin") {
    md <- md %>% dplyr::group_by(dataset) %>%
      dplyr::mutate(mean = mean(log10(`cfin_CV_VI` + 1), na.rm = TRUE),
                    sd = sd(log10(`cfin_CV_VI` + 1), na.rm = TRUE),
                    anomaly = (log10(`cfin_CV_VI` + 1) - mean) / sd) %>%
      dplyr::ungroup()
  } else if (species == "ctyp") {
    md <- md %>% dplyr::group_by(dataset) %>%
      dplyr::mutate(mean = mean(log10(`ctyp_total` + 1), na.rm = TRUE),
                    sd = sd(log10(`ctyp_total` + 1), na.rm = TRUE),
                    anomaly = (log10(`ctyp_total` + 1) - mean) / sd) %>%
      dplyr::ungroup()
  } else if (species == "pcal") {
    md <- md %>% dplyr::group_by(dataset) %>%
      dplyr::mutate(mean = mean(log10(`pcal_total` + 1), na.rm = TRUE),
                    sd = sd(log10(`pcal_total` + 1), na.rm = TRUE),
                    anomaly = (log10(`pcal_total` + 1) - mean) / sd) %>%
      dplyr::ungroup()
  }
  
  # -------- Take the log of count data --------
  if (species == "cfin") {
    md$abund <- as.data.frame(md[paste0(species, "_CV_VI")])$cfin_CV_VI
  } else if (species == "ctyp") {
    md$abund <- as.data.frame(md[paste0(species, "_total")])$ctyp_total
  } else if (species == "pcal") {
    md$abund <- as.data.frame(md[paste0(species, "_total")])$pcal_total
  }
  
  # -------- Exclude NAs and select columns --------
  md <- md %>% dplyr::select(lat, lon, year, month, abund, wind, fetch, jday, chl, int_chl, bots, bott, sss, sst, sst_grad, lag_sst, uv, bat, dist, slope) %>%
    as.data.frame() %>%
    na.exclude() %>%
    dplyr::mutate(season = if_else(month %in% c(1:3), 1,
                                   if_else(month %in% c(4:6), 2,
                                           if_else(month %in% c(7:9), 3, 4))))
  
  # -------- Take log of bathymetry and chlorophyll --------
  md$chl <- log(abs(md$chl))
  md$int_chl <- log(abs(md$int_chl))
  md$bat <- log(abs(md$bat))
  
  # -------- Load world map data --------
  worldmap <- ggplot2::map_data("world")
  
  modelOptions <- BIOMOD_ModelingOptions(GAM = list(k = 4))
  
  for (j in 1:12) {
    
    # -------- Isolate month data --------
    trainingData <- md %>% dplyr::filter(month == j)
    # -------- Select presence or absence based on right whale feeding threshold --------
    trainingData$pa <- if_else(trainingData$abund < threshold, 0, 1)
    
    if (nrow(trainingData) < 100 | length(unique(trainingData$pa)) != 2) {
      print("skipped")
      next
    }
    
    # Isolate binary presence/absence data
    trainingPA <- as.data.frame(trainingData[,'pa'])
    # Isolate presence/absence coordiantes
    trainingXY <- trainingData[,c('lon', 'lat')]
    # Isolate environmental covariates
    # Select variables based on month
    trainingCovars <- as.data.frame(trainingData[, env_covars])
    
    # Format data for use in Biomod2 modelling function
    biomodData <- BIOMOD_FormatingData(resp.var = trainingPA,
                                       expl.var = trainingCovars,
                                       resp.xy = trainingXY,
                                       resp.name = paste0(species, version))
    
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
    
    # ---------------------- SAVE EVALUATIONS & VARIABLE CONTRIBUTION ----------------------
    # Retrieves model evaluations
    #'@param obj <BIOMOD.models.out> model produced using BIOMOD_Modeling
    #'@param as.data.frame <boolean> if TRUE, function returns evaluations as a dataframe
    modelEvals <- get_evaluations(obj = modelOut, as.data.frame = TRUE)
    # Saves model evaluations to a csv file in the results directory
    # Allows for easy analysis
    write.csv(modelEvals, file = file.path(fp_out, species, version, "Biomod", "Evals", paste0("evals.csv")), row.names = TRUE)
    
    # Retrieves variable contribution
    #'@param obj <BIOMOD.models.out> model produced using BIOMOD_Modeling
    #'@param as.data.frame <boolean> if TRUE, function returns variable contributions as a dataframe
    varImportance <- get_variables_importance(obj = modelOut, as.data.frame = TRUE)
    # Saves model evaluations to a csv file in the results directory
    # Allows for easy analysis
    write.csv(varImportance, file = file.path(fp_out, species, version, "Biomod", "Evals", paste0("var_importance_", i, "_", j, ".csv")), row.names = TRUE)
    
    
    # # GBM response curve; check biological plausibility
    # response.plot2(models  = BIOMOD_LoadModels(modelOut, models='GBM'),
    #                Data = get_formal_data(modelOut,'expl.var'), 
    #                show.variables= get_formal_data(modelOut,'expl.var.names'),
    #                do.bivariate = FALSE,
    #                fixed.var.metric = 'median',
    #                col = c("blue", "red", "green"),
    #                save.file = "png",
    #                name = file.path(fp_out, species, version, "Biomod", "Plots", paste0("BRT_response_curve_", i, "_", j, ".png")),
    #                legend = FALSE,
    #                main = "BRT Response Curves",
    #                data_species = get_formal_data(modelOut,'resp.var'))
    # 
    # # GAM response curve; check biological plausibility
    # response.plot2(models  = BIOMOD_LoadModels(modelOut, models='GAM'),
    #                Data = get_formal_data(modelOut,'expl.var'), 
    #                show.variables= get_formal_data(modelOut,'expl.var.names'),
    #                do.bivariate = FALSE,
    #                fixed.var.metric = 'median',
    #                col = c("blue", "red", "green"),
    #                save.file = "png",
    #                name = file.path(fp_out, species, version, "Biomod", "Plots", paste0("GAM_response_curve_", i, "_", j, ".png")),
    #                legend = FALSE,
    #                main = "GAM Response Curves",
    #                data_species = get_formal_data(modelOut,'resp.var'))
    # 
    #---------------------- BUILD ENSEMBLE MODEL ----------------------
    biomodEM <- BIOMOD_EnsembleModeling(
      modeling.output = modelOut,
      chosen.models = 'all',
      em.by = 'all',
      eval.metric = c('ROC'),
      eval.metric.quality.threshold = c(0.7),
      prob.mean = TRUE,
      prob.cv = TRUE,
      prob.ci = TRUE,
      prob.ci.alpha = 0.05,
      prob.median = TRUE,
      committee.averaging = TRUE,
      prob.mean.weight = TRUE,
      prob.mean.weight.decay = 'proportional')
    
    #Retrieves ensemble model evaluations
    ensembleEvals <- get_evaluations(obj = biomodEM, as.data.frame = TRUE)
    # Saves ensemble model evaluations to a csv file in the results directory
    # Allows for easy analysis
    write.csv(ensembleEvals, file = file.path(fp_out, species, version, "Biomod", "Evals", paste0("ensemble_evals_", i, "_", j, ".csv")), row.names = TRUE)
    
    # -------- Load environmental covariates for projection --------
    covars <- get_climatology(fp_covars = fp_covars, env_covars = env_covars, month = j)
    
    # -------- Project ensemble model --------
    biomodProjFuture <- BIOMOD_Projection(modeling.output = modelOut,
                                          new.env = covars,
                                          proj.name = 'all',
                                          selected.models = 'all',
                                          binary.meth = 'ROC',
                                          compress = 'xz',
                                          build.clamping.mask = FALSE,
                                          output.format = '.grd')
    
    # Build ensemble forecast
    myBiomodEF <- BIOMOD_EnsembleForecasting(EM.output = biomodEM,
                                             projection.output = biomodProjFuture,
                                             proj.name = paste0("Ensemble_", j))
    
    # Load ensemble forecast as raster
    # Divide by 1000 to convert probabilities to percentages
    ensemble_proj_raster <- raster(file.path(paste0(species, version), paste0("proj_Ensemble_", j), paste0("proj_Ensemble_", j, "_", species, version, "_ensemble.grd"))) %>%
      `/`(1000)
    
    crs(ensemble_proj_raster) <- '+init=epsg:4121 +proj=longlat +ellps=GRS80 +datum=GGRS87 +no_defs +towgs84=-199.87,74.79,246.62'
    
    # Save projection raster as data frame with xy coords and no NAs
    ensemble_proj_df <- as.data.frame(ensemble_proj_raster, xy = TRUE, na.rm = TRUE)
    
    # Assign column names
    names(ensemble_proj_df) <- c('pred', 'x', 'y')
    
    # -------- Plot projection --------
    ggplot() + 
      # Add projection data
      geom_tile(data = ensemble_proj_df, aes(x, y, fill = pred)) +
      # Add projection color gradient and label
      scale_fill_gradientn(colors = inferno(500), limits = c(0,1), na.value = "white") +
      labs(x = "", 
           y = "",
           fill = "Probability") +
      # Add world map data
      geom_polygon(data = worldmap, aes(long, lat, group = group), fill = NA, colour = "gray43") +
      coord_quickmap(xlim = c(round(min(ensemble_proj_df$x)), round(max(ensemble_proj_df$x))), 
                     ylim = c(round(min(ensemble_proj_df$y)), round(max(ensemble_proj_df$y))),
                     expand = TRUE) +
      # Remove grid lines
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      # Save plot to hard drive
      ggsave(filename = file.path(fp_out, species, version, "Biomod", "Plots", paste0("ensemble_proj_", j, ".png")), width = 7, height = 7)
    
    
    # -------- Extract predicted values --------
    month_md$ensemble_pred <- raster::extract(ensemble_proj_raster, month_md$abund)
    
    # -------- Convert NAs to zeros --------
    month_md$ensemble_pred[is.na(month_md$ensemble_pred)] <- 0
    month_md$abund[is.na(month_md$abund)] <- 0
    
    # -------- Unlist variables --------
    month_md$abund <- unlist(month_md$abund)
    month_md$ensemble_pred <- unlist(month_md$ensemble_pred)
    
    readr::write_csv(month_md %>% dplyr::select(abund, ensemble_pred), file.path(fp_out, species, version, "Biomod", "Projections", paste0("ensemble_abund_vs_pred_", i, "_", j, ".csv")))
    
    # -------- Plot actual vs. predicted values --------
    ggplot(data = month_md, aes(x = log10(abund + 1), y = ensemble_pred)) +
      geom_point() +
      ylim(c(0, 1)) +
      ggsave(filename = file.path(fp_out, species, version, "Biomod", "Plots", paste0("ensemble_actualvspred_", j, ".png")))
    
    # ------- Project GAMS -------
    # Create vector all runs of model algorithm for projection
    select_models <- c()
    for (k in 1:10) {
      select_models[k] <- paste0(species, version, "_AllData_RUN", k, "_GAM")
    }
    
    gamProj <- BIOMOD_Projection(modeling.output = modelOut,
                                 new.env = covars,
                                 proj.name = paste0("GAM_", j),
                                 selected.models = select_models,
                                 binary.meth = 'ROC',
                                 compress = 'xz',
                                 build.clamping.mask = TRUE,
                                 output.format = '.grd')
    
    
    # Load ensemble forecast as raster
    # Divide by 1000 to convert probabilities to percentages
    gam_proj_raster <- raster(file.path(paste0(species, version), paste0("proj_GAM_", j), paste0('proj_GAM_', j, "_", species, version, '.grd'))) %>%
      `/`(1000)
    
    crs(gam_proj_raster) <- '+init=epsg:4121 +proj=longlat +ellps=GRS80 +datum=GGRS87 +no_defs +towgs84=-199.87,74.79,246.62'
    
    # Save projection raster as data frame with xy coords and no NAs
    gam_proj_df <- as.data.frame(gam_proj_raster, xy = TRUE, na.rm = TRUE)
    
    # Assign column names
    names(gam_proj_df) <- c('pred', 'x', 'y')
    
    # -------- Plot projection --------
    ggplot() + 
      # Add projection data
      geom_tile(data = gam_proj_df, aes(x, y, fill = pred)) +
      # Add projection color gradient and label
      scale_fill_gradientn(colors = inferno(500), limits = c(0,1), na.value = "white") +
      labs(x = "", 
           y = "") +
      # Add world map data
      geom_polygon(data = worldmap, aes(long, lat, group = group), fill = NA, colour = "gray43") +
      coord_quickmap(xlim = c(round(min(gam_proj_df$x)), round(max(gam_proj_df$x))), 
                     ylim = c(round(min(gam_proj_df$y)), round(max(gam_proj_df$y))),
                     expand = TRUE) +
      # Remove grid lines
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      # Save plot to hard drive
      ggsave(filename = file.path(fp_out, species, version, "Biomod", "Plots", paste0("gam_proj_", j, ".png")), width = 7, height = 7)    
    
    # -------- Extract predicted values --------
    month_md$gam_pred <- raster::extract(gam_proj_raster, month_md$abund)
    
    # -------- Convert NAs to zeros --------
    month_md$gam_pred[is.na(month_md$gam_pred)] <- 0
    month_md$abund[is.na(month_md$abund)] <- 0
    
    # -------- Unlist variables --------
    month_md$abund <- unlist(month_md$abund)
    month_md$gam_pred <- unlist(month_md$gam_pred)
    
    readr::write_csv(month_md %>% dplyr::select(abund, gam_pred), file.path(fp_out, species, version, "Biomod", "Projections", paste0("gam_abund_vs_pred_", j, ".csv")))
    
    # -------- Plot actual vs. predicted values --------
    ggplot(data = month_md, aes(x = log10(abund + 1), y = gam_pred)) +
      geom_point() +
      ylim(c(0, 1)) +
      ggsave(filename = file.path(fp_out, species, version, "Biomod", "Plots", paste0("gam_actualvspred_", j, ".png")))
    
    # ------- Project BRT -------
    # Create vector all runs of model algorithm for projection
    select_models <- c()
    for (k in 1:10) {
      select_models[k] <- paste0(species, version, "_AllData_RUN", k, "_GBM")
    }
    
    brtProj <- BIOMOD_Projection(modeling.output = modelOut,
                                 new.env = covars,
                                 proj.name = paste0("GBM_", j),
                                 selected.models = select_models,
                                 binary.meth = 'ROC',
                                 compress = 'xz',
                                 build.clamping.mask = TRUE,
                                 output.format = '.grd')
    
    
    # Load ensemble forecast as raster
    # Divide by 1000 to convert probabilities to percentages
    brt_proj_raster <- raster(file.path(paste0(species, version), paste0('proj_GBM_', j), paste0('proj_GBM_', j, "_", species, version, '.grd'))) %>%
      `/`(1000)
    
    crs(brt_proj_raster) <- '+init=epsg:4121 +proj=longlat +ellps=GRS80 +datum=GGRS87 +no_defs +towgs84=-199.87,74.79,246.62'
    
    # Save projection raster as data frame with xy coords and no NAs
    brt_proj_df <- as.data.frame(brt_proj_raster, xy = TRUE, na.rm = TRUE)
    
    # Assign column names
    names(brt_proj_df) <- c('pred', 'x', 'y')
    
    # -------- Plot projection --------
    ggplot() + 
      # Add projection data
      geom_tile(data = brt_proj_df, aes(x, y, fill = pred)) +
      # Add projection color gradient and label
      scale_fill_gradientn(colors = inferno(500), limits = c(0,1), na.value = "white") +
      labs(x = "", 
           y = "") +
      # Add world map data
      geom_polygon(data = worldmap, aes(long, lat, group = group), fill = NA, colour = "gray43") +
      coord_quickmap(xlim = c(round(min(brt_proj_df$x)), round(max(brt_proj_df$x))), 
                     ylim = c(round(min(brt_proj_df$y)), round(max(brt_proj_df$y))),
                     expand = TRUE) +
      # Remove grid lines
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      # Save plot to hard drive
      ggsave(filename = file.path(fp_out, species, version, "Biomod", "Plots", paste0("brt_proj_", j, ".png")), width = 7, height = 7)    
    
    # -------- Extract predicted values --------
    month_md$brt_pred <- raster::extract(brt_proj_raster, month_md$abund)
    
    # -------- Convert NAs to zeros --------
    month_md$brt_pred[is.na(month_md$brt_pred)] <- 0
    month_md$abund[is.na(month_md$abund)] <- 0
    
    # -------- Unlist variables --------
    month_md$abund <- unlist(month_md$abund)
    month_md$brt_pred <- unlist(month_md$brt_pred)
    
    readr::write_csv(month_md %>% dplyr::select(abund, brt_pred), file.path(fp_out, species, version, "Biomod", "Projections", paste0("brt_abund_vs_pred_", j, ".csv")))
    
    # -------- Plot actual vs. predicted values --------
    ggplot(data = month_md, aes(x = log10(abund + 1), y = brt_pred)) +
      geom_point() +
      ylim(c(0, 1)) +
      ggsave(filename = file.path(fp_out, species, version, "Biomod", "Plots", paste0("brt_actualvspred_", j, ".png")))
    
    # ------- Project RF -------
    # Create vector all runs of model algorithm for projection
    select_models <- c()
    for (k in 1:10) {
      select_models[k] <- paste0(species, version, "_AllData_RUN", k, "_RF")
    }
    
    rfProj <- BIOMOD_Projection(modeling.output = modelOut,
                                new.env = covars,
                                proj.name = paste0("RF_", j),
                                selected.models = select_models,
                                binary.meth = 'ROC',
                                compress = 'xz',
                                build.clamping.mask = TRUE,
                                output.format = '.grd')
    
    
    # Load ensemble forecast as raster
    # Divide by 1000 to convert probabilities to percentages
    rf_proj_raster <- raster(file.path(paste0(species, version), paste0('proj_RF_', j), paste0('proj_RF_', j, "_", species, version, '.grd'))) %>%
      `/`(1000)
    
    crs(rf_proj_raster) <- '+init=epsg:4121 +proj=longlat +ellps=GRS80 +datum=GGRS87 +no_defs +towgs84=-199.87,74.79,246.62'
    
    # Save projection raster as data frame with xy coords and no NAs
    rf_proj_df <- as.data.frame(rf_proj_raster, xy = TRUE, na.rm = TRUE)
    
    # Assign column names
    names(rf_proj_df) <- c('pred', 'x', 'y')
    
    # -------- Plot projection --------
    ggplot() + 
      # Add projection data
      geom_tile(data = rf_proj_df, aes(x, y, fill = pred)) +
      # Add projection color gradient and label
      scale_fill_gradientn(colors = inferno(500), limits = c(0,1), na.value = "white") +
      labs(x = "", 
           y = "") +
      # Add world map data
      geom_polygon(data = worldmap, aes(long, lat, group = group), fill = NA, colour = "gray43") +
      coord_quickmap(xlim = c(round(min(rf_proj_df$x)), round(max(rf_proj_df$x))), 
                     ylim = c(round(min(rf_proj_df$y)), round(max(rf_proj_df$y))),
                     expand = TRUE) +
      # Remove grid lines
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      # Save plot to hard drive
      ggsave(filename = file.path(fp_out, species, version, "Biomod", "Plots", paste0("rf_proj_", j, ".png")), width = 7, height = 7)    
    
    
    # -------- Extract predicted values --------
    month_md$rf_pred <- raster::extract(rf_proj_raster, month_md$abund)
    
    # -------- Convert NAs to zeros --------
    month_md$rf_pred[is.na(month_md$rf_pred)] <- 0
    month_md$abund[is.na(month_md$abund)] <- 0
    
    # -------- Unlist variables --------
    month_md$abund <- unlist(month_md$abund)
    month_md$rf_pred <- unlist(month_md$rf_pred)
    
    readr::write_csv(month_md %>% dplyr::select(abund, rf_pred), file.path(fp_out, species, version, "Biomod", "Projections", paste0("rf_abund_vs_pred_", j, ".csv")))
    
    # -------- Plot actual vs. predicted values --------
    ggplot(data = month_md, aes(x = log10(abund + 1), y = rf_pred)) +
      geom_point() +
      ylim(c(0, 1)) +
      ggsave(filename = file.path(fp_out, species, version, "Biomod", "Plots", paste0("rf_actualvspred_", j, ".png")))
    
    
  }
}





