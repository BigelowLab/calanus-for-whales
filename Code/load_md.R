
# Source data formatting function
source("./calanus-for-whales/Code/format_model_data.R")
# Source data binding function
source("./calanus_data/Code/bind_years.R")

load_md <- function(version, fp_md, biomod_dataset, fp_covars, env_covars, 
                    years, fp_out, species,
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
    md$abund <- as.data.frame(md[paste0(species, "_CIV")])$cfin_CIV + 
      as.data.frame(md[paste0(species, "_CV_VI")])$cfin_CV_VI
  } else if (species == "ctyp") {
    md$abund <- as.data.frame(md[paste0(species, "_total")])$ctyp_total
  } else if (species == "pcal") {
    md$abund <- as.data.frame(md[paste0(species, "_total")])$pcal_total
  }
  
  # -------- Exclude NAs and select columns --------
  md <- md %>%
    as.data.frame() %>%
    na.exclude() %>%
    dplyr::mutate(season = if_else(month %in% c(1:3), 1,
                                   if_else(month %in% c(4:6), 2,
                                           if_else(month %in% c(7:9), 3, 4))))
  
  # -------- Take log of bathymetry and chlorophyll --------
  md$chl <- log(abs(md$chl))
  md$int_chl <- log(abs(md$int_chl))
  md$bat <- log(abs(md$bat))
  
  return(md)
  
}