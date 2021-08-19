
library(biomod2)

load_biomod_data <- function(md, j, threshold, env_covars, species, version) {
  
  
  # -------- Isolate month data --------
  trainingData <- md %>% dplyr::filter(month == j)
  # -------- Select presence or absence based on right whale feeding threshold --------
  trainingData$pa <- if_else(trainingData$abund < threshold, 0, 1)
  
  if (nrow(trainingData) < 100 | length(unique(trainingData$pa)) != 2) {
    print("skipped")
    skipped = TRUE
  } else {
    skipped = FALSE
  }
  
  if (skipped) {
    return(NULL)
  } else {
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
    
    return(biomodData)
  }
}