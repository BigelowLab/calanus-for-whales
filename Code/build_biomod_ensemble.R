library(biomod2)

build_biomod_ensemble <- function(modelOut, j, fp_out, species, version) {
  
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
  write.csv(ensembleEvals, file = file.path(fp_out, species, version, "Biomod", "Evals", paste0("ensemble_evals_", j, ".csv")), row.names = TRUE)
  
  return(biomodEM)
}