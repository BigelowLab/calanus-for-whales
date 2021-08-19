library(biomod2)

# Function to retrieve and save biomod evaluations and variable contributions
compile_biomod_evals <- function(modelOut, j, fp_out, species, version) {
  
  # ---------------------- SAVE EVALUATIONS & VARIABLE CONTRIBUTION ----------------------
  # Retrieves model evaluations
  #'@param obj <BIOMOD.models.out> model produced using BIOMOD_Modeling
  #'@param as.data.frame <boolean> if TRUE, function returns evaluations as a dataframe
  modelEvals <- get_evaluations(obj = modelOut, as.data.frame = TRUE)
  # Saves model evaluations to a csv file in the results directory
  # Allows for easy analysis
  write.csv(modelEvals, file = file.path(fp_out, species, version, "Biomod", "Evals", paste0("evals_", j, ".csv")), row.names = TRUE)
  
  # Retrieves variable contribution
  #'@param obj <BIOMOD.models.out> model produced using BIOMOD_Modeling
  #'@param as.data.frame <boolean> if TRUE, function returns variable contributions as a dataframe
  varImportance <- get_variables_importance(obj = modelOut, as.data.frame = TRUE)
  # Saves model evaluations to a csv file in the results directory
  # Allows for easy analysis
  write.csv(varImportance, file = file.path(fp_out, species, version, "Biomod", "Evals", paste0("var_importance_", j, ".csv")), row.names = TRUE)
  
}