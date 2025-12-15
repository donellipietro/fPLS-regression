## room for solutions
results_evaluation <- list()

## load results_evaluation if available
if (file.exists(paste(path_batch, "batch_", i, "_results_evaluation.RData", sep = ""))) {
  load(paste(path_batch, "batch_", i, "_results_evaluation.RData", sep = ""))
}

### Model MV-PLS ----
model_MV_PLS <- NULL
file_model <- paste(path_batch, "batch_", i, "_fitted_model_MV_PLS.RData", sep = "")
if (file.exists(file_model) && !FORCE_FIT) {
  if(FORCE_EVALUATE) {
    cat("- Loading fitted MV-PLS ... \n")
    load(file_model)
  }
} else if ("MV_PLS" %in% names_models) {
  cat("- Fitting MV-PLS ... ")
  
  ## model fit
  start.time <- Sys.time()
  model_MV_PLS <- MV_PLS_wrapped(data, n_comp = n_comp, center = FALSE, mode = mode_MV)
  end.time <- Sys.time()
  cat(paste("finished after", end.time - start.time, attr(end.time - start.time, "units"), "\n"))
  
  ## add flags
  model_MV_PLS$model_traits$is_functional <- FALSE
  model_MV_PLS$model_traits$has_interpolator <- FALSE
  
  ## add execution time to results
  model_MV_PLS$results$execution_time <- end.time - start.time
  
  ## save fitted model
  save(
    index_batch = i,
    model_MV_PLS,
    file = file_model
  )
}
if(!is.null(model_MV_PLS)) {
  
  ## model evaluation
  results_evaluation$MV_PLS <- evaluate_results(model_MV_PLS, generated_data)$indexes
  
  ## clean the workspace
  rm(model_MV_PLS)
}


### Model fPLS no calibration ----
model_fPLS_off <- NULL
file_model <- paste(path_batch, "batch_", i, "_fitted_model_fPLS_off.RData", sep = "")
if (file.exists(file_model) && !FORCE_FIT) {
  if(FORCE_EVALUATE) {
    cat("- Loading fitted fPLS sequential (no calibration) ... \n")
    load(file_model)
  }
} else if ("fPLS_off" %in% names_models) {
  cat("- Fitting fPLS sequential (no calibration) ... ")
  
  ## model fit
  start.time <- Sys.time()
  model_fPLS_off <- fdaPDE2::fPLS(
    data = data,
    center = FALSE,
    mode = mode_fun
  )
  model_fPLS_off$fit(lambda_fixed, n_comp = n_comp)
  end.time <- Sys.time()
  cat(paste("finished after", end.time - start.time, attr(end.time - start.time, "units"), "\n"))
  
  ## add flag
  model_fPLS_off$model_traits$is_functional <- TRUE
  model_fPLS_off$model_traits$has_interpolator <- TRUE
  
  ## add execution time to results
  model_fPLS_off$results$execution_time <- end.time - start.time
  
  ## evaluating fields at locations
  model_fPLS_off <- evaluate_fPLS_at_locations(model_fPLS_off, n_comp, mesh, nodes_HR)
  
  ## save fitted model
  save(
    index_batch = i,
    model_fPLS_off,
    file = file_model
  )
  
}
if(!is.null(model_fPLS_off)) {
  ## model evaluation
  results_evaluation$fPLS_off <- evaluate_results(model_fPLS_off, generated_data)$indexes
  
  ## clean the workspace
  rm(model_fPLS_off)
}


### Model fPLS (gcv calibration) ----
model_fPLS_gcv <- NULL
file_model <- paste(path_batch, "batch_", i, "_fitted_model_fPLS_gcv.RData", sep = "")
if (file.exists(file_model) && !FORCE_FIT) {
  if(FORCE_EVALUATE) {
    cat("- Loading fitted fPLS sequential (gcv calibration) ... \n")
    load(file_model)
  }
} else if ("fPLS_gcv" %in% names_models) {
  cat("- Fitting fPLS sequential (gcv calibration) ... ")
  
  ## model fit
  start.time <- Sys.time()
  model_fPLS_gcv <- fdaPDE2::fPLS(
    data = data,
    center = FALSE,
    solver = sequential(),
    mode = mode_fun
  )
  model_fPLS_gcv$fit(calibrator = gcv(lambda = lambda_grid, seed = seed), n_comp = n_comp)
  end.time <- Sys.time()
  cat(paste("finished after", end.time - start.time, attr(end.time - start.time, "units"), "\n"))
  
  ## add flag
  model_fPLS_gcv$model_traits$is_functional <- TRUE
  model_fPLS_gcv$model_traits$has_interpolator <- TRUE
  
  ## add execution time to results
  model_fPLS_gcv$results$execution_time <- end.time - start.time
  
  ## evaluating fields at locations
  model_fPLS_gcv <- evaluate_fPLS_at_locations(model_fPLS_gcv, n_comp, mesh, nodes_HR)
  
  ## save fitted model
  save(
    index_batch = i,
    model_fPLS_gcv,
    file = file_model
  )
}
if(!is.null(model_fPLS_gcv)) {
  ## model evaluation
  results_evaluation$fPLS_gcv <- evaluate_results(model_fPLS_gcv, generated_data)$indexes
  
  ## clean the workspace
  rm(model_fPLS_gcv)
}


### save results evaluation ----
if("generated_data" %in% ls())
results_evaluation$NSR_X <- generated_data$NSR_X

save(
  ## batch index
  index_batch = i,
  ## results
  results_evaluation,
  ## path
  file = paste(path_batch, "batch_", i, "_results_evaluation.RData", sep = "")
)

## clean the workspace
if("generated_data" %in% ls()) rm(generated_data)
rm(results_evaluation)

cat(paste("- Batch", i, "completed.\n"))
