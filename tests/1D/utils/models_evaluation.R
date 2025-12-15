# performances evaluation ----

evaluate_fPLS <- function(model, n_comp) {
  # reconstruction
  model$results$X_hat_HR <- list()
  model$results$X_hat_locs <- list()
  for(h in 1:n_comp) {
    model$results$X_hat_HR[[h]] <- model$results$X_hat[[h]] %*% t(Psi_HR)
    model$results$X_hat_locs[[h]] <- model$results$X_hat[[h]] %*% t(Psi_locs)
  }
  # Beta
  model$results$Beta_hat_HR <- list()
  model$results$Beta_hat_locs <- list()
  if(model$MODE == "fPLS-R") {
    for(h in 1:n_comp) {
      model$results$Beta_hat_HR[[h]] <- Psi_HR %*% model$results$Beta_hat[[h]]
      model$results$Beta_hat_locs[[h]] <- Psi_locs %*% model$results$Beta_hat[[h]]
    }
  }
  # directions & loadings
  model$results$X_space_directions_HR <- Psi_HR %*% model$results$X_space_directions
  model$results$X_space_directions_locs <- Psi_locs %*% model$results$X_space_directions
  model$results$X_loadings_HR <- Psi_HR %*% model$results$X_loadings
  model$results$X_loadings_locs <- Psi_locs %*% model$results$X_loadings
  
  return(model)
}

evaluate_B_fPLS <- function(model, n_comp) {
  
  basisobj <- model$basisobj
  Psi_HR <- eval.basis(grid_HR, basisobj)
  Psi_locs <- eval.basis(locations, basisobj)
  
  # reconstruction
  model$results$X_hat_HR <- list()
  model$results$X_hat_locs <- list()
  for(h in 1:n_comp) {
    model$results$X_hat_HR[[h]] <- model$results$X_hat[[h]] %*% t(Psi_HR)
    model$results$X_hat_locs[[h]] <- model$results$X_hat[[h]] %*% t(Psi_locs)
  }
  # Beta
  model$results$Beta_hat_HR <- list()
  model$results$Beta_hat_locs <- list()
  if(model$MODE == "fPLS-R") {
    for(h in 1:n_comp) {
      model$results$Beta_hat_HR[[h]] <- Psi_HR %*% model$results$Beta_hat[[h]]
      model$results$Beta_hat_locs[[h]] <- Psi_locs %*% model$results$Beta_hat[[h]]
    }
  }
  # directions & loadings
  model$results$X_space_directions_HR <- Psi_HR %*% model$results$X_space_directions
  model$results$X_space_directions_locs <- Psi_locs %*% model$results$X_space_directions
  model$results$X_loadings_HR <- Psi_HR %*% model$results$X_loadings
  model$results$X_loadings_locs <- Psi_locs %*% model$results$X_loadings
  
  return(model)
}

RMSE_nf <- function(x, y) {
  if(any(is.na(x))){
    return(NA)
  }
  x <- x/norm_l2(x)
  y <- y/norm_l2(y)
  if(RMSE(x-y) > RMSE(x+y)) {
    return(RMSE(x+y))
  } else {
    return(RMSE(x-y))
  }
}

evaluate_results <- function(model, generated_data) {
  
  ## number of computed components
  n_comp <- generated_data$dimensions$n_comp
  
  ## room for results
  rmse <- list()
  irmse <- list()
  
  ## execution time
  execution_time <- model$results$execution_time
  
  ## RMSE Y and X (at locations) ----

  ## data reconstruction
  for(h in 1:n_comp) {
    norm <- 1 # RMSE(generated_data$expected_results$Y_true)
    rmse$Y_reconstruction[h] <- RMSE(model$results$Y_hat[[h]] - generated_data$expected_results$Y_true) / norm
    norm <- 1 # RMSE(generated_data$expected_results$X_true_locs)
    rmse$X_reconstruction_locs[h] <- RMSE(model$results$X_hat_locs[[h]] - generated_data$expected_results$X_true_locs) / norm
  }
  
  ## Beta
  if(model$MODE == "fPLS-R" || model$MODE == "PLS-R") {
    for(h in 1:n_comp) {
      norm <- 1 # RMSE(generated_data$expected_results$Beta_locs)
      rmse$Beta_locs[h] <- RMSE_nf(model$results$Beta_hat_locs[[h]], generated_data$expected_results$Beta_true_locs) / norm
    }
  }
  
  ## directions, loadings & scores
  for (h in 1:n_comp) {
    ## directions
    rmse$Y_space_directions[h] <- RMSE_nf(model$results$Y_space_directions[, h], generated_data$expected_results$Y_space_directions_true[, h])
    rmse$X_space_directions_locs[h] <- RMSE_nf(model$results$X_space_directions_locs[, h], generated_data$expected_results$X_space_directions_true_locs[, h])
    ## latent scores
    rmse$X_latent_scores[h] <- RMSE_nf(model$results$X_latent_scores[, h], generated_data$expected_results$X_latent_scores_true[, h])
    ## loadings
    rmse$Y_loadings[h] <- RMSE_nf(model$results$Y_loadings[, h], generated_data$expected_results$Y_loadings_true[, h])
    rmse$X_loadings_locs[h] <- RMSE_nf(model$results$X_loadings_locs[, h], generated_data$expected_results$X_loadings_true_locs[, h])
  }
  
  ## RMSE at HR (if possible) ----

  if (model$model_traits$has_interpolator) {

    ## data reconstruction
    for(h in 1:n_comp) {
      norm <- 1 # RMSE(generated_data$expected_results$X_true)
      rmse$X_reconstruction[h] <- RMSE(model$results$X_hat_HR[[h]] - generated_data$expected_results$X_true) / norm
    }

    ## Beta
    if(model$MODE == "fPLS-R" || model$MODE == "PLS-R") {
      for(h in 1:n_comp) {
        norm <- 1 # RMSE(generated_data$expected_results$Beta)
        rmse$Beta[h] <- RMSE_nf(model$results$Beta_hat_HR[[h]], generated_data$expected_results$Beta_true) / norm
      }
    }

    ## directions & loadings
    for (h in 1:n_comp) {
      ## directions
      rmse$X_space_directions[h] <- RMSE_nf(model$results$X_space_directions_HR[, h], generated_data$expected_results$X_space_directions_true[, h])
      ## loadings
      rmse$X_loadings[h] <- RMSE_nf(model$results$X_loadings_HR[, h], generated_data$expected_results$X_loadings_true[, h])
    }

  }
  
  return(list(
    indexes = list(
      execution_time = execution_time,
      rmse = rmse,
      irmse = irmse
    )
  ))
}
