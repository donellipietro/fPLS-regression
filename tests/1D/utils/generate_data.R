## loadings_true generator
fourier_basis <- function(locs, i) {
  k <- c(1, 1, 2, 2)[i]
  a <- c(0, 1, 0, 1)[i]
  return((1-a)*sqrt(2)*cos(2*pi*k*locs) + a*sqrt(2)*sin(2*pi*k*locs))
}

build_psi <- function(locations, grid_HR){
  # Create a logical matrix of dimension length(locations) x length(grid_HR),
  # where TRUE indicates that locations[i] == grid_HR[j].
  # Then convert the logical matrix (TRUE/FALSE) to numeric (1/0).
  psi_matrix <- outer(locations, grid_HR, FUN = "==") + 0
  
  return(psi_matrix)
}

## data generator
generate_data <- function(grid_HR, locs = NULL,
                          x_basis_true_generator = NULL,
                          n_stat_units = 50,
                          n_comp = 3,
                          score_dist = "unif",
                          NSR_X_lc = 0.5,
                          NSR_Y = 0.5,
                          seed = 1,
                          VERBOSE = FALSE) {
  
  # defaults & dimensions -------------------------------------------------------
  
  ## set defaults
  if (is.null(x_basis_true_generator)) {
    x_basis_true_generator <- fourier_basis
  }
  G_true <- matrix(c(1,1,1,1), ncol = n_comp)
  
  ## HR grid
  Psi_FEM <- build_psi(locations, grid_HR)
  # spline_basis_HR <- create.bspline.basis(rangeval = range(grid_HR), breaks = grid_HR)
  n_knots_HR <- length(grid_HR)
  
  ## dimensions
  n_locs <- length(locs)
  n_resp <- nrow(G_true)
  
  
  # pre-generation --------------------------------------------------------------
  
  ## generating the x basis functions
  F_true <- matrix(0, nrow = n_knots_HR, ncol = n_comp)
  # F_true_locs <- matrix(0, nrow = n_locs, ncol = n_comp)
  for (i in 1:n_comp) {
    F_true[, i] <- x_basis_true_generator(grid_HR, i)
    # F_true_locs[, i] <- x_basis_true_generator(locs, i)
  }
  
  ## scores sd
  sigma_S <- exp(seq(0, log(0.25), length = n_comp))
  
  ## sample scores
  set.seed(seed)
  if(score_dist == "norm"){
    Sampled <- mvrnorm(n_stat_units, mu = rep(0, n_comp), diag(sigma_S^2))
  } else if(score_dist == "unif") {
    Sampled <- NULL
    for(i in 1:n_comp) {
      set.seed(seed*100000 + i)
      Sampled <- cbind(Sampled, runif(n_stat_units, min = -2*sigma_S[i], max = 2*sigma_S[i]))
    }
  }
  
  S_true <- scale(Sampled, scale = FALSE)
  
  ## X true data
  X_true <- S_true %*% t(F_true)
  
  ## Y true data
  Y_true <- S_true %*% t(G_true)
  
  
  # PLS -------------------------------------------------------------------------
  
  results_PLS <- PLS(Y_true, X_true, n_comp, center = FALSE, mode = "PLS-R")
  if(norm(Y_true - results_PLS$Y_hat[[n_comp]], "I") > 1e-9 || norm(X_true - results_PLS$X_hat[[n_comp]], "I") > 1e-9){
    stop("The reconstruction error is too big!")
  }
  
  ## external
  Beta_true <- results_PLS$Beta_hat[[n_comp]]
  
  ## internal
  V_true <- results_PLS$Y_space_directions
  W_true <- results_PLS$X_space_directions
  T_true <- results_PLS$X_latent_scores
  C_true <- results_PLS$Y_loadings
  R_true <- results_PLS$X_loadings
  
  
  # locations ------------------------------------------------------------------
  
  ## external
  X_true_locs <- X_true %*% t(Psi_FEM)
  Beta_true_locs <- Psi_FEM %*% Beta_true
  
  ## internal
  R_true_locs <- Psi_FEM %*% R_true
  W_true_locs <- Psi_FEM %*% W_true
  
  
  # noise ----------------------------------------------------------------------
  
  X_true_locs_lc <- T_true[, n_comp] %*% t(R_true_locs[, n_comp])
  
  # NSR <- noise_power / signal_power
  # noise_power <- NSR * signal_power 
  # => sigma_noise <- sqrt(noise_power) <- sqrt(NSR * signal_power)
  
  ## X
  signal_X_power <- mean(X_true_locs^2)
  signal_X_locs_lc_power <- mean(X_true_locs_lc^2)
  sigma_X_noise <- sqrt(NSR_X_lc * signal_X_locs_lc_power)
  NSR_X <- sigma_X_noise^2 / signal_X_power
  
  ## Y
  Y_powers <- colMeans(Y_true^2)
  sigma_Y_noise <- sqrt(NSR_Y * Y_powers)
  
  ## sampling
  set.seed(seed*10000)
  Sampled <- mvrnorm(n_stat_units, mu = rep(0, n_locs + n_resp),
                     diag(c(rep(sigma_X_noise^2, n_locs), sigma_Y_noise^2)))
  EE <- scale(as.matrix(Sampled[, 1:(n_locs)], ncol = n_locs), scale = FALSE)
  FF <- scale(as.matrix(Sampled[, (n_locs + 1):(n_locs + n_resp)], ncol = n_resp), scale = FALSE)
  
  
  # observed data --------------------------------------------------------------
  
  X <- X_true_locs + EE
  Y <- Y_true + FF
  
  return(list(
    ## dimensions
    dimensions = list(
      n_stat_units = n_stat_units,
      n_comp = n_comp,
      n_locs = n_locs,
      n_resp = n_resp
    ),
    ## data
    data = list(
      Y = Y,
      X = X
    ),
    expected_results = list(
      ## expected results: reconstruction
      Y_true = Y_true,
      X_true = X_true,
      X_true_locs = X_true_locs,
      ## expected results: decomposition
      Beta_true = Beta_true,
      Beta_true_locs = Beta_true_locs,
      Y_space_directions_true = V_true,
      X_space_directions_true = W_true,
      X_space_directions_true_locs = W_true_locs,
      X_latent_scores_true = T_true,
      X_loadings_true = R_true,
      X_loadings_true_locs = R_true_locs,
      Y_loadings_true = C_true
    ),
    ## computed
    sigma_S = sigma_S,
    sigma_Y_noise = sigma_Y_noise,
    sigma_X_noise = sigma_X_noise,
    NSR_X = NSR_X
  ))
}
