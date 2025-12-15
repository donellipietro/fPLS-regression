
MV_PLS_wrapped <- function(data, n_comp = 4, center = TRUE, mode = "PLS-R"){
  
  model_PLS <- PLS(data$Y, data$X, n_comp = n_comp, center = center, mode = mode)
  
  return(list(
    model = NULL,
    MODE = mode,
    results = list(
      ## centering
      Y_mean = model_PLS$Y_mean,
      X_mean_locs = model_PLS$X_mean,
      X_mean = NULL,
      ## space directions
      Y_space_directions = model_PLS$Y_space_directions,
      X_space_directions_locs = model_PLS$X_space_directions,
      X_space_directions = NULL,
      ## latent scores
      Y_latent_scores = model_PLS$Y_latent_scores,
      X_latent_scores = model_PLS$X_latent_scores,
      ## loadings
      Y_loadings = model_PLS$Y_loadings,
      X_loadings_locs = model_PLS$X_loadings,
      X_loadings = NULL,
      ## reconstruction
      Y_hat = model_PLS$Y_hat,
      X_hat_locs = model_PLS$X_hat,
      X_hat = NULL,
      ## Beta
      Beta_hat_locs = model_PLS$Beta_hat,
      Beta_hat = NULL
    )
  ))
}


B_fPLS_wrapped <- function(data, n_comp = 4, center = FALSE, mode = "fPLS-R", presmoothing = FALSE, penalized = FALSE, lambda_grid = NULL){
  
  if(mode != "fPLS-R") stop("Error, this mode is not implemented!")
  
  ## data
  Y <- data$Y
  X <- data$X
  
  ## nodes and locations
  nodes <- as.vector(data$domain$nodes)
  locations <- as.vector(data$locations)
  
  ## dimensions
  n_stat_units <- nrow(data$Y)
  n_resps <- ncol(data$Y)
  n_nodes <- length(nodes)
  n_locs <- length(locations)
  
  # ## spline utils
  basisobj <- create.bspline.basis(rangeval = range(nodes), breaks = nodes)
  
  Psi <- eval.basis(locations, basisobj)
  R0 <- fda::bsplinepen(basisobj = basisobj, Lfdobj = 0)
  
  if(center) stop("Not implemented, you need to center the data by yourself!")
  Xc <- X
  Yc <- Y
  
  model_B_fPLS <- B_fPLS(Yc, Xc, nodes, locations, presmoothing, penalized, lambda_grid, basisobj, Psi, R0)
  TT <- model_B_fPLS$TT
  W <- model_B_fPLS$W
  V <- model_B_fPLS$V
  R <- model_B_fPLS$R
  C <- model_B_fPLS$C
  Beta_hat <- model_B_fPLS$Beta_hat
  
  ## at locations
  R_locs <- Psi %*% R
  W_locs <- Psi %*% W
  
  ## at nodes
  Psi_nodes <- eval.basis(as.vector(domain$nodes), basisobj)
  R_nodes <- Psi_nodes %*% R
  W_nodes <- Psi_nodes %*% W
  
  # results
  Y_hat_B_fPLS  <- list()
  X_hat_B_fPLS  <- list()
  X_hat_locs_B_fPLS  <- list()
  X_hat_nodes_B_fPLS  <- list()
  B_hat_B_fPLS  <- list()
  B_hat_locs_B_fPLS  <- list()
  B_hat_nodes_B_fPLS  <- list()
  for(h in 1:n_comp) {
    Y_hat_B_fPLS[[h]] <- as.matrix(model_B_fPLS$fitted.values[,,h])
    X_hat_B_fPLS[[h]] <- TT[,1:h] %*% t(R[,1:h])
    X_hat_locs_B_fPLS[[h]] <- TT[,1:h] %*% t(R_locs[,1:h])
    X_hat_nodes_B_fPLS[[h]] <- TT[,1:h] %*% t(R_nodes[,1:h])
    B_hat_B_fPLS[[h]] <- as.matrix(Beta_hat[,,h])
    B_hat_locs_B_fPLS[[h]] <- Psi %*% as.matrix(Beta_hat[,,h])
    B_hat_nodes_B_fPLS[[h]] <- Psi_nodes %*% as.matrix(Beta_hat[,,h])
  }
  
  return(list(
    MODE = mode,
    basisobj = basisobj,
    results = list(
      ## centering
      Y_mean = NULL,
      X_mean_locs = NULL,
      X_mean = NULL,
      ## space directions
      Y_space_directions = t(rep(1, n_comp)),
      X_space_directions_locs = as.matrix(W_locs),
      X_space_directions_nodes = as.matrix(W_nodes),
      X_space_directions = as.matrix(W),
      ## latent scores
      Y_latent_scores = as.matrix(V),
      X_latent_scores = as.matrix(TT),
      ## loadings
      Y_loadings =  as.matrix(C),
      X_loadings_locs = as.matrix(R_locs),
      X_loadings_nodes = as.matrix(R_nodes),
      X_loadings = as.matrix(R),
      ## reconstruction
      Y_hat = Y_hat_B_fPLS,
      X_hat_locs = X_hat_locs_B_fPLS,
      X_hat_nodes = X_hat_nodes_B_fPLS,
      X_hat = X_hat_B_fPLS,
      ## Beta
      Beta_hat_locs = B_hat_locs_B_fPLS,
      Beta_hat_nodes = B_hat_nodes_B_fPLS,
      Beta_hat = B_hat_B_fPLS
    )
  ))
  
}


B_fPLS <- function(Y, X, nodes, locations, presmoothing, penalized, lambda_grid, basisobj, Psi, R0) {
  
  ## dimensions
  n_stat_units <- nrow(data$Y)
  n_resps <- ncol(data$Y)
  n_nodes <- length(nodes)
  n_locs <- length(locations)
  
  ## penalty
  dorder <- 2
  P <- fda::bsplinepen(basisobj = basisobj, Lfdobj = dorder)
  
  ## basis representation
  if(is.null(lambda_grid) || !presmoothing) {
    bsplineX <- fda::Data2fd(locations, y = Matrix::t(X), basisobj = basisobj)
    A <- Matrix::t(bsplineX[["coefs"]])
    best_lambda <- 0
  } else {
    
    # Select the overall best lambda
    lambdas <- lapply(1:nrow(X), function(i) {
      y <- X[i, ]
      
      # Apply smooth.basis with a grid of lambda
      gcv_results <- sapply(lambda_grid$space, function(lambda) {
        fdPar_obj <- fdPar(basisobj, Lfdobj = 2, lambda = lambda)  # Lfdobj = 2 for roughness penalty
        fit <- smooth.basis(locations, y, fdPar_obj)
        sum(fit$gcv)  # Return the total GCV score
      })
      
      # Select the lambda that minimizes GCV
      best_lambda <- lambda_grid$space[which.min(gcv_results)]
      
      # Save the best lambda to the vector
      best_lambda
      
    })
    
    # select the overall optimal lambda
    log10_lambdas <- log10(unlist(lambdas))
    density_estimate <- density(log10_lambdas)
    mode_value <- density_estimate$x[which.max(density_estimate$y)]
    
    best_lambda <- 10^mode_value
    
    # Plot the density and highlight the mode
    # plot(density_estimate, main = "Density and Mode", xlab = "Value", ylab = "Density")
    # abline(v = mode_value, col = "red", lwd = 2, lty = 2)  # Add a vertical line for the mode
    # legend("topleft", legend = paste("Mode =", round(mode_value, 2)), col = "red", lty = 2, lwd = 2)
    
    A <- t(solve(t(Psi) %*% Psi + best_lambda*P, t(Psi) %*% t(X)))
  }
  
  if(penalized){
    tuning <- tune_lambda_cv(Y, A, lambda_grid$space, n_comp, K = 5, R0, P)
    returned <- B_fPLS_core(Y, A, tuning$best_lambda, n_comp, R0, P)
    mvpls_model <- returned$model
    L <- returned$L
    inv_t_L <- returned$inv_t_L
  } else {
    returned <- B_fPLS_core(Y, A, 0, n_comp, R0, P)
    mvpls_model <- returned$model
    L <- returned$L
    inv_t_L <- returned$inv_t_L
  }
  
  ## Get MV-model components:
  TT <- as.matrix(mvpls_model[["scores"]])
  V <-  as.matrix(mvpls_model[["Yscores"]])
  R_tilde <- as.matrix(mvpls_model[["loadings"]] )
  C <-  as.matrix(mvpls_model[["Yloadings"]])
  W_tilde <- as.matrix(mvpls_model[["loading.weights"]] )
  
  ## coefficient function:
  Beta_hat <- array(NA, dim = c(n_nodes+2, n_resps, n_comp))
  for (h in 1:n_comp) {
    tmp_tilde <- W_tilde[, 1:h] %*% solve(t(R_tilde[, 1:h]) %*% W_tilde[, 1:h])
    beta_coeff_h <- inv_t_L %*% tmp_tilde %*% C[, 1:h]
    Beta_hat[ , , h] <- solve(t(Psi) %*% Psi + best_lambda*P, R0 %*% beta_coeff_h)
  }
  
  ## back-transform X loadings
  R <- solve(R0) %*% L %*% R_tilde
  
  ## back-transform X directions
  W <- array(NA, dim = c(n_nodes+2, n_comp))
  # W[, 1] <- solve(t(Psi) %*% Psi) %*% R0 %*% inv_t_L %*% W_tilde[, 1]
  # for (h in 2:n_comp) {
  #   tmp <- W[, 1:(h-1)] %*% solve(t(R[, 1:(h-1)]) %*% W[, 1:(h-1)])
  #   W[, h] <- ((diag(n_nodes+2) + tmp %*% t(R[, 1:(h-1)])) %*% solve(t(Psi) %*% Psi) %*% R0 %*% inv_t_L - tmp %*% t(R_tilde[, 1:(h-1)])) %*% W_tilde[, h]
  # }
  
  return(list(
    TT = TT,
    W = W,
    V = V,
    R = R,
    C = C,
    best_lambda = best_lambda,
    Beta_hat = Beta_hat,
    fitted.values = mvpls_model$fitted.values
  ))
}

B_fPLS_core <- function(Y, A, lambda, n_comp, R0, P) {
  ## transformation matrix
  LtL <- R0 + lambda*P
  L <- expm::sqrtm(LtL)
  inv_t_L <- Matrix::solve(L)
  
  ## data transformation
  data_pls <- A %*% R0 %*% inv_t_L
  
  # PLS model:
  mvpls_model <- pls::plsr(Y ~ data_pls,
                           ncomp =  n_comp,
                           method = "oscorespls",
                           center = FALSE,
                           scale = FALSE)
  return(list(
    model = mvpls_model,
    L = L,
    inv_t_L = inv_t_L
  ))
}


tune_lambda_cv <- function(Y, A, lambda_grid, n_comp, K = 5, R0, P) {
  # Check if the number of rows in Y and A are the same
  if (nrow(Y) != nrow(A)) stop("Y and A must have the same number of rows")
  
  # Create folds
  n <- nrow(Y)
  folds <- sample(rep(1:K, length.out = n))
  
  # Initialize a vector to store the cross-validated errors for each lambda
  cv_errors <- numeric(length(lambda_grid))
  
  # Loop over each lambda
  for (l in seq_along(lambda_grid)) {
    lambda <- lambda_grid[l]
    fold_errors <- numeric(K)  # Store errors for each fold
    
    # Perform K-fold CV
    for (k in 1:K) {
      # Split into training and validation sets
      train_idx <- which(folds != k)
      val_idx <- which(folds == k)
      
      Y_train <- Y[train_idx, , drop = FALSE]
      A_train <- A[train_idx, , drop = FALSE]
      Y_val <- Y[val_idx, , drop = FALSE]
      A_val <- A[val_idx, , drop = FALSE]
      
      # Fit the model on the training data
      model <- B_fPLS_core(Y_train, A_train, lambda, n_comp, R0, P)$model
      
      # Extract the fitted values for the validation set
      Y_pred <- predict(model, newdata = A_val, n_comp = 4)
      
      # Compute the validation error (e.g., mean squared error)
      fold_errors[k] <- mean((Y_val - Y_pred[,,4])^2)
    }
    
    # Compute the average cross-validated error for this lambda
    cv_errors[l] <- mean(fold_errors)
  }
  
  # Find the best lambda
  best_lambda <- lambda_grid[which.min(cv_errors)]
  
  # Return results
  return(list(
    best_lambda = best_lambda,
    cv_errors = cv_errors,
    lambda_grid = lambda_grid
  ))
}


old_B_fPLS_wrapped <- function(data, n_comp = 4, center = TRUE, mode = "fPLS-R", lambda_grid = 0){
  
  if(mode != "fPLS-R") stop("Error, this mode is not implemented!")
  
  nodes <- data$domain$nodes
  
  N <- nrow(data$Y)
  k_folds <- 5
  
  if(ncol(nodes) == 1) {
    method = "fpls_bs"
    
    # Ruppert's law:  nbasis = nbreaks + norder - 2  and norder = degree + 1
    n_breaks <- min(round(nrow(nodes)/4), 40)
    n_basis <- n_breaks + (3+1) - 2
    
    # B-spline basis:
    basisobj <- create.bspline.basis(rangeval = range(nodes),
                                     nbasis = n_basis)
    
    # Evaluate basis functions:
    Psi <- fda::eval.basis(evalarg = nodes[,1],
                           basisobj = basisobj,
                           Lfdobj=0, returnMatrix=TRUE)
    R0 <- NULL
    
  } else if(ncol(nodes) == 2) {
    method = "fpls_tps"
    
    basisobj <- NULL
    
    ## initialization
    n_basis <- 10
    gam_fit <- mgcv::gam(data$X[1, ] ~ s(nodes[ , 1],
                                         nodes[ , 2],
                                         bs = "tp",
                                         k = n_basis_tps))
    # Evaluate basis functions:
    Psi <- stats::model.matrix(gam_fit)
    # Matrix of inner products (mass):
    R0 <- matrix(NA, nrow = ncol(Psi), ncol = ncol(Psi))
    # Numerical approx. of the inner products:
    for (ii in 1:nrow(R0)) {
      for (jj in ii:ncol(R0)) {
        df <- as.data.frame(nodes)
        df$z = as.numeric(Psi[, ii]*Psi[, jj])
        R0[ii,jj] <-  penR1FPLS:::getVolume(df)
      }
    }
    R0[lower.tri(R0)] <- R0[upper.tri(R0)]
  }
  
  ## fit
  cv_pTPS <- cv_unique_par(X = data$X,
                           Y = data$Y,
                           center = center,
                           nodes = nodes,
                           argvals = data$locations,
                           nbasis = n_basis,
                           penalty_vec = lambda_grid,
                           basisobj = basisobj,
                           ncomp = n_comp,
                           folds = k_folds,
                           method = method,
                           verbose = FALSE,
                           stripped = FALSE,
                           R0 = R0)
  model_B_fPLS <- cv_pTPS$final_model
  
  # results
  Y_hat_B_fPLS  <- list()
  X_hat_B_fPLS  <- list()
  B_hat_B_fPLS  <- list()
  X_mean_B_fPLS  <- model_B_fPLS$X_mean
  Y_mean_B_fPLS  <- model_B_fPLS$Y_mean
  C_hat_B_fPLS  <- as.matrix(model_B_fPLS$C)
  TT <- as.matrix(model_B_fPLS$TT)
  for(h in 1:n_comp) {
    Y_hat_B_fPLS[[h]] <- as.matrix(model_B_fPLS$fitted.values[,,h])
    X_hat_B_fPLS[[h]] <- TT[,1:h] %*% t(C_hat_B_fPLS[,1:h]) + rep(1, N) %*% t(X_mean_B_fPLS)
    B_hat_B_fPLS[[h]] <- as.matrix(model_B_fPLS$coefficient_function[,,h])
  }
  
  return(list(
    model = model_B_fPLS,
    MODE = mode,
    Psi = Psi,
    R0 = R0,
    results = list(
      ## centering
      Y_mean = model_B_fPLS$Y_mean,
      X_mean_locs = model_B_fPLS$X_mean,
      X_mean = NULL,
      ## space directions
      Y_space_directions = t(rep(1, n_comp)),
      X_space_directions_locs = as.matrix(model_B_fPLS$W),
      X_space_directions = NULL,
      ## latent scores
      Y_latent_scores = as.matrix(model_B_fPLS$V),
      X_latent_scores = as.matrix(model_B_fPLS$TT),
      ## loadings
      Y_loadings = t(as.numeric(model_B_fPLS$D)),
      X_loadings_locs = as.matrix(model_B_fPLS$C),
      X_loadings = NULL,
      ## reconstruction
      Y_hat = Y_hat_B_fPLS,
      X_hat_locs = X_hat_B_fPLS,
      X_hat = NULL,
      ## Beta
      Beta_hat_locs = B_hat_B_fPLS,
      Beta_hat = NULL
    )
  ))
  
}



