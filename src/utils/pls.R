PLS <- function(Y, X, n_comp = 4, center = TRUE, mode = "PLS-R") {
  
  if (mode != "PLS-A" & mode != "PLS-SB" & mode != "PLS-R") {
    stop("Error: mode must be either 'PLS-A', 'PLS-SB', or 'PLS-R'")
  }
  
  #### Initialization  ####
  
  # Centering
  if(center) {
    X_c = scale(X, scale = FALSE)
    X_mean = attr(X_c, "scaled:center")
    Y_c <- scale(Y, scale = FALSE)
    Y_mean = attr(Y_c, "scaled:center")
  } else {
    X_c = X
    X_mean = rep(0, ncol(X))
    Y_c = Y
    Y_mean = rep(0, ncol(Y))
  }
  
  # Extracting dimensions
  N <- nrow(X_c) #N
  L <- ncol(Y_c) #L
  S <- ncol(X_c) #S
  
  # X and Y directions
  X_space_directions <- matrix(0, nrow = S, ncol = n_comp)
  Y_space_directions <- matrix(0, nrow = L, ncol = n_comp)
  # X and Y scores/components
  X_latent_scores <- matrix(0, nrow = N, ncol = n_comp)
  Y_latent_scores <- matrix(0, nrow = N, ncol = n_comp)
  # X and Y loadings
  X_loadings <- matrix(0, nrow = S, ncol = n_comp)
  Y_loadings <- matrix(0, nrow = L, ncol = n_comp)
  
  # B matrix
  B <- matrix(0, nrow = n_comp, ncol = n_comp)
  
  # Mean values in matrix form
  X_MEAN = matrix(X_mean, nrow = N, ncol = S, byrow = TRUE)
  Y_MEAN = matrix(Y_mean, nrow = N, ncol = L,  byrow = TRUE)
  
  # Matrices for intermediate steps
  XX <- list()
  YY <- list()
  
  # X and Y residual matrices
  EE <- X_c
  FF <- Y_c
  XX[[1]] <- EE
  YY[[1]] <- FF
  
  ## room for results
  Y_hat <- list()
  X_hat <- list()
  Beta_hat <- list()
  
  #### Main loop ####
  
  for(h in 1:n_comp){
    
    ### Step 1: covariance maximization ### 
    
    # Weight extraction
    C.YX <- t(FF) %*% EE
    SVD <- svd(C.YX, nu = 1, nv = 1)
    X_space_directions[,h] <- SVD$v  
    Y_space_directions[,h] <- SVD$u 
    
    # Score/component computation
    X_latent_scores[,h] <- EE %*% X_space_directions[,h] 
    Y_latent_scores[,h] <- FF %*% Y_space_directions[,h] 
    
    ### Step 2-3: Loading computation ### 
    
    # Loading computation
    tt = sum(X_latent_scores[,h]^2)
    uu = sum(Y_latent_scores[,h]^2)
    X_loadings[,h] <- t(EE) %*% X_latent_scores[,h] / (tt)
    if(mode == "PLS-R") {
      Y_loadings[,h] <- t(FF) %*% X_latent_scores[,h] / (tt)
    } else if(mode == "PLS-A") {
      Y_loadings[,h] <- t(FF) %*% Y_latent_scores[,h] / (uu)
    } else if(mode == "PLS-SB") {
      # Loadings are set to be the directions
      X_loadings[,h] <- X_space_directions[,h] 
      Y_loadings[,h] <- Y_space_directions[,h] 
    }
    
    ### Step 4: Regression (only in PLS_12) + Deflation ### 
    
    # X-deflation
    EE <- EE - X_latent_scores[,h] %*% t(X_loadings[,h])
    
    # Y-deflation
    if(mode=="PLS-R"){
      B[h,h] <- t(Y_latent_scores[,h]) %*% X_latent_scores[,h] / (tt) 
      FF <- FF - X_latent_scores[,h] %*% (B[h,h] * t(Y_space_directions[,h]))
    }else{
      FF <- FF - Y_latent_scores[,h] %*% t(Y_loadings[,h])  
    }
    
    # Saving intermediate results
    XX[[h+1]] = EE
    YY[[h+1]] = FF
    
    # X estimate
    X_hat[[h]] <- X_latent_scores[,1:h] %*% t(X_loadings[,1:h]) + X_MEAN
    
    # Y (and Beta) estimates
    if(mode == "PLS-R"){
      W_h <- X_space_directions[,1:h]
      R_h <- X_loadings[,1:h]
      V_h <- Y_space_directions[,1:h]
      if(L == 1) V_h <- t(V_h)
      Beta_hat[[h]] <- W_h %*% solve(t(R_h) %*% W_h , B[1:h,1:h] %*% t(V_h))
      Y_hat[[h]] <- X_c %*% Beta_hat[[h]] + Y_MEAN
    }else{
      Y_hat[[h]] <- Y_latent_scores[,1:h] %*% t(Y_loadings[,1:h]) + Y_MEAN
    }
  }
  
  return(list(X_space_directions = X_space_directions,
              Y_space_directions = Y_space_directions,
              X_latent_scores = X_latent_scores,
              Y_latent_scores = Y_latent_scores,
              X_loadings = X_loadings,
              Y_loadings = Y_loadings,
              Beta_hat = Beta_hat,
              X_mean = X_mean,
              Y_mean = Y_mean,
              Y_hat = Y_hat,
              X_hat = X_hat))
}