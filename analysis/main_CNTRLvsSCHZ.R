# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Application Brain: Analysis %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
graphics.off()


## libraries ----

# fdaPDE
library(fdaPDE)
library(fdaPDE2)

# classification
library(pROC)

# matrix
library(Matrix)


## sources ----

source("src/utils/directories.R")
source("src/utils/domain_and_locations.R")
source("src/utils/meshes.R")
source("src/utils/plots.R")
source("src/utils/wrappers.R")
source("src/utils/errors.R")
source("src/utils/pls.R")
source("analysis/utils/vtu.R")


## global variables ----

# directories
directory.images <- "images/brain/"
path_results <- "results/brain/"
mkdir(c(path_results, directory.images))

# code flow control
RUN <- list()
RUN[["Mean Estimation - ColMean"]] <- FALSE
RUN[["Mean Estimation - fMean"]] <- FALSE
RUN[["MV-PLS"]] <- FALSE
RUN[["fPLS"]] <- FALSE
RUN[["fPLS - GCV"]] <- FALSE


## data ----

# Mesh
FEMbasis <- readRDS("data/brain/FEMbasis.Rds")
mesh <- FEMbasis$mesh

# Domain
mesh_data <- list(
  nodes = mesh$nodes,
  edges = mesh$faces,
  elements = mesh$tetrahedrons,
  neigh = mesh$neighbors,
  boundary = mesh$nodesmarkers
)
domain <- fdaPDE2::Mesh(mesh_data)

# Brain connectivity maps
X <- readRDS("data/brain/X_imputed.rds")


sub1 <- FEM(as.numeric(X[1, ]), FEMbasis)
mkdir(paste0(path_results, "data/"))
write.vtu(sub1, paste(path_results, "data/subject_1.vtu", sep = ""))


# Responses
Ycl <- readRDS("data/brain/Ycl.rds")
Ycl[Ycl == "CONTROL"] <- "CTRL" 
Ycl <- factor(Ycl)
levels(Ycl)

# Names
names.groups <- c("Control", "Schz")

# Colors
colors <- c("#206400", '#AE017E')
names(colors) <- levels(Ycl)

# Y01 <- readRDS("data/Y01.rds")
# Yp <- readRDS("data/Yp.rds")


## support variables ----

# Indexes groups
index.CONTROL <- Ycl == "CTRL"
index.SCHZ <- Ycl == "SCHZ"

# Check
any(!(index.CONTROL | index.SCHZ))
any((index.CONTROL & index.SCHZ))

# Response re-scaling
p1 <- sum(index.CONTROL)/length(Ycl)
p2 <- sum(index.SCHZ)/length(Ycl)
Yp <- rep(0, length(Ycl))
Yp[index.CONTROL] <- sqrt(p2/p1)
Yp[index.SCHZ] <- -sqrt(p1/p2)

# Ycl: CONTROL                      vs    SCHZ
# Yp:  -sqrt(P{CONTROL}/P{SCHZ})    vs    sqrt(P{SCHZ}/P{CONTROL})


## mean estimation ----

### column-wise mean ----

ColMean <- function(X, domain = NULL, lambdas = NULL, verbose = FALSE){
  
  
  fitted <- fitted_nodes <- colMeans(X)
  
  return(list(fitted = fitted,
              fitted_nodes = fitted_nodes))
  
}

name <- "ColMean"
method <- ColMean

if(RUN[["Mean Estimation - ColMean"]]){
  source("analysis/templates/mean_estimation.R")
}


### functional mean ----

fMean <- function(X, domain = NULL, lambdas = NULL, verbose = FALSE){
  
  data_centering <- fdaPDE2::functional_data(
    domain = domain,
    X = X
  )
  
  model <- fdaPDE2::fCentering(data_centering, calibrator = gcv(lambdas, seed = 0))
  model$fit()
  
  fitted <- fitted_nodes <- model$results$mean
  
  return(list(fitted = fitted,
              fitted_nodes = fitted_nodes))
  
}

name <- "fMean"
method <- fMean

lambdas <- fdaPDE2::hyperparameters(10^seq(-2, -2, by = 0.1))

if(RUN[["Mean Estimation - fMean"]]){
  source("analysis/templates/mean_estimation.R")
}


## sensitivity analysis (on lambda) ----


### data ----

## load mean
load(paste("results/brain/mean_estimation/mean_estimation_fMean.RData", sep = ""))

# center data

name_test <- "overall_mean"
Xc <- X - rep(1, nrow(X)) %*% t(means$overall$fitted)

name_test <- "group_means"
Xc <- X
Xc[index.CONTROL] <- Xc[index.CONTROL] - rep(1, sum(index.CONTROL)) %*% t(means$CONTROL$fitted)
Xc[index.SCHZ] <- Xc[index.SCHZ] - rep(1, sum(index.SCHZ)) %*% t(means$SCHZ$fitted)

# name_test <- "controll_mean"
# Xc <- X - rep(1, nrow(X)) %*% t(means$CONTROL$fitted)


### MV-PLS ----

if(RUN[["MV-PLS"]]){
  results_PLS <- PLS(matrix(Yp, ncol = 1), Xc, 3, center = FALSE)
  save(results_PLS, file = paste(path_results, "results_PLS.RData", sep = ""))
}
load( paste(path_results, "results_PLS.RData", sep = ""))


### fPLS ----

# directories
path_results <- paste(path_results, "sensitivity_analysis/", sep = "")
mkdir(c(path_results))

# directories
path_results <- paste(path_results, name_test, "/", sep = "")
mkdir(c(path_results))

lambda_vect <- 10^seq(-9,2, by = 0.2)
save(mesh_data, Yp, Xc, lambda_vect, file = paste(path_results, "data.RData", sep = ""))
if(RUN[["fPLS"]]){
  system(paste("bash run_analysis_parallel.sh", "sensitivity_analysis.R", name_test, length(lambda_vect)))
}

### post-processing ----

mass <- Diagonal(x = rowSums(fdaPDE:::CPP_get.FEM.Mass.Matrix(FEMbasis)))
stiff <- fdaPDE:::CPP_get.FEM.Stiff.Matrix(FEMbasis)
P <- stiff %*% solve(mass, stiff)
save(P, file = "data/brain/P.RData")
load("data/brain/P.RData")

compute_smoothness <- function(f){
  as.numeric(t(f) %*% P %*% f)
}

## load data
Y1_all <- NULL
Y2_all <- NULL
Y3_all <- NULL
optimal_coords_1 <- NULL
optimal_coords_2 <- NULL
optimal_coords_3 <- NULL
smoothness <- NULL
for(i in 1:length(lambda_vect)){
  
  ## load data
  load(paste(path_results, paste("results_fPLS_", sprintf("%0*d", 3, i),".RData", sep=""), sep = ""))
  Y1_all <- rbind(Y1_all, t(results_fPLS$Y_hat[[1]]))
  Y2_all <- rbind(Y2_all, t(results_fPLS$Y_hat[[2]]))
  Y3_all <- rbind(Y3_all, t(results_fPLS$Y_hat[[3]]))
  
  # compute optimal threshold
  roc_curve <- roc(Ycl, as.vector(results_fPLS$Y_hat[[1]]))
  optimal_coords_1 <- rbind(optimal_coords_1, coords(roc_curve, "best", ret = c("threshold", "sensitivity", "specificity")))
  
  # compute optimal threshold
  roc_curve <- roc(Ycl, as.vector(results_fPLS$Y_hat[[2]]))
  optimal_coords_2 <- rbind(optimal_coords_2, coords(roc_curve, "best", ret = c("threshold", "sensitivity", "specificity")))
  
  # compute optimal threshold
  roc_curve <- roc(Ycl, as.vector(results_fPLS$Y_hat[[3]]))
  optimal_coords_3 <- rbind(optimal_coords_3, coords(roc_curve, "best", ret = c("threshold", "sensitivity", "specificity")))
  
  # # smoothness
  # smoothness <- rbind(
  #   smoothness, 
  #   data.frame(
  #     r1 = compute_smoothness(results_fPLS$X_loadings[,1]),
  #     r2 = compute_smoothness(results_fPLS$X_loadings[,2]),
  #     r3 = compute_smoothness(results_fPLS$X_loadings[,3])
  #   )        
  # )
  
}

## plots
plot_lambda_dep <- function(Y_all, optimal_coords, h) {
  par(mfrow = c(1,2))
  matplot(NaN, xlim = c(1e-10, 1e2), ylim = range(Y_all), log = "x", main = "CONTROL", ylab = "Y_hat", xaxt="n")
  axis(1, at=10^((-10):2), labels= c("MV", 10^((-9):2)), las = 2)
  grid()
  matplot(lambda_vect, Y_all[, Ycl == "CONTROL"], col = colors["CONTROL"], type = "l", add = TRUE)
  matpoints(rep(1e-10, sum(Ycl == "CONTROL")), results_PLS$Y_hat[[h]][Ycl == "CONTROL", 1], pch = 20, col = colors["CONTROL"])
  matplot(lambda_vect, optimal_coords$threshold, col = "black", lwd = 3, type = "l", add = TRUE)
  matplot(NaN, xlim = c(1e-10, 1e2), ylim = range(Y_all), log = "x", main = "SCHZ", ylab = "Y_hat", xaxt="n")
  axis(1, at=10^((-10):2), labels= c("MV", 10^((-9):2)), las = 2)
  grid()
  matplot(lambda_vect, Y_all[, Ycl == "SCHZ"], col = colors["SCHZ"], type = "l", add = TRUE)
  matpoints(rep(1e-10, sum(Ycl == "SCHZ")), results_PLS$Y_hat[[h]][Ycl == "SCHZ", 1], pch = 20, col = colors["SCHZ"])
  matplot(lambda_vect, optimal_coords$threshold, col = "black", lwd = 3, type = "l", add = TRUE)
  par(mfrow = c(1,1))
  plot(NaN, xlim = range(lambda_vect), ylim = c(0,1), log = "x", xlab = "", ylab = "", main = "Metrics", xaxt="n")
  grid()
  axis(1, at=10^((-9):2), labels= c(10^((-9):2)))
  points(lambda_vect, optimal_coords$sensitivity, type = "l", lwd = 2, col = "darkgreen")
  points(lambda_vect, optimal_coords$specificity, type = "l", lwd = 2, col = "orange")
  legend("bottomleft", c("Sensitivity", "Specificity"), col = c("darkgreen", "orange"), lty = c(1,1), lwd = 2)
}
plot_lambda_dep(Y1_all, optimal_coords_1, 1)
plot_lambda_dep(Y2_all, optimal_coords_2, 2)
plot_lambda_dep(Y3_all, optimal_coords_3, 3)

# plot(NaN, xlim = range(lambda_vect), ylim = range(smoothness), log = "x", xlab = "", ylab = "", main = "Smoothness", xaxt="n")
# grid()
# axis(1, at=10^((-9):2), labels= c(10^((-9):2)))
# points(lambda_vect, smoothness[,1], type = "l", lwd = 2, col = "purple")
# points(lambda_vect, smoothness[,2], type = "l", lwd = 2, col = "blue")
# points(lambda_vect, smoothness[,3], type = "l", lwd = 2, col = "lightblue")
# legend("bottomleft", c("r1", "r2", "r3"), col = c("purple", "blue", "lightblue"), lty = c(1,1,1), lwd = 2)

## decomposition ----


index_selected <- 1 # 41
type <- "multivariate"

## load data
load(paste(path_results, paste("results_fPLS_", sprintf("%0*d", 3, index_selected),".RData", sep=""), sep = ""))
optimal_lambda <- lambda_vect[index_selected]
optimal_treshold <- optimal_coords_3$threshold[index_selected]


if(RUN[["fPLS - GCV"]]){
  domain <- fdaPDE2::Mesh(mesh_data)
  data <- fdaPDE2::functional_regression_data(
    domain = domain,
    Y = Yp,
    X = Xc
  )
  
  lambda_grid <- fdaPDE2::hyperparameters(10^seq(-3, -1, by = 0.2))
  model_fPLS <- fdaPDE2::fPLS(data, center = FALSE)
  model_fPLS$fit(gcv(lambda_grid, seed = 0), n_comp = 4)
  results_fPLS <- model_fPLS$results
  save(results_fPLS, file = paste(path_results, "results_fPLS_gcv.RData", sep = ""))
}

load(paste(path_results, "results_fPLS_gcv.RData", sep = ""))
type <- "fPLS_gcv"
Y_hat <- results_fPLS$Y_hat[[4]]
X_hat <- results_fPLS$X_hat
roc_curve <- roc(Ycl, as.vector(Y_hat))
optimal_coords_3 <- coords(roc_curve, "best", ret = c("threshold", "sensitivity", "specificity"))
optimal_treshold <- optimal_coords_3$threshold


# Explained X
res_list <- vector("list", length = nrow(X) * 3)
index <- 1
for (i in 1:nrow(X)) {
  for (h in 1:4) {
    res_list[[index]] <- data.frame(
      Group = factor(h),
      SubGroup = factor(as.character(Ycl)[i]),
      Score = 1- RMSE(Xc[i, ] - X_hat[[h]][i, ]) / RMSE(Xc[i, ])
    )
    index <- index + 1
  }
}
res <- do.call(rbind, res_list)
colors_new <- setNames(colors, unique(res$SubGroup))
ggplot(res, aes(x = Group, y = Score, fill = SubGroup)) +  
  geom_boxplot() +
  labs(x = "", y = "") +
  scale_fill_manual(values = colors_new) +
  standard_plot_settings()


df <- data.frame(
  s = Yp-Y_hat[, 1],
  group = Ycl
)
ggplot(df, aes(x = s, color = group, fill = group)) +
  geom_density(alpha = 0.5) +
  # match your colors
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  # same x-range, same y-limit as in your base R plot
  # coord_cartesian(xlim = range_all, ylim = c(0, 1)) +
  labs(x = "", y = "") +
  theme_minimal() + 
  geom_vline(xintercept = optimal_treshold, linetype = "dashed")


library(plotly)
generate_basis <- function(M) {
  # Input: A vector v from R^d
  # Output: An orthonormal basis for R^d containing v
  
  v <- M[, 1]
  d <- length(v) # Dimension of the space
  basis <- matrix(0, nrow = d, ncol = d) # To store the basis vectors
  
  # Normalize the initial vector
  basis[, 1] <- v / sqrt(sum(v^2))
  
  # Generate additional vectors orthogonal to the initial vector
  for (i in 2:d) {
    # Create a random vector
    new_vec <- M[, i]
    
    # Orthogonalize against all previous basis vectors
    for (j in 1:(i - 1)) {
      new_vec <- new_vec - sum(basis[, j] * new_vec) * basis[, j]
    }
    
    # Normalize the new vector
    basis[, i] <- new_vec / sqrt(sum(new_vec^2))
  }
  
  return(basis)
}
plot.latent_space <- function(latent_scores, plane = NULL) {
  
  # Support vectors
  # index.support_vectors <- svm_model$index
  # support_vectors <- latent_scores[index.support_vectors, ]
  # support_vectors$group <- "Support vectors"
  # data <- rbind(latent_scores, support_vectors)
  
  # Scatter-plot latent scores
  fig <- plot_ly(latent_scores,
                 x = ~t1, y = ~t2, z = ~t3,
                 color = ~group, 
                 colors = c("#206400", "#AE017E"), 
                 marker = list(size = 5)) %>%
    add_markers() %>%
    layout(scene = list(
      xaxis = list(title = bquote(t[1])),
      yaxis = list(title = bquote(t[2])),
      zaxis = list(title = bquote(t[3]))
    ),
    title = "Latent space")
  
  if(!is.null(plane)){
    
    # Decision boundary coefficients
    intercept <- plane[1]
    a <- plane[2:4]
    
    # Mesh-grid for the decision boundary plane
    t1_range <- seq(min(latent_scores$t1), max(latent_scores$t1), length.out = 20)
    t2_range <- seq(min(latent_scores$t2), max(latent_scores$t2), length.out = 20)
    meshgrid <- expand.grid(t1 = t1_range, t2 = t2_range)
    meshgrid$t3 <-  -(a[1] * meshgrid$t1 + a[2] * meshgrid$t2 + intercept) / a[3]
    
    
    # Plot the decision boundary plane
    fig <- fig %>% add_trace(data = meshgrid,
                             x = ~t1, y = ~t2, z = ~t3,
                             type = "mesh3d",
                             intensity = rep(0.5, nrow(meshgrid)), # Single color mapped via intensity
                             opacity = 0.5,
                             colorscale = list(c(0, "grey"), c(1, "grey")), # Uniform grey colorscale
                             showscale = FALSE,
                             name = 'Boundary plane',
                             inherit = FALSE)
  }
  
  return(fig)
}

## Y loadings
C <- results_fPLS$Y_loadings[,1:3]

## Latent space
TT <- data.frame(results_fPLS$X_latent_scores)[,1:3]
colnames(TT) <- c("t1", "t2", "t3")
TT$group[index.CONTROL] <- 'Control'
TT$group[index.SCHZ] <- 'Schz'
TT$group <- as.factor(TT$group)

plot.latent_space(TT, c(-optimal_treshold, C))


## Transformed latent space
Q <- generate_basis(cbind(C, c(0, 1, 0), c(0, 0, 1)))
SS <- as.matrix(TT[,1:3]) %*% Q + rep(1, nrow(X)) %*% t(c(- optimal_treshold / (t(Q) %*% C)[1], 0, 0))
colnames(SS) <- c("s1", "s2", "s3")
pairs(TT[,1:3], col = colors[Ycl])
pairs(SS[,1:3], col = colors[Ycl], asp = 1)
# Q <- generate_basis(cbind(t(C), c(0, 1, 0, 0), c(0, 0, 1, 0), c(0, 0, 0, 1)))
# SS <- as.matrix(TT[,1:4]) %*% Q + rep(1, nrow(X)) %*% t(c(- optimal_treshold / (t(Q) %*% t(C))[1], 0, 0, 0))
# colnames(SS) <- c("s1", "s2", "s3")
# pairs(TT[,1:4], col = colors[Ycl])
# pairs(SS[,1:4], col = colors[Ycl], asp = 1)

plot(NaN, xlim = range(SS[, 1]), ylim = range(SS[, 2]), xlab = "s1", ylab = "s2", main = "Latent space")
grid()
abline(h = 0, lty = 2)
abline(v = 0, lty = 2)
points(SS[, 1], SS[, 3], col = colors[Ycl])
points(mean(SS[index.CONTROL, 1]), mean(SS[index.CONTROL, 3]), col = colors["CTRL"], pch = 17, cex = 3)
points(mean(SS[index.SCHZ, 1]), mean(SS[index.SCHZ, 3]), col = colors["SCHZ"], pch = 17, cex = 3)





library(GGally)
# Using ggpairs:
pairs <- ggpairs(
  data = cbind(data.frame(SS), group = factor(Ycl)),
  columns = c("s1", "s2", "s3"),  # which columns to use
  mapping = aes(color = group, fill = group),  
  diag = list(continuous = wrap("densityDiag", alpha = 0.5)),  
  lower = list(continuous = wrap("points", alpha = 0.7, size = 2)),
  upper = list(continuous = wrap("points", alpha = 0.7, size = 2)),
  xlab = NULL, ylab = NULL, columnLabels  = NULL
) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  theme_bw() +
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank()
  ) 

ggsave(pairs, filename = paste(path_results, type, "/pairs.pdf", sep = ""),
       width = 8, height = 8)



## X loadings
R <- results_fPLS$X_loadings[,1:3]
R_tilde <- R %*% Q

for(h in 1:4) {
  FEMobj <- FEM(as.numeric(R_tilde[,h]/norm_l2(R_tilde[,h])), FEMbasis)
  write.vtu(FEMobj, paste(path_results, type, "/r", h, ".vtu", sep = ""))
}

for(h in 1:4) {
  FEMobj <- FEM(as.numeric(results_fPLS$Beta_hat[[h]]), FEMbasis)
  write.vtu(FEMobj, paste(path_results, type, "/beta_", h, ".vtu", sep = ""))
}

i = 2

prova <- FEM(as.numeric(Xc[i, ]), FEMbasis)
write.vtu(prova, paste(path_results, "data.vtu", sep = ""))

prova <- FEM(as.numeric(X_hat[[h]][i, ]), FEMbasis)
write.vtu(prova, paste(path_results, "fitted.vtu", sep = ""))

prova <- FEM(as.numeric(Xc[i, ] - X_hat[[h]][i, ]), FEMbasis)
write.vtu(prova, paste(path_results, "residual.vtu", sep = ""))
