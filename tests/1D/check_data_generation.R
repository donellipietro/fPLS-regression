# % %%%%%%%%%%%%%%%%%%% %
# % % Test: fPLS - 1D % %
# % %%%%%%%%%%%%%%%%%%% %

rm(list = ls())
graphics.off()

## global variables ----

test_suite <- "1D"
TEST_SUITE <- "fPLS - 1D"


## libraries ----

## algebraic utils
suppressMessages(library(pracma))

## statistical utils
suppressMessages(library(MASS))

## visualization
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(viridis))
suppressMessages(library(stringr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))

## sampling
suppressMessages(library(sf))
suppressMessages(library(sp))
suppressMessages(library(raster))

## fdaPDE
suppressMessages(library(fda))
suppressMessages(library(fdaPDE))
suppressMessages(library(fdaPDE2))

## sources ----
source("src/utils/directories.R")
source("src/utils/domain_and_locations.R")
source("src/utils/meshes.R")
source("src/utils/plots.R")
source("src/utils/wrappers.R")
source("src/utils/errors.R")
source("src/utils/pls.R")
source(paste("tests/", test_suite, "/utils/generate_data.R", sep = ""))


## directories ----
path_images <- paste("images/", test_suite, "/", sep = "") 
mkdir(path_images)
path_images <- paste(path_images, "check_data_generation", "/", sep = "") 
mkdir(path_images)


## options ----

name_mesh <- "unit_segment"
n_knots <- 101
n_locs <- 200
n_stat_units <- 100
locs_eq_nodes <- FALSE
n_comp <- 4
NSR_X_lc <- 3
NSR_Y <- 0.1


## domain & locations ----

## domain
domain <- fdaPDE2::Mesh(unit_segment(n_knots))
domain_boundary <- c(0,1)
spline_basis <- create.bspline.basis(rangeval = domain_boundary, breaks = domain$nodes)

## HR grid
n_knots_HR <- 601
grid_HR <- seq(0,1, length = n_knots_HR)
Psi_HR <- eval.basis(grid_HR, spline_basis)

## locations
set.seed(0)
locations <- seq(0,1, length = n_locs) + rnorm(n_locs, 0, 0.25/n_locs)
locations <- sort(unique(c(0,1, closest(grid_HR, locations))))
Psi_locs <- eval.basis(locations, spline_basis)

## data ----

set.seed(0)
generated_data <- generate_data(
  grid_HR, locations,
  n_stat_units = n_stat_units, 
  n_comp = n_comp,
  score_dist = "norm",
  NSR_X_lc = NSR_X_lc,
  NSR_Y = NSR_Y,
  seed = 1
)
data <- functional_regression_data(
  domain = domain,
  locations = locations,
  Y = generated_data$data$Y,
  X = generated_data$data$X
)


## data exploration

## Y, X

pdf(file = paste(path_images, "Y.pdf", sep = ""), width = 5, height = 5)
par(mar = c(0, 0, 0, 0))
plot(NaN, asp = 1, ylim = range(generated_data$data$Y), xlim = range(generated_data$expected_results$Y_true), 
     xlab = "", ylab = "", xaxt = "n", yaxt = "n") # , main = "Y vs Y_true", )
grid()
abline(0 ,1, col = "darkgrey", lty = 2, lwd = 2)
points(generated_data$expected_results$Y_true, generated_data$data$Y, pch = 20, cex = 1.2)
# abline(a = -2*generated_data$sigma_Y_noise, b = 1, col = "red", lty = 2)
# abline(a = -1*generated_data$sigma_Y_noise, b = 1, col = "red", lty = 2)
# abline(a = 1*generated_data$sigma_Y_noise, b = 1, col = "red", lty = 2)
# abline(a = 2*generated_data$sigma_Y_noise, b = 1, col = "red", lty = 2)
dev.off()


pdf(file = paste(path_images, "X_points.pdf", sep = ""), width = 11, height = 5)
par(mar = c(0, 0, 0, 0))
matplot(NaN, NaN, type = "l", lty = 1, ylim = range(generated_data$data$X[1:4, ]), xlim = range(grid_HR), 
        xlab = "", ylab = "", xaxt = "n", yaxt = "n", lwd = 2, col = brewer.pal(5, "RdPu")[2:5]) # , main = "X true", )
matpoints(locations, t(generated_data$data$X)[, 1:4], lwd = 2, pch = 20, lty = 2, col = brewer.pal(5, "RdPu")[2:5])
grid()
dev.off()

pdf(file = paste(path_images, "X_interpolant.pdf", sep = ""), width = 11, height = 5)
par(mar = c(0, 0, 0, 0))
matplot(NaN, NaN, type = "l", lty = 1, ylim = range(generated_data$data$X[1:4, ]), xlim = range(grid_HR), 
        xlab = "", ylab = "", xaxt = "n", yaxt = "n", lwd = 2, col = brewer.pal(5, "RdPu")[2:5]) # , main = "X true", )
matpoints(locations, t(generated_data$data$X)[, 1:4], lwd = 2, type = "b", pch = 20, lty = 2, col = brewer.pal(5, "RdPu")[2:5])
grid()
dev.off()


pdf(file = paste(path_images, "X.pdf", sep = ""), width = 11, height = 5)
par(mar = c(0, 0, 0, 0))
matplot(NaN, NaN, type = "l", lty = 1, ylim = range(generated_data$data$X[1:4, ]), xlim = range(grid_HR), 
        xlab = "", ylab = "", xaxt = "n", yaxt = "n", lwd = 2, col = brewer.pal(5, "RdPu")[2:5]) # , main = "X true", )
matpoints(grid_HR, t(generated_data$expected_results$X_true)[, 1:4], lwd = 2, type = "l", pch = 20, lty = 2, col = brewer.pal(5, "RdPu")[2:5])
grid()
dev.off()

f_hat <- NULL
for(i in 1:4){

  ## data init
  data <- spatial_data(
    domain = domain,
    observations = t(generated_data$data$X)[, i],
    locations = locations
  )
  
  ## model init
  model <- SRPDE(
    data = data,
    penalty = simple_laplacian_penalty(3, "SPLINE"),
    calibrator = gcv()
  )
  
  ## model fit
  lambda <- hyperparameters(space = 1)
  model$fit()
  
  f_hat <- cbind(f_hat, Psi_HR %*% model$results$f)

}


pdf(file = paste(path_images, "X_oversmoothing.pdf", sep = ""), width = 11, height = 5)
par(mar = c(0, 0, 0, 0))
matplot(grid_HR, f_hat, type = "l", lty = 1, ylim = range(generated_data$data$X[1:4, ]), 
        xlab = "", ylab = "", xaxt = "n", yaxt = "n", lwd = 2, col = brewer.pal(5, "RdPu")[2:5]) # , main = "X true", )
matlines(locations, t(generated_data$data$X)[, 1:4], lwd = 2, type = "b", pch = 20, lty = 2, col = brewer.pal(5, "RdPu")[2:5])
grid()
dev.off()


## Beta

title <- "Beta"
Beta_true <- generated_data$expected_results$Beta_true
Beta_locs <- generated_data$expected_results$Beta_true_locs
plot_list <- list()
for (h in 1:1) {
  plot_list[[h]] <- plot.curve_points(
    locations, Beta_locs[, h], size = 1
  ) + standard_plot_settings_fields()
}
for (h in 1:1) {
  plot_list[[1 + h]] <- plot.curve(
    grid_HR, Beta_true[, h]
  ) + standard_plot_settings_fields()
}
plot <- arrangeGrob(grobs = plot_list, nrow = 1)
plot <- labled_plots_grid(
  plot,
  title = title,
  width = 5,
  height = 5
)
grid.arrange(plot)

## X directions
title <- "X space directions"
labels_cols <- as.list(paste("w", 1:n_comp, sep = ""))
X_space_directions_true <- generated_data$expected_results$X_space_directions_true
X_space_directions_true_locs <- generated_data$expected_results$X_space_directions_true_locs
plot_list <- list()
for (h in 1:n_comp) {
  plot_list[[h]] <- plot.curve_points(
    locations, X_space_directions_true_locs[, h], size = 1
  ) + standard_plot_settings_fields()
}
for (h in 1:n_comp) {
  plot_list[[n_comp + h]] <- plot.curve(
    grid_HR, X_space_directions_true[, h]
  ) + standard_plot_settings_fields()
}
plot <- arrangeGrob(grobs = plot_list, nrow = 2)
plot <- labled_plots_grid(
  plot,
  title = title,
  labels_cols = labels_cols,
  width = 5,
  height = 5
)
grid.arrange(plot)

## X loadings
title <- "X loadings"
labels_cols <- as.list(paste("f", 1:n_comp, sep = ""))
X_loadings_true <- generated_data$expected_results$X_loadings_true
X_loadings_true_locs <- generated_data$expected_results$X_loadings_true_locs
plot_list <- list()
for (h in 1:n_comp) {
  plot_list[[h]] <- plot.curve_points(
    locations, X_loadings_true_locs[, h], size = 1
  ) +
    standard_plot_settings_fields()
}
for (h in 1:n_comp) {
  plot_list[[n_comp + h]] <- plot.curve(
    grid_HR, X_loadings_true[, h]
  ) +
    standard_plot_settings_fields()
}
plot <- arrangeGrob(grobs = plot_list, nrow = 2)
plot <- labled_plots_grid(
  plot,
  title = title,
  labels_cols = labels_cols,
  width = 5,
  height = 5
)
grid.arrange(plot)


## MV-PLS
model_MV_PLS <- MV_PLS_wrapped(data, n_comp = n_comp, center = FALSE, mode = "PLS-R")
Beta_MV_PLS <- model_MV_PLS$results$Beta_hat_locs[[n_comp]]
W_MV_PLS <- model_MV_PLS$results$X_space_directions_locs
R_MV_PLS <- model_MV_PLS$result$X_loadings_locs
T_MV_PLS <- model_MV_PLS$result$X_latent_scores


pdf(file = paste(path_images, "MV-PLS.pdf", sep = ""), height = 6, width = 5)

Beta_true_HR <- generated_data$expected_results$Beta_true
norm <-  norm_l2(Beta_true_HR)
Beta_true_HR <- Beta_true_HR / norm
Beta_MV_PLS <- Beta_MV_PLS / norm_l2(Beta_MV_PLS) * norm_l2(generated_data$expected_results$Beta_true_locs) / norm
plot(NA, xlim = c(0,1), ylim = range(c(Beta_MV_PLS, Beta_true_HR)), xlab = "", ylab = "", main = "Beta")
grid()
points(locations, Beta_MV_PLS)
points(grid_HR, Beta_true_HR, type = "l", col = "green")

X_space_directions_true_HR <- generated_data$expected_results$X_space_directions_true
for(h in 1:n_comp) {
  norm_X_latent_scores_true <- norm_l2(generated_data$expected_results$X_latent_scores_true[, h])
  norm_T_MV_PLS <- norm_l2(T_MV_PLS[, h])
  sign <- ifelse(norm_l2(generated_data$expected_results$X_latent_scores_true[, h] - T_MV_PLS[, h]) > norm_l2(generated_data$expected_results$X_latent_scores_true[, h] + T_MV_PLS[, h]), -1, 1)
  plot(NA, xlim = c(0,1), ylim = range(c(W_MV_PLS[, h]*norm_T_MV_PLS, X_space_directions_true_HR[, h]*norm_X_latent_scores_true)), xlab = "", ylab = "", main = paste("w", h, sep = ""))
  grid()
  points(grid_HR, X_space_directions_true_HR[, h]*norm_X_latent_scores_true, type = "l", col = "green")
  points(locations, sign*W_MV_PLS[, h]*norm_T_MV_PLS)
}

X_loadings_true_HR <- generated_data$expected_results$X_loadings_true
for(h in 1:n_comp) {
  norm_X_latent_scores_true <- norm_l2(generated_data$expected_results$X_latent_scores_true[, h])
  norm_T_MV_PLS <- norm_l2(T_MV_PLS[, h])
  sign <- ifelse(norm_l2(generated_data$expected_results$X_latent_scores_true[, h] - T_MV_PLS[, h]) > norm_l2(generated_data$expected_results$X_latent_scores_true[, h] + T_MV_PLS[, h]), -1, 1)
  plot(NA, xlim = c(0,1), ylim = range(c(R_MV_PLS[, h]*norm_T_MV_PLS, X_loadings_true_HR[, h]*norm_X_latent_scores_true)), xlab = "", ylab = "", main = paste("r", h, sep = ""))
  grid()
  points(grid_HR, X_loadings_true_HR[, h]*norm_X_latent_scores_true, type = "l", col = "green")
  points(locations, sign*R_MV_PLS[, h]*norm_T_MV_PLS)
}

dev.off()


## fPLS
lambda_grid <- fdaPDE2::hyperparameters(10^seq(-9, 6, by = 0.1))
model_fPLS <- fdaPDE2::fPLS(
  data = data,
  center = FALSE,
  solver = sequential(),
  mode = "fPLS-R",
  penalty = simple_laplacian_penalty(3, "SPLINE")
)
model_fPLS$fit(calibrator = gcv(lambda = lambda_grid, seed = 0), n_comp = n_comp)
Beta_fPLS_HR <- Psi_HR %*% model_fPLS$results$Beta_hat[[n_comp]]
W_fPLS_HR <- Psi_HR %*% model_fPLS$results$X_space_directions
R_fPLS_HR <- Psi_HR %*% model_fPLS$results$X_loadings
T_fPLS <- model_fPLS$results$X_latent_scores


pdf(file = paste(path_images, "fPLS.pdf", sep = ""), height = 6, width = 5)

Beta_true_HR <- generated_data$expected_results$Beta_true
Beta_true_HR <- Beta_true_HR / norm_l2(Beta_true_HR)
Beta_fPLS_HR <- Beta_fPLS_HR / norm_l2(Beta_fPLS_HR)
plot(NA, xlim = c(0,1), ylim = range(c(Beta_fPLS_HR, Beta_true_HR)), xlab = "", ylab = "", main = "Beta")
grid()
points(grid_HR, Beta_fPLS_HR, type = "l")
points(grid_HR, Beta_true_HR, type = "l", col = "green")

X_space_directions_true_HR <- generated_data$expected_results$X_space_directions_true
for(h in 1:n_comp) {
  norm_X_latent_scores_true <- norm_l2(generated_data$expected_results$X_latent_scores_true[, h])
  norm_T_fPLS <- norm_l2(T_fPLS[, h])
  sign <- ifelse(norm_l2(generated_data$expected_results$X_latent_scores_true[, h] - T_fPLS[, h]) > norm_l2(generated_data$expected_results$X_latent_scores_true[, h] + T_fPLS[, h]), -1, 1)
  plot(NA, xlim = c(0,1), ylim = range(c(W_fPLS_HR[, h]*norm_T_fPLS, X_space_directions_true_HR[, h]*norm_X_latent_scores_true)), xlab = "", ylab = "", main = paste("w", h, sep = ""))
  grid()
  points(grid_HR, X_space_directions_true_HR[, h]*norm_X_latent_scores_true, type = "l", col = "green")
  points(grid_HR, sign*W_fPLS_HR[, h]*norm_T_fPLS, type = "l")
}

X_loadings_true_HR <- generated_data$expected_results$X_loadings_true
for(h in 1:n_comp) {
  norm_X_latent_scores_true <- norm_l2(generated_data$expected_results$X_latent_scores_true[, h])
  norm_T_fPLS <- norm_l2(T_fPLS[, h])
  sign <- ifelse(norm_l2(generated_data$expected_results$X_latent_scores_true[, h] - T_fPLS[, h]) > norm_l2(generated_data$expected_results$X_latent_scores_true[, h] + T_fPLS[, h]), -1, 1)
  plot(NA, xlim = c(0,1), ylim = range(c(R_fPLS_HR[, h]*norm_T_fPLS, X_loadings_true_HR[, h]*norm_X_latent_scores_true)), xlab = "", ylab = "", main = paste("r", h, sep = ""))
  grid()
  points(grid_HR, X_loadings_true_HR[, h]*norm_X_latent_scores_true, type = "l", col = "green")
  points(grid_HR, sign*R_fPLS_HR[, h]*norm_T_fPLS, type = "l")
}

dev.off()


## B-fPLS
lambda_grid <- fdaPDE2::hyperparameters(10^seq(-9, 6, by = 0.1))
model_B_fPLS <- B_fPLS_wrapped(
  data = data,
  center = FALSE,
  mode = "fPLS-R",
  presmoothing = FALSE,
  penalized = FALSE
)
Psi_HR_B <- eval.basis(grid_HR, model_B_fPLS$basisobj)
Beta_B_fPLS_HR <- Psi_HR_B %*% model_B_fPLS$results$Beta_hat[[n_comp]]
W_B_fPLS_HR <- Psi_HR_B %*% model_B_fPLS$results$X_space_directions
R_B_fPLS_HR <- Psi_HR_B %*% model_B_fPLS$results$X_loadings
T_B_fPLS <- model_B_fPLS$results$X_latent_scores


pdf(file = paste(path_images, "B-fPLS.pdf", sep = ""), height = 6, width = 5)

Beta_true_HR <- generated_data$expected_results$Beta_true
Beta_true_HR <- Beta_true_HR / norm_l2(Beta_true_HR)
Beta_B_fPLS_HR <- Beta_B_fPLS_HR / norm_l2(Beta_B_fPLS_HR)
plot(NA, xlim = c(0,1), ylim = range(c(Beta_B_fPLS_HR, Beta_true_HR)), xlab = "", ylab = "", main = "Beta")
grid()
points(grid_HR, Beta_B_fPLS_HR, type = "l")
points(grid_HR, Beta_true_HR, type = "l", col = "green")

# for(h in 1:n_comp) {
#   X_loadings_true_HR <- eval.fd(grid_HR, Data2fd(as.vector(domain$nodes), generated_data$expected_results$X_space_directions_true, spline_basis))
#   norm_X_latent_scores_true <- norm_l2(generated_data$expected_results$X_latent_scores_true[, h])
#   norm_T_B_fPLS <- norm_l2(T_B_fPLS[, h])
#   sign <- ifelse(norm_l2(generated_data$expected_results$X_latent_scores_true[, h] - T_B_fPLS[, h]) > norm_l2(generated_data$expected_results$X_latent_scores_true[, h] + T_B_fPLS[, h]), -1, 1)
#   plot(NA, xlim = c(0,1), ylim = range(c(W_B_fPLS_HR[, h]*norm_T_B_fPLS, X_space_directions_true_HR[, h]*norm_X_latent_scores_true)), xlab = "", ylab = "", main = paste("r", h, sep = ""))
#   points(grid_HR, X_space_directions_true_HR[, h]*norm_X_latent_scores_true, type = "l", col = "green")
#   points(grid_HR, 100*sign*W_B_fPLS_HR[, h]*norm_T_B_fPLS, type = "l")
# }

X_loadings_true_HR <- generated_data$expected_results$X_loadings_true
for(h in 1:n_comp) {
  norm_X_latent_scores_true <- norm_l2(generated_data$expected_results$X_latent_scores_true[, h])
  norm_T_B_fPLS <- norm_l2(T_B_fPLS[, h])
  sign <- ifelse(norm_l2(generated_data$expected_results$X_latent_scores_true[, h] - T_B_fPLS[, h]) > norm_l2(generated_data$expected_results$X_latent_scores_true[, h] + T_B_fPLS[, h]), -1, 1)
  plot(NA, xlim = c(0,1), ylim = range(c(R_B_fPLS_HR[, h]*norm_T_B_fPLS, X_loadings_true_HR[, h]*norm_X_latent_scores_true)), xlab = "", ylab = "", main = paste("r", h, sep = ""))
  grid()
  points(grid_HR, X_loadings_true_HR[, h]*norm_X_latent_scores_true, type = "l", col = "green")
  points(grid_HR, sign*R_B_fPLS_HR[, h]*norm_T_B_fPLS, type = "l")
}

dev.off()


## B-fPLS-smooth
lambda_grid <- fdaPDE2::hyperparameters(10^seq(-9, 6, by = 0.1))
model_B_fPLS_smooth <- B_fPLS_wrapped(
  data = data,
  center = FALSE,
  mode = "fPLS-R",
  presmoothing = TRUE,
  penalized = FALSE,
  lambda_grid = lambda_grid
)
Psi_HR_B <- eval.basis(grid_HR, model_B_fPLS$basisobj)
Beta_B_fPLS_smooth_HR <- Psi_HR_B %*% model_B_fPLS_smooth$results$Beta_hat[[n_comp]]
W_B_fPLS_smooth_HR <- Psi_HR_B %*% model_B_fPLS_smooth$results$X_space_directions
R_B_fPLS_smooth_HR <- Psi_HR_B %*% model_B_fPLS_smooth$results$X_loadings
T_B_fPLS_smooth <- model_B_fPLS_smooth$results$X_latent_scores


pdf(file = paste(path_images, "B-fPLS-smooth.pdf", sep = ""), height = 6, width = 5)

Beta_true_HR <- generated_data$expected_results$Beta_true
Beta_true_HR <- Beta_true_HR / norm_l2(Beta_true_HR)
Beta_B_fPLS_smooth_HR <- Beta_B_fPLS_smooth_HR / norm_l2(Beta_B_fPLS_smooth_HR)
plot(NA, xlim = c(0,1), ylim = range(c(Beta_B_fPLS_smooth_HR, Beta_true_HR)), xlab = "", ylab = "", main = "Beta")
grid()
points(grid_HR, Beta_B_fPLS_smooth_HR, type = "l")
points(grid_HR, Beta_true_HR, type = "l", col = "green")

# for(h in 1:n_comp) {
#   X_loadings_true_HR <- eval.fd(grid_HR, Data2fd(as.vector(domain$nodes), generated_data$expected_results$X_space_directions_true, spline_basis))
#   norm_X_latent_scores_true <- norm_l2(generated_data$expected_results$X_latent_scores_true[, h])
#   norm_T_B_fPLS <- norm_l2(T_B_fPLS[, h])
#   sign <- ifelse(norm_l2(generated_data$expected_results$X_latent_scores_true[, h] - T_B_fPLS[, h]) > norm_l2(generated_data$expected_results$X_latent_scores_true[, h] + T_B_fPLS[, h]), -1, 1)
#   plot(NA, xlim = c(0,1), ylim = range(c(W_B_fPLS_HR[, h]*norm_T_B_fPLS, X_space_directions_true_HR[, h]*norm_X_latent_scores_true)), xlab = "", ylab = "", main = paste("r", h, sep = ""))
#   points(grid_HR, X_space_directions_true_HR[, h]*norm_X_latent_scores_true, type = "l", col = "green")
#   points(grid_HR, 100*sign*W_B_fPLS_HR[, h]*norm_T_B_fPLS, type = "l")
# }

X_loadings_true_HR <- generated_data$expected_results$X_loadings_true
for(h in 1:n_comp) {
  norm_X_latent_scores_true <- norm_l2(generated_data$expected_results$X_latent_scores_true[, h])
  norm_T_B_fPLS_smooth <- norm_l2(T_B_fPLS_smooth[, h])
  sign <- ifelse(norm_l2(generated_data$expected_results$X_latent_scores_true[, h] - T_B_fPLS_smooth[, h]) > norm_l2(generated_data$expected_results$X_latent_scores_true[, h] + T_B_fPLS_smooth[, h]), -1, 1)
  plot(NA, xlim = c(0,1), ylim = range(c(R_B_fPLS_smooth_HR[, h]*norm_T_B_fPLS_smooth, X_loadings_true_HR[, h]*norm_X_latent_scores_true)), xlab = "", ylab = "", main = paste("r", h, sep = ""))
  grid()
  points(grid_HR, X_loadings_true_HR[, h]*norm_X_latent_scores_true, type = "l", col = "green")
  points(grid_HR, sign*R_B_fPLS_smooth_HR[, h]*norm_T_B_fPLS_smooth, type = "l")
}

dev.off()


## PB-fPLS
lambda_grid <- fdaPDE2::hyperparameters(10^seq(-9, 6, by = 0.1))
model_PB_fPLS <- B_fPLS_wrapped(
  data = data,
  center = FALSE,
  mode = "fPLS-R",
  presmoothing = FALSE,
  penalized = TRUE,
  lambda_grid = lambda_grid
)
Psi_HR_B <- eval.basis(grid_HR, model_PB_fPLS$basisobj)
Beta_PB_fPLS_HR <- Psi_HR_B %*% model_PB_fPLS$results$Beta_hat[[n_comp]]
W_PB_fPLS_HR <- Psi_HR_B %*% model_PB_fPLS$results$X_space_directions
R_PB_fPLS_HR <- Psi_HR_B %*% model_PB_fPLS$results$X_loadings
T_PB_fPLS <- model_PB_fPLS$results$X_latent_scores


pdf(file = paste(path_images, "PB-fPLS.pdf", sep = ""), height = 6, width = 5)

Beta_true_HR <- generated_data$expected_results$Beta_true
Beta_true_HR <- Beta_true_HR / norm_l2(Beta_true_HR)
Beta_PB_fPLS_HR <- Beta_PB_fPLS_HR / norm_l2(Beta_PB_fPLS_HR)
plot(NA, xlim = c(0,1), ylim = range(c(Beta_PB_fPLS_HR, Beta_true_HR)), xlab = "", ylab = "", main = "Beta")
grid()
points(grid_HR, Beta_PB_fPLS_HR, type = "l")
points(grid_HR, Beta_true_HR, type = "l", col = "green")

# for(h in 1:n_comp) {
#   X_loadings_true_HR <- eval.fd(grid_HR, Data2fd(as.vector(domain$nodes), generated_data$expected_results$X_space_directions_true, spline_basis))
#   norm_X_latent_scores_true <- norm_l2(generated_data$expected_results$X_latent_scores_true[, h])
#   norm_T_PB_fPLS <- norm_l2(T_PB_fPLS[, h])
#   sign <- ifelse(norm_l2(generated_data$expected_results$X_latent_scores_true[, h] - T_PB_fPLS[, h]) > norm_l2(generated_data$expected_results$X_latent_scores_true[, h] + T_PB_fPLS[, h]), -1, 1)
#   plot(NA, xlim = c(0,1), ylim = range(c(W_PB_fPLS_HR[, h]*norm_T_PB_fPLS, X_space_directions_true_HR[, h]*norm_X_latent_scores_true)), xlab = "", ylab = "", main = paste("r", h, sep = ""))
#   points(grid_HR, X_space_directions_true_HR[, h]*norm_X_latent_scores_true, type = "l", col = "green")
#   points(grid_HR, 100*sign*W_PB_fPLS_HR[, h]*norm_T_PB_fPLS, type = "l")
# }

X_loadings_true_HR <- generated_data$expected_results$X_loadings_true
for(h in 1:n_comp) {
  norm_X_latent_scores_true <- norm_l2(generated_data$expected_results$X_latent_scores_true[, h])
  norm_T_PB_fPLS <- norm_l2(T_PB_fPLS[, h])
  sign <- ifelse(norm_l2(generated_data$expected_results$X_latent_scores_true[, h] - T_PB_fPLS[, h]) > norm_l2(generated_data$expected_results$X_latent_scores_true[, h] + T_PB_fPLS[, h]), -1, 1)
  plot(NA, xlim = c(0,1), ylim = range(c(R_PB_fPLS_HR[, h]*norm_T_PB_fPLS, X_loadings_true_HR[, h]*norm_X_latent_scores_true)), xlab = "", ylab = "", main = paste("r", h, sep = ""))
  grid()
  points(grid_HR, X_loadings_true_HR[, h]*norm_X_latent_scores_true, type = "l", col = "green")
  points(grid_HR, sign*R_PB_fPLS_HR[, h]*norm_T_PB_fPLS, type = "l")
}

dev.off()

pdf(file = paste(path_images, "Beta_", n_knots, ".pdf", sep = ""), height = 5, width = 8)

par(mai = c(0, 0, 0, 0))
# range(c(Beta_MV_PLS, Beta_PB_fPLS_HR))
colors <- c(brewer.pal(3, "Greys")[2], brewer.pal(3, "Blues")[3], brewer.pal(9, "RdPu")[c(5, 7, 9)])
plot(NaN, xlim = c(0,1), ylim = c(-0.08557198, 0.11170637), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
grid()
points(locations, Beta_MV_PLS, col = colors[1], type = "b", lty = 2, pch = 20)
points(grid_HR, Beta_B_fPLS_HR, col = colors[3], type = "l", pch = 20, lwd = 2)
points(grid_HR, Beta_B_fPLS_smooth_HR, col = colors[4], type = "l", pch = 20, lwd = 2)
points(grid_HR, Beta_PB_fPLS_HR, col = colors[5], type = "l", pch = 20, lwd = 2)
points(grid_HR, Beta_fPLS_HR, col = colors[2], type = "l", pch = 20, lwd = 2)
points(grid_HR, Beta_true_HR, col = "black", type = "l", lty = 2, pch = 20, lwd = 1)

dev.off()


