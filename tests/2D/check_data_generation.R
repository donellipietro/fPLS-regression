# % %%%%%%%%%%%%%%%%%%% %
# % % Test: fPLS - 2D % %
# % %%%%%%%%%%%%%%%%%%% %

rm(list = ls())
graphics.off()

## global variables ----

test_suite <- "2D"
TEST_SUITE <- "fPLS - 2D"


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

plot_qualitative_check <- function(GT, Est1, Est2, title, label, at_locs) {
  K <- ncol(GT)
  labels_cols <- as.list(paste(label, " ", 1:K, sep = ""))
  plot_list <- list()
  
  colormap <- c(
    "#053061", # (R=0.0196078, G=0.188235,  B=0.380392)
    "#2166AC", # (R=0.129412,  G=0.4,       B=0.67451)
    "#4393C3", # (R=0.262745,  G=0.576471,  B=0.764706)
    "#92C5DE", # (R=0.572549,  G=0.772549,  B=0.870588)
    "#D1E5F0", # (R=0.819608,  G=0.898039,  B=0.941176)
    "#F7F7F7", # (R=0.968627,  G=0.968627,  B=0.968627)
    "#FDDBC7", # (R=0.992157,  G=0.858824,  B=0.780392)
    "#F4A582", # (R=0.956863,  G=0.647059,  B=0.509804)
    "#D6604D", # (R=0.839216,  G=0.376471,  B=0.301961)
    "#B2182B", # (R=0.698039,  G=0.0941176, B=0.168627)
    "#67001F"  # (R=0.403922,  G=0,         B=0.121569)
  )
  
  
  for (k in 1:K) {
    limits <- range(GT[, k])
    plot_list[[k]] <- plot.field_tile_gradient(
      locations, GT[, k],
      limits = limits,
      LEGEND = FALSE,
      ISOLINES = TRUE,
      colormap = colormap,
      boundary = domain_boundary,
    ) + standard_plot_settings_fields()
  }
  for (k in 1:K) {
    limits <- range(GT[, k])
    plot_list[[K + k]] <- plot.field_tile_gradient(
      locations, Est1[, k],
      limits = limits,
      LEGEND = FALSE,
      ISOLINES = TRUE,
      colormap = colormap,
      boundary = domain_boundary
    ) + standard_plot_settings_fields()
  }
  for (k in 1:K) {
    limits <- range(GT[, k])
    plot_list[[2*K + k]] <- plot.field_tile_gradient(
      locations, Est2[, k],
      limits = limits,
      LEGEND = FALSE,
      ISOLINES = TRUE,
      colormap = colormap,
      boundary = domain_boundary,
    ) + standard_plot_settings_fields()
  }
  plot <- arrangeGrob(grobs = plot_list, nrow = 2, as.table = FALSE)
  plot <- labled_plots_grid(
    plot,
    title = title,
    labels_cols = c("True", "MV-PLS", "R1-fPLS"),
    labels_rows = labels_cols,
    width = 5,
    height = 5
  )
  grid.arrange(plot)
}

## directories ----

path_images <- paste("images/", test_suite, "/", sep = "") 
mkdir(path_images)
path_images <- paste(path_images, "check_data_generation", "/", sep = "") 
mkdir(path_images)


## options ----

name_mesh <- "unit_square"
n_nodes <- 30*30
n_locs <- 1600
n_stat_units <- 100
locs_eq_nodes <- FALSE
n_comp <- 4
NSR_X_lc <- 3
NSR_Y <- 0.1


## domain & locations ----

## domain
generated_domain <- generate_domain(name_mesh, n_nodes)
domain <- generated_domain$domain
domain_boundary <- generated_domain$domain_boundary
loadings_true_generator <- generated_domain$loadings_true_generator
mesh <- generated_domain$mesh

## HR mesh
n_nodes_HR <- 79*79
generated_domain <- generate_domain(name_mesh, n_nodes_HR)
mesh_HR <- generated_domain$mesh
nodes_HR <- mesh_HR$nodes

## locations
generated_domain <- generate_domain(name_mesh, 40*40)
locations <- generated_domain$mesh$nodes

plot(nodes_HR, asp = 1, pch = 20)
points(locations, col = "red", pch = 20)



## data ----

generated_data <- generate_data(
  mesh_HR, locations,
  n_stat_units = n_stat_units, 
  n_comp = n_comp,
  score_dist = "norm",
  NSR_X_lc = NSR_X_lc,
  NSR_Y = NSR_Y,
  seed = 7
)
data <- functional_regression_data(
  domain = domain,
  locations = locations,
  Y = generated_data$data$Y,
  X = generated_data$data$X
)
n_resps <- ncol(generated_data$data$Y)
T_true <- generated_data$expected_results$X_latent_scores_true

## data exploration
par(mfrow = c(1,2))

plot(NaN, asp = 1, xlab = "Y_true", ylab = "Y", main = "Y1 vs Y1_true", ylim = range(generated_data$data$Y[,1]), xlim = range(generated_data$expected_results$Y_true[,1]))
grid()
points(generated_data$expected_results$Y_true[,1], generated_data$data$Y[,1])
abline(0 ,1, col = "red")
abline(a = -2*generated_data$sigma_Y_noises[1], b = 1, col = "red", lty = 2)
abline(a = -1*generated_data$sigma_Y_noises[1], b = 1, col = "red", lty = 2)
abline(a = 1*generated_data$sigma_Y_noises[1], b = 1, col = "red", lty = 2)
abline(a = 2*generated_data$sigma_Y_noises[1], b = 1, col = "red", lty = 2)

plot(NaN, asp = 1, xlab = "Y_true", ylab = "Y", main = "Y2 vs Y2_true", ylim = range(generated_data$data$Y[,2]), xlim = range(generated_data$expected_results$Y_true[,2]))
grid()
points(generated_data$expected_results$Y_true[,2], generated_data$data$Y[,2])
abline(0 ,1, col = "red")
abline(a = -2*generated_data$sigma_Y_noises[2], b = 1, col = "red", lty = 2)
abline(a = -1*generated_data$sigma_Y_noises[2], b = 1, col = "red", lty = 2)
abline(a = 1*generated_data$sigma_Y_noises[2], b = 1, col = "red", lty = 2)
abline(a = 2*generated_data$sigma_Y_noises[2], b = 1, col = "red", lty = 2)

## view generators
title <- "X generators"
labels_cols <- as.list(paste("f", 1:n_comp, sep = ""))
plot_list <- list()
for (h in 1:n_comp) {
  plot_list[[h]] <- plot.field_tile(
    nodes_HR, cube_eigenfunction(nodes_HR, h),
    LEGEND = FALSE,
    boundary = domain_boundary,
  ) +
    standard_plot_settings_fields()
}
for (h in 1:n_comp) {
  plot_list[[n_comp + h]] <- plot.field_points(
    locations, cube_eigenfunction(locations, h),
    LEGEND = FALSE,
    boundary = domain_boundary, size = 1.5
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


## view X space directions
title <- "X space directions"
labels_cols <- as.list(paste("w", 1:n_comp, sep = ""))
X_space_directions_true <- generated_data$expected_results$X_space_directions_true
X_space_directions_true_locs <- generated_data$expected_results$X_space_directions_true_locs
plot_list <- list()
for (h in 1:n_comp) {
  plot_list[[h]] <- plot.field_tile(
    nodes_HR, X_space_directions_true[, h],
    LEGEND = FALSE,
    boundary = domain_boundary
  ) + standard_plot_settings_fields()
}
for (h in 1:n_comp) {
  plot_list[[n_comp + h]] <- plot.field_points(
    locations, X_space_directions_true_locs[, h],
    LEGEND = FALSE,
    boundary = domain_boundary, size = 1.5
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

## view loadings
title <- "X loadings"
labels_cols <- as.list(paste("r", 1:n_comp, sep = ""))
X_loadings_true <- generated_data$expected_results$X_loadings_true
X_loadings_true_locs <- generated_data$expected_results$X_loadings_true_locs
plot_list <- list()
for (h in 1:n_comp) {
  plot_list[[h]] <- plot.field_tile(
    nodes_HR, X_loadings_true[, h],
    LEGEND = FALSE,
    boundary = domain_boundary,
  ) +
    standard_plot_settings_fields()
}
for (h in 1:n_comp) {
  plot_list[[n_comp + h]] <- plot.field_points(
    locations, X_loadings_true_locs[, h],
    LEGEND = FALSE,
    boundary = domain_boundary, size = 1.5
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


## view beta
title <- "Beta"
labels_cols <- as.list(paste("Beta", 1:n_resps, sep = ""))
Beta_true <- generated_data$expected_results$Beta_true
Beta_true_locs <- generated_data$expected_results$Beta_true_locs
plot_list <- list()
for (l in 1:n_resps) {
  limits <- range(Beta_true[, l])
  plot_list[[l]] <- plot.field_tile(
    nodes_HR, Beta_true[, l],
    limits = limits,
    LEGEND = FALSE,
    boundary = domain_boundary,
  ) +
    standard_plot_settings_fields()
}
for (l in 1:n_resps) {
  limits <- range(Beta_true[, l])
  plot_list[[n_resps + l]] <- plot.field_points(
    locations, Beta_true_locs[, l],
    limits = limits,
    LEGEND = FALSE,
    boundary = domain_boundary, size = 1.5
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

unique_seed=14
## data samples
title <- "X samples"
labels_cols <- list("True", "True at locations", "Observed")
set.seed(unique_seed)
n_samples <- 2
sampled_units<-sample(1:n_stat_units, n_samples,replace=FALSE)
sampled_collapsed=paste(as.character(sampled_units),collapse= " ")
labels_rows <- as.list(paste("x ",sampled_units[1:3], sep = ""))
plot_list <- list()
for (i in 1:length(sampled_units)) {
  curr_index=sampled_units[i]
  limits <- range(generated_data$data$X[curr_index, ])
  plot_list[[(i - 1) * 3 + 1]] <- plot.field_tile(nodes_HR, generated_data$expected_results$X_true[curr_index, ], boundary = domain_boundary, limits = limits, ISOLINES = TRUE) + standard_plot_settings_fields()
  plot_list[[(i - 1) * 3 + 2]] <- plot.field_points(locations, generated_data$expected_results$X_true_locs[curr_index, ], boundary = domain_boundary, limits = limits, size = 1.5) + standard_plot_settings_fields()
  plot_list[[(i - 1) * 3 + 3]] <- plot.field_points(locations, generated_data$data$X[curr_index, ], boundary = domain_boundary, limits = limits, size = 1.5) + standard_plot_settings_fields()
}
plot <- arrangeGrob(grobs = plot_list, nrow = n_samples)
plot <- labled_plots_grid(plot, title, labels_cols)
grid.arrange(plot)


#### Comparative plots

## MV-PLS
model_MV_PLS <- MV_PLS_wrapped(data, n_comp = n_comp, center = FALSE, mode = "PLS-R")
Beta_MV_PLS <- model_MV_PLS$results$Beta_hat_locs[[n_comp]]
W_MV_PLS <- model_MV_PLS$results$X_space_directions_locs
R_MV_PLS <- model_MV_PLS$result$X_loadings_locs
T_MV_PLS <- model_MV_PLS$result$X_latent_scores

## fPLS
lambda_grid <- fdaPDE2::hyperparameters(10^seq(-9, 2, by = 0.2))
model_fPLS <- fdaPDE2::fPLS(
  data = data,
  center = FALSE,
  solver = sequential(),
  mode = "fPLS-R"
)
model_fPLS$fit(calibrator = gcv(lambda = lambda_grid, seed = 0), n_comp = n_comp)
Beta_fPLS_locs <- evaluate_field(locations, model_fPLS$results$Beta_hat[[n_comp]] , mesh)
Beta_fPLS <- evaluate_field(nodes_HR, model_fPLS$results$Beta_hat[[n_comp]], mesh)   
W_fPLS <- evaluate_field(nodes_HR, model_fPLS$results$X_space_directions, mesh)   
R_fPLS <- evaluate_field(nodes_HR, model_fPLS$results$X_loadings, mesh)   
T_fPLS <- model_fPLS$results$X_latent_scores



Beta_true_locs <- generated_data$expected_results$Beta_true_locs
for(l in 1:n_resps) {
  Beta_true_locs[,l] <- Beta_true_locs[,l] / norm_l2(Beta_true_locs[,l])
  Beta_MV_PLS[,l] <- Beta_MV_PLS[,l] / norm_l2(Beta_MV_PLS[,l])
  Beta_fPLS_locs[,l] <- Beta_fPLS_locs[,l] / norm_l2(Beta_fPLS_locs[,l])
}

pdf(file = paste(path_images, "MV-PLS-beta.pdf"), height = 6, width = 8)
plot_qualitative_check(Beta_true_locs, Beta_MV_PLS, Beta_fPLS_locs, title = "Beta", label = "Beta", at_locs = TRUE) 
dev.off()


######


R_true <- generated_data$expected_results$X_loadings_true
for(h in 1:n_comp) {
  R_true[,h] <- R_true[,h] / norm_l2(generated_data$expected_results$X_loadings_true_locs[,h])
  T_true[,h] <- T_true[,h] / norm_l2(T_true[,h])
  R_MV_PLS[,h] <- R_MV_PLS[,h] / norm_l2(R_MV_PLS[,h])
  T_MV_PLS[,h] <- T_MV_PLS[,h] / norm_l2(T_MV_PLS[,h])
  sign <- ifelse(norm_l2(T_true[,h] - T_MV_PLS[,h]) > norm_l2(T_true[,h] + T_MV_PLS[,h]), -1, 1)
  R_MV_PLS[,h] <- sign * R_MV_PLS[,h]
}
# pdf(file = paste(path_images, "MV-PLS-directions.pdf"), height = 5, width = 10)
plot_qualitative_check(R_true, R_MV_PLS, title = "X space directions", label = "r", at_locs = TRUE) 
# dev.off()

W_true <- generated_data$expected_results$X_space_directions_true
for(h in 1:n_comp) {
  W_true[,h] <- W_true[,h] / norm_l2(generated_data$expected_results$X_space_directions_true_locs[,h])
  T_true[,h] <- T_true[,h] / norm_l2(T_true[,h])
  W_MV_PLS[,h] <- W_MV_PLS[,h] / norm_l2(W_MV_PLS[,h])
  T_MV_PLS[,h] <- T_MV_PLS[,h] / norm_l2(T_MV_PLS[,h])
  sign <- ifelse(norm_l2(T_true[,h] - T_MV_PLS[,h]) > norm_l2(T_true[,h] + T_MV_PLS[,h]), -1, 1)
  W_MV_PLS[,h] <- sign * W_MV_PLS[,h]
}
# pdf(file = paste(path_images, "MV-PLS-loadings.pdf"), height = 5, width = 10)
plot_qualitative_check(W_true, W_MV_PLS, title = "X loadings", label = "w", at_locs = TRUE) 
# dev.off()







# pdf(file = paste(path_images, "fPLS-beta.pdf"), height = 5, width = 5)
plot_qualitative_check(Beta_true, Beta_fPLS, title = "Beta", label = "Beta", at_locs = FALSE)
# dev.off()

R_true <- generated_data$expected_results$X_loadings_true
for(h in 1:n_comp) {
  R_true[,h] <- R_true[,h] / norm_l2(R_true[,h])
  T_true[,h] <- T_true[,h] / norm_l2(T_true[,h])
  R_fPLS[,h] <- R_fPLS[,h] / norm_l2(R_fPLS[,h])
  T_fPLS[,h] <- T_fPLS[,h] / norm_l2(T_fPLS[,h])
  sign <- ifelse(norm_l2(T_true[,h] - T_fPLS[,h]) > norm_l2(T_true[,h] + T_fPLS[,h]), -1, 1)
  R_fPLS[,h] <- sign * R_fPLS[,h]
}
# pdf(file = paste(path_images, "fPLS-directions.pdf"), height = 5, width = 10)
plot_qualitative_check(R_true, R_fPLS, title = "X space directions", label = "r", at_locs = FALSE) 
# dev.off()

W_true <- generated_data$expected_results$X_space_directions_true
for(h in 1:n_comp) {
  W_true[,h] <- W_true[,h] / norm_l2(W_true[,h])
  T_true[,h] <- T_true[,h] / norm_l2(T_true[,h])
  W_fPLS[,h] <- W_fPLS[,h] / norm_l2(W_fPLS[,h])
  T_fPLS[,h] <- T_fPLS[,h] / norm_l2(T_fPLS[,h])
  sign <- ifelse(norm_l2(T_true[,h] - T_fPLS[,h]) > norm_l2(T_true[,h] + T_fPLS[,h]), -1, 1)
  W_fPLS[,h] <- sign * W_fPLS[,h]
}
# pdf(file = paste(path_images, "fPLS-loadings.pdf"), height = 5, width = 10)
plot_qualitative_check(W_true, W_fPLS, title = "X loadings", label = "w", at_locs = FALSE) 
# dev.off()




