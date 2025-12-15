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
path_images <- paste(path_images, "check_population_convergence", "/", sep = "") 
mkdir(path_images)


## options ----

name_mesh <- "unit_segment"
n_knots <- 11
n_stat_units <- 200
locs_eq_nodes <- FALSE
n_comp <- 4
NSR_X_lc <- 1
NSR_Y <- 0.01


## domain & locations ----

## domain
domain <- fdaPDE2::Mesh(unit_segment(n_knots))
domain_boundary <- c(0,1)

## HR grid
n_knots_HR <- 200
grid_HR <- seq(0,1, length = n_knots_HR)

## locations
locations <- grid_HR


## data ----

Beta_list_GT <- list()
W_list_GT <- list()
R_list_GT <- list()

Beta_list_MV_PLS <- list()
W_list_MV_PLS <- list()
R_list_MV_PLS <- list()

Beta_list_fPLS <- list()
W_list_fPLS <- list()
R_list_fPLS <- list()


n_stat_units_vect <- round(10^(seq(2, 6, by = 0.25)))
for(i in 1:length(n_stat_units_vect)) {
  
  ## generate data
  generated_data <- generate_data(
    grid_HR, locations,
    n_stat_units = n_stat_units_vect[i], 
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
  
  ## GT
  Beta_list_GT[[i]] <- generated_data$expected_results$Beta_true
  W_list_GT[[i]] <- generated_data$expected_results$X_space_directions_true
  R_list_GT[[i]] <- generated_data$expected_results$X_loadings_true
  
  ## MV-PLS
  model_MV_PLS <- MV_PLS_wrapped(data, n_comp = n_comp, center = FALSE, mode = "PLS-R")
  Beta_list_MV_PLS[[i]] <- model_MV_PLS$results$Beta_hat_locs[[n_comp]]
  W_list_MV_PLS[[i]] <- model_MV_PLS$results$X_space_directions_locs
  R_list_MV_PLS[[i]] <- model_MV_PLS$result$X_loadings_locs
  
  ## fPLS
  lambda_grid <- fdaPDE2::hyperparameters(10^seq(-9, 6, by = 0.2))
  model_fPLS <- fdaPDE2::fPLS(
    data = data,
    center = FALSE,
    solver = sequential(),
    mode = "fPLS-R",
    penalty = simple_laplacian_penalty(3, "SPLINE")
  )
  model_fPLS$fit(calibrator = gcv(lambda = lambda_grid, seed = 0), n_comp = n_comp)
  Beta_list_fPLS[[i]] <- model_fPLS$evaluate(model_fPLS$results$Beta_hat[[n_comp]])
  W_list_fPLS[[i]] <- model_fPLS$evaluate(model_fPLS$results$X_space_directions)
  R_list_fPLS[[i]] <- model_fPLS$evaluate(model_fPLS$results$X_loadings)
}


plot_population_convergence <- function(Beta_list, R_list, W_list, population = TRUE) {
  
  n <- length(n_stat_units_vect)
  
  grayscale_colors <- gray.colors(length(n_stat_units_vect), start = 0.8, end = 0)
  
  for(i in 1:n) {
    Beta_list[[i]] <- Beta_list[[i]]/norm_l2(Beta_list[[i]])
    Beta_list_GT[[i]] <- Beta_list_GT[[i]]/norm_l2(Beta_list_GT[[i]])
    for(h in 1:n_comp) {
      R_list[[i]][,h] <- R_list[[i]][,h]/norm_l2(R_list[[i]][,h])
      R_list_GT[[i]][,h] <- R_list_GT[[i]][,h]/norm_l2(R_list_GT[[i]][,h])
      W_list[[i]][,h] <- W_list[[i]][,h]/norm_l2(W_list[[i]][,h])
      W_list_GT[[i]][,h] <- W_list_GT[[i]][,h]/norm_l2(W_list_GT[[i]][,h])
    }
  }
  
  
  
  par(mfrow = c(1, 2))
  error_PGT <- c()
  error_GT <- c()
  plot(NaN, xlim = c(0,1), ylim = range(Beta_list), xlab = "", ylab = "", main = "Beta")
  grid()
  for(i in 1:n) {
    points(as.vector(locations), Beta_list[[i]], type = "l", col = grayscale_colors[i])
    error_GT[i] <- norm_l2(Beta_list[[i]] - Beta_list_GT[[i]])
    error_PGT[i] <- norm_l2(Beta_list[[i]] - Beta_list_GT[[n]])
  }
  points(as.vector(grid_HR), Beta_list_GT[[n]], type = "l", col = "red")
  ## errors
  if(population){ ylim <- range(c(error_PGT, error_GT[-n]))
  } else { ylim <- range(error_PGT[-n]) }
  plot(NaN, log = "xy", ylab = "", xlab = "# stat. units", main = "Error",
       xlim = range(n_stat_units_vect), ylim = ylim)
  grid()
  points(n_stat_units_vect[-n], error_PGT[-n], col = "red", pch = 2)
  if(population){ points(n_stat_units_vect, error_GT)}
  
  
  for(h in 1:n_comp){
    par(mfrow = c(1, 2))
    error_PGT <- c()
    error_GT <- c()
    plot(NaN, xlim = c(0,1), ylim = range(W_list), xlab = "", ylab = "", main = paste("w", h, sep = ""))
    grid()
    for(i in 1:n) {
      if(norm_l2(W_list[[i]][,h] + W_list_GT[[n]][,h]) < norm_l2(W_list[[i]][,h] - W_list_GT[[n]][,h]))
        W_list[[i]][, h] = -W_list[[i]][, h]
      points(as.vector(locations), W_list[[i]][, h], type = "l", col = grayscale_colors[i])
      error_GT[i] <- norm_l2(W_list[[i]][,h] - W_list_GT[[i]][,h])
      error_PGT[i] <- norm_l2(W_list[[i]][,h] - W_list_GT[[n]][,h])
    }
    points(as.vector(grid_HR), W_list_GT[[n]][,h], type = "l", col = "red")
    ## errors
    if(population){ ylim <- range(c(error_PGT, error_GT[-n]))
    } else { ylim <- range(error_PGT[-n]) }
    plot(NaN, log = "xy", ylab = "", xlab = "# stat. units", main = "Error",
         xlim = range(n_stat_units_vect), ylim = ylim)
    grid()
    points(n_stat_units_vect, error_PGT, col = "red", pch = 2)
    if(population){ points(n_stat_units_vect, error_GT)}
  }
  
  for(h in 1:n_comp){
    par(mfrow = c(1, 2))
    error_PGT <- c()
    error_GT <- c()
    plot(NaN, xlim = c(0,1), ylim = range(R_list), xlab = "", ylab = "", main = paste("r", h, sep = ""))
    grid()
    for(i in 1:n) {
      if(norm_l2(R_list[[i]][,h] + R_list_GT[[n]][,h]) < norm_l2(R_list[[i]][,h] - R_list_GT[[n]][,h]))
        R_list[[i]][, h] = -R_list[[i]][, h]
      points(as.vector(locations), R_list[[i]][, h], type = "l", col = grayscale_colors[i])
      error_GT[i] <- norm_l2(R_list[[i]][,h] - R_list_GT[[i]][,h])
      error_PGT[i] <- norm_l2(R_list[[i]][,h] - R_list_GT[[n]][,h])
    }
    points(as.vector(grid_HR), R_list_GT[[n]][,h], type = "l", col = "red")
    ## errors
    if(population){ ylim <- range(c(error_PGT, error_GT[-n]))
    } else { ylim <- range(error_PGT[-n]) }
    plot(NaN, log = "xy", ylab = "", xlab = "# stat. units", main = "Error",
         xlim = range(n_stat_units_vect), ylim = ylim)
    grid()
    points(n_stat_units_vect, error_PGT, col = "red", pch = 2)
    if(population){ points(n_stat_units_vect, error_GT)}
  }
  
}

## GT
pdf(file = paste(path_images, "GT.pdf"), height = 5, width = 10)
plot_population_convergence(Beta_list_GT, R_list_GT, W_list_GT, population = FALSE)
dev.off()

## MV-PLS
pdf(file = paste(path_images, "MV_PLS.pdf"), height = 5, width = 10)
plot_population_convergence(Beta_list_MV_PLS, R_list_MV_PLS, W_list_MV_PLS)
dev.off()

## fPLS
pdf(file = paste(path_images, "fPLS.pdf"), height = 5, width = 10)
plot_population_convergence(Beta_list_fPLS, R_list_fPLS, W_list_fPLS)
dev.off()
