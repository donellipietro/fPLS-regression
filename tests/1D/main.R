# % %%%%%%%%%%%%%%%%%%% %
# % % Test: fPLS - 1D % %
# % %%%%%%%%%%%%%%%%%%% %

rm(list = ls())
graphics.off()

## global variables ----

test_suite <- "1D"
TEST_SUITE <- "fPLS - 1D"

mode_MV <- "PLS-R"
mode_fun <- "fPLS-R"


# libraries ----

## fda
suppressMessages(library(fda))
suppressMessages(library(fdaPDE))
suppressMessages(library(fdaPDE2))

## algebraic utils
suppressMessages(library(pracma))

## statistical utils
suppressMessages(library(MASS))

## data visualization
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(viridis))
suppressMessages(library(stringr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))

## json
suppressMessages(library(jsonlite))

## sampling
suppressMessages(library(sf))
suppressMessages(library(sp))
suppressMessages(library(raster))


# sources ----

## general functions
source("src/utils/cat.R")
source("src/utils/directories.R")
source("src/utils/meshes.R")
source("src/utils/domain_and_locations.R")
source("src/utils/wrappers.R")
source("src/utils/errors.R")
source("src/utils/results_management.R")
source("src/utils/plots.R")
source("src/utils/PLS.R")

## test specific functions
source(paste("tests/", test_suite, "/utils/generate_data.R", sep = ""))
source(paste("tests/", test_suite, "/utils/models_evaluation.R", sep = ""))


# paths ----

path_results <- paste("results/", test_suite, "/", sep = "")
path_images <- paste("images/", test_suite, "/", sep = "")
mkdir(c(path_results, path_images))

path_queue <- paste("queue/", sep = "")
path_logs <- paste("logs/", sep = "")
mkdir(c(path_logs))


# options ----

## force testing even if a fit is already available
FORCE_FIT <- FALSE
FORCE_EVALUATE <- FALSE

## execution flow modifiers
RUN <- list()
RUN$tests <- TRUE
RUN$analysis <- TRUE
RUN$quantitative_analysis <- TRUE
RUN$qualitative_analysis <- FALSE

## global variables
RSTUDIO <- FALSE


## calibration parameters ----
seed <- 0 # for gcv calibration procedure
lambda_fixed <- fdaPDE2::hyperparameters(1e-9)
lambda_grid <- fdaPDE2::hyperparameters(10^seq(-9, 2, by = 0.1))


## visualization options ----

## colors used in the plots
colors <- c(brewer.pal(3, "Greys")[3], brewer.pal(3, "Blues")[2:3], brewer.pal(3, "Purples")[1:3])

## names and labels
names_models <- c("MV_PLS", "fPLS_off", "fPLS_gcv", "B_fPLS", "B_fPLS_smooth", "PB_fPLS")
lables_models <- c("MV-PLS", "fPLS (no calibration)", "fPLS (gcv calibration)", "B-fPLS", "B-fPLS smoothing splines", "PB-fPLS")

## resolution of the high resolution grid
n_knots_HR <- 500


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Test: fPLS - 1D ----
cat.script_title(paste("Test:", TEST_SUITE))


## options ----
cat.section_title("Options")

## check arguments passed by terminal, set default if not present
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  RSTUDIO <- TRUE
  source(paste("tests/", test_suite, "/utils/generate_options.R", sep = ""))
  args[1] <- "test0"
  generate_options(args[1], path_queue)
  args[2] <- sort(list.files(path_queue))[1]
}

## read arguments provided
name_main_test <- args[1]
file_options <- args[2]

## load options
parsed_json <- fromJSON(paste(path_queue, file_options, sep = ""))
name_test <- parsed_json$test$name_test
name_mesh <- parsed_json$mesh$name_mesh
n_nodes <- parsed_json$dimensions$n_nodes
n_locs <- parsed_json$dimensions$n_locs
n_stat_units <- parsed_json$dimensions$n_stat_units
n_comp <- parsed_json$dimensions$n_comp
n_reps <- parsed_json$dimensions$n_reps
locs_eq_nodes <- parsed_json$data$locs_eq_nodes
NSR_X_lc <- parsed_json$noise$NSR_X_lc
NSR_Y <- parsed_json$noise$NSR_Y

## log file
file_log <- paste(path_logs, "log_",name_test,".txt", sep = "")
if(!RSTUDIO){
  sink(file_log, append = TRUE)
}

## options visualization
cat.json(parsed_json)

## create results directory
path_results <- paste(path_results, name_main_test, "/", sep = "")
mkdir(path_results)
path_results <- paste(path_results, name_test, "/", sep = "")
mkdir(path_results)

## create images directory
path_images <- paste(path_images, name_main_test, "/", sep = "")
mkdir(path_images)
path_images <- paste(path_images, "single_tests/", sep = "")
mkdir(path_images)


## test ----
cat.section_title("Test")

## domain
domain <- fdaPDE2::Mesh(unit_segment(n_nodes))
domain_boundary <- c(0,1)
spline_basis <- create.bspline.basis(rangeval = domain_boundary, breaks = domain$nodes)

## HR grid
grid_HR <- seq(0,1, length = n_knots_HR)
Psi_HR <- eval.basis(grid_HR, spline_basis)

## locations
if(locs_eq_nodes) {
  locations <- as.vector(domain$nodes)
  Psi_locs <- eval.basis(locations, spline_basis)
} else {
  locations <- seq(0,1, length = n_locs) + rnorm(n_locs, 0, 0.25/n_locs)
  locations <- closest(grid_HR, locations)
  Psi_locs <- eval.basis(locations, spline_basis)
}


## fit the models n_reps times
if (RUN$tests) {
  for (i in 1:n_reps) {
    ## message
    cat(paste("\nBatch ", i, ":\n", sep = ""))
    
    ## create batch directory
    path_batch <- paste(path_results, "batch_", i, "/", sep = "")
    mkdir(path_batch)
    
    ## generate data if necessary
    file_model_vect <- paste(path_batch, "batch_", i, "_fitted_model_", names_models, ".RData", sep = "")
    if (any(!file.exists(file_model_vect)) || FORCE_FIT || FORCE_EVALUATE) {
      ## generate data
      cat("- Data generation\n")
      generated_data <- generate_data(
        grid_HR, locations,
        n_stat_units = n_stat_units, 
        n_comp = n_comp,
        score_dist = "norm",
        NSR_X_lc = NSR_X_lc,
        NSR_Y = NSR_Y,
        seed = i
      )
      ## assembly functional data
      data <- functional_regression_data(
        domain = domain,
        locations = locations,
        Y = generated_data$data$Y,
        X = generated_data$data$X
      )
      if(name_main_test == "test3") {
        data_B_fPLS <- data
      } else {
        data_B_fPLS <- functional_regression_data(
          domain = list(nodes = seq(0, 1, length = 11)), # min(round(n_locs/4), 40) + 2
          locations = locations,
          Y = generated_data$data$Y,
          X = generated_data$data$X
        )
      }
    }
    
    ## fit models
    source(paste("tests/", test_suite, "/templates/fit_and_evaluate.R", sep = ""))
  }
} else { ## end run tests
  cat("Skipped, relying on the saved results!\n")
}


## results analysis ----
cat.section_title("Results analysis")

## load saved results
if (RUN$analysis) {
  source(paste("tests/", test_suite, "/templates/load_results.R", sep = ""))
}


### quantitative analysis ----
cat.subsection_title("Quantitative analysis")


if (RUN$quantitative_analysis) {
  ## open a pdf where to save plots (quantitative analysis)
  pdf(file = paste(path_images, name_test, "_quantitative.pdf", sep = ""))
  
  ## plots
  source(paste("tests/", test_suite, "/templates/plot_quantitative_results.R", sep = ""))
  
  ## close pdf (quantitative analysis)
  dev.off()
}


### qualitative analysis ----
cat.subsection_title("Qualitative analysis")

## plots
if (RUN$qualitative_analysis) {
  ## open a pdf where to save plots (qualitative analysis)
  if (!RSTUDIO) {
    pdf(file = paste(path_images, name_test, "_qualitative.pdf", sep = ""), width = 10, height = 7)
  }
  
  ## plot
  source(paste("tests/", test_suite, "/templates/plot_qualitative_results.R", sep = ""))
  
  ## close pdf (qualitative analysis)
  dev.off()
}



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


## remove options file
file.remove(paste(path_queue, file_options, sep = ""))

cat("\n\n")

if(!RSTUDIO){
  sink()
}
