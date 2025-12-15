# % %%%%%%%%%%%%%%%%%%% %
# % % Test: fPLS - 2D % %
# % %%%%%%%%%%%%%%%%%%% %

rm(list = ls())
graphics.off()

## global variables ----

test_suite <- "2D"
TEST_SUITE <- "fPLS - 2D"

mode_MV <- "PLS-R"
mode_fun <- "fPLS-R"


# libraries ----

## fda
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
FORCE_EVALUATE <- TRUE

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
colors <- c(brewer.pal(3, "Greys")[3], brewer.pal(3, "Blues")[2:3])

## names and labels
names_models <- c("MV_PLS", "fPLS_off", "fPLS_gcv")
lables_models <- c("MV-PLS", "fPLS (no calibration)", "fPLS (gcv calibration)")

## resolution of the high resolution grid
n_nodes_HR <- 80*80


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
  args[1] <- "test1"
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
generated_domain <- generate_domain(name_mesh, n_nodes)
domain <- generated_domain$domain
domain_boundary <- generated_domain$domain_boundary
loadings_true_generator <- generated_domain$loadings_true_generator
mesh <- generated_domain$mesh

## HR mesh
generated_domain <- generate_domain(name_mesh, n_nodes_HR)
mesh_HR <- generated_domain$mesh
nodes_HR <- mesh_HR$nodes

## locations
locations <- generate_locations(generated_domain, locs_eq_nodes, n_locs)
locations <- closest_2D(nodes_HR, locations)

plot(nodes_HR, asp = 1, pch = 20)
points(locations, col = "red", pch = 20)

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
        mesh_HR, locations,
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
