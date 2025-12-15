# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Application Brain: Analysis %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
graphics.off()

## logs
if(!dir.exists("logs")){ dir.create("logs") }

## check arguments passed by terminal, set default if not present
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  cat("Test defaulted to overall_mean\n")
  args[1] <- "overall_mean"
  cat("Index defaulted to 1\n")
  args[2] <- 1
}

name_test <- args[1]
index <- as.numeric(args[2])

sink(paste("logs/log_", name_test, "_", index ,".txt", sep = ""))

cat(paste("The test is", name_test ,"\n"))
cat(paste("The index is", index ,"\n"))

## libraries ----

# fdaPDE
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
source("analysis/utils/vtu.R")


## data ----
path_results <- paste("results/brain/sensitivity_analysis/", name_test, "/", sep = "")
load(paste(path_results, "data.RData", sep = ""))

## domain
domain <- fdaPDE2::Mesh(mesh_data)

## functional data
data <- fdaPDE2::functional_regression_data(
  domain = domain, 
  Y = Yp,
  X = Xc
)


## analysis ----

lambda <- fdaPDE2::hyperparameters(lambda_vect[index])
model_fPLS <- fdaPDE2::fPLS(data, center = FALSE)
model_fPLS$fit(lambda)
results_fPLS <- model_fPLS$results
save(results_fPLS, file = paste(path_results, paste("results_fPLS_", sprintf("%0*d", 3, index),".RData", sep=""), sep = ""))

cat("\nCompleted!")
sink()


