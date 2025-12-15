# Mean estimation ----

## Directory ----

directory.results_path <- paste(directory.results, "mean_estimation/", sep = "")
if (!file.exists(directory.results_path)){
  dir.create(directory.results_path)
}


## Mean estimation ----

means <- list()

# Overall mean
means[["overall"]] <- method(X, domain, lambdas)
FEMobj <- FEM(as.numeric(means[["overall"]]$fitted_nodes), FEMbasis)
write.vtu(FEMobj, paste(directory.results_path, "mean_", name, "_overall.vtu", sep = ""))

# CONTROL mean
means[["CONTROL"]] <- method(X[index.CONTROL,], domain, lambdas)
FEMobj <- FEM(as.numeric(means[["CONTROL"]]$fitted_nodes), FEMbasis)
write.vtu(FEMobj, paste(directory.results_path, "mean_", name, "_CONTROL.vtu", sep = ""))

# SCHZ mean
means[["SCHZ"]] <- method(X[index.SCHZ,], domain, lambdas)
FEMobj <- FEM(as.numeric(means[["SCHZ"]]$fitted_nodes), FEMbasis)
write.vtu(FEMobj, paste(directory.results_path, "mean_", name, "_SCHZ.vtu", sep = ""))


## Save results ----

save(means,
     file = paste(directory.results_path, "mean_estimation_", name, ".RData", sep = ""))


## Clean ----

rm(means, FEMobj) 