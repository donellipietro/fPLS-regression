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


## prerequisite ----
source(paste("tests/", test_suite, "/utils/generate_options.R", sep = ""))


## libraries ----

## json
suppressMessages(library(jsonlite))

## algebraic utils
suppressMessages(library(pracma))

## data visualization
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(viridis))
suppressMessages(library(stringr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))

## sources ----
source("src/utils/directories.R")
source("src/utils/results_management.R")
source("src/utils/plots.R")


## paths ----
path_options <- paste("queue/", sep = "")
path_results <- paste("results/", test_suite, "/", sep = "")
path_images <- paste("images/", test_suite, "/", sep = "")
file_log <- "log.txt"


## options ----

## colors used in the plots
colors <- c(brewer.pal(3, "Greys")[3], brewer.pal(3, "Blues")[2:3], brewer.pal(9, "RdPu")[c(5, 7, 9)])

## names and labels
names_models <- c("MV_PLS", "fPLS_off", "fPLS_gcv",
                  "B_fPLS", "B_fPLS_smooth", "PB_fPLS")
lables_models <- c("MV-PLS", "R1-fPLS (no calibration)", "R1-fPLS (GCV Calibration)",
                   "B-fPLS (Regression Splines)", "B-fPLS (Smoothing Splines)", "PB-fPLS")


selected <- c(1, 4, 5, 6, 3)
colors <- colors[selected]
names_models <- names_models[selected]
lables_models <- lables_models[selected]

## load data ----

## check arguments passed by terminal
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  args[1] <- "test3"
}

## main test name
name_main_test <- args[1]
cat(paste("\nTest selected:", name_main_test, "\n"))
path_results <- paste(path_results, name_main_test, "/", sep = "")
path_images <- paste(path_images, name_main_test, "/", sep = "")
mkdir(path_images)

## figure size
figure_width <- 20
figure_height <- 10

## generate options
generate_options(name_main_test, path_options)

## list of available tests
file_test_vect <- sort(list.files(path_options))

## room for solutions
names_columns <- c("Group", "n_nodes", "n_locs", "n_stat_units", "NSR_Y", "NSR_X", names_models)
empty_df <- data.frame(matrix(NaN, nrow = 0, ncol = length(names_columns)))
colnames(empty_df) <- names_columns
times <- empty_df
rmses <- list()
irmses <- list()
angles <- list()


for (file_test in file_test_vect) {
  ## load specs
  file_json <- paste(path_options, file_test, sep = "")
  parsed_json <- fromJSON(file_json)
  name_test <- parsed_json$test$name_test
  n_nodes <- parsed_json$dimensions$n_nodes
  n_locs <- parsed_json$dimensions$n_locs
  n_stat_units <- parsed_json$dimensions$n_stat_units
  NSR_Y <- parsed_json$noise$NSR_Y
  NSR_X <- parsed_json$noise$NSR_X_lc
  n_reps <- parsed_json$dimensions$n_reps
  
  cat(paste("\nTest ", name_test, ":\n", sep = ""))
  
  ## load batches
  for (i in 1:n_reps) {
    ## laod batch and log if not present
    sink(file_log, append = TRUE)
    tryCatch(
      {
        path_batch <- paste(path_results, name_test, "/", "batch_", i, "/", sep = "")
        load(paste(path_batch, "batch_", i, "_results_evaluation.RData", sep = ""))
      },
      error = function(e) {
        cat(paste("Error in test ", name_test, " - batch ", i, ": ", conditionMessage(e), "\n", sep = ""))
      }
    )
    sink()
    
    ## times
    times <- add_results(
      times,
      c(
        list(n_nodes = n_nodes, n_locs = n_locs, n_stat_units = n_stat_units, NSR_Y = NSR_Y, NSR_X = NSR_X),
        extract_new_results(results_evaluation, names_models, "execution_time")
      ),
      names_columns
    )
    
    ## rmse
    names <- c("Y_reconstruction", "X_reconstruction_locs", "X_reconstruction")
    if(mode_MV == "PLS-R") {
      names <- c(names, "Beta_locs", "Beta")
    }
    for (name in names) {
      rmses[[name]] <- add_results(
        rmses[[name]],
        c(
          list(n_nodes = n_nodes, n_locs = n_locs, n_stat_units = n_stat_units, NSR_Y = NSR_Y, NSR_X = NSR_X),
          extract_new_results(results_evaluation, names_models, c("rmse", name))
        ),
        names_columns
      )
    }
    
    names <- c("Y_space_directions", "X_space_directions_locs",
               "Y_latent_scores", "X_latent_scores",
               "Y_loadings", "X_loadings_locs")
    for (name in names) {
      rmses[[name]] <- add_results(
        rmses[[name]],
        c(
          list(n_nodes = n_nodes, n_locs = n_locs, n_stat_units = n_stat_units, NSR_Y = NSR_Y, NSR_X = NSR_X),
          extract_new_results(results_evaluation, names_models, c("rmse", name))
        ),
        names_columns
      )
    }
    
    cat(paste("- Batch", i, "loaded\n"))
  }
  
  ## remove option file
  file.remove(file_json)
}


## analysis ----

## names
if(name_main_test == "test1") {
  name_aggregation_option_vect <- c("n_stat_units", "NSR_Y", "NSR_X")
  name_group_vect <- c("N", "NSR_Y", "NSR_X")
  name_group_labels_vect <- c(
    "number of statistical units (N)",
    "NSR Y", "NSR X last comp"
  )
} else if(name_main_test == "test2") {
  name_aggregation_option_vect <- c("n_locs", "NSR_Y", "NSR_X")
  name_group_vect <- c("S", "NSR_Y", "NSR_X")
  name_group_labels_vect <- c(
    "number of locations (S)",
    "NSR Y", "NSR X last comp"
  )
} else if(name_main_test == "test3") {
  name_aggregation_option_vect <- c("n_nodes", "NSR_Y", "NSR_X")
  name_group_vect <- c("K", "NSR_Y", "NSR_X")
  name_group_labels_vect <- c(
    "number of nodes (K)",
    "NSR Y", "NSR X last comp"
  )
} 

### time complexity ----

## open a pdf where to save the plots
pdf(paste(path_images, "time_complexity.pdf", sep = ""), width = figure_width, height = figure_height)

## data and titles
data_plot <- times
values_name <- "Time [seconds]"
title_vect <- paste(
  "Execution times w.r.t the",
  name_group_labels_vect
)

## options
options_grid <- list()
for (name_ao in name_aggregation_option_vect) {
  options_grid[[name_ao]] <- unique(data_plot[, name_ao])
}


for (i in c(1)) {
  boxplot_list <- list()
  
  name_aggregation_option <- name_aggregation_option_vect[i]
  group_name <- name_group_vect[i]
  title <- title_vect[i]
  
  options_grid_selected <- options_grid
  options_grid_selected[[name_aggregation_option]] <- NULL
  names_options_selected <- names(options_grid_selected)
  labels_options_selected <- name_group_vect[-i]
  mg <- meshgrid(options_grid_selected[[1]], options_grid_selected[[2]])
  combinations_options <- data.frame(as.vector(mg[[1]]), as.vector(mg[[2]]))
  colnames(combinations_options) <- names_options_selected
  
  limits <- c(0, max(data_plot[, names_models]))
  range <- range(times[, names_models])
  
  labels_rows <- paste(labels_options_selected[1], "=", options_grid_selected[[1]])
  labels_cols <- paste(labels_options_selected[2], "=", options_grid_selected[[2]])
  
  for (j in 1:nrow(combinations_options)) {
    ## data preparation
    data_plot_trimmed <- data_plot[
      data_plot[, names_options_selected[1]] == combinations_options[j, 1] &
        data_plot[, names_options_selected[2]] == combinations_options[j, 2],
      c(name_aggregation_option, names_models)
    ]
    colnames(data_plot_trimmed) <- c("Group", names_models)
    indexes <- which(!is.nan(colSums(data_plot_trimmed[, names_models])))
    data_plot_aggregated <- aggregate(. ~ Group, data = data_plot_trimmed, FUN = median)
    
    ## plots
    boxplot_list[[j]] <- plot.grouped_boxplots(
      data_plot_trimmed[, c("Group", names(indexes))],
      values_name = NULL, # values_name,
      group_name = group_name,
      subgroup_name = "Approches",
      subgroup_labels = lables_models[indexes],
      subgroup_colors = colors[indexes],
      limits = limits,
      LEGEND = FALSE
    ) + standard_plot_settings()
  }
  boxplot <- arrangeGrob(grobs = boxplot_list, ncol = length(labels_cols))
  boxplot <- labled_plots_grid(boxplot, title, labels_cols, labels_rows, 9, 7)
  grid.arrange(boxplot)
}

dev.off()


### overall quantitative results ----


## open a pdf where to save the plots
pdf(paste(path_images, "overall_quantitative_results.pdf", sep = ""), width = figure_width, height = figure_height)

## Y reconstruction RMSE at locations
data_plot <- rmses[["Y_reconstruction"]][rmses[["Y_reconstruction"]]$Group == "4", ]
title_vect <- paste(
  "Y reconstruction RMSE w.r.t the",
  name_group_labels_vect
)
limits <- NULL
source(paste("tests/", test_suite, "/templates/plot_overall.R", sep = ""))

## X reconstruction RMSE at locations
data_plot <- rmses[["X_reconstruction_locs"]][rmses[["X_reconstruction_locs"]]$Group == "4", ]
title_vect <- paste(
  "X reconstruction RMSE at locations w.r.t the",
  name_group_labels_vect
)
limits <- NULL
source(paste("tests/", test_suite, "/templates/plot_overall.R", sep = ""))

## Beta RMSE at locations
if(mode_MV == "PLS-R") {
  ## Beta 1
  data_plot <- rmses[["Beta_locs"]][rmses[["Beta_locs"]]$Group == "4", ]
  title_vect <- paste(
    "Beta RMSE at locations w.r.t the",
    name_group_labels_vect
  )
  limits <- NULL
  source(paste("tests/", test_suite, "/templates/plot_overall.R", sep = ""))
}

## X reconstruction RMSE at locations
data_plot <- rmses[["X_reconstruction"]][rmses[["X_reconstruction"]]$Group == "4", ]
title_vect <- paste(
  "X reconstruction RMSE w.r.t the",
  name_group_labels_vect
)
limits <- NULL
source(paste("tests/", test_suite, "/templates/plot_overall.R", sep = ""))

## Beta RMSE at locations
if(mode_MV == "PLS-R") {
  ## Beta 1
  data_plot <- rmses[["Beta"]][rmses[["Beta"]]$Group == "4", ]
  title_vect <- paste(
    "Beta RMSE w.r.t the",
    name_group_labels_vect
  )
  limits <- NULL
  source(paste("tests/", test_suite, "/templates/plot_overall.R", sep = ""))
}

dev.off()

### X quantitative results ----

## open a pdf where to save the plots
pdf(paste(path_images, "X_quantitative_results.pdf", sep = ""), width = figure_width, height = figure_height)

name_vect <- c("1", "2", "3", "4")
for (k in 1:4) {
  data_plot <- rmses[["X_space_directions_locs"]][rmses[["X_space_directions_locs"]]$Group == k, ]
  title_vect <- paste(
    "X space directions at locations RMSE component", name_vect[k], "w.r.t the",
    name_group_labels_vect
  )
  limits <- NULL
  source(paste("tests/", test_suite, "/templates/plot_overall.R", sep = ""))
}
for (k in 1:4) {
  data_plot <- rmses[["X_latent_scores"]][rmses[["X_latent_scores"]]$Group == k, ]
  title_vect <- paste(
    "X latent scores RMSE component", name_vect[k], "w.r.t the",
    name_group_labels_vect
  )
  limits <- NULL
  source(paste("tests/", test_suite, "/templates/plot_overall.R", sep = ""))
}
for (k in 1:4) {
  data_plot <- rmses[["X_loadings_locs"]][rmses[["X_loadings_locs"]]$Group == k, ]
  title_vect <- paste(
    "X loadings at locations RMSE component", name_vect[k], "w.r.t the",
    name_group_labels_vect
  )
  limits <- NULL
  source(paste("tests/", test_suite, "/templates/plot_overall.R", sep = ""))
}

dev.off()


### Y quantitative results ----

## open a pdf where to save the plots
pdf(paste(path_images, "Y_quantitative_results.pdf", sep = ""), width = figure_width, height = figure_height)

name_vect <- c("1", "2", "3", "4")
for (k in 1:4) {
  data_plot <- rmses[["Y_space_directions"]][rmses[["Y_space_directions"]]$Group == k, ]
  title_vect <- paste(
    "Y space directions RMSE component", name_vect[k], "w.r.t the",
    name_group_labels_vect
  )
  limits <- NULL
  source(paste("tests/", test_suite, "/templates/plot_overall.R", sep = ""))
}
for (k in 1:4) {
  data_plot <- rmses[["Y_loadings"]][rmses[["Y_loadings"]]$Group == k, ]
  title_vect <- paste(
    "Y loadings RMSE component", name_vect[k], "w.r.t the",
    name_group_labels_vect
  )
  limits <- NULL
  source(paste("tests/", test_suite, "/templates/plot_overall.R", sep = ""))
}

dev.off()


## Plot paper

standard_plot_settings <-function() {
  standard_plot_settings <- theme_bw() +
    theme(
      text = element_text(size = 18),
      plot.title = element_text(
        color = "black",
        face = "bold",
        size = 14,
        hjust = 0.5,
        vjust = 1
      )
    )
}

if(name_main_test == "test3") {

  pdf(paste(path_images, "paper.pdf", sep = ""), width = 8, height = 5)
  
  index_NSR_Y <- rmses$Y_reconstruction$NSR_Y == 0.1
  index_NSR_X <- rmses$Y_reconstruction$NSR_X == 3
  index_group <- rmses$Y_reconstruction$Group == 4
  data <- rmses$Y_reconstruction[index_NSR_Y & index_NSR_X & index_group, c(2,7, 10:12, 9)]
  data_colnames <- colnames(data)
  data_colnames[1] <- "Group"
  colnames(data) <- data_colnames
  plot.grouped_boxplots(data, group_name = "",
                        subgroup_colors = colors[c(1, 4:6, 3)], values_name = "",
                        subgroup_labels = lables_models[c(1, 4:6, 3)],
                        LEGEND = FALSE) + standard_plot_settings()
  
  
  index_NSR_Y <- rmses$X_reconstruction_locs$NSR_Y == 0.1
  index_NSR_X <- rmses$X_reconstruction_locs$NSR_X == 3
  data <- rmses$X_reconstruction_locs[index_NSR_Y & index_NSR_X & index_group, c(2,7, 10:12, 9)]
  data_colnames <- colnames(data)
  data_colnames[1] <- "Group"
  colnames(data) <- data_colnames
  plot.grouped_boxplots(data, group_name = "",
                        subgroup_colors = colors[c(1, 4:6, 3)], values_name = "",
                        subgroup_labels = lables_models[c(1, 4:6, 3)],
                        # limits = c(0.05, 0.13),
                        LEGEND = FALSE) + standard_plot_settings()
  
  
  index_NSR_Y <- rmses$Beta_locs$NSR_Y == 0.1
  index_NSR_X <- rmses$Beta_locs$NSR_X == 3
  data <- rmses$Beta_locs[index_NSR_Y & index_NSR_X & index_group, c(2,7, 10:12, 9)]
  data_colnames <- colnames(data)
  data_colnames[1] <- "Group"
  colnames(data) <- data_colnames
  plot.grouped_boxplots(data, group_name = "",
                        subgroup_colors = colors[c(1, 4:6, 3)], values_name = "",
                        subgroup_labels = lables_models[c(1, 4:6, 3)],
                        # limits = c(0, 0.031),
                        LEGEND = FALSE) + standard_plot_settings()
  
  dev.off()

} else if(name_main_test == "test1") {
  
  pdf(paste(path_images, "paper.pdf", sep = ""), width = 8, height = 5)
  
  index_NSR_Y <- rmses$Y_reconstruction$NSR_Y == 0.1
  index_NSR_X <- rmses$Y_reconstruction$NSR_X == 3
  index_group <- rmses$Y_reconstruction$Group == 4
  data <- rmses$Y_reconstruction[index_NSR_Y & index_NSR_X & index_group, c(4, 7, 10:12, 9)]
  data_colnames <- colnames(data)
  data_colnames[1] <- "Group"
  colnames(data) <- data_colnames
  plot.grouped_boxplots(data, group_name = "",
                        subgroup_colors = colors[c(1, 4:6, 3)], values_name = "",
                        subgroup_labels = lables_models[c(1, 4:6, 3)],
                        LEGEND = FALSE) + standard_plot_settings()
  
  
  index_NSR_Y <- rmses$X_reconstruction_locs$NSR_Y == 0.1
  index_NSR_X <- rmses$X_reconstruction_locs$NSR_X == 3
  data <- rmses$X_reconstruction_locs[index_NSR_Y & index_NSR_X & index_group, c(4 ,7, 10:12, 9)]
  data_colnames <- colnames(data)
  data_colnames[1] <- "Group"
  colnames(data) <- data_colnames
  plot.grouped_boxplots(data, group_name = "",
                        subgroup_colors = colors[c(1, 4:6, 3)], values_name = "",
                        subgroup_labels = lables_models[c(1, 4:6, 3)],
                        LEGEND = FALSE) + standard_plot_settings()
  
  
  index_NSR_Y <- rmses$Beta_locs$NSR_Y == 0.1
  index_NSR_X <- rmses$Beta_locs$NSR_X == 3
  data <- rmses$Beta_locs[index_NSR_Y & index_NSR_X & index_group, c(4, 7, 10:12, 9)]
  data_colnames <- colnames(data)
  data_colnames[1] <- "Group"
  colnames(data) <- data_colnames
  plot.grouped_boxplots(data, group_name = "",
                        subgroup_colors = colors[c(1, 4:6, 3)], values_name = "",
                        subgroup_labels = lables_models[c(1, 4:6, 3)],
                        LEGEND = FALSE) + standard_plot_settings()
  
  dev.off()
  
} else if(name_main_test == "test2") {
  
  pdf(paste(path_images, "paper.pdf", sep = ""), width = 8, height = 5)
  
  index_NSR_Y <- rmses$Y_reconstruction$NSR_Y == 0.1
  index_NSR_X <- rmses$Y_reconstruction$NSR_X == 3
  index_group <- rmses$Y_reconstruction$Group == 4
  data <- rmses$Y_reconstruction[index_NSR_Y & index_NSR_X & index_group, c(3, 7, 10:12, 9)]
  data_colnames <- colnames(data)
  data_colnames[1] <- "Group"
  colnames(data) <- data_colnames
  plot.grouped_boxplots(data, group_name = "",
                        subgroup_colors = colors[c(1, 4:6, 3)], values_name = "",
                        subgroup_labels = lables_models[c(1, 4:6, 3)],
                        LEGEND = FALSE) + standard_plot_settings()
  
  
  index_NSR_Y <- rmses$X_reconstruction_locs$NSR_Y == 0.1
  index_NSR_X <- rmses$X_reconstruction_locs$NSR_X == 3
  data <- rmses$X_reconstruction_locs[index_NSR_Y & index_NSR_X & index_group, c(3 ,7, 10:12, 9)]
  data_colnames <- colnames(data)
  data_colnames[1] <- "Group"
  colnames(data) <- data_colnames
  plot.grouped_boxplots(data, group_name = "",
                        subgroup_colors = colors[c(1, 4:6, 3)], values_name = "",
                        subgroup_labels = lables_models[c(1, 4:6, 3)],
                        LEGEND = FALSE) + standard_plot_settings()
  
  
  index_NSR_Y <- rmses$Beta_locs$NSR_Y == 0.1
  index_NSR_X <- rmses$Beta_locs$NSR_X == 3
  data <- rmses$Beta_locs[index_NSR_Y & index_NSR_X & index_group, c(3, 7, 10:12, 9)]
  data_colnames <- colnames(data)
  data_colnames[1] <- "Group"
  colnames(data) <- data_colnames
  plot.grouped_boxplots(data, group_name = "",
                        subgroup_colors = colors[c(1, 4:6, 3)], values_name = "",
                        subgroup_labels = lables_models[c(1, 4:6, 3)],
                        LEGEND = FALSE) + standard_plot_settings()
  
  dev.off()
  
}



