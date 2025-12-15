

## time
indexes <- which(!is.nan(colSums(times[, names_models])))
plot <- plot.grouped_boxplots(
  times[, c("Group", names(indexes))],
  values_name = "Time [seconds]",
  group_name = "",
  group_labels = "",
  subgroup_name = "Approaches",
  subgroup_labels = lables_models[indexes],
  subgroup_colors = colors[indexes]
) + standard_plot_settings() + ggtitle("Time")
print(plot)

## RMSE
name_beta <- c()
title_beta <- c()
if(mode_MV == "PLS-R" || mode_fun == "fPLS-R") {
  name_beta <- c("Beta1_locs", "Beta2_locs")
  title_beta <- c("Beta1_locs", "Beta2_locs")
}
names <- c("Y_reconstruction", "X_reconstruction_locs", name_beta,
           "Y_space_directions", "X_space_directions_locs",
           "X_latent_scores",
           "Y_loadings", "X_loadings_locs")
titles <- c("Y_reconstruction", "X_reconstruction_locs", title_beta,
            "Y_space_directions", "X_space_directions_locs",
            "X_latent_scores",
            "Y_loadings", "X_loadings_locs")
for (i in 1:length(names)) {
  name <- names[i]
  title <- titles[i]
  indexes <- which(!is.nan(colSums(rmses[[name]][, names_models])))
  plot <- plot.grouped_boxplots(
    rmses[[name]][, c("Group", names(indexes))],
    values_name = "RMSE",
    group_name = "",
    group_labels = "",
    subgroup_name = "Approaches",
    subgroup_labels = lables_models[indexes],
    subgroup_colors = colors[indexes]
  ) + standard_plot_settings() + ggtitle(title)
  print(plot)
}