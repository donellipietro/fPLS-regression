## options
options_grid <- list()
for(name_ao in name_aggregation_option_vect) {
  options_grid[[name_ao]] <- unique(data_plot[, name_ao])
}

## plots
for(i in c(1)) {
  
  boxplot_list <- list()
  plot_list <- list()
  
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
  
  if(is.null(limits)) {
    limits <- c(0, max(data_plot[, names_models], na.rm = TRUE))
  }
  
  labels_rows <- paste(labels_options_selected[1], "=", options_grid_selected[[1]])
  labels_cols <- paste(labels_options_selected[2], "=", options_grid_selected[[2]])
  
  for(j in 1:nrow(combinations_options)) {
    
    ## data preparation
    data_plot_trimmed <- data_plot[data_plot[,names_options_selected[1]] == combinations_options[j, 1] & 
                                     data_plot[,names_options_selected[2]] == combinations_options[j, 2],
                                   c(name_aggregation_option, names_models)]
    colnames(data_plot_trimmed) <- c("Group", names_models)
    indexes <- which(!is.nan(colSums(data_plot_trimmed[, names_models])))
    # data_plot_aggregated <- aggregate(. ~ Group, data = data_plot_trimmed[, c("Group", names(indexes))], FUN = median)
    
    ## plots
    boxplot_list[[j]] <- plot.grouped_boxplots(
      data_plot_trimmed[, c("Group", names(indexes))],
      values_name = NULL,
      group_name = group_name,
      subgroup_name = "Approches",
      subgroup_labels = lables_models[indexes],
      subgroup_colors = colors[indexes],
      limits = limits,
      LEGEND = FALSE
    ) + standard_plot_settings()
  }
  boxplot <- arrangeGrob(grobs = boxplot_list, ncol = length(labels_cols))
  boxplot <- labled_plots_grid(boxplot, title, labels_cols, labels_rows, 7, 7)
  grid.arrange(boxplot)
  
}