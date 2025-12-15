# utils ----

## settings

standard_plot_settings <- function() {
  standard_plot_settings <- theme_bw() +
    theme(
      text = element_text(size = 12),
      plot.title = element_text(
        color = "black",
        face = "bold",
        size = 14,
        hjust = 0.5,
        vjust = 1
      )
    )
}

standard_plot_settings_fields <- function() {
  standard_plot_settings_fiends <- theme_minimal() +
    theme(
      text = element_text(size = 12),
      plot.title = element_text(
        color = "black",
        face = "bold",
        size = 14,
        hjust = 0.5,
        vjust = 1
      ),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(angle = 45, hjust = 1)
    )
}


# points ----

plot.points <- function(locations, group = NULL, boundary = NULL,
                        group_name = "Groups", group_colors = NULL, group_labels = NULL,
                        size = 1, LEGEND = FALSE) {
  ## assemble data
  if (is.null(group)) {
    group <- rep(0, nrow(locations))
  }
  data <- data.frame(locations, group)
  colnames(data) <- c("x", "y", "Group")
  
  ## grouping
  group_levels <- unique(data$Group)
  if (is.null(group_labels)) {
    group_labels <- group_levels
  }
  
  ## add association between subgroup_labels and subgroup_colors
  if (is.null(group_colors)) {
    if (length(group_labels) > 1) {
      group_colors <- rainbow(length(group_labels))
    } else {
      group_colors <- c("#000000")
    }
  }
  names(group_colors) <- group_labels
  
  ## refactor categorical variables
  data <- data %>%
    mutate(Group = factor(Group, levels = group_levels, labels = group_labels))
  
  ## plot
  plot <- ggplot() +
    geom_point(data = data, aes(x = x, y = y, color = Group), size = size) +
    coord_fixed() +
    scale_color_manual(
      name = group_name,
      values = group_colors,
      labels = group_labels
    )
  
  if (!is.null(boundary)) {
    plot <- plot +
      geom_polygon(data = fortify(boundary), aes(x = long, y = lat), fill = "transparent", color = "black", linewidth = 1)
  }
  
  ## add a legend if required
  if (LEGEND) {
    # plot <- plot + guides(fill=guide_legend(nrow=1, byrow=TRUE))
  } else {
    plot <- plot + guides(color = "none")
  }
  
  return(plot)
}


## plot curves ----

plot.curve <- function(locations, f, true = NULL,  limits = NULL){
  
  data <- data.frame(locations, y = f)
  colnames(data) <- c("x", "y")
  
  plot <- ggplot() +
    geom_line(data = data, aes(x = x, y = y)) # +
  # coord_fixed()
  
  if(!is.null(true)) {
    data_true <- data.frame(locations, y = true)
    colnames(data_true) <- c("x", "y")
    plot <- plot +
      geom_line(data = data_true, aes(x = x, y = y), linetype = "dotted")
  }
  
  if(!is.null(limits)) {
    plot <- plot + ylim(limits[1], limits[2])
  }
  
  return(plot)
}

plot.curve_points <- function(locations, f, true = NULL, size = 1, limits = NULL){
  
  data <- data.frame(locations, y = f)
  colnames(data) <- c("x", "y")
  
  plot <- ggplot() +
    geom_point(data = data, aes(x = x, y = y), size = size) # +
  # coord_fixed()
  
  if(!is.null(true)) {
    data_true <- data.frame(locations, y = true)
    colnames(data_true) <- c("x", "y")
    plot <- plot +
      geom_line(data = data_true, aes(x = x, y = y))
  }
  
  if(!is.null(limits)) {
    plot <- plot + ylim(limits[1], limits[2])
  }
  
  return(plot)
}


# fields ----

## plot field: points version

plot.field_points <- function(locations, f, boundary = NULL,
                              size = 1, limits = NULL, colormap = "D",
                              discrete = FALSE, LEGEND = FALSE) {
  
  if(is.null(f)) {
    return(ggplot() + theme_void())
  }
  
  data <- data.frame(locations, value = f)
  colnames(data) <- c("x", "y", "value")
  
  plot <- ggplot() +
    geom_point(data = data, aes(x = x, y = y, color = value), size = size) +
    coord_fixed()
  
  if (!discrete) {
    if (is.null(limits)) {
      plot <- plot + scale_color_viridis(
        option = colormap
      )
    } else {
      plot <- plot + scale_color_viridis(
        option = colormap,
        limits = limits
      )
    }
  } else {
    plot <- plot + scale_color_viridis_d(
      option = colormap
    )
  }
  
  if (!is.null(boundary)) {
    plot <- plot +
      geom_polygon(data = fortify(boundary), aes(x = long, y = lat), fill = "transparent", color = "black", linewidth = 1)
  }
  
  ## add a legend if required
  if (LEGEND) {
    # plot <- plot + guides(fill=guide_legend(nrow=1, byrow=TRUE))
  } else {
    plot <- plot + guides(color = "none")
  }
  
  return(plot)
}

## plot field: tile version

plot.field_tile <- function(nodes, f, boundary = NULL,
                            limits = NULL, breaks = NULL, colormap = "D",
                            discrete = FALSE, ISOLINES = FALSE, LEGEND = FALSE) {
  
  if(is.null(f)) {
    return(ggplot() + theme_void())
  }
  
  data <- data.frame(nodes, value = f)
  colnames(data) <- c("x", "y", "value")
  
  data <- na.omit(data)
  
  plot <- ggplot() +
    geom_tile(data = data, aes(x = x, y = y, fill = value)) +
    coord_fixed()
  
  if(!is.null(breaks) || ISOLINES) {
    color = "black"
    limits_real <- range(data$value)
    if(is.null(breaks)) {
      breaks <- seq(limits_real[1], limits_real[2], length = 10)
    }
    breaks_initial <- breaks
    h = breaks[2]-breaks[1]
    if(limits_real[1] < min(breaks)) {
      breaks <- c(sort(seq(min(breaks), limits_real[1]-h, by = -h)[-1]), breaks)
    }
    if(limits_real[2] > max(breaks)) {
      breaks <- c(breaks, seq(max(breaks), limits_real[2]+h, by = h)[-1])
    }
    if(length(breaks) > 2*length(breaks_initial) ) {
      breaks <- breaks_initial
      color = "red"
    }
    plot <- plot +
      geom_contour(data = data, aes(x = x, y = y, z = value), color = color, breaks = breaks)
  }
  
  if (!discrete) {
    if(!is.null(breaks)) {
      h = breaks[2]-breaks[1]
      limits <- limits + c(-h, h)
    }
    if (is.null(limits)) {
      plot <- plot + scale_fill_viridis(
        option = colormap # ,
        # labels = scales::scientific_format()
      )
    } else {
      plot <- plot + scale_fill_viridis(
        option = colormap,
        limits = limits # ,
        # labels = scales::scientific_format()
      )
    }
  } else {
    if(length(colormap)==1) {
      plot <- plot + scale_fill_viridis_d(
        option = colormap # ,
        # labels = scales::scientific_format()
      )
    } else {
      plot <- plot + scale_fill_gradientn(
        colours = colormap # ,
        # labels = scales::scientific_format()
      )
    }
  }
  
  if (!is.null(boundary)) {
    plot <- plot +
      geom_polygon(data = fortify(boundary), aes(x = long, y = lat), fill = "transparent", color = "black", linewidth = 1)
  }
  
  ## add a legend if required
  if (LEGEND) {
    # plot <- plot + guides(fill=guide_legend(nrow=1, byrow=TRUE))
  } else {
    plot <- plot + guides(fill = "none")
  }
  
  return(plot)
}



plot.field_tile_gradient <- function(nodes, f, boundary = NULL,
                                     limits = NULL, breaks = NULL, colormap = "D",
                                     discrete = FALSE, ISOLINES = FALSE, LEGEND = FALSE) {
  
  if(is.null(f)) {
    return(ggplot() + theme_void())
  }
  
  data <- data.frame(nodes, value = f)
  colnames(data) <- c("x", "y", "value")
  
  data <- na.omit(data)
  
  plot <- ggplot() +
    geom_tile(data = data, aes(x = x, y = y, fill = value)) +
    coord_fixed()
  
  if(!is.null(breaks) || ISOLINES) {
    color = "black"
    limits_real <- range(data$value)
    if(is.null(breaks)) {
      breaks <- seq(limits_real[1], limits_real[2], length = 10)
    }
    breaks_initial <- breaks
    h = breaks[2]-breaks[1]
    if(limits_real[1] < min(breaks)) {
      breaks <- c(sort(seq(min(breaks), limits_real[1]-h, by = -h)[-1]), breaks)
    }
    if(limits_real[2] > max(breaks)) {
      breaks <- c(breaks, seq(max(breaks), limits_real[2]+h, by = h)[-1])
    }
    if(length(breaks) > 2*length(breaks_initial) ) {
      breaks <- breaks_initial
      color = "red"
    }
    plot <- plot +
      geom_contour(data = data, aes(x = x, y = y, z = value), color = color, breaks = breaks)
  }
  
  if (!discrete) {
    if(!is.null(breaks)) {
      h = breaks[2]-breaks[1]
      limits <- limits + c(-h, h)
    }
    if (is.null(limits)) {
      plot <- plot + scale_fill_gradientn(
        colours = colormap #,
        # labels = scales::scientific_format()
      )
    } else {
      plot <- plot + scale_fill_gradientn(
        colours = colormap,
        limits = limits # ,
        # labels = scales::scientific_format()
      )
    }
  } else {
    plot <- plot + scale_fill_gradientn(
      colours = colormap # ,
      # labels = scales::scientific_format()
    )
  }
  
  if (!is.null(boundary)) {
    plot <- plot +
      geom_polygon(data = fortify(boundary), aes(x = long, y = lat), fill = "transparent", color = "black", linewidth = 1)
  }
  
  ## add a legend if required
  if (LEGEND) {
    # plot <- plot + guides(fill=guide_legend(nrow=1, byrow=TRUE))
  } else {
    plot <- plot + guides(fill = "none")
  }
  
  return(plot)
}


# box-plots ----

## grouped box-plot

# Input format:
# Group | Model1 | ... | ModelN
# Output:
# box-plot gouped by Group with levels {Model1, ..., ModelN}

plot.grouped_boxplots <- function(data,
                                  ## groups options
                                  group_name = "Components", group_labels = NULL,
                                  ## subgroup options
                                  subgroup_name = "Models", subgroup_labels = NULL, subgroup_colors = NULL,
                                  values_name = "Score",
                                  ## limits
                                  limits = NULL,
                                  ## options
                                  DIVIDERS = TRUE, LEGEND = TRUE) {
  ## data integrety check
  if (!("Group" %in% names(data))) stop("The dataframe must contain a column named 'Group'")
  
  ## data reformat
  data <- data %>%
    pivot_longer(cols = -Group, names_to = "SubGroup", values_to = "Score")
  
  ## extract names
  groups_levels <- unique(data$Group)
  subgroup_levels <- unique(data$SubGroup)
  if (is.null(group_labels)) {
    group_labels <- groups_levels
  }
  if (is.null(subgroup_labels)) {
    subgroup_labels <- subgroup_levels
  }
  
  ## add association between subgroup_labels and subgroup_colors
  if (is.null(subgroup_colors)) {
    subgroup_colors <- rainbow(length(subgroup_labels))
  }
  names(subgroup_colors) <- subgroup_labels
  
  ## refactor categorical variables
  data <- data %>%
    mutate(Group = factor(Group, levels = groups_levels, labels = group_labels)) %>%
    mutate(SubGroup = factor(SubGroup, levels = subgroup_levels, labels = subgroup_labels))
  
  ## grouped box-plot
  plot <- ggplot(data = data, aes(x = Group, y = Score, fill = SubGroup)) +
    geom_boxplot(na.rm = TRUE) +
    labs(x = group_name, y = values_name) +
    scale_fill_manual(
      name = subgroup_name,
      values = subgroup_colors,
      labels = subgroup_labels
    )
  
  if (!is.null(limits)) {
    plot <- plot + scale_y_continuous(
      limits = limits
    )
  }
  
  ## add groups divider if required
  if (DIVIDERS) {
    if (length(group_labels) > 1) {
      plot <- plot +
        geom_vline(
          xintercept = seq(1.5, length(unique(group_labels)) - 0.5, 1),
          lwd = 0.2, colour = "grey"
        )
    }
  }
  
  ## add a legend if required
  if (LEGEND) {
    # plot <- plot + guides(fill=guide_legend(nrow=1, byrow=TRUE))
  } else {
    plot <- plot + guides(fill = "none")
  }
  
  return(plot)
}

## multiple lines

# Input format:
# x | y1 | ... | yN

plot.multiple_lines <- function(data,
                                ## x axis options
                                x_name = "Components",
                                x_breaks = NULL,
                                ## subgroup options
                                subgroup_name = "Models", subgroup_labels = NULL, subgroup_colors = NULL,
                                values_name = "Score",
                                ## limits
                                limits = NULL,
                                ## options
                                LOGX = FALSE,
                                LOGY = FALSE,
                                LOGLOG = FALSE,
                                NORMALIZED = FALSE,
                                LEGEND = TRUE) {
  ## rename first column
  columns_names <- colnames(data)
  columns_names[1] <- "x"
  colnames(data) <- columns_names
  
  if(NORMALIZED) {
    for(name in columns_names[-1])
      data[, name] <- data[, name] / min(data[, name])
  }
  
  if(LOGLOG == TRUE) {
    LOGX <- TRUE
    LOGY <- TRUE
  }
  
  if(is.logical(x_breaks) && x_breaks) {
    x_breaks <- unique(data$x)
  }
  
  ## data reformat
  data <- data %>%
    pivot_longer(cols = -x, names_to = "SubGroup", values_to = "Score")
  
  ## extract names
  subgroup_levels <- unique(data$SubGroup)
  if (is.null(subgroup_labels)) {
    subgroup_labels <- subgroup_levels
  }
  
  ## add association between subgroup_labels and subgroup_colors
  if (is.null(subgroup_colors)) {
    subgroup_colors <- rainbow(length(subgroup_labels))
  }
  names(subgroup_colors) <- subgroup_labels
  
  ## refactor categorical variables
  data <- data %>%
    mutate(SubGroup = factor(SubGroup, levels = subgroup_levels, labels = subgroup_labels))
  
  ## grouped lines
  plot <- ggplot(data = data, aes(x = x, y = Score, color = SubGroup)) +
    geom_line(linewidth = 1) +
    labs(x = x_name, y = values_name) +
    scale_color_manual(
      name = subgroup_name,
      values = subgroup_colors,
      labels = subgroup_labels
    )
  
  ## LOGLOG scale and ref
  if (LOGLOG) {
    if(NORMALIZED) {
      x <- seq(min(data$x), max(data$x), length = 10)
      plot <- plot +
        geom_line(data = data.frame(x = x, y = x/x[1]), aes(x = x, y = y), linetype = "dashed", color = "grey", linewidth = 0.3) +
        geom_line(data = data.frame(x = x, y = (x/x[1])^2), aes(x = x, y = y), linetype = "dashed", color = "grey", linewidth = 0.3) +
        geom_line(data = data.frame(x = x, y = (x/x[1])^3), aes(x = x, y = y), linetype = "dashed", color = "grey", linewidth = 0.3)
    }
  }
  if (LOGX) {
    if(!is.null(x_breaks)) {
      plot <- plot +
        scale_x_log10(breaks = x_breaks)
    } else {
      plot <- plot +
        scale_x_log10()
    }
  } else {
    if(!is.null(x_breaks)) {
      plot <- plot +
        scale_x_continuous(breaks = x_breaks)
    } else {
      plot <- plot +
        scale_x_continuous()
    }
  }
  if(LOGY) {
    if (!is.null(limits)) {
      plot <- plot + scale_y_log10(
        limits = limits
      )
    } else {
      plot <- plot + scale_y_log10()
    }
  } else {
    if (!is.null(limits)) {
      plot <- plot + scale_y_continuous(
        limits = limits
      )
    }
  }
  
  ## add a legend if required
  if (LEGEND) {
    # plot <- plot + guides(fill=guide_legend(nrow=1, byrow=TRUE))
  } else {
    plot <- plot + guides(color = "none")
  }
  
  return(plot)
}

## grids of plots

labled_plots_grid <- function(plot, title = NULL, labels_cols = NULL, labels_rows = NULL, height = 4, width = 4) {
  
  ## grid dimensions
  n_row <- max(plot$layout$t)
  n_col <- max(plot$layout$l)
  
  ## labels columns
  if(!is.null(labels_cols)) {
    labels_grobs_cols <- list()
    for(col in 1:length(labels_cols)) {
      label_cols <- labels_cols[[col]]
      labels_grobs_cols[[col]] <- textGrob(label_cols, gp = gpar(fontsize = 12, fontface = 'bold'))
    }
    labels_grobs_cols <- arrangeGrob(grobs = labels_grobs_cols, nrow = 1)
  }
  
  ## labels rows
  add <- 0
  if(!is.null(labels_rows)) {
    labels_grobs_rows <- list()
    if(!is.null(labels_cols)) {
      add <- 1
      labels_grobs_rows[[1]] <- textGrob(" ", gp = gpar(fontsize = 12, fontface = 'bold'),)
    }
    for(row in 1:length(labels_rows)+add) {
      label_row <- labels_rows[[row-add]]
      labels_grobs_rows[[row]] <- textGrob(label_row, gp = gpar(fontsize = 12, fontface = 'bold'), rot = 90)
    }
    labels_grobs_rows <- arrangeGrob(grobs = labels_grobs_rows, ncol = 1, heights = c(1, rep(height, n_row)))
  }
  
  ## title
  if(!is.null(title)) {
    title_grob <- textGrob(title, gp = gpar(fontsize = 14, fontface = 'bold'))
  }
  
  ## add labels
  if(!is.null(labels_cols)) {
    plot <- arrangeGrob(labels_grobs_cols, plot, heights = c(1, height * n_row))
  }
  if(!is.null(labels_rows)) {
    plot <- arrangeGrob(labels_grobs_rows, plot, widths = c(1, width * n_col))
  }
  
  ## add title
  if(!is.null(title)) {
    plot <- arrangeGrob(title_grob, plot, heights = c(1, add + height * n_row))
  }
  
  return(plot)
}
