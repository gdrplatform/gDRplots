#' Draw interactive Metric Distribution plot
#'
#' Build an interactive plot of metrics distribution.
#'
#' Draws a box- or violin plot, depending on user input, with optionally added points.
#' Outliers are omitted from the boxplot.
#' The color scheme, determined by \code{var_col}, applies to both box/violin and points.
#'
#' Special rules apply for GR50, IC50 and EC50.
#' These variables are plotted in log-scale.
#' The argument \code{axis_limits} determines whether the Y axis limits will be
#' fixed at 1e-4 to 100 (fixed), derived from data (free), or whether one
#' limit will be fixed and the other derived from the data (upper and lower).
#'
#' Note that in the app axes are editable and the user can alter the range by dragging the axis.
#'
#' @param data a data table containing metrics for various conditions
#'             (cell lines treated with drugs)
#' @param var_x,var_y,var_col,var_grp variables of \code{data} to assign to plot elements
#' @param type character string defining type of plot: box or violin
#' @param points character string specifying whether to add points on top of boxes/violins
#' @param axis_limits character string; fix or release limits of Y axis, see \code{Details}
#' @param show_tick_labels boolean; whether or not to show X axis labels
#'
#' @return a plotly object
#'
#' @seealso \code{MetricDistribution}
#'
#' @keywords plugin_plot
#' @export
#'
plotlyMD <- function(data, 
                     var_x, 
                     var_y, 
                     var_grp, 
                     var_col, 
                     type, 
                     points = "off",
                     axis_limits = c("fixed", "upper", "lower", "free"),
                     show_tick_labels) {
  
  checkmate::assert_data_table(data)
  checkmate::assert_string(var_x)
  checkmate::assert_choice(var_x, names(data))
  checkmate::assert_string(var_y)
  checkmate::assert_choice(var_y, names(data))
  checkmate::assert_string(var_grp)
  checkmate::assert_choice(var_grp, c("none", names(data)))
  checkmate::assert_string(var_col)
  checkmate::assert_choice(var_col, c("none", names(data)))
  checkmate::assert_string(type)
  checkmate::assert_choice(type, choices = c("box", "violin"))
  checkmate::assert_string(points)
  checkmate::assert_choice(points, choices = c("on", "off"))
  axis_limits <- match.arg(axis_limits)
  
  if (is.element(var_y, get_metrics_to_transform())) {
    title_x <- paste(var_y, "[&mu;M]")
    title_x <- sub("50", "<sub>50</sub>", title_x)
    type_y <- "log"
    range_y <- switch(axis_limits,
                      "fixed" = log10(c(1e-4, 100) * c(1 / 2, 2)),
                      "upper" = log10(c(NA, 100) * c(1 / 2, 2)),
                      "lower" = log10(c(1e-4, NA) * c(1 / 2, 2)),
                      "free" = c(NA, NA))
    ticks_y <- 10 ^ (-4:4)
    horizontal <- 1
  } else {
    title_x <- var_y
    type_y <- "linear"
    shift_yaxis <- 0.2
    range_y <-
      c(get_axis_min(data[[var_y]], shift_yaxis),
        get_axis_max(data[[var_y]], shift_yaxis))
    ticks_y <- switch(
      var_y,
      "E max" = seq(-2.5, 2.5, shift_yaxis),
      "GR max" = seq(-5, 5, shift_yaxis),
      seq(-2, 2, shift_yaxis)
    )
    horizontal <- 0
  }
  title_y <- switch(var_grp,
                    "none" = var_x,
                    sprintf("%s / %s", var_grp, var_x))
  title_plot <- paste("Distribution of", title_x)
  line_horizontal <- list(type = "line",
                          line = list(width = 1, color = "#A6A6A6"), layer = "below",
                          x0 = 0, x1 = 1, xref = "paper", y0 = horizontal, y1 = horizontal)
  
  # add labels
  data$label <- buildLabel(data, "distribution")
  
  # initialize plot
  plot_base <- plotly::plot_ly()
  # add traces
  if (var_grp == "none" && var_col == "none") { ## no coloring, single trace
    plot_populated <- plotly::add_trace(plot_base,
                                        x = data[[var_x]],
                                        y = data[[var_y]],
                                        text = data[["label"]],
                                        type = type)
    plot_populated <- plotly::layout(plot_populated, xaxis = list(title = title_y))
  } else if (var_grp != "none" && var_col == "none") { ## grouping no colouring
    # add traces
    cols <- c(var_x, var_grp)
    plot_base <- plotly::add_trace(plot_base,
                                   x = t(data[, cols, with = FALSE]),
                                   y = data[[var_y]],
                                   text = data[["label"]],
                                   color = I("#1d6cab"),
                                   type = type)
    if (show_tick_labels) {
      plot_populated <- plotly::layout(plot_base, xaxis = list(type = "multicategory", title = title_y),
                                       showlegend = FALSE)
    } else {
      plot_populated <- plotly::layout(plot_base,
                                       xaxis = list(type = "multicategory", showdividers = FALSE,
                                                    title = list(text = title_y, standoff = 380L)),
                                       showlegend = FALSE)
    }
    
    
  } else { ## with coloring, multiple traces
    # unique values of coloring variables
    colors <- unique(data[[var_col]])
    # prepare an ad hoc color palette TODO this is a temporary solution
    color_palette <- paletteBrew(length(colors), "Accent")
    # add traces
    for (col in seq_along(colors)) {
      # subset data of one color
      data_subset <- data[data[[var_col]] == colors[col], ]
      # add trace
      cols <- c(var_x, var_col)
      plot_base <- plotly::add_trace(plot_base,
                                     x = t(data_subset[, cols, with = FALSE]),
                                     y = data_subset[[var_y]],
                                     text = data_subset[["label"]],
                                     color = I(color_palette[col]),
                                     name = colors[col],
                                     type = type)
    }
    if (show_tick_labels) {
      plot_populated <- plotly::layout(plot_base, xaxis = list(type = "multicategory",
                                                               title = sprintf("%s / %s", var_col, var_x)))
    } else {
      plot_populated <- plotly::layout(plot_base,
                                       xaxis = list(type = "multicategory",
                                                    showdividers = FALSE,
                                                    title = list(text = sprintf("%s / %s", var_col, var_x),
                                                                 standoff = 380L)))
    }
  }
  
  # TODO disabling violin plot until I figure out how to properly draw violins in these log scales
  if (type == "violin" && is.element(var_y, get_metrics_to_transform())) {
    plot_base <- plotly::plotly_empty(type = "violin")
    viol_message <-
      paste0("Violin plots for metrics displayed in log scale<br>",
             "are temporarily unavailable.<br><br>",
             "Enjoy the box plots.")
    plot_populated <- plotly::layout(plot_base,
                                     annotations = list(showarrow = FALSE,
                                                        text = viol_message,
                                                        font = list(size = 20, color = "#808B96")))
  }
  # TODO end
  
  # constant style options
  style_options <- list(marker = list(symbol = "circle-open"),
                        pointpos = 0, jitter = 0.5)
  # variable style options
  point_option <- switch(type,
                         "box" = list(boxpoints = switch(points, "on" = "all", "off" = "none"),
                                      whiskerwidth = 0.2),
                         "violin" = list(points = switch(points, "on" = "all", "off" = FALSE)))
  hover_option <- switch(type,
                         "box" = list(hoverinfo = "text+points+boxes+color"),
                         "violin" = list(hoverinfo = "text+points+color"))
  # collate
  style_arguments <- c(list(p = plot_populated), style_options, point_option, hover_option)
  # apply
  plot_styled <- do.call(plotly::style, style_arguments)
  
  # remove y-range limits for violin plots
  if (type == "violin") {
    range_y <- ""
  }
  plot_laidout <- plotly::layout(plot_styled,
                                 title = list(text = title_plot, font = list(size = 15)),
                                 xaxis = list(showticklabels = show_tick_labels),
                                 yaxis = list(title = title_x, showgrid = FALSE,
                                              type = type_y, range = range_y, tickvals = ticks_y,
                                              exponentformat = "e", zeroline = FALSE),
                                 shapes = line_horizontal)
  # modify config
  plot_final <- gDR_plotly_config(plot_laidout,
                                  edits = get_plotly_edits(),
                                  showAxisRangeEntryBoxes = FALSE)
  
  return(plot_final)
}

#' get min value for the axis based on the data to be visualized
#'
#' get min value for the axis based on the data to be visualized
#'
#' @param data numeric vector
#' @param nearest_shift single number nearest value to which round (floor)
#'
#' @keywords internal
#' @return single number
#'
get_axis_min <- function(data, nearest_shift = 0.1) {
  checkmate::assert_numeric(data)
  checkmate::assert_number(nearest_shift)
  
  my_min <- min(data) - 0.01 * diff(range(data))
  floor(my_min / nearest_shift) * nearest_shift
}

#' get max value for the axis based on the data to be visualized
#'
#' get max value for the axis based on the data to be visualized
#'
#' @param data numeric vector
#' @param nearest_shift single number nearest value to which round (ceiling)
#'
#' @keywords internal
#' @return single number
#'
get_axis_max <- function(data, nearest_shift = 0.1) {
  checkmate::assert_numeric(data)
  checkmate::assert_number(nearest_shift)
  
  my_max <- max(data) + 0.01 * diff(range(data))
  ceiling(my_max / nearest_shift) * nearest_shift
}
