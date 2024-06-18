#' Prepare data for plotting Metric Contrast
#'
#' Select columns and reshape data, bringing it to a format expected by \code{plotlyMC}.
#'
#' If \code{variable} is one of the \code{transformed_metrics} (stored in global variable),
#' it will be log transformed. Infinite values that arise here will be changed to NAs,
#' which will be dropped in the app server before plotting. The app will display a message
#' if that happens
#'
#' @param data a data table containing growth metrics
#' @param choices a list such as one returned by the \code{selectGroups} module,
#'                of variables and their values that will be plotted
#' @param variable name of column to be plotted
#' @param var_col name of column used for coloring dots
#'
#' @return A \code{data.table} where the column containing the two primary variable values
#'         has separated in two, and the column containing the secondary variable
#'         values remained.
#'         If \code{data} contains co-treatment information, it will be preserved.
#'
#' @keywords plugin_plot
#' @seealso \code{plotlyMC} \code{MetricContrast}
#'
#' @export
#'
prepareDataMC <- function(data, 
                          choices, 
                          variable, 
                          var_col) {
  
  # get prettified versions of selected identifiers
  pidfs <- gDRutils::get_prettified_identifiers(simplify = TRUE)
  cell_name <- pidfs[["cellline_name"]]
  drug2_name <- pidfs[["drug_name2"]]
  concentration2_name <- pidfs[["concentration2"]]
  
  # check inputs
  checkmate::assert_string(variable)
  
  checkmate::assert_list(choices)
  checkmate::assert_true(length(choices) == 4)
  checkmate::assert_true(identical(names(choices)[1:2], c("primary", "secondary")))
  checkmate::assert_true(length(choices[[choices$primary]]) == 2)
  checkmate::assert_true(length(choices[[choices$secondary]]) > 0)
  
  checkmate::assert_data_table(data)
  checkmate::assert_choice(choices[["primary"]], names(data))
  checkmate::assert_choice(choices[["secondary"]], names(data))
  checkmate::assert_choice(variable, names(data))
  checkmate::assert_choice(var_col, c(names(data), "none"))
  
  # there are subsets with no co-drug data but still with drug_name_2 and concentration2 columns
  # let's remove columns with codrug data in such cases
  if (!gDRutils::has_valid_codrug_data(data)) {
    data <- gDRutils::remove_codrug_data(data)
  }
  
  # prepare variable combinations as data table
  conditions <- data.table::CJ(choices[[choices$primary]], choices[[choices$secondary]])
  data.table::setnames(conditions, c(choices$primary, choices$secondary))
  # subset data by joining conditions, including co-treatment identifiers
  vars_cotreatment <- intersect(c(drug2_name, concentration2_name), names(data))
  vars_keep <- intersect(names(data),
                         c(names(conditions), vars_cotreatment, var_col, variable))
  data_long <- data[conditions, on = names(conditions)][, vars_keep, with = FALSE]
  # prepare formula for long-to-wide reshaping
  # protect against space in the variable names
  sec_converted <- gsub("(.*)", "`\\1`", choices$secondary)
  prim_converted <- gsub("(.*)", "`\\1`", choices$primary)
  if (var_col != "none") {
    var_col_converted <- gsub("(.*)", "`\\1`", var_col)
    form <- paste(paste(sec_converted, var_col_converted, sep = " + "), prim_converted, sep = " ~ ")
  } else {
    form <- paste(paste(sec_converted, sep = " + "), prim_converted, sep = " ~ ")
  }
  if (length(vars_cotreatment) > 0) {
    # protect against space in the variable names
    vars_cotreatment <- gsub("(.*)", "`\\1`", vars_cotreatment)
    form <- paste(paste(vars_cotreatment, collapse = " + "), form, sep = " + ")
  }
  form <- stats::as.formula(form)
  # reshape data so that primary variable is split and secondary variable is kept together
  data_wide <- data.table::dcast(data_long, form, value.var = variable)
  # check output
  checkmate::assert(
    checkmate::check_data_table(data_long),
    checkmate::check_false(choices$primary %in% names(data_wide)),
    checkmate::check_choice(choices$secondary, choices = names(data_wide)),
    all(vapply(choices[[choices$primary]], 
               function(x) checkmate::check_choice(x, names(data_wide)), logical(1))),
    combine = "and"
  )
  
  return(data_wide)
}



#' Draw interactive Metric Contrast plot
#'
#' Build a an intereactive plot for comparing a metric between two groups.
#'
#' Draws a scatter plot, where one metric is compared between different conditions.
#' Depending on the variables assigned to the axes (one drug or one cell line each),
#' the points correspond to the values of the selected metric for those cell lines / drugs.
#' If cell lines are put on axes, drugs will be color coded and vice versa.
#'
#' Special rules apply for GR50, IC50 and EC50.
#' These variables are plotted in log-scale.
#' The argument \code{axis_limits} determines whether the Y axis limits will be
#' fixed at 1e-4 to 100 (fixed), derived from data (free), or whether one
#' limit will be fixed and the other derived from the data (upper and lower).
#'
#' Note that in the app axes are editable and the user can alter the range by dragging the axis.
#'
#' On request from shiny app an identity line (\code{y = x}) will be drawn.
#' Likewise, a Spearman correlation test can be run on the data and its results can be
#' printed on the plot. Unless \code{metric} is one of the *50s, a linear fit will also be plotted.
#'
#' The plot is generated with \code{\link[plotly]{plot_ly}} and has some interactivity.
#'
#' @param data a data table containing metrics for various conditions
#'             (cell lines treated with drugs)
#' @param var_x,var_y character strings assigning variables to axes
#' @param var_txt character string assigning variable to tooltip
#' @param var_col character string assigning variable to coloring
#' @param metric character string representing name of growth metric that is plotted
#' @param identity logical flag specifying whether to add an identity line
#' @param correlation logical flag specifying whether to print statistics
#'                    of a Spearman correlation test; a linear fit is also added for
#'                    the metrics that are plotted on linear axes
#' @param axis_limits character string; fix or release limits of Y axis, see \code{Details}
#' @param source character string representing a match of this string's value 
#'     with the source argument in \code{event_data()}
#' @param with_labelR logical flag whether to enabale \code{plotlyLabelR} functionality
#'
#' @return a plotly object
#'
#' @keywords plugin_plot
#' @seealso \code{prepareDataMC} \code{MetricContrast}
#'
#' @export
#'
plotlyMC <- function(data, 
                     var_x, 
                     var_y, 
                     var_col, 
                     var_txt,
                     metric = "growth metric", 
                     identity = FALSE, 
                     correlation = FALSE,
                     axis_limits = c("fixed", "upper", "lower", "free"),
                     source = "metricContrast",
                     with_labelR = FALSE) {
  
  checkmate::assert_data_table(data)
  checkmate::assert_string(var_x)
  checkmate::assert_string(var_y)
  checkmate::assert_string(var_txt)
  checkmate::assert_string(metric)
  checkmate::assert_flag(identity)
  checkmate::assert_flag(correlation)
  checkmate::assert_choice(var_x, names(data))
  checkmate::assert_choice(var_y, names(data))
  checkmate::assert_choice(var_col, c(names(data), "none"))
  checkmate::assert_choice(var_txt, names(data))
  axis_limits <- match.arg(axis_limits)
  checkmate::assert_string(source)
  checkmate::assert_flag(with_labelR)
  
  # get prettified versions of selected identifiers
  pidfs <- gDRutils::get_prettified_identifiers(simplify = TRUE)
  concentration2_name <- pidfs[["concentration2"]]
  
  if (nrow(data) == 0) {
    return(plotly::plotly_empty())
  }
  
  # check for presence of the single codrug data
  has_codrug_data <- gDRutils::has_single_codrug_data(names(data))
  
  # should the legend be shown ?
  do_show_legend <- do_show_legend(var_col, data)
  
  # prepare axis options and plot title
  plot_title <- sprintf("Contrast of %s\nbetween %s and %s", metric, var_x, var_y)
  if (is.element(metric, get_metrics_to_transform())) {
    axis_type <- "log"
    axis_ticks <- 10 ^ (-4:4)
    data_limits <- switch(axis_limits,
                          "fixed" = c(1e-4, 100) * c(1 / 2, 2),
                          "upper" = c(NA, 100) * c(1 / 2, 2),
                          "lower" = c(1e-4, NA) * c(1 / 2, 2),
                          "free" = c(NA, NA))
    plot_limits <- log10(data_limits)
    plot_title <- sub("50", "<sub>50</sub> [&mu;M]", plot_title)
  } else {
    axis_type <- "linear"
    axis_ticks <- switch(metric, "E max" = -10:10 / 4, "GR max" = -10:10 / 2, -10:10 / 5)
    data_limits <- switch(metric,
                          "E max" = c(0, 1.1),
                          "GR max" = c(-1, 1.1),
                          c(-0.1, max(data[, c(var_x, var_y), with = FALSE], na.rm = TRUE)))
    data_limits[2] <- max(data_limits[2], 0.1, na.rm = TRUE)
    plot_limits <- data_limits * c(0.9, 1.1) ^ sign(data_limits)
  }
  
  # prepare label for tooltips
  data$label <- build_label(data, "contrast")
  
  # adding info about color variable
  if (var_col != "none") {
    data$label <- paste0(data$label, "\n", var_col, ": ", data[[var_col]])
  }
  
  # build plot
  if (has_codrug_data && var_col == "none") {
    
    val_for_symbols <- as.factor(data[[concentration2_name]])
    # Plot symbols for "Concentration_2
    plot_base <- 
      plotly::plot_ly(source = source, customdata = data[[var_txt]],
                      x = data[[var_x]], y = data[[var_y]], text = data[["label"]],
                      type = "scatter", mode = "markers",
                      hoverinfo = "text", showlegend = do_show_legend,
                      symbol = val_for_symbols,
                      symbols = .get_symbol_list()[seq_along(levels(val_for_symbols))],
                      marker = list(color = "#5A5A5A", size = 11))

  } else if (has_codrug_data && var_col != "none") {
    val_for_symbols <- as.factor(data[[concentration2_name]])
    
    plot_base <- 
      plotly::plot_ly(source = source, customdata = data[[var_txt]], colors = "Set1", 
                      type = "scatter", mode = "markers")
    # Add legend for "Concentration_2"
    plot_base <- 
      plotly::add_trace(plot_base,
                        x = data[[var_x]], y = data[[var_y]],
                        symbol = val_for_symbols,
                        symbols = .get_symbol_list()[seq_along(levels(val_for_symbols))],
                        hoverinfo = "none", showlegend = TRUE, legendgroup = concentration2_name,
                        marker = list(color = "#5A5A5A", size = 0.01))
    # Add legend for var_col
    plot_base <- 
      plotly::add_trace(plot_base,
                        x = data[[var_x]], y = data[[var_y]],
                        color = data[[var_col]],
                        hoverinfo = "none", showlegend = TRUE, legendgroup = var_col,
                        marker = list(size = 0.01))
    # Plot symbols for "Concentration_2"
    plot_base <- 
      plotly::add_trace(plot_base,
                        x = data[[var_x]], y = data[[var_y]], text = data[["label"]],
                        color = data[[var_col]],
                        hoverinfo = "text", showlegend = FALSE,
                        symbol = val_for_symbols,
                        symbols = .get_symbol_list()[seq_along(levels(val_for_symbols))],
                        marker = list(size = 11))
    # Add double legend
    plot_base <- 
      plotly::layout(plot_base,
                     showlegend = do_show_legend,
                     legend = list(
                       orientation = "v",
                       itemsizing = "constant", itemclick = FALSE, tracegroupgap = 30,
                       title = get_legend_title(var_col,
                                                has_codrug_data = has_codrug_data))
      )
    
  } else {
    plot_base <- plotly::plot_ly(x = data[[var_x]], y = data[[var_y]], text = data[["label"]],
                                 color = data[[var_col]], colors = "Set1",
                                 type = "scatter", mode = "markers",
                                 hoverinfo = "text", showlegend = do_show_legend,
                                 source = source, customdata = data[[var_txt]])
  }
  
  # adding plotlyLabelR functionality # nolint start
  if (with_labelR) {
    ## enable adjustable labels
    plot_base <- plotlyLabelR::enable_adjustable_labels(plot_base)
    ## register events with the plot
    plot_base <- plotly::event_register(plot_base, event = "plotly_click")
    # meta + click or alt + click
    plot_base <- plotlyLabelR::lbr_event_register(plot_base, event = "plotly_metaclick")
    # meta + click on label
    plot_base <- plotlyLabelR::lbr_event_register(plot_base, event = "plotly_removemetaclick") 
    # meta + doubleclick
    plot_base <- plotlyLabelR::lbr_event_register(plot_base, event = "plotly_doublemetaclick") 
  } # nolint end
  
  # options for correlation information
  if (correlation) {
    # PART ONE
    
    # plot linear fit
    model_formula <-
      stats::reformulate(sprintf("`%s`", var_x), sprintf("`%s`", var_y))
    # build linear model and add fitted values to curated data
    if (is.na(data_limits[1])) {
      data_limits[1] <- min(c(data[[var_x]], data[[var_y]]))
    }
    
    c_data <-
      data[, c(var_x, var_y), with = FALSE]
    c_data_limits <- data_limits
    if (is.element(metric, get_metrics_to_transform())) {
      c_data <- log10(c_data)
      c_data_limits <- log10(c_data_limits)
    }
    
    linear_fit <- stats::lm(model_formula, data = c_data)
    # use model to predict response within range of plot
    c_predicted <-
      data.table::data.table(V = seq(
        from = c_data_limits[1],
        to = c_data_limits[2],
        length.out = 1000
      ))
    data.table::setnames(c_predicted, var_x)
    c_predicted[[var_y]] <-
      stats::predict(linear_fit, c_predicted)
    predicted <-
      if (is.element(metric, get_metrics_to_transform())) {
        10 ^ c_predicted
      } else {
        c_predicted
      }
    
    # add prediction to plot
    plot_base <-
      plotly::add_lines(plot_base,
                        inherit = FALSE,
                        x = predicted[[var_x]],
                        y = predicted[[var_y]],
                        color = I("black"),
                        showlegend = FALSE,
                        hoverinfo = "none"
      )
    # PART TWO
    # add statistics
    test_dummy <- list(method = "Spearman's rank correlation rho", estimate = NA, p.value = NA)
    cor_test <- tryCatch(
      stats::cor.test(formula = stats::reformulate(sprintf("`%s` + `%s`", var_y, var_x)), 
                      data = data, method = "spearman"),
      error = function(e) return(test_dummy))
    test_info1 <- cor_test$method
    padding_x <- if (nchar(var_x) > 15) {
      "\n  "
    } else {
      ""
    }
    padding_y <- if (nchar(var_y) > 15) {
      "\n  "
    } else {
      ""
    }
    test_info2 <- sprintf("data:\t%s%s%s and %s%s", padding_x, var_x, padding_x, padding_y, var_y)
    test_info3 <- if (!is.na(cor_test$estimate)) {
      sprintf("  rho: \t%.3f", cor_test$estimate)
    } else {
      "  unable to quantify correlation"
    }
    test_info4 <- if (!is.na(cor_test$p.value)) {
      sprintf("  p-value: \t%.3f", cor_test$estimate)
    } else {
      "  unable to determine p-value"
    }
    test_info <- c(test_info1, test_info2, test_info3, test_info4)
    test_info <- paste(test_info, collapse = "\n")
    # location of test info
    x_corr <- data_limits[2]
    y_corr <- data_limits[1]
    
    plot_base <-
      plotly::add_text(plot_base,
                       x = x_corr, y = y_corr, text = test_info,
                       type = "scatter", mode = "text",
                       textposition = "top left",
                       textfont = list(size = 11, color = plotly::toRGB("grey50")),
                       showlegend = FALSE, hoverinfo = "none")
  }
  
  # add layout
  plot_laidout <-
    plotly::layout(
      plot_base,
      title = list(text = plot_title, font = list(size = 15)),
      xaxis = list(title = var_x, showgrid = FALSE, range = plot_limits,
                   type = axis_type, tickvals = axis_ticks, exponentformat = "e",
                   showline = TRUE, zeroline = TRUE),
      yaxis = list(title = var_y, showgrid = FALSE, range = plot_limits,
                   type = axis_type, tickvals = axis_ticks, exponentformat = "e",
                   scaleanchor  = "x",
                   showline = TRUE, zeroline = TRUE),
      legend = list(orientation = "v",
                    title = get_legend_title(var_col, has_codrug_data = has_codrug_data))
    )
  
  # options for identity line
  if (identity) {
    # draw the identity line for the observed range of
    # 'x' and 'y' variables in the 'data'
    data_range_v <- range(c(data[[var_x]], data[[var_y]]))
    data_range <- diff(data_range_v)
    
    # set minimum to -0.01 if lowest value is > -0.01
    # to make identity line always crossing (0,0)
    data_vals <-
      c(min(c(-0.01, data_range_v[1] - 0.1 * data_range)),
        data_range_v[2] + 0.1 * data_range)
    
    plot_laidout <- plotly::add_segments(
      p = plot_laidout,
      x = data_vals[1],
      xend = data_vals[2],
      y = data_vals[1],
      yend = data_vals[2],
      line = list(dash = "dash", color = "rgba(220,220,220,1)"),
      showlegend = FALSE,
      hoverinfo = "none"
    )
  }
  
  # modify config
  plot_final <- gDR_plotly_config(plot_laidout,
                                  edits = get_plotly_edits(),
                                  showAxisRangeEntryBoxes = FALSE)
  
  return(plot_final)
  
}

#' List of symbols for plotly
#'
#' @description Since the markers will be colored, it is important that they are shapes with a area.
#' \code{plotly} has own supported list and markers with area are type \code{-dot},
#' but not \code{-open-dot}
#'
#' vals <- plotly::schema(FALSE)$traces$scatter$attributes$marker$symbol$values
#' vals <- vals[!grepl("-open-dot", vals)]
#' vals <- grep("-dot$", vals, value = TRUE)
#'
#' Note: item "hash-dot" was removed from list as one without area.
#'
#' @seealso [Working with symbols](https://plotly-r.com/working-with-symbols)
#'
#' @keywords internal
.get_symbol_list <- function() {
  c("circle-dot", "square-dot", "diamond-dot", "cross-dot", "x-dot", "triangle-up-dot",
    "triangle-down-dot", "triangle-left-dot", "triangle-right-dot", "triangle-ne-dot", "triangle-se-dot",
    "triangle-sw-dot", "triangle-nw-dot", "pentagon-dot", "hexagon-dot", "hexagon2-dot", "octagon-dot",
    "star-dot", "hexagram-dot", "star-triangle-up-dot", "star-triangle-down-dot", "star-square-dot",
    "star-diamond-dot", "diamond-tall-dot", "diamond-wide-dot"
  )
}
