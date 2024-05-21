#' prepare curve plotting data
#'
#' Simulate and plot dose response curves for all existing combinations of Cell Line and Drug.
#'
#' @param metrics a \code{data.table} with drug response metrics taken from a \code{SummarizedExperiment}
#' @param range_x range of concentration for which to predict response; numeric vector of length 2
#' @param density number of concentrations for which to predict response; single numeric
#' 
#' @examples
#' SE <- gDRutils::get_synthetic_data("small")[[1]]
#' dt <- gDRcomponents::convert_se_assay_to_custom_dt(SE, "Metrics")
#' prepareCurves(dt)
#'
#' @return A data.table with predicted response for the given concentration range.
#'
#' Prediction is made for both GR value and Relative Viability,
#' across a \code{density}-long sequence of concentrations within \code{range_x}.
#' The concentrations are generated with \code{logSeq}.
#'
#' @keywords drawResponseOverview
#' @export
#'
prepareCurves <- function(metrics, 
                          range_x = c(1e-3, 50e+0), 
                          density = 100) {
  
  # get prettified identifiers
  pidfs <- gDRutils::get_prettified_identifiers(simplify = TRUE)
  drug_id <- pidfs[["drug"]]
  drug_name <- pidfs[["drug_name"]]
  cell_name <- pidfs[["cellline_name"]]
  drug2_name <- pidfs[["drug_name2"]]
  concentration_name <- pidfs[["concentration"]]
  concentration2_name <- pidfs[["concentration2"]]
  duration <- pidfs[["duration"]]
  data_source <- pidfs[["data_source"]]
  
  checkmate::assert_data_table(metrics)
  checkmate::assert_names(names(metrics), must.include = c(cell_name, drug_name))
  vars_fit_params <- c("GR Inf", "GR 0", "GEC50", "h GR", "E Inf", "E0", "EC50", "h RV")
  checkmate::assert_names(names(metrics), must.include = vars_fit_params)
  checkmate::assert_numeric(range_x, len = 2)
  checkmate::assert_number(density, lower = 1, finite = TRUE)
  
  # prepare list of grouping variables, including co-treatment data
  vars_wish_list <- c(cell_name, drug_name, drug2_name,
                      concentration2_name, duration, data_source)
  vars_id <- intersect(vars_wish_list, names(metrics))
  
  # set concentration range
  concentrations <- logSeq(range_x[1], range_x[2], density)
  
  # safety check 
  # we are expecting single row in 'data' for selected combination of vars  in (vars_id)
  if (NROW(metrics) != nrow(unique(metrics[, vars_id, with = FALSE]))) {
    stop(
      sprintf(
        "sth wrong with the data model: selected combination of vars: '%s' is not unique in the data",
        toString(vars_id)
      )
    )
  }
  
  prediction <- metrics[, .(
    Concentration = concentrations,
    `GR value` = gDRutils::predict_efficacy_from_conc(
      c = concentrations,
      x_inf = `GR Inf`,
      x_0 = `GR 0`,
      ec50 = `GEC50`,
      h = `h GR`
    ),
    `Relative Viability` = gDRutils::predict_efficacy_from_conc(
      c = concentrations,
      x_inf = `E Inf`,
      x_0 = `E0`,
      ec50 = `EC50`,
      h = `h RV`
    )
  ),
  by = vars_id]
  
  data.table::setnames(prediction, "Concentration", concentration_name)
  
  # this column will be used to separate curves of individual co-tratments
  if (all(is.element(vars_wish_list, names(prediction)))) {
    # TODO GDR-2513 # nolint start
    # prediction[["Treatment"]] <-
    #   paste(prediction[[cell_name]], prediction[[drug_name]], prediction[[drug2_name]],
    #         prediction[[concentration2_name]], sep = "__")
    # nolint end
  } else {
    prediction[["Treatment"]] <- paste(prediction[[cell_name]], prediction[[drug_name]], sep = "__")
  }
  # this column will be used for highlighting curves coming from the same primary treatment
  prediction[["Highlight"]] <- paste(prediction[[cell_name]], prediction[[drug_name]], sep = "__")
  
  return(prediction)
}


#' Prepare coordinates of extra features
#'
#' @inheritParams prepareCurves
#' 
#' @examples
#' SE <- gDRutils::get_synthetic_data("small")[[1]]
#' dt <- gDRcomponents::convert_se_assay_to_custom_dt(SE, "Metrics")
#' prepareExtras(dt)
#' 
#' @return
#' List of length 2.
#' The first item is a \code{data.table} with coordinates of points
#' passed to \code{plotly::add_markers}.
#' The second is a named list of line descriptions
#' passed the \code{shapes} argument of \code{plotly::layout}.
#'
#' @keywords drawResponseOverview
#' @export
#'
prepareExtras <- function(metrics, range_x = c(1e-3, 50e+0)) {
  
  # get prettified identifiers
  pidfs <- gDRutils::get_prettified_identifiers(simplify = TRUE)
  drug_name <- pidfs[["drug_name"]]
  cell_name <- pidfs[["cellline_name"]]
  drug2_name <- pidfs[["drug_name2"]]
  duration <- pidfs[["duration"]]
  concentration_name <- pidfs[["concentration"]]
  concentration2_name <- pidfs[["concentration2"]]
  checkmate::assert_data_table(metrics, min.rows = 1, min.cols = 1)
  checkmate::assert_names(names(metrics), must.include = c(cell_name, drug_name))
  vars_fit_params <- c("GR Inf", "GR 0", "GEC50", "h GR", "E Inf", "E0", "EC50", "h RV")
  checkmate::assert_names(names(metrics), must.include = vars_fit_params)
  # prepare list of grouping variables, including co-treatment data
  vars_wish_list <- c(cell_name, drug_name, drug2_name,
                      concentration2_name, duration)
  varsID <- intersect(vars_wish_list, names(metrics))
  checkmate::assert_numeric(range_x, len = 2)
  
  data <- data.table::as.data.table(metrics)
  
  # extract columns of interest and find response for GR50 and IS50;
  points <-
    data[,
         .(GR50,
           `GR value` = gDRutils::predict_efficacy_from_conc(
             c = `GR50`, x_inf = `GR Inf`, x_0 = `GR 0`, ec50 = `GEC50`, h = `h GR`),
           IC50,
           `Relative Viability` = gDRutils::predict_efficacy_from_conc(
             c = `IC50`, x_inf = `E Inf`, x_0 = `E0`, ec50 = `EC50`, h = `h RV`),
           `GR Max`, `E Max`),
         by = varsID]
  # reshape
  points_long <- data.table::melt(points,
                                  id.vars = varsID,
                                  measure.vars = list(Concentration = c("GR50", "IC50"),
                                                      Response = c("GR value", "Relative Viability"),
                                                      MaxEffectiveness = c("GR Max", "E Max")),
                                  variable.name = "metric", variable.factor = FALSE)
  data.table::setnames(points_long, "Concentration", concentration_name)
  points_long[, metric := rep(c("GR value", "Relative Viability"), each = .N / 2)]
  # trace vertical lines from axes to points
  verticals <- points_long[, .(x0 = get(concentration_name), x1 = get(concentration_name),
                               y0 = ifelse(metric == "GR value", -1, 0), y1 = Response),
                           by = c(varsID, "metric")]
  verticals <- lapply(split(verticals, by = c(varsID, "metric"), keep.by = FALSE), as.list)
  # trace horizontal lines from axes to bottoms of curves
  horizontals <- points_long[, .(x0 = range_x[1], x1 = range_x[2], y0 = MaxEffectiveness, y1 = MaxEffectiveness),
                             by = c(varsID, "metric")]
  horizontals <- lapply(split(horizontals, by = c(varsID, "metric"), keep.by = FALSE), as.list)
  # collate all lines
  lines <- c(verticals, horizontals)
  # add common line format
  line_properties <- list(type = "line", line = list(width = 1, dash = 3), opacity = 0.2, layer =  "below")
  lines <- lapply(lines, function(x) c(x, line_properties))
  
  ans <- list(points = points_long, lines = lines)
  return(ans)
}



#' plot reactive response curves
#'
#' Plot all curves in a given data set on a reactive plot.
#'
#' @param curve_data a data.table containing coordinates for plotting predicted dose response curves
#' @param var_y which variable to plot: GR value or Relative Viability
#' @param range_x numeric vector of length 2, specifying range of X axis
#' @param plot_width,plot_height numeric; dimensions of plot
#' 
#' @examples
#' SE <- gDRutils::get_synthetic_data("small")[[1]]
#' dt <- gDRcomponents::convert_se_assay_to_custom_dt(SE, "Metrics")
#' prepared_curves <- prepareCurves(dt)
#' plotlyRCAll(prepared_curves, "GR value")
#'
#' @return A plotly object with highlights. The plot contains only response curves.
#'
#' The app will tracks clicks and create a secondary plot with a more detailed view of the selected curves.
#'
#' @keywords drawResponseOverview
#' @export
#'
plotlyRCAll <- function(curve_data, 
                        var_y, 
                        range_x = c(1e-3, 50e+0),
                        plot_width = 400L,
                        plot_height = 300L) {
  
  checkmate::assert_numeric(plot_width, lower = 1)
  checkmate::assert_numeric(plot_height, lower = 1)
  checkmate::assert_numeric(range_x, len = 2)
  
  # get prettified identifiers
  pidfs <- gDRutils::get_prettified_identifiers(simplify = TRUE)
  drug_name <- pidfs[["drug_name"]]
  cell_name <- pidfs[["cellline_name"]]
  concentration_name <- pidfs[["concentration"]]
  
  checkmate::assert_data_table(curve_data, min.rows = 1, min.cols = 1)
  checkmate::assert_names(names(curve_data), must.include = c(cell_name, drug_name, "Treatment", "Highlight"))
  checkmate::assert_string(var_y)
  checkmate::assert_choice(var_y, choices = c("GR value", "Relative Viability"))
  
  # drug/cell line combination (for exception handling)
  comb_name <- paste0(curve_data[[cell_name]][1], " x ", curve_data[[drug_name]][1])
  
  # add a label column
  curve_data$label <- buildLabel(curve_data, view = "grid")
  # build plot title
  ## this assumes that either Cell Line or Drug has only one value!
  cols <- c(cell_name, drug_name)
  var_long <- gDRcomponents::getLongest(curve_data[, cols, with = FALSE])
  var_short <- if (var_long == cell_name) {
    drug_name
  } else if (var_long == drug_name) {
    cell_name
  } 
  
  plot_title <- sprintf("Dose response curves for %s: %s", var_short, 
                        paste(unique(curve_data[[var_short]]), collapse  = ", "))
  # rename columns of interest
  curve_data <- data.table::copy(curve_data)
  data.table::setnames(curve_data, old = var_y, new = "var_y")
  
  # establish plot limits
  range_y <- switch(var_y,
                    "GR value" = c(-1, 1.1),
                    "Relative Viability" = c(0, 1.1))
  # line properties
  line_horizontal_top <- list(
    type = "line", line = list(width = 1, color = "#A6A6A6"), layer = "below",
    x0 = range_x[1], x1 = range_x[2], y0 = 1, y1 = 1)
  line_horizontal_mid <- list(
    type = "line", line = list(width = 1, color = "#A6A6A6", dash = "dash"), layer = "below",
    x0 = range_x[1], x1 = range_x[2], y0 = 0.5, y1 = 0.5)
  line_horizontal_bot <- list(
    type = "line", line = list(width = 1, color = "black"), layer = "below",
    x0 = range_x[1], x1 = range_x[2], y0 = 0, y1 = 0)
  line_vertical <- list(
    type = "line", line = list(width = 1), layer = "below",
    x0 = 1, x1 = 1, y0 = range_y[1], y1 = range_y[2])
  
  # positions of axis ticks
  ticks_x <- 10 ^ (-10:10)
  ticks_y <- -2:4 / 2
  
  # drop data points that would not fit in the plotting area
  # show exception_data if all data has been dropped
  curve_data <- curve_data[data.table::between(curve_data[["var_y"]], range_y[1], range_y[2]), ]
  if (!NROW(curve_data)) {
    # tmp solution - start
    exception_data <- gDRimport::get_exception_data(36)
    txt_msg <- sprintf(exception_data$sprintf_text, comb_name)
    txt_msg <- gsub("\\. |: ", "\n", txt_msg)
    txt_msg <- gsub(" or ", " \nor ", txt_msg)
    txt_msg <- paste(exception_data$title, "\n\n", txt_msg)
    # tmp solution - end
    plt_err <- plotly::layout(
      p = plotly::plot_ly(mode = "text", type = "scatter"),
      xaxis = list(visible = FALSE), yaxis = list(visible = FALSE),
      annotations = list(text = txt_msg, align = "justify",
                         showarrow = FALSE, font = list(size = 12, color = "darkred")))
    
    return(gDRcomponents::gDR_plotly_config(plt_err))
  }
  
  # add highlight key
  plot_data <- plotly::highlight_key(curve_data, ~ Highlight)
  # build plot
  plot_base <- plotly::plot_ly(plot_data,
                               x = ~ get(concentration_name), y = ~ var_y,
                               split = ~ Treatment, text = ~ label,
                               type = "scatter", mode = "lines",
                               line = list(width = 1), color = I("cornflowerblue"),
                               showlegend = FALSE, hoverinfo = "text", source = "curvePlot",
                               width = plot_width, height = plot_height)
  # add layout options
  plot_laidout <- plotly::layout(plot_base,
                                 title = list(text = plot_title, font = list(size = 14)),
                                 shapes = list(line_horizontal_top, line_horizontal_mid, line_horizontal_bot,
                                               line_vertical),
                                 xaxis = list(title = sprintf("%s [&mu;M]", concentration_name),
                                              range = log10(range_x), type = "log",
                                              showgrid = FALSE, zeroline = FALSE,
                                              tickvals = ticks_x, exponentformat = "e"),
                                 yaxis = list(title = var_y,
                                              range = range_y,
                                              showgrid = FALSE, zeroline = FALSE,
                                              tickvals = ticks_y))
  # add highlights
  plot_highlighted <- plotly::highlight(
    plot_laidout, on = "plotly_hover", off = "plotly_doubleclick",
    persistent = FALSE, color = "red", opacityDim = 0.4)
  # modify config
  plot_final <- gDRcomponents::gDR_plotly_config(plot_highlighted,
                                                 edits = gDRcomponents::get_plotly_edits(),
                                                 showAxisRangeEntryBoxes = FALSE)
  
  return(plot_final)
}


#' plot drug dose response
#'
#' Make an interactive plot of cell fitness against a range of drug concentrations.
#'
#' Any combination of four layers can be drawn, depending on \code{layers}:
#' all measurement values with "observations",
#' average GR value for each drug concentration with "average",
#' error bars (standard deviation) with "error",
#' and the fitted response curve with "curve" (the default).
#'
#' The fifth layer, "extras", highlights points corresponding to GR50/IC50 and GR max/E max.
#' This quasi-layer is composed of points, which are a trace, lines, which are shapes,
#' and labels, which are also a trace. The labels serve as axis ticks for the vertical and horizontal lines.
#' The information necessary to draw the layer is supplied by the \code{extras} argument,
#' see \code{\link{prepareExtras}}.
#' It is a list of length 2,
#' the first item being a \code{data.table} with coordinates of points
#' passed to \code{plotly::add_markers},
#' and the second being a list of line descriptions
#' passed the \code{shapes} argument of \code{plotly::layout}.
#'
#' Curve data is supplied in a separate data frame.
#'
#' @section Recursive plotting:
#' If the data set contains multiple co-drugs, this function is called recursively
#' so that data for each co-drug is plotted on a separate panel. Panel are arranged vertically.
#'
#' @param data a data.table with cell fitness data taken from a \code{SummarizedExperiment}
#' @param var_y variable to plot on Y axis
#' @param layers character vector of layers to plot, see \code{Details}
#' @param curves a data.table with drug response predicted from growth metrics
#' @param range_x numeric vector of length 2, specifying range of X axis
#' @param extras a list of extra features that comprise the fifth layer, see \code{Details}
#' @param plot_width,plot_height numeric; dimensions of plot
#' 
#' @examples
#' SE <- gDRutils::get_synthetic_data("small")[[1]]
#' dt <- gDRcomponents::convert_se_assay_to_custom_dt(SE, "Metrics")
#' prepared_curves <- prepareCurves(dt)
#' plotlyRCSelected(
#'   prepared_curves, 
#'   "GR value",
#'   layers = c("curve", "average", "error"),
#'   curves = prepared_curves,
#'   extras = prepareExtras(dt)
#' )
#' 
#'
#' @return a plotly object
#'
#' @keywords internal
#' @seealso \code{replaceValues}
#'
#' @keywords drawResponseOverview
#' @export
#'
plotlyRCSelected <- function(data,
                             var_y,
                             layers,
                             curves,
                             range_x = c(1e-3, 50e+0),
                             extras,
                             plot_width = 400L,
                             plot_height = 300L) {
  
  checkmate::assert_numeric(plot_width, lower = 1)
  checkmate::assert_numeric(plot_height, lower = 1)
  checkmate::assert_numeric(range_x, len = 2)
  
  # get prettified identifiers
  pidfs <- gDRutils::get_prettified_identifiers(simplify = TRUE)
  drug_name <- pidfs[["drug_name"]]
  cell_name <- pidfs[["cellline_name"]]
  concentration_name <- pidfs[["concentration"]]
  drug2_name <- pidfs[["drug_name2"]]
  concentration2_name <- pidfs[["concentration2"]]
  untreated_tag <- pidfs[["untreated_tag"]]
  
  checkmate::assert_data_table(data, min.rows = 1, min.cols = 1)
  checkmate::assert_string(var_y)
  checkmate::assert_choice(var_y, names(data))
  checkmate::assert_character(layers, null.ok = TRUE)
  lapply(layers, checkmate::assert_choice,
         choices = c("curve", "average", "error", "observations", "extras"), null.ok = TRUE)
  checkmate::assert_data_table(curves, min.rows = 1, min.cols = 1)
  checkmate::assert_names(names(curves), must.include = var_y)
  checkmate::assert_list(extras)
  checkmate::assert_named(extras)
  checkmate::assert_names(names(extras), must.include = c("points", "lines"))
  
  # fast end due to lack of data
  if (is.null(layers)) {
    return(gDRcomponents::gDR_plotly_config(plotly::plotly_empty(type = "scatter", mode = "markers")))
  }
  if (all(is.na(data[[var_y]]))) {
    return(gDRcomponents::gDR_plotly_config(plotly::plotly_empty(type = "scatter", mode = "markers")))
  }
  
  # drug/cell line combination (for exception handling)
  comb_name <-
    paste0(data[[cell_name]][1], " x ", data[[drug_name]][1])
  # TODO GDR-2513 # nolint start
  ### recursive incursion
  # vars_cotreatment <- intersect(c(drug2_name, concentration2_name), names(data))
  # if (length(vars_cotreatment) > 0) {
  #   co_drugs <- setdiff(unique(stats::na.omit(data[[drug2_name]])), untreated_tag)
  #   #  and there is more than one drug
  #   if (length(co_drugs) > 1) {
  #     # replace "untreated" with respective co-drug name
  #     base_untreated_tag <- intersect(untreated_tag, data[[drug2_name]])
  #     data <-
  #       replaceValues(data, drug2_name, base_untreated_tag, co_drugs)
  #     curves <-
  #       replaceValues(curves, drug2_name, base_untreated_tag, co_drugs)
  #     
  #     # split data on co-drugs
  #     data_split <- split(data, data[[drug2_name]])
  #     curves_split <- split(curves, curves[[drug2_name]])
  #     # plot each subset
  #     plot_list <- mapply(FUN = plotlyRCSelected, data = data_split, curves = curves_split,
  #                         MoreArgs = list(var_y = var_y, layers = layers, range_x = range_x, extras = extras),
  #                         SIMPLIFY = FALSE, USE.NAMES = FALSE)
  #     # determine final plot dimensions (all panels in one column)
  #     plot_height <- length(plot_list) * plot_height
  #     
  #     # collate into one subplot
  #     plot_sub <- plotly::subplot(
  #       plot_list, nrows = length(plot_list),
  #       margin = 0.05, titleY = TRUE, titleX = TRUE)
  #     plot_laidout <- plotly::layout(plot_sub, width = plot_width, height = plot_height)
  #     plot_final <- gDRcomponents::gDR_plotly_config(plot_laidout,
  #                                                    edits = gDRcomponents::get_plotly_edits(),
  #                                                    showAxisRangeEntryBoxes = FALSE)
  #     
  #     return(plot_final)
  #   }
  # }
  ### end recursive incursion
  # nolint end
  
  # determine variable to color by
  cols <- c(cell_name, drug_name)
  var_col <- gDRcomponents::getLongest(data[, cols, with = FALSE])
  var_not_col <- if (var_col == cell_name)  {
    drug_name
  } else if (var_col == drug_name)  {
    cell_name
  }
  
  # rename columns of interest
  data <- data.table::copy(data)
  var_std <- paste("Std", var_y)
  data.table::setnames(data, old = c(var_y, var_std, var_col, var_not_col),
                       new = c("var_y", "var_y.err", "var_col", "var_not_col"), skip_absent = TRUE)
  # add label column for tooltips
  data$label <- buildLabel(data, view = "curve")
  
  # build plot title
  plot_title <- sprintf("Drug dose response for %s: %s", var_not_col, 
                        paste(unique(stats::na.omit(data)[["var_not_col"]]), collapse = ", "))
  
  # determine limits of axes
  range_y <- switch(var_y,
                    "GR value" = c(-1, 1.1),
                    "Relative Viability" = c(0, 1.1))
  
  # line properties
  lines <- NULL # reset
  line_horizontal_top <- list(
    type = "line", line = list(width = 1, color = "#A6A6A6"), layer = "below",
    x0 = range_x[1], x1 = range_x[2], y0 = 1, y1 = 1)
  line_horizontal_mid <- list(
    type = "line", line = list(width = 1, color = "#A6A6A6", dash = "dash"), layer = "below",
    x0 = range_x[1], x1 = range_x[2], y0 = 0.5, y1 = 0.5)
  line_horizontal_bot <- list(
    type = "line", line = list(width = 1, color = "black"), layer = "below",
    x0 = range_x[1], x1 = range_x[2], y0 = 0, y1 = 0)
  line_vertical <- list(
    type = "line", line = list(width = 1), layer = "below",
    x0 = 1, x1 = 1, y0 = range_y[1], y1 = range_y[2])
  lines <- list(line_horizontal_top, line_horizontal_mid, line_horizontal_bot, line_vertical)
  # add lines to *50 points
  if ("extras" %in% layers) {
    lines_extra <- extras$lines[grepl(var_y, names(extras$lines))]
    names(lines_extra) <- NULL
    lines <- c(lines, lines_extra)
  }
  
  # positions of axis ticks
  ticks_x <- 10 ^ (-10:10)
  ticks_y <- -2:4 / 2
  # axis options
  axis_options_x <- list(title = sprintf("%s [&mu;M]", concentration_name), type = "log",
                         showgrid = FALSE, zeroline = FALSE, range = log10(range_x),
                         tickvals = ticks_x, exponentformat = "e")
  axis_options_y <- list(title = var_y,
                         showgrid = FALSE, zeroline = FALSE, range = range_y, tickvals = ticks_y)
  
  # drop data points that would not fit in the plotting area
  # show exception_data if all data has been dropped
  data <- data[data.table::between(data[["var_y"]], range_y[1], range_y[2]), ]
  if (!NROW(data)) {
    # tmp solution - start
    exception_data <- gDRimport::get_exception_data(36)
    txt_msg <- sprintf(exception_data$sprintf_text, comb_name)
    txt_msg <- gsub("\\. |: ", "\n", txt_msg)
    txt_msg <- gsub(" or ", " \nor ", txt_msg)
    txt_msg <- paste(exception_data$title, "\n\n", txt_msg)
    # tmp solution - end
    plt_err <- plotly::layout(
      p = plotly::plot_ly(mode = "text", type = "scatter"),
      xaxis = list(visible = FALSE), yaxis = list(visible = FALSE),
      annotations = list(text = txt_msg, align = "justify",
                         showarrow = FALSE, font = list(size = 12, color = "darkred")))
    
    return(gDRcomponents::gDR_plotly_config(plt_err))
  }
  
  # format curve data
  if (is.element("curve", layers)) {
    data_curves <- data.table::as.data.table(curves)
    data.table::setnames(data_curves, old = c(var_y, var_col, var_not_col),
                         new = c("var_y", "var_col", "var_not_col"))
    # add column for tooltips
    data_curves$label <- buildLabel(data_curves, view = "curve")
  }
  
  # initialize plot
  plot_base <- plotly::plot_ly(hoverinfo = "text", showlegend = FALSE, colors = "Set1",
                               width = plot_width, height = plot_height)
  # add layers
  if ("curve" %in% layers) {
    plot_base <- plotly::add_lines(
      plot_base,
      x = data_curves[[concentration_name]], y = data_curves[["var_y"]],
      color = data_curves[["var_col"]], text = data_curves[["label"]],
      split = data_curves[["Treatment"]], line = list(width = 1))
  }
  if ("observations" %in% layers) {
    plot_base <- plotly::add_markers(
      plot_base,
      x = data[[concentration_name]], y = data[["var_y"]],
      color = data[["var_col"]], text = data[["label"]],
      size = I(5))
  }
  if ("average" %in% layers) {
    plot_base <- plotly::add_markers(
      plot_base,
      x = data[[concentration_name]], y = data[["var_y"]],
      color = data[["var_col"]], text = data[["label"]])
  }
  
  if ("error" %in% layers) {
    plot_base <- plotly::add_markers(
      plot_base,
      x = data[[concentration_name]], y = data[["var_y"]],
      color = data[["var_col"]], text = data[["label"]],
      error_y = list(array = data[["var_y.err"]],
                     thickness = 1, width = 3), marker = list(opacity = 0))
  }
  
  if ("extras" %in% layers) {
    data_extra <- extras$points[metric == var_y, ]
    # add points for half response
    plot_base <- plotly::add_markers(
      plot_base,
      x = data_extra[[concentration_name]], y = data_extra[["Response"]],
      size = I(5), color = I("black"))
    # add ticks for half response
    plot_base <- plotly::add_trace(plot_base,
                                   type = "scatter", mode = "text+markers",
                                   x = range_x[1] * 1.17, y = data_extra[["MaxEffectiveness"]],
                                   text = round(data_extra[["MaxEffectiveness"]], 2),
                                   textposition = "middle right", color = I("gray50"),
                                   marker = list(symbol = 25, size = 8))
    # add ticks for max effectiveness
    plot_base <- plotly::add_trace(plot_base,
                                   type = "scatter", mode = "text+markers",
                                   x = data_extra[[concentration_name]], y = range_y[1] + 0.033,
                                   text = signif(data_extra[[concentration_name]], 2),
                                   textposition = "top center", color = I("gray50"),
                                   marker = list(symbol = 26, size = 8))
  }
  
  # apply layout
  plot_laidout <- plotly::layout(plot_base,
                                 title = list(text = plot_title, font = list(size = 14)),
                                 xaxis = axis_options_x,
                                 yaxis = axis_options_y,
                                 shapes = lines)
  # modify config
  plot_final <- gDRcomponents::gDR_plotly_config(plot_laidout,
                                                 edits = gDRcomponents::get_plotly_edits(),
                                                 showAxisRangeEntryBoxes = FALSE)
  
  return(plot_final)
}

#' create a log-sequence
#'
#' Create a sequence of numbers growing in log-domain.
#'
#' The result is a numeric vector of length \code{length}.
#' Differences between items are constant in logarithmic domain
#' and therefore geometrically increase in linear domain.
#'
#' @param start,end numeric, lower and upper margins of the sequence
#' @param length integer, resulting sequence length
#' 
#' @keywords internal
#' @return A numeric vector, see \code{Details}.
#'
#' @seealso \code{\link{prepareCurves}}
#'
logSeq <- function(start, end, length) {
  
  checkmate::assert_number(start, lower = 0, finite = TRUE)
  checkmate::assert_number(end, lower = 0, finite = TRUE)
  checkmate::assert_number(length, lower = 1, finite = TRUE)
  
  limits <- c(start, end)
  limits_log <- log10(limits)
  sequence_log <- seq(from = limits_log[1], limits_log[2], length.out = length)
  sequence <- 10 ^ sequence_log
  return(sequence)
}

