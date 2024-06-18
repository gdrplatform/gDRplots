#' Draw interactive Response Grid panel
#'
#' Build a panel for the drug response grid.
#'
#' @param data a data.table of curve coordinates
#' @param var_y name of metric to plot: GR value or Relative Viability, given as character string
#' @param range_x numeric vector of length 2 that specifies the limits of the X axis
#' @param title character string defining the column, whose content will form the panel's title;
#'              if set to NA (the default), no title is displayed
#'
#' @return A \code{plotly} object.
#'
#' @keywords plugin_plot
#'
#' @seealso \code{\link[gDRplots]{prepareCurves}}
#'
#' @export
#'
plotlyRGPanel <- function(data, var_y, range_x = c(1e-3, 50e+0), title = NA) {
  
  checkmate::assert_data_table(data)
  checkmate::assert_string(var_y)
  checkmate::assert_choice(var_y, choices = c("GR value", "Relative Viability"))
  checkmate::assert_numeric(range_x, len = 2, finite = TRUE, any.missing = FALSE)
  checkmate::assert_string(title, na.ok = TRUE)
  
  pidfs <- gDRutils::get_prettified_identifiers(simplify = TRUE)
  # add a label column
  data$label <- build_label(data, "grid")
  # determine plot title
  plot_title <- unique(data[[title]])
  # rename column of interest
  data.table::setnames(data, old = var_y, new = "var_y")
  
  # set range of Y axis
  range_y <- switch(var_y,
                    "GR value" = c(-1, 1.1),
                    "Relative Viability" = c(0, 1.1))
  # line properties
  line_horizontal_top <- list(type = "line",
                              line = list(width = 1, color = "#A6A6A6"), layer = "below",
                              x0 = range_x[1], x1 = range_x[2], y0 = 1, y1 = 1)
  line_horizontal_mid <- list(type = "line",
                              line = list(width = 1, color = "#A6A6A6", dash = "dash"), layer = "below",
                              x0 = range_x[1], x1 = range_x[2], y0 = 0.5, y1 = 0.5)
  line_horizontal_bot <- list(type = "line",
                              line = list(width = 1, color = "#000000"), layer = "below",
                              x0 = range_x[1], x1 = range_x[2], y0 = 0, y1 = 0)
  line_vertical <- list(type = "line", line = list(width = 1), layer = "below",
                        x0 = 1, x1 = 1, y0 = range_y[1], y1 = range_y[2])
  
  # drop data points that would not fit in the plotting area
  data <- data[data.table::between(data[["var_y"]], range_y[1], range_y[2]), ]
  
  # build plots
  # 'Treatment' is not and identifier, it's a column defined in `gDRplots::prepareCurves`
  plot_base <- plotly::plot_ly(x = data[[pidfs[["concentration"]]]],
                               y = data[["var_y"]], 
                               split = data[["Treatment"]],
                               type = "scatter", mode = "lines",
                               text = data[["label"]], hoverinfo = "text",
                               showlegend = FALSE, color = I("grey50"), opacity = 0.25,
                               line = list(width = 1),
                               width = 240, height = 175,
                               source = "source")
  if (nrow(data) == 0) {
    plot_base <- plotly::plotly_empty()
  }
  
  # add layout
  plot_laidout <- plotly::layout(plot_base,
                                 title = list(text = plot_title,
                                              xref = "paper", x = 0.02, xanchor = "left",
                                              yref = "container", y = 0.98, yanchor = "top",
                                              font = list(size = 13)),
                                 shapes = list(line_horizontal_top, line_horizontal_mid, line_horizontal_bot,
                                               line_vertical),
                                 xaxis = list(title = list(
                                   text = sprintf("log10(%s)", pidfs[["concentration"]]),
                                   font = list(size = 11)
                                 ),
                                 range = log10(range_x), type = "log",
                                 zeroline = FALSE, showgrid = FALSE,
                                 fixedrange = TRUE, showticklabels = FALSE),
                                 yaxis = list(title = list(text = var_y, font = list(size = 11)),
                                              range = range_y,
                                              zeroline = FALSE, showgrid = FALSE,
                                              fixedrange = TRUE, showticklabels = FALSE),
                                 margin = list(l = 0, r = 0, t = 15, b = 0))
  # remove bar
  plot_final <- plotly::config(plot_laidout, displayModeBar = FALSE,
                               edits = get_plotly_edits(),
                               showAxisRangeEntryBoxes = FALSE)
  
  return(plot_final)
}
