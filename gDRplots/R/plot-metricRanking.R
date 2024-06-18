#' Draw interactive Metric Ranking plot
#'
#' Create a barplot that ranks drugs or cell lines across a growth metric.
#'
#' Special rules apply for GR50, IC50 and EC50.
#' These variables are plotted in log-scale.
#' The argument \code{axis_limits} determines whether the Y axis limits will be
#' fixed at 1e-4 to 100 (fixed), derived from data (free), or whether one
#' limit will be fixed and the other derived from the data (upper and lower).
#'
#' Note that in the app axes are editable and the user can alter the range by dragging the axis.
#'
#' @param data \code{data.table} prepared by \code{prepareDataMR}
#' @param var_x character string; name of variable to put on X axis
#' @param var_y character string; name of variable to put on Y axis
#' @param var_col character string; name of variable to color bars by
#' @param var_grp character string; name of variable to group values by
#' @param title name of variable to construct plot title from
#' @param show_tick_labels boolean; whether or not to show X axis labels
#' @param axis_limits character string; fix or release limits of Y axis, see \code{Details}
#'
#' @return a \code{plotly} object
#'
#' @keywords plugin_plot
#'
#' @export
#'
plotly_metric_ranking <- function(data, 
                                  var_x, 
                                  var_y, 
                                  var_col, 
                                  var_grp, 
                                  title, 
                                  show_tick_labels,
                                  axis_limits = c("fixed", "upper", "lower", "free")) {
  
  checkmate::assert_data_table(data)
  checkmate::assert_string(var_x)
  checkmate::assert_choice(var_x, choices = names(data))
  checkmate::assert_string(var_y)
  checkmate::assert_choice(var_y, choices = names(data))
  checkmate::assert_string(var_col)
  checkmate::assert_choice(var_col, choices = c("none", names(data)))
  checkmate::assert_string(var_grp)
  checkmate::assert_choice(var_grp, choices = c("none", names(data)))
  checkmate::assert_string(title)
  checkmate::assert_choice(title, choices = c("none", names(data)))
  axis_limits <- match.arg(axis_limits)
  
  if (nrow(data) == 0) {
    return(plotly::plotly_empty(type = "bar"))
  }
  
  # get prettified identifiers
  pidfs <- gDRutils::get_prettified_identifiers(simplify = TRUE)
  drug_name2 <- pidfs[["drug_name2"]]
  concentration2 <- pidfs[["concentration2"]]
  untreated_tags <- pidfs[["untreated_tag"]]
  
  # there are subsets with no co-drug data but still with drug_name_2 and concentration2 columns
  # let's remove columns with codrug data in such cases
  if (!gDRutils::has_valid_codrug_data(data)) {
    data <- gDRutils::remove_codrug_data(data)
  }
  
  # check for presence of the single codrug data
  has_codrug_data <- gDRutils::has_single_codrug_data(names(data))
  
  ### recursive incursion
  vars_cotreatment <- intersect(c(drug_name2, concentration2), names(data))
  
  if (has_codrug_data)  {
    co_drugs <- setdiff(unique(data[[drug_name2]]), untreated_tags)
    #  if co-drugs contains both one of `untreated_tags` and codrugs not being `untreated_tags`
    # replace former with the latter using `gDRplots::replace_values`
    # TODO: propose a more robust approach to handle more than two untreated tags
    if (length(co_drugs) > 1 && any(untreated_tags %in% data[[drug_name2]])) {
      data <-
        if (untreated_tags[1] %in% data[[drug_name2]]) {
          replace_values(data, drug_name2, untreated_tags[1], co_drugs)
        } else if (untreated_tags[2] %in% data[[drug_name2]]) {
          replace_values(data, drug_name2, untreated_tags[2], co_drugs)
        } else {
          stop(sprintf("unsupported value for 'untreated_tags:%s", untreated_tags))
        }
      # split data on co-drugs
      data_split <- split(data, data[[drug_name2]])
      # plot each subset
      plot_list <- lapply(
        data_split,
        plotlyMR,
        var_x = var_x, var_y = var_y, var_col = var_col, var_grp = var_grp,
        title = title,
        show_tick_labels = show_tick_labels
      )
      # determine final plot dimensions
      # collate into one subplot
      plot_sub <- plotly::subplot(plot_list, nrows = length(plot_list),
                                  margin = 0.05, titleX = TRUE)
      # modify config
      plot_final <- gDR_plotly_config(plot_sub,
                                      edits = get_plotly_edits(),
                                      showAxisRangeEntryBoxes = FALSE)
      
      return(plot_final)
    }
  }
  ### end recursive incursion
  
  # axis options and titles
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
    range_y <- switch(var_y,
                      "E max" = c(0, 1.1),
                      "GR max" = c(-1, 1.1),
                      c(-0.1, NA))
    ticks_y <- switch(var_y,
                      "E max" = -10:10 / 4,
                      "GR max" = -10:10 / 2,
                      -10:10 / 5)
    horizontal <- 0
  }
  titlePlot <- sprintf("Ranking %s for %s: %s", title_x, title, unique(data[[title]]))
  line_horizontal <- list(type = "line",
                          line = list(width = 1, color = "#A6A6A6"), layer = "below",
                          x0 = 0, x1 = 1, xref = "paper", y0 = horizontal, y1 = horizontal)
  
  # determine font size for X axis labels as plotly omits some of them if they don"t fit
  font_size <- floor(225 / max(nchar(unique(data[[var_x]]))))
  font_size <- if (font_size > 15) {
    15
  } else if (font_size < 9) {
    9
  }
  # determine bottom margin based on the labels" length
  longest_label <- max(nchar(unique(data[[var_x]])))
  margin_bottom <- longest_label * font_size / 2
  
  # order by grouping variable and ascending metric values
  if (var_grp == "none") {
    data.table::setorderv(data, var_y)
  } else {
    data.table::setorderv(data, c(var_grp, var_y))
  }
  # convert to factor to maintain bar order
  data[[var_x]] <- factor(data[[var_x]], levels = unique(data[[var_x]]))
  
  # add label column
  data$label <- build_label(data, "ranking")
  # modify column Concentration_2 for proper sorting of bars
  if (is.element(concentration2, names(data))) {
    data[[concentration2]] <- as.character(data[[concentration2]])
    data[[concentration2]] <- sub("(^\\d\\..*)", "0\\1", data[[concentration2]])
    data[[concentration2]] <- as.factor(data[[concentration2]])
    levels(data[[concentration2]]) <- sub("^(0+)(\\d.*)", "\\2", levels(data[[concentration2]]))
  } else {
    data[[concentration2]] <- NULL
  }
  
  # should the legend be shown ?
  do_show_legend <- do_show_legend(var_col, data)
  
  # base plot
  plot_base <- plotly::plot_ly(
    x = data[[var_x]], y = data[[var_y]],
    color = data[[var_col]],
    split = data[[concentration2]],
    text = data[["label"]], hoverinfo = "text", textposition = "none",
    type = "bar",
    showlegend = do_show_legend)
  # add layout
  plot_laidout <- plotly::layout(
    plot_base,
    title = list(text = titlePlot),
    xaxis = list(title = var_x, tickangle = -90, tickfont = list(size = font_size),
                 # condition showing tick labels on whether all labels will fit
                 showticklabels = show_tick_labels),
    yaxis = list(title = title_x, showgrid = FALSE,
                 type = type_y, range = range_y, tickvals = ticks_y,
                 exponentformat = "e", zeroline = FALSE),
    shapes = line_horizontal,
    margin = list(b = margin_bottom),
    legend = list(orientation = "v",
                  title = get_legend_title(var_col, has_codrug_data = has_codrug_data))
  )
  # modify config
  plot_final <- gDR_plotly_config(plot_laidout,
                                  edits = get_plotly_edits(),
                                  showAxisRangeEntryBoxes = FALSE)
  
  return(plot_final)
}
