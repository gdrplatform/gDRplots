#' Plot box plots for metric for single-agent data grouped by selected variable
#' 
#' @param dt_metrics data.table representing data from the \code{Metrics} assay,
#'    outputted by \code{gDRutils::convert_se_assay_to_dt(se, "Metrics")}
#'    and single-agent \code{SummarizedExperiment}
#' @param group_var string name of group variable; one of: \code{"CellLineName"} or \code{"DrugName"}
#' @param normalization_type string with normalization types to be selected
#'                           one of: "GR" ("GRvalue") or "RV" ("RelativeViability")
#' @param metric string name of the metric;
#'    one of: "xc50" ("GR50" or "IC50" - respectively depending on \code{normalization_type}), 
#'    "x_max" ("GR Max" or "E Max") or "x_mean" ("GR Mean" or "RV Mean");
#'    but the values from any numeric column can be displayed.
#' @param fit_source string source name for metrics
#' @param grouped_flag a logical flag whether the boxplots should be grouped and 
#'    colored by \code{Tissue} for \code{group_var} set as \code{"CellLineName"} 
#'    and \code{drug_moa} - for \code{"DrugName"}
#' @param colored_pts_flag a logical flag whether the points should be colored by grouped variable -
#'    for \code{group_var} equal \code{"CellLineName"} points will be colored by \code{"DrugName"}
#'    and similarly vice versa
#' @param colors_vec character vector with colors (name or hex value) to color boxplots; 
#'    for \code{grouped_flag} set as FALSE only first from vector will be used.
#' @param with_inf a logical flag indicating whether infinite values should be shown on boxplots
#' 
#' @return \code{ggplot} object containing boxplots for selected single-agent metric 
#'    grouped by \code{group_var}
#' 
#' @keywords single-agent_plots
#' 
#' @author Janina Smoła \email{janina.smola@@contractors.roche.com}
#' 
#' @examples
#' mae <- gDRutils::get_synthetic_data("combo_matrix")
#' se <- mae[[gDRutils::get_supported_experiments("sa")]]
#' 
#' dt_metrics <- gDRutils::convert_se_assay_to_dt(se, "Metrics")
#' 
#' plot_boxplot_metric_sa(dt_metrics,
#'                        group_var = "CellLineName")
#' 
#' plot_boxplot_metric_sa(dt_metrics,
#'                        group_var = "DrugName")
#' 
#' plot_boxplot_metric_sa(dt_metrics,
#'                        group_var = "DrugName",
#'                        metric = "x_AOC_range",
#'                        colors_vec = "grey",
#'                        colored_pts_flag = TRUE)
#' 
#' plot_boxplot_metric_sa(dt_metrics,
#'                        group_var = "CellLineName",
#'                        normalization_type = "RV",
#'                        metric = "x_max",
#'                        colors_vec = c("gold", "darkorange", "darkcyan", "darkblue"),
#'                        grouped_flag = TRUE)
#' 
#' @export
plot_boxplot_metric_sa <- function(
    dt_metrics,
    group_var,
    normalization_type = "GR",
    metric = "xc50",
    fit_source = "gDR",
    grouped_flag = FALSE,
    colored_pts_flag = FALSE,
    colors_vec = NULL,
    with_inf = FALSE
) {
  
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  tissue <- gDRutils::get_env_identifiers("cellline_tissue")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  drug_MOA <- gDRutils::get_env_identifiers("drug_moa")
  
  checkmate::assert_data_table(dt_metrics)
  checkmate::assert_choice(group_var, choices = c(cellline_name, drug_name))
  
  # check input data
  if (group_var == cellline_name) {
    point_var <- drug_name
    col_var <- tissue
  } else if (group_var == drug_name) {
    point_var <- cellline_name
    col_var <- drug_MOA
  }
  
  checkmate::assert_choice(normalization_type, choices = c("GR", "RV"))
  numeric_columns <- names(dt_metrics)[vapply(dt_metrics, is.numeric, logical(1))]
  checkmate::assert_choice(metric, choices = numeric_columns)
  checkmate::assert_names(names(dt_metrics), 
                          must.include = c(cellline_name, drug_name, col_var, metric))
  checkmate::assert_string(fit_source, null.ok = TRUE)
  checkmate::assert_flag(grouped_flag)
  checkmate::assert_flag(colored_pts_flag)
  checkmate::assert_character(colors_vec, null.ok = TRUE)
  checkmate::assert_flag(with_inf)
  boxplot_fill <- 
    gDRutils::get_settings_from_json("BOXPLOT_FILL",
                                     system.file(package = "gDRplots", "settings.json"))
  hline_color <- 
    gDRutils::get_settings_from_json("HLINE_COLOR",
                                     system.file(package = "gDRplots", "settings.json"))
  jitter_poinst_color <- 
    gDRutils::get_settings_from_json("JITTER_POINST_COLOR",
                                     system.file(package = "gDRplots", "settings.json"))
  edge_color <- 
    gDRutils::get_settings_from_json("EDGE_COLOR",
                                     system.file(package = "gDRplots", "settings.json"))
  
  if (grouped_flag && colored_pts_flag) {
    message("Please, choose only one coloring option: either `grouped_flag` or `colored_pts_flag`.")
  }
  
  # filter data for normalization type
  filter_expr <- substitute(normalization_type == norm_type & fit_source == fit_src,
                            list(norm_type = normalization_type, fit_src = fit_source))
  dt_met <- dt_metrics[eval(filter_expr)]
  dt_met <- dt_met[, c(group_var, point_var, col_var, metric), with = FALSE]
  
  if (metric == "xc50") {
    dt_met[, (metric) := log10(get(metric))] 
  }
  
  # handle -Inf (NA will be not shown on boxplots when with_inf = FALSE)
  if (!with_inf) {
    dt_met[is.infinite(get(metric)), (metric) := NA] 
  }
  
  group_label <- if (point_var == cellline_name) {
    "celllines"
  } else {
    "drugs"
  }
  plt_title <- sprintf("Number of unique %s: %s", group_label, NROW(unique(dt_met[[point_var]])))
  
  if (grouped_flag) {
    data.table::setorderv(dt_met, col_var)
    dt_met[[group_var]] <- factor(dt_met[[group_var]], levels = unique(dt_met[[group_var]]))
    
    fill_colors <- if (is.null(colors_vec) || !all(vapply(colors_vec, is_valid_color, logical(1)))) {
      get_qual_colors(NROW(unique(dt_met[[col_var]])))
    } else if (NROW(colors_vec) != NROW(unique(dt_met[[col_var]]))) {
      grDevices::colorRampPalette(colors_vec)(NROW(unique(dt_met[[col_var]])))
    } else {
      colors_vec
    }
    names(fill_colors) <- unique(dt_met[[col_var]])
    
    plt <- 
      ggplot2::ggplot(data = dt_met,
                      mapping = ggplot2::aes(x = get(group_var), y = get(metric))) +
      ggplot2::geom_hline(yintercept = 0, color = hline_color, linetype = "solid") +
      ggplot2::geom_point(ggplot2::aes(fill = get(col_var)), size = -1, alpha = 0.25, na.rm = TRUE) +
      ggplot2::geom_boxplot(ggplot2::aes(fill = get(col_var)), 
                            color = edge_color, alpha = 0.25, staplewidth = 0.5,
                            na.rm = TRUE, outliers = FALSE, show.legend = FALSE) +
      ggplot2::scale_fill_manual(name = col_var, values = fill_colors) +
      ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(shape = 22, size = 10)))
    
  } else {
    data.table::setorderv(dt_met, group_var)
    dt_met[[group_var]] <- factor(dt_met[[group_var]])
    
    fill_color <- if (is.null(colors_vec) || !all(vapply(colors_vec, is_valid_color, logical(1)))) {
      boxplot_fill
    } else {
      colors_vec[1]
    }
    
    plt <- 
      ggplot2::ggplot(data = dt_met,
                      mapping = ggplot2::aes(x = get(group_var), y = get(metric))) +
      ggplot2::geom_hline(yintercept = 0, color = hline_color, linetype = "solid") +
      ggplot2::geom_boxplot(fill = fill_color, 
                            color = edge_color, alpha = 0.25, staplewidth = 0.5,
                            na.rm = TRUE, outliers = FALSE) +
      ggplot2::theme(legend.position = "none")
  }
  
  # validation for coloring points
  if (colored_pts_flag) {
    if (NROW(unique(dt_met[[point_var]])) >= 10) colored_pts_flag <- FALSE
  }
  # jitter point
  if (colored_pts_flag) {
    color_points <- get_qual_colors(NROW(unique(dt_met[[point_var]])))
    names(color_points) <- unique(dt_met[[point_var]])
    
    plt <- plt +
      ggplot2::geom_jitter(mapping = ggplot2::aes(color = get(point_var)), size = 2, alpha = 0.75,
                           width = 0.2, height = 0, na.rm = TRUE) +
      ggplot2::scale_color_manual(name = drug_name, values = color_points) +
      ggplot2::guides(color = ggplot2::guide_legend(title = point_var))
  } else {
    plt <- plt +
      ggplot2::geom_jitter(color = jitter_poinst_color, alpha = 0.75,
                           width = 0.2, height = 0, na.rm = TRUE)
  }
  
  # final
  plt <- plt +
    ggplot2::labs(title = plt_title,
                  y = get_hm_title(metric, normalization_type), 
                  x = "") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8, angle = 90, vjust = 1, hjust = 1),
                   axis.text.y = ggplot2::element_text(size = 8),
                   plot.title = ggplot2::element_text(size = 10),
                   panel.grid.minor = ggplot2::element_blank())
  
  return(plt)
}

#' Plot box plots for metric for single-agent data grouped by cell line names
#' 
#' @inheritParams plot_boxplot_metric_sa
#' 
#' @return \code{ggplot} object containing boxplots for selected single-agent metric 
#'    grouped by cellline names
#' 
#' @keywords single-agent_plots
#' 
#' @author Janina Smoła \email{janina.smola@@contractors.roche.com}
#' 
#' @examples
#' mae <- gDRutils::get_synthetic_data("combo_matrix")
#' se <- mae[[gDRutils::get_supported_experiments("sa")]]
#' 
#' dt_metrics <- gDRutils::convert_se_assay_to_dt(se, "Metrics")
#' 
#' plot_boxplot_metric_sa_by_CLs(dt_metrics)
#' 
#' plot_boxplot_metric_sa_by_CLs(dt_metrics,
#'                               grouped_flag = TRUE)
#'                               
#' plot_boxplot_metric_sa_by_CLs(dt_metrics,
#'                               with_inf = TRUE)
#'                               
#' plot_boxplot_metric_sa_by_CLs(dt_metrics,
#'                               metric = "x_AOC",
#'                               colors_vec = "gold")
#' 
#' plot_boxplot_metric_sa_by_CLs(
#'   dt_metrics,
#'   metric = "x_max",
#'   colored_pts_flag = TRUE)
#' 
#' plot_boxplot_metric_sa_by_CLs(
#'   dt_metrics,
#'   metric = "x_mean",
#'   grouped_flag = TRUE,
#'   colors_vec = c("deeppink", "darkcyan", "orange", "darkblue"))
#' 
#' @export
plot_boxplot_metric_sa_by_CLs <- function(
    dt_metrics,
    normalization_type = "GR",
    metric = "xc50",
    fit_source = "gDR",
    grouped_flag = FALSE,
    colored_pts_flag = FALSE,
    colors_vec = NULL,
    with_inf = FALSE
) {
  
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  tissue <- gDRutils::get_env_identifiers("cellline_tissue")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  
  checkmate::assert_data_table(dt_metrics)
  checkmate::assert_choice(normalization_type, choices = c("GR", "RV"))
  numeric_columns <- names(dt_metrics)[vapply(dt_metrics, is.numeric, logical(1))]
  checkmate::assert_choice(metric, choices = numeric_columns)
  checkmate::assert_names(names(dt_metrics), 
                          must.include = c(cellline_name, tissue, drug_name, metric))
  checkmate::assert_string(fit_source, null.ok = TRUE)
  checkmate::assert_flag(grouped_flag)
  checkmate::assert_flag(colored_pts_flag)
  checkmate::assert_character(colors_vec, null.ok = TRUE)
  checkmate::assert_flag(with_inf)
  
  plt <- 
    plot_boxplot_metric_sa(
      dt_metrics = dt_metrics,
      group_var = cellline_name,
      normalization_type = normalization_type,
      metric = metric,
      fit_source = fit_source,
      grouped_flag = grouped_flag,
      colored_pts_flag = colored_pts_flag,
      colors_vec = colors_vec,
      with_inf = with_inf
    )
  
  return(plt)
}

#' Plot box plots for metric for single-agent data grouped by drug names
#' 
#' @inheritParams plot_boxplot_metric_sa
#' 
#' @return \code{ggplot} object containing boxplots for selected single-agent metric
#'    grouped by drug names
#' 
#' @keywords single-agent_plots
#' 
#' @author Janina Smoła \email{janina.smola@@contractors.roche.com}
#' 
#' @examples
#' mae <- gDRutils::get_synthetic_data("combo_matrix")
#' se <- mae[[gDRutils::get_supported_experiments("sa")]]
#' 
#' dt_metrics <- gDRutils::convert_se_assay_to_dt(se, "Metrics")
#' 
#' plot_boxplot_metric_sa_by_drugs(dt_metrics,
#'                                 normalization_type = "RV")
#' 
#' plot_boxplot_metric_sa_by_drugs(dt_metrics,
#'                                 normalization_type = "RV",
#'                                 grouped_flag = TRUE)
#'                               
#' plot_boxplot_metric_sa_by_drugs(dt_metrics,
#'                                 metric = "x_AOC",
#'                                 colors_vec = "gold")
#' 
#' plot_boxplot_metric_sa_by_drugs(
#'   dt_metrics,
#'   metric = "x_max",
#'   colored_pts_flag = TRUE)
#' 
#' plot_boxplot_metric_sa_by_drugs(
#'   dt_metrics,
#'   metric = "x_mean",
#'   grouped_flag = TRUE,
#'   colors_vec = c("deeppink", "darkcyan", "orange", "darkblue", "limegreen"))
#' 
#' @export
plot_boxplot_metric_sa_by_drugs <- function(
    dt_metrics,
    normalization_type = "GR",
    metric = "xc50",
    fit_source = "gDR",
    grouped_flag = FALSE,
    colored_pts_flag = FALSE,
    colors_vec = NULL,
    with_inf = FALSE
) {
  
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  drug_MOA <- gDRutils::get_env_identifiers("drug_moa")
  
  checkmate::assert_data_table(dt_metrics)
  checkmate::assert_choice(normalization_type, choices = c("GR", "RV"))
  numeric_columns <- names(dt_metrics)[vapply(dt_metrics, is.numeric, logical(1))]
  checkmate::assert_choice(metric, choices = numeric_columns)
  checkmate::assert_names(names(dt_metrics), 
                          must.include = c(cellline_name, drug_name, drug_MOA, metric))
  checkmate::assert_string(fit_source, null.ok = TRUE)
  checkmate::assert_flag(grouped_flag)
  checkmate::assert_flag(colored_pts_flag)
  checkmate::assert_character(colors_vec, null.ok = TRUE)
  checkmate::assert_flag(with_inf)
  
  plt <- 
    plot_boxplot_metric_sa(
      dt_metrics = dt_metrics,
      group_var = drug_name,
      normalization_type = normalization_type,
      metric = metric,
      fit_source = fit_source,
      grouped_flag = grouped_flag,
      colored_pts_flag = colored_pts_flag,
      colors_vec = colors_vec,
      with_inf = with_inf
    )
  
  return(plt)
}

#' Plot box plots for metric for single-agent data grouped by selected variable
#' 
#' @inheritParams plot_boxplot_metric_sa
#' 
#' @param selection_var string name of selected main variable - one value from column 
#'    \code{"CellLineName"} or \code{"DrugName"}
#' @param selection_name string name of selected variable value from column \code{selection_var} 
#'    to filter data for plotting;
#' @param group_var string name of group variable; not numeric variablee from \code{dt_metrics} 
#'    different than \code{selection_var} and not containing unique values for each row;
#' @param group_names character vector with names to subset from column \code{group_var};
#'    if \code{NULL} then all values will be plotted
#' @param named_n_bottom number of n-bottom \code{metric} values to be labeled on the plot
#'    for \code{group_var} equal \code{"DrugName"} points will be labeled by \code{"CellLineName"}
#'    and similarly vice versa
#' @param grouped_flag logical flag whether the boxplots should be colored by \code{group_var}
#'   
#' @return \code{ggplot} object containing boxplots for selected single-agent metric 
#'    grouped by selected variable
#' 
#' @keywords single-agent_plots
#' 
#' @author Janina Smoła \email{janina.smola@@contractors.roche.com}
#' 
#' @examples
#' mae <- gDRutils::get_synthetic_data("combo_matrix")
#' se <- mae[[gDRutils::get_supported_experiments("sa")]]
#' 
#' dt_metrics <- gDRutils::convert_se_assay_to_dt(se, "Metrics")
#' invisible(dt_metrics[, Tissue_grp := data.table::fifelse(Tissue == "tissue_w", 
#'                                                          "tissue_w", 
#'                                                          "tissue_other")])
#' 
#' plot_boxplot_metric_sa_by_grp(dt_metrics,
#'                               selection_var = "DrugName",
#'                               selection_name = "drug_001",
#'                               group_var = "Tissue")
#' 
#' plot_boxplot_metric_sa_by_grp(dt_metrics,
#'                               selection_var = "DrugName",
#'                               selection_name = "drug_001",
#'                               group_var = "Tissue_grp")
#' 
#' plot_boxplot_metric_sa_by_grp(dt_metrics,
#'                               selection_var = "DrugName",
#'                               selection_name = "drug_001",
#'                               group_var = "Tissue",
#'                               with_inf = TRUE)
#' 
#' plot_boxplot_metric_sa_by_grp(dt_metrics,
#'                               selection_var = "DrugName",
#'                               selection_name = "drug_001",
#'                               group_var = "Tissue",
#'                               named_n_bottom = 0)
#' 
#' plot_boxplot_metric_sa_by_grp(dt_metrics,
#'                               selection_var = "DrugName",
#'                               selection_name = "drug_001",
#'                               group_var = "Tissue",
#'                               grouped_flag = TRUE)
#'                               
#' plot_boxplot_metric_sa_by_grp(dt_metrics,
#'                               selection_var = "DrugName",
#'                               selection_name = "drug_001",
#'                               group_var = "Tissue",
#'                               colors_vec = c("darkblue", "deeppink"))
#' 
#' @export
plot_boxplot_metric_sa_by_grp <- function(
    dt_metrics,
    selection_var,
    selection_name,
    group_var,
    group_names = NULL,
    normalization_type = "GR",
    metric = "xc50",
    fit_source = "gDR",
    named_n_bottom = 5,
    grouped_flag = FALSE,
    colors_vec = NULL,
    with_inf = FALSE
) {
  
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  numeric_columns <- names(dt_metrics)[vapply(dt_metrics, is.numeric, logical(1))]
  
  checkmate::assert_data_table(dt_metrics)
  checkmate::assert_choice(selection_var, choices = c(cellline_name, drug_name))
  
  # point data
  if (selection_var == cellline_name) {
    point_var <- drug_name
  } else if (selection_var == drug_name) {
    point_var <- cellline_name
  }
  
  checkmate::assert_string(selection_name)
  checkmate::assert_choice(selection_name, choices = unique(dt_metrics[[selection_var]]))
  checkmate::assert_string(group_var)
  checkmate::assert_choice(
    group_var, 
    choices = names(dt_metrics)[!names(dt_metrics) %in% c(numeric_columns, cellline_name, drug_name)])
  # TODO add validation for number of levels >1 and !=NROW(dt_metrics)
  checkmate::assert_choice(group_names, choices = unique(dt_metrics[[group_var]]), null.ok = TRUE)
  checkmate::assert_choice(normalization_type, choices = c("GR", "RV"))
  checkmate::assert_choice(metric, choices = numeric_columns)
  checkmate::assert_names(names(dt_metrics), 
                          must.include = c(cellline_name, drug_name, group_var, metric))
  checkmate::assert_string(fit_source, null.ok = TRUE)
  checkmate::assert_number(named_n_bottom, lower = 0)
  checkmate::assert_flag(grouped_flag)
  checkmate::assert_character(colors_vec, null.ok = TRUE)
  checkmate::assert_flag(with_inf)
  boxplot_fill <- 
    gDRutils::get_settings_from_json("BOXPLOT_FILL",
                                     system.file(package = "gDRplots", "settings.json"))
  hline_color <- 
    gDRutils::get_settings_from_json("HLINE_COLOR",
                                     system.file(package = "gDRplots", "settings.json"))
  jitter_poinst_color <- 
    gDRutils::get_settings_from_json("JITTER_POINST_COLOR",
                                     system.file(package = "gDRplots", "settings.json"))
  edge_color <- 
    gDRutils::get_settings_from_json("EDGE_COLOR",
                                     system.file(package = "gDRplots", "settings.json"))
  
  # filter data for normalization type
  filter_expr <- substitute(normalization_type == norm_type & fit_source == fit_src,
                            list(norm_type = normalization_type, fit_src = fit_source))
  dt_met <- dt_metrics[eval(filter_expr)]
  # filter by selection
  dt_met <- dt_met[get(selection_var) == selection_name, ]
  # select min required data for plotting
  dt_met <- dt_met[, c(group_var, point_var, metric), with = FALSE]
 
  
  # coloring points by rank
  data.table::setorderv(dt_met, cols = metric)
  dt_met[, `:=`(is_bottom = FALSE, label = "")]
  if (named_n_bottom > 0) {
    named_n_bottom <- min(named_n_bottom, NROW(dt_met)) # deal with less than n bigger than table
    dt_met <- 
      dt_met[order(get(metric)), ][seq_len(named_n_bottom), `:=`(is_bottom = TRUE, label = get(point_var))]
  }
  
  if (metric == "xc50") {
    dt_met[, (metric) := log10(get(metric))] 
  }
  
  # handle -Inf (NA will be not shown on boxplots when with_inf = FALSE)
  if (!with_inf) {
    dt_met[is.infinite(get(metric)), (metric) := NA] 
  }
  
  # update group (it depends on user choice for `group_names` and `selection_name`)
  group_names <- if (is.null(group_names)) {
    sort(unique(dt_met[[group_var]]))
  } else {
    sort(intersect(unique(dt_met[[group_var]]), group_names))
  }
  
  plt_title <- sprintf("%s for %s by %s", 
                       gDRplots::get_hm_title(metric, normalization_type), selection_name, group_var)
  
  dt_met_lbl <- data.table::copy(dt_met)[!is.na(get(metric)), ]
  data.table::setorderv(dt_met, group_var)
  dt_met[[group_var]] <- factor(dt_met[[group_var]], levels = group_names)
  
  if (grouped_flag) {
    fill_colors <- if (is.null(colors_vec) || !all(vapply(colors_vec, is_valid_color, logical(1)))) {
      get_qual_colors(NROW(unique(dt_met[[group_var]])))
    } else if (NROW(colors_vec) != NROW(unique(dt_met[[group_var]]))) {
      grDevices::colorRampPalette(colors_vec)(NROW(unique(dt_met[[group_var]])))
    } else {
      colors_vec
    }
    names(fill_colors) <- unique(dt_met[[group_var]])
    
    plt <- 
      ggplot2::ggplot(data = dt_met,
                      mapping = ggplot2::aes(x = get(group_var), 
                                             y = get(metric))) +
      ggplot2::geom_hline(yintercept = 0, color = hline_color, linetype = "solid") +
      ggplot2::geom_boxplot(ggplot2::aes(fill = get(group_var)), 
                            color = edge_color, alpha = 0.25, staplewidth = 0.5,
                            na.rm = TRUE, outliers = FALSE, show.legend = FALSE) +
      ggplot2::scale_fill_manual(name = group_var, values = fill_colors) +
      ggplot2::guides(fill = "none")
    
  } else {
    fill_color <- if (is.null(colors_vec) || !all(vapply(colors_vec, is_valid_color, logical(1)))) {
      boxplot_fill
    } else {
      colors_vec[1]
    }
    
    plt <- 
      ggplot2::ggplot(data = dt_met,
                      mapping = ggplot2::aes(x = get(group_var), 
                                             y = get(metric))) +
      ggplot2::geom_hline(yintercept = 0, color = hline_color, linetype = "solid") +
      ggplot2::geom_boxplot(fill = fill_color, 
                            color = edge_color, alpha = 0.25, staplewidth = 0.5,
                            na.rm = TRUE, outliers = FALSE)
  }
  
  # adding lbls and points colored when bottom
  if (named_n_bottom > 0) {
    plt <- plt +
      ggrepel::geom_text_repel(data = dt_met_lbl,
                               mapping = ggplot2::aes(x = get(group_var), 
                                                      y = get(metric),
                                                      label = label),
                               size = 3, max.overlaps = 20, show.legend = FALSE) +
      ggplot2::geom_hline(yintercept = max(dt_met[is_bottom == TRUE][[metric]]), 
                          color = "red", linetype = "dashed")
  }
  plt <- plt +
    ggplot2::geom_jitter(mapping = ggplot2::aes(color = is_bottom), 
                         size = 2, alpha = 0.75,
                         width = 0.2, height = 0, na.rm = TRUE, 
                         show.legend = (named_n_bottom > 0)) +
    ggplot2::scale_color_manual(values = c("TRUE" = "red", "FALSE" = jitter_poinst_color))
  
  # final
  plt <- plt +
    ggplot2::labs(title = plt_title,
                  y = get_hm_title(metric, normalization_type), 
                  x = "",
                  color = sprintf("Bottom %s", named_n_bottom)) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8, angle = 90, vjust = 1, hjust = 1),
                   axis.text.y = ggplot2::element_text(size = 8),
                   plot.title = ggplot2::element_text(size = 10),
                   panel.grid.minor = ggplot2::element_blank())
  
  return(plt)
} 


#' Plot box plots for metric for combo data grouped by selected variable
#' 
#' @inheritParams plot_boxplot_metric_sa
#' 
#' @param dt_scores \code{data.table} representing data from the \code{scores} assay,
#'   outputted by \code{gDRutils::convert_se_assay_to_dt(se, "scores")}
#'   and combo \code{SummarizedExperiment}
#' @param group_var string name of group variable; one of: \code{"CellLineName"} or \code{"DrugName"};
#'   for \code{group_var} set as  \code{"DrugName"} points will be grouped by drug 
#'   combinations \code{"DrugName"} x \code{"DrugName_2"}
#' @param metric string name of the combo metric;
#'   one of: "hsa_score"("Bliss Excess GR" or "Bliss Excess RV" - respectively 
#'   depending on \code{normalization_type}), "bliss_score" ("Bliss Score GR" or "Bliss Score RV")
#' 
#' @return \code{ggplot} object containing boxplots for selected combo metric 
#'   grouped by \code{group_var}
#' 
#' @keywords combo_plots
#' 
#' @author Janina Smoła \email{janina.smola@@contractors.roche.com}
#' 
#' @examples
#' mae <- gDRutils::get_synthetic_data("combo_matrix")
#' se <- mae[[gDRutils::get_supported_experiments("combo")]]
#' 
#' dt_scores <- gDRutils::convert_se_assay_to_dt(se = se,
#'                                               assay_name = "scores")
#' 
#' plot_boxplot_metric_combo(dt_scores,
#'                           group_var = "DrugName")
#' 
#' plot_boxplot_metric_combo(dt_scores,
#'                           group_var = "DrugName", 
#'                           colored_pts_flag = TRUE,
#'                           colors_vec = "grey")
#' 
#' plot_boxplot_metric_combo(dt_scores,
#'                           group_var = "CellLineName",
#'                           metric = "hsa_score",
#'                           normalization_type = "RV",
#'                           grouped_flag = TRUE)
#' 
#' plot_boxplot_metric_combo(
#'   dt_scores,
#'   group_var = "CellLineName",
#'   metric = "hsa_score",
#'   grouped_flag = TRUE,
#'   colors_vec = c("deeppink", "darkblue"))
#' 
#' @export
plot_boxplot_metric_combo <- function(
    dt_scores,
    group_var,
    normalization_type = "GR",
    metric = "hsa_score",
    fit_source = "gDR",
    grouped_flag = FALSE,
    colored_pts_flag = FALSE,
    colors_vec = NULL
) {
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  tissue <- gDRutils::get_env_identifiers("cellline_tissue")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  drug_name_2 <- gDRutils::get_env_identifiers("drug_name2")
  
  checkmate::assert_data_table(dt_scores)
  checkmate::assert_choice(group_var, choices = c(cellline_name, drug_name))
  
  # check input data
  if (group_var == cellline_name) {
    point_var <- "DrugCombination"
    col_var <- tissue
  } else if (group_var == drug_name) {
    point_var <- cellline_name
    group_var <- "DrugCombination"
    if (grouped_flag) message("Coloring box by group is not available for this scenario.")
    grouped_flag <- FALSE
    col_var <- NULL
  }
  
  checkmate::assert_choice(normalization_type, choices = c("GR", "RV"))
  checkmate::assert_choice(metric, choices = c("hsa_score", "bliss_score"))
  checkmate::assert_names(names(dt_scores), 
                          must.include = c(cellline_name, drug_name, drug_name_2, col_var, metric))
  checkmate::assert_string(fit_source, null.ok = TRUE)
  checkmate::assert_flag(grouped_flag)
  checkmate::assert_flag(colored_pts_flag)
  checkmate::assert_character(colors_vec, null.ok = TRUE)
  boxplot_fill <- 
    gDRutils::get_settings_from_json("BOXPLOT_FILL",
                                     system.file(package = "gDRplots", "settings.json"))
  hline_color <- 
    gDRutils::get_settings_from_json("HLINE_COLOR",
                                     system.file(package = "gDRplots", "settings.json"))
  jitter_poinst_color <- 
    gDRutils::get_settings_from_json("JITTER_POINST_COLOR",
                                     system.file(package = "gDRplots", "settings.json"))
  edge_color <- 
    gDRutils::get_settings_from_json("EDGE_COLOR",
                                     system.file(package = "gDRplots", "settings.json"))
  
  if (grouped_flag && colored_pts_flag) {
    message("Please, choose only one coloring option: either `grouped_flag` or `colored_pts_flag`.")
  }
  
  # filter data for normalization type
  filter_expr <- substitute(normalization_type == norm_type & fit_source == fit_src,
                            list(norm_type = normalization_type, fit_src = fit_source))
  dt_sco_norm <- dt_scores[eval(filter_expr)]
  
  dt_sco <- dt_sco_norm[, c(cellline_name, drug_name, drug_name_2, col_var, metric), with = FALSE]
  dt_sco$DrugCombination <-
    paste(dt_sco[[drug_name]], "x", dt_sco[[drug_name_2]])
  
  group_label <- if (point_var == cellline_name) {
    "celllines"
  } else {
    "drug combinations"
  }
  plt_title <- sprintf("Number of unique %s: %s", group_label, NROW(unique(dt_sco[[point_var]])))
  
  if (grouped_flag) {
    data.table::setorderv(dt_sco, col_var)
    dt_sco[[group_var]] <- factor(dt_sco[[group_var]], levels = unique(dt_sco[[group_var]]))
    
    fill_colors <- if (is.null(colors_vec) || !all(vapply(colors_vec, is_valid_color, logical(1)))) {
      get_qual_colors(NROW(unique(dt_sco[[col_var]])))
    } else if (NROW(colors_vec) != NROW(unique(dt_sco[[col_var]]))) {
      grDevices::colorRampPalette(colors_vec)(NROW(unique(dt_sco[[col_var]])))
    } else {
      colors_vec
    }
    names(fill_colors) <- unique(dt_sco[[col_var]])
    
    plt <-
      ggplot2::ggplot(data = dt_sco,
                      mapping = ggplot2::aes(x = get(group_var), y = get(metric))) +
      ggplot2::geom_hline(yintercept = 0, color = hline_color, linetype = "solid") +
      ggplot2::geom_point(ggplot2::aes(fill = get(col_var)), size = -1, alpha = 0.25) +
      ggplot2::geom_boxplot(ggplot2::aes(fill = get(col_var)),
                            color = edge_color, alpha = 0.25, staplewidth = 0.5,
                            na.rm = TRUE, outliers = FALSE, show.legend = FALSE) +
      ggplot2::scale_fill_manual(name = col_var, values = fill_colors) +
      ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(shape = 22, size = 10)))
    
  } else {
    fill_color <- if (is.null(colors_vec) || !all(vapply(colors_vec, is_valid_color, logical(1)))) {
      boxplot_fill
    } else {
      colors_vec[1]
    }
    
    plt <- 
      ggplot2::ggplot(data = dt_sco,
                      mapping = ggplot2::aes(x = get(group_var), y = get(metric))) +
      ggplot2::geom_hline(yintercept = 0, color = hline_color, linetype = "solid") +
      ggplot2::geom_boxplot(fill = fill_color, 
                            color = edge_color, alpha = 0.25, 
                            outliers = FALSE, staplewidth = 0.5) +
      ggplot2::theme(legend.position = "none")
  }
  
  # validation for coloring points
  if (colored_pts_flag) {
    if (NROW(unique(dt_sco[[point_var]])) >= 10) colored_pts_flag <- FALSE
  }
  # jitter point
  if (colored_pts_flag) {
    color_points <- get_qual_colors(NROW(unique(dt_sco[[point_var]])))
    names(color_points) <- unique(dt_sco[[point_var]])
    
    plt <- plt +
      ggplot2::geom_jitter(mapping = ggplot2::aes(color = get(point_var)), size = 2, alpha = 0.75,
                           width = 0.2, height = 0, na.rm = TRUE) +
      ggplot2::scale_color_manual(name = drug_name, values = color_points) +
      ggplot2::guides(color = ggplot2::guide_legend(title = point_var))
  } else {
    plt <- plt +
      ggplot2::geom_jitter(color = jitter_poinst_color, alpha = 0.75,
                           width = 0.2, height = 0, na.rm = TRUE)
  }
  
  # final
  plt <- plt +
    ggplot2::labs(title = plt_title,
                  y = get_hm_title(metric, normalization_type), 
                  x = "") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8, angle = 90, vjust = 1, hjust = 1),
                   axis.text.y = ggplot2::element_text(size = 8),
                   plot.title = ggplot2::element_text(size = 10),
                   panel.grid.minor = ggplot2::element_blank())
  
  return(plt)
}

#' Plot box plots for metric for combo data grouped by cell line names
#' 
#' @inheritParams plot_boxplot_metric_combo
#' 
#' @return \code{ggplot} object containing boxplots for selected combo metric 
#'   grouped by cell line names
#' 
#' @keywords combo_plots
#' 
#' @author Janina Smoła \email{janina.smola@@contractors.roche.com}
#' 
#' @examples
#' mae <- gDRutils::get_synthetic_data("combo_matrix")
#' se <- mae[[gDRutils::get_supported_experiments("combo")]]
#' 
#' dt_scores <- gDRutils::convert_se_assay_to_dt(se = se,
#'                                               assay_name = "scores")
#' 
#' plot_boxplot_metric_combo_by_CLs(dt_scores)
#' 
#' plot_boxplot_metric_combo_by_CLs(dt_scores,
#'                                  normalization_type = "RV",
#'                                  grouped_flag = TRUE)
#' 
#' plot_boxplot_metric_combo_by_CLs(dt_scores,
#'                                  metric = "bliss_score",
#'                                  colors_vec = "gold")
#' 
#' plot_boxplot_metric_combo_by_CLs(
#'   dt_scores,
#'   metric = "bliss_score",
#'   colored_pts_flag = TRUE)
#' 
#' plot_boxplot_metric_combo_by_CLs(
#'   dt_scores,
#'   metric = "hsa_score",
#'   normalization_type = "RV",
#'   grouped_flag = TRUE,
#'   colors_vec = c("tomato", "darkgreen", "orange", "darkblue"))
#' 
#' @export
plot_boxplot_metric_combo_by_CLs <- function(
    dt_scores,
    normalization_type = "GR",
    metric = "hsa_score",
    fit_source = "gDR",
    grouped_flag = FALSE,
    colored_pts_flag = FALSE,
    colors_vec = NULL
) {
  
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  tissue <- gDRutils::get_env_identifiers("cellline_tissue")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  drug_name_2 <- gDRutils::get_env_identifiers("drug_name2")
  
  checkmate::assert_data_table(dt_scores)
  checkmate::assert_choice(normalization_type, choices = c("GR", "RV"))
  checkmate::assert_choice(metric, choices = c("hsa_score", "bliss_score"))
  checkmate::assert_names(names(dt_scores), 
                          must.include = c(cellline_name, tissue, drug_name, drug_name_2, metric))
  checkmate::assert_string(fit_source, null.ok = TRUE)
  checkmate::assert_flag(grouped_flag)
  checkmate::assert_flag(colored_pts_flag)
  checkmate::assert_character(colors_vec, null.ok = TRUE)
  boxplot_fill <- 
    gDRutils::get_settings_from_json("BOXPLOT_FILL",
                                     system.file(package = "gDRplots", "settings.json"))
  hline_color <- 
    gDRutils::get_settings_from_json("HLINE_COLOR",
                                     system.file(package = "gDRplots", "settings.json"))
  jitter_poinst_color <- 
    gDRutils::get_settings_from_json("JITTER_POINST_COLOR",
                                     system.file(package = "gDRplots", "settings.json"))
  edge_color <- 
    gDRutils::get_settings_from_json("EDGE_COLOR",
                                     system.file(package = "gDRplots", "settings.json"))
  
  plt <- 
    plot_boxplot_metric_combo(
      dt_scores = dt_scores,
      group_var = cellline_name,
      normalization_type = normalization_type,
      metric = metric,
      fit_source = fit_source,
      grouped_flag = grouped_flag,
      colored_pts_flag = colored_pts_flag,
      colors_vec = colors_vec
    )
  
  return(plt)
}

#' Plot box plots for metric for combo data grouped by drug names
#' 
#' @inheritParams plot_boxplot_metric_combo_by_CLs
#' 
#' @return \code{ggplot} object containing boxplots for selected combo metric
#'    grouped by drug names
#' 
#' @keywords combo_plots
#' 
#' @author Janina Smoła \email{janina.smola@@contractors.roche.com}
#' 
#' @examples
#' mae <- gDRutils::get_synthetic_data("combo_matrix")
#' se <- mae[[gDRutils::get_supported_experiments("combo")]]
#' 
#' dt_scores <- gDRutils::convert_se_assay_to_dt(se = se,
#'                                               assay_name = "scores")
#' 
#' plot_boxplot_metric_combo_by_drugs(dt_scores)
#' 
#' plot_boxplot_metric_combo_by_drugs(dt_scores,
#'                                    normalization_type = "RV",
#'                                    colors_vec = "gold")
#' 
#' plot_boxplot_metric_combo_by_drugs(dt_scores,
#'                                    metric = "bliss_score",
#'                                    normalization_type = "RV",
#'                                    colored_pts_flag = TRUE)
#' 
#' @export
plot_boxplot_metric_combo_by_drugs <- function(
    dt_scores,
    normalization_type = "GR",
    metric = "hsa_score",
    fit_source = "gDR",
    colored_pts_flag = FALSE,
    colors_vec = NULL
) {
  
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  drug_name_2 <- gDRutils::get_env_identifiers("drug_name2")
  
  checkmate::assert_data_table(dt_scores)
  checkmate::assert_choice(normalization_type, choices = c("GR", "RV"))
  checkmate::assert_choice(metric, choices = c("hsa_score", "bliss_score"))
  checkmate::assert_names(names(dt_scores), 
                          must.include = c(cellline_name, drug_name, drug_name_2, metric))
  checkmate::assert_string(fit_source, null.ok = TRUE)
  checkmate::assert_flag(colored_pts_flag)
  checkmate::assert_character(colors_vec, null.ok = TRUE)
  boxplot_fill <- 
    gDRutils::get_settings_from_json("BOXPLOT_FILL",
                                     system.file(package = "gDRplots", "settings.json"))
  hline_color <- 
    gDRutils::get_settings_from_json("HLINE_COLOR",
                                     system.file(package = "gDRplots", "settings.json"))
  jitter_poinst_color <- 
    gDRutils::get_settings_from_json("JITTER_POINST_COLOR",
                                     system.file(package = "gDRplots", "settings.json"))
  edge_color <- 
    gDRutils::get_settings_from_json("EDGE_COLOR",
                                     system.file(package = "gDRplots", "settings.json"))
  
  plt <- 
    plot_boxplot_metric_combo(
      dt_scores = dt_scores,
      group_var = drug_name,
      normalization_type = normalization_type,
      metric = metric,
      fit_source = fit_source,
      grouped_flag = FALSE,
      colored_pts_flag = colored_pts_flag,
      colors_vec = colors_vec
    )
  
  return(plt)
}

#' Plot box plots for metric for combo data grouped by selected variable
#' 
#' @inheritParams plot_boxplot_metric_combo
#' 
#' @param selection_var string name of selected main variable - one value from column 
#'    \code{"CellLineName"} or \code{"DrugName"}
#' @param selection_name string name of selected variable value from column \code{selection_var} 
#'    to filter data for plotting;
#' @param group_var string name of group variable; not numeric variablee from \code{dt_metrics} 
#'    different than \code{selection_var} and not containing unique values for each row;
#' @param group_names character vector with names to subset from column \code{group_var};
#'    if \code{NULL} then all values will be plotted
#' @param named_n_bottom number of n-bottom \code{metric} values to be labeled on the plot
#'    for \code{group_var} equal \code{"DrugName"} points will be labeled by \code{"CellLineName"}
#'    and similarly vice versa
#' @param grouped_flag logical flag whether the boxplots should be colored by \code{group_var}
#' 
#' @return \code{ggplot} object containing boxplots for selected combo metric 
#'   grouped by selected variable
#' 
#' @keywords combo_plots
#' 
#' @author Janina Smoła \email{janina.smola@@contractors.roche.com}
#' 
#' @examples
#' mae <- gDRutils::get_synthetic_data("combo_matrix")
#' se <- mae[[gDRutils::get_supported_experiments("combo")]]
#' 
#' dt_scores <- gDRutils::convert_se_assay_to_dt(se = se,
#'                                               assay_name = "scores")
#' invisible(dt_scores[, Tissue_grp := data.table::fifelse(Tissue == "tissue_w", 
#'                                                         "tissue_w", 
#'                                                         "tissue_other")])
#' 
#' plot_boxplot_metric_combo_by_grp(dt_scores,
#'                                  selection_var = "DrugName",
#'                                  selection_name = c("drug_001", "drug_021"),
#'                                  group_var = "Tissue")
#' 
#' plot_boxplot_metric_combo_by_grp(dt_scores,
#'                                  selection_var = "DrugName",
#'                                  selection_name = c("drug_001", "drug_021"),
#'                                  group_var = "Tissue_grp",
#'                                  grouped_flag = TRUE,
#'                                  colors_vec = c("darkblue", "deeppink"))
#' 
#' plot_boxplot_metric_combo_by_grp(dt_scores,
#'                                  selection_var = "DrugName",
#'                                  selection_name = c("drug_001", "drug_021"),
#'                                  group_var = "Tissue",
#'                                  grouped_flag = TRUE)
#' 
#' plot_boxplot_metric_combo_by_grp(dt_scores,
#'                                  selection_var = "DrugName",
#'                                  selection_name = c("drug_001", "drug_021"),
#'                                  group_var = "Tissue",
#'                                  colors_vec = c("darkblue", "orange"))
#' 
#' @export
plot_boxplot_metric_combo_by_grp <- function(
    dt_scores,
    selection_var,
    selection_name,
    group_var,
    group_names = NULL,
    normalization_type = "GR",
    metric = "hsa_score",
    fit_source = "gDR",
    named_n_bottom = 5,
    grouped_flag = FALSE,
    colors_vec = NULL
) {
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  drug_name_2 <- gDRutils::get_env_identifiers("drug_name2")
  numeric_columns <- names(dt_scores)[vapply(dt_scores, is.numeric, logical(1))]
  
  checkmate::assert_data_table(dt_scores)
  checkmate::assert_choice(selection_var, choices = c(cellline_name, drug_name))
  
  # check input data
  if (selection_var == cellline_name) {
    checkmate::assert_string(selection_name)
    checkmate::assert_choice(selection_name, choices = unique(dt_scores[[selection_var]]))
    point_var <- "DrugCombination"
  } else if (selection_var == drug_name) {
    checkmate::assert_character(selection_name, len = 2, any.missing = FALSE)
    checkmate::assert_true(
      NROW(dt_scores[(get(drug_name) == selection_name[1] & get(drug_name_2) == selection_name[2]) |
                       (get(drug_name) == selection_name[2] & get(drug_name_2) == selection_name[1])]) > 1)
    point_var <- cellline_name
  }
  checkmate::assert_string(group_var)
  checkmate::assert_choice(
    group_var, 
    choices = names(dt_scores)[!names(dt_scores) %in% c(numeric_columns, cellline_name, drug_name, drug_name_2)])
  # TODO add validation for number of levels >1 and !=NROW(dt_scores)
  checkmate::assert_choice(group_names, choices = unique(dt_scores[[group_var]]), null.ok = TRUE)
  checkmate::assert_choice(normalization_type, choices = c("GR", "RV"))
  checkmate::assert_choice(metric, choices = c("hsa_score", "bliss_score"))
  checkmate::assert_names(names(dt_scores), 
                          must.include = c(cellline_name, drug_name, drug_name_2, group_var, metric))
  checkmate::assert_string(fit_source, null.ok = TRUE)
  checkmate::assert_flag(grouped_flag)
  checkmate::assert_character(colors_vec, null.ok = TRUE)
  boxplot_fill <- 
    gDRutils::get_settings_from_json("BOXPLOT_FILL",
                                     system.file(package = "gDRplots", "settings.json"))
  hline_color <- 
    gDRutils::get_settings_from_json("HLINE_COLOR",
                                     system.file(package = "gDRplots", "settings.json"))
  jitter_poinst_color <- 
    gDRutils::get_settings_from_json("JITTER_POINST_COLOR",
                                     system.file(package = "gDRplots", "settings.json"))
  edge_color <- 
    gDRutils::get_settings_from_json("EDGE_COLOR",
                                     system.file(package = "gDRplots", "settings.json"))
  
  # filter data for normalization type
  filter_expr <- substitute(normalization_type == norm_type & fit_source == fit_src,
                            list(norm_type = normalization_type, fit_src = fit_source))
  dt_sco <- dt_scores[eval(filter_expr)]
  # filter by selection
  if (selection_var == cellline_name) {
    dt_sco <- dt_sco[get(selection_var) == selection_name, ]
  } else {
    dt_sco <-
      dt_sco[(get(drug_name) == selection_name[1] & get(drug_name_2) == selection_name[2]) |
               (get(drug_name) == selection_name[2] & get(drug_name_2) == selection_name[1])]
  }
  dt_sco$DrugCombination <-
    paste(dt_sco[[drug_name]], "x", dt_sco[[drug_name_2]])
  # select min required data for plotting
  dt_sco <- dt_sco[, c(group_var, point_var, metric), with = FALSE]
  
  # coloring points by rank
  data.table::setorderv(dt_sco, cols = metric)
  dt_sco[, `:=`(is_bottom = FALSE, label = "")]
  if (named_n_bottom > 0) {
    named_n_bottom <- min(named_n_bottom, NROW(dt_sco)) # deal with less than n bigger than table
    dt_sco <- 
      dt_sco[order(get(metric)), ][seq_len(named_n_bottom), `:=`(is_bottom = TRUE, label = get(point_var))]
  }
  
  # update group (it depends on user choice for `group_names` and `selection_name`)
  group_names <- if (is.null(group_names)) {
    sort(unique(dt_sco[[group_var]]))
  } else {
    sort(intersect(unique(dt_sco[[group_var]]), group_names))
  }
  
  plt_title <- sprintf("%s for %s by %s", 
                       gDRplots::get_hm_title(metric, normalization_type), 
                       paste(selection_name, collapse = " + "), 
                       group_var)
  
  dt_sco_lbl <- data.table::copy(dt_sco)[!is.na(get(metric)), ]
  data.table::setorderv(dt_sco, group_var)
  dt_sco[[group_var]] <- factor(dt_sco[[group_var]], levels = group_names)
  
  
  if (grouped_flag) {
    fill_colors <- if (is.null(colors_vec) || !all(vapply(colors_vec, is_valid_color, logical(1)))) {
      get_qual_colors(NROW(unique(dt_sco[[group_var]])))
    } else if (NROW(colors_vec) != NROW(unique(dt_sco[[group_var]]))) {
      grDevices::colorRampPalette(colors_vec)(NROW(unique(dt_sco[[group_var]])))
    } else {
      colors_vec
    }
    names(fill_colors) <- unique(dt_sco[[group_var]])
    
    plt <-
      ggplot2::ggplot(data = dt_sco,
                      mapping = ggplot2::aes(x = get(group_var), 
                                             y = get(metric))) +
      ggplot2::geom_hline(yintercept = 0, color = hline_color, linetype = "solid") +
      ggplot2::geom_boxplot(ggplot2::aes(fill = get(group_var)),
                            color = edge_color, alpha = 0.25, staplewidth = 0.5,
                            na.rm = TRUE, outliers = FALSE, show.legend = FALSE) +
      ggplot2::scale_fill_manual(name = group_var, values = fill_colors) +
      ggplot2::guides(fill = "none")
    
  } else {
    fill_color <- if (is.null(colors_vec) || !all(vapply(colors_vec, is_valid_color, logical(1)))) {
      boxplot_fill
    } else {
      colors_vec[1]
    }
    
    plt <- 
      ggplot2::ggplot(data = dt_sco,
                      mapping = ggplot2::aes(x = get(group_var), y = get(metric))) +
      ggplot2::geom_hline(yintercept = 0, color = hline_color, linetype = "solid") +
      ggplot2::geom_boxplot(fill = fill_color, 
                            color = edge_color, alpha = 0.25, 
                            outliers = FALSE, staplewidth = 0.5) 
  }
  
  # adding lbls and points colored when bottom
  if (named_n_bottom > 0) {
    plt <- plt +
      ggrepel::geom_text_repel(data = dt_sco_lbl,
                               mapping = ggplot2::aes(x = get(group_var), 
                                                      y = get(metric),
                                                      label = label),
                               size = 3, max.overlaps = 20, show.legend = FALSE) +
      ggplot2::geom_hline(yintercept = max(dt_sco[is_bottom == TRUE][[metric]]), 
                          color = "red", linetype = "dashed")
  }
  
  plt <- plt +
    ggplot2::geom_jitter(mapping = ggplot2::aes(color = is_bottom), 
                         size = 2, alpha = 0.75,
                         width = 0.2, height = 0, na.rm = TRUE, 
                         show.legend = (named_n_bottom > 0)) +
    ggplot2::scale_color_manual(values = c("TRUE" = "red", "FALSE" = jitter_poinst_color))
  
  # final
  plt <- plt +
    ggplot2::labs(title = plt_title,
                  y = get_hm_title(metric, normalization_type), 
                  x = "",
                  color = sprintf("Bottom %s", named_n_bottom)) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8, angle = 90, vjust = 1, hjust = 1),
                   axis.text.y = ggplot2::element_text(size = 8),
                   plot.title = ggplot2::element_text(size = 10),
                   panel.grid.minor = ggplot2::element_blank())
  
  return(plt) 
}
