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
    message("Please, choose only one coloring options: or `grouped_flag` or `colored_pts_flag`.")
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
  
  plt_title <- sprintf("Number of unique %s: %s", 
                       ifelse(point_var == cellline_name, "celllines", "drugs"),
                       NROW(unique(dt_met[[point_var]])))
  
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
      ggplot2::geom_jitter(mapping = ggplot2::aes(color = get(point_var)), size = 2,
                           width = 0.2, height = 0, na.rm = TRUE) +
      ggplot2::scale_color_manual(name = drug_name, values = color_points) +
      ggplot2::guides(color = ggplot2::guide_legend(title = point_var))
  } else {
    plt <- plt +
      ggplot2::geom_jitter(color = jitter_poinst_color, width = 0.2, height = 0, na.rm = TRUE)
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
#'   grouped_flag = TRUE)
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
#'   grouped_flag = TRUE)
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
#' @param dt_scores \code{data.table} representing data from the \code{scores} assay,
#'   outputted by \code{gDRutils::convert_se_assay_to_dt(se, "scores")}
#'   and combo \code{SummarizedExperiment}
#' @param metric string name of the combo metric;
#'   one of: "hsa_score"("Bliss Excess GR" or "Bliss Excess RV" - respectively 
#'   depending on \code{normalization_type}), "bliss_score" ("Bliss Score GR" or "Bliss Score RV")
#' 
#' @return \code{ggplot} object containing boxplots for selected combo metric 
#'   grouped by \code{group_var}
#' 
#' @keywords combo_plots
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
#'   colors_vec = c("deeppink", "darkcyan", "orange", "darkblue"))
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
    if(grouped_flag) message("Coloring box by group is not available for this scenario.")
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
    message("Please, choose only one coloring options: or `grouped_flag` or `colored_pts_flag`.")
  }
  
  # filter data for normalization type
  filter_expr <- substitute(normalization_type == norm_type & fit_source == fit_src,
                            list(norm_type = normalization_type, fit_src = fit_source))
  dt_sco_norm <- dt_scores[eval(filter_expr)]
  
  dt_sco <- dt_sco_norm[, c(cellline_name, drug_name, drug_name_2, col_var, metric), with = FALSE]
  dt_sco$DrugCombination <-
    paste(dt_sco[[drug_name]], "x", dt_sco[[drug_name_2]])
  
  plt_title <- sprintf("Number of unique %s: %s", 
                       ifelse(point_var == cellline_name, "celllines", "drug combinations"),
                       NROW(unique(dt_sco[[point_var]])))
  
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
      ggplot2::geom_jitter(mapping = ggplot2::aes(color = get(point_var)), size = 2,
                           width = 0.2, height = 0, na.rm = TRUE) +
      ggplot2::scale_color_manual(name = drug_name, values = color_points) +
      ggplot2::guides(color = ggplot2::guide_legend(title = point_var))
  } else {
    plt <- plt +
      ggplot2::geom_jitter(color = jitter_poinst_color, width = 0.2, height = 0, na.rm = TRUE)
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
#'   colors_vec = c("deeppink", "darkcyan", "orange", "darkblue"))
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
