#' Plot box plots for metric for single-agent data grouped by cell line names
#' 
#' @param dt_metrics data.table representing data from the \code{Metrics} assay,
#'    outputted by \code{gDRutils::convert_se_assay_to_dt(se, "Metrics")}
#'    and single-agent \code{SummarizedExperiment}
#' @param normalization_type string with normalization types to be selected
#'                           one of: "GR" ("GRvalue") or "RV" ("RelativeViability")
#' @param metric string name of the metric;
#'    one of: "xc50" ("GR50" or "IC50" - respectively depending on \code{normalization_type}), 
#'    "x_max" ("GR Max" or "E Max") or "x_mean" ("GR Mean" or "RV Mean");
#'    but the values from any numeric colum can be displayed.
#' @param fit_source string source name for metrics
#' @param grouped_flag a logical flag whether the boxplots should be grouped and 
#'    colored by \code{Tissue}
#' @param colors_vec character vector with colors (name or hex value) to color boxplots
#' 
#' @return \code{ggplot} object containing boxplots for selected single-agent grouped by cellline names
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
    colors_vec = NULL
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
  checkmate::assert_character(colors_vec, null.ok = TRUE)
  
  # filter data for normalization type
  filter_expr <- substitute(normalization_type == norm_type & fit_source == fit_src,
                            list(norm_type = normalization_type, fit_src = fit_source))
  dt_met <- dt_metrics[eval(filter_expr)]
  dt_met <- dt_met[, c(cellline_name, tissue, drug_name, metric), with = FALSE]
  # handle -Inf (NA will be not shown on boxplots)
  dt_met[[metric]] <- 
    vapply(dt_met[[metric]], function(x) ifelse(is.infinite(x), NA, x), numeric(1))
  
  plt_title <- sprintf("Number of unique drugs: %s", NROW(unique(dt_met[[drug_name]])))
  
  if (grouped_flag) {
    data.table::setorderv(dt_met, tissue)
    dt_met[[cellline_name]] <- factor(dt_met[[cellline_name]], levels = unique(dt_met[[cellline_name]]))

    fill_colors <- if (is.null(colors_vec) || !all(vapply(colors_vec, is_valid_color, logical(1)))) {
      get_qual_colors(NROW(unique(dt_met[[tissue]])))
    } else if (NROW(colors_vec) != NROW(unique(dt_met[[tissue]]))) {
      grDevices::colorRampPalette(colors_vec)(NROW(unique(dt_met[[tissue]])))
    } else {
      colors_vec
    }
    names(fill_colors) <- unique(dt_met[[tissue]])
    
    plt <- 
      ggplot2::ggplot(data = dt_met,
                      mapping = ggplot2::aes(x = get(cellline_name), y = get(metric))) +
      ggplot2::geom_hline(yintercept = 0, color = "#B3B3B3", linetype = "solid") +
      ggplot2::geom_point(ggplot2::aes(fill = get(tissue)), size = -1, alpha = 0.25, na.rm = TRUE) +
      ggplot2::geom_boxplot(ggplot2::aes(fill = get(tissue)), 
                            color = "#A9A9A9", alpha = 0.25, show.legend = FALSE) +
      ggplot2::scale_fill_manual(name = tissue, values = fill_colors) +
      ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(shape = 22, size = 10)))
    
  } else {
    data.table::setorderv(dt_met, cellline_name)
    dt_met[[cellline_name]] <- factor(dt_met[[cellline_name]])
    
    fill_color <- if (is.null(colors_vec) || !all(vapply(colors_vec, is_valid_color, logical(1)))) {
      "#A6CEE3"
    } else {
      colors_vec[1]
    }
    
    plt <- 
      ggplot2::ggplot(data = dt_met,
                      mapping = ggplot2::aes(x = get(cellline_name), y = get(metric))) +
      ggplot2::geom_hline(yintercept = 0, color = "#B3B3B3", linetype = "solid") +
      ggplot2::geom_boxplot(fill = fill_color, color = "#A9A9A9", alpha = 0.25, na.rm = TRUE) +
      ggplot2::theme(legend.position = "none")
  }
  
  # final
  plt <- plt +
    ggplot2::geom_jitter(width = 0.2, height = 0, color = "#4C4C4C", na.rm = TRUE) +
    ggplot2::labs(title = plt_title,
                  y = sprintf("%s for %s", metric, normalization_type), 
                  x = "") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8, angle = 90, vjust = 1, hjust = 1),
                   axis.text.y = ggplot2::element_text(size = 8),
                   plot.title = ggplot2::element_text(size = 10),
                   panel.grid.minor = ggplot2::element_blank())
  
  return(plt)
}

#' Plot box plots for metric for single-agent data grouped by drug names
#' 
#' @inheritParams plot_boxplot_metric_sa_by_CLs
#' 
#' @return \code{ggplot} object containing boxplots for selected single-agent grouped by drug names
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
    colors_vec = NULL
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
  checkmate::assert_character(colors_vec, null.ok = TRUE)
  
  # filter data for normalization type
  filter_expr <- substitute(normalization_type == norm_type & fit_source == fit_src,
                            list(norm_type = normalization_type, fit_src = fit_source))
  dt_met <- dt_metrics[eval(filter_expr)]
  dt_met <- dt_met[, c(cellline_name, drug_name, drug_MOA, metric), with = FALSE]

  # handle -Inf (NA will be not shown on boxplots)
  dt_met[[metric]] <- 
    vapply(dt_met[[metric]], function(x) ifelse(is.infinite(x), NA, x), numeric(1))
  
  plt_title <- sprintf("Number of unique celllines: %s", NROW(unique(dt_met[[cellline_name]])))
  
  if (grouped_flag) {
    data.table::setorderv(dt_met, drug_MOA)
    dt_met[[drug_name]] <- factor(dt_met[[drug_name]], levels = unique(dt_met[[drug_name]]))

    fill_colors <- if (is.null(colors_vec) || !all(vapply(colors_vec, is_valid_color, logical(1)))) {
      get_qual_colors(NROW(unique(dt_met[[drug_MOA]])))
    } else if (NROW(colors_vec) != NROW(unique(dt_met[[drug_MOA]]))) {
      grDevices::colorRampPalette(colors_vec)(NROW(unique(dt_met[[drug_MOA]])))
    } else {
      colors_vec
    }
    names(fill_colors) <- unique(dt_met[[drug_MOA]])
    
    plt <- 
      ggplot2::ggplot(data = dt_met,
                      mapping = ggplot2::aes(x = get(drug_name), y = get(metric))) +
      ggplot2::geom_hline(yintercept = 0, color = "#B3B3B3", linetype = "solid") +
      ggplot2::geom_point(ggplot2::aes(fill = get(drug_MOA)), size = -1, alpha = 0.25, na.rm = TRUE) +
      ggplot2::geom_boxplot(ggplot2::aes(fill = get(drug_MOA)), 
                            color = "#A9A9A9", alpha = 0.25, show.legend = FALSE, na.rm = TRUE) +
      ggplot2::scale_fill_manual(name = drug_MOA, values = fill_colors) +
      ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(shape = 22, size = 10)))
    
  } else {
    data.table::setorderv(dt_met, drug_name)
    dt_met[[drug_name]] <- factor(dt_met[[drug_name]])

    fill_color <- if (is.null(colors_vec) || !all(vapply(colors_vec, is_valid_color, logical(1)))) {
      "#A6CEE3"
    } else {
      colors_vec[1]
    }
    
    plt <- 
      ggplot2::ggplot(data = dt_met,
                      mapping = ggplot2::aes(x = get(drug_name), y = get(metric))) +
      ggplot2::geom_hline(yintercept = 0, color = "#B3B3B3", linetype = "solid") +
      ggplot2::geom_boxplot(fill = fill_color, color = "#A9A9A9", alpha = 0.25, na.rm = TRUE) +
      ggplot2::theme(legend.position = "none")
  }
  
  # final
  plt <- plt +
    ggplot2::geom_jitter(width = 0.2, height = 0, color = "#4C4C4C", na.rm = TRUE) +
    ggplot2::labs(title = plt_title,
                  y = sprintf("%s for %s", metric, normalization_type), 
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
#' @inheritParams plot_boxplot_metric_sa_by_CLs
#' @param dt_scores \code{data.table} representing data from the \code{scores} assay,
#'   outputted by \code{gDRutils::convert_se_assay_to_dt(se, "scores")}
#'   and combo \code{SummarizedExperiment}
#' @param metric string name of the combo metric;
#'   one of: "hsa_score"("Bliss Excess GR" or "Bliss Excess RV" - respectively 
#'   depending on \code{normalization_type}), "bliss_score" ("Bliss Score GR" or "Bliss Score RV")
#' 
#' @return \code{ggplot} object containing boxplots for selected combo grouped by drug names
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
#'   grouped_flag = TRUE)
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
  checkmate::assert_character(colors_vec, null.ok = TRUE)
  
  # filter data for normalization type
  filter_expr <- substitute(normalization_type == norm_type & fit_source == fit_src,
                            list(norm_type = normalization_type, fit_src = fit_source))
  dt_sco_norm <- dt_scores[eval(filter_expr)]
  
  dt_sco <- dt_sco_norm[, c(cellline_name, tissue, drug_name, drug_name_2, metric), with = FALSE]
  dt_sco$DrugCombination <-
    paste(dt_sco[[drug_name]], "x", dt_sco[[drug_name_2]])
  
  plt_title <- sprintf("Number of unique drug combinations: %s", NROW(unique(dt_sco[["DrugCombination"]])))
  
  if (grouped_flag) {
    data.table::setorderv(dt_sco, tissue)
    dt_sco[[cellline_name]] <- factor(dt_sco[[cellline_name]], levels = unique(dt_sco[[cellline_name]]))
    
    fill_colors <- if (is.null(colors_vec) || !all(vapply(colors_vec, is_valid_color, logical(1)))) {
      get_qual_colors(NROW(unique(dt_sco[[tissue]])))
    } else if (NROW(colors_vec) != NROW(unique(dt_sco[[tissue]]))) {
      grDevices::colorRampPalette(colors_vec)(NROW(unique(dt_sco[[tissue]])))
    } else {
      colors_vec
    }
    names(fill_colors) <- unique(dt_sco[[tissue]])
    
    plt <-
      ggplot2::ggplot(data = dt_sco,
                      mapping = ggplot2::aes(x = get(cellline_name), y = get(metric))) +
      ggplot2::geom_hline(yintercept = 0, color = "#B3B3B3", linetype = "solid") +
      ggplot2::geom_point(ggplot2::aes(fill = get(tissue)), size = -1, alpha = 0.25) +
      ggplot2::geom_boxplot(ggplot2::aes(fill = get(tissue)),
                            color = "#A9A9A9", alpha = 0.25, show.legend = FALSE) +
      ggplot2::scale_fill_manual(name = tissue, values = fill_colors) +
      ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(shape = 22, size = 10)))
    
  } else {
    fill_color <- if (is.null(colors_vec) || !all(vapply(colors_vec, is_valid_color, logical(1)))) {
      "#A6CEE3"
    } else {
      colors_vec[1]
    }
    
    plt <- 
      ggplot2::ggplot(data = dt_sco,
                      mapping = ggplot2::aes(x = get(cellline_name), y = get(metric))) +
      ggplot2::geom_hline(yintercept = 0, color = "#B3B3B3", linetype = "solid") +
      ggplot2::geom_boxplot(fill = fill_color, color = "#A9A9A9", alpha = 0.25) +
      ggplot2::theme(legend.position = "none")
  }
  
  # final
  plt <- plt +
    ggplot2::geom_jitter(width = 0.2, height = 0, color = "#4C4C4C") +
    ggplot2::labs(title = plt_title,
                  y = sprintf("%s for %s", metric, normalization_type), 
                  x = "") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8, angle = 90, vjust = 1, hjust = 1),
                   axis.text.y = ggplot2::element_text(size = 8),
                   plot.title = ggplot2::element_text(size = 10),
                   panel.grid.minor = ggplot2::element_blank())
  
  return(plt)
}

#' Plot box plots for metric for combo data grouped by drug names
#' 
#' @inheritParams plot_boxplot_metric_combo_by_CLs
#' 
#' @return \code{ggplo} object containing boxplots for selected combo grouped by drug names
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
#' @export
plot_boxplot_metric_combo_by_drugs <- function(
    dt_scores,
    normalization_type = "GR",
    metric = "hsa_score",
    fit_source = "gDR",
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
  checkmate::assert_character(colors_vec, null.ok = TRUE)
  
  # filter data for normalization type
  filter_expr <- substitute(normalization_type == norm_type & fit_source == fit_src,
                            list(norm_type = normalization_type, fit_src = fit_source))
  dt_sco_norm <- dt_scores[eval(filter_expr)]
  
  dt_sco <- dt_sco_norm[, c(cellline_name, drug_name, drug_name_2, metric), with = FALSE]
  dt_sco$DrugCombination <-
    paste(dt_sco[[drug_name]], "x", dt_sco[[drug_name_2]])
  
  plt_title <- sprintf("Number of unique celllines: %s", NROW(unique(dt_sco[[cellline_name]])))
  
  fill_color <- if (is.null(colors_vec) || !all(vapply(colors_vec, is_valid_color, logical(1)))) {
    "#A6CEE3"
  } else {
    colors_vec[1]
  }
  
  plt <- 
    ggplot2::ggplot(data = dt_sco,
                    mapping = ggplot2::aes(x = DrugCombination, y = get(metric))) +
    ggplot2::geom_hline(yintercept = 0, color = "#B3B3B3", linetype = "solid") +
    ggplot2::geom_boxplot(fill = fill_color, color = "#A9A9A9", alpha = 0.25) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::geom_jitter(width = 0.2, height = 0, color = "#4C4C4C") +
    ggplot2::labs(title = plt_title,
                  y = sprintf("%s for %s", metric, normalization_type), 
                  x = "") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8, angle = 90, vjust = 1, hjust = 1),
                   axis.text.y = ggplot2::element_text(size = 8),
                   plot.title = ggplot2::element_text(size = 10),
                   panel.grid.minor = ggplot2::element_blank())
  
  return(plt)
}
