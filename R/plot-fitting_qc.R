#' Lollipop plot for metric single-agent data to control quality of the data
#'
#' @param dt_assay data.table representing data from the \code{Metrics} assay,
#'    outputted by \code{gDRutils::convert_se_assay_to_dt(se, "Metrics")}
#'    and single-agent \code{SummarizedExperiment}
#' @param cl_name string cell line name to be plotted (Cell Line Name)
#' @param metric string with variable name to be plotted; it has to be in \code{dt_assay}
#' @param normalization_type string with normalization_types to be selected
#'                           one of: "GR" ("GRvalue") or "RV" ("RelativeViability")
#' @param with_table logical whether table with metric values should be shown next to the plot
#'
#' @return \code{ggplot} object containing lollipop plot with stat value for each drug
#'
#' @keywords QC_plot
#'
#' @author Bartosz Czech \email{czech.bartosz@@external.gene.com}
#'
#' @examples
#' mae <- gDRutils::get_synthetic_data("small")
#' se <- mae[[gDRutils::get_supported_experiments("sa")]]
#'
#' dt_metrics <- gDRutils::convert_se_assay_to_dt(se, "Metrics")
#' cl_name <- dt_metrics[["CellLineName"]][1]
#'
#' plot_var_stat_qc(dt_assay = dt_metrics,
#'                  cl_name = cl_name)
#'
#' plot_var_stat_qc(dt_assay = dt_metrics,
#'                  cl_name = cl_name,
#'                  metric = "r2",
#'                  normalization_type = "RV")
#'
#' plot_var_stat_qc(dt_assay = dt_metrics,
#'                  cl_name = cl_name,
#'                  metric = "x_AOC",
#'                  normalization_type = "RV",
#'                  with_table = TRUE)
#'
#' @export
plot_var_stat_qc <- function(dt_assay,
                             cl_name,
                             metric = "r2",
                             normalization_type = "GR",
                             with_table = FALSE) {

  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  clid <- gDRutils::get_env_identifiers("cellline")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  gnumber <- gDRutils::get_env_identifiers("drug")

  checkmate::assert_data_table(dt_assay)
  checkmate::assert_string(cl_name)
  checkmate::assert_choice(cl_name, choices = unique(dt_assay[[cellline_name]]))
  checkmate::assert_choice(metric, choices = names(dt_assay))
  checkmate::assert_choice(normalization_type, choices = c("GR", "RV"))
  checkmate::assert_flag(with_table)
  hline_color <-
    .get_setting("HLINE_COLOR")

  cl_clid <- unique(dt_assay[get(cellline_name) == cl_name, clid])

  # filter data for normalization type
  filter_expr <- substitute(normalization_type == norm_type, list(norm_type = normalization_type))
  dt_assay <- dt_assay[eval(filter_expr)]

  tab_subplot <- dt_assay[get(cellline_name) == cl_name, ]

  plt_title <- sprintf("%s (%s)", cl_name, cl_clid)
  color_palette <- get_qual_colors(NROW(unique(tab_subplot[[drug_name]])))

  plt <- ggplot2::ggplot(tab_subplot,
                         ggplot2::aes(x = get(drug_name), y = !!rlang::sym(metric))) +
    ggplot2::geom_hline(yintercept = 1, color = hline_color, linetype = "dashed") +
    ggplot2::geom_hline(yintercept = 0, color = hline_color, linetype = "solid") +
    ggplot2::geom_segment(
      ggplot2::aes(x = get(drug_name), xend = get(drug_name), y = 0, yend = !!rlang::sym(metric))) +
    ggplot2::geom_point(ggplot2::aes(fill = get(drug_name), color = get(drug_name)),
                        alpha = 0.75, size = 5, shape = 21, stroke = 1) +
    ggplot2::theme_minimal() +
    ggplot2::scale_fill_manual(values = color_palette) +
    ggplot2::scale_color_manual(values = color_palette) +
    ggplot2::labs(y = sprintf("%s for %s", metric, normalization_type), x = drug_name, title = plt_title) +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1))

  if (with_table) {
    tab_dt <- tab_subplot[, .SD, .SDcols = c(drug_name, metric)][order(get(metric))]
    tab_plt <- .table_to_ggplot(tab_dt)
    plt <- patchwork::wrap_plots(plt, tab_plt, widths = c(2, 1))
  }

  return(plt)
}


#' Visualization for the quality control of the fitting for single-agent data
#'
#' @param dt_assay data.table representing data from the \code{Metrics} assay,
#'    outputted by \code{gDRutils::convert_se_assay_to_dt(se, "Metrics")}
#'    and single-agent \code{SummarizedExperiment}
#' @param cl_name string cell line name to be plotted (Cell Line Name)
#' @param normalization_type string with normalization_types to be selected
#'                           one of: "GR" ("GRvalue") or "RV" ("RelativeViability")
#'
#' @return \code{ggplot} object containing panel with lollipop plots with r2 and rss values for each drug
#'
#' @keywords QC_plot
#'
#' @author Bartosz Czech \email{czech.bartosz@@external.gene.com}
#'
#' @examples
#' mae <- gDRutils::get_synthetic_data("small")
#' se <- mae[[gDRutils::get_supported_experiments("sa")]]
#'
#' dt_metrics <- gDRutils::convert_se_assay_to_dt(se, "Metrics")
#' cl_name <- dt_metrics[["CellLineName"]][1]
#'
#' plot_fitting_acc(dt_assay = dt_metrics,
#'                  cl_name = cl_name,
#'                  normalization_type = "RV")
#' @export
plot_fitting_acc <- function(dt_assay,
                             cl_name,
                             normalization_type = "GR") {

  cellline_name <- gDRutils::get_env_identifiers("cellline_name")

  checkmate::assert_data_table(dt_assay)
  checkmate::assert_string(cl_name)
  checkmate::assert_choice(cl_name, choices = unique(dt_assay[[cellline_name]]))
  checkmate::assert_choice(normalization_type, choices = c("GR", "RV"))

  if (all(is.na(dt_assay[get(cellline_name) == cl_name, "r2"]))) {
    warning(sprintf("Missing data for %s in %s normalization type.", cl_name, normalization_type))
    return(ggplot2::ggplot() + ggplot2::theme_void())
  }

  r2 <- plot_var_stat_qc(dt_assay,
                         cl_name,
                         metric = "r2",
                         normalization_type = normalization_type,
                         with_table = FALSE) +
    ggplot2::labs(x = NULL) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank())

  r2 <- r2 +
    ggplot2::geom_text(data = subset(r2$data, p_value < 0.05),
                       ggplot2::aes(label = ifelse(p_value < 0.001, "***", ifelse(p_value < 0.01, "**", "*")),
                                    y = 1.01),
                       position = ggplot2::position_dodge(0.5),
                       size = 4,
                       vjust = 0) +
    ggplot2::labs(caption = "*** p < 0.001, ** p < 0.01, * p < 0.05")


  rss <- plot_var_stat_qc(dt_assay,
                          cl_name,
                          metric = "rss",
                          normalization_type = normalization_type,
                          with_table = FALSE) +
    ggplot2::labs(title = NULL)

  rss$layers <- rss$layers[-1]

  # Combine plots vertically (one on top of the other)
  combined_plot <- patchwork::wrap_plots(r2, rss, ncol = 1, heights = c(0.7, 0.7))
  metric_cols <- c("r2", "rss")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  data2table <- r2$data[, c(drug_name, metric_cols), with = FALSE]
  data.table::setorder(data2table, "rss")
  tab_plt <- .table_to_ggplot(data2table, base_size = 8)
  patchwork::wrap_plots(combined_plot,
                       tab_plt,
                       widths = c(2, 1))
}


#' Plot heatmap of mapping controls to treated for single-agent and combo data to control quality of the data
#'
#' @param dt_treat data.table representation of the data in \code{RawTreated} assay,
#'    outputted by \code{gDRutils::convert_se_assay_to_dt(se, "RawTreated")}
#'    and \code{SummarizedExperiment} with chosen data type: single-agent or combo
#' @param dt_controls data.table representation of the data in \code{Controls} assay,
#'    outputted by \code{gDRutils::convert_se_assay_to_dt(se, "Controls")}
#'    and \code{SummarizedExperiment} with chosen data type: single-agent or combo
#'
#' @return \code{pheatmap} object containing hetamap of mapping controls to treated
#'
#' @keywords QC_plot
#'
#' @author Bartosz Czech \email{czech.bartosz@@external.gene.com}
#'
#' @examples
#' mae <- gDRutils::get_synthetic_data("small")
#' se <- mae[[gDRutils::get_supported_experiments("sa")]]
#'
#' dt_treat <- gDRutils::convert_se_assay_to_dt(se, "RawTreated")
#' dt_controls <- gDRutils::convert_se_assay_to_dt(se, "Controls")
#'
#' heatmap_control_mapping_qc(dt_treat = dt_treat,
#'                            dt_controls = dt_controls)
#'
#' heatmap_control_mapping_qc(dt_treat = dt_treat[1:1350, ],
#'                            dt_controls = dt_controls)
#'
#' dt_treat_NA <- dt_treat[-c(1:135, 270:405),]
#' heatmap_control_mapping_qc(dt_treat = dt_treat_NA,
#'                            dt_controls = dt_controls)
#'
#' dt_controls_NA <- dt_controls[-c(1:305, 611:763, 1221:1750),]
#' heatmap_control_mapping_qc(dt_treat = dt_treat,
#'                            dt_controls = dt_controls_NA)
#'
#' heatmap_control_mapping_qc(dt_treat = dt_treat,
#'                            dt_controls = dt_controls[1:3660, ])
#'
#'
#' @export
heatmap_control_mapping_qc <- function(dt_treat,
                                       dt_controls) {
  checkmate::assert_data_table(dt_treat)
  checkmate::assert_data_table(dt_controls)
  qc_heatmap_palette <-
    .get_setting("QC_HEATMAP_PALETTE")

  # calculate the frequency of each (rID, cID) combination in Controls
  frequency <- dt_controls[, .N, by = .(rId, cId)]
  # merge the frequency with the Treated data.table
  result <- merge(unique(dt_treat[, c("rId", "cId")]), frequency, by = c("rId", "cId"), all.x = TRUE)

  # Convert the result to a matrix format suitable for pheatmap
  result_matrix <- data.table::dcast(result, rId ~ cId, value.var = "N")
  rownames <- result_matrix$rId
  result_matrix <- as.matrix(result_matrix[, !("rId"), with = FALSE])
  rownames(result_matrix) <- rownames

  # Replace 0 with NA to use na_col for red color
  result_matrix[result_matrix == 0] <- NA

  # Generate the breaks for integers
  maxval <- max(result_matrix, na.rm = TRUE)
  minval <- min(c(0, min(result_matrix, na.rm = TRUE)))

  breaks <- seq(from = minval, to = maxval, by = 1)
  hm_color_palette <- grDevices::colorRampPalette(qc_heatmap_palette)(NROW(breaks) - 1)

  unique_values <- unique(stats::na.omit(as.vector(result_matrix)))

  # renaming rows and columns
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  drug_name_2 <- gDRutils::get_env_identifiers("drug_name2")
  if (all(colnames(result_matrix) %in% unique(dt_treat$cId))) {
    dict <- unique(dt_treat[, .SD, .SDcols = c("cId", cellline_name)])[order(colnames(result_matrix))]
    colnames(result_matrix) <- dict[[cellline_name]]
  }
  if (all(rownames(result_matrix) %in% unique(dt_treat$rId))) {
    if (drug_name_2 %in% names(dt_treat)) {
      # combo data
      dict <-
        unique(dt_treat[, .SD, .SDcols = c("rId", drug_name, drug_name_2)])[order(rownames(result_matrix))]
      dict$drug_lbl <- paste(dict[[drug_name]], dict[[drug_name_2]], sep = " x ")
      rownames(result_matrix) <- dict[["drug_lbl"]]
    } else {
      # single-agent data
      dict <- unique(dt_treat[, .SD, .SDcols = c("rId", drug_name)])[order(rownames(result_matrix))]
      rownames(result_matrix) <- dict[[drug_name]]
    }
  }

  # Generate the heatmap
  pheatmap::pheatmap(result_matrix,
                     color = hm_color_palette,
                     breaks = breaks,
                     na_col = "red",
                     main = "Counts of mapped controls",
                     cluster_cols = FALSE,
                     cluster_rows = FALSE,
                     angle_col = 90,
                     legend_breaks = unique_values,
                     legend_labels = unique_values)
}

#' Render a data.table as a text-based ggplot table
#'
#' @param dt data.table to render
#' @param base_size numeric base font size (default 10)
#' @return ggplot object displaying the table
#'
#' @examples
#' dt <- data.table::data.table(drug = c("A", "B"), r2 = c(0.9, 0.8))
#' .table_to_ggplot(dt)
#'
#' @keywords internal
.table_to_ggplot <- function(dt, base_size = 10) {
  checkmate::assert_data_table(dt, min.rows = 1)
  checkmate::assert_number(base_size, lower = 1)
  dt <- data.table::copy(dt)
  cols <- names(dt)
  dt[, .row := .N - seq_len(.N) + 1L]
  header_y <- max(dt$.row) + 1L
  labels <- lapply(cols, function(col) {
    vals <- as.character(dt[[col]])
    vals <- ifelse(is.na(vals), "", vals)
    data.table::data.table(
      x = match(col, cols),
      y = c(header_y, dt$.row),
      label = c(col, vals),
      fontface = c("bold", rep("plain", NROW(dt)))
    )
  })
  lbl_dt <- data.table::rbindlist(labels)
  ggplot2::ggplot(lbl_dt, ggplot2::aes(x = x, y = y, label = label)) +
    ggplot2::geom_text(ggplot2::aes(fontface = fontface),
                       size = base_size / 3, hjust = 0.5) +
    ggplot2::theme_void() +
    ggplot2::coord_cartesian(clip = "off")
}
