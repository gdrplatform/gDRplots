#' Plot drug response curves for single-agent data
#'
#' @param dt_metrics data.table representing data from the \code{Metrics} assay,
#'    outputted by \code{gDRutils::convert_se_assay_to_dt(se, "Metrics")}
#' @param dt_average data.table representing data from the \code{Averaged} assay,
#'    outputted by \code{gDRutils::convert_se_assay_to_dt(se, "Averaged")}
#' @param grouping string name of dim of se to plot; it has to be opppsoite to this used to filter se
#'    (that is when rownames(se) == 1 it has to be "cId", otherwise - "rId")
#' @param group_names character vector with names to subset from se (the same dim as \code{grouping});
#'    if \code{NULL} then all values will be plotted
#' @param normalization_type string with normalization_types to be selected
#'                           one of: "GR" ("GRvalue") or "RV" ("RelativeViability")
#' @param colormap character vector with colors for \code{group_names} - name or hex value
#' @param plot_averaged_flag logical flag whether plot points with average values 
#' @param plot_fit_flag logical flag whether plot points with fitted values 
#'
#' @return plot with dose-response curves
#'
#' @keywords single-agent_plots
#' @examples
#' mae <- gDRutils::get_synthetic_data("small")
#' se <- mae[[gDRutils::get_supported_experiments("sa")]]
#' iR <- rownames(se)[1]
#' grouping <- "CellLineName"
#' dt_metrics <- gDRutils::convert_se_assay_to_dt(se[iR], "Metrics")
#' dt_average <- gDRutils::convert_se_assay_to_dt(se[iR], "Averaged")
#' group_names <- unique(dt_metrics[[grouping]])[1:3]
#' 
#' plot_dose_response_sa(dt_metrics = dt_metrics, 
#'                       dt_average = dt_average, 
#'                       grouping = grouping,
#'                       group_names = group_names)
#'             
#' iC <- colnames(se)[1]                    
#' grouping <- "DrugName"
#' dt_metrics <- gDRutils::convert_se_assay_to_dt(se[, iC], "Metrics")
#' dt_average <- gDRutils::convert_se_assay_to_dt(se[, iC], "Averaged")
#' group_names <- unique(dt_metrics[[grouping]])[1:3]
#' 
#' plot_dose_response_sa(dt_metrics = dt_metrics, 
#'                       dt_average = dt_average, 
#'                       grouping = grouping,
#'                       group_names = group_names)
#'                       
#' @export
plot_dose_response_sa <- function(dt_metrics, 
                                  dt_average, 
                                  grouping, 
                                  group_names = NULL, 
                                  normalization_type = "GR", 
                                  colormap = NULL, 
                                  plot_averaged_flag = TRUE, 
                                  plot_fit_flag = TRUE) {
  
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  clid <- gDRutils::get_env_identifiers("cellline")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  gnumber <- gDRutils::get_env_identifiers("drug")
  
  checkmate::expect_data_table(dt_metrics)
  checkmate::expect_data_table(dt_average)
  checkmate::expect_choice(grouping, choices = c(cellline_name, drug_name))
  checkmate::expect_character(group_names, null.ok = TRUE)
  checkmate::expect_choice(normalization_type, choices = c("GR", "RV"))
  checkmate::expect_character(colormap, null.ok = TRUE)
  checkmate::expect_flag(plot_averaged_flag)
  checkmate::expect_flag(plot_fit_flag)
  
  # check input data
  if (grouping == cellline_name) {
    stopifnot("grouping` does not fit to `dt_metrics` and `dt_average`" =
                (NROW(unique(dt_metrics[[drug_name]])) == 1 && NROW(unique(dt_average[[drug_name]])) == 1))
  } else if (grouping == drug_name) {
    stopifnot("grouping` does not fit to `dt_metrics` and `dt_average`" =
                (NROW(unique(dt_metrics[[cellline_name]])) == 1 && NROW(unique(dt_average[[cellline_name]])) == 1))
  }
  stopifnot("empty plot was selected" = any(plot_averaged_flag, plot_fit_flag))
  
  # filter data for normalization type
  data.table::setkeyv(dt_metrics, "normalization_type")
  dt_met_norm <- dt_metrics[normalization_type]
  data.table::setkey(dt_met_norm, NULL)
  data.table::setkeyv(dt_average, "normalization_type")
  dt_avg_norm <- dt_average[normalization_type]
  data.table::setkey(dt_avg_norm, NULL)
  
  # variables
  conc <- gDRutils::get_env_identifiers("concentration")
  
  # prep value ranges for plot
  data_range <- c(min(min(dt_avg_norm$x), 0) - 0.05, max(max(dt_avg_norm$x), 1) + 0.05)
  min_conc <- min(dt_avg_norm[dt_avg_norm[[conc]] > 0, ][[conc]])
  max_conc <- max(dt_avg_norm[[conc]])
  conc_range <- 0.5 * c(floor(2 * log10(min_conc) - 0.5), ceiling(2 * log10(max(max_conc)) + 0.3)) 
  # remove conc = 0
  dt_avg_norm[[conc]][dt_avg_norm[[conc]] == 0] <- min_conc / 10
  
  # group
  if (is.null(group_names)) group_names <- unique(dt_met_norm[[grouping]])
  
  # prep fitted data
  sel_conc <- 10 ^ (seq(conc_range[1], conc_range[2], 0.05))
  dt_fit <- data.table::data.table()
  
  for (icol in group_names) {
    sel_metrics <- dt_met_norm[get(grouping) == icol, ]
    dt_fit <- rbind(dt_fit, 
                    cbind(sel_metrics[, grouping, with = FALSE],
                          data.table::data.table(
                            conc_col = sel_conc,
                            x = gDRutils::predict_efficacy_from_conc(sel_conc, 
                                                                     sel_metrics$x_inf, 
                                                                     sel_metrics$x_0, 
                                                                     sel_metrics$ec50, 
                                                                     sel_metrics$h)
                          )
                    )
    )
  }
  dt_fit <- dt_fit[!is.na(x)]
  data.table::setnames(dt_fit, "conc_col", conc)
  
  # colors
  if (is.null(colormap) || !all(vapply(colormap, is_valid_color, logical(1)))) {
    colormap <- grDevices::colorRampPalette(c("coral", "chartreuse", "darkblue"))(NROW(group_names))
    
  } else if (NROW(colormap) != NROW(group_names)) {
    colormap <- grDevices::colorRampPalette(colormap)(NROW(group_names))
  }
  names(colormap) <- group_names
  color_values <- colormap[names(colormap) %in% unique(dt_fit[[grouping]])]
  
  # prep averaged data
  dt_avg <- dt_avg_norm[get(grouping) %in% group_names, ]
  
  # levels
  dt_avg$grouping <- factor(dt_avg[[grouping]], levels = group_names)
  dt_fit$grouping <- factor(dt_fit[[grouping]], levels = group_names)
  
  plt_title <- sprintf(
    "%s (%s)",
    ifelse(grouping == cellline_name, unique(dt_metrics[[drug_name]]), unique(dt_metrics[[cellline_name]])),
    ifelse(grouping == cellline_name, unique(dt_metrics[[gnumber]]), unique(dt_metrics[[clid]]))
  )
  
  # final plot
  plt <- 
    ggplot2::ggplot(mapping = ggplot2::aes(x = log10(get(conc)), y = x, color = grouping, group = grouping)) +
    ggplot2::geom_hline(yintercept = c(0, 1), color = "#555555") +
    ggplot2::scale_color_manual(values = color_values,
                                name = ifelse(grouping == "cId", "Cell Line", "Drug")) +
    ggplot2::coord_cartesian(xlim = conc_range, ylim = data_range) +
    ggplot2::scale_x_continuous(breaks = -5:2, labels = c("1e-5", "1e-4", 10 ^ (-3:2))) +
    ggplot2::xlab(bquote(.(conc) ~ "[" ~ mu * M ~ "]")) +
    ggplot2::ylab(sprintf("log10(%s)", normalization_type)) +
    ggplot2::ggtitle(plt_title) +
    ggplot2::theme_bw()
  
  if (plot_averaged_flag) {
    plt <- plt + ggplot2::geom_point(data = dt_avg)
  }
  
  if (plot_fit_flag) {
    plt <- plt + ggplot2::geom_line(data = dt_fit)
  }
  
  # define legend
  plt <- plt + 
    ggplot2::guides(colour = ggplot2::guide_legend(position = "left"))
  
  return(plt)
}


#' Plot drug response curves for single-agent data for selected call lines and drugs
#' 
#' @inheritParams plot_dose_response_sa
#' @param cellline_name_vec character vector with cell line to be plotted (Cell Line Name)
#' @param drug_name_vec character vector with cell line to be plotted (Drug Name)
#'    
#' @return list of plots with dose-response curves
#' 
#' @keywords single-agent_plots
#' @examples
#' mae <- gDRutils::get_synthetic_data("small")
#' se <- mae[[gDRutils::get_supported_experiments("sa")]]
#' dt_metrics <- gDRutils::convert_se_assay_to_dt(se, "Metrics")
#' dt_average <- gDRutils::convert_se_assay_to_dt(se, "Averaged")
#' cellline_name_vec <- unique(dt_metrics[["CellLineName"]])[2:5]
#' drug_name_vec <- unique(dt_metrics[["DrugName"]])[5:7]
#' 
#' plot_dose_response_sa_by_CLs(dt_metrics = dt_metrics, 
#'                              dt_average = dt_average, 
#'                              cellline_name_vec = cellline_name_vec, 
#'                              drug_name_vec = drug_name_vec, 
#'                              normalization_type = "RV", 
#'                              colormap = c("#B9D3EE", "#FF6347", "#C2F970"))
#' 
#' @export
plot_dose_response_sa_by_CLs <- function(dt_metrics, 
                                         dt_average,
                                         cellline_name_vec = NULL, 
                                         drug_name_vec = NULL,
                                         normalization_type = "GR", 
                                         colormap = NULL, 
                                         plot_averaged_flag = TRUE, 
                                         plot_fit_flag = TRUE) {
  
  checkmate::expect_data_table(dt_metrics)
  checkmate::expect_data_table(dt_average)
  checkmate::expect_character(cellline_name_vec, null.ok = TRUE)
  checkmate::expect_character(drug_name_vec, null.ok = TRUE)
  checkmate::expect_choice(normalization_type, choices = c("GR", "RV"))
  checkmate::expect_character(colormap, null.ok = TRUE)
  checkmate::expect_flag(plot_averaged_flag)
  checkmate::expect_flag(plot_fit_flag)
  
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  clid <- gDRutils::get_env_identifiers("cellline")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  gnumber <- gDRutils::get_env_identifiers("drug")
  
  available_drugs <- unique(dt_metrics[[drug_name]])
  if (is.null(drug_name_vec) || all(!drug_name_vec %in% available_drugs)) {
    drug_name_vec  <- available_drugs
  } else if (!all(drug_name_vec %in% available_drugs)) {
    drug_name_vec <- drug_name_vec[drug_name_vec  %in% available_drugs]
  }  
  
  available_cellline <- unique(dt_metrics[[cellline_name]])
  if (is.null(cellline_name_vec) || all(!cellline_name_vec %in% available_cellline)) {
    cellline_name_vec <- available_cellline
  } else if (!all(cellline_name_vec %in% available_cellline)) {
    cellline_name_vec <- cellline_name_vec[cellline_name_vec  %in% available_cellline]
  } 
  
  plt_list <- list()
  for (iR in drug_name_vec) {
    
    # subset data
    dt_metrics_subset <- dt_metrics[get(drug_name) == iR & get(cellline_name) %in% cellline_name_vec]
    dt_average_subset <- dt_average[get(drug_name) == iR & get(cellline_name) %in% cellline_name_vec]
    
    plt_title <- sprintf("%s (%s)", iR, unique(dt_metrics[get(drug_name) == iR, ][[gnumber]]))
    
    plt <- 
      plot_dose_response_sa(dt_metrics = dt_metrics_subset, 
                            dt_average = dt_average_subset, 
                            grouping = cellline_name,
                            group_names = cellline_name_vec,
                            normalization_type = normalization_type,
                            colormap = colormap,
                            plot_averaged_flag = plot_averaged_flag,
                            plot_fit_flag = plot_fit_flag)
    
    plt_list[[plt_title]] <- plt
  }
  
  return(plt_list)
}

#' Plot drug response curves for single-agent data for selected call lines and drugs
#' 
#' @inheritParams plot_dose_response_sa
#' @param cellline_name_vec character vector with cell line to be plotted (Cell Line Name)
#' @param drug_name_vec character vector with cell line to be plotted (Drug Name)
#'    
#' @return list of plots with dose-response curves
#' 
#' @keywords single-agent_plots
#' @examples
#' mae <- gDRutils::get_synthetic_data("small")
#' se <- mae[[gDRutils::get_supported_experiments("sa")]]
#' 
#' dt_metrics <- gDRutils::convert_se_assay_to_dt(se, "Metrics")
#' dt_average <- gDRutils::convert_se_assay_to_dt(se, "Averaged")
#' cellline_name_vec <- unique(dt_metrics[["CellLineName"]])[2:5]
#' drug_name_vec <- unique(dt_metrics[["DrugName"]])[5:7]
#' 
#' plot_dose_response_sa_by_drugs(dt_metrics = dt_metrics, 
#'                                dt_average = dt_average, 
#'                                cellline_name_vec = cellline_name_vec, 
#'                                drug_name_vec = drug_name_vec, 
#'                                normalization_type = "RV", 
#'                                colormap = c("#B9D3EE", "#FF6347", "#C2F970"))
#' 
#' @export
plot_dose_response_sa_by_drugs <- function(dt_metrics, 
                                           dt_average,
                                           cellline_name_vec = NULL, 
                                           drug_name_vec = NULL,
                                           normalization_type = "GR", 
                                           colormap = NULL, 
                                           plot_averaged_flag = TRUE, 
                                           plot_fit_flag = TRUE) {
  
  checkmate::expect_data_table(dt_metrics)
  checkmate::expect_data_table(dt_average)
  checkmate::expect_character(cellline_name_vec, null.ok = TRUE)
  checkmate::expect_character(drug_name_vec, null.ok = TRUE)
  checkmate::expect_choice(normalization_type, choices = c("GR", "RV"))
  checkmate::expect_character(colormap, null.ok = TRUE)
  checkmate::expect_flag(plot_averaged_flag)
  checkmate::expect_flag(plot_fit_flag)
  
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  clid <- gDRutils::get_env_identifiers("cellline")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  gnumber <- gDRutils::get_env_identifiers("drug")
  
  available_cellline <- unique(dt_metrics[[cellline_name]])
  if (is.null(cellline_name_vec) || all(!cellline_name_vec %in% available_cellline)) {
    cellline_name_vec <- available_cellline
  } else if (!all(cellline_name_vec %in% available_cellline)) {
    cellline_name_vec <- cellline_name_vec[cellline_name_vec  %in% available_cellline]
  }  
  
  available_drugs <- unique(dt_metrics[[drug_name]])
  if (is.null(drug_name_vec) || all(!drug_name_vec %in% available_drugs)) {
    drug_name_vec  <- available_drugs
  } else if (!all(drug_name_vec %in% available_drugs)) {
    drug_name_vec <- drug_name_vec[drug_name_vec  %in% available_drugs]
  }  
  
  plt_list <- list()
  for (iC in cellline_name_vec) {
    
    # subset data
    dt_metrics_subset <- dt_metrics[get(cellline_name) == iC & get(drug_name) %in% drug_name_vec]
    dt_average_subset <- dt_average[get(cellline_name) == iC & get(drug_name) %in% drug_name_vec]
    
    plt_title <- sprintf("%s (%s)", iC, unique(dt_metrics[get(cellline_name) == iC, ][[clid]]))
    
    plt <- 
      plot_dose_response_sa(dt_metrics = dt_metrics_subset, 
                            dt_average = dt_average_subset, 
                            grouping = drug_name,
                            group_names = drug_name_vec,
                            normalization_type = normalization_type,
                            colormap = colormap,
                            plot_averaged_flag = plot_averaged_flag,
                            plot_fit_flag = plot_fit_flag)
    
    plt_list[[plt_title]] <- plt
  }
  
  return(plt_list)
}

#' Plot drug response curves for single-agent data to control quality of the data
#'
#' @param dt_metrics data.table representing data from the \code{Metrics} assay,
#'    outputted by \code{gDRutils::convert_se_assay_to_dt(se, "Metrics")}
#' @param dt_average data.table representing data from the \code{Averaged} assay,
#'    outputted by \code{gDRutils::convert_se_assay_to_dt(se, "Averaged")}
#' @param cl_name string cell line name to be plotted (Cell Line Name)
#' @param d_name string vector with drug name to be plotted (Drug Name)
#' @param normalization_type string with normalization_types to be selected
#'                           one of: "GR" ("GRvalue") or "RV" ("RelativeViability")
#' @param fit_source string source name for metrics
#'
#' @return plot with dose-response curves
#'
#' @keywords QC_plot
#' @examples
#' mae <- gDRutils::get_synthetic_data("small")
#' se <- mae[[gDRutils::get_supported_experiments("sa")]]
#'
#' dt_metrics <- gDRutils::convert_se_assay_to_dt(se, "Metrics")
#' dt_average <- gDRutils::convert_se_assay_to_dt(se, "Averaged")
#' cl_name <- dt_metrics[["CellLineName"]][1]
#' d_name <- dt_metrics[["DrugName"]][1]
#' 
#' plot_dose_response_sa_qc(dt_metrics = dt_metrics,
#'                          dt_average = dt_average,
#'                          cl_name = cl_name,
#'                          d_name = d_name)
#' 
#' @export
plot_dose_response_sa_qc <- function(dt_metrics, 
                                     dt_average, 
                                     cl_name, 
                                     d_name,
                                     normalization_type = "GR",
                                     fit_source = "gDR") {
  
  checkmate::expect_data_table(dt_metrics)
  checkmate::expect_data_table(dt_average)
  checkmate::expect_string(cl_name)
  checkmate::expect_string(d_name)
  checkmate::expect_choice(normalization_type, choices = c("GR", "RV"))
  checkmate::assert_string(fit_source, null.ok = TRUE)
  
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  clid <- gDRutils::get_env_identifiers("cellline")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  gnumber <- gDRutils::get_env_identifiers("drug")
  conc <- gDRutils::get_env_identifiers("concentration")
  
  checkmate::expect_choice(cl_name, choices = dt_metrics[[cellline_name]])
  checkmate::expect_choice(cl_name, choices = dt_average[[cellline_name]])
  checkmate::expect_choice(d_name, choices = dt_metrics[[drug_name]])
  checkmate::expect_choice(d_name, choices = dt_average[[drug_name]])

  # filter data for normalization_type
  data.table::setkeyv(dt_metrics, "normalization_type")
  dt_metrics <- dt_metrics[normalization_type]
  data.table::setkey(dt_metrics, NULL)
  data.table::setkeyv(dt_average, "normalization_type")
  dt_average <- dt_average[normalization_type]
  data.table::setkey(dt_average, NULL)
  
  # filter data for min required data
  dt_metrics <-
    dt_metrics[, .SD, .SDcols = c(drug_name, gnumber, cellline_name, clid, "x_inf", "x_0", "ec50", "h")]
  dt_average <- 
    dt_average[, .SD, .SDcols = c(drug_name, gnumber, cellline_name, clid, conc, "x", "x_std")]
  
  selected_combination <- data.table::data.table(cellline_name = cl_name,
                                                 drug_name = d_name)
  data.table::setnames(selected_combination, c("cellline_name", "drug_name"), c(cellline_name, drug_name))
  
  # tab_plots
  dt_average_plot <- dt_average[selected_combination, on = c(cellline_name, drug_name)]
  dt_metrics_plot <- dt_metrics[selected_combination, on = c(cellline_name, drug_name)]
  
  if (NROW(dt_metrics_plot) > 0) {
    min_conc <- min(dt_average_plot[get(conc) != 0][[conc]])
    max_conc <- max(dt_average_plot[[conc]])
    sampled_conc <- gDRplots::create_log_seq(min_conc, max_conc, 50) 
    fitted_curve_sampled <- gDRutils::predict_efficacy_from_conc(sampled_conc,
                                                                 dt_metrics_plot$x_inf,
                                                                 dt_metrics_plot$x_0,
                                                                 dt_metrics_plot$ec50,
                                                                 dt_metrics_plot$h)
    dt_reconstructed_fit <- data.table::data.table(
      Concentration = sampled_conc,
      x = fitted_curve_sampled
    )
    
    # set min and max values for y 
    ymin <- min(c(0, min(dt_average_plot$x)))
    ymax <- max(c(1.2, max(dt_average_plot$x)))
    
    plt_title <- sprintf("%s (%s)", dt_metrics_plot[[drug_name]], dt_metrics_plot[[gnumber]])
    
    # plot
    plt <- 
      ggplot2::ggplot() +
      ggplot2::geom_errorbar(
        data = dt_average_plot, 
        ggplot2::aes(x = get(conc), y = x,  ymin = x - x_std, ymax = x + x_std, color = "Errors Bar"), 
        width = 0.1) + 
      ggplot2::geom_line(
        data = dt_average_plot, 
        ggplot2::aes(x = get(conc), y = x, color = "Averaged Data"), 
        linetype = "dashed") +
      ggplot2::geom_line(
        data = dt_reconstructed_fit, 
        ggplot2::aes(x = get(conc), y = x, color = "Fitted Curve")) +
      ggplot2::geom_hline(yintercept = 0, color = "#A9A9A9") +
      ggplot2::scale_x_continuous(trans = "log10") +
      ggplot2::scale_y_continuous(lim = c(ymin, ymax)) +
      ggplot2::xlab(bquote(.(conc) ~ "[" ~ mu * M ~ "]")) +
      ggplot2::ylab(normalization_type) + 
      ggplot2::ggtitle(plt_title) +
      ggplot2::labs(color = "") +
      ggplot2::theme_bw() +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) +
      ggplot2::scale_color_manual(values = c("Averaged Data" = "black", 
                                             "Fitted Curve" = "red", 
                                             "Errors Bar" = "#A9A9A9"))
  } else {
    txt_err <- sprintf(
      "Dose response curve \nfor Drug Name: %s (%s) and CellLine: %s (%s) \n could not be calculated.",
      dt_metrics_plot[[drug_name]], 
      dt_metrics_plot[[gnumber]], 
      dt_metrics_plot[[cellline_name]], 
      dt_metrics_plot[[clid]])
    plt <- 
      ggplot2::ggplot() +
      ggplot2::geom_text(ggplot2::aes(x = 0, y = 0, label = txt_err),
                         color = "darkred", size = 5) + 
      ggplot2::theme_void()
  }
  
  return(plt)
}

#' Plot panel with drug response curves for single-agent data to control quality of the data
#'
#' @inheritParams plot_dose_response_sa_qc
#' @param d_names character vector with drug names to be plotted (Drug Name); 
#'    if NULL - all available drugs will be plotted
#'
#' @return panle with plot with dose-response curves for selected cell line by drugs
#'
#' @keywords QC_plot
#' @examples
#' mae <- gDRutils::get_synthetic_data("small")
#' se <- mae[[1]]
#' 
#' dt_metrics <- gDRutils::convert_se_assay_to_dt(se, "Metrics")
#' dt_average <- gDRutils::convert_se_assay_to_dt(se, "Averaged")
#' cl_name <- dt_metrics[["CellLineName"]][1]
#' d_names <- unique(dt_metrics[["DrugName"]])[1:3]
#' 
#' plot_dose_response_sa_qc_panel(dt_metrics = dt_metrics,
#'                                dt_average = dt_average,
#'                                cl_name = cl_name)
#' 
#' plot_dose_response_sa_qc_panel(dt_metrics = dt_metrics,
#'                                dt_average = dt_average,
#'                                cl_name = cl_name,
#'                                d_names = d_names)
#' 
#' @export
plot_dose_response_sa_qc_panel <- function(dt_metrics, 
                                           dt_average, 
                                           cl_name, 
                                           d_names = NULL,
                                           normalization_type = "GR",
                                           fit_source = "gDR") {
  
  checkmate::expect_data_table(dt_metrics)
  checkmate::expect_data_table(dt_average)
  checkmate::expect_string(cl_name)
  checkmate::expect_character(d_names, null.ok = TRUE)
  checkmate::expect_choice(normalization_type, choices = c("GR", "RV"))
  checkmate::assert_string(fit_source, null.ok = TRUE)
  
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  clid <- gDRutils::get_env_identifiers("cellline")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  gnumber <- gDRutils::get_env_identifiers("drug")
  conc <- gDRutils::get_env_identifiers("concentration")
  
  checkmate::expect_choice(cl_name, choices = dt_metrics[[cellline_name]])
  checkmate::expect_choice(cl_name, choices = dt_average[[cellline_name]])
  
  available_drugs <- unique(dt_metrics[[drug_name]])
  if (is.null(d_names) || all(!d_names %in% available_drugs)) {
    d_names  <- available_drugs
  } else if (!all(d_names %in% available_drugs)) {
    d_names <- drug_name[drug_name %in% available_drugs]
  } 
  
  ls_drug <- list(d_name = d_names)
  
  # list of plots for each drug
  ls_plt <- purrr::pmap(ls_drug, 
                        gDRplots::plot_dose_response_sa_qc,
                        dt_metrics = dt_metrics,
                        dt_average = dt_average,
                        cl_name = cl_name,
                        normalization_type = normalization_type,
                        fit_source = fit_source)
  
  # panel title
  cl_clid <- unique(dt_metrics[get(cellline_name) == cl_name, ][[clid]]) 
  panel_title <- sprintf("%s (%s)", cl_name, cl_clid)
  
  # final panel
  panel <- ggpubr::annotate_figure(
    ggpubr::ggarrange(plotlist = ls_plt, common.legend = TRUE),
    top = panel_title) + 
    ggpubr::bgcolor("white") + ggpubr::border("white")
  
  return(panel)
}
