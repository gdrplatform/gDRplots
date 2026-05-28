#' Plot drug response curves for single-agent data
#'
#' @note inspired by the \code{grob_SA} function written by Marc Hafner
#'
#' @param dt_metrics data.table representing data from the \code{Metrics} assay,
#'    outputted by \code{gDRutils::convert_se_assay_to_dt(se, "Metrics")}
#'    and single-agent \code{SummarizedExperiment}
#' @param dt_average data.table representing data from the \code{Averaged} assay,
#'    outputted by \code{gDRutils::convert_se_assay_to_dt(se, "Averaged")}
#'    and single-agent \code{SummarizedExperiment}
#' @param selection_name string name of selected main variable - one value from column
#'    \code{"CellLineName"} or \code{"DrugName"}
#' @param group_var string name of group variable; one of: \code{"CellLineName"} or \code{"DrugName"}
#' @param group_names character vector with names to subset from se (the same dim as \code{group_var});
#'    if \code{NULL} then all values will be plotted
#' @param normalization_type string with normalization_types to be selected
#'                           one of: "GR" ("GRvalue") or "RV" ("RelativeViability")
#' @param colors_vec character vector with colors for \code{group_names} - name or hex value
#' @param plot_averaged_flag logical flag whether plot points with average values
#' @param plot_fit_flag logical flag whether plot points with fitted values
#' @param fit_source string source name for metrics
#'
#' @return \code{ggplot} object containing plot of dose-response curves
#'
#' @keywords single-agent_plots
#' @examples
#' mae <- gDRutils::get_synthetic_data("small")
#' se <- mae[[gDRutils::get_supported_experiments("sa")]]
#' selected_drug <- "drug_002"
#' group_var <- "CellLineName"
#' dt_metrics <- gDRutils::convert_se_assay_to_dt(se, "Metrics")
#' dt_average <- gDRutils::convert_se_assay_to_dt(se, "Averaged")
#' celline_names <- unique(dt_metrics[[group_var]])[1:3]
#'
#' plot_dose_response_sa(dt_metrics = dt_metrics,
#'                       dt_average = dt_average,
#'                       selection_name = selected_drug,
#'                       group_var = group_var,
#'                       group_names = celline_names)
#'
#' plot_dose_response_sa(dt_metrics = dt_metrics,
#'                       dt_average = NULL,
#'                       selection_name = selected_drug,
#'                       group_var = group_var,
#'                       group_names = celline_names)
#'
#' selected_cellline <- "cellline_HB"
#' group_var <- "DrugName"
#' dt_metrics <- gDRutils::convert_se_assay_to_dt(se, "Metrics")
#' dt_average <- gDRutils::convert_se_assay_to_dt(se, "Averaged")
#' group_names <- unique(dt_metrics[[group_var]])[1:3]
#'
#' plot_dose_response_sa(dt_metrics = dt_metrics,
#'                       dt_average = dt_average,
#'                       selection_name = selected_cellline,
#'                       group_var = group_var,
#'                       group_names = group_names)
#'
#' @export
plot_dose_response_sa <- function(dt_metrics,
                                  dt_average,
                                  selection_name,
                                  group_var,
                                  group_names = NULL,
                                  normalization_type = "GR",
                                  colors_vec = NULL,
                                  plot_averaged_flag = TRUE,
                                  plot_fit_flag = TRUE,
                                  fit_source = "gDR") {

  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  clid <- gDRutils::get_env_identifiers("cellline")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  gnumber <- gDRutils::get_env_identifiers("drug")
  conc <- gDRutils::get_env_identifiers("concentration")
  zero_conc_scaling_factor <-
    gDRutils::get_settings_from_json("ZERO_CONC_SCALING_FACTOR",
                                     system.file(package = "gDRplots", "settings.json"))
  hline_color <-
    gDRutils::get_settings_from_json("HLINE_COLOR",
                                     system.file(package = "gDRplots", "settings.json"))

  checkmate::assert_data_table(dt_metrics)
  checkmate::assert_data_table(dt_average, null.ok = TRUE)
  checkmate::assert_string(selection_name)
  checkmate::assert_choice(group_var, choices = c(cellline_name, drug_name))
  checkmate::assert_character(group_names, null.ok = TRUE)
  checkmate::assert_choice(normalization_type, choices = c("GR", "RV"))
  checkmate::assert_character(colors_vec, null.ok = TRUE)
  checkmate::assert_flag(plot_averaged_flag)
  checkmate::assert_flag(plot_fit_flag)
  checkmate::assert_string(fit_source, null.ok = TRUE)

  # check input data
  if (group_var == cellline_name) {
    main_var <- drug_name
  } else if (group_var == drug_name) {
    main_var <- cellline_name
  }
  stopifnot("Empty plot was selected" = any(plot_averaged_flag, plot_fit_flag))

  # filter data for normalization type
  filter_expr <- substitute(normalization_type == norm_type & fit_source == fit_src,
                            list(norm_type = normalization_type, fit_src = fit_source))
  dt_met_norm <- dt_metrics[eval(filter_expr)]
  dt_avg_norm <- dt_average[eval(filter_expr)]

  # filter data for selected main variable
  dt_met <- dt_met_norm[get(main_var) == selection_name, ]
  dt_avg <- if (is.null(dt_avg_norm)) {
    NULL
  } else {
    dt_avg_norm[get(main_var) == selection_name, ]
  }

  # update group (it depends on user choice for `group_names` and `selection_name`)
  group_names <- if (is.null(group_names)) {
    unique(c(unique(dt_met[[group_var]]), unique(dt_avg[[group_var]])))
  } else {
    intersect(unique(c(unique(dt_met[[group_var]]), unique(dt_avg[[group_var]]))),
              group_names)
  }

  # filter data
  dt_met <- dt_met[get(group_var) %in% group_names, ][!is.na(x_inf) & !is.na(x_0) & !is.na(ec50) & !is.na(h)]
  if (!is.null(dt_avg)) dt_avg <- dt_avg[get(group_var) %in% group_names, ][!is.na(x), ]

  # if conc is NA
  if (all(is.na(unique(dt_avg[[conc]])))) dt_avg <- NULL
  if (!NROW(dt_met)) dt_met <- NULL

  # plot title
  if (NROW(dt_met) == 0 && NROW(dt_avg) == 0) {
    plt_title <- selection_name
  } else {
    dt_src <- if (NROW(dt_met) == 0) dt_avg else dt_met
    dt_src <- unique(dt_src[get(main_var) == selection_name, c(drug_name, gnumber, cellline_name, clid), with = FALSE])

    if (group_var == cellline_name) {
      title_name <- unique(dt_src[[drug_name]])
      title_id <- unique(dt_src[[gnumber]])
    } else {
      title_name <- unique(dt_src[[cellline_name]])
      title_id <- unique(dt_src[[clid]])
    }
    plt_title <- sprintf("%s (%s)", title_name, title_id)
  }

  if (NROW(dt_met) == 0 && NROW(dt_avg) == 0) {
    plt <-
      ggplot2::ggplot() +
      ggplot2::theme(aspect.ratio = 1)
  } else {
    # prep value ranges for x-axis
    if (is.null(dt_avg)) {
      min_conc <- 1e-3
      max_conc <- max(10 ^ dt_met[["maxlog10Concentration"]])
      if (is.na(max_conc)) max_conc <- 30
    } else {
      min_conc <- min(dt_avg[dt_avg[[conc]] > 0, ][[conc]], na.rm = TRUE)
      max_conc <- max(dt_avg[[conc]], na.rm = TRUE)
      # handle conc = 0
      dt_avg[[conc]][dt_avg[[conc]] == 0] <- min_conc / zero_conc_scaling_factor
    }
    conc_range <- 0.5 * c(floor(2 * log10(min_conc) - 0.5), ceiling(2 * log10(max(max_conc)) + 0.3))

    # prep fitted data
    sel_conc <- 10 ^ (seq(conc_range[1], conc_range[2], 0.05))

    if (!NROW(dt_met)) {
      dt_fit <- NULL
    } else {
      dt_fit <- data.table::data.table()
      for (grp_nm in group_names) {
        sel_metrics <- dt_met[get(group_var) == grp_nm, ]
        if (NROW(sel_metrics) == 0) next
        dt_fit <- rbind(dt_fit,
                        cbind(sel_metrics[, group_var, with = FALSE],
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
      if (!NROW(dt_fit)) {
        dt_fit <- NULL
      } else {
        dt_fit <- dt_fit[!is.na(x)]
        data.table::setnames(dt_fit, "conc_col", conc)
      }
    }

    # prep value ranges for y-axis
    min_val <- min(c(dt_avg$x, dt_fit$x, 0), na.rm = TRUE) - 0.05
    max_val <- max(c(dt_avg$x, dt_fit$x, 0, 1), na.rm = TRUE) + 0.05
    data_range <- c(min_val, max_val)

    # prep color palette
    color_values <- if (is.null(colors_vec) || !all(vapply(colors_vec, is_valid_color, logical(1)))) {
      get_qual_colors(NROW(group_names))
    } else if (NROW(colors_vec) != NROW(group_names)) {
      grDevices::colorRampPalette(colors_vec)(NROW(group_names))
    } else {
      colors_vec
    }
    names(color_values) <- group_names

    # levels
    if (!is.null(dt_avg)) dt_avg$group_var <- factor(dt_avg[[group_var]], levels = group_names)
    if (!is.null(dt_fit)) dt_fit$group_var <- factor(dt_fit[[group_var]], levels = group_names)

    # final plot
    plt <-
      ggplot2::ggplot(mapping = ggplot2::aes(x = log10(get(conc)), y = x, color = group_var, group = group_var)) +
      ggplot2::geom_hline(yintercept = c(-1, 0, 1), color = hline_color) +
      ggplot2::geom_hline(yintercept = 0.5, color = hline_color, linetype = "dashed") +
      ggplot2::scale_color_manual(values = color_values,
                                  name = ifelse(group_var == cellline_name, "Cell Line", "Drug")) +
      ggplot2::coord_cartesian(xlim = conc_range, ylim = data_range) +
      ggplot2::scale_x_continuous(breaks = -5:2, labels = c("1e-5", "1e-4", 10 ^ (-3:2))) +
      ggplot2::scale_y_continuous(breaks = seq(-1, 1, by = 0.25))


    if (plot_averaged_flag && !is.null(dt_avg)) {
      plt <- plt + ggplot2::geom_point(data = dt_avg, na.rm = TRUE)
    }

    if (plot_fit_flag && !is.null(dt_fit)) {
      plt <- plt + ggplot2::geom_line(data = dt_fit, na.rm = TRUE)
    }

    # define legend
    plt <- plt +
      ggplot2::guides(color = ggplot2::guide_legend(position = "left"))
  }

  plt <- plt +
    ggplot2::xlab(bquote(.(conc) ~ "[" ~ mu * M ~ "]")) +
    ggplot2::ylab(normalization_type) +
    ggplot2::ggtitle(plt_title) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
                   axis.text.y = ggplot2::element_text(size = 8),
                   plot.title = ggplot2::element_text(size = 10),
                   aspect.ratio = 1)

  return(plt)
}


#' Plot drug response curves for single-agent data for selected call lines and drugs
#'
#' @inheritParams plot_dose_response_sa
#' @param cellline_name_vec character vector with cell line names to be plotted (Cell Line Name)
#' @param drug_name_vec character vector with drug names to be plotted (Drug Name)
#'
#' @return list of \code{ggplot objects} containing plots of dose-response curves
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
#'                              colors_vec = c("#00008B", "#FF6347", "#4CBB17"))
#'
#' @export
plot_dose_response_sa_by_CLs <- function(dt_metrics,
                                         dt_average,
                                         cellline_name_vec = NULL,
                                         drug_name_vec = NULL,
                                         normalization_type = "GR",
                                         colors_vec = NULL,
                                         plot_averaged_flag = TRUE,
                                         plot_fit_flag = TRUE,
                                         fit_source = "gDR") {

  checkmate::assert_data_table(dt_metrics)
  checkmate::assert_data_table(dt_average)
  checkmate::assert_character(cellline_name_vec, null.ok = TRUE)
  checkmate::assert_character(drug_name_vec, null.ok = TRUE)
  checkmate::assert_choice(normalization_type, choices = c("GR", "RV"))
  checkmate::assert_character(colors_vec, null.ok = TRUE)
  checkmate::assert_flag(plot_averaged_flag)
  checkmate::assert_flag(plot_fit_flag)
  checkmate::assert_string(fit_source, null.ok = TRUE)

  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  drug_name <- gDRutils::get_env_identifiers("drug_name")

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
  for (d_name in drug_name_vec) {

    plt_list[[d_name]] <-
      plot_dose_response_sa(dt_metrics = dt_metrics,
                            dt_average = dt_average,
                            selection_name = d_name,
                            group_var = cellline_name,
                            group_names = cellline_name_vec,
                            normalization_type = normalization_type,
                            colors_vec = colors_vec,
                            plot_averaged_flag = plot_averaged_flag,
                            plot_fit_flag = plot_fit_flag)
  }

  return(plt_list)
}

#' Plot drug response curves for single-agent data for selected call lines and drugs
#'
#' @inheritParams plot_dose_response_sa
#' @param cellline_name_vec character vector with cell line names to be plotted (Cell Line Name)
#' @param drug_name_vec character vector with drug names to be plotted (Drug Name)
#'
#' @return list of \code{ggplot objects} containing plots of dose-response curves
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
#'                                colors_vec = c("#00008B", "#FF6347", "#4CBB17"))
#'
#' @export
plot_dose_response_sa_by_drugs <- function(dt_metrics,
                                           dt_average,
                                           cellline_name_vec = NULL,
                                           drug_name_vec = NULL,
                                           normalization_type = "GR",
                                           colors_vec = NULL,
                                           plot_averaged_flag = TRUE,
                                           plot_fit_flag = TRUE,
                                           fit_source = "gDR") {

  checkmate::assert_data_table(dt_metrics)
  checkmate::assert_data_table(dt_average)
  checkmate::assert_character(cellline_name_vec, null.ok = TRUE)
  checkmate::assert_character(drug_name_vec, null.ok = TRUE)
  checkmate::assert_choice(normalization_type, choices = c("GR", "RV"))
  checkmate::assert_character(colors_vec, null.ok = TRUE)
  checkmate::assert_flag(plot_averaged_flag)
  checkmate::assert_flag(plot_fit_flag)
  checkmate::assert_string(fit_source, null.ok = TRUE)

  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  drug_name <- gDRutils::get_env_identifiers("drug_name")

  # filter data for normalization_type and fit_source
  filter_expr <- substitute(normalization_type == norm_type & fit_source == fit_src,
                            list(norm_type = normalization_type, fit_src = fit_source))
  dt_metrics <- dt_metrics[eval(filter_expr)]
  dt_average <- dt_average[eval(filter_expr)]

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
  for (cl_name in cellline_name_vec) {

    plt_list[[cl_name]] <-
      plot_dose_response_sa(dt_metrics = dt_metrics,
                            dt_average = dt_average,
                            selection_name = cl_name,
                            group_var = drug_name,
                            group_names = drug_name_vec,
                            normalization_type = normalization_type,
                            colors_vec = colors_vec,
                            plot_averaged_flag = plot_averaged_flag,
                            plot_fit_flag = plot_fit_flag)
  }

  return(plt_list)
}

#' Plot drug response curves for single-agent data to control quality of the data
#'
#' @param dt_metrics data.table representing data from the \code{Metrics} assay,
#'    outputted by \code{gDRutils::convert_se_assay_to_dt(se, "Metrics")}
#'    and single-agent \code{SummarizedExperiment}
#' @param dt_average data.table representing data from the \code{Averaged} assay,
#'    outputted by \code{gDRutils::convert_se_assay_to_dt(se, "Averaged")}
#'    and single-agent \code{SummarizedExperiment}
#' @param cl_name string cell line name to be plotted (Cell Line Name)
#' @param d_name string vector with drug name to be plotted (Drug Name)
#' @param normalization_type string with normalization_types to be selected
#'                           one of: "GR" ("GRvalue") or "RV" ("RelativeViability")
#' @param fit_source string source name for metrics
#'
#' @return \code{ggplot} object containing plot of dose-response curves (observed and fitted values)
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

  checkmate::assert_data_table(dt_metrics)
  checkmate::assert_data_table(dt_average)
  checkmate::assert_string(cl_name)
  checkmate::assert_string(d_name)
  checkmate::assert_choice(normalization_type, choices = c("GR", "RV"))
  checkmate::assert_string(fit_source, null.ok = TRUE)

  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  clid <- gDRutils::get_env_identifiers("cellline")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  gnumber <- gDRutils::get_env_identifiers("drug")
  conc <- gDRutils::get_env_identifiers("concentration")

  checkmate::assert_choice(cl_name, choices = dt_metrics[[cellline_name]])
  checkmate::assert_choice(cl_name, choices = dt_average[[cellline_name]])
  checkmate::assert_choice(d_name, choices = dt_metrics[[drug_name]])
  checkmate::assert_choice(d_name, choices = dt_average[[drug_name]])
  hline_color <-
    gDRutils::get_settings_from_json("HLINE_COLOR",
                                     system.file(package = "gDRplots", "settings.json"))

  # filter data for normalization_type and fit_source
  filter_expr <- substitute(normalization_type == norm_type & fit_source == fit_src,
                            list(norm_type = normalization_type, fit_src = fit_source))
  dt_metrics <- dt_metrics[eval(filter_expr)]
  dt_average <- dt_average[eval(filter_expr)]

  # filter data for min required data
  dt_metrics <-
    dt_metrics[, c(drug_name, gnumber, cellline_name, clid, "x_inf", "x_0", "ec50", "h"), with = FALSE]
  dt_average <-
    dt_average[, c(drug_name, gnumber, cellline_name, clid, conc, "x", "x_std"), with = FALSE]

  selected_combination <- data.table::data.table(cellline_name = cl_name,
                                                 drug_name = d_name)
  data.table::setnames(selected_combination, c("cellline_name", "drug_name"), c(cellline_name, drug_name))

  # tab_plots
  dt_average_plot <- dt_average[selected_combination, on = c(cellline_name, drug_name), nomatch = NULL]
  dt_metrics_plot <- dt_metrics[selected_combination, on = c(cellline_name, drug_name), nomatch = NULL]

  if (NROW(dt_metrics_plot) > 0) {

    if (NROW(stats::na.omit(dt_metrics_plot)) > 0 && NROW(unique(dt_average_plot[[conc]])) > 1) {

      min_conc <- min(dt_average_plot[get(conc) != 0][[conc]])
      max_conc <- max(dt_average_plot[[conc]])
      sampled_conc <- create_log_seq(min_conc, max_conc, 50)
      fitted_curve_sampled <- gDRutils::predict_efficacy_from_conc(sampled_conc,
                                                                   dt_metrics_plot$x_inf,
                                                                   dt_metrics_plot$x_0,
                                                                   dt_metrics_plot$ec50,
                                                                   dt_metrics_plot$h)
      dt_reconstructed_fit <- data.table::data.table(Concentration = sampled_conc,
                                                     x = fitted_curve_sampled)

    } else {
      dt_reconstructed_fit <- data.table::data.table(Concentration = numeric(),
                                                     x = numeric())
    }

    # set min and max values for y
    ymin <- min(c(0, min(dt_average_plot$x)))
    ymax <- max(c(1.2, max(dt_average_plot$x)))

    legend_title <- d_name
    plt_title <- cl_name

    # plot
    plt <-
      ggplot2::ggplot() +
      ggplot2::geom_hline(yintercept = c(0, 1), color = hline_color) +
      ggplot2::geom_line(
        data = dt_reconstructed_fit,
        ggplot2::aes(x = get(conc), y = x, color = "Fitted Curve", group = "Fitted Curve")) +
      ggplot2::geom_errorbar(
        data = dt_average_plot,
        ggplot2::aes(x = get(conc), y = x,  ymin = x - x_std, ymax = x + x_std,
                     color = "Errors Bar"),
        width = 0.1, position = ggplot2::position_dodge(0.1)) +
      ggplot2::geom_point(
        data = dt_average_plot,
        ggplot2::aes(x = get(conc), y = x, color = "Averaged Data")) +
      ggplot2::geom_line(
        data = dt_average_plot,
        ggplot2::aes(x = get(conc), y = x, color = "Averaged Data", group = "Averaged Data"),
        linetype = "longdash") +
      ggplot2::scale_x_log10(oob = scales::squish_infinite) +
      ggplot2::scale_y_continuous(lim = c(ymin, ymax)) +
      ggplot2::xlab(bquote(.(conc) ~ "[" ~ mu * M ~ "]")) +
      ggplot2::ylab(normalization_type) +
      ggplot2::ggtitle(plt_title) +
      ggplot2::labs(color = legend_title) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
                     axis.text.y = ggplot2::element_text(size = 8),
                     plot.title = ggplot2::element_text(size = 10),
                     legend.title = ggplot2::element_text(size = 10),
                     panel.grid.minor = ggplot2::element_blank(),
                     aspect.ratio = 1) +
      ggplot2::scale_color_manual(values = c("Errors Bar" = "black",
                                             "Averaged Data" = "black",
                                             "Fitted Curve" = "red"))
  } else {
    plt <-
      ggplot2::ggplot() +
      ggplot2::labs(x = bquote(.(conc) ~ "[" ~ mu * M ~ "]"),
                    y = normalization_type,
                    title = cl_name) +
      ggplot2::theme_bw() +
      ggplot2::theme(aspect.ratio = 1)
  }

  return(plt)
}

#' Plot panel with drug response curves for single-agent data to control quality of the data
#'
#' @inheritParams plot_dose_response_sa_qc
#' @param d_names character vector with drug names to be plotted (Drug Name);
#'    if NULL - all available drugs will be plotted
#'
#' @return \code{ggplot} object with panel with plots of dose-response curves for selected
#'    cell line by drugs (observed and fitted values)
#'
#' @keywords QC_plot
#' @examples
#' mae <- gDRutils::get_synthetic_data("small")
#' se <- mae[[1]]
#'
#' dt_metrics <- gDRutils::convert_se_assay_to_dt(se, "Metrics")
#' dt_average <- gDRutils::convert_se_assay_to_dt(se, "Averaged")
#' cl_name <- dt_metrics[["CellLineName"]][1]
#' d_names <- unique(dt_metrics[["DrugName"]])[1:5]
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

  checkmate::assert_data_table(dt_metrics)
  checkmate::assert_data_table(dt_average)
  checkmate::assert_string(cl_name)
  checkmate::assert_character(d_names, null.ok = TRUE)
  checkmate::assert_choice(normalization_type, choices = c("GR", "RV"))
  checkmate::assert_string(fit_source, null.ok = TRUE)

  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  clid <- gDRutils::get_env_identifiers("cellline")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  gnumber <- gDRutils::get_env_identifiers("drug")
  conc <- gDRutils::get_env_identifiers("concentration")

  checkmate::assert_choice(cl_name, choices = dt_metrics[[cellline_name]])
  checkmate::assert_choice(cl_name, choices = dt_average[[cellline_name]])
  hline_color <-
    gDRutils::get_settings_from_json("HLINE_COLOR",
                                     system.file(package = "gDRplots", "settings.json"))

  # filter data for normalization_type and fit_source
  filter_expr <- substitute(normalization_type == norm_type & fit_source == fit_src,
                            list(norm_type = normalization_type, fit_src = fit_source))
  dt_metrics <- dt_metrics[eval(filter_expr)]
  dt_average <- dt_average[eval(filter_expr)]

  available_drugs <- unique(dt_metrics[get(cellline_name) %in% cl_name, ][[drug_name]])
  if (is.null(d_names) || all(!d_names %in% available_drugs)) {
    d_names  <- available_drugs
  } else if (!all(d_names %in% available_drugs)) {
    d_names <- d_names[d_names %in% available_drugs]
  }

  # filter data for min required data
  dt_metrics <-
    dt_metrics[, c(drug_name, gnumber, cellline_name, clid, "x_inf", "x_0", "ec50", "h"), with = FALSE]
  dt_average <-
    dt_average[, c(drug_name, gnumber, cellline_name, clid, conc, "x", "x_std"), with = FALSE]

  selected_combination <- data.table::data.table(cellline_name = cl_name,
                                                 drug_name = d_names)
  data.table::setnames(selected_combination, c("cellline_name", "drug_name"), c(cellline_name, drug_name))

  # tab_plots
  dt_average_plot <- dt_average[selected_combination, on = c(cellline_name, drug_name)]
  dt_metrics_plot <- dt_metrics[selected_combination, on = c(cellline_name, drug_name)]

  dt_reconstructed_fit <- data.table::data.table(drug_name = character(),
                                                 conc = numeric(),
                                                 x = numeric())

  for (d_name in d_names) {
    dt_met_plot <- dt_metrics_plot[get(drug_name) == d_name, ]
    dt_avg_plot <- dt_average_plot[get(drug_name) == d_name, ]

    if (NROW(stats::na.omit(dt_met_plot)) > 0 && NROW(unique(dt_avg_plot[[conc]])) > 1) {
      min_conc <- min(dt_avg_plot[get(conc) != 0][[conc]])
      max_conc <- max(dt_avg_plot[[conc]])
      sampled_conc <- create_log_seq(min_conc, max_conc, 50)

      fitted_curve_sampled <- gDRutils::predict_efficacy_from_conc(sampled_conc,
                                                                   dt_met_plot$x_inf,
                                                                   dt_met_plot$x_0,
                                                                   dt_met_plot$ec50,
                                                                   dt_met_plot$h)
      dt_fit <- data.table::data.table(
        drug_name = d_name,
        conc = sampled_conc,
        x = fitted_curve_sampled
      )

      dt_reconstructed_fit <- rbind(dt_reconstructed_fit, dt_fit)
    }
  }
  data.table::setnames(dt_reconstructed_fit, c("drug_name", "conc"), c(drug_name, conc))

  # panel title
  cl_clid <- unique(dt_metrics[get(cellline_name) == cl_name, ][[clid]])
  panel_title <- sprintf("%s (%s)", cl_name, cl_clid)

  # set min and max values for y
  ymin <- min(c(0, min(dt_average_plot$x)))
  ymax <- max(c(1.2, max(dt_average_plot$x)))

  plt <-
    ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = c(0, 1), color = hline_color) +
    ggplot2::geom_line(
      data = dt_reconstructed_fit,
      ggplot2::aes(x = get(conc), y = x, color = "Fitted Curve")) +
    ggplot2::geom_errorbar(
      data = dt_average_plot,
      ggplot2::aes(x = get(conc), y = x,  ymin = x - x_std, ymax = x + x_std, color = "Errors Bar"),
      width = 0.1, position = ggplot2::position_dodge(0.1)) +
    ggplot2::geom_point(
      data = dt_average_plot,
      ggplot2::aes(x = get(conc), y = x, color = "Averaged Data"),
      shape = 20) +
    ggplot2::scale_x_log10(oob = scales::squish_infinite) +
    ggplot2::scale_y_continuous(lim = c(ymin, ymax)) +
    ggplot2::xlab(bquote(.(conc) ~ "[" ~ mu * M ~ "]")) +
    ggplot2::ylab(normalization_type) +
    ggplot2::ggtitle(panel_title) +
    ggplot2::labs(color = "") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
                   axis.text.y = ggplot2::element_text(size = 8),
                   plot.title = ggplot2::element_text(size = 10),
                   panel.grid.minor = ggplot2::element_blank(),
                   aspect.ratio = 1) +
    ggplot2::scale_color_manual(values = c("Errors Bar" = "black",
                                           "Averaged Data" = "black",
                                           "Fitted Curve" = "red")) +
    ggplot2::facet_wrap(~get(drug_name)) +
    ggplot2::theme(strip.background = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(face = "bold", hjust = 0),
                   legend.position = "top",
                   plot.title = ggplot2::element_text(size = 12, hjust = 0.5))

  return(plt)
}
