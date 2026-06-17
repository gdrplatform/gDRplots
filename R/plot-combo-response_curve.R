#' Plot drug response curves for combo data
#'
#' @param dt_average data.table representing data from the \code{Averaged} assay,
#'    outputted by \code{gDRutils::convert_se_assay_to_dt(se, "Averaged")}
#'    and combo \code{SummarizedExperiment}
#' @param drug1_name string with drug name to be plotted (identifiers \code{DrugName})
#' @param drug2_name string with co-drug name to be plotted (identifiers \code{DrugName_2})
#' @param cl_name string with cell line name to be plotted (identifiers \code{CellLineName})
#' @param normalization_type string with normalization_types to be selected
#'                           one of: "GR" ("GRvalue") or "RV" ("RelativeViability")
#' @param colors_vec character vector with colors for \code{group_names} - name or hex value
#'    note that the first color will be assigned to the min value of \code{Concentration_2},
#'    and the last one - to the max of \code{Concentration_2}; the default is the orange-black palette
#' @param split_by_conc split_by_conc logical flag indicating whether curves
#'    for \code{Concentration_2} should be plotted on a single plot or separately
#'
#' @return \code{ggplot} object containing plot with dose-response curves for combo data
#'
#' @keywords combo_plots
#' @examples
#' mae <- gDRutils::get_synthetic_data("combo_matrix")
#' se <- mae[[gDRutils::get_supported_experiments("combo")]]
#' dt_average <- gDRutils::convert_se_assay_to_dt(se, "Averaged")
#'
#' cl_name <- "cellline_BC"
#' drug1_name <- "drug_011"
#' drug2_name <- "drug_021"
#'
#' plot_dose_response_combo(dt_average = dt_average,
#'                          drug1_name = drug1_name,
#'                          drug2_name = drug2_name,
#'                          cl_name = cl_name)
#'
#' plot_dose_response_combo(dt_average = dt_average,
#'                          drug1_name = drug1_name,
#'                          drug2_name = drug2_name,
#'                          cl_name = cl_name,
#'                          split_by_conc = TRUE)
#'
#' plot_dose_response_combo(dt_average = dt_average,
#'                          drug1_name = drug1_name,
#'                          drug2_name = drug2_name,
#'                          cl_name = cl_name,
#'                          normalization_type = "RV",
#'                          colors_vec = c("lightblue", "darkblue"))
#'
#' @export
plot_dose_response_combo <- function(dt_average,
                                     drug1_name,
                                     drug2_name,
                                     cl_name,
                                     normalization_type = "GR",
                                     colors_vec = NULL,
                                     split_by_conc = FALSE) {

  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  clid <- gDRutils::get_env_identifiers("cellline")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  drug_name_2 <- gDRutils::get_env_identifiers("drug_name2")
  conc <- gDRutils::get_env_identifiers("concentration")
  conc_2 <- gDRutils::get_env_identifiers("concentration2")
  zero_conc_scaling_factor <-
    gDRutils::get_settings_from_json("ZERO_CONC_SCALING_FACTOR",
                                     system.file(package = "gDRplots", "settings.json"))

  checkmate::assert_data_table(dt_average)
  checkmate::assert_choice(drug1_name, choices = unique(dt_average[[drug_name]]))
  checkmate::assert_choice(drug2_name, choices = unique(dt_average[[drug_name_2]]))
  checkmate::assert_string(cl_name)
  checkmate::assert_choice(cl_name, choices = unique(dt_average[[cellline_name]]))
  checkmate::assert_choice(normalization_type, choices = c("GR", "RV"))
  checkmate::assert_character(colors_vec, null.ok = TRUE)
  checkmate::assert_flag(split_by_conc)
  hline_color <-
    gDRutils::get_settings_from_json("HLINE_COLOR",
                                     system.file(package = "gDRplots", "settings.json"))

  # check input data
  drugs_combination <-
    unique(dt_average[get(cellline_name) == cl_name, .SD, .SDcols = c(cellline_name, drug_name, drug_name_2)])
  stopifnot("combination of drugs and cell line does not exist" =
              any(drug2_name %in% drugs_combination[[drug_name_2]],
                  drug1_name %in% drugs_combination[[drug_name]]))

  # plt title
  cl_clid <- unique(dt_average[get(cellline_name) == cl_name, ][[clid]])
  plt_title <- sprintf("%s (%s)", cl_name, cl_clid)

  # filter data for normalization type
  filter_expr <- substitute(normalization_type == norm_type, list(norm_type = normalization_type))
  dt_avg <- dt_average[eval(filter_expr)]

  # and selected cell line
  required_cols <- c(cellline_name, drug_name, drug_name_2, conc, conc_2, "x")
  dt_avg <- dt_avg[get(cellline_name) == cl_name, ][, .SD, .SDcols = required_cols]

  # filter data for combination cell line (drug x drug2)
  selected_combination <-
    drugs_combination[get(drug_name) == drug1_name & get(drug_name_2) == drug2_name, ]

  dt_avg <- dt_avg[selected_combination, on = c(cellline_name, drug_name, drug_name_2)]

  # Ensure the drug with more dose levels is on the x-axis
  n_conc_1 <- sum(unique(dt_avg[[conc]]) > 0)
  n_conc_2 <- sum(unique(dt_avg[[conc_2]]) > 0)
  if (n_conc_2 > n_conc_1) {
    tmp <- dt_avg[[conc]]
    data.table::set(dt_avg, j = conc, value = dt_avg[[conc_2]])
    data.table::set(dt_avg, j = conc_2, value = tmp)
    tmp_name <- drug1_name
    drug1_name <- drug2_name
    drug2_name <- tmp_name
  }

  dt_avg[[conc_2]] <- factor(dt_avg[[conc_2]],
                             levels = sort(unique(dt_avg[[conc_2]])),
                             labels = round_to_unique_string(sort(unique(dt_avg[[conc_2]]))))

  # handle conc = 0
  min_conc <- min(dt_avg[dt_avg[[conc]] > 0, ][[conc]])
  dt_avg[[conc]][dt_avg[[conc]] == 0] <- min_conc / zero_conc_scaling_factor

  # colors
  ls_conc_2 <- unique(dt_avg[[conc_2]])
  if (is.null(colors_vec) || !all(vapply(colors_vec, is_valid_color, logical(1)))) {
    colormap <- .get_combo_curves_colors(ls_conc_2)
  } else {
    colormap <- grDevices::colorRampPalette(colors_vec)(NROW(ls_conc_2))
    names(colormap) <- levels(ls_conc_2)
  }

  # set min and max values for y
  ymin <- min(c(0, min(dt_avg$x)))
  ymax <- max(c(1.2, max(dt_avg$x)))

  # final plot
  plt <-
    ggplot2::ggplot(dt_avg,
                    ggplot2::aes(x = get(conc), y = x, color = get(conc_2), group = get(conc_2))) +
    ggplot2::geom_hline(yintercept = c(0, 1), color = hline_color) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::scale_y_continuous(lim = c(ymin, ymax)) +
    ggplot2::scale_x_log10(oob = scales::squish_infinite) +
    ggplot2::scale_color_manual(values = colormap,
                                labels = levels(ls_conc_2)) +
    ggplot2::xlab(bquote(.(drug1_name) ~ "[" ~ mu * M ~ "]")) +
    ggplot2::ylab(normalization_type) +
    ggplot2::ggtitle(plt_title) +
    ggplot2::labs(color = bquote(.(drug2_name) ~ "[" ~ mu * M ~ "]")) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 8),
      plot.title = ggplot2::element_text(size = 10),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "left",
      aspect.ratio = 1)

  # split
  if (split_by_conc) {
    lbls <- sprintf("%.4f", as.numeric(levels(ls_conc_2)))
    names(lbls) <- levels(ls_conc_2)

    plt <- plt +
      ggplot2::facet_wrap(~get(conc_2), labeller = ggplot2::as_labeller(lbls)) +
      ggplot2::theme(
        legend.position = "none",
        strip.text = ggplot2::element_text(face = "bold"))
  }

  return(plt)
}

#' Plot panel with drug response curves for single-agent data to control quality of the data
#'
#' @inheritParams plot_dose_response_combo
#' @param d_names character vector with drug names to be plotted (Drug Name);
#'    if NULL - all available drugs will be plotted
#'
#' @return \code{ggplot} object containing panel with plot with dose-response curves
#'    for selected cell line by drugs
#'
#' @keywords combo_plots
#' @examples
#' mae <- gDRutils::get_synthetic_data("combo_matrix")
#' se <- mae[[gDRutils::get_supported_experiments("combo")]]
#' dt_average <- gDRutils::convert_se_assay_to_dt(se, "Averaged")
#'
#' cl_name <- "cellline_IB"
#'
#' plot_dose_response_combo_panel(dt_average = dt_average,
#'                                cl_name = cl_name)
#'
#' d_names <- c("drug_001", "drug_002")
#' plot_dose_response_combo_panel(dt_average = dt_average,
#'                                cl_name = cl_name,
#'                                d_names = d_names)
#'
#' @export
plot_dose_response_combo_panel <- function(dt_average,
                                           cl_name,
                                           d_names = NULL,
                                           normalization_type = "GR",
                                           colors_vec = NULL) {

  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  clid <- gDRutils::get_env_identifiers("cellline")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  drug_name_2 <- gDRutils::get_env_identifiers("drug_name2")
  conc <- gDRutils::get_env_identifiers("concentration")
  conc_2 <- gDRutils::get_env_identifiers("concentration2")
  zero_conc_scaling_factor <-
    gDRutils::get_settings_from_json("ZERO_CONC_SCALING_FACTOR",
                                     system.file(package = "gDRplots", "settings.json"))
  hline_color <-
    gDRutils::get_settings_from_json("HLINE_COLOR",
                                     system.file(package = "gDRplots", "settings.json"))

  checkmate::assert_data_table(dt_average)
  checkmate::assert_string(cl_name)
  checkmate::assert_choice(cl_name, choices = unique(dt_average[[cellline_name]]))
  checkmate::assert_character(d_names, null.ok = TRUE)
  checkmate::assert_choice(normalization_type, choices = c("GR", "RV"))
  checkmate::assert_character(colors_vec, null.ok = TRUE)


  available_drugs <- unique(dt_average[[drug_name]])
  if (is.null(d_names) || all(!d_names %in% available_drugs)) {
    d_names  <- available_drugs
  } else if (!all(d_names %in% available_drugs)) {
    d_names <- d_names[d_names %in% available_drugs]
  }

  # check input data
  selected_combination <-
    unique(dt_average[get(cellline_name) == cl_name, .SD, .SDcols = c(cellline_name, drug_name, drug_name_2)])
  selected_combination <- selected_combination[get(drug_name) %in% d_names, ]
  stopifnot("combination of drugs and cell line does not exist" =
              any(d_names %in% selected_combination[[drug_name]]))

  # panel title
  cl_clid <- unique(dt_average[get(cellline_name) == cl_name, ][[clid]])
  panel_title <- sprintf("%s (%s)", cl_name, cl_clid)

  # filter data for normalization type
  filter_expr <- substitute(normalization_type == norm_type, list(norm_type = normalization_type))
  dt_avg <- dt_average[eval(filter_expr)]

  # and selected cell line
  required_cols <- c(cellline_name, drug_name, drug_name_2, conc, conc_2, "x")
  dt_avg <- dt_avg[, .SD, .SDcols = required_cols]

  # filter data for combination cell line (drug x drug2)
  dt_avg <- dt_avg[selected_combination, on = c(cellline_name, drug_name, drug_name_2)]

  # Ensure the drug with more dose levels is on the x-axis for each pair
  pairs <- unique(dt_avg[, .SD, .SDcols = c(drug_name, drug_name_2)])
  for (i in seq_len(NROW(pairs))) {
    d1 <- pairs[[drug_name]][i]
    d2 <- pairs[[drug_name_2]][i]
    idx <- which(dt_avg[[drug_name]] == d1 & dt_avg[[drug_name_2]] == d2)
    n_conc_1 <- sum(unique(dt_avg[[conc]][idx]) > 0)
    n_conc_2 <- sum(unique(dt_avg[[conc_2]][idx]) > 0)
    if (n_conc_2 > n_conc_1) {
      old_conc <- dt_avg[[conc]][idx]
      old_conc_2 <- dt_avg[[conc_2]][idx]
      old_d <- dt_avg[[drug_name]][idx]
      old_d2 <- dt_avg[[drug_name_2]][idx]
      data.table::set(dt_avg, idx, conc, old_conc_2)
      data.table::set(dt_avg, idx, conc_2, old_conc)
      data.table::set(dt_avg, idx, drug_name, old_d2)
      data.table::set(dt_avg, idx, drug_name_2, old_d)
    }
  }

  dt_avg[[conc_2]] <- factor(dt_avg[[conc_2]],
                             levels = sort(unique(dt_avg[[conc_2]])),
                             labels = round_to_unique_string(sort(unique(dt_avg[[conc_2]]))))

  # handle conc = 0
  min_conc <- min(dt_avg[dt_avg[[conc]] > 0, ][[conc]])
  dt_avg[[conc]][dt_avg[[conc]] == 0] <- min_conc / zero_conc_scaling_factor

  # colors
  ls_conc_2 <- unique(dt_avg[[conc_2]])
  if (is.null(colors_vec) || !all(vapply(colors_vec, is_valid_color, logical(1)))) {
    colormap <- .get_combo_curves_colors(ls_conc_2)
  } else if (NROW(colors_vec) != NROW(ls_conc_2)) {
    colormap <- grDevices::colorRampPalette(colors_vec)(NROW(ls_conc_2))
    names(colormap) <- levels(ls_conc_2)
  }

  # set min and max values for y
  ymin <- min(c(0, min(dt_avg$x)))
  ymax <- max(c(1.2, max(dt_avg$x)))

  # final plot
  plt <-
    ggplot2::ggplot(dt_avg,
                    ggplot2::aes(x = get(conc), y = x, color = get(conc_2), group = get(conc_2))) +
    ggplot2::geom_hline(yintercept = c(0, 1), color = hline_color) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::scale_y_continuous(lim = c(ymin, ymax)) +
    ggplot2::scale_x_log10(oob = scales::squish_infinite) +
    ggplot2::scale_color_manual(values = colormap,
                                labels = levels(ls_conc_2)) +
    ggplot2::xlab(bquote(~ "Concentration [" ~ mu * M ~ "]")) +
    ggplot2::ylab(normalization_type) +
    ggplot2::ggtitle(panel_title) +
    ggplot2::labs(color = bquote(~ "Co-treatment [" ~ mu * M ~ "]")) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 8, angle = 45, vjust = 1, hjust = 1, margin = ggplot2::margin()),
      axis.text.y = ggplot2::element_text(size = 8, margin = ggplot2::margin()),
      plot.title = ggplot2::element_text(size = 10, margin = ggplot2::margin()),
      panel.grid.minor = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 8, margin = ggplot2::margin()),
      legend.title = ggplot2::element_text(size = 8, margin = ggplot2::margin()),
      aspect.ratio = 1) +
    ggplot2::facet_wrap(~get(drug_name) + get(drug_name_2)) +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 8, face = "bold", hjust = 0, margin = ggplot2::margin()),
      legend.position = "left",
      plot.title = ggplot2::element_text(size = 12, hjust = 0.5, margin = ggplot2::margin()))

  return(plt)
}

#' Get color palette for the dose response curves for combination data
#'
#' @param ls_conc_2 factor vector with values for \code{Concentration_2}
#'
#' @return gDR palette for Concentration_2 given in \code{ls_conc_2}
#'
#' @keywords internal
#' @examples
#' \dontrun{
#' ls_conc <- factor(c("0.001", "0.01", "1"))
#' .get_combo_curves_colors(ls_conc)
#' }
.get_combo_curves_colors <- function(ls_conc_2) {
  checkmate::assert_factor(ls_conc_2)

  # order iso level
  ls_conc_2 <- ls_conc_2[order(as.numeric(ls_conc_2))]

  conc_2_colors <-
    grDevices::colorRampPalette(
      gDRutils::get_settings_from_json("COMBO_CURVES_PALETTE",
                                       system.file(package = "gDRplots", "settings.json"))
    )(2 * NROW(ls_conc_2))[2 * seq_along(ls_conc_2)]
  names(conc_2_colors) <- levels(ls_conc_2)

  conc_2_colors
}
