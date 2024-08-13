#' Plot drug response curves for combo data
#'
#' @param dt_average data.table representing data from the \code{Averaged} assay,
#'    outputted by \code{gDRutils::convert_se_assay_to_dt(se, "Averaged")}
#' @param drug1_name string with drug name to be plotted (identifiers \code{DrugName})
#' @param drug2_name string with co-drug name to be plotted (identifiers \code{DrugName_2})
#' @param cl_name string with cell line name to be plotted (identifiers \code{CellLineName})
#' @param normalization_type string with normalization_types to be selected
#'                           one of: "GR" ("GRvalue") or "RV" ("RelativeViability")
#' @param colormap character vector with colors for \code{group_names} - name or hex value
#' @param split_by_conc split_by_conc logical flag indicating whether curves 
#'    for \code{Concentration_2} should be plotted on a single plot or separately
#'
#' @return plot with dose-response curves for combo data
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
#'                          colormap = c("orange", "darkred"))                     
#'                       
#' @export
plot_dose_response_combo <- function(dt_average,
                                     drug1_name,
                                     drug2_name,
                                     cl_name,
                                     normalization_type = "GR",
                                     colormap = NULL, 
                                     split_by_conc = FALSE) {
  
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  clid <- gDRutils::get_env_identifiers("cellline")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  drug_name_2 <- gDRutils::get_env_identifiers("drug_name2")
  conc <- gDRutils::get_env_identifiers("concentration")
  conc_2 <- gDRutils::get_env_identifiers("concentration2")
  
  checkmate::assert_data_table(dt_average)
  checkmate::assert_choice(drug1_name, choices = unique(dt_average[[drug_name]]))
  checkmate::assert_choice(drug2_name, choices = unique(dt_average[[drug_name_2]]))
  checkmate::assert_string(cl_name)
  checkmate::assert_choice(cl_name, choices = unique(dt_average[[cellline_name]]))
  checkmate::assert_choice(normalization_type, choices = c("GR", "RV"))
  checkmate::assert_character(colormap, null.ok = TRUE)
  checkmate::assert_flag(split_by_conc)
  
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
  
  dt_avg[[conc_2]] <- factor(dt_avg[[conc_2]])
  
  # handle conc = 0
  min_conc <- min(dt_avg[dt_avg[[conc]] > 0, ][[conc]])
  dt_avg[[conc]][dt_avg[[conc]] == 0] <- min_conc / 100
  
  # colors
  ls_conc_2 <- unique(dt_avg[[conc_2]])
  if (is.null(colormap) || !all(vapply(colormap, is_valid_color, logical(1)))) {
    number_of_color <- NROW(ls_conc_2)
    colormap <-
      rev(colorspace::sequential_hcl(number_of_color + 1, palette = "viridis")[seq_along(ls_conc_2)])
  } else if (NROW(colormap) != NROW(ls_conc_2)) {
    colormap <- grDevices::colorRampPalette(colormap)(NROW(ls_conc_2))
  }
  names(colormap) <- levels(ls_conc_2)
  
  # set min and max values for y
  ymin <- min(c(0, min(dt_avg$x)))
  ymax <- max(c(1.2, max(dt_avg$x)))
  
  # final plot
  plt <-
    ggplot2::ggplot(dt_avg,
                    ggplot2::aes(x = get(conc), y = x, color = get(conc_2), group = get(conc_2))) +
    ggplot2::geom_hline(yintercept = 0, color = "#A9A9A9") +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::scale_y_continuous(lim = c(ymin, ymax)) +
    ggplot2::scale_x_continuous(trans = "log10") +
    ggplot2::scale_color_manual(values = colormap, 
                                labels = sprintf("%.4f", as.numeric(levels(ls_conc_2)))) +
    ggplot2::xlab(bquote(.(drug1_name) ~ "[" ~ mu * M ~ "]")) +
    ggplot2::ylab(sprintf("log10(%s)", normalization_type)) +
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
#' @return panel with plot with dose-response curves for selected cell line by drugs
#'
#' @keywords QC_plot
#' @examples
#' mae <- gDRutils::get_synthetic_data("combo_matrix")
#' se <- mae[[gDRutils::get_supported_experiments("combo")]]
#' dt_average <- gDRutils::convert_se_assay_to_dt(se, "Averaged")
#' 
#' cl_name <- "cellline_IB"
#' 
#' plot_dose_response_combo_qc_panel(dt_average = dt_average,
#'                                   cl_name = cl_name)
#' 
#' d_names <- c("drug_001", "drug_002")
#' plot_dose_response_combo_qc_panel(dt_average = dt_average,
#'                                   cl_name = cl_name,
#'                                   d_names = d_names)
#' 
#' @export
plot_dose_response_combo_qc_panel <- function(dt_average,
                                              cl_name,
                                              d_names = NULL,
                                              normalization_type = "GR",
                                              colormap = NULL) {
  
  
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  clid <- gDRutils::get_env_identifiers("cellline")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  drug_name_2 <- gDRutils::get_env_identifiers("drug_name2")
  conc <- gDRutils::get_env_identifiers("concentration")
  conc_2 <- gDRutils::get_env_identifiers("concentration2")
  
  checkmate::assert_data_table(dt_average)
  checkmate::assert_string(cl_name)
  checkmate::assert_choice(cl_name, choices = unique(dt_average[[cellline_name]]))
  checkmate::assert_character(d_names, null.ok = TRUE)
  checkmate::assert_choice(normalization_type, choices = c("GR", "RV"))
  checkmate::assert_character(colormap, null.ok = TRUE)
  
  available_drugs <- unique(dt_average[[drug_name]])
  if (is.null(d_names) || all(!d_names %in% available_drugs)) {
    d_names  <- available_drugs
  } else if (!all(d_names %in% available_drugs)) {
    d_names <- drug_name[drug_name %in% available_drugs]
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
  
  dt_avg[[conc_2]] <- factor(dt_avg[[conc_2]])
  
  # handle conc = 0
  min_conc <- min(dt_avg[dt_avg[[conc]] > 0, ][[conc]])
  dt_avg[[conc]][dt_avg[[conc]] == 0] <- min_conc / 100
  
  # colors
  ls_conc_2 <- unique(dt_avg[[conc_2]])
  if (is.null(colormap) || !all(vapply(colormap, is_valid_color, logical(1)))) {
    number_of_color <- NROW(ls_conc_2)
    colormap <-
      rev(colorspace::sequential_hcl(number_of_color + 1, palette = "viridis")[seq_along(ls_conc_2)])
  } else if (NROW(colormap) != NROW(ls_conc_2)) {
    colormap <- grDevices::colorRampPalette(colormap)(NROW(ls_conc_2))
  }
  names(colormap) <- levels(ls_conc_2)
  
  # set min and max values for y
  ymin <- min(c(0, min(dt_avg$x)))
  ymax <- max(c(1.2, max(dt_avg$x)))
  
  # final plot
  plt <-
    ggplot2::ggplot(dt_avg,
                    ggplot2::aes(x = get(conc), y = x, color = get(conc_2), group = get(conc_2))) +
    ggplot2::geom_hline(yintercept = 0, color = "#A9A9A9") +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::scale_y_continuous(lim = c(ymin, ymax)) +
    ggplot2::scale_x_continuous(trans = "log10") +
    ggplot2::scale_color_manual(values = colormap, 
                                labels = sprintf("%.4f", as.numeric(levels(ls_conc_2)))) +
    ggplot2::xlab(bquote(~ "Concentration of Drug [" ~ mu * M ~ "]")) +
    ggplot2::ylab(sprintf("log10(%s)", normalization_type)) +
    ggplot2::ggtitle(panel_title) +
    ggplot2::labs(color = bquote(~ "Codrug [" ~ mu * M ~ "]")) +
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
      plot.title = ggplot2::element_text(size = 12, hjust = 0.5, margin = ggplot2::margin()),
      panel.spacing = ggplot2::unit(0.5, "pt"))
  
  return(plt)
}
