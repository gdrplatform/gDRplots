#' Plot heatmaps of fitted values for combination metrics data
#'
#' @param dt_excess data.table representing data from the \code{excess} assay,
#'    outputted by \code{gDRutils::convert_se_assay_to_dt(se, "excess")}
#'    and combo \code{SummarizedExperiment}
#' @param dt_isobolograms data.table representing data from the \code{isobolograms} assay,
#'    outputted by \code{gDRutils::convert_se_assay_to_dt(se, "isobolograms")}
#'    and combo \code{SummarizedExperiment}
#' @param drug1_name string with drug name to be plotted (identifiers \code{DrugName})
#' @param drug2_name string with co-drug name to be plotted (identifiers \code{DrugName_2})
#' @param cl_name string with cell line to be plotted (identifiers \code{CellLineName})
#' @param metric string name of combo metric;
#'   one of: "smooth"("Smooth GR" or "Smooth RV" - respectively depending on \code{normalization_type}), 
#'   "hsa_excess"("Bliss Excess GR" or "Bliss Excess RV")
#'   or "bliss_excess" ("Bliss Excess GR" or "Bliss Excess RV")
#' @param normalization_type string with normalization_types to be selected
#'                           one of: "GR" ("GRvalue") or "RV" ("RelativeViability")
#' @param iso_levels character vector with isobologram levels to be selected
#' @param colors_vec_smooth character vector of colors (valid names or hex codes) used in the heatmap
#'    for smooth values; the default is the dark purple-light grey palette
#' @param colors_vec_excess character vector of colors (valid name or hex codes) used in the heatmap
#'    for excess values; the default is the blue-light grey-red color scale
#' @param limit numeric vector of length two providing limits of the scale. 
#'    Use NA to refer to the existing minimum or maximum
#' @param no_breaks numeric number of breaks on scale
#' @param swap_axes logical flag indicating whether to swap the axes with drugs of the heatmap
#' @param show_values logical flag indicating whether to show values of the metric on the heatmap
#'
#' @return \code{ggplot} object containing heatmap with the values for the selected combo metric 
#'    for selected drugs and cell line with selected isoline (when chosen)
#' 
#' @keywords combo_plots
#' @examples
#' cl_name <- "cellline_BC"
#' drug1_name <- "drug_001"
#' drug2_name <- "drug_026"
#' 
#' mae <- gDRutils::get_synthetic_data("combo_matrix")
#' se <- mae[[gDRutils::get_supported_experiments("combo")]]
#' dt_excess <- gDRutils::convert_se_assay_to_dt(se, "excess")
#' dt_isobolograms <- gDRutils::convert_se_assay_to_dt(se, "isobolograms")
#' 
#' heatmap_combo_metrics(dt_excess,
#'                       dt_isobolograms,
#'                       drug1_name, drug2_name,
#'                       cl_name,
#'                       normalization_type = "GR")
#'                       
#' heatmap_combo_metrics(dt_excess,
#'                       dt_isobolograms,
#'                       drug1_name, drug2_name,
#'                       cl_name,
#'                       normalization_type = "GR",
#'                       swap_axes = TRUE)
#'         
#' heatmap_combo_metrics(dt_excess,
#'                       dt_isobolograms,
#'                       drug1_name, drug2_name,
#'                       cl_name,
#'                       metric = "hsa_excess",
#'                       normalization_type = "GR",
#'                       limit = c(NA, 0.1))
#'                       
#' cl_name <- "cellline_JE"
#' drug1_name <- "drug_011"
#' drug2_name <- "drug_026"
#' 
#' heatmap_combo_metrics(dt_excess,
#'                       dt_isobolograms,
#'                       drug1_name, drug2_name,
#'                       cl_name,
#'                       metric = "bliss_excess",
#'                       normalization_type = "RV",
#'                       iso_levels = "0.5")
#'                       
#' heatmap_combo_metrics(dt_excess,
#'                       dt_isobolograms,
#'                       drug1_name, drug2_name,
#'                       cl_name,
#'                       metric = "bliss_excess",
#'                       normalization_type = "RV",
#'                       iso_levels = NULL)
#'
#' @export
heatmap_combo_metrics <- function(
    dt_excess,
    dt_isobolograms = NULL,
    drug1_name,
    drug2_name,
    cl_name,
    metric = "smooth",
    normalization_type = "GR",
    iso_levels =  c("0.25", "0.5", "0.75"),
    colors_vec_smooth = NULL,
    colors_vec_excess = NULL,
    limit = NULL,
    no_breaks = 50,
    swap_axes = FALSE,
    show_values = FALSE) {
  
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  clid <- gDRutils::get_env_identifiers("cellline")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  drug_name_2 <- gDRutils::get_env_identifiers("drug_name2")
  conc <- gDRutils::get_env_identifiers("concentration")
  conc_2 <- gDRutils::get_env_identifiers("concentration2")
  duration <- gDRutils::get_env_identifiers("duration")
  
  checkmate::assert_data_table(dt_excess)
  checkmate::assert_string(drug1_name)
  checkmate::assert_choice(drug1_name, choices = dt_excess[[drug_name]])
  checkmate::assert_string(drug2_name)
  checkmate::assert_choice(drug2_name, choices = dt_excess[[drug_name_2]])
  checkmate::assert_string(cl_name)
  checkmate::assert_choice(cl_name, choices = dt_excess[[cellline_name]])
  checkmate::assert_choice(metric, names(gDRutils::get_combo_excess_field_names()))
  checkmate::assert_choice(normalization_type, choices = c("GR", "RV"))
  checkmate::assert_data_table(dt_isobolograms, null.ok = TRUE)
  if (!is.null(dt_isobolograms)) {
    checkmate::assert_character(iso_levels, null.ok = TRUE)
    checkmate::assert_numeric(as.numeric(iso_levels))
    checkmate::assert_names(names(dt_isobolograms), must.include = "iso_level")
  }
  checkmate::assert_numeric(limit, len = 2, null.ok = TRUE)
  checkmate::assert_int(no_breaks, lower = 2)
  checkmate::assert_flag(swap_axes)
  checkmate::assert_flag(show_values)
  
  # data filtering and processing
  filter_expr <- substitute(normalization_type == norm_type, list(norm_type = normalization_type))
  dt_excess <- dt_excess[eval(filter_expr)]
  if (!is.null(dt_isobolograms) && !is.null(iso_levels)) {
    dt_isobolograms <- dt_isobolograms[eval(filter_expr)]
    dt_isobolograms <- dt_isobolograms[get(cellline_name) == cl_name & get(drug_name) == drug1_name &
                                         get(drug_name_2) == drug2_name]
    dt_isobolograms <- dt_isobolograms[iso_level %in% iso_levels, ]
  }
  
  dt_excess <- dt_excess[get(cellline_name) == cl_name & get(drug_name) == drug1_name & get(drug_name_2) == drug2_name]
  
  # check if isolines are available and adjust plotting logic accordingly
  available_iso_lvl <- if (!is.null(dt_isobolograms) && !is.null(iso_levels)) {
    unique(dt_isobolograms[["iso_level"]])
  } else {
    NULL
  }
  
  iso_colors <- if (!is.null(available_iso_lvl)) {
    .get_iso_colors(available_iso_lvl)
  } else {
    NULL
  }
  
  # title
  main_title <- sprintf("%s (%s)",
                        cl_name,
                        unique(dt_excess[get(cellline_name) == cl_name][[clid]]))
  # legend
  legend_title_iso <- "Iso Levels"
  label_prefix <- if (normalization_type == "GR") {
    "GR"
  } else {
    "IC"
  }
  legend_lbl_iso <- sprintf("%s%s", label_prefix, 100 - 100 * as.numeric(available_iso_lvl))
  
  # prep hm color palette
  hm_color_palette_smooth <-
    if (is.null(colors_vec_smooth)  || !all(vapply(colors_vec_smooth, is_valid_color, logical(1)))) {
      .get_smooth_palette(no_breaks)
    } else {
      grDevices::colorRampPalette(colors_vec_smooth)(no_breaks + 1)
    }
  
  hm_color_palette_excess <- 
    if (is.null(colors_vec_excess) || !all(vapply(colors_vec_excess, is_valid_color, logical(1)))) {
      .get_excess_palette(no_breaks)
    } else {
      grDevices::colorRampPalette(colors_vec_excess)(no_breaks + 1)
    }
  
  # plot title
  plt_title <- sprintf("%s for %s, T=%sh",
                       gDRutils::prettify_flat_metrics(x = metric, human_readable = TRUE),
                       normalization_type,
                       unique(dt_excess[get(cellline_name) == cl_name][[duration]]))
  # plot data
  dt_ <- dt_excess[, c(conc, conc_2, metric), with = FALSE]
  
  x_axis_drug <- if (swap_axes) {
    drug1_name
  } else {
    drug2_name
  }
  
  y_axis_drug <- if (swap_axes) {
    drug2_name
  } else {
    drug1_name
  }
  
  x_axis_lab <- sprintf("%s [\U00B5M]", x_axis_drug)
  y_axis_lab <- sprintf("%s [\U00B5M]", y_axis_drug)
  
  if (!NROW(dt_) > 1) { # co-dilution input data is like: (conc = 0, conc_2 = 0, metric = 1)
    plt <- 
      ggplot2::ggplot() +
      ggplot2::labs(x = x_axis_lab,
                    y = y_axis_lab,
                    title = plt_title) +
      ggplot2::theme_bw() +
      ggplot2::theme(aspect.ratio = 1)
  } else {
    dt_[[metric]] <- pmin(1.1, dt_[[metric]])
    
    conc_y <- if (swap_axes) dt_[[conc_2]] else dt_[[conc]]
    conc_x <- if (swap_axes) dt_[[conc]] else dt_[[conc_2]]
    
    dt_$pos_y <- transform_log_conc(conc_y)
    dt_$pos_x <- transform_log_conc(conc_x)
    
    lbl_y <- sprintf("%.2g", gDRutils::round_concentration(sort(unique(conc_y))))
    mrk_y <- sort(unique(dt_$pos_y))
    
    lbl_x <- sprintf("%.2g", gDRutils::round_concentration(sort(unique(conc_x))))
    mrk_x <- sort(unique(dt_$pos_x))
    
    tile_height <- .get_tile_size(mrk_y)
    tile_width <- .get_tile_size(mrk_x)
    
    # prep hm color palette
    hm_color_palette <- if (metric == "smooth") {
      hm_color_palette_smooth
    } else {
      hm_color_palette_excess
    }
    
    # legend title
    legend_title_fill <- sprintf("%s %s",
                                 gDRutils::prettify_flat_metrics(x = metric, human_readable = TRUE),
                                 normalization_type)
    
    # prep limits
    limit_fill <- if (is.null(limit)) {
      prep_hm_limits(dt_[[metric]],   
                     metric = metric,
                     normalization_type = normalization_type,
                     symmetric = metric != "smooth")
    } else {
      limit
    }
    
    # base plot
    plt <-
      ggplot2::ggplot(dt_, ggplot2::aes(x = pos_x, y = pos_y)) +
      ggplot2::geom_tile(ggplot2::aes(fill = get(metric)), height = tile_height, width = tile_width) +
      ggplot2::labs(x = x_axis_lab,
                    y = y_axis_lab,
                    title = plt_title,
                    fill = legend_title_fill)  +
      ggplot2::scale_x_continuous(breaks = mrk_x,
                                  labels = lbl_x,
                                  expand = c(0, 0)) +
      ggplot2::scale_y_continuous(breaks = mrk_y,
                                  labels = lbl_y,
                                  expand = c(0, 0)) +
      ggplot2::scale_fill_gradientn(colors = hm_color_palette,
                                    limit = limit_fill,
                                    labels = function(x) sprintf("%.2f", x),
                                    na.value = "lightgrey") + 
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
                     axis.text.y = ggplot2::element_text(size = 8),
                     plot.title = ggplot2::element_text(size = 10),
                     panel.grid = ggplot2::element_blank(),
                     aspect.ratio = 1)
    
    if (show_values) {
      plt <-  plt + 
        ggplot2::geom_text(
          ggplot2::aes(label = ifelse(is.na(get(metric)), "", sprintf("%.2f", get(metric)))),
          size = 2,
          color = "black")
    }
    
    # add isoline
    if (NROW(available_iso_lvl)) { # isobolograms as lines
      if (all(available_iso_lvl %in% c("0.25", "0.5", "0.75"))) {
        # friendly for user with color vision deficiency
        plt <- plt +
          ggplot2::geom_path(data = dt_isobolograms, linewidth = 1,
                             ggplot2::aes(x = if (swap_axes) pos_y else pos_x,
                                          y = if (swap_axes) pos_x else pos_y,
                                          color = iso_level,
                                          linetype = iso_level)) +
          ggplot2::scale_color_manual(values = iso_colors[available_iso_lvl],
                                      breaks = available_iso_lvl,
                                      labels = legend_lbl_iso) +
          ggplot2::scale_linetype_manual(values = c("solid", "twodash", "dashed"),
                                         breaks = available_iso_lvl,
                                         labels = legend_lbl_iso) +
          ggplot2::theme(legend.key.width = ggplot2::unit(3, "line")) +
          ggplot2::labs(color = legend_title_iso,
                        linetype = legend_title_iso)
      } else {
        plt <- plt +
          ggplot2::geom_path(data = dt_isobolograms, linewidth = 1,
                             ggplot2::aes(x = if (swap_axes) pos_y else pos_x,
                                          y = if (swap_axes) pos_x else pos_y,
                                          color = iso_level)) +
          ggplot2::scale_color_manual(values = iso_colors[available_iso_lvl],
                                      breaks = available_iso_lvl,
                                      labels = legend_lbl_iso) +
          ggplot2::labs(color = legend_title_iso)
      }
    }
  }
  return(plt)
}

#' Plot panel of heatmaps of fitted values for combination metrics data
#'
#' @inheritParams heatmap_combo_metrics
#' @param as_list logical flag whether return list of plot or panel
#' @param one_row_panel logical flag whether return panel 2x2 (containing heatmaps for combination metrics and CI plot)
#'    or 3x1 (containing only heatmaps for combination metrics); it is working only for \code{as_list = TRUE}
#'
#' @return \code{ggplot} object containing panel with heatmaps with value for excess assays for 
#'    selected drugs and cell line with selected isoline and comparison of iso levels
#'    or list of \code{ggplot} object containing these plots.
#'    
#' @keywords combo_plots
#' @examples
#' cl_name <- "cellline_BC"
#' drug1_name <- "drug_001"
#' drug2_name <- "drug_026"
#' 
#' mae <- gDRutils::get_synthetic_data("combo_matrix")
#' se <- mae[[gDRutils::get_supported_experiments("combo")]]
#' dt_excess <- gDRutils::convert_se_assay_to_dt(se, "excess")
#' dt_isobolograms <- gDRutils::convert_se_assay_to_dt(se, "isobolograms")
#' 
#' heatmap_combo_metrics_panel(dt_excess,
#'                             dt_isobolograms,
#'                             drug1_name, drug2_name,
#'                             cl_name,
#'                             normalization_type = "GR")
#'                       
#' cl_name <- "cellline_JE"
#' drug1_name <- "drug_011"
#' drug2_name <- "drug_026"
#' 
#' heatmap_combo_metrics_panel(dt_excess,
#'                             dt_isobolograms,
#'                             drug1_name, drug2_name,
#'                             cl_name,
#'                             normalization_type = "RV",
#'                             iso_levels = "0.5",
#'                             as_list = TRUE)
#'                             
#' heatmap_combo_metrics_panel(dt_excess,
#'                             dt_isobolograms,
#'                             drug1_name, drug2_name,
#'                             cl_name,
#'                             normalization_type = "RV",
#'                             iso_levels = "0.5",
#'                             one_row_panel = TRUE)
#' 
#' heatmap_combo_metrics_panel(dt_excess,
#'                             dt_isobolograms,
#'                             drug1_name, drug2_name,
#'                             cl_name,
#'                             normalization_type = "RV",
#'                             iso_levels = NULL,
#'                             as_list = FALSE,
#'                             swap_axes = FALSE)
#' 
#' heatmap_combo_metrics_panel(dt_excess,
#'                             dt_isobolograms,
#'                             drug1_name, drug2_name,
#'                             cl_name,
#'                             normalization_type = "RV",
#'                             iso_levels = NULL,
#'                             as_list = FALSE,
#'                             swap_axes = TRUE)
#'
#' @export
heatmap_combo_metrics_panel <- function(
    dt_excess,
    dt_isobolograms = NULL,
    drug1_name,
    drug2_name,
    cl_name,
    normalization_type = "GR",
    iso_levels = c("0.25", "0.5", "0.75"),
    colors_vec_smooth = NULL,
    colors_vec_excess = NULL,
    no_breaks = 50,
    as_list = FALSE,
    one_row_panel = FALSE,
    swap_axes = FALSE,
    show_values = FALSE) {
  
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  clid <- gDRutils::get_env_identifiers("cellline")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  drug_name_2 <- gDRutils::get_env_identifiers("drug_name2")
  conc <- gDRutils::get_env_identifiers("concentration")
  conc_2 <- gDRutils::get_env_identifiers("concentration2")
  duration <- gDRutils::get_env_identifiers("duration")
  mx_names <- names(gDRutils::get_combo_excess_field_names())
  
  checkmate::assert_data_table(dt_excess)
  checkmate::assert_string(drug1_name)
  checkmate::assert_choice(drug1_name, choices = dt_excess[[drug_name]])
  checkmate::assert_string(drug2_name)
  checkmate::assert_choice(drug2_name, choices = dt_excess[[drug_name_2]])
  checkmate::assert_string(cl_name)
  checkmate::assert_choice(cl_name, choices = dt_excess[[cellline_name]])
  checkmate::assert_choice(normalization_type, choices = c("GR", "RV"))
  checkmate::assert_data_table(dt_isobolograms, null.ok = TRUE)
  if (!is.null(dt_isobolograms)) {
    checkmate::assert_character(iso_levels, null.ok = TRUE)
    checkmate::assert_numeric(as.numeric(iso_levels))
    checkmate::assert_names(names(dt_isobolograms), must.include = "iso_level")
  }
  checkmate::assert_int(no_breaks, lower = 2)
  checkmate::assert_flag(as_list)
  checkmate::assert_flag(one_row_panel)
  checkmate::assert_flag(swap_axes)
  checkmate::assert_flag(show_values)
  hline_color <- 
    gDRutils::get_settings_from_json("HLINE_COLOR",
                                     system.file(package = "gDRplots", "settings.json"))
  
  # data filtering and processing
  filter_expr <- substitute(normalization_type == norm_type, list(norm_type = normalization_type))
  dt_excess <- dt_excess[eval(filter_expr)]
  if (!is.null(dt_isobolograms) && !is.null(iso_levels)) {
    dt_isobolograms <- dt_isobolograms[eval(filter_expr)]
    dt_isobolograms <- dt_isobolograms[get(cellline_name) == cl_name & get(drug_name) == drug1_name &
                                         get(drug_name_2) == drug2_name]
    dt_isobolograms <- dt_isobolograms[iso_level %in% iso_levels, ]
  }
  
  dt_excess <- dt_excess[get(cellline_name) == cl_name & get(drug_name) == drug1_name & get(drug_name_2) == drug2_name]
  
  # check if isolines are available and adjust plotting logic accordingly
  available_iso_lvl <- if (!is.null(dt_isobolograms) && !is.null(iso_levels)) {
    unique(dt_isobolograms[["iso_level"]])
  } else {
    NULL
  }
  
  iso_colors <- if (!is.null(available_iso_lvl)) {
    .get_iso_colors(available_iso_lvl)
  } else {
    NULL
  }
  
  # title
  main_title <- sprintf("%s (%s)",
                        cl_name,
                        unique(dt_excess[get(cellline_name) == cl_name][[clid]]))
  # legend
  legend_title_iso <- "Iso Levels"
  label_prefix <- if (normalization_type == "GR") {
    "GR"
  } else {
    "IC"
  }
  legend_lbl_iso <- sprintf("%s%s", label_prefix, 100 - 100 * as.numeric(available_iso_lvl))
  
  # prep hm color palette
  hm_color_palette_smooth <-
    if (is.null(colors_vec_smooth)  || !all(vapply(colors_vec_smooth, is_valid_color, logical(1)))) {
      .get_smooth_palette(no_breaks)
    } else {
      grDevices::colorRampPalette(colors_vec_smooth)(no_breaks + 1)
    }
  
  hm_color_palette_excess <- 
    if (is.null(colors_vec_excess) || !all(vapply(colors_vec_excess, is_valid_color, logical(1)))) {
      .get_excess_palette(no_breaks)
    } else {
      grDevices::colorRampPalette(colors_vec_excess)(no_breaks + 1)
    }
  
  # plots
  mx_plts <- lapply(mx_names, function(mx_name) {
    # plot title
    plt_title <- sprintf("%s for %s, T=%sh",
                         gDRutils::prettify_flat_metrics(x = mx_name, human_readable = TRUE),
                         normalization_type,
                         unique(dt_excess[get(cellline_name) == cl_name][[duration]]))
    if (as_list) plt_title <- paste(main_title, plt_title, sep = " : ")
    
    dt_ <- dt_excess[, c(conc, conc_2, mx_name), with = FALSE]
    
    x_axis_drug <- if (swap_axes) {
      drug1_name
    } else {
      drug2_name
    }
    
    y_axis_drug <- if (swap_axes) {
      drug2_name
    } else {
      drug1_name
    }
    
    x_axis_lab <- sprintf("%s [\U00B5M]", x_axis_drug)
    y_axis_lab <- sprintf("%s [\U00B5M]", y_axis_drug)
    
    if (!NROW(dt_) > 1 || # co-dilution input data is like: (conc = 0, conc_2 = 0, mx_name = 1)
        all(is.na(dt_[get(conc) != 0 & get(conc_2) != 0][[mx_name]]))) { # lack of smooth & excess data
      plt <- 
        ggplot2::ggplot() +
        ggplot2::labs(x = x_axis_lab,
                      y = y_axis_lab,
                      title = plt_title) +
        ggplot2::theme_bw() +
        ggplot2::theme(aspect.ratio = 1)
    } else {
      dt_[[mx_name]] <- pmin(1.1, dt_[[mx_name]])
      
      conc_y <- if (swap_axes) dt_[[conc_2]] else dt_[[conc]]
      conc_x <- if (swap_axes) dt_[[conc]] else dt_[[conc_2]]
      
      dt_$pos_y <- transform_log_conc(conc_y)
      dt_$pos_x <- transform_log_conc(conc_x)
      
      lbl_y <- sprintf("%.2g", gDRutils::round_concentration(sort(unique(conc_y))))
      mrk_y <- sort(unique(dt_$pos_y))
      
      lbl_x <- sprintf("%.2g", gDRutils::round_concentration(sort(unique(conc_x))))
      mrk_x <- sort(unique(dt_$pos_x))
      
      tile_height <- .get_tile_size(mrk_y)
      tile_width <- .get_tile_size(mrk_x)
      
      # prep hm color palette
      hm_color_palette <- if (mx_name == "smooth") {
        hm_color_palette_smooth
      } else {
        hm_color_palette_excess
      }
      
      # legend title
      legend_title_fill <- sprintf("%s %s",
                                   gDRutils::prettify_flat_metrics(x = mx_name, human_readable = TRUE),
                                   normalization_type)
      
      # prep limits
      limits <- prep_hm_limits(dt_[[mx_name]],   
                               metric = mx_name,
                               normalization_type = normalization_type,
                               symmetric = mx_name != "smooth")
      
      # base plot
      plt <-
        ggplot2::ggplot(dt_, ggplot2::aes(x = pos_x, y = pos_y)) +
        ggplot2::geom_tile(ggplot2::aes(fill = get(mx_name)), height = tile_height, width = tile_width) +
        ggplot2::labs(x = x_axis_lab,
                      y = y_axis_lab,
                      title = plt_title,
                      fill = legend_title_fill)  +
        ggplot2::scale_x_continuous(breaks = mrk_x,
                                    labels = lbl_x,
                                    expand = c(0, 0)) +
        ggplot2::scale_y_continuous(breaks = mrk_y,
                                    labels = lbl_y,
                                    expand = c(0, 0)) +
        ggplot2::scale_fill_gradientn(colors = hm_color_palette,
                                      limit = limits,
                                      labels = function(x) sprintf("%.2f", x),
                                      na.value = "lightgrey") + 
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
                       axis.text.y = ggplot2::element_text(size = 8),
                       plot.title = ggplot2::element_text(size = 10),
                       panel.grid = ggplot2::element_blank(),
                       aspect.ratio = 1)
      
      if (show_values) {
        plt <- plt + 
          ggplot2::geom_text(
            ggplot2::aes(label = ifelse(is.na(get(mx_name)), "", sprintf("%.2f", get(mx_name)))),
            size = 2,
            color = "black")
      }
      
      # add isoline
      if (NROW(available_iso_lvl)) { # isobolograms as lines
        if (all(available_iso_lvl %in% c("0.25", "0.5", "0.75"))) {
          # friendly for user with color vision deficiency
          plt <- plt +
            ggplot2::geom_path(data = dt_isobolograms, linewidth = 1,
                               ggplot2::aes(x = if (swap_axes) pos_y else pos_x,
                                            y = if (swap_axes) pos_x else pos_y,
                                            color = iso_level,
                                            linetype = iso_level)) +
            ggplot2::scale_color_manual(values = iso_colors[available_iso_lvl],
                                        breaks = available_iso_lvl,
                                        labels = legend_lbl_iso) +
            ggplot2::scale_linetype_manual(values = c("solid", "twodash", "dashed"),
                                           breaks = available_iso_lvl,
                                           labels = legend_lbl_iso) +
            ggplot2::theme(legend.key.width = ggplot2::unit(3, "line")) +
            ggplot2::labs(color = legend_title_iso,
                          linetype = legend_title_iso)
        } else {
          plt <- plt +
            ggplot2::geom_path(data = dt_isobolograms, linewidth = 1,
                               ggplot2::aes(x = if (swap_axes) pos_y else pos_x,
                                            y = if (swap_axes) pos_x else pos_y,
                                            color = iso_level)) +
            ggplot2::scale_color_manual(values = iso_colors[available_iso_lvl],
                                        breaks = available_iso_lvl,
                                        labels = legend_lbl_iso) +
            ggplot2::labs(color = legend_title_iso)
        }
      }
    }
    plt
  })
  names(mx_plts) <- mx_names
  
  # isobolograms across range of concentration ratios
  if (NROW(available_iso_lvl)) { # isobolograms as lines
    plt_title <- sprintf("%s for %s, T=%sh",
                         gDRutils::prettify_flat_metrics(x = "smooth", human_readable = TRUE),
                         normalization_type,
                         unique(dt_excess[get(cellline_name) == cl_name][[duration]]))
    
    if (as_list) plt_title <- paste(main_title, plt_title, sep = " : ")
    # base plot
    plt_iso_compare <-
      ggplot2::ggplot(mapping = ggplot2::aes(x = log10_ratio_conc, y = log2_CI)) +
      ggplot2::geom_line(
        data = data.table::data.table(log10_ratio_conc = c(-2, 2), log2_CI = c(0, 0))) +
      ggplot2::geom_hline(yintercept = 0, color = hline_color)
    
    if (all(available_iso_lvl %in% c("0.25", "0.5", "0.75"))) {
      # friendly for user with color vision deficiency
      plt_iso_compare <- plt_iso_compare +
        ggplot2::geom_path(data = dt_isobolograms, linewidth = 0.5,
                           ggplot2::aes(x = log10_ratio_conc, y = log2_CI, color = iso_level, linetype = iso_level)) +
        ggplot2::scale_color_manual(values = iso_colors[available_iso_lvl],
                                    breaks = available_iso_lvl,
                                    labels = legend_lbl_iso) +
        ggplot2::scale_linetype_manual(values = c("solid", "twodash", "dashed"),
                                       breaks = available_iso_lvl,
                                       labels = legend_lbl_iso) +
        ggplot2::theme(legend.key.width = ggplot2::unit(3, "line")) +
        ggplot2::labs(color = legend_title_iso,
                      linetype = legend_title_iso)
    } else {
      plt_iso_compare <- plt_iso_compare +
        ggplot2::geom_path(data = dt_isobolograms, linewidth = 0.5,
                           ggplot2::aes(x = log10_ratio_conc, y = log2_CI, color = iso_level)) +
        ggplot2::scale_color_manual(values = iso_colors[available_iso_lvl],
                                    breaks = available_iso_lvl,
                                    labels = legend_lbl_iso) +
        ggplot2::labs(color = legend_title_iso)
    }
    
    # add x and y scales
    plt_iso_compare <- plt_iso_compare +
      ggplot2::scale_y_continuous(breaks = -5:4, 
                                  labels = c(paste0("1/", 2 ^ (5:1)), 2 ^ (0:4))) +
      ggplot2::scale_x_continuous(breaks = -3:3, 
                                  labels = c(paste0("1/", 10 ^ (3:1)), 10 ^ (0:3))) +
      ggplot2::coord_cartesian(ylim = c(-5, 4)) +
      ggplot2::labs(y = "CI",
                    x = paste(drug2_name, "/", drug1_name, "ratio"),
                    title = plt_title) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
                     axis.text.y = ggplot2::element_text(size = 8),
                     plot.title = ggplot2::element_text(size = 10),
                     panel.grid.minor = ggplot2::element_blank(),
                     aspect.ratio = 1)
    
    # final plots
    ls_plts <- append(mx_plts, list(iso_compare = plt_iso_compare))
  } else {
    ls_plts <- mx_plts
  }
  
  final_plot <- if (as_list) {
    ls_plts
  } else if (one_row_panel) {
    # build panel 3x1
    ggpubr::ggarrange(
      plotlist = list(
        ls_plts[["smooth"]] + ggplot2::guides(linetype = "none", color = "none"),
        ls_plts[["hsa_excess"]] + ggplot2::labs(fill = "Excess"),
        ls_plts[["bliss_excess"]] + ggplot2::labs(fill = "Excess")
      ),
      ncol = 3, common.legend = FALSE, legend = "left")
  } else {
    # build panel 2x2
    ggpubr::annotate_figure(
      ggpubr::ggarrange(
        ggpubr::ggarrange(
          plotlist = list(
            ls_plts[["smooth"]] + ggplot2::guides(linetype = "none", color = "none"),
            ls_plts[["iso_compare"]] + ggplot2::guides(linetype = "none", color = "none")
          ),
          ncol = 2, common.legend = TRUE, legend = "left"),
        ggpubr::ggarrange(
          plotlist = list(
            ls_plts[["hsa_excess"]] + ggplot2::labs(fill = "Excess"),
            ls_plts[["bliss_excess"]] + ggplot2::labs(fill = "Excess")
          ),
          ncol = 2, common.legend = TRUE, legend = "left"),
        common.legend = TRUE, nrow = 2),
      top = main_title) +
      ggpubr::bgcolor("white") + ggpubr::border("white")
  }
  # final
  final_plot
}

#' Plot line plot of combination index
#'
#' @inheritParams heatmap_combo_metrics_panel
#' @param colors_vec_iso character vector of colors (valid name or hex) used for the isolines; 
#'     the default is the dark red-orange palette
#'
#' @return \code{ggplot} object containing combination index plot at different ratios of the two drugs
#'    
#' @keywords combo_plots
#' @examples
#' cl_name <- "cellline_BC"
#' drug1_name <- "drug_001"
#' drug2_name <- "drug_026"
#' 
#' mae <- gDRutils::get_synthetic_data("combo_matrix")
#' se <- mae[[gDRutils::get_supported_experiments("combo")]]
#' dt_excess <- gDRutils::convert_se_assay_to_dt(se, "excess")
#' dt_isobolograms <- gDRutils::convert_se_assay_to_dt(se, "isobolograms")
#' 
#' plot_combination_index(dt_excess,
#'                        dt_isobolograms,
#'                        drug1_name, drug2_name,
#'                        cl_name,
#'                        normalization_type = "GR")
#'                        
#' plot_combination_index(dt_excess,
#'                        dt_isobolograms,
#'                        drug1_name, drug2_name,
#'                        cl_name,
#'                        normalization_type = "RV",                       
#'                        colors_vec_iso = c("darkblue", "darkcyan"))                       
#' 
#' cl_name <- "cellline_JE"
#' drug1_name <- "drug_011"
#' drug2_name <- "drug_026"
#' 
#' plot_combination_index(dt_excess,
#'                        dt_isobolograms,
#'                        drug1_name, drug2_name,
#'                        cl_name,
#'                        normalization_type = "RV",
#'                        iso_levels = "0.5")
#' 
#' @export
plot_combination_index <- function(
    dt_excess,
    dt_isobolograms,
    drug1_name,
    drug2_name,
    cl_name,
    normalization_type = "GR",
    iso_levels =  c("0.25", "0.5", "0.75"),
    colors_vec_iso = NULL) {
  
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  drug_name_2 <- gDRutils::get_env_identifiers("drug_name2")
  duration <- gDRutils::get_env_identifiers("duration")
  
  checkmate::assert_data_table(dt_excess, null.ok = TRUE)
  checkmate::assert_data_table(dt_isobolograms)
  checkmate::assert_string(drug1_name)
  checkmate::assert_choice(drug1_name, choices = dt_isobolograms[[drug_name]])
  checkmate::assert_string(drug2_name)
  checkmate::assert_choice(drug2_name, choices = dt_isobolograms[[drug_name_2]])
  checkmate::assert_string(cl_name)
  checkmate::assert_choice(cl_name, choices = dt_isobolograms[[cellline_name]])
  checkmate::assert_choice(normalization_type, choices = c("GR", "RV"))
  checkmate::assert_character(iso_levels)
  checkmate::assert_numeric(as.numeric(iso_levels))
  checkmate::assert_names(names(dt_isobolograms), must.include = "iso_level")
  checkmate::assert_character(colors_vec_iso, null.ok = TRUE)
  hline_color <- 
    gDRutils::get_settings_from_json("HLINE_COLOR",
                                     system.file(package = "gDRplots", "settings.json"))
  
  # data filtering and processing
  filter_expr <- substitute(normalization_type == norm_type, list(norm_type = normalization_type))
  
  dt_isobolograms <- dt_isobolograms[eval(filter_expr)]
  dt_isobolograms <- dt_isobolograms[get(cellline_name) == cl_name & get(drug_name) == drug1_name &
                                       get(drug_name_2) == drug2_name]
  dt_isobolograms <- dt_isobolograms[iso_level %in% iso_levels, ]
  available_iso_lvl <- unique(dt_isobolograms[["iso_level"]])
  
  if (!is.null(dt_excess)) {
    dt_excess <- dt_excess[eval(filter_expr)]
    dt_excess <- dt_excess[get(cellline_name) == cl_name & get(drug_name) == drug1_name & 
                             get(drug_name_2) == drug2_name]
  }
  
  # title
  plt_title <- if (!is.null(dt_excess)) {
    sprintf("%s for %s, T=%sh",
            gDRutils::prettify_flat_metrics(x = "smooth", human_readable = TRUE),
            normalization_type,
            unique(dt_excess[get(cellline_name) == cl_name][[duration]]))
  } else {
    sprintf("%s for %s",
            gDRutils::prettify_flat_metrics(x = "smooth", human_readable = TRUE),
            normalization_type)
  }
  
  # base plot
  plt <-
    ggplot2::ggplot(mapping = ggplot2::aes(x = log10_ratio_conc, y = log2_CI)) +
    ggplot2::geom_line(
      data = data.table::data.table(log10_ratio_conc = c(-2, 2), 
                                    log2_CI = c(0, 0))) +
    ggplot2::geom_hline(yintercept = 0, color = hline_color)
  
  # check if isolines are available and adjust plotting logic accordingly
  if (NROW(available_iso_lvl) > 0) {
    legend_title_iso <- "Iso Levels"
    label_prefix <- if (normalization_type == "GR") {
      "GR"
    } else {
      "IC"
    }
    legend_lbl_iso <- sprintf("%s%s", label_prefix, 100 - 100 * as.numeric(available_iso_lvl))
    
    iso_colors <- 
      if (is.null(colors_vec_iso) || !all(vapply(colors_vec_iso, is_valid_color, logical(1)))) {
        .get_iso_colors(available_iso_lvl)
      } else {
        ls_ <- grDevices::colorRampPalette(colors_vec_iso)(NROW(available_iso_lvl))
        names(ls_) <- available_iso_lvl
        ls_ 
      }
    
    if (all(available_iso_lvl %in% c("0.25", "0.5", "0.75"))) {
      # friendly for user with color vision deficiency
      plt <- plt +
        ggplot2::geom_path(data = dt_isobolograms, linewidth = 0.5,
                           ggplot2::aes(x = log10_ratio_conc, 
                                        y = log2_CI, 
                                        color = iso_level, 
                                        linetype = iso_level)) +
        ggplot2::scale_color_manual(values = iso_colors[available_iso_lvl],
                                    breaks = available_iso_lvl,
                                    labels = legend_lbl_iso,
                                    name = legend_title_iso) +
        ggplot2::scale_linetype_manual(values = c("solid", "twodash", "dashed"),
                                       breaks = available_iso_lvl,
                                       labels = legend_lbl_iso,
                                       name = legend_title_iso) +
        ggplot2::theme(legend.key.width = ggplot2::unit(3, "line"))
    } else {
      plt <- plt +
        ggplot2::geom_path(data = dt_isobolograms, linewidth = 0.5,
                           ggplot2::aes(x = log10_ratio_conc, 
                                        y = log2_CI, 
                                        color = iso_level)) +
        ggplot2::scale_color_manual(values = iso_colors[available_iso_lvl],
                                    breaks = available_iso_lvl,
                                    labels = legend_lbl_iso,
                                    name = legend_title_iso)
    }
  }
  
  # add x and y scales
  plt <- plt +
    ggplot2::scale_y_continuous(breaks = -5:4, 
                                labels = c(paste0("1/", 2 ^ (5:1)), 2 ^ (0:4))) +
    ggplot2::scale_x_continuous(breaks = -3:3, 
                                labels = c(paste0("1/", 10 ^ (3:1)), 10 ^ (0:3))) +
    ggplot2::coord_cartesian(ylim = c(-5, 4)) +
    ggplot2::labs(y = "CI",
                  x = paste(drug2_name, "/", drug1_name, "ratio"),
                  title = plt_title) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
                   axis.text.y = ggplot2::element_text(size = 8),
                   plot.title = ggplot2::element_text(size = 10),
                   panel.grid.minor = ggplot2::element_blank(),
                   aspect.ratio = 1)
  
  # final
  return(plt)
}


#' Plot heatmaps of averaged values for combination data
#'
#' @param dt_excess data.table representing data from the \code{excess} assay,
#'    outputted by \code{gDRutils::convert_se_assay_to_dt(se, "excess")}
#'    and combo \code{SummarizedExperiment}
#' @param dt_isobolograms data.table representing data from the \code{isobolograms} assay,
#'    outputted by \code{gDRutils::convert_se_assay_to_dt(se, "isobolograms")}
#'    and combo \code{SummarizedExperiment}
#' @param drug1_name string with drug name to be plotted (identifiers \code{DrugName})
#' @param drug2_name string with co-drug name to be plotted (identifiers \code{DrugName_2})
#' @param cl_name string with cell line to be plotted (identifiers \code{CellLineName})
#' @param normalization_type string with normalization_types to be selected
#'                           one of: "GR" ("GRvalue") or "RV" ("RelativeViability")
#' @param metric string name of the combo metric;
#'   one of: "smooth" ("Smooth GR" or "Smooth RV" - respectively depending on \code{normalization_type})
#'   "hsa_excess" ("Bliss Excess GR" or "Bliss Excess RV") or "bliss_excess" 
#'   ("Bliss Excess GR" or "Bliss Excess RV")
#' @param iso_levels character vector with  isobologram levels to be selected;
#'     when \code{NULL} - no isolines will be displayed
#' @param colors_vec character vector of colors (valid name or hex) used in heatmap; 
#'     the default depends on \code{metric}: for "smooth" - the dark purple-light grey palette
#'     and for "hsa_excess" and "bliss_excess" - the blue-light grey-red color scale
#' @param colors_vec_iso character vector of colors (valid name or hex) used for the isolines; 
#'     the default is the dark red-orange palette
#' @param no_breaks numeric number of breaks on scale
#' @param swap_axes logical flag whether to swap the axes with drugs of the heatmap
#'
#' @return \code{ggplot} object containing heatmap for fitted values and reference data 
#'    for isobolograms for selected drug and co-drug and selected cell line
#'    
#' @keywords combo_plots
#' @examples
#' cl_name <- "cellline_BC"
#' drug1_name <- "drug_001"
#' drug2_name <- "drug_026"
#' 
#' mae <- gDRutils::get_synthetic_data("combo_matrix")
#' se <- mae[[gDRutils::get_supported_experiments("combo")]]
#' dt_excess <- gDRutils::convert_se_assay_to_dt(se, "excess")
#' dt_isobolograms <- gDRutils::convert_se_assay_to_dt(se, "isobolograms")
#' 
#' heatmap_combo_with_isoref(dt_excess,
#'                           dt_isobolograms,
#'                           drug1_name, drug2_name,
#'                           cl_name)
#'                           
#' heatmap_combo_with_isoref(dt_excess,
#'                           dt_isobolograms,
#'                           drug1_name, drug2_name,
#'                           cl_name,
#'                           metric = "hsa_excess",
#'                           iso_levels = c("-0.2", "0.2"))
#'                           
#' heatmap_combo_with_isoref(dt_excess,
#'                           dt_isobolograms,
#'                           drug1_name, drug2_name,
#'                           cl_name,
#'                           iso_levels = NULL,
#'                           colors_vec = c("darkcyan", "snow", "darkorange"))                      
#'  
#' heatmap_combo_with_isoref(dt_excess,
#'                           dt_isobolograms,
#'                           drug1_name, drug2_name,
#'                           cl_name,
#'                           normalization_type = "RV",
#'                           iso_levels = c("0.25", "0.75"),
#'                           swap_axes = FALSE)
#'                                                    
#' heatmap_combo_with_isoref(dt_excess,
#'                           dt_isobolograms,
#'                           drug1_name, drug2_name,
#'                           cl_name,
#'                           normalization_type = "RV",
#'                           iso_levels = c("0.25", "0.75"),
#'                           swap_axes = TRUE)
#'
#' heatmap_combo_with_isoref(dt_excess,
#'                           dt_isobolograms,
#'                           drug1_name, drug2_name,
#'                           cl_name,
#'                           metric = "hsa_excess",
#'                           iso_levels = c("0.25", "0.75"),
#'                           colors_vec_iso = c("0.25" = "darkcyan",
#'                                              "0.75" = "darkblue"),
#'                           swap_axes = FALSE)
#'                           
#' @export
heatmap_combo_with_isoref <- function(
    dt_excess,
    dt_isobolograms,
    drug1_name,
    drug2_name,
    cl_name,
    normalization_type = "GR",
    metric = "smooth",
    iso_levels = "0.5",
    colors_vec = NULL,
    colors_vec_iso = NULL,
    no_breaks = 50,
    swap_axes = FALSE) {
  
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  clid <- gDRutils::get_env_identifiers("cellline")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  drug_name_2 <- gDRutils::get_env_identifiers("drug_name2")
  conc <- gDRutils::get_env_identifiers("concentration")
  conc_2 <- gDRutils::get_env_identifiers("concentration2")
  
  checkmate::assert_data_table(dt_excess)
  checkmate::assert_data_table(dt_isobolograms)
  checkmate::assert_string(drug1_name)
  checkmate::assert_choice(drug1_name, choices = dt_excess[[drug_name]])
  checkmate::assert_string(drug2_name)
  checkmate::assert_choice(drug2_name, choices = dt_excess[[drug_name_2]])
  checkmate::assert_string(cl_name)
  checkmate::assert_choice(cl_name, choices = dt_excess[[cellline_name]])
  checkmate::assert_choice(normalization_type, choices = c("GR", "RV"))
  checkmate::assert_choice(metric, choices = names(gDRutils::get_combo_excess_field_names()))
  checkmate::assert_character(colors_vec_iso, null.ok = TRUE)
  if (!is.null(iso_levels)) checkmate::assert_numeric(as.numeric(iso_levels))
  checkmate::assert_character(colors_vec, null.ok = TRUE)
  checkmate::assert_character(iso_levels, null.ok = TRUE)
  checkmate::assert_int(no_breaks, lower = 2)
  checkmate::assert_flag(swap_axes)
  
  # filter data for normalization type
  filter_expr <- substitute(normalization_type == norm_type, list(norm_type = normalization_type))
  dt_excess <- dt_excess[eval(filter_expr)]
  dt_isobolograms <- dt_isobolograms[eval(filter_expr)]
  
  # filter data for combination cell line (drug x drug2)
  dt_excess <-
    dt_excess[get(cellline_name) == cl_name & get(drug_name) == drug1_name & get(drug_name_2) == drug2_name]
  dt_isobolograms <-
    dt_isobolograms[get(cellline_name) == cl_name & get(drug_name) == drug1_name & get(drug_name_2) == drug2_name]
  
  # prep hm color palette
  hm_color_palette <- 
    if (is.null(colors_vec) || !all(vapply(colors_vec, is_valid_color, logical(1)))) {
      if (metric == "smooth") {
        .get_smooth_palette(no_breaks)
      } else {
        .get_excess_palette(no_breaks)
      }
    } else {
      grDevices::colorRampPalette(colors_vec)(no_breaks + 1)
    }
  
  # panel title
  cl_clid <- unique(dt_excess[get(cellline_name) == cl_name, ][[clid]])
  plt_title <- sprintf("%s (%s)", cl_name, cl_clid)
  
  # prep plot data
  dt_ <- dt_excess[, c(conc, conc_2, metric), with = FALSE]
  
  x_axis_drug <- if (swap_axes) {
    drug1_name
  } else {
    drug2_name
  }
  
  y_axis_drug <- if (swap_axes) {
    drug2_name
  } else {
    drug1_name
  }
  
  x_axis_lab <- sprintf("%s [\U00B5M]", x_axis_drug)
  y_axis_lab <- sprintf("%s [\U00B5M]", y_axis_drug)
  
  if (!NROW(dt_) > 1) { # co-dilution input data is like: (conc = 0, conc_2 = 0, metric = 1)
    plt <- 
      ggplot2::ggplot() +
      ggplot2::labs(x = x_axis_lab,
                    y = y_axis_lab,
                    title = plt_title) +
      ggplot2::theme_bw() +
      ggplot2::theme(aspect.ratio = 1)
  } else {
    dt_[[metric]] <- pmin(1.1, dt_[[metric]])
    
    conc_y <- if (swap_axes) dt_[[conc_2]] else dt_[[conc]]
    conc_x <- if (swap_axes) dt_[[conc]] else dt_[[conc_2]]
    
    dt_$pos_y <- transform_log_conc(conc_y)
    dt_$pos_x <- transform_log_conc(conc_x)
    
    lbl_y <- sprintf("%.2g", gDRutils::round_concentration(sort(unique(conc_y))))
    mrk_y <- sort(unique(dt_$pos_y))
    
    lbl_x <- sprintf("%.2g", gDRutils::round_concentration(sort(unique(conc_x))))
    mrk_x <- sort(unique(dt_$pos_x))
    
    tile_height <- .get_tile_size(mrk_y)
    tile_width <- .get_tile_size(mrk_x)
    
    range_x <- c(min(mrk_x) - 0.65 * tile_width, max(mrk_x) + 0.65 * tile_width)
    range_y <- c(min(mrk_y) - 0.65 * tile_height, max(mrk_y) + 0.65 * tile_height)
    
    range_xy <- c(min(range_x[1], range_y[1]), max(range_x[2], range_y[2]))
    
    # legend title
    legend_title_fill <- sprintf("%s %s",
                                 gDRutils::prettify_flat_metrics(x = metric, human_readable = TRUE),
                                 normalization_type)
    
    # prep limits
    limits <- prep_hm_limits(num_vec = dt_[[metric]],   
                             metric = metric,
                             normalization_type = normalization_type,
                             symmetric = (metric != "smooth"))
    
    # base plot
    plt <-
      ggplot2::ggplot(dt_, ggplot2::aes(x = pos_x, y = pos_y)) +
      ggplot2::geom_tile(ggplot2::aes(fill = get(metric)), 
                         height = tile_height, width = tile_width, alpha = 0.90) +
      ggplot2::labs(x = x_axis_lab,
                    y = y_axis_lab,
                    title = plt_title,
                    fill = legend_title_fill) +
      ggplot2::scale_x_continuous(breaks = mrk_x,
                                  labels = lbl_x,
                                  expand = c(0, 0)) +
      ggplot2::scale_y_continuous(breaks = mrk_y,
                                  labels = lbl_y,
                                  expand = c(0, 0)) +
      ggplot2::scale_fill_gradientn(colors = hm_color_palette,
                                    limit = limits,
                                    labels = function(x) sprintf("%.2f", x),
                                    na.value = "lightgrey") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
                     axis.text.y = ggplot2::element_text(size = 8),
                     plot.title = ggplot2::element_text(size = 10),
                     panel.grid = ggplot2::element_blank(),
                     legend.key.width = ggplot2::unit(2, "line"),
                     aspect.ratio = 1)
    
    # plot isobologram
    if (!is.null(dt_isobolograms$iso_level) && !is.null(iso_levels)) { # add isolines - if there are such data
      # iso level availability
      dt_isobolograms <- dt_isobolograms[iso_level %in% iso_levels, ]
      available_iso_lvl <- unique(dt_isobolograms[["iso_level"]])
      iso_levels <- iso_levels[iso_levels %in% available_iso_lvl]
      
      if (NROW(iso_levels)) {
        # order iso level
        iso_levels <- iso_levels[order(as.numeric(iso_levels))]
        
        req_cols <- c(cellline_name, drug_name, drug_name_2, gDRutils::get_header("iso_position"))
        dt_iso <- 
          dt_isobolograms[iso_level %in% iso_levels, .SD, .SDcols = req_cols]
        
        # colors for isoline
        iso_colors <- 
          if (is.null(colors_vec_iso) || !all(vapply(colors_vec_iso, is_valid_color, logical(1)))) {
            .get_iso_colors(iso_levels)
          } else if (all(iso_levels %in% names(colors_vec_iso))) {
            colors_vec_iso[iso_levels]
          } else {
            ls_ <- grDevices::colorRampPalette(colors_vec_iso)(NROW(available_iso_lvl))
            names(ls_) <- iso_levels
            ls_ 
          }
        
        # plot
        label_prefix <- if (normalization_type == "GR") {
          "GR"
        } else {
          "IC"
        }
        
        iso_label <- sprintf("%s%s", label_prefix, 100 - 100 * as.numeric(available_iso_lvl))
        names(iso_label) <- available_iso_lvl
        
        tab_measured <- dt_iso[, .SD, .SDcols = -c("pos_x_ref", "pos_y_ref")]
        tab_measured[, iso_source := "measured"]
        tab_expected <- dt_iso[, .SD, .SDcols = -c("pos_x", "pos_y")]
        tab_expected[, iso_source := "expected"]
        data.table::setnames(tab_expected, old = c("pos_x_ref", "pos_y_ref"), new = c("pos_x", "pos_y"))
        
        tab_isoline <- rbind(tab_measured, tab_expected)
        # adjust isoline range to heatmap
        tab_isoline <-
          tab_isoline[data.table::between(pos_x, range_xy[1], range_xy[2]) &
                        data.table::between(pos_y, range_xy[1], range_xy[2]), ]
        
        
        if (NROW(iso_levels) == 1) {
          plt <- plt +
            ggplot2::geom_path(data = tab_isoline,
                               ggplot2::aes(x = if (swap_axes) pos_y else pos_x,
                                            y = if (swap_axes) pos_x else pos_y,
                                            linetype = iso_source),
                               linewidth = 1, color = iso_colors) +
            ggplot2::scale_linetype_manual(values = c("measured" = "solid", "expected" = "twodash")) +
            ggplot2::labs(linetype = as.character(iso_label))
        } else {
          plt <- plt +
            ggplot2::geom_path(data = tab_isoline,
                               ggplot2::aes(x = if (swap_axes) pos_y else pos_x,
                                            y = if (swap_axes) pos_x else pos_y,
                                            linetype = iso_source,
                                            color = iso_level),
                               linewidth = 1) +
            ggplot2::scale_linetype_manual(values = c("measured" = "solid", "expected" = "twodash")) +
            ggplot2::scale_color_manual(values = iso_colors,
                                        breaks = names(iso_colors),
                                        labels = iso_label) +
            ggplot2::labs(linetype = normalization_type,
                          color = "Iso Levels")
        }
      }
    }
    # final plot
    plt <- plt +
      ggplot2::guides(fill = ggplot2::guide_colorbar(order = 1),
                      linetype = ggplot2::guide_legend(order = 2), 
                      color = ggplot2::guide_legend(order = 3))
  }
  
  return(plt)
}

#' Plot panel of heatmaps with fitted and reference data for isobolograms
#' 
#' @inheritParams heatmap_combo_with_isoref
#' @param cl_names character vector with cell line names to be plotted (Cell Line Name);
#'    if \code{NULL} - all available cell lines will be plotted
#'
#' @return \code{ggplot} object containing panel with heatmaps for fitted values and reference data 
#'    for isobolograms for selected drug and co-drug and all selected cell line
#'    
#' @keywords combo_plots
#' @examples
#' cl_names <-
#'   c("cellline_AA", "cellline_EA", "cellline_IB",
#'   "cellline_MC", "cellline_BC", "cellline_FD")
#' 
#' drug1_name <- "drug_001"
#' drug2_name <- "drug_026"
#' 
#' mae <- gDRutils::get_synthetic_data("combo_matrix")
#' se <- mae[[gDRutils::get_supported_experiments("combo")]]
#' dt_excess <- gDRutils::convert_se_assay_to_dt(se, "excess")
#' dt_isobolograms <- gDRutils::convert_se_assay_to_dt(se, "isobolograms")
#' 
#' heatmap_combo_with_isoref_panel(dt_excess,
#'                                 dt_isobolograms,
#'                                 drug1_name, drug2_name,
#'                                 cl_names)
#' 
#' dt_excess_2 <- data.table::copy(dt_excess)
#' invisible(dt_excess_2[CellLineName %in% cl_names[1:2],
#'                       Concentration := Concentration / 10])
#' 
#' heatmap_combo_with_isoref_panel(dt_excess_2,
#'                                 dt_isobolograms,
#'                                 drug1_name, drug2_name,
#'                                 cl_names,
#'                                 iso_levels = c("0.25", "0.75"),
#'                                 colors_vec = c("darkcyan", "snow", "coral"))
#' 
#' @export
heatmap_combo_with_isoref_panel <- function(
    dt_excess,
    dt_isobolograms,
    drug1_name,
    drug2_name,
    cl_names,
    normalization_type = "GR",
    metric = "smooth",
    iso_levels = "0.5",
    colors_vec = NULL,
    no_breaks = 50,
    swap_axes = FALSE) {
  
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  gnumber <- gDRutils::get_env_identifiers("drug")
  drug_name_2 <- gDRutils::get_env_identifiers("drug_name2")
  gnumber_2 <- gDRutils::get_env_identifiers("drug2")
  conc <- gDRutils::get_env_identifiers("concentration")
  conc_2 <- gDRutils::get_env_identifiers("concentration2")
  
  checkmate::assert_data_table(dt_excess)
  checkmate::assert_data_table(dt_isobolograms)
  checkmate::assert_string(drug1_name)
  checkmate::assert_choice(drug1_name, choices = dt_excess[[drug_name]])
  checkmate::assert_string(drug2_name)
  checkmate::assert_choice(drug2_name, choices = dt_excess[[drug_name_2]])
  checkmate::assert_character(cl_names, null.ok = TRUE)
  checkmate::assert_choice(normalization_type, choices = c("GR", "RV"))
  checkmate::assert_choice(metric, choices = names(gDRutils::get_combo_excess_field_names()))
  checkmate::assert_character(iso_levels, null.ok = TRUE)
  if (!is.null(iso_levels)) {
    stopifnot("`iso_levels` must be a valid numeric value" = 
                all(vapply(iso_levels, function(i) grepl("^[-]*0\\.?[0-9]*$", i), logical(1))))
  }
  checkmate::assert_character(colors_vec, null.ok = TRUE)
  checkmate::assert_int(no_breaks, lower = 2)
  checkmate::assert_flag(swap_axes)
  
  available_cls <- unique(dt_excess[[cellline_name]])
  if (is.null(cl_names) || all(!cl_names %in% available_cls)) {
    cl_names  <- available_cls
  } else if (!all(cl_names %in% available_cls)) {
    cl_names <- cl_names[cl_names %in% available_cls]
  }
  
  # filter data for normalization type
  filter_expr <- substitute(normalization_type == norm_type, list(norm_type = normalization_type))
  dt_excess <- dt_excess[eval(filter_expr)]
  
  # filter data for combination cell line (drug x drug2)
  dt_excess <-
    dt_excess[get(cellline_name) %in% cl_names & get(drug_name) == drug1_name & get(drug_name_2) == drug2_name]
  
  # check whether concentrations are  common for all cell line
  ls_vec_conc <- lapply(cl_names, function(cl_nm) {
    unique(dt_excess[get(cellline_name) == cl_nm, ][[conc]])
  })
  ls_vec_conc_2 <- lapply(cl_names, function(cl_nm) {
    unique(dt_excess[get(cellline_name) == cl_nm, ][[conc_2]])
  })
  
  panel <- if ("independent" %in% c(.get_combo_panel_type(ls_vec_conc), 
                                    .get_combo_panel_type(ls_vec_conc_2))) {
    heatmap_combo_with_isoref_panel_independent(
      dt_excess = dt_excess,
      dt_isobolograms = dt_isobolograms,
      drug1_name = drug1_name,
      drug2_name = drug2_name,
      cl_names = cl_names,
      normalization_type = normalization_type,
      metric = metric,
      iso_levels = iso_levels,
      colors_vec = colors_vec,
      no_breaks = no_breaks,
      swap_axes = swap_axes)
  } else {
    heatmap_combo_with_isoref_panel_common(
      dt_excess = dt_excess,
      dt_isobolograms = dt_isobolograms,
      drug1_name = drug1_name,
      drug2_name = drug2_name,
      cl_names = cl_names,
      normalization_type = normalization_type,
      metric = metric,
      iso_levels = iso_levels,
      colors_vec = colors_vec,
      no_breaks = no_breaks,
      swap_axes = swap_axes)
  }
  
  # final
  return(panel)
}


#' #' Plot panel of heatmaps with fitted and reference data for isobolograms
#' 
#' @inheritParams heatmap_combo_with_isoref_panel
#'    
#' @return \code{ggplot} object containing panel with heatmaps for fitted values and reference data 
#'    for isobolograms for selected drug and co-drug by list of cell lines
#'    
#' @keywords combo_plots
#' @examples
#' cl_names <-
#'   c("cellline_AA", "cellline_EA", "cellline_IB", 
#'   "cellline_MC", "cellline_BC", "cellline_FD")
#' 
#' drug1_name <- "drug_001"
#' drug2_name <- "drug_026"
#' 
#' mae <- gDRutils::get_synthetic_data("combo_matrix")
#' se <- mae[[gDRutils::get_supported_experiments("combo")]]
#' dt_excess <- gDRutils::convert_se_assay_to_dt(se, "excess")
#' dt_isobolograms <- gDRutils::convert_se_assay_to_dt(se, "isobolograms")
#' 
#' heatmap_combo_with_isoref_panel_common(dt_excess,
#'                                        dt_isobolograms,
#'                                        drug1_name, drug2_name,
#'                                        cl_names)
#' 
#' heatmap_combo_with_isoref_panel_common(dt_excess,
#'                                        dt_isobolograms,
#'                                        drug1_name, drug2_name,
#'                                        cl_names,
#'                                        iso_levels = c("0.25", "0.5"))
#' 
#' heatmap_combo_with_isoref_panel_common(dt_excess,
#'                                        dt_isobolograms,
#'                                        drug1_name, drug2_name,
#'                                        cl_names,
#'                                        metric = "hsa_excess",
#'                                        iso_levels = c("-0.25", "0.25"))
#' 
#' heatmap_combo_with_isoref_panel_common(dt_excess,
#'                                        dt_isobolograms,
#'                                        drug1_name, drug2_name,
#'                                        cl_names,
#'                                        normalization_type = "RV",
#'                                        iso_levels = NULL,
#'                                        colors_vec = c("darkcyan", "snow", "darkorange"),
#'                                        swap_axes = FALSE)
#' 
#' heatmap_combo_with_isoref_panel_common(dt_excess,
#'                                        dt_isobolograms,
#'                                        drug1_name, drug2_name,
#'                                        cl_names,
#'                                        normalization_type = "RV",
#'                                        iso_levels = NULL,
#'                                        colors_vec = c("darkcyan", "snow", "darkorange"),
#'                                        swap_axes = TRUE)
#' 
#' heatmap_combo_with_isoref_panel_common(dt_excess,
#'                                        dt_isobolograms,
#'                                        drug1_name, drug2_name,
#'                                        cl_names,
#'                                        metric = "hsa_excess",
#'                                        iso_levels = NULL,
#'                                        swap_axes = FALSE)
#' 
#' @export
heatmap_combo_with_isoref_panel_common <- function(
    dt_excess,
    dt_isobolograms,
    drug1_name,
    drug2_name,
    cl_names,
    normalization_type = "GR",
    metric = "smooth",
    iso_levels = "0.5",
    colors_vec = NULL,
    no_breaks = 50,
    swap_axes = FALSE) {
  
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  gnumber <- gDRutils::get_env_identifiers("drug")
  drug_name_2 <- gDRutils::get_env_identifiers("drug_name2")
  gnumber_2 <- gDRutils::get_env_identifiers("drug2")
  conc <- gDRutils::get_env_identifiers("concentration")
  conc_2 <- gDRutils::get_env_identifiers("concentration2")
  
  checkmate::assert_data_table(dt_excess)
  checkmate::assert_data_table(dt_isobolograms)
  checkmate::assert_string(drug1_name)
  checkmate::assert_choice(drug1_name, choices = dt_excess[[drug_name]])
  checkmate::assert_string(drug2_name)
  checkmate::assert_choice(drug2_name, choices = dt_excess[[drug_name_2]])
  checkmate::assert_character(cl_names, null.ok = TRUE)
  checkmate::assert_choice(normalization_type, choices = c("GR", "RV"))
  checkmate::assert_choice(metric, choices = names(gDRutils::get_combo_excess_field_names()))
  checkmate::assert_character(iso_levels, null.ok = TRUE)
  if (!is.null(iso_levels)) {
    stopifnot("`iso_levels` must be a valid numeric value" = 
                all(vapply(iso_levels, function(i) grepl("^[-]*0\\.?[0-9]*$", i), logical(1))))
  }
  checkmate::assert_character(colors_vec, null.ok = TRUE)
  checkmate::assert_int(no_breaks, lower = 2)
  checkmate::assert_flag(swap_axes)
  
  available_cls <- unique(dt_excess[[cellline_name]])
  if (is.null(cl_names) || all(!cl_names %in% available_cls)) {
    cl_names  <- available_cls
  } else if (!all(cl_names %in% available_cls)) {
    cl_names <- cl_names[cl_names %in% available_cls]
  }
  
  # panel title
  panel_title <- sprintf("%s (%s) x %s (%s)",
                         drug1_name,
                         unique(dt_excess[get(drug_name) == drug1_name, ][[gnumber]]),
                         drug2_name,
                         unique(dt_excess[get(drug_name_2) == drug2_name, ][[gnumber_2]]))
  
  # filter data for normalization type
  filter_expr <- substitute(normalization_type == norm_type, list(norm_type = normalization_type))
  dt_excess <- dt_excess[eval(filter_expr)]
  dt_isobolograms <- dt_isobolograms[eval(filter_expr)]
  
  # filter data for combination cell line (drug x drug2)
  selected_combination <-
    unique(dt_excess[get(cellline_name) %in% cl_names & get(drug_name) == drug1_name & get(drug_name_2) == drug2_name, 
                     .SD, .SDcols = c(cellline_name, drug_name, drug_name_2)])
  
  dt_excess <-
    dt_excess[selected_combination, on = c(cellline_name, drug_name, drug_name_2)]
  dt_isobolograms <-
    dt_isobolograms[selected_combination, on = c(cellline_name, drug_name, drug_name_2)]
  
  # check for overlap of concentration
  ls_vec_conc <- lapply(cl_names, function(cl_nm) {
    unique(dt_excess[get(cellline_name) == cl_nm, ][[conc]])
  })
  ls_vec_conc_2 <- lapply(cl_names, function(cl_nm) {
    unique(dt_excess[get(cellline_name) == cl_nm, ][[conc_2]])
  })
  if (.get_combo_panel_type(ls_vec_conc) != "common") {
    stop("Concentration values for drug 1 are not common for all selected cell lines.
          Consider using `heatmap_combo_with_isoref_panel_independent` function.")
  }
  if (.get_combo_panel_type(ls_vec_conc_2) != "common") {
    stop("Concentration values for drug 2 are not common for all selected cell lines.
          Consider using `heatmap_combo_with_isoref_panel_independent` function.")
  }
  
  # prep hm color palette
  hm_color_palette <- 
    if (is.null(colors_vec) || !all(vapply(colors_vec, is_valid_color, logical(1)))) {
      if (metric == "smooth") {
        .get_smooth_palette(no_breaks)
      } else {
        .get_excess_palette(no_breaks)
      }
    } else {
      grDevices::colorRampPalette(colors_vec)(no_breaks + 1)
    }
  
  # prep panel elements
  dt_all <- dt_excess[, c(cellline_name, conc, conc_2, metric), with = FALSE]
  
  conc_y <- if (swap_axes) {
    conc_2
  } else {
    conc
  }
  
  conc_x <- if (swap_axes) {
    conc
  } else {
    conc_2
  }
  
  dt_tile <- dt_all[get(cellline_name) %in% cl_names, ][, 
                                                        `:=`(
                                                          metric = pmin(1.1, get(metric)),
                                                          pos_y = transform_log_conc(get(conc_y)),
                                                          pos_x = transform_log_conc(get(conc_x))
                                                        ), 
                                                        by = cellline_name
  ][, .SD, .SDcols = -metric]
  data.table::setnames(dt_tile, "metric", metric)
  
  # tiles positioning 
  dt_tile$pos_y <- transform_log_conc(dt_tile[[conc_y]])
  dt_tile$pos_x <- transform_log_conc(dt_tile[[conc_x]])
  
  lbl_y <- sprintf("%.2g", gDRutils::round_concentration(sort(unique(dt_tile[[conc_y]]))))
  mrk_y <- sort(unique(dt_tile$pos_y))
  
  lbl_x <- sprintf("%.2g", gDRutils::round_concentration(sort(unique(dt_tile[[conc_x]]))))
  mrk_x <- sort(unique(dt_tile$pos_x))
  
  tile_height <- .get_tile_size(mrk_y)
  tile_width <- .get_tile_size(mrk_x)
  
  # plot range
  range_x <- c(min(mrk_x) - 0.65 * tile_width, max(mrk_x) + 0.65 * tile_width)
  range_y <- c(min(mrk_y) - 0.65 * tile_height, max(mrk_y) + 0.65 * tile_height)
  
  range_xy <- c(min(range_x[1], range_y[1]), max(range_x[2], range_y[2]))
  
  # prep limits
  limits <- prep_hm_limits(dt_tile[[metric]],   
                           metric = metric,
                           normalization_type = normalization_type,
                           symmetric = (metric != "smooth"))
  # legend title
  legend_title_fill <- sprintf("%s %s",
                               gDRutils::prettify_flat_metrics(x = metric, human_readable = TRUE),
                               normalization_type)
  # base plot
  x_axis_drug <- if (swap_axes) {
    drug1_name
  } else {
    drug2_name
  }
  
  y_axis_drug <- if (swap_axes) {
    drug2_name
  } else {
    drug1_name
  }
  
  x_axis_lab <- sprintf("%s [\U00B5M]", x_axis_drug)
  y_axis_lab <- sprintf("%s [\U00B5M]", y_axis_drug)
  
  plt <-
    ggplot2::ggplot(dt_tile,
                    ggplot2::aes(x = pos_x, y = pos_y)) +
    ggplot2::geom_tile(ggplot2::aes(fill = get(metric)), 
                       height = tile_height, width = tile_width, alpha = 0.90) +
    ggplot2::labs(x = x_axis_lab,
                  y = y_axis_lab,
                  title = panel_title,
                  fill = legend_title_fill) +
    ggplot2::scale_x_continuous(breaks = mrk_x,
                                labels = lbl_x,
                                expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = mrk_y,
                                labels = lbl_y,
                                expand = c(0, 0)) +
    ggplot2::scale_fill_gradientn(colors = hm_color_palette,
                                  limit = limits,
                                  labels = function(x) sprintf("%.2f", x),
                                  na.value = "lightgrey") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
                   axis.text.y = ggplot2::element_text(size = 8),
                   plot.title = ggplot2::element_text(size = 10),
                   panel.grid = ggplot2::element_blank(),
                   legend.key.width = ggplot2::unit(2, "line"),
                   aspect.ratio = 1)
  
  # isoline data
  if (!is.null(dt_isobolograms$iso_level) && !is.null(iso_levels)) {
    # iso level availability
    dt_isobolograms <- dt_isobolograms[iso_level %in% iso_levels, ]
    available_iso_lvl <- unique(dt_isobolograms[["iso_level"]])
    iso_levels <- iso_levels[iso_levels %in% available_iso_lvl]
    
    if (NROW(iso_levels)) {
      # order iso level
      iso_levels <- iso_levels[order(as.numeric(iso_levels))]
      
      req_cols <- c(cellline_name, drug_name, drug_name_2, gDRutils::get_header("iso_position"))
      dt_iso <- 
        dt_isobolograms[iso_level %in% iso_levels, .SD, .SDcols = req_cols]
      
      # colors for isoline
      iso_colors <- .get_iso_colors(iso_levels)
      
      # plot
      label_prefix <- if (normalization_type == "GR") {
        "GR"
      } else {
        "IC"
      }
      iso_label <- sprintf("%s%s", label_prefix, 100 - 100 * as.numeric(iso_levels))
      names(iso_label) <- iso_levels
      
      tab_measured <- dt_iso[, .SD, .SDcols = -c("pos_x_ref", "pos_y_ref")]
      tab_measured[, iso_source := "measured"]
      tab_expected <- dt_iso[, .SD, .SDcols = -c("pos_x", "pos_y")]
      tab_expected[, iso_source := "expected"]
      data.table::setnames(tab_expected, old = c("pos_x_ref", "pos_y_ref"), new = c("pos_x", "pos_y"))
      
      tab_isoline <- rbind(tab_measured, tab_expected)
      # adjust isoline range to heatmap
      tab_isoline <- 
        tab_isoline[data.table::between(pos_x, range_xy[1], range_xy[2]) & 
                      data.table::between(pos_y, range_xy[1], range_xy[2]), ]
      
      plt <- plt +
        ggplot2::geom_path(data = tab_isoline,
                           ggplot2::aes(x = if (swap_axes) pos_y else pos_x, 
                                        y = if (swap_axes) pos_x else pos_y,
                                        linetype = iso_source, color = iso_level),
                           linewidth = 1) +
        ggplot2::scale_linetype_manual(values = c("measured" = "solid", "expected" = "twodash")) +
        ggplot2::scale_color_manual(values = iso_colors,
                                    breaks = names(iso_label),
                                    labels = iso_label) +
        ggplot2::labs(color = "Iso Levels",
                      linetype = normalization_type)
    }
  }
  
  # final plot
  plt <- plt +
    ggplot2::facet_wrap(~get(cellline_name)) +
    ggplot2::guides(fill = ggplot2::guide_colorbar(order = 1),
                    linetype = ggplot2::guide_legend(order = 2), 
                    color = ggplot2::guide_legend(order = 3)) +
    ggplot2::theme(
      legend.position = "left", 
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 10, face = "bold", hjust = 0, margin = ggplot2::margin()))
  
  return(plt)
}


#' Plot panel of heatmaps with fitted and reference data for isobolograms
#' 
#' This function is dedicated to cases in which given cell lines are exposed to drugs of different concentrations 
#' and have almost no or no common values.
#'
#' @inheritParams heatmap_combo_with_isoref_panel
#'    
#' @return \code{ggplot} object containing panel with heatmaps for fitted values and reference data 
#'    for isobolograms for selected drug and co-drug by list of cell lines
#'    
#' @keywords combo_plots
#' @examples
#' cl_names <-
#'   c("cellline_AA", "cellline_EA", "cellline_IB",
#'   "cellline_MC", "cellline_BC", "cellline_FD")
#' 
#' drug1_name <- "drug_001"
#' drug2_name <- "drug_026"
#' 
#' mae <- gDRutils::get_synthetic_data("combo_matrix")
#' se <- mae[[gDRutils::get_supported_experiments("combo")]]
#' dt_excess <- gDRutils::convert_se_assay_to_dt(se, "excess")
#' dt_isobolograms <- gDRutils::convert_se_assay_to_dt(se, "isobolograms")
#' 
#' heatmap_combo_with_isoref_panel_independent(dt_excess,
#'                                             dt_isobolograms,
#'                                             drug1_name, drug2_name,
#'                                             cl_names)
#' 
#' heatmap_combo_with_isoref_panel_independent(dt_excess,
#'                                             dt_isobolograms,
#'                                             drug1_name, drug2_name,
#'                                             cl_names,
#'                                             iso_levels = c("-0.25", "0.25"))
#'                                             
#' heatmap_combo_with_isoref_panel_independent(
#'   dt_excess,
#'   dt_isobolograms,
#'   drug1_name, drug2_name,
#'   cl_names = c("cellline_FD", "cellline_MC", "cellline_AA", "cellline_EA"),
#'   iso_levels =  c("-0.25", "-0.05", "0.2", "0.65"))
#' 
#' heatmap_combo_with_isoref_panel_independent(dt_excess,
#'                                             dt_isobolograms,
#'                                             drug1_name, drug2_name,
#'                                             cl_names,
#'                                             normalization_type = "RV",
#'                                             iso_levels = NULL,
#'                                             colors_vec = c("darkcyan", "snow", "darkorange"),
#'                                             swap_axes = TRUE)
#' 
#' @export
heatmap_combo_with_isoref_panel_independent <- function(
    dt_excess,
    dt_isobolograms,
    drug1_name,
    drug2_name,
    cl_names,
    normalization_type = "GR",
    metric = "smooth",
    iso_levels = "0.5",
    colors_vec = NULL,
    no_breaks = 50,
    swap_axes = FALSE) {
  
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  gnumber <- gDRutils::get_env_identifiers("drug")
  drug_name_2 <- gDRutils::get_env_identifiers("drug_name2")
  gnumber_2 <- gDRutils::get_env_identifiers("drug2")
  conc <- gDRutils::get_env_identifiers("concentration")
  conc_2 <- gDRutils::get_env_identifiers("concentration2")
  
  checkmate::assert_data_table(dt_excess)
  checkmate::assert_data_table(dt_isobolograms)
  checkmate::assert_string(drug1_name)
  checkmate::assert_choice(drug1_name, choices = dt_excess[[drug_name]])
  checkmate::assert_string(drug2_name)
  checkmate::assert_choice(drug2_name, choices = dt_excess[[drug_name_2]])
  checkmate::assert_character(cl_names, null.ok = TRUE)
  checkmate::assert_choice(normalization_type, choices = c("GR", "RV"))
  checkmate::assert_choice(metric, choices = names(gDRutils::get_combo_excess_field_names()))
  checkmate::assert_character(iso_levels, null.ok = TRUE)
  if (!is.null(iso_levels)) {
    stopifnot("`iso_levels` must be a valid numeric value" = 
                all(vapply(iso_levels, function(i) grepl("^[-]*0\\.?[0-9]*$", i), logical(1))))
  }
  checkmate::assert_character(colors_vec, null.ok = TRUE)
  checkmate::assert_int(no_breaks, lower = 2)
  checkmate::assert_flag(swap_axes)
  
  available_cls <- unique(dt_excess[[cellline_name]])
  if (is.null(cl_names) || all(!cl_names %in% available_cls)) {
    cl_names  <- available_cls
  } else if (!all(cl_names %in% available_cls)) {
    cl_names <- cl_names[cl_names %in% available_cls]
  }
  
  # filter data for normalization type
  filter_expr <- substitute(normalization_type == norm_type, list(norm_type = normalization_type))
  dt_excess <- dt_excess[eval(filter_expr)]
  dt_isobolograms <- dt_isobolograms[eval(filter_expr)]
  
  # filter data for combination cell line (drug x drug2)
  selected_combination <-
    unique(dt_excess[get(cellline_name) %in% cl_names & get(drug_name) == drug1_name & get(drug_name_2) == drug2_name, 
                     .SD, .SDcols = c(cellline_name, drug_name, drug_name_2)])
  
  dt_excess <-
    dt_excess[selected_combination, on = c(cellline_name, drug_name, drug_name_2)]
  dt_isobolograms <-
    dt_isobolograms[selected_combination, on = c(cellline_name, drug_name, drug_name_2)]
  
  # removing cell line with no data
  cl_names_with_data <- cl_names[cl_names %in% unique(dt_excess[get(cellline_name) %in% cl_names][[cellline_name]])]
  
  # panel title
  panel_title <- sprintf("%s (%s) x %s (%s)",
                         drug1_name,
                         unique(dt_excess[get(drug_name) == drug1_name, ][[gnumber]]),
                         drug2_name,
                         unique(dt_excess[get(drug_name_2) == drug2_name, ][[gnumber_2]]))
  
  # set consisten colors accros panels
  colors_vec_iso <- .get_iso_colors(sort(unique(dt_isobolograms[iso_level %in% iso_levels][["iso_level"]])))
  
  plt_list <- lapply(cl_names_with_data, function(cl_nm) {
    plt <- heatmap_combo_with_isoref(
      dt_excess = dt_excess,
      dt_isobolograms = dt_isobolograms,
      drug1_name = drug1_name,
      drug2_name = drug2_name,
      cl_name = cl_nm,
      normalization_type = normalization_type,
      metric = metric,
      iso_levels = iso_levels,
      colors_vec = colors_vec,
      colors_vec_iso = colors_vec_iso,
      no_breaks = no_breaks,
      swap_axes = swap_axes)
  })
  names(plt_list) <- cl_names_with_data
  
  # find the maximum legend
  if (!is.null(iso_levels) && any(iso_levels %in% unique(dt_isobolograms$iso_level))) {
    dt_num_iso <-
      unique(dt_isobolograms[iso_level %in% iso_levels, .SD, .SDcols = c(cellline_name, "iso_level")])
    lbl_legend <- dt_num_iso[, .N, by = cellline_name][order(N)][N == max(N), get(cellline_name)][[1]]
  } else {
    lbl_legend <- names(plt_list)[1]
  }
  
  panel <- ggpubr::annotate_figure(
    ggpubr::ggarrange(plotlist = plt_list, widths = c(1, 1),
                      common.legend = TRUE, legend.grob = ggpubr::get_legend(plt_list[[lbl_legend]]),
                      legend = "left"),
    top = panel_title)
  
  # final panel
  return(panel)
}

#' Calculate limit for combo heatmap with gDR assumptions
#'
#' @param num_vec numeric vector
#' @param metric string name of the combo metric;
#'   one of: "smooth" ("Smooth GR" or "Smooth RV" - respectively depending on \code{normalization_type})
#'   "hsa_excess" ("Bliss Excess GR" or "Bliss Excess RV") or "bliss_excess" 
#'   ("Bliss Excess GR" or "Bliss Excess RV")
#' @param normalization_type string with normalization_types to be selected
#'                           one of: "GR" ("GRvalue") or "RV" ("RelativeViability")
#' @param symmetric logical indicating if limits should be symmetric around 0
#'
#' @return capped limits (min and max) for given numeric vector
#'
#' @keywords internal
#' @examples
#' \dontrun{
#' vec <- c(-0.1, -0.3, 0, 0.5, Inf, NA)
#' prep_hm_limits(vec)
#' prep_hm_limits(vec, metric = "hsa_excess", symmetric = TRUE)
#' }
#' 
prep_hm_limits <- function(num_vec,
                           metric = "smooth",
                           normalization_type = "GR",
                           symmetric = FALSE) {
  
  checkmate::assert_numeric(num_vec)
  checkmate::assert_choice(metric, choices = names(gDRutils::get_combo_excess_field_names()))
  checkmate::assert_choice(normalization_type, choices = c("GR", "RV"))
  checkmate::assert_logical(symmetric)
  
  vec_range <- range(num_vec, na.rm = TRUE, finite = TRUE)
  min_data <- min(vec_range)
  max_data <- max(vec_range)
  
  max_val <- if (metric == "smooth") {
    if (max_data > 1) {
      max_data
    } else {
      1
    }
  } else {
    if (max_data < 0.25) {
      0.25
    } else {
      max_data
    }
  }
  
  min_val <- if (metric == "smooth") {
    if (normalization_type == "GR") {
      min(0, min_data)
    } else {
      0
    }
  } else {
    if (min_data > -0.25) {
      -0.25
    } else {
      min_data
    }
  }
  
  if (symmetric) {
    max_abs_val <- max(abs(c(min_val, max_val)))
    min_val <- -max_abs_val
    max_val <- max_abs_val
  }
  
  return(c(min_val, max_val))
}

#' Transform concentrations values with log10
#'
#' @param conc_vec numeric vector with concentration values
#'
#' @return numeric vector with log10 concentration values; log10 for concentration equal 0 (-Inf)
#'    is replaced with one step less in the dose dilution
#'
#' @keywords internal
#' @examples
#' \dontrun{
#' vec <- c(0, 0.003, 0.01, 0.03)
#' transform_log_conc(vec)
#' }
#' 
transform_log_conc <- function(conc_vec) {
  checkmate::assert_numeric(conc_vec, lower = 0, any.missing = FALSE, finite = TRUE)
  
  log_values <- log10(conc_vec)
  # replace the -Inf value coming from the 0 dose with one step less in the dose dilution
  idx_inf <- (conc_vec == 0)
  doses <- sort(unique(log_values[!is.infinite(log_values)]))
  zero_replacement <- doses[1] + (doses[1] - doses[2]) 
  if (is.na(zero_replacement)) zero_replacement <- doses[1] - 0.5 # only two unique conc and one is 0
  log_values[idx_inf] <- zero_replacement
  # final
  log_values
}

#' Calculate the size of tiles based on pos_x/pos_y values
#'
#' Since \code{ggplot2::geom_tile} uses the center of the tile and its size (x, y, width, height),
#' x and y are given as pos_x and pos_y, it is required to calculate the width and height.
#'
#' @param pos_vec numeric vector with pos_x or pos_y values
#'
#' @return size of tile in \code{ggplot2::geom_tile}
#'
#' @keywords internal
.get_tile_size <- function(pos_vec) {
  checkmate::assert_numeric(pos_vec)
  
  diff_ <- sort(unique(diff(sort(unique(round(pos_vec, 4))))), decreasing = TRUE)
  
  tile_size <- if (NROW(diff_) > 1) {
    max(diff_)
  } else if (NROW(diff_) == 1) {
    diff_
  } else { 
    0.5
  }
  tile_size
}

#' Get color palette for the isobologram levels
#' 
#' @param iso_levels character vector with isobologram levels
#' 
#' @return gDR palette for isoline given in \code{iso_levels}
#' 
#' @keywords internal
#' @examples
#' \dontrun{
#' ls_iso_lvl <- c("0.25", "0.5", "0.75")
#' .get_iso_colors(ls_iso_lvl)
#' }
.get_iso_colors <- function(iso_levels) {
  checkmate::assert_character(iso_levels)
  
  # order iso level
  iso_levels <- iso_levels[order(as.numeric(iso_levels))]
  
  iso_colors <- 
    grDevices::colorRampPalette(
      gDRutils::get_settings_from_json("ISOLINE_PALETTE",
                                       system.file(package = "gDRplots", "settings.json"))
    )(2 * NROW(iso_levels))[2 * seq_along(iso_levels)]
  names(iso_colors) <- iso_levels
  
  iso_colors
}


#' Get color palette for the smooth values
#' 
#' @param no_breaks numeric number of breaks on scale
#' 
#' @return gDR palette for smooth values with given \code{no_breaks}
#' 
#' @keywords internal
#' @examples
#' \dontrun{
#' .get_smooth_palette(25)
#' }
.get_smooth_palette <- function(no_breaks) {
  checkmate::assert_int(no_breaks, lower = 2)
  
  grDevices::colorRampPalette(
    gDRutils::get_settings_from_json("SMOOTH_PALETTE",
                                     system.file(package = "gDRplots", "settings.json"))
  )(no_breaks)
}

#' Get color palette for the excess values
#' 
#' @param no_breaks numeric number of breaks on scale
#' 
#' @return gDR palette for excess values with given \code{no_breaks}
#' 
#' @keywords internal
#' @examples
#' \dontrun{
#' .get_excess_palette(20)
#' }
.get_excess_palette <- function(no_breaks) {
  checkmate::assert_int(no_breaks, lower = 2)
  
  grDevices::colorRampPalette(
    gDRutils::get_settings_from_json("EXCESS_PALETTE",
                                     system.file(package = "gDRplots", "settings.json"))
  )(no_breaks)
}


#' Check type of concentrations set for combination of drug per cell line
#' 
#' This function checks if the concentration vectors per cell line have common part or not.
#' It is required to decide which function to use for plotting the heatmaps panel (the set of heatmaps
#' with combo metrics plot for combination of selected drug with selected codrug, each heatmap per one cell line):
#' \code{\link{heatmap_combo_with_isoref_panel_common}} that plot heatmaps with shared axes
#' or \code{\link{heatmap_combo_with_isoref_panel_independent}} that plot heatmaps independently.
#' 
#' Possible combinations (0 is not taken into account as it is always present):
#' \itemize{
#'   \item \code{common}
#'     \itemize{
#'       \item all vectors have common part and start and end conc are the same\cr
#'             0, 0.03, 0.1, 0.3, 1\cr
#'       \item all vectors have common part and start and end conc are the same, but there are some gap\cr
#'             0, 0.003, 0.01, 0.03, 0.1, 0.3, 1\cr
#'             0, 0.003, 0.01, 0.03, ___, 0.3, 1\cr
#'          }
#'   \item \code{independent}
#'     \itemize{
#'       \item all vectors have common part and start conc are the same, but end conc is different (no gap)\cr
#'             0, 0.003, 0.01, 0.03, 0.1, 0.3, __\cr
#'             0, 0.003, 0.01, 0.03, 0.1, 0.3, 1\cr
#'       \item no common part between vectors\cr
#'             0, 0.003, 0.01, 0.03, __, __, __\cr
#'             0, ____,  ____, ____, 0.1, 0.3, 1\cr
#'       \item all vectors have common part but start and end conc are different (shifted range)\cr
#'             0, 0.003, 0.01, 0.03, 0.1, 0.3, 1, __, __\cr
#'             0, ____,  ____, 0.03, 0.1, 0.3, 1, 3, 10\cr
#'          }
#' }
#' 
#' @param ls_vec_conc a list with vectors with concentration per cell line
#'
#' @return a string decribing type of concentration list - one of:
#' \itemize{
#'   \item \code{common} vectors in the input list have common part; heatmaps in the panel
#'   should be created jointly with function \code{\link{heatmap_combo_with_isoref_panel_common}}
#'   \item \code{independent} vectors in the input list do not have common part; heatmaps in the panel
#'   should be created independently with function \code{\link{heatmap_combo_with_isoref_panel_independent}}
#' }
#'
#' @keywords internal
#' 
#' @author Janina Smoła \email{janina.smola@@contractors.roche.com}
#' 
#' @examples
#' \dontrun{
#' ls_conc <- list(c(0, 0.003, 0.01, 0.03), c(0, 0.003, 0.01, 0.03, 0.1))
#' .get_combo_panel_type(ls_conc)
#' }
#' 
.get_combo_panel_type <- function(ls_vec_conc) {
  checkmate::assert_list(ls_vec_conc)
  stopifnot("Must be a list with numeric vectors." = all(
    vapply(ls_vec_conc, function(x) is.numeric(x), logical(1))
  ))
  
  # find unique vectors
  ls_vec_conc <- ls_vec_conc[!duplicated(lapply(ls_vec_conc, sort))]
  ls_vec_conc <- ls_vec_conc[vapply(ls_vec_conc, function(x) NROW(x) > 0, logical(1))]
  
  if (NROW(ls_vec_conc) == 1) {
    # there is only one concentration set
    return("common")
  } else {
    # clean list 
    digits <- max(vapply(ls_vec_conc, function(x) { 
      max(nchar(as.character(x))) }, numeric(1)))
    # clean unique vectors
    ls_vec_conc_clean <- lapply(ls_vec_conc, function(x) {
      x <- x[!is.na(x) & x > 0] # remove NAs & 0
      x <- round(as.numeric(x), digits) # handle floating point differences
      x
    })
    # common part of vectors
    common_conc <- Reduce(intersect, ls_vec_conc_clean) 
    # all posible values and start&end conditions
    all_conc <- unique(unlist(ls_vec_conc_clean))
    ls_range_conc <- lapply(ls_vec_conc_clean, function(x) which(all_conc %in% x))
    start_is_same <- NROW(unique(vapply(ls_vec_conc_clean, function(x) min(x), numeric(1)))) == 1
    end_is_same <- NROW(unique(vapply(ls_vec_conc_clean, function(x) max(x), numeric(1)))) == 1
    
    if (NROW(common_conc) == 0 || 
        (all(vapply(ls_vec_conc_clean, function(x) NROW(setdiff(x, common_conc)) > 0, logical(1))) &&
         !(start_is_same && end_is_same))) {
      # each vector is fully independent or the shift between is too big
      return("independent")
    } else {
      return("common")
    }
  }
}
