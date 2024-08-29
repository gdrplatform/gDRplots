#' Plot heatmaps of fitted values for combination metrics data
#'
#' @param dt_excess data.table data.table representing data from the \code{excess} assay,
#'    outputted by \code{gDRutils::convert_se_assay_to_dt(se, "excess")}
#' @param dt_isobolograms data.table data.table representing data from the \code{isobolograms} assay,
#'    outputted by \code{gDRutils::convert_se_assay_to_dt(se, "isobolograms")}
#' @param drug1_name string with drug name to be plotted (identifiers \code{DrugName})
#' @param drug2_name string with co-drug name to be plotted (identifiers \code{DrugName_2})
#' @param cl_name string with cell line to be plotted (identifiers \code{CellLineName})
#' @param normalization_type string with normalization_types to be selected
#'                           one of: "GR" ("GRvalue") or "RV" ("RelativeViability")
#' @param iso_levels character vector with  isobologram levels to be selected
#' @param colors_vec_smooth character vector of colors (valid names or hex codes) used in the heatmap
#'    for smooth values; the default is the viridis palette
#' @param colors_vec_excess character vector of colors (valid name or hex codes) used in the heatmap
#'    for excess values; the default is a blue-light grey-red color scale
#' @param no_breaks numeric number of breaks on scale
#' @param as_panel logical flag whether return list of plot or panel
#'
#' @return list or panel with heatmaps with value for excess assays for selected drugs and cell line with
#'    selected isoline and comparison of iso levels
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
#' cl_name <- "cellline_JE"
#' drug1_name <- "drug_011"
#' drug2_name <- "drug_026"
#' 
#' heatmap_combo_metrics(dt_excess,
#'                       dt_isobolograms,
#'                       drug1_name, drug2_name,
#'                       cl_name,
#'                       normalization_type = "RV",
#'                       iso_levels = "0.5",
#'                       as_panel = FALSE)
#'
#' @export
heatmap_combo_metrics <- function(
    dt_excess,
    dt_isobolograms,
    drug1_name,
    drug2_name,
    cl_name,
    normalization_type = "GR",
    iso_levels =  c("0.25", "0.5", "0.75"),
    colors_vec_smooth = NULL,
    colors_vec_excess = NULL,
    no_breaks = 50,
    as_panel = TRUE) {
  
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  clid <- gDRutils::get_env_identifiers("cellline")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  drug_name_2 <- gDRutils::get_env_identifiers("drug_name2")
  conc <- gDRutils::get_env_identifiers("concentration")
  conc_2 <- gDRutils::get_env_identifiers("concentration2")
  duration <- gDRutils::get_env_identifiers("duration")
  mx_names <- names(gDRutils::get_combo_excess_field_names())
  
  checkmate::assert_data_table(dt_excess)
  checkmate::assert_data_table(dt_isobolograms)
  checkmate::assert_string(drug1_name)
  checkmate::assert_choice(drug1_name, choices = dt_excess[[drug_name]])
  checkmate::assert_string(drug2_name)
  checkmate::assert_choice(drug2_name, choices = dt_excess[[drug_name_2]])
  checkmate::assert_string(cl_name)
  checkmate::assert_choice(cl_name, choices = dt_excess[[cellline_name]])
  checkmate::assert_choice(normalization_type, choices = c("GR", "RV"))
  checkmate::assert_character(iso_levels, null.ok = TRUE)
  if (!is.null(iso_levels)) checkmate::assert_numeric(as.numeric(iso_levels))
  stopifnot("Must be a valid color name" = all(vapply(colors_vec_smooth, is_valid_color, logical(1))))
  stopifnot("Must be a valid color name" = all(vapply(colors_vec_excess, is_valid_color, logical(1))))
  checkmate::assert_int(no_breaks, lower = 2)
  checkmate::assert_flag(as_panel)
  
  # filter data for normalization type
  filter_expr <- substitute(normalization_type == norm_type, list(norm_type = normalization_type))
  dt_excess <- dt_excess[eval(filter_expr)]
  dt_isobolograms <- dt_isobolograms[eval(filter_expr)]
  
  # filter data for combination cell line (drug x drug2)
  dt_excess <-
    dt_excess[get(cellline_name) == cl_name & get(drug_name) == drug1_name & get(drug_name_2) == drug2_name]
  dt_isobolograms <-
    dt_isobolograms[get(cellline_name) == cl_name & get(drug_name) == drug1_name & get(drug_name_2) == drug2_name]
  
  # isoline data
  if (!is.null(dt_isobolograms$iso_level)) {
    dt_isobolograms <- dt_isobolograms[iso_level %in% iso_levels, ]
  }
  available_iso_lvl <- unique(dt_isobolograms[["iso_level"]])
  iso_colors <- get_iso_colors(normalization_type)[available_iso_lvl]
  
  # title
  main_title <- sprintf("%s (%s)",
                        cl_name,
                        unique(dt_excess[get(cellline_name) == cl_name][[clid]]))
  # legend
  legend_title_iso <- "Iso Levels"
  legend_lbl_iso <- NULL # due to NSE notes in R CMD check
  legend_lbl_iso <- paste0(ifelse(normalization_type == "GR", "GR", "IC"),
                           100 - 100 * as.numeric(available_iso_lvl))
  
  # prep hm color palette
  hm_color_palette_smooth <- if (is.null(colors_vec_smooth)) {
    grDevices::colorRampPalette(
      c("#510046", "#b3009a", "#e400c4", "#F2F2F2"))(no_breaks + 1)
  } else {
    grDevices::colorRampPalette(colors_vec_smooth)(no_breaks + 1)
  }
  
  hm_color_palette_excess <- if (is.null(colors_vec_excess)) {
    grDevices::colorRampPalette(
      c("royalblue3", "royalblue1", "grey95", "grey95", "firebrick1", "firebrick3"))(no_breaks + 1)
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
    if (!as_panel) plt_title <- paste(main_title, plt_title, sep = " : ")
    
    dt_ <- dt_excess[!is.na(get(mx_name)), c(conc, conc_2, mx_name), with = FALSE]
    
    if (!NROW(dt_) > 1) { # co-dilution input data is like: (conc = 0, conc_2 = 0, mx_name = 1)
      plt <- 
        ggplot2::ggplot() +
        ggplot2::labs(x = bquote(.(drug2_name) ~ "[" ~ mu * M ~ "]"),
                      y = bquote(.(drug1_name) ~ "[" ~ mu * M ~ "]"),
                      title = plt_title) +
        ggplot2::theme_bw() +
        ggplot2::theme(aspect.ratio = 1)
    } else {
      dt_[[mx_name]] <- pmin(1.1, dt_[[mx_name]])
      dt_$pos_y <- transform_log_conc(dt_[[conc]])
      dt_$pos_x <- transform_log_conc(dt_[[conc_2]])
      
      ls_axes <- gDRutils::define_matrix_grid_positions(dt_[[conc]], dt_[[conc_2]])
      drug1_axis <- ls_axes$axis_1
      drug2_axis <- ls_axes$axis_2
      tile_height <- .get_tile_size(drug1_axis$pos_y)
      tile_width <- .get_tile_size(drug2_axis$pos_x)
      
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
                               normalization_type = normalization_type)
      
      # base plot
      plt <-
        ggplot2::ggplot(dt_, ggplot2::aes(x = pos_x, y = pos_y)) +
        ggplot2::geom_tile(ggplot2::aes(fill = get(mx_name)), height = tile_height, width = tile_width) +
        ggplot2::labs(x = bquote(.(drug2_name) ~ "[" ~ mu * M ~ "]"),
                      y = bquote(.(drug1_name) ~ "[" ~ mu * M ~ "]"),
                      title = plt_title,
                      fill = legend_title_fill)  +
        ggplot2::scale_x_continuous(breaks = drug2_axis$pos_x,
                                    labels = drug2_axis$marks_x,
                                    expand = c(0, 0)) +
        ggplot2::scale_y_continuous(breaks = drug1_axis$pos_y,
                                    labels = drug1_axis$marks_y,
                                    expand = c(0, 0)) +
        ggplot2::scale_fill_gradientn(colors = hm_color_palette,
                                      limit = limits,
                                      labels = function(x) sprintf("%.2f", x)) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
                       axis.text.y = ggplot2::element_text(size = 8),
                       plot.title = ggplot2::element_text(size = 10),
                       panel.grid = ggplot2::element_blank(),
                       aspect.ratio = 1)
      
      # add isoline
      if (NROW(available_iso_lvl)) { # isobolograms as lines
        if (all(available_iso_lvl %in% c("0.25", "0.5", "0.75"))) {
          # friendly for user with color vision deficiency
          plt <- plt +
            ggplot2::geom_path(data = dt_isobolograms, linewidth = 1,
                               ggplot2::aes(x = pos_x, y = pos_y, color = iso_level, linetype = iso_level)) +
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
            ggplot2::geom_path(data = dt_isobolograms, linewidth = 1,
                               ggplot2::aes(x = pos_x, y = pos_y, color = iso_level)) +
            ggplot2::scale_color_manual(values = iso_colors[available_iso_lvl],
                                        breaks = available_iso_lvl,
                                        labels = legend_lbl_iso,
                                        name = legend_title_iso)
        }
      }
    }
    plt
  })
  names(mx_plts) <- mx_names
  
  # isobolograms across range of concentration ratios
  plt_title <- sprintf("%s for %s, T=%sh",
                       gDRutils::prettify_flat_metrics(x = "smooth", human_readable = TRUE),
                       normalization_type,
                       unique(dt_excess[get(cellline_name) == cl_name][[duration]]))
  
  if (!as_panel) plt_title <- paste(main_title, plt_title, sep = " : ")
  # base plot
  plt_iso_compare <-
    ggplot2::ggplot(mapping = ggplot2::aes(x = log10_ratio_conc, y = log2_CI)) +
    ggplot2::geom_line(
      data = data.table::data.table(log10_ratio_conc = c(-2, 2), log2_CI = c(0, 0))) +
    ggplot2::geom_hline(yintercept = 0, color = "#A9A9A9")
  
  # add isoline
  if (NROW(available_iso_lvl)) { # isobolograms as lines
    if (all(available_iso_lvl %in% c("0.25", "0.5", "0.75"))) {
      # friendly for user with color vision deficiency
      plt_iso_compare <- plt_iso_compare +
        ggplot2::geom_path(data = dt_isobolograms, linewidth = 0.5,
                           ggplot2::aes(x = log10_ratio_conc, y = log2_CI, color = iso_level, linetype = iso_level)) +
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
      plt_iso_compare <- plt_iso_compare +
        ggplot2::geom_path(data = dt_isobolograms, linewidth = 0.5,
                           ggplot2::aes(x = log10_ratio_conc, y = log2_CI, color = iso_level))
      ggplot2::scale_color_manual(values = iso_colors[available_iso_lvl],
                                  breaks = available_iso_lvl,
                                  labels = legend_lbl_iso,
                                  name = legend_title_iso)
    }
  }
  
  # add x and y scales
  plt_iso_compare <- plt_iso_compare +
    ggplot2::scale_y_continuous(breaks = -5:4, labels = c(paste0("1/", 2 ^ (5:1)), 2 ^ (0:4))) +
    ggplot2::scale_x_continuous(breaks = -3:3, labels = c(paste0("1/", 10 ^ (3:1)), 10 ^ (0:3))) +
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
  
  final_plot <- if (as_panel) {
    ggpubr::annotate_figure(
      ggpubr::ggarrange(
        ggpubr::ggarrange(
          plotlist = list(
            ls_plts[["smooth"]] + ggplot2::guides(linetype = "none", color = "none"),
            ls_plts[["iso_compare"]] + ggplot2::guides(linetype = "none", color = "none")
          ),
          ncol = 2, common.legend = TRUE, legend = "right"),
        ggpubr::ggarrange(
          plotlist = list(
            ls_plts[["hsa_excess"]] + ggplot2::labs(fill = "Excess"),
            ls_plts[["bliss_excess"]] + ggplot2::labs(fill = "Excess")
          ),
          ncol = 2, common.legend = TRUE, legend = "right"),
        common.legend = TRUE, nrow = 2),
      top = main_title) +
      ggpubr::bgcolor("white") + ggpubr::border("white")
  } else {
    ls_plts
  }
  return(final_plot)
}

#' Plot heatmaps of averaged values for combination data
#'
#' @param dt_excess data.table representing data from the \code{excess} assay,
#'    outputted by \code{gDRutils::convert_se_assay_to_dt(se, "excess")}
#' @param dt_isobolograms data.table representing data from the \code{isobolograms} assay,
#'    outputted by \code{gDRutils::convert_se_assay_to_dt(se, "isobolograms")}
#' @param drug1_name string with drug name to be plotted (identifiers \code{DrugName})
#' @param drug2_name string with co-drug name to be plotted (identifiers \code{DrugName_2})
#' @param cl_name string with cell line to be plotted (identifiers \code{CellLineName})
#' @param normalization_type string with normalization_types to be selected
#'                           one of: "GR" ("GRvalue") or "RV" ("RelativeViability")
#' @param iso_levels character vector with  isobologram levels to be selected;
#'     when \code{NULL} - no isolines will be displayed
#' @param colors_vec character vector of colors (valid name or hex) used in heatmap
#' @param no_breaks numeric number of breaks on scale
#'
#' @return list or panel with heatmaps with values for excess assays for selected drugs and cell line with
#'    selected isoline and comparison of iso levels
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
#'                           iso_levels = NULL,
#'                           colors_vec = c("darkcyan", "snow", "darkorange"))                      
#'                           
#' heatmap_combo_with_isoref(dt_excess,
#'                           dt_isobolograms,
#'                           drug1_name, drug2_name,
#'                           cl_name,
#'                           normalization_type = "RV",
#'                           iso_levels = c("0.25", "0.75"))
#'
#' @export
heatmap_combo_with_isoref <- function(
    dt_excess,
    dt_isobolograms,
    drug1_name,
    drug2_name,
    cl_name,
    normalization_type = "GR",
    iso_levels = "0.5",
    colors_vec = NULL,
    no_breaks = 50) {
  
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
  checkmate::assert_character(iso_levels, null.ok = TRUE)
  if (!is.null(iso_levels)) checkmate::assert_numeric(as.numeric(iso_levels))
  checkmate::assert_character(colors_vec, null.ok = TRUE)
  if (!is.null(colors_vec)) {
    stopifnot("Must be a valid color name" = all(vapply(colors_vec, is_valid_color, logical(1))))
  }
  checkmate::assert_int(no_breaks, lower = 2)
  
  # filter data for normalization type
  filter_expr <- substitute(normalization_type == norm_type, list(norm_type = normalization_type))
  dt_excess <- dt_excess[eval(filter_expr)]
  dt_isobolograms <- dt_isobolograms[eval(filter_expr)]
  
  # filter data for combination cell line (drug x drug2)
  dt_excess <-
    dt_excess[get(cellline_name) == cl_name & get(drug_name) == drug1_name & get(drug_name_2) == drug2_name]
  dt_isobolograms <-
    dt_isobolograms[get(cellline_name) == cl_name & get(drug_name) == drug1_name & get(drug_name_2) == drug2_name]
  
  # isoline data
  if (!is.null(dt_isobolograms$iso_level)) {
    dt_isobolograms <- dt_isobolograms[iso_level %in% iso_levels, ]
  }
  available_iso_lvl <- unique(dt_isobolograms[["iso_level"]])
  iso_colors <- get_iso_colors()[available_iso_lvl]
  
  # prep hm color palette
  hm_color_palette <- if (is.null(colors_vec)) {
    colorspace::sequential_hcl(no_breaks + 1, palette = "viridis")
  } else {
    grDevices::colorRampPalette(colors_vec)(no_breaks + 1)
  }
  
  # panel title
  cl_clid <- unique(dt_excess[get(cellline_name) == cl_name, ][[clid]])
  plt_title <- sprintf("%s (%s)", cl_name, cl_clid)
  
  # prep plot data
  mx_name <- "smooth"
  dt_ <- dt_excess[, c(conc, conc_2, mx_name), with = FALSE]
  # correction of NA for conc = 0 or conc_2 = 0
  dt_[(get(conc) == 0 | get(conc_2) == 0) & is.na(get(mx_name))] <- 0
  
  if (!NROW(dt_) > 1) { # co-dilution input data is like: (conc = 0, conc_2 = 0, mx_name = 1)
    plt <- 
      ggplot2::ggplot() +
      ggplot2::labs(x = bquote(.(drug2_name) ~ "[" ~ mu * M ~ "]"),
                    y = bquote(.(drug1_name) ~ "[" ~ mu * M ~ "]"),
                    title = plt_title) +
      ggplot2::theme_bw() +
      ggplot2::theme(aspect.ratio = 1)
  } else {
    dt_[[mx_name]] <- pmin(1.1, dt_[[mx_name]])
    dt_$pos_y <- transform_log_conc(dt_[[conc]])
    dt_$pos_x <- transform_log_conc(dt_[[conc_2]])
    
    ls_axes <- gDRutils::define_matrix_grid_positions(dt_[[conc]], dt_[[conc_2]])
    drug1_axis <- ls_axes$axis_1
    drug2_axis <- ls_axes$axis_2
    tile_height <- .get_tile_size(drug1_axis$pos_y)
    tile_width <- .get_tile_size(drug2_axis$pos_x)
    
    range_x <- c(min(drug2_axis$pos_x), max(drug2_axis$pos_x) + 0.5 * tile_width)
    range_y <- c(min(drug1_axis$pos_y), max(drug1_axis$pos_y) + 0.5 * tile_height)
    
    # legend title
    legend_title_fill <- sprintf("%s %s",
                                 gDRutils::prettify_flat_metrics(x = mx_name, human_readable = TRUE),
                                 normalization_type)
    
    # prep limits
    limits <- prep_hm_limits(dt_[[mx_name]],   
                             metric = mx_name,
                             normalization_type = normalization_type)
    
    # base plot
    plt <-
      ggplot2::ggplot(dt_, ggplot2::aes(x = pos_x, y = pos_y)) +
      ggplot2::geom_tile(ggplot2::aes(fill = get(mx_name), ), 
                         height = tile_height, width = tile_width, alpha = 0.90) +
      ggplot2::labs(x = bquote(.(drug2_name) ~ "[" ~ mu * M ~ "]"),
                    y = bquote(.(drug1_name) ~ "[" ~ mu * M ~ "]"),
                    title = plt_title,
                    fill = legend_title_fill) +
      ggplot2::scale_fill_gradientn(colors = hm_color_palette,
                                    limit = limits,
                                    labels = function(x) sprintf("%.2f", x))
    
    # plot isobologram
    if (NROW(available_iso_lvl)) { # add isolines - if there are such data
      iso_label <- sprintf("%s%s",
                           ifelse(normalization_type == "GR", "GR", "IC"),
                           100 - 100 * as.numeric(available_iso_lvl))
      names(iso_label) <- available_iso_lvl
      
      iso_source <- NULL # due to NSE notes in R CMD check
      tab_measured <- dt_isobolograms[, .SD, .SDcols = -c("pos_x_ref", "pos_y_ref")]
      tab_measured[, iso_source := "measured"]
      tab_expected <- dt_isobolograms[, .SD, .SDcols = -c("pos_x", "pos_y")]
      tab_expected[, iso_source := "expected"]
      data.table::setnames(tab_expected, old = c("pos_x_ref", "pos_y_ref"), new = c("pos_x", "pos_y"))
      
      tab_isoline <- rbind(tab_measured, tab_expected)
      # adjust isoline range to heatmap
      tab_isoline <- 
        tab_isoline[data.table::between(pos_x, range_x[1], range_x[2]) & 
                      data.table::between(pos_y, range_x[1], range_x[2]), ]
      
      
      if (NROW(available_iso_lvl) == 1) {
        plt <- plt +
          ggplot2::geom_path(data = tab_isoline,
                             ggplot2::aes(x = pos_x, y = pos_y, linetype = iso_source),
                             linewidth = 1, color = "red") +
          ggplot2::scale_linetype_manual(values = c("measured" = "solid", "expected" = "dashed"),
                                         name = iso_label)
      } else {
        iso_colors <-
          grDevices::colorRampPalette(c("red", "darkred"))(2 * NROW(available_iso_lvl))[seq_along(available_iso_lvl) * 2] # nolint
        names(iso_colors) <- available_iso_lvl
        plt <- plt +
          ggplot2::geom_path(data = tab_isoline,
                             ggplot2::aes(x = pos_x, y = pos_y, linetype = iso_source, color = iso_level),
                             linewidth = 1) +
          ggplot2::scale_linetype_manual(values = c("measured" = "solid", "expected" = "dashed"),
                                         name = normalization_type) +
          ggplot2::scale_color_manual(values = iso_colors,
                                      label = iso_label,
                                      breaks = available_iso_lvl,
                                      name = "Iso Levels")
      }
    }
    
    # final plot
    plt <- plt +
      ggplot2::scale_x_continuous(breaks = drug2_axis$pos_x,
                                  labels = drug2_axis$marks_x,
                                  expand = c(0, 0)) +
      ggplot2::scale_y_continuous(breaks = drug1_axis$pos_y,
                                  labels = drug1_axis$marks_y,
                                  expand = c(0, 0)) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
                     axis.text.y = ggplot2::element_text(size = 8),
                     plot.title = ggplot2::element_text(size = 10),
                     panel.grid.minor = ggplot2::element_blank(),
                     legend.key.width = ggplot2::unit(2, "line"),
                     aspect.ratio = 1)
  }
  
  return(plt)
}


#' Plot panel of heatmaps with fitted and reference data for isobolograms
#' to control quality of the data
#'
#' @inheritParams heatmap_combo_with_isoref
#' @param cl_names character vector with cell line names to be plotted (Cell Line Name);
#'    if \code{NULL} - all available cell lines will be plotted
#'    
#' @return panel with heatmaps for fitted values and reference data for isobolograms
#'    for selected drug and co-drug by cell line names
#'
#' @keywords QC_plot
#' @examples
#' cl_names <- 
#'   c("cellline_AA", "cellline_EA", "cellline_IB", "cellline_MC", "cellline_BC", "cellline_FD")
#' 
#' drug1_name <- "drug_001"
#' drug2_name <- "drug_026"
#' 
#' mae <- gDRutils::get_synthetic_data("combo_matrix")
#' se <- mae[[gDRutils::get_supported_experiments("combo")]]
#' dt_excess <- gDRutils::convert_se_assay_to_dt(se, "excess")
#' dt_isobolograms <- gDRutils::convert_se_assay_to_dt(se, "isobolograms")
#' 
#' heatmap_combo_with_isoref_qc_panel(dt_excess,
#'                                    dt_isobolograms,
#'                                    drug1_name, drug2_name,
#'                                    cl_names)
#' 
#' heatmap_combo_with_isoref_qc_panel(dt_excess,
#'                                    dt_isobolograms,
#'                                    drug1_name, drug2_name,
#'                                    cl_names,
#'                                    iso_levels = c("0.25", "0.5"))
#'                                    
#' heatmap_combo_with_isoref_qc_panel(dt_excess,
#'                                    dt_isobolograms,
#'                                    drug1_name, drug2_name,
#'                                    cl_names,
#'                                    normalization_type = "RV",
#'                                    iso_levels = NULL,
#'                                    colors_vec = c("darkcyan", "snow", "darkorange"))
#' 
#' @export
heatmap_combo_with_isoref_qc_panel <- function(
    dt_excess,
    dt_isobolograms,
    drug1_name,
    drug2_name,
    cl_names,
    normalization_type = "GR",
    iso_levels = "0.5",
    colors_vec = NULL,
    no_breaks = 50) {
  
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
  checkmate::assert_character(iso_levels, null.ok = TRUE)
  if (!is.null(iso_levels)) {
    stopifnot("`iso_levels` must be a valid numeric value" = 
                all(vapply(iso_levels, function(i) grepl("^0\\.?[0-9]*$", i), logical(1))))
  }
  checkmate::assert_character(colors_vec, null.ok = TRUE)
  if (!is.null(colors_vec)) {
    stopifnot("`colors_vec` must be a valid color name" = all(vapply(colors_vec, is_valid_color, logical(1))))
  }
  checkmate::assert_int(no_breaks, lower = 2)
  
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
  selecteted_combination <-
    unique(dt_excess[get(cellline_name) %in% cl_names & get(drug_name) == drug1_name & get(drug_name_2) == drug2_name, 
                     .SD, .SDcols = c(cellline_name, drug_name, drug_name_2)])
  
  dt_excess <-
    dt_excess[selecteted_combination, on = c(cellline_name, drug_name, drug_name_2)]
  dt_isobolograms <-
    dt_isobolograms[selecteted_combination, on = c(cellline_name, drug_name, drug_name_2)]
  
  # prep hm color palette
  hm_color_palette <- if (is.null(colors_vec)) {
    colorspace::sequential_hcl(no_breaks + 1, palette = "viridis")
  } else {
    grDevices::colorRampPalette(colors_vec)(no_breaks + 1)
  }
  
  # prep panel elements
  mx_name <- "smooth"
  # prep plot data
  dt_all <- dt_excess[, c(cellline_name, conc, conc_2, mx_name), with = FALSE]
  # correction of NA for conc = 0 or conc_2 = 0
  dt_all[(get(conc) == 0 | get(conc_2) == 0) & is.na(get(mx_name))] <- 0
  
  # prep data for heatmat
  dt_tile <- dt_all[get(cellline_name) %in% cl_names, ][, 
                                                        `:=`(
                                                          mx_name = pmin(1.1, get(mx_name)),
                                                          pos_y = transform_log_conc(get(conc)),
                                                          pos_x = transform_log_conc(get(conc_2))
                                                        ), 
                                                        by = cellline_name
  ][, .SD, .SDcols = -mx_name]
  data.table::setnames(dt_tile, "mx_name", mx_name)
  
  # tiles positioning 
  ls_axes_all <- gDRutils::define_matrix_grid_positions(dt_tile[[conc]], dt_all[[conc_2]])
  drug1_axis_all <- ls_axes_all$axis_1
  drug2_axis_all <- ls_axes_all$axis_2
  tile_height <- .get_tile_size(drug1_axis_all$pos_y)
  tile_width <- .get_tile_size(drug2_axis_all$pos_x)
  
  range_x <- c(min(drug2_axis_all$pos_x), max(drug2_axis_all$pos_x) + tile_width)
  range_y <- c(min(drug1_axis_all$pos_y), max(drug1_axis_all$pos_y) + tile_height)
  
  # prep limits
  limits <- prep_hm_limits(dt_tile[[mx_name]],   
                           metric = mx_name,
                           normalization_type = normalization_type)
  # legend title
  legend_title_fill <- sprintf("%s %s",
                               gDRutils::prettify_flat_metrics(x = mx_name, human_readable = TRUE),
                               normalization_type)
  # base plot
  plt <-
    ggplot2::ggplot(dt_tile, ggplot2::aes(x = pos_x, y = pos_y)) +
    ggplot2::geom_tile(ggplot2::aes(fill = get(mx_name), ), 
                       height = tile_height, width = tile_width, alpha = 0.90) +
    ggplot2::labs(x = bquote(.(drug2_name) ~ "[" ~ mu * M ~ "]"),
                  y = bquote(.(drug1_name) ~ "[" ~ mu * M ~ "]"),
                  title = panel_title,
                  fill = legend_title_fill) +
    ggplot2::scale_fill_gradientn(colors = hm_color_palette,
                                  limit = limits,
                                  labels = function(x) sprintf("%.2f", x))
  
  # isoline data
  if (!is.null(dt_isobolograms$iso_level) && !is.null(iso_levels)) {
    # iso level availability
    available_iso_lvl <- unique(dt_isobolograms[["iso_level"]])
    iso_levels <- iso_levels[iso_levels %in% available_iso_lvl]
    
    if (NROW(iso_levels)) {
      # order iso level
      iso_levels <- iso_levels[order(as.numeric(iso_levels))]
      
      req_cols <- c(cellline_name, drug_name, drug_name_2, gDRutils::get_header("iso_position"))
      dt_iso <- 
        dt_isobolograms[iso_level %in% iso_levels, .SD, .SDcols = req_cols]
      
      # colors for isoline
      iso_colors <- if (NROW(iso_levels) == 1) {
        "red"
      } else {
        grDevices::colorRampPalette(c("red", "darkred"))(2 * NROW(iso_levels))[2 * seq_along(iso_levels)] # nolint
      }
      names(iso_colors) <- iso_levels 
      
      # plot
      iso_label <- sprintf("%s%s",
                           ifelse(normalization_type == "GR", "GR", "IC"),
                           100 - 100 * as.numeric(iso_levels))
      names(iso_label) <- iso_levels
      
      iso_source <- NULL # due to NSE notes in R CMD check
      tab_measured <- dt_iso[, .SD, .SDcols = -c("pos_x_ref", "pos_y_ref")]
      tab_measured[, iso_source := "measured"]
      tab_expected <- dt_iso[, .SD, .SDcols = -c("pos_x", "pos_y")]
      tab_expected[, iso_source := "expected"]
      data.table::setnames(tab_expected, old = c("pos_x_ref", "pos_y_ref"), new = c("pos_x", "pos_y"))
      
      tab_isoline <- rbind(tab_measured, tab_expected)
      # adjust isoline range to heatmap
      tab_isoline <- 
        tab_isoline[data.table::between(pos_x, range_x[1], range_x[2]) & 
                      data.table::between(pos_y, range_x[1], range_x[2]), ]
      
      plt <- plt +
        ggplot2::geom_path(data = tab_isoline,
                           ggplot2::aes(x = pos_x, y = pos_y, linetype = iso_source, color = iso_level),
                           linewidth = 1) +
        ggplot2::scale_linetype_manual(values = c("measured" = "solid", "expected" = "dashed"),
                                       name = normalization_type) +
        ggplot2::scale_color_manual(values = iso_colors,
                                    label = iso_label,
                                    name = "Iso Levels")
    }
  }
  
  # final plot
  plt <- plt +
    ggplot2::scale_x_continuous(breaks = drug2_axis_all$pos_x,
                                labels = drug2_axis_all$marks_x,
                                expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = drug1_axis_all$pos_y,
                                labels = drug1_axis_all$marks_y,
                                expand = c(0, 0)) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
                   axis.text.y = ggplot2::element_text(size = 8),
                   plot.title = ggplot2::element_text(size = 10),
                   panel.grid.minor = ggplot2::element_blank(),
                   legend.key.width = ggplot2::unit(2, "line"),
                   legend.title = ggplot2::element_text(size = 8),
                   aspect.ratio = 1) +
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

#' Calculate limit for combo heatmap with gDR assumptions
#'
#' @param num_vec numeric vector
#' @param metric string name of combo exccess metric;
#'    one of: "smooth", "hsa_excess", "bliss_excess"
#' @param normalization_type string with normalization_types to be selected
#'                           one of: "GR" ("GRvalue") or "RV" ("RelativeViability")
#'
#' @return capped limits (min and max) for given numeric vector
#'
#' @keywords internal
#' @examples
#' \dontrun{
#' vec <- c(-0.1, -0.3, 0, 0.5, Inf, NA)
#' prep_hm_limits(vec)
#' prep_hm_limits(vec, metric = "hsa_excess")
#' }
#' 
prep_hm_limits <- function(num_vec,
                           metric = "smooth",
                           normalization_type = "GR") {
  
  checkmate::assert_numeric(num_vec)
  checkmate::assert_choice(metric, choices = names(gDRutils::get_combo_excess_field_names()))
  checkmate::assert_choice(normalization_type, choices = c("GR", "RV"))
  
  vec_range <- range(num_vec, na.rm = TRUE, finite = TRUE)
  min_data <- min(vec_range)
  max_data <- max(vec_range)
  
  max_val <- if (metric == "smooth") {
    ifelse(max_data > 1, max_data, 1)
  } else {
    ifelse(max_data < 0.25, 0.25, max_data)
  }
  
  min_val <- if (metric == "smooth") {
    ifelse(normalization_type == "GR", min(0, min_data), 0)
  } else {
    ifelse(min_data > -0.25, -0.25, min_data)
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
  
  return(log_values)
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
  
  diff_ <- sort(unique(diff(pos_vec)), decreasing = TRUE)
  
  tile_size <- if (NROW(diff_) > 1) {
    diff_[2] # 1st in related to conc = 0 and and is not conclusive
  } else if (NROW(diff_) == 1) {
    diff_
  } else { 
    0.5
  }
  tile_size
}
