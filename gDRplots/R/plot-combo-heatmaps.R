#' Plot fitted values heatmaps for combo data for combo metrics
#'
#' @param dt_excess data.table representation of the data in \code{excess} assay
#'    output from \code{gDRutils::convert_se_assay_to_dt(se, "excess")}
#' @param dt_isobolograms data.table representation of the data in \code{isobolograms} assay
#'    output from \code{gDRutils::convert_se_assay_to_dt(se, "isobolograms")}
#' @param drug1_name string with drug name to be plotted (identifiers \code{DrugName})
#' @param drug2_name string with co-drug name to be plotted (identifiers \code{DrugName_2})
#' @param cl_name string with cell line to be plotted (identifiers \code{CellLineName}) 
#' @param normalization_type string with normalization_types to be selected
#'                           one of: "GR" ("GRvalue") or "RV" ("RelativeViability")
#' @param iso_levels character vector with  isobologram levels to be selected
#' @param colors_vec_smooth character vector of colors (valid name or hex) used in heatmap 
#'    for smooth values; as default will be used viridis pallette
#' @param colors_vec_excess character vector of colors (valid name or hex) used in heatmap 
#'    for excess values; as default will be used blue - light grey - red color scale
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
  gnumber <- gDRutils::get_env_identifiers("drug")
  drug_name_2 <- gDRutils::get_env_identifiers("drug_name2")
  gnumber_2 <- gDRutils::get_env_identifiers("drug2")
  conc <- gDRutils::get_env_identifiers("concentration")
  conc_2 <- gDRutils::get_env_identifiers("concentration2")
  duration <- gDRutils::get_env_identifiers("duration")
  mx_names <- names(gDRutils::get_combo_excess_field_names())
  
  checkmate::expect_data_table(dt_excess)
  checkmate::expect_data_table(dt_isobolograms)
  checkmate::assert_string(drug1_name)
  checkmate::assert_choice(drug1_name, choices = dt_excess[[drug_name]])
  checkmate::assert_string(drug2_name)
  checkmate::assert_choice(drug2_name, choices = dt_excess[[drug_name_2]])
  checkmate::assert_string(cl_name)
  checkmate::assert_choice(cl_name, choices = dt_excess[[cellline_name]])
  checkmate::assert_choice(normalization_type, choices = c("GR", "RV"))
  checkmate::assert_character(iso_levels, null.ok = TRUE)
  if (!is.null(iso_levels)) checkmate::assert_numeric(as.numeric(iso_levels))
  stopifnot("Must be valid color name" = all(vapply(colors_vec_smooth, gDRplots::is_valid_color, logical(1))))
  stopifnot("Must be valid color name" = all(vapply(colors_vec_excess, gDRplots::is_valid_color, logical(1))))
  checkmate::assert_int(no_breaks, lower = 2)
  checkmate::assert_flag(as_panel)
  
  # filter data for normalization type
  data.table::setkeyv(dt_excess, "normalization_type")
  dt_excess <- dt_excess[normalization_type]
  data.table::setkey(dt_excess, NULL)
  
  data.table::setkeyv(dt_isobolograms, "normalization_type")
  dt_isobolograms <- dt_isobolograms[normalization_type]
  data.table::setkey(dt_isobolograms, NULL)
  
  # filter data for combination cell line (drug x drug2)
  dt_excess <- 
    dt_excess[get(cellline_name) == cl_name & get(drug_name) == drug1_name & get(drug_name_2) == drug2_name]
  dt_isobolograms <- 
    dt_isobolograms[get(cellline_name) == cl_name & get(drug_name) == drug1_name & get(drug_name_2) == drug2_name]
  
  # colors for isoline
  iso_colors <- gDRutils::get_iso_colors()[iso_levels]
  dt_isobolograms <- dt_isobolograms[iso_level %in% iso_levels, ]
  avialable_iso <- unique(dt_isobolograms$iso_level)
  
  # title 
  main_title <- sprintf("%s (%s)",
                        cl_name,
                        unique(dt_excess[get(cellline_name) == cl_name][[clid]]))
  # legend
  legend_title_iso <- "Iso Levels"
  legend_lbl_iso <- NULL # due to NSE notes in R CMD check
  legend_lbl_iso <- paste0(ifelse(normalization_type == "GR", "GR", "IC"),
                           100 - 100 * as.numeric(avialable_iso))
  
  # prep hm color palette
  hm_color_palette_smooth <- if (is.null(colors_vec_smooth)) {
    colorspace::sequential_hcl(no_breaks + 1, palette = "viridis")
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
    dt_ <- dt_excess[, c(conc, conc_2, mx_name), with = FALSE]
    # correction of NA for conc = 0 ir conc_2 = 0
    dt_[(get(conc) == 0 | get(conc_2) == 0) & is.na(get(mx_name))] <- 0
    
    dt_[[mx_name]] <- pmin(1.1, dt_[[mx_name]])
    dt_$pos_y <- log10(dt_[[conc]])
    dt_$pos_x <- log10(dt_[[conc_2]])
    
    ls_axes <- gDRutils::define_matrix_grid_positions(dt_[[conc]], dt_[[conc_2]])
    drug1_axis <- ls_axes$axis_1
    drug2_axis <- ls_axes$axis_2
    tile_height <- diff(drug1_axis$pos_y[3:4])
    tile_width <- diff(drug2_axis$pos_x[3:4])
    
    # prep hm color palette
    hm_color_palette <- if (mx_name == "smooth") {
      hm_color_palette_smooth
    } else {
      hm_color_palette_excess
    }
    
    # plot title
    plt_title <- sprintf("%s for %s, T=%sh",
                         gDRutils::prettify_flat_metrics(x = mx_name, human_readable = TRUE),
                         normalization_type,
                         unique(dt_excess[get(cellline_name) == cl_name][[duration]]))
    if (!as_panel) plt_title <- paste(main_title, plt_title, sep = " : ")
    
    legend_title_fill <- sprintf("%s %s",
                                 gDRutils::prettify_flat_metrics(x = mx_name, human_readable = TRUE),
                                 normalization_type)
    
    # prep limits
    limits <- prep_hm_limits(dt_[[mx_name]])
    
    # base plot
    plt <- 
      ggplot2::ggplot(dt_, ggplot2::aes(x = pos_x, y = pos_y)) +
      ggplot2::geom_tile(ggplot2::aes(fill = get(mx_name)), height = tile_height, width = tile_width) +
      ggplot2::labs(x = bquote(.(drug2_name) ~ "[" ~ mu * M ~ "]"),
                    y = bquote(.(drug1_name) ~ "[" ~ mu * M ~ "]"),
                    title = plt_title,
                    fill = legend_title_fill) +
      ggplot2::theme_bw() + 
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = 9, angle = 45, vjust = 1, hjust = 1),
                     axis.text.y = ggplot2::element_text(size = 9),
                     plot.title = ggplot2::element_text(size = 11),
                     panel.grid.minor = ggplot2::element_blank()) +
      ggplot2::scale_x_continuous(breaks = drug2_axis$pos_x, 
                                  labels = drug2_axis$marks_x,
                                  expand = c(0, 0)) +
      ggplot2::scale_y_continuous(breaks = drug1_axis$pos_y, 
                                  labels = drug1_axis$marks_y,
                                  expand = c(0, 0)) + 
      ggplot2::scale_fill_gradientn(colors = hm_color_palette,
                                    limit = limits)
    
    # add isoline
    if (NROW(avialable_iso)) { # isobolograms as lines
      if (all(avialable_iso %in% c("0.25", "0.5", "0.75"))) {
        # friendly for user with color vision deficiency
        plt <- plt +
          ggplot2::geom_path(data = dt_isobolograms, linewidth = 1,
                             ggplot2::aes(x = pos_x, y = pos_y, color = iso_level, linetype = iso_level)) +
          ggplot2::scale_color_manual(values = iso_colors[avialable_iso],
                                      breaks = avialable_iso,
                                      labels = legend_lbl_iso,
                                      name = legend_title_iso) +
          ggplot2::scale_linetype_manual(values = c("solid", "twodash", "dashed"),
                                         breaks = avialable_iso,
                                         labels = legend_lbl_iso,
                                         name = legend_title_iso) +
          ggplot2::theme(legend.key.width = ggplot2::unit(3, "line"))
        
      } else {
        plt <- plt +
          ggplot2::geom_path(data = dt_isobolograms, linewidth = 1,
                             ggplot2::aes(x = pos_x, y = pos_y, color = iso_level)) +
          ggplot2::scale_color_manual(values = iso_colors[avialable_iso],
                                      breaks = avialable_iso,
                                      labels = legend_lbl_iso,
                                      name = legend_title_iso)
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
  legend_title <- ifelse(normalization_type == "GR", "GR", "IC")
  
  if (!as_panel) plt_title <- paste(main_title, plt_title, sep = " : ")
  # base plot
  plt_iso_compare <- 
    ggplot2::ggplot(mapping = ggplot2::aes(x = log10_ratio_conc, y = log2_CI)) +
    ggplot2::geom_line(
      data = data.table::data.table(log10_ratio_conc = c(-2, 2), log2_CI = c(0, 0))) +
    ggplot2::geom_hline(yintercept = 0, color = "#A9A9A9")
  
  # add isoline
  if (all(avialable_iso %in% c("0.25", "0.5", "0.75"))) {
    # friendly for user with color vision deficiency
    plt_iso_compare <- plt_iso_compare +
      ggplot2::geom_path(data = dt_isobolograms, linewidth = 0.5,
                         ggplot2::aes(x = log10_ratio_conc, y = log2_CI, color = iso_level, linetype = iso_level)) +
      ggplot2::scale_color_manual(values = iso_colors[avialable_iso],
                                  breaks = avialable_iso,
                                  labels = legend_lbl_iso,
                                  name = legend_title_iso) +
      ggplot2::scale_linetype_manual(values = c("solid", "twodash", "dashed"),
                                     breaks = avialable_iso,
                                     labels = legend_lbl_iso,
                                     name = legend_title_iso) +
      ggplot2::theme(legend.key.width = ggplot2::unit(3, "line"))
    
  } else {
    plt_iso_compare <- plt_iso_compare +
      ggplot2::geom_path(data = dt_isobolograms, linewidth = 0.5,
                         ggplot2::aes(x = log10_ratio_conc, y = log2_CI, color = iso_level))
    ggplot2::scale_color_manual(values = iso_colors[avialable_iso],
                                breaks = avialable_iso,
                                labels = legend_lbl_iso,
                                name = legend_title_iso)
  }
  
  # add x and y scales
  plt_iso_compare <- plt_iso_compare +
    ggplot2::scale_y_continuous(breaks = -5:4, labels = c(paste0("1/", 2 ^ (5:1)), 2 ^ (0:4))) +
    ggplot2::scale_x_continuous(breaks = -3:3, labels = c(paste0("1/", 10 ^ (3:1)), 10 ^ (0:3))) +
    ggplot2::coord_cartesian(ylim = c(-5, 4)) + 
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 9, angle = 45, vjust = 1, hjust = 1),
                   axis.text.y = ggplot2::element_text(size = 9),
                   plot.title = ggplot2::element_text(size = 11),
                   panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::labs(y = "CI",
                  x = paste(drug2_name, "/", drug1_name, "ratio"),
                  title = plt_title)
  
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
            ls_plts[["hsa_excess"]] + ggplot2::guides(fill = ggplot2::guide_legend(title = "Excess")),
            ls_plts[["bliss_excess"]] + ggplot2::guides(fill = ggplot2::guide_legend(title = "Excess"))
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

#' Plot averaged values heatmaps for combo data
#'
#' @param dt_excess data.table representation of the data in \code{excess} assay
#'    output from \code{gDRutils::convert_se_assay_to_dt(se, "excess")}
#' @param dt_isobolograms data.table representation of the data in \code{isobolograms} assay
#'    output from \code{gDRutils::convert_se_assay_to_dt(se, "isobolograms")}
#' @param drug1_name string with drug name to be plotted (identifiers \code{DrugName})
#' @param drug2_name string with co-drug name to be plotted (identifiers \code{DrugName_2})
#' @param cl_name string with cell line to be plotted (identifiers \code{CellLineName}) 
#' @param normalization_type string with normalization_types to be selected
#'                           one of: "GR" ("GRvalue") or "RV" ("RelativeViability")
#' @param iso_levels character vector with  isobologram levels to be selected
#' @param colors_vec character vector of colors (valid name or hex) used in heatmap
#' @param no_breaks numeric number of breaks on scale
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
  
  checkmate::expect_data_table(dt_excess)
  checkmate::expect_data_table(dt_isobolograms)
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
    stopifnot("Must be valid color name" = all(vapply(colors_vec, gDRplots::is_valid_color, logical(1))))
  }
  checkmate::assert_int(no_breaks, lower = 2)
  
  # filter data for normalization type
  data.table::setkeyv(dt_excess, "normalization_type")
  dt_excess <- dt_excess[normalization_type]
  data.table::setkey(dt_excess, NULL)
  
  data.table::setkeyv(dt_isobolograms, "normalization_type")
  dt_isobolograms <- dt_isobolograms[normalization_type]
  data.table::setkey(dt_isobolograms, NULL)
  
  # filter data for combination cell line (drug x drug2)
  dt_excess <- 
    dt_excess[get(cellline_name) == cl_name & get(drug_name) == drug1_name & get(drug_name_2) == drug2_name]
  dt_isobolograms <- 
    dt_isobolograms[get(cellline_name) == cl_name & get(drug_name) == drug1_name & get(drug_name_2) == drug2_name]
  
  dt_isobolograms <- dt_isobolograms[iso_level %in% iso_levels, ]
  available_iso_lvl <- unique(dt_isobolograms[["iso_level"]])
  
  # prep hm color palette
  hm_color_palette <- if (is.null(colors_vec)) {
    colorspace::sequential_hcl(no_breaks + 1, palette = "viridis")
  } else {
    grDevices::colorRampPalette(colors_vec)(no_breaks + 1)
  }
  
  # prep plot data
  mx_name <- "smooth" 
  dt_ <- dt_excess[, c(conc, conc_2, mx_name), with = FALSE]
  # correction of NA for conc = 0 ir conc_2 = 0
  dt_[(get(conc) == 0 | get(conc_2) == 0) & is.na(get(mx_name))] <- 0
  
  dt_[[mx_name]] <- pmin(1.1, dt_[[mx_name]])
  dt_$pos_y <- log10(dt_[[conc]])
  dt_$pos_x <- log10(dt_[[conc_2]])
  
  ls_axes <- gDRutils::define_matrix_grid_positions(dt_[[conc]], dt_[[conc_2]])
  drug1_axis <- ls_axes$axis_1
  drug2_axis <- ls_axes$axis_2
  tile_height <- diff(drug1_axis$pos_y[3:4])
  tile_width <- diff(drug2_axis$pos_x[3:4])
  
  # panel title
  cl_clid <- unique(dt_excess[get(cellline_name) == cl_name, ][[clid]]) 
  plt_title <- sprintf("%s (%s)", cl_name, cl_clid)
  
  legend_title_fill <- sprintf("%s %s",
                               gDRutils::prettify_flat_metrics(x = mx_name, human_readable = TRUE),
                               normalization_type)

  # prep limits
  limits <- prep_hm_limits(dt_[[mx_name]])
  
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
                                  limit = limits)
  
  # plot isobologram
  if (NROW(dt_isobolograms)) { # add isolines - if there are such data
    iso_label <- sprintf("%s%s", 
                         ifelse(normalization_type == "GR", "GR", "IC"), 
                         100 - 100 * as.numeric(available_iso_lvl))
    names(iso_label) <- available_iso_lvl
    
    iso_source <- NULL # due to NSE notes in R CMD check
    tab_measured <- dt_isobolograms[, .SD, .SDcols = -c("pos_x_ref", "pos_y_ref")][, iso_source := "measured"]
    tab_expected <- dt_isobolograms[, .SD, .SDcols = -c("pos_x", "pos_y")][, iso_source := "expected"]
    data.table::setnames(tab_expected, old = c("pos_x_ref", "pos_y_ref"), new = c("pos_x", "pos_y"))
    
    tab_isoline <- rbind(tab_measured, tab_expected)
    
    if (NROW(available_iso_lvl) == 1) {
      plt <- plt +
        ggplot2::geom_path(data = tab_isoline, 
                           ggplot2::aes(x = pos_x, y = pos_y, linetype = iso_source), 
                           linewidth = 1, color = "red") +
        ggplot2::scale_linetype_manual(values = c("measured" = "solid", "expected" = "dashed"),
                                       name = iso_label) 
    } else {
      iso_colors <- 
        grDevices::colorRampPalette(c("red", "darkred"))(2 * NROW(available_iso_lvl))[seq_along(available_iso_lvl) * 2]
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
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 9, angle = 45, vjust = 1, hjust = 1),
                   axis.text.y = ggplot2::element_text(size = 9),
                   plot.title = ggplot2::element_text(size = 11),
                   panel.grid.minor = ggplot2::element_blank(),
                   legend.key.width = ggplot2::unit(2, "line"))
  
  
  return(plt)
}


#' Plot panel with thefitted values heatmaps and references data for isolobograms 
#' to control quality of the data
#'
#' @inheritParams heatmap_combo_with_isoref
#' @param cl_names character vector with cell line names names to be plotted (Cell Line NAme); 
#'    if NULL - all available cell liene will be plotted
#'    
#' @return panel with heatmaps for fitted values and references data for isolobograms 
#'    for selected drug and co-drug by cell line names
#'
#' @keywords QC_plot
#' @examples
#' drug1_name <- "drug_001"
#' drug2_name <- "drug_026"
#' 
#' mae <- gDRutils::get_synthetic_data("combo_matrix")
#' se <- mae[[gDRutils::get_supported_experiments("combo")]]
#' dt_excess <- gDRutils::convert_se_assay_to_dt(se, "excess")
#' dt_isobolograms <- gDRutils::convert_se_assay_to_dt(se, "isobolograms")
#' 
#' cl_names <- unique(dt_excess[["CellLineName"]])[1:4]
#' 
#' heatmap_combo_with_isoref_qc_panel(dt_excess,
#'                                    dt_isobolograms,
#'                                    drug1_name, drug2_name,
#'                                    cl_names,
#'                                    normalization_type = "GR")
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
    colors_vec = NULL) {
  
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  clid <- gDRutils::get_env_identifiers("cellline")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  gnumber <- gDRutils::get_env_identifiers("drug")
  drug_name_2 <- gDRutils::get_env_identifiers("drug_name2")
  gnumber_2 <- gDRutils::get_env_identifiers("drug2")
  
  checkmate::expect_data_table(dt_excess)
  checkmate::expect_data_table(dt_isobolograms)
  checkmate::assert_string(drug1_name)
  checkmate::assert_choice(drug1_name, choices = dt_excess[[drug_name]])
  checkmate::assert_string(drug2_name)
  checkmate::assert_choice(drug2_name, choices = dt_excess[[drug_name_2]])
  checkmate::assert_character(cl_names, null.ok = TRUE)
  checkmate::assert_choice(normalization_type, choices = c("GR", "RV"))
  checkmate::assert_character(iso_levels, null.ok = TRUE)
  if (!is.null(iso_levels)) checkmate::assert_numeric(as.numeric(iso_levels))
  checkmate::assert_character(colors_vec, null.ok = TRUE)
  if (!is.null(colors_vec)) {
    stopifnot("Must be valid color name" = all(vapply(colors_vec, gDRplots::is_valid_color, logical(1))))
  }
  
  available_cls <- unique(dt_excess[[cellline_name]])
  if (is.null(cl_names) || all(!cl_names %in% available_cls)) {
    cl_names  <- available_cls
  } else if (!all(cl_names %in% available_cls)) {
    cl_names <- cl_names[cl_names %in% available_cls]
  } 
  
  ls_celllines <- list(cl_name = cl_names)
  
  # panel title
  panel_title <- sprintf("%s (%s) x %s (%s)",
                         drug1_name,
                         unique(dt_excess[get(drug_name) == drug1_name, ][[gnumber]]),
                         drug2_name,
                         unique(dt_excess[get(drug_name_2) == drug2_name, ][[gnumber_2]]))
  
  # list of plots for each drug
  ls_plt <- purrr::pmap(ls_celllines, 
                        gDRplots::heatmap_combo_with_isoref,
                        dt_excess = dt_excess,
                        dt_isobolograms = dt_isobolograms,
                        drug1_name = drug1_name,
                        drug2_name = drug2_name,
                        normalization_type = normalization_type,
                        iso_levels = iso_levels,
                        colors_vec = colors_vec)
  
  names(ls_plt) <- cl_names
  
  # final panel
  panel <- ggpubr::annotate_figure(
    ggpubr::ggarrange(plotlist = ls_plt, common.legend = TRUE, legend = "right"),
    top = panel_title) + 
    ggpubr::bgcolor("white") + ggpubr::border("white")
  
  return(panel)
}




#' Calculate limit 
#'
#' @param num_vec numeric vector
#' @param lower_cap numeric lower capping value
#' @param upper_cap numeric upper capping value
#'
#' @return limit for given numeric vector
#'
#' @keywords internal
#' @examples
#' \dontrun{
#' vec <- c(-5, -0.3, 0, Inf, NA)
#' prep_hm_limits(vec)
#' }
#' 
prep_hm_limits <- function(num_vec,
                           lower_cap = -0.25,
                           upper_cap = 0.25
) {
  checkmate::assert_numeric(num_vec)
  checkmate::assert_number(lower_cap)
  checkmate::assert_number(upper_cap)
  
  vec_range <- range(num_vec, na.rm = TRUE, finite = TRUE)
  
  min_ <- min(c(lower_cap, min(vec_range)))
  max_ <- max(c(upper_cap, max(vec_range))) 

  return(c(min_, max_))
}
