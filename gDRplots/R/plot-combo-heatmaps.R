#' Plot averaged values heatmaps for combo data
#'
#' @param se combination \code{SummarizedExperiment} object holding raw and/or processed dose-response 
#'    data in its assays for one cell line
#' @param drug1_name string with drug name to be plotted (identifiers \code{DrugName})
#' @param drug2_name string with co-drug name to be plotted (identifiers \code{DrugName_2})
#' @param cl_name string with cell line to be plotted (identifiers \code{CellLineName}) 
#' @param normalization_type string with normalization_types to be selected
#'                           one of: "GR" ("GRvalue") or "RV" ("RelativeViability")
#' @param iso_levels character vectore with  isobologram levels to be selected
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
#' se <- mae[["combination"]]
#' 
#' heatmap_combo_metrics(se, 
#'                       drug1_name, drug2_name, 
#'                       cl_name, 
#'                       normalization_type = "GR")
#'                       
#' heatmap_combo_metrics(se, 
#'                       drug1_name, drug2_name, 
#'                       cl_name, 
#'                       normalization_type = "RV", 
#'                       as_panel = FALSE)
#'
#' @export
heatmap_combo_metrics <- function(se,
                                  drug1_name,
                                  drug2_name,
                                  cl_name,
                                  normalization_type = "GR",
                                  iso_levels =  c("0.25", "0.5", "0.75"),
                                  as_panel = TRUE) {
  
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  drug_name2 <- gDRutils::get_env_identifiers("drug_name2")
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  
  checkmate::assert_class(se, "SummarizedExperiment")
  checkmate::assert_string(drug1_name)
  checkmate::assert_choice(drug1_name, choices = SummarizedExperiment::rowData(se)[[drug_name]])
  checkmate::assert_string(drug2_name)
  checkmate::assert_choice(drug2_name, choices = SummarizedExperiment::rowData(se)[[drug_name2]])
  checkmate::assert_string(cl_name)
  checkmate::assert_choice(cl_name, choices = SummarizedExperiment::colData(se)[[cellline_name]])
  checkmate::assert_character(iso_levels)
  checkmate::assert_numeric(as.numeric(iso_levels))
  checkmate::assert_choice(normalization_type, choices = c("GR", "RV"))
  
  selected_col <- 
    SummarizedExperiment::colData(se)[SummarizedExperiment::colData(se)$CellLineName == cl_name, ]
  selected_row <- 
    SummarizedExperiment::rowData(se)[SummarizedExperiment::rowData(se)$DrugName == drug1_name & 
                                        SummarizedExperiment::rowData(se)$DrugName_2 == drug2_name, ]
  sel_se <- se[rownames(selected_row), rownames(selected_col)]
  
  stopifnot("`se` has to by filterd to one drug and one cell line" = 
              (NCOL(sel_se) == 1 && NROW(sel_se) == 1)) 
  
  dt_excess <- data.table::as.data.table(
    BumpyMatrix::unsplitAsDataFrame(SummarizedExperiment::assay(sel_se, "excess")))
  data.table::setkeyv(dt_excess, "normalization_type")
  dt_excess <- dt_excess[normalization_type]
  data.table::setkey(dt_excess, NULL)
  
  dt_isobolograms <- data.table::as.data.table(
    BumpyMatrix::unsplitAsDataFrame(SummarizedExperiment::assay(sel_se, "isobolograms")))
  data.table::setkeyv(dt_isobolograms, "normalization_type")
  dt_isobolograms <- dt_isobolograms[normalization_type]
  data.table::setkey(dt_isobolograms, NULL)

  iso_colors <- gDRutils::get_iso_colors()[iso_levels]
  dt_isobolograms <- 
    dt_isobolograms[iso_level %in% iso_levels, ]
  avialable_iso <- unique(dt_isobolograms$iso_level)
  
  
  # variables
  conc <- gDRutils::get_env_identifiers("concentration")
  conc_2 <- gDRutils::get_env_identifiers("concentration2")
  duration <- gDRutils::get_env_identifiers("duration")
  clid <- gDRutils::get_env_identifiers("cellline")
  mx_names <- names(gDRutils::get_combo_excess_field_names())
  
  main_title <- sprintf("%s (%s)",
                        cl_name,
                        selected_col[[clid]])
  # plots
  mx_plts <- lapply(mx_names, function(mx_name) {
    dt_ <- dt_excess[, c(conc, conc_2, mx_name), with = FALSE]
    dt_[[mx_name]] <- pmin(1.1, dt_[[mx_name]])
    dt_$pos_y <- log10(dt_[[conc]])
    dt_$pos_x <- log10(dt_[[conc_2]])
    
    ls_axes <- gDRutils::define_matrix_grid_positions(dt_[[conc]], dt_[[conc_2]])
    drug1_axis <- ls_axes$axis_1
    drug2_axis <- ls_axes$axis_2
    tile_height <- diff(drug1_axis$pos_y[3:4])
    tile_width <- diff(drug2_axis$pos_x[3:4])
    
    plt_title <- sprintf("%s for %s, T=%sh",
                         gDRutils::get_combo_excess_field_names()[[mx_name]],
                         normalization_type,
                         selected_row[[duration]])
    if (!as_panel) plt_title <- paste(main_title, plt_title, sep = " : ")
    
    # base plot
    plt <- 
      ggplot2::ggplot(dt_, ggplot2::aes(x = pos_x, y = pos_y)) +
      ggplot2::geom_tile(ggplot2::aes(fill = get(mx_name)), height = tile_height, width = tile_width) +
      ggplot2::labs(x = bquote(.(drug2_name) ~ "[" ~ mu * M ~ "]"),
                    y = bquote(.(drug1_name) ~ "[" ~ mu * M ~ "]"),
                    title = plt_title,
                    fill = gDRutils::get_combo_excess_field_names()[[mx_name]]) +
      ggplot2::theme_bw() + 
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = 9, angle = 45, vjust = 1, hjust = 1),
                     axis.text.y = ggplot2::element_text(size = 9),
                     plot.title = ggplot2::element_text(size = 11)) +
      ggplot2::scale_x_continuous(breaks = drug2_axis$pos_x, labels = drug2_axis$marks_x,
                                  expand = c(0, 0)) +
      ggplot2::scale_y_continuous(breaks = drug1_axis$pos_y, labels = drug1_axis$marks_y,
                                  expand = c(0, 0)) +
      ggplot2::scale_shape_discrete(name = paste0(ifelse(normalization_type == "GR", "GR", "IC"), "50"))
    
    # add color scale
    if (normalization_type == "GR") {
      plt <- plt +
        ggplot2::scale_fill_gradientn(
          colors = c("black", "#b06000", "#c07700", "white"),
          values = c(0, 0.59 / 1.7, 0.61 / 1.7, 1), 
          limits = c(-0.6, 1.1),
          name = "GR",
          oob = scales::squish)
    } else {
      plt <- plt +
        ggplot2::scale_fill_gradientn(
          colors = c("#440000", "#ff5500", "white"), values = c(0, 0.4, 1),
          limits = c(0, 1.1), 
          name = "RV")
    }

    if (length(avialable_iso) > 0) { # three isobolograms as lines avialable
      plt <- plt +
        ggplot2::geom_path(data = dt_isobolograms, linewidth = 0.5,
                           ggplot2::aes(x = pos_x, y = pos_y, color = iso_level)) +
        ggplot2::scale_color_manual(values = iso_colors[avialable_iso],
                                    breaks = names(iso_colors[avialable_iso]),
                                    labels = paste0(ifelse(normalization_type == "GR", "GR", "IC"),
                                                    100 - 100 * as.numeric(avialable_iso)),
                                    name = "Iso Levels")
    }
    return(plt)
  })
  names(mx_plts) <- mx_names
  
  # isobolograms across range of concentration ratios
  plt_title <- sprintf("%s for %s, T=%sh", 
                       gDRutils::get_combo_excess_field_names()[["smooth"]],
                       normalization_type,
                       selected_row[[duration]]) 
  if (!as_panel) plt_title <- paste(main_title, plt_title, sep = " : ")
  # base plot
  plt_iso_compare <- 
    ggplot2::ggplot(mapping = ggplot2::aes(x = log10_ratio_conc, y = log2_CI)) +
    ggplot2::geom_line(
      data = data.table::data.table(log10_ratio_conc = c(-2, 2), log2_CI = c(0, 0))) +
    ggplot2::geom_hline(yintercept = 0)
  
  # add isoline
  plt_iso_compare <- plt_iso_compare +
    ggplot2::geom_path(data = dt_isobolograms, linewidth = 0.5,
                       ggplot2::aes(x = log10_ratio_conc, y = log2_CI, color = iso_level)) + 
    ggplot2::scale_color_manual(values = iso_colors[avialable_iso],
                                name = ifelse(normalization_type == "GR", "GR", "IC")) +
    ggplot2::labs(color = "Iso levels") +
    ggplot2::scale_y_continuous(breaks = -5:4, labels = c(paste0("1/", 2 ^ (5:1)), 2 ^ (0:4))) +
    ggplot2::scale_x_continuous(breaks = -3:3, labels = c(paste0("1/", 10 ^ (3:1)), 10 ^ (0:3))) +
    ggplot2::coord_cartesian(ylim = c(-5, 4)) + 
    ggplot2::theme_bw() +
    ggplot2::labs(y = "CI",
                  x = paste(drug2_name, "/", drug1_name, "ratio"),
                  title = plt_title,
                  color = ifelse(normalization_type == "GR", "GR", "IC"))
  
  # final plots
  ls_plts <- append(mx_plts, list(iso_compare = plt_iso_compare))
  
  final_plot <- if (as_panel) {
    ggpubr::annotate_figure(
      ggpubr::ggarrange(plotlist = ls_plts, nrow = 2, ncol = 2, 
                        common.legend = TRUE, legend = "right"),
      top = main_title)
  } else {
    ls_plts
  }
  return(final_plot)
}
