#' Plot averaged values heatmaps for combo data
#'
#' @param se combination \code{SummarizedExperiment} object holding raw and/or processed dose-response 
#'    data in its assays for one cell line
#' @param drug1_name string with cell line to be plotted (identifiers \code{DrugName})
#' @param drug2_name string with cell line to be plotted (identifiers \code{DrugName_2})
#' @param cl_name string with cell line to be plotted (identifiers \code{CellLineName}) 
#' @param normalization_type string with normalization_types to be selected
#'                           one of: "GR" ("GRvalue") or "RV" ("RelativeViability")
#'
#' @return list of heatmaps with value for excess assays for selected drugs and cell line with
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
#' gDR_combo_plot(se, drug1_name, drug2_name, cl_name, normalization_type = "GR")
#'
#' @export
gDR_combo_plot <- function(se,
                           drug1_name,
                           drug2_name,
                           cl_name,
                           normalization_type = "GR") {
  
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
  all_iso <- unique(dt_isobolograms$iso_level)
  iso_colors <- gDRutils::get_iso_colors()[all_iso]
  
  # variables
  conc <- gDRutils::get_env_identifiers("concentration")
  conc_2 <- gDRutils::get_env_identifiers("concentration2")
  duration <- gDRutils::get_env_identifiers("duration")
  clid <- gDRutils::get_env_identifiers("cellline")
  mx_names <- names(gDRutils::get_combo_excess_field_names())
  
  # plots
  mx_plts <- lapply(mx_names, function(mx_name) {
    dt_ <- dt_excess[, c(conc, conc_2, mx_name), with = FALSE]
    dt_[[mx_name]] <- pmin(1.1, dt_[[mx_name]])
    dt_$pos_y <- log10(dt_[[conc]])
    dt_$pos_x <- log10(dt_[[conc_2]])
    
    ls_axes <- gDRcore:::define_matrix_grid_positions(dt_[[conc]], dt_[[conc_2]])
    drug1_axis <- ls_axes$axis_1
    drug2_axis <- ls_axes$axis_2
    tile_height <- diff(drug1_axis$pos_y[3:4])
    tile_width <- diff(drug2_axis$pos_x[3:4])
    
    plt_title <- sprintf("%s (%s) : %s for %s, T=%sh",
                         cl_name,
                         selected_col[[clid]], 
                         gDRutils::get_combo_excess_field_names()[[mx_name]],
                         normalization_type,
                         selected_row[[duration]])
    
    # base plot
    plt <- 
      ggplot2::ggplot(dt_, ggplot2::aes(x = pos_x, y = pos_y)) +
      ggplot2::geom_tile(ggplot2::aes(fill = get(mx_name)), height = tile_height, width = tile_width) +
      # TODO
      # ggplot2::geom_abline(slope = 1, intercept = log10(ref_x50["conc_1"]/ref_x50["conc_2"])) + 
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
    if (!(mx_name %in% c("hsa_excess", "bliss_excess"))) { # heatmaps with readout values
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
    } else { # bliss/hsa excess matrix
      plt <- plt +
        ggplot2::scale_fill_gradientn(
          colors = c("black", "#ffffaa", "white", "white", "#aaffff", "blue"),
          values = c(0, 0.35, 0.48, 0.51, 0.65, 1), limits = c(-0.6, 0.6),
          name = gDRutils::get_combo_excess_field_names()[[mx_name]],
          oob = scales::squish)
    }
    
    # TODO add selected iso level
    # if ("0.5" %in% all_iso) { # points of the isobologram at GR/IC50
    #   plt <- plt + ggplot2::geom_point(data = data.table::as.data.table(
    #     all_iso[[max(1, which(all_iso == "0.5"))]]$dt_iso),
    #     ggplot2::aes(shape = fit_type), show.legend = FALSE)
    # }
    
    if (length(all_iso) > 0) { # three isobolograms as lines
      select_iso <- all_iso[
        sort(unique(c(max(1, which(all_iso == "0.5")),
                      which(all_iso %in% ifelse(normalization_type[c(1, 1, 1)] == "GR",
                                                c(0.5, 0.25, 0), c(0.75, 0.5, 0.25)))
        )), TRUE)]
      
      plt <- plt +
        ggplot2::geom_path(data = dt_isobolograms[iso_level %in% select_iso, ], 
                           linewidth = 0.5,
                           ggplot2::aes(x = pos_x, y = pos_y, color = iso_level)) +
        ggplot2::scale_color_manual(values = iso_colors[select_iso],
                                    breaks = names(iso_colors[select_iso]),
                                    labels = paste0(ifelse(normalization_type == "GR", "GR", "IC"),
                                                    100 - 100 * as.numeric(select_iso)),
                                    name = "Iso Levels")
    }
    return(plt)
  })
  names(mx_plts) <- mx_names
  
  # isobolograms across range of concentration ratios
  plt_title <- sprintf("%s (%s) : %s for %s, T=%sh",
                       cl_name,
                       selected_col[[clid]], 
                       gDRutils::get_combo_excess_field_names()[["smooth"]],
                       normalization_type,
                       selected_row[[duration]]) 
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
    ggplot2::scale_color_manual(values = iso_colors[all_iso],
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
  
  return(
    append(mx_plts, list(iso_compare = plt_iso_compare))
  )
}
