#' Volcano plot with association
#'
#' @param dt_assoc \code{data.table} with the calculated linear association between DepMap and metrics
#'     outputted by \code{kaleidoscope::calc_assoc}
#' @param feature_info string describing name of associated feature/metadata from DepMap
#' @param condition_info string describing experiment condition 
#'     (preferred: \code{"DrugName"}_\code{"Gnumber"}_\code{"drug_moa"}_\code{"Duration"})
#' @param q_cutoff numeric cutoff to identify statistically significant correlations
#' @param named_p_top numeric value for p-top statistically significant correlations to be labeled on the plot
#' @param max_N numeric value for limit maximum number of non-statistically significant points to plot; 
#'     for default \code{NULL} all points will be plotted.
#'
#' @return a volcano plot with association
#' @keywords prism_plots
#' 
#' @export
plot_volcano_assoc <- function(dt_assoc,
                               feature_info,
                               condition_info = NULL,
                               q_cutoff = 0.05,
                               named_p_top = 10,
                               max_N = NULL) {
  
  checkmate::assert_data_table(dt_assoc)
  checkmate::assert_string(feature_info)
  checkmate::assert_string(condition_info, null.ok = TRUE)
  checkmate::assert_number(q_cutoff, lower = 0, upper = 1)
  checkmate::assert_number(named_p_top, lower = 0)
  checkmate::assert_number(max_N, lower = 10, null.ok = TRUE)
  
  x_lbl <- "rho"
  y_lbl <- "neglog_q_value"
  stat_sig <- q_value <- label <- feature <- NULL # due to NSE notes in R CMD check
  
  checkmate::assert_names(names(dt_assoc), must.include = c(x_lbl, "q_value", "feature"))
  
  plt_title <- sprintf("%s__%s", unique(dt_assoc[["response"]]), feature_info)
  
  if (all(is.na(dt_assoc[["q_value"]]))) {
    # empty plot
    plt <- 
      ggplot2::ggplot() + 
      ggplot2::labs(title = paste(plt_title, ": all NAs")) +
      ggplot2::theme_bw()
  } else {
    tab_plot <- data.table::setorderv(data.table::copy(dt_assoc), cols = "q_value")
    
    # prep column with statistically significant
    tab_plot[, stat_sig :=  data.table::fifelse(q_value <= q_cutoff, "yes", "no")]
    
    # downsample non-statistically significant dots
    if (!is.null(max_N) && NROW(tab_plot[stat_sig == "no", ]) > max_N) {
      idx_no <- which(tab_plot$stat_sig == "no")
      tab_plot <- rbind(tab_plot[stat_sig == "yes", ], 
                        tab_plot[sample(x = idx_no, size = max_N, replace = FALSE), ])
    }
    
    # add labels to named_p_top correlations
    data.table::setorderv(tab_plot, cols = "q_value")
    named_p_top_act <- min(named_p_top, NROW(tab_plot)) # deal with less than p-entries
    tab_plot[["label"]] <- ""
    tab_plot[1:named_p_top_act, label := feature]
    
    # prep data
    tab_plot[[y_lbl]] <- -log10(tab_plot$q_value)
    
    # volcano plot
    plt <- 
      ggplot2::ggplot(
        data = tab_plot,
        mapping = ggplot2::aes(x = get(x_lbl), y = get(y_lbl), label = label, color = stat_sig)) +
      ggplot2::geom_point() +
      ggplot2::scale_x_continuous(trans = "identity", name = x_lbl) +
      ggplot2::scale_y_continuous(trans = "identity", name = y_lbl) +
      ggplot2::scale_color_manual(values = list(yes = "black", no = "#A9A9A9"),
                                  name = "Statistically Significant") +
      ggrepel::geom_text_repel(size = 4, show.legend = FALSE) +
      ggplot2::labs(title = plt_title, subtitle = condition_info) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "bottom", 
                     legend.title = ggplot2::element_text(vjust = 0.5, hjust = 1))
  }
  
  return(plt)
}

#' Plot scatter with correlation
#'
#' @param dt_response \code{data.table} with experimental response data (rows are samples) for one metric
#'  outputted by one of functions: \code{\link[gDRplots]{prep_dt_response_metric_sa}},
#'  \code{\link[gDRplots]{prep_dt_response_dose_sa}}, \code{\link[gDRplots]{prep_dt_response_scores}}
#'  or \code{\link[gDRplots]{prep_dt_response_metric_diff}}, 
#' @param dt_depmap \code{data.table} with dependent variables data loaded from DepMap - for one
#'    feature or one metadata; (rows are samples, columns are features or meta). 
#' @param selected_feat string with name of selected feature from \code{dt_depmap}
#' @param selected_feat_meta_col string name of feature column in DepMap
#'
#' @return a scatter plot with correlation
#' @keywords prism_plots
#' 
#' @export
plot_scatter_with_corr <- function(dt_response,
                                   dt_depmap, 
                                   selected_feat,
                                   selected_feat_meta_col = NULL) {
  
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  
  checkmate::assert_data_table(dt_response)
  checkmate::assert_data_table(dt_depmap)
  checkmate::assert_names(names(dt_depmap), must.include = c("CCLEName", selected_feat))
  checkmate::assert_string(selected_feat)
  checkmate::assert_string(selected_feat_meta_col, null.ok = TRUE)
  
  selected_metric <- setdiff(names(dt_response), 
                             c(cellline_name, "rId", "cId"))
  stopifnot("Provide `dt_response` for one metric." = NROW(selected_metric) == 1)
  
  CCLEName <- NULL # due to NSE notes in R CMD check
  
  # prep table with data to plot
  X_dt <- dt_depmap[, c("CCLEName", selected_feat), with = FALSE]
  Y_dt <- dt_response[, c(cellline_name, selected_metric), with = FALSE]
  tab_plot <- Y_dt[X_dt, on = .(CellLineName = CCLEName), nomatch = NULL]
  # remove NA
  tab_plot <- stats::na.omit(tab_plot)
  
  # re-calculate correlation
  c <- stats::cor(tab_plot[[selected_feat]], tab_plot[[selected_metric]], 
                  method = "pearson", use = "pairwise.complete.obs") 
  # calculate slope (b1 = r * (sd(y) / sd(x)))
  slope <- c * (stats::sd(tab_plot[[selected_metric]]) / stats::sd(tab_plot[[selected_feat]])) 
  # calculate intercept (b0 = mean(y) - b1 * mean(x))
  intercept <- mean(tab_plot[[selected_metric]]) - slope * mean(tab_plot[[selected_feat]])
  
  # plot title
  plt_subtitle <- 
    sprintf("corr=%2.2f, slope=%2.2f, intercept=%2.2f", c, slope, intercept)
  
  plt <-        
    ggplot2::ggplot(
      data = tab_plot,
      mapping =  ggplot2::aes(x = get(selected_feat), y = get(selected_metric), label = get(cellline_name))) +
    ggplot2::geom_point() +
    ggplot2::scale_x_continuous(trans = "identity", name = selected_feat) +
    ggplot2::scale_y_continuous(trans = "identity", name = selected_metric) +
    # ggrepel::geom_text_repel(size = 2) + # nolint
    ggplot2::geom_abline(intercept = intercept, slope = slope, color = "red") +   
    ggplot2::labs(title = selected_feat_meta_col, subtitle = plt_subtitle) +
    ggplot2::theme_bw()
  
  return(plt)
}


#' Plot boxplot for metric values grouped by metadata from DepMap
#'
#' @param dt_response \code{data.table} with experimental response data (rows are samples) for one metric
#'  outputted by one of functions: \code{\link[gDRplots]{prep_dt_response_metric_sa}},
#'  \code{\link[gDRplots]{prep_dt_response_dose_sa}}, \code{\link[gDRplots]{prep_dt_response_scores}}
#'  or \code{\link[gDRplots]{prep_dt_response_metric_diff}}, 
#' @param dt_depmap_lng \code{data.table} with dependent variables data loaded from DepMap - for one
#'    metadata in long format - column with cellline names and second with values for meta.
#' @param selected_meta string with name of selected meta data from \code{dt_depmap}
#' @param with_1_item_grp logical flag whether to show group with only one item
#' @param max_x_lbl_length numeric value for max character number of x-axis label
#'
#' @return a boxplot
#' 
#' @keywords prism_plots
#' 
#' @export
plot_boxplot_meta <- function(dt_response,
                              dt_depmap_lng, 
                              selected_meta,
                              with_1_item_grp = TRUE,
                              max_x_lbl_length = 60) {
  
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  
  checkmate::assert_data_table(dt_response)
  checkmate::assert_data_table(dt_depmap_lng)
  checkmate::assert_string(selected_meta)
  checkmate::assert_names(names(dt_depmap_lng), must.include = c("CCLEName", selected_meta))
  checkmate::assert_flag(with_1_item_grp)
  checkmate::assert_int(max_x_lbl_length, lower = 5)
  
  selected_metric <- setdiff(names(dt_response), 
                             c(cellline_name, "rId", "cId"))
  stopifnot("Provide `dt_response` for one metric." = NROW(selected_metric) == 1)
  
  stopifnot("There is no data in `dt_depmap_lng` (all `selected_meta` is NA)." = !all(is.na(dt_depmap_lng[[selected_meta]]))) # nolint
  stopifnot(
    "It seems that `dt_depmap_lng` has too many categories for `selected_meta` - try use `plot_scatter_with_corr()`." = 
      all(
        NROW(unique(dt_depmap_lng[!is.na(get(selected_meta)), ][[selected_meta]])) < NROW(dt_depmap_lng[!is.na(get(selected_meta)), ]), # nolint
        NROW(unique(dt_depmap_lng[get(selected_meta) != "", ][[selected_meta]])) < NROW(dt_depmap_lng[get(selected_meta) != "", ]) # nolint
      )
  )
  
  N <- CCLEName <- NULL # due to NSE notes in R CMD check
  
  # prep table with data to plot
  X_dt <- dt_depmap_lng[, c("CCLEName", selected_meta), with = FALSE]
  if (is.numeric(X_dt[[selected_meta]])) {
    X_dt[[selected_meta]] <- as.factor(X_dt[[selected_meta]])
  }
  Y_dt <- dt_response[, c(cellline_name, selected_metric), with = FALSE]
  tab_plot <- Y_dt[X_dt, on = .(CellLineName = CCLEName), nomatch = NULL]
  # remove NA
  tab_plot <- stats::na.omit(tab_plot)
  
  # plot without group with only on item
  if (!with_1_item_grp) {
    multi_item_grp <- tab_plot[, .N, by = selected_meta][N > 1, ][[selected_meta]]
    tab_plot <- tab_plot[get(selected_meta) %in% multi_item_grp, ]
  }
  
  # some labels may be too long to see the boxes 
  if (is.character(tab_plot[[selected_meta]]) && 
      any(nchar(unique(tab_plot[[selected_meta]])) > max_x_lbl_length)) {
    too_long_lbl <- which(nchar(tab_plot[[selected_meta]]) > max_x_lbl_length)
    tab_plot[too_long_lbl, ][[selected_meta]] <- 
      paste0(substr(tab_plot[too_long_lbl, ][[selected_meta]], 1, max_x_lbl_length - 3), "...")
  }
  
  # final plt
  plt <-        
    ggplot2::ggplot(
      data = tab_plot,
      mapping =  ggplot2::aes(x = get(selected_meta), y = get(selected_metric))) +
    ggplot2::geom_hline(yintercept = 1, color = "#B3B3B3", linetype = "dashed") +
    ggplot2::geom_hline(yintercept = 0, color = "#B3B3B3", linetype = "solid") +
    ggplot2::geom_boxplot(fill = "#A6CEE3", color = "#A9A9A9", alpha = 0.25) +
    ggplot2::geom_jitter(width = 0.2, height = 0, color = "#4C4C4C") + 
    ggplot2::labs(y = selected_metric, x = "", title = selected_meta) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 90, vjust = 1, hjust = 1))
  
  return(plt)
}

#' Plot panel with volcano plot and scatter plots (top 4) features from DepMap
#'
#' @param dt_response \code{data.table} with experimental response data (rows are samples) for one metric
#'  outputted by one of functions: \code{\link[gDRplots]{prep_dt_response_metric_sa}},
#'  \code{\link[gDRplots]{prep_dt_response_dose_sa}}, \code{\link[gDRplots]{prep_dt_response_scores}}
#'  or \code{\link[gDRplots]{prep_dt_response_metric_diff}}, 
#' @param dt_depmap \code{data.table} with dependent variables data load from DepMap.
#'  (rows are samples, columns are meta);  
#'  outputted by one of \code{\link[gDRplots]{prep_dt_depmap_feat}}
#' @param selected_metric string name of metric in \code{dt_response}
#' @param selected_feat string with name of selected meta in \code{dt_depmap}
#'
#' @return a panel with volcano plot and scatter plots with correlation for feature set
#' 
#' @keywords prism_plots
#' 
#' @export
plot_volcano_corr_panel <- function(dt_response,
                                    dt_depmap,
                                    selected_metric,  
                                    selected_feat) {
  
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  
  checkmate::assert_data_table(dt_response)
  checkmate::assert_data_table(dt_depmap)
  checkmate::assert_string(selected_metric)
  checkmate::assert_string(selected_feat)
  checkmate::assert_names(names(dt_response), must.include = c(cellline_name, selected_metric))
  checkmate::assert_names(names(dt_depmap), must.include = "CCLEName")

  # TODO add validation for NROW(intersect(dt_response[[cellline_name]],  dt_depmap[["CCLEName"]])) == 0
  
  # plot data
  dt_response_ <- dt_response[, c("rId", "cId", cellline_name, selected_metric), with = FALSE]
  
  obj_assoc <- gDRplots::prep_dt_assoc(dt_response = dt_response_,
                                       dt_depmap = dt_depmap,
                                       selected_feat_meta_col = selected_feat)
  # volcano plot
  plt_vol <- plot_volcano_assoc(dt_assoc = obj_assoc[["dt_assoc"]],
                                feature_info = selected_feat)
  
  # scatter plot with corr
  top_4 <- data.table::setorderv(obj_assoc[["dt_assoc"]], cols = "q_value")[["feature"]][1:4]
  ls_plt_corr <- lapply(top_4, function(top_feat) {
    gDRplots::plot_scatter_with_corr(dt_response = dt_response_,
                                     dt_depmap = dt_depmap,
                                     selected_feat = top_feat,
                                     selected_feat_meta_col = selected_feat)
  })
  
  # final panel
  panel <- ggpubr::annotate_figure(
    ggpubr::ggarrange(plotlist = list(plt_vol,
                                      ggpubr::ggarrange(plotlist = ls_plt_corr)), 
                      widths = c(1, 1)),
    top = obj_assoc[["condition_info"]])
  return(panel)
}

#' Plot panel with volcano plot and boxplot for metric values grouped by metadata from DepMap
#'
#' @param dt_response \code{data.table} with experimental response data (rows are samples) for one metric
#'  outputted by one of functions: \code{\link[gDRplots]{prep_dt_response_metric_sa}},
#'  \code{\link[gDRplots]{prep_dt_response_dose_sa}}, \code{\link[gDRplots]{prep_dt_response_scores}}
#'  or \code{\link[gDRplots]{prep_dt_response_metric_diff}}, 
#' @param dt_depmap \code{data.table} with dependent variables data load from DepMap.
#'  (rows are samples, columns are meta);  
#'  outputted by one of \code{\link[gDRplots]{prep_dt_depmap_meta}}
#' @param selected_metric string name of metric in \code{dt_response}
#' @param selected_meta string with name of selected meta in \code{dt_depmap}
#'
#' @return a panel with volcano plot and boxplot for mata
#' 
#' @keywords prism_plots
#' 
#' @export
plot_volcano_box_panel <- function(dt_response,
                                   dt_depmap,
                                   selected_metric,  
                                   selected_meta) {
  
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  
  checkmate::assert_data_table(dt_response)
  checkmate::assert_data_table(dt_depmap)
  checkmate::assert_string(selected_metric)
  checkmate::assert_string(selected_meta)
  checkmate::assert_names(names(dt_response), must.include = c(cellline_name, selected_metric))
  checkmate::assert_names(names(dt_depmap), must.include = "CCLEName")
  
  value <- NULL # due to NSE notes in R CMD check
  
  # TODO add validation for NROW(intersect(dt_response[[cellline_name]],  dt_depmap[["CCLEName"]])) == 0
  
  # plot data
  dt_response_ <- dt_response[, c("rId", "cId", cellline_name, selected_metric), with = FALSE]
  
  obj_assoc <- gDRplots::prep_dt_assoc(dt_response = dt_response_,
                                       dt_depmap = dt_depmap,
                                       selected_feat_meta_col = selected_meta)
  dt_depmap_lng <- 
    data.table::melt(dt_depmap, 
                     id.vars = c("ModelID", "CCLEName"), variable.name = selected_meta
    )[value == 1, ][, value := NULL]
  
  # volcano plot
  plt_vol <- plot_volcano_assoc(dt_assoc = obj_assoc[["dt_assoc"]],
                                feature_info = selected_meta)
  # boxplot
  plt_box <- plot_boxplot_meta(dt_response = dt_response_,
                               dt_depmap_lng = dt_depmap_lng,
                               selected_meta = selected_meta)
  
  # final panel
  panel <- ggpubr::annotate_figure(
    ggpubr::ggarrange(plotlist = list(plt_vol, plt_box), widths = c(1, 1)),
    top = obj_assoc[["condition_info"]])
  return(panel)
}