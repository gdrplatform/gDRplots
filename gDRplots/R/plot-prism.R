#' Volcano plot with association
#'
#' @param dt_assoc \code{data.table} with the calculated linear association between depmap and metrics
#'     outputted by \code{kaleidoscope::calc_assoc}
#' @param condition_txt string describing experiment condition 
#'     (preferred: \code{"DrugName"}_\code{"Gnumber"}_\code{"drug_moa"}_\code{"Duration"})
#' @param correlation_txt string describing association 
#'     (preferred: \code{metric}_\code{depmap feature of metdata})
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
                               condition_txt,
                               correlation_txt,
                               q_cutoff = 0.05,
                               named_p_top = 10,
                               max_N = NULL) {
  
  checkmate::assert_data_table(dt_assoc)
  checkmate::assert_string(condition_txt)
  checkmate::assert_string(correlation_txt)
  checkmate::assert_number(q_cutoff, lower = 0, upper = 1)
  checkmate::assert_number(named_p_top, lower = 0)
  checkmate::assert_number(max_N, lower = 10, null.ok = TRUE)
  
  x_lbl <- "rho"
  y_lbl <- "neglog_q_value"
  stat_sig <- q_value <- label <- feature <- NULL # due to NSE notes in R CMD check
  
  checkmate::assert_names(names(dt_assoc), must.include = c(x_lbl, y_lbl, "q_value", "feature"))
  
  if (all(is.na(dt_assoc[["q_value"]]))) {
    # empty plot
    plt <- 
      ggplot2::ggplot() + 
      ggplot2::labs(title = paste(correlation_txt, ": all NAs")) +
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
      ggplot2::ggplot(data = tab_plot,
                      mapping = ggplot2::aes(x = get(x_lbl), y = get(y_lbl), label = label, color = stat_sig)) +
      ggplot2::geom_point() +
      ggplot2::scale_x_continuous(trans = "identity", name = x_lbl) +
      ggplot2::scale_y_continuous(trans = "identity", name = y_lbl) +
      ggplot2::scale_color_manual(values = list(yes = "black", no = "#A9A9A9"),
                                  name = "Statistically\nSignificant") +
      ggrepel::geom_text_repel(size = 4, show.legend = FALSE) +
      ggplot2::labs(title = sprintf("%s (%s)", correlation_txt, condition_txt)) +
      ggplot2::theme_bw()
  }
  
  return(plt)
}

#' Plot scatter with correlation
#'
#' @param dt_response \code{data.table} with experimental response data (rows are samples) for one metric
#'  outputted by one of functions: \code{\link[gDRplots]{.prep_dt_response_metric_sa}},
#'  \code{\link[gDRplots]{.prep_dt_response_dose_sa}}, \code{\link[gDRplots]{.prep_dt_response_scores}}
#'  or \code{\link[gDRplots]{.prep_dt_response_metric_diff}}, 
#' @param dt_depmap \code{data.table} with dependent variables data loaded from DepMap - for one
#'    feature or one metadata; (rows are samples, columns are features or meta). 
#' @param selected_feat string with name of selected feature from \code{dt_depmap}
#'
#' @return a scatter plot with correlation
#' @keywords prism_plots
#' 
#' @export
plot_scatter_with_corr <- function(dt_response,
                                   dt_depmap, 
                                   selected_feat) {
  
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  clid <- gDRutils::get_env_identifiers("cellline")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  gnumber <- gDRutils::get_env_identifiers("drug")
  drug_name_2 <- gDRutils::get_env_identifiers("drug_name2")
  gnumber_2 <- gDRutils::get_env_identifiers("drug2")
  
  checkmate::assert_data_table(dt_response)
  checkmate::assert_data_table(dt_depmap)
  checkmate::assert_names(names(dt_depmap), must.include = c("CCLEName", selected_feat))
  checkmate::assert_string(selected_feat)
  
  selected_metric <- setdiff(names(dt_response), 
                             c(cellline_name, clid, drug_name, gnumber, drug_name_2, gnumber_2))
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
  feature_type <- "placeholder"
  plt_title <- 
    sprintf("%s\n corr=%2.2f, slope=%2.2f, intercept=%2.2f", 
            feature_type, c, slope, intercept)
  
  plt <-        
    ggplot2::ggplot(
      data = tab_plot,
      mapping =  ggplot2::aes(x = get(selected_feat), y = get(selected_metric), label = get(cellline_name))) +
    ggplot2::geom_point() +
    ggplot2::scale_x_continuous(trans = "identity", name = selected_feat) +
    ggplot2::scale_y_continuous(trans = "identity", name = selected_metric) +
    ggrepel::geom_text_repel(size = 2) +
    ggplot2::geom_abline(intercept = intercept, slope = slope, color = "red") +   
    ggplot2::labs(title = plt_title) +
    ggplot2::theme_bw()
  
  return(plt)
}


#' Plot boxplot for metric values grouped by metadata from DepMap
#'
#' @param dt_response \code{data.table} with experimental response data (rows are samples) for one metric
#'  outputted by one of functions: \code{\link[gDRplots]{.prep_dt_response_metric_sa}},
#'  \code{\link[gDRplots]{.prep_dt_response_dose_sa}}, \code{\link[gDRplots]{.prep_dt_response_scores}}
#'  or \code{\link[gDRplots]{.prep_dt_response_metric_diff}}, 
#' @param dt_depmap \code{data.table} with dependent variables data loaded from DepMap - for one
#'    metadata; (rows are samples, columns are features or meta). 
#' @param selected_meta string with name of selected meta data from \code{dt_depmap}
#'
#' @return a boxplot
#' @keywords prism_plots
#' 
#' @export
plot_boxplot_meta <- function(dt_response,
                              dt_depmap, 
                              selected_meta) {
  
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  clid <- gDRutils::get_env_identifiers("cellline")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  gnumber <- gDRutils::get_env_identifiers("drug")
  drug_name_2 <- gDRutils::get_env_identifiers("drug_name2")
  gnumber_2 <- gDRutils::get_env_identifiers("drug2")
  
  checkmate::assert_data_table(dt_response)
  checkmate::assert_data_table(dt_depmap)
  checkmate::assert_names(names(dt_depmap), must.include = c("CCLEName", selected_meta))
  checkmate::assert_string(selected_meta)
  
  selected_metric <- setdiff(names(dt_response), 
                             c(cellline_name, clid, drug_name, gnumber, drug_name_2, gnumber_2))
  stopifnot("Provide `dt_response` for one metric" = NROW(selected_metric) == 1)
  
  stopifnot("There is no data in `dt_depmap`" = !all(is.na(dt_depmap[[selected_meta]])))
  stopifnot("It seems that `dt_depmap` has too many categories - try use `plot_scatter_with_corr()`" = 
              NROW(unique(dt_depmap[!is.na(get(selected_meta)), ])) <  NROW(dt_depmap[!is.na(get(selected_meta)), ]))
  
  CCLEName <- NULL # due to NSE notes in R CMD check
  
  # prep table with data to plot
  X_dt <- dt_depmap[, c("CCLEName", selected_meta), with = FALSE]
  Y_dt <- dt_response[, c(cellline_name, selected_metric), with = FALSE]
  tab_plot <- Y_dt[X_dt, on = .(CellLineName = CCLEName), nomatch = NULL]
  # remove NA
  tab_plot <- stats::na.omit(tab_plot)

  plt <-        
    ggplot2::ggplot(
      data = tab_plot,
      mapping =  ggplot2::aes(x = get(selected_meta), y = get(selected_metric))) +
    ggplot2::geom_hline(yintercept = 1, color = "#B3B3B3", linetype = "dashed") +
    ggplot2::geom_hline(yintercept = 0, color = "#B3B3B3", linetype = "solid") +
    ggplot2::geom_boxplot(fill = "#A6CEE3", color = "#A9A9A9", alpha = 0.25) +
    ggplot2::geom_jitter(width = 0.2, height = 0, color = "#4C4C4C") + 
    ggplot2::labs(y = selected_metric, x = selected_meta, title = selected_meta) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 90, vjust = 1, hjust = 1))
  
  return(plt)
}