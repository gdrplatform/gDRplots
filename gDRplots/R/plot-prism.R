#' Volcano plot with association
#'
#' @param dt_assoc \code{data.table} with the calculated linear association between depmap and metrics
#'     outputted by \code{kaleidoscope::calc_assoc}
#' @param condition_txt string describing experiment condition 
#'     (prefered: \code{"DrugName"}_\code{"Gnumber"}_\code{"drug_moa"}_|code{"Duration"})
#' @param correlation_txt string describing association 
#'     (prefered: \code{metric}_\code{depmap feature of metdata})
#' @param q_cutoff numeric cutoff to identify statistically significant correlations
#' @param named_p_top numeric value for p-top statistically significant correlations to be labeled on the plot
#' @param max_N numeric value for limit maximum number of non-statistically significant points to plot; 
#'     for default \code{NULL} all points will be plotted.
#'
#' @return a volcano plot with 
#' @keywords prism_plots
#' 
#' @export
plot_volcano_assoc <- function(dt_assoc,
                               condition_txt,
                               correlation_txt,
                               q_cutoff = 0.05,
                               named_p_top = 10,
                               max_N = NULL
) {
  
  checkmate::assert_data_table(dt_assoc)
  checkmate::assert_string(condition_txt)
  checkmate::assert_string(correlation_txt)
  checkmate::assert_number(q_cutoff, lower = 0, upper = 1)
  checkmate::assert_number(named_p_top, lower = 0)
  checkmate::assert_number(max_N, lower = 10, null.ok = TRUE)
  
  x_lbl <- "rho"
  y_lbl <- "neglog_q_value"
  
  checkmate::assert_names(names(dt_assoc), must.include = c(x_lbl, y_lbl, "q_value", "feature"))
  
  if (all(is.na(dt_assoc[["q_value"]]))) {
    # empty plot
    plt <- 
      ggplot2::ggplot() + 
      ggplot2::labs(title = paste(correlation_txt, ': all NAs')) +
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
      ggplot2::scale_x_continuous(trans = 'identity', name = x_lbl) +
      ggplot2::scale_y_continuous(trans = 'identity', name = y_lbl) +
      ggplot2::scale_color_manual(values = list(yes = "black", no = "#A9A9A9"),
                                  name = "Statistically\nSignificant") +
      ggrepel::geom_text_repel(size = 4, show.legend = FALSE) +
      ggplot2::labs(title = sprintf("%s (%s)", correlation_txt, condition_txt)) +
      ggplot2::theme_bw()
  }
  
  return(plt)
}
