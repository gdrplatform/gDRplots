#' Volcano plot with association
#'
#' @param dt_assoc \code{data.table} with the calculated linear association between DepMap and metrics
#'     outputted by \code{calc_assoc}
#' @param selected_feat_meta_col string describing the name of the associated feature/metadata from DepMap
#' @param selected_metric string describing the name of the selected metric used the association calculation
#' @param condition_info string describing experiment condition 
#'     (preferred: \code{"DrugName"}_\code{"Gnumber"}_\code{"drug_moa"}_\code{"Duration"})
#' @param alpha numeric cutoff to identify statistically significant correlations
#' @param named_p_top numeric value for p-top statistically significant correlations 
#'     to be labeled on the plot
#' @param max_N numeric value for limit the maximum number of non-statistically 
#'     significant points to plot; for default \code{NULL} all points will be plotted.
#'
#' @return \code{ggplot} object containing a volcano plot with association
#' @keywords prism_plots
#' 
#' @export
plot_volcano_assoc <- function(dt_assoc,
                               selected_feat_meta_col,
                               selected_metric,
                               condition_info = NULL,
                               alpha = 0.05,
                               named_p_top = 10,
                               max_N = NULL) {
  
  checkmate::assert_data_table(dt_assoc)
  checkmate::assert_string(selected_feat_meta_col)
  checkmate::assert_string(selected_metric)
  checkmate::assert_string(condition_info, null.ok = TRUE)
  checkmate::assert_number(alpha, lower = 0, upper = 1)
  checkmate::assert_number(named_p_top, lower = 0)
  checkmate::assert_number(max_N, lower = 10, null.ok = TRUE)
  
  x_lbl <- "rho"
  y_lbl <- "neglog_q_value"
  
  checkmate::assert_names(names(dt_assoc), must.include = c(x_lbl, "q_value", "feature"))
  
  plt_title <- sprintf("%s__%s", selected_metric, selected_feat_meta_col)
  
  if (NROW(dt_assoc) == 0 || all(is.na(dt_assoc[["q_value"]]))) {
    # empty plot
    plt <- 
      ggplot2::ggplot() + 
      ggplot2::labs(title = paste(plt_title, ": all NAs"),
                    subtitle = condition_info,
                    x = x_lbl,
                    y = y_lbl) +
      ggplot2::theme_bw()
  } else {
    tab_plot <- data.table::setorderv(data.table::copy(dt_assoc), cols = "q_value")
    
    # prep column with statistically significant
    tab_plot[, stat_sig :=  data.table::fifelse(q_value <= alpha, "yes", "no")]
    
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
    tab_plot[seq_len(named_p_top_act), label := feature]
    
    # prep data
    tab_plot[[y_lbl]] <- -log10(tab_plot$q_value)
    
    # volcano plot
    plt <- 
      ggplot2::ggplot(data = tab_plot,
                      mapping = ggplot2::aes(x = get(x_lbl), 
                                             y = get(y_lbl), 
                                             label = label, 
                                             color = stat_sig)) +
      ggplot2::geom_point() +
      ggplot2::scale_color_manual(values = list(yes = "black", no = "darkgray"),
                                  name = "Statistically Significant") +
      ggrepel::geom_text_repel(size = 3, max.overlaps = 20, show.legend = FALSE) +
      ggplot2::labs(title = plt_title, 
                    subtitle = condition_info,
                    x = x_lbl,
                    y = y_lbl) +
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
#' @param dt_depmap \code{data.table} with dependent variables data loaded from DepMap - for one feature;
#'  (rows are samples, columns are features). 
#'  outputted by \code{\link[gDRplots]{prep_dt_depmap_feat}}
#' @param selected_feat string with name of selected feature from \code{dt_depmap}
#' @param selected_feat_meta_col string with the name of a feature column in DepMap 
#'  (will be used as a plot title)
#'
#' @return \code{ggplot} object containing a scatter plot with correlation
#' @keywords prism_plots
#' 
#' @export
plot_scatter_with_corr <- function(dt_response,
                                   dt_depmap, 
                                   selected_feat,
                                   selected_feat_meta_col = NULL) {
  
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  
  checkmate::assert_data_table(dt_response)
  checkmate::assert_data_table(dt_depmap)
  checkmate::assert_string(selected_feat)
  checkmate::assert_names(names(dt_depmap), must.include = c("CCLEName", selected_feat))
  checkmate::assert_string(selected_feat_meta_col, null.ok = TRUE)
  
  selected_metric <- setdiff(names(dt_response), 
                             c(cellline_name, "rId", "cId"))
  stopifnot("Provide `dt_response` for one metric." = NROW(selected_metric) == 1)
  
  # prep table with data to plot
  X_dt <- dt_depmap[, c("CCLEName", selected_feat), with = FALSE]
  Y_dt <- dt_response[, c(cellline_name, selected_metric), with = FALSE]
  tab_plot <- Y_dt[X_dt, on = .(CellLineName = CCLEName), nomatch = NULL]
  # remove NA
  tab_plot <- stats::na.omit(tab_plot)
  
  if (NROW(tab_plot) == 0) {
    plt <- 
      ggplot2::ggplot() + 
      ggplot2::labs(title = paste(selected_feat_meta_col, ": all NAs"),
                    x = selected_feat, 
                    y = selected_metric,
                    caption = unique(dt_response$rId)) +
      ggplot2::theme_bw()
  } else {
    # re-calculate correlation, slope and intercept
    # add label for points driving the correlation
    fit <- stats::lm(get(selected_metric) ~ get(selected_feat), tab_plot)
    intercept <- stats::coef(fit)[1]
    slope <- stats::coef(fit)[2]
    r_squared <- summary(fit)$r.squared
    correlation <- sqrt(r_squared)
    
    dist_cooks <- sort(stats::cooks.distance(fit), decreasing = TRUE)
    top_driving_corr <- as.numeric(names(dist_cooks)[seq_len(5)])
    tab_plot$label <- ""
    tab_plot[top_driving_corr, ]$label <- tab_plot[top_driving_corr, ][[cellline_name]] 
    tab_plot$col <- "no"
    tab_plot[top_driving_corr, ]$col <- "yes"
    
    # plot title
    plt_subtitle <- 
      sprintf("corr=%2.2f, slope=%2.2f, intercept=%2.2f", correlation, slope, intercept)
    
    plt <-        
      ggplot2::ggplot(data = tab_plot,
                      mapping =  ggplot2::aes(x = get(selected_feat), 
                                              y = get(selected_metric), 
                                              label = label, 
                                              color = col)) +
      ggplot2::geom_point(shape = 21, fill = "black", size = 1, stroke = 1) +
      ggrepel::geom_text_repel(size = 3, max.overlaps = 20, color = "black") +
      ggplot2::geom_abline(intercept = intercept, slope = slope, color = "red") +   
      ggplot2::labs(title = selected_feat_meta_col, 
                    subtitle = plt_subtitle, 
                    x = selected_feat, 
                    y = selected_metric,
                    caption = unique(dt_response$rId)) +
      ggplot2::theme_bw() +
      ggplot2::guides(color = "none") +
      ggplot2::scale_color_manual(values = c(yes = "red", no = "black"))
  } 
  return(plt)
}


#' Plot panel with scatter with correlation
#' 
#' @inheritParams plot_scatter_with_corr
#' @param selected_feats character vector with names of selected features from \code{dt_depmap}
#'
#' @return \code{ggplot} object containing panel of scatter plot with correlation for selected features
#' @keywords prism_plots
#' 
#' @export
plot_scatter_with_corr_panel <- function(dt_response,
                                         dt_depmap, 
                                         selected_feats,
                                         selected_feat_meta_col = NULL) {
  
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  
  checkmate::assert_data_table(dt_response)
  checkmate::assert_data_table(dt_depmap)
  checkmate::assert_character(selected_feats)
  checkmate::assert_names(names(dt_depmap), must.include = "CCLEName")
  checkmate::assert_string(selected_feat_meta_col, null.ok = TRUE)
  
  selected_metric <- setdiff(names(dt_response), 
                             c(cellline_name, "rId", "cId"))
  stopifnot("Provide `dt_response` for one metric." = NROW(selected_metric) == 1)
  
  if (all(is.na(selected_feats))) {
    plt <- 
      ggplot2::ggplot() + 
      ggplot2::labs(title = paste(selected_feat_meta_col, ": all NAs"),
                    x = "", 
                    y = selected_metric) +
      ggplot2::theme_bw()
  } else {
    available_feats <- setdiff(names(dt_depmap), c("ModelID", "CCLEName"))
    # prep table with data to plot
    X_dt <- dt_depmap[, c("CCLEName", selected_feats[selected_feats %in% available_feats]), with = FALSE]
    Y_dt <- dt_response[, c(cellline_name, selected_metric), with = FALSE]
    tab_plot <- Y_dt[X_dt, on = .(CellLineName = CCLEName), nomatch = NULL]
    
    tab_plot_all <- data.table::data.table()
    
    if (sum(is.na(selected_feats)) > 1) {
      tmp <- data.table::data.table(selected_feats)
      tmp[, N := seq_len(.N), by = selected_feats]
      tmp[, selected_feats := ifelse(is.na(selected_feats), sprintf("%s_%s", selected_feats, N), selected_feats)]
      selected_feats <- tmp$selected_feats
    }
    
    feat_lbl_levels <- selected_feats
    
    for (selected_feat in selected_feats) {
      
      if (selected_feat %chin% available_feats) {
        # remove NA
        tab_plot_sel <- stats::na.omit(tab_plot[, c(cellline_name, selected_metric, selected_feat), with = FALSE])
        
        if (NROW(tab_plot_sel) > 0) { 
          # re-calculate correlation, slope and intercept
          fit <- stats::lm(get(selected_metric) ~ get(selected_feat), tab_plot_sel)
          intercept <- stats::coef(fit)[1]
          slope <- stats::coef(fit)[2]
          r_squared <- summary(fit)$r.squared
          correlation <- sqrt(r_squared)
          # add label for points driving the correlation
          dist_cooks <- sort(stats::cooks.distance(fit), decreasing = TRUE)
          top_driving_corr <- as.numeric(names(dist_cooks)[seq_len(5)])
          tab_plot_sel$label <- ""
          tab_plot_sel[top_driving_corr, ]$label <- tab_plot_sel[top_driving_corr, ][[cellline_name]] 
          tab_plot_sel$col <- "no"
          tab_plot_sel[top_driving_corr, ]$col <- "yes"
          tab_plot_sel$feat_lbl <- selected_feat
          data.table::setnames(tab_plot_sel, selected_feat, "feat_val")
        } else {
          # dummy data when all data is NA
          feat_lbl <- paste(selected_feat, ": all NAs")
          feat_lbl_levels[which(selected_feats == selected_feat)] <- feat_lbl
          
          tab_plot_sel <- data.table::data.table(
            cellline_name = "",
            selected_metric = 0,
            feat_val = 0,
            label = "",
            col = "NA",
            feat_lbl = feat_lbl
          )
          data.table::setnames(tab_plot_sel, 
                               old = c("cellline_name", "selected_metric"), 
                               new = c(cellline_name, selected_metric))
        }
      } else {
        # dummy data required for faceting
        feat_lbl <- paste(selected_feat, ": all NAs")
        feat_lbl_levels[which(selected_feats == selected_feat)] <- feat_lbl
        
        tab_plot_sel <- data.table::data.table(
          cellline_name = "",
          selected_metric = 0,
          feat_val = 0,
          label = "",
          col = "NA",
          feat_lbl = feat_lbl
        )
        data.table::setnames(tab_plot_sel, 
                             old = c("cellline_name", "selected_metric"), 
                             new = c(cellline_name, selected_metric))
      }
      tab_plot_all <- rbind(tab_plot_all, tab_plot_sel)
    }
    # order vis as in selected_feats
    tab_plot_all$feat_lbl <- factor(tab_plot_all$feat_lbl, levels = feat_lbl_levels)
    
    plt <-
      ggplot2::ggplot(data = tab_plot_all,
                      mapping = ggplot2::aes(x = feat_val, 
                                             y = get(selected_metric), 
                                             label = label, 
                                             color = col)) +
      ggplot2::geom_point(mapping = ggplot2::aes(alpha = col), 
                          fill = "black", shape = 21, size = 1, stroke = 1) +
      ggrepel::geom_text_repel(size = 3, max.overlaps = 20, color = "black") +
      ggplot2::labs(title = selected_feat_meta_col, 
                    x = "", 
                    y = selected_metric,
                    caption = unique(dt_response$rId)) +
      ggplot2::theme_bw() +
      ggplot2::guides(color = "none") +
      ggplot2::scale_color_manual(values = c(yes = "red", no = "black", "NA" = "black")) +
      ggplot2::scale_alpha_manual(values = c(yes = 1, no = 1, "NA" = 0)) +
      ggplot2::facet_wrap(~feat_lbl, scales = "free") +
      ggplot2::geom_smooth(ggplot2::aes(x = feat_val, 
                                        y = get(selected_metric)), 
                           color = "red",
                           formula = y ~ x, method = "lm", se = FALSE, inherit.aes = FALSE) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8),
                     axis.text.y = ggplot2::element_text(size = 8),
                     plot.title = ggplot2::element_text(size = 12),
                     panel.grid.minor = ggplot2::element_blank(), 
                     aspect.ratio = 1,
                     strip.background = ggplot2::element_blank(),
                     strip.text = ggplot2::element_text(size = 10, face = "bold", hjust = 0, 
                                                        margin = ggplot2::margin()),
                     legend.position = "none")
  }
  return(plt)
}


#' Plot boxplot for categorical features
#'
#' @param dt_response \code{data.table} with experimental response data (rows are samples) for one metric
#'  outputted by one of functions: \code{\link[gDRplots]{prep_dt_response_metric_sa}},
#'  \code{\link[gDRplots]{prep_dt_response_dose_sa}}, \code{\link[gDRplots]{prep_dt_response_scores}}
#'  or \code{\link[gDRplots]{prep_dt_response_metric_diff}}, 
#' @param dt_depmap \code{data.table} with dependent variables data loaded from DepMap - for one feature;
#'  (rows are samples, columns are features). 
#'  outputted by \code{\link[gDRplots]{prep_dt_depmap_feat}}
#' @param selected_feat string with name of selected feature from \code{dt_depmap}
#' @param selected_feat_meta_col string with the name of a feature column in DepMap (will be used as a plot title)
#'  that has 0-1 values only (categorical character but without relation one-to-one for ids and feature)
#'
#' @return \code{ggplot} object containing boxplots for variable levels
#' @keywords prism_plots
#' 
#' @export
plot_boxplot_num <- function(dt_response,
                             dt_depmap, 
                             selected_feat,
                             selected_feat_meta_col = NULL) {
  
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  
  checkmate::assert_data_table(dt_response)
  checkmate::assert_data_table(dt_depmap)
  checkmate::assert_string(selected_feat)
  checkmate::assert_names(names(dt_depmap), must.include = c("CCLEName", selected_feat))
  checkmate::assert_string(selected_feat_meta_col, null.ok = TRUE)
  boxplot_fill <- 
    gDRutils::get_settings_from_json("BOXPLOT_FILL",
                                     system.file(package = "gDRplots", "settings.json"))
  hline_color <- 
    gDRutils::get_settings_from_json("HLINE_COLOR",
                                     system.file(package = "gDRplots", "settings.json"))
  hline_color <- 
    gDRutils::get_settings_from_json("HLINE_COLOR",
                                     system.file(package = "gDRplots", "settings.json"))
  jitter_poinst_color <- 
    gDRutils::get_settings_from_json("JITTER_POINST_COLOR",
                                     system.file(package = "gDRplots", "settings.json"))
  edge_color <- 
    gDRutils::get_settings_from_json("EDGE_COLOR",
                                     system.file(package = "gDRplots", "settings.json"))
  
  selected_metric <- setdiff(names(dt_response), 
                             c(cellline_name, "rId", "cId"))
  stopifnot("Provide `dt_response` for one metric." = NROW(selected_metric) == 1)
  
  # prep table with data to plot
  X_dt <- dt_depmap[, c("CCLEName", selected_feat), with = FALSE]
  Y_dt <- dt_response[, c(cellline_name, selected_metric), with = FALSE]
  tab_plot <- Y_dt[X_dt, on = .(CellLineName = CCLEName), nomatch = NULL]
  
  # remove NA
  tab_plot <- stats::na.omit(tab_plot)
  
  if (NROW(tab_plot) == 0) {
    plt <- 
      ggplot2::ggplot() + 
      ggplot2::labs(title = paste(selected_feat_meta_col, ": all NAs"),
                    x = selected_feat, 
                    y = selected_metric,
                    caption = unique(dt_response$rId)) +
      ggplot2::theme_bw()
  } else {
    # prep value ranges for y-axis
    range_y <- range(tab_plot[!is.infinite(get(selected_metric)), ][[selected_metric]])
    min_val <- min(c(range_y[1], 0), na.rm = TRUE) - 0.05 * (range_y[2] - range_y[1]) 
    
    tab_plot[[selected_feat]] <- factor(tab_plot[[selected_feat]])
    
    # prep the number of items in each category
    tab_count <- tab_plot[, .N, by = selected_feat][, xlbl := sprintf("%s (%s)", get(selected_feat), N)]
    xlbl <- tab_count$xlbl
    names(xlbl) <- tab_count[[selected_feat]]
    
    plt <- 
      ggplot2::ggplot(data = tab_plot,
                      mapping =  ggplot2::aes(x = get(selected_feat), 
                                              y = get(selected_metric))) +
      ggplot2::geom_hline(yintercept = 0, color = hline_color, linetype = "solid") +
      ggplot2::geom_boxplot(fill = boxplot_fill, color = edge_color, alpha = 0.25, staplewidth = 0.5, 
                            outliers = FALSE, na.rm = TRUE) +
      ggplot2::geom_jitter(width = 0.2, height = 0, color = jitter_poinst_color, na.rm = TRUE) + 
      ggplot2::scale_y_continuous(limits = c(min_val, NA)) + 
      ggplot2::labs(title = selected_feat_meta_col,
                    x = selected_feat,
                    y = selected_metric, 
                    caption = unique(dt_response$rId)) +
      ggplot2::theme_bw() +
      ggplot2::scale_x_discrete(labels = xlbl, drop = FALSE) +
      ggplot2::theme(legend.position = "none",
                     axis.text.x = ggplot2::element_text(size = 8))
  }
  return(plt)
}


#' Plot panel with boxplots
#' 
#' @inheritParams plot_boxplot_num
#' @param selected_feats character vector with names of selected features from \code{dt_depmap}
#' @param ncol number of plot column in panel
#'
#' @return \code{ggplot} object containing boxplots for variable levels (0-1)
#' @keywords prism_plots
#' 
#' @export
plot_boxplot_num_panel <- function(dt_response,
                                   dt_depmap, 
                                   selected_feats,
                                   selected_feat_meta_col = NULL,
                                   ncol = NULL) {
  
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  
  checkmate::assert_data_table(dt_response)
  checkmate::assert_data_table(dt_depmap)
  checkmate::assert_character(selected_feats)
  checkmate::assert_names(names(dt_depmap), must.include = "CCLEName")
  checkmate::assert_string(selected_feat_meta_col, null.ok = TRUE)
  checkmate::assert_int(ncol, lower = 1, null.ok = TRUE)
  boxplot_fill <- 
    gDRutils::get_settings_from_json("BOXPLOT_FILL",
                                     system.file(package = "gDRplots", "settings.json"))
  hline_color <- 
    gDRutils::get_settings_from_json("HLINE_COLOR",
                                     system.file(package = "gDRplots", "settings.json"))
  jitter_poinst_color <- 
    gDRutils::get_settings_from_json("JITTER_POINST_COLOR",
                                     system.file(package = "gDRplots", "settings.json"))
  edge_color <- 
    gDRutils::get_settings_from_json("EDGE_COLOR",
                                     system.file(package = "gDRplots", "settings.json"))
  
  selected_metric <- setdiff(names(dt_response), 
                             c(cellline_name, "rId", "cId"))
  stopifnot("Provide `dt_response` for one metric." = NROW(selected_metric) == 1)
  
  # prep table with data to plot
  available_feats <- setdiff(names(dt_depmap), c("ModelID", "CCLEName"))
  X_dt <- dt_depmap[, c("CCLEName", selected_feats[selected_feats %in% available_feats]), with = FALSE]
  Y_dt <- dt_response[, c(cellline_name, selected_metric), with = FALSE]
  tab_plot <- Y_dt[X_dt, on = .(CellLineName = CCLEName), nomatch = NULL]
  
  if (all(is.na(selected_feats)) || 
      all(is.na(tab_plot[[selected_metric]])) ||
      all(vapply(selected_feats[selected_feats %in% available_feats], 
                 function(nm) all(is.na(tab_plot[[nm]])), FUN.VALUE = logical(1)))) {
    plt <- 
      ggplot2::ggplot() + 
      ggplot2::labs(title = paste(selected_feat_meta_col, ": all NAs"),
                    x = "", 
                    y = selected_metric) +
      ggplot2::theme_bw()
  } else {
    # prep dummy values for x-axis
    ls_lbl <- unique(as.vector(as.matrix(tab_plot[, -c(cellline_name, selected_metric), with = FALSE])))
    dummy_feat_val <- ls_lbl[!is.na(ls_lbl)]
    if (all(is.na(dummy_feat_val))) dummy_feat_val <- 0
    
    tab_plot_all <- data.table::data.table()
    
    if (sum(is.na(selected_feats)) > 1) {
      tmp <- data.table::data.table(selected_feats)
      tmp[, N := seq_len(.N), by = selected_feats]
      tmp[, selected_feats := data.table::fifelse(
        is.na(selected_feats), sprintf("%s_%s", selected_feats, N), selected_feats)]
      selected_feats <- tmp$selected_feats
    }
    
    feat_lbl_levels <- selected_feats
    
    for (selected_feat in selected_feats) {
      
      if (selected_feat %chin% available_feats) {
        # remove NA
        tab_plot_tmp <- stats::na.omit(tab_plot[, c(cellline_name, selected_metric, selected_feat), with = FALSE])
        
        if (NROW(tab_plot_tmp) > 0) { 
          # add label 
          tab_plot_tmp$feat_lbl <- selected_feat
          data.table::setnames(tab_plot_tmp, selected_feat, "feat_val")
        } else {
          # dummy data when all data is NA
          feat_lbl <- paste(selected_feat, ": all NAs")
          feat_lbl_levels[which(selected_feats == selected_feat)] <- feat_lbl
          
          tab_plot_tmp <- data.table::data.table(
            cellline_name = "",
            selected_metric = NA,
            feat_val = dummy_feat_val,
            feat_lbl = feat_lbl
          )
          data.table::setnames(tab_plot_tmp, 
                               old = c("cellline_name", "selected_metric"), 
                               new = c(cellline_name, selected_metric))
        }
      } else {
        # dummy data required for faceting
        feat_lbl <- paste(selected_feat, ": all NAs")
        feat_lbl_levels[which(selected_feats == selected_feat)] <- feat_lbl
        
        tab_plot_tmp <- data.table::data.table(
          cellline_name = "",
          selected_metric = NA,
          feat_val = dummy_feat_val,
          feat_lbl = feat_lbl
        )
        data.table::setnames(tab_plot_tmp, 
                             old = c("cellline_name", "selected_metric"), 
                             new = c(cellline_name, selected_metric))
      }
      tab_plot_all <- rbind(tab_plot_all, tab_plot_tmp)
    }
    
    # prep value ranges for y-axis
    range_y <- range(tab_plot_all[!is.infinite(get(selected_metric)), ][[selected_metric]])
    min_val <- min(c(range_y[1], 0), na.rm = TRUE) - 0.05 * (range_y[2] - range_y[1]) 
    
    # order vis as in selected_feats
    tab_plot_all$feat_lbl <- factor(tab_plot_all$feat_lbl, levels = feat_lbl_levels)
    tab_plot_all$feat_val <- factor(tab_plot_all$feat_val)
    
    # prep the number of items in each category
    tab_count_all <- tab_plot_all[!is.na(get(selected_metric)), .N, by = c("feat_val", "feat_lbl")]
    
    # fill lacking labels
    lack_feat_lbl <- feat_lbl_levels[!feat_lbl_levels %in% tab_count_all$feat_lbl]
    if (NROW(lack_feat_lbl) > 0) {
      tab_count_all <- rbind(tab_count_all,
                             expand.grid(feat_val = unique(tab_count_all$feat_val),
                                         feat_lbl = lack_feat_lbl,
                                         N = 0))
    }
    xlbl <- NULL # fix check
    tab_count_all[, xlbl := sprintf("%s (%s)", feat_val, N)]
    
    tab_plot_all <- merge(tab_plot_all, tab_count_all, by = c("feat_val", "feat_lbl"), all.x = TRUE)
    
    plt <- 
      ggplot2::ggplot(data = tab_plot_all,
                      mapping =  ggplot2::aes(x = xlbl, 
                                              y = get(selected_metric))) +
      ggplot2::geom_hline(yintercept = 0, color = hline_color, linetype = "solid") +
      ggplot2::geom_boxplot(fill = boxplot_fill, color = edge_color, alpha = 0.25, staplewidth = 0.5, 
                            outliers = FALSE, na.rm = TRUE) +
      ggplot2::geom_jitter(width = 0.2, height = 0, color = jitter_poinst_color, na.rm = TRUE) + 
      ggplot2::scale_y_continuous(limits = c(min_val, NA)) + 
      ggplot2::labs(title = selected_feat_meta_col,
                    x = "",
                    y = selected_metric, 
                    caption = unique(dt_response$rId)) +
      ggplot2::theme_bw() +
      ggplot2::scale_x_discrete(drop = FALSE) +
      ggplot2::facet_wrap(~feat_lbl, scales = "free", ncol = ncol) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8),
                     axis.text.y = ggplot2::element_text(size = 8),
                     plot.title = ggplot2::element_text(size = 12),
                     panel.grid.minor = ggplot2::element_blank(), 
                     aspect.ratio = 1,
                     strip.background = ggplot2::element_blank(),
                     strip.text = ggplot2::element_text(size = 10, face = "bold", hjust = 0, 
                                                        margin = ggplot2::margin()),
                     legend.position = "none")
  }
  return(plt)
}

#' Plot boxplot for metric values grouped by metadata from DepMap
#'
#' @param dt_response \code{data.table} with experimental response data (rows are samples) for one metric
#'  outputted by one of functions: \code{\link[gDRplots]{prep_dt_response_metric_sa}},
#'  \code{\link[gDRplots]{prep_dt_response_dose_sa}}, \code{\link[gDRplots]{prep_dt_response_scores}}
#'  or \code{\link[gDRplots]{prep_dt_response_metric_diff}}, 
#' @param dt_depmap \code{data.table} with dependent variables data load from DepMap - for one metadata;
#'  (rows are samples, columns are metadata levels);  
#'  outputted by \code{\link[gDRplots]{prep_dt_depmap_meta}}
#' @param selected_feat_meta_col string with the name of the selected metadata from \code{dt_depmap}
#'  (will be used as a plot title)
#' @param with_1_item_grp logical flag indicating whether to show group with only one item
#' @param max_x_lbl_length numeric value for the maximum number of characters in the x-axis label
#' @param with_inf a logical flag indicating whether infinite values should be shown on boxplots
#'
#' @return \code{ggplot} object containing boxplots for variable levels
#' 
#' @keywords prism_plots
#' 
#' @export
plot_boxplot_meta <- function(dt_response,
                              dt_depmap, 
                              selected_feat_meta_col,
                              with_1_item_grp = TRUE,
                              max_x_lbl_length = 60,
                              with_inf = FALSE) {
  
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  
  checkmate::assert_data_table(dt_response)
  checkmate::assert_data_table(dt_depmap)
  checkmate::assert_string(selected_feat_meta_col)
  checkmate::assert_names(names(dt_depmap), must.include = "CCLEName")
  checkmate::assert_flag(with_1_item_grp)
  checkmate::assert_number(max_x_lbl_length, lower = 5)
  checkmate::assert_flag(with_inf)
  boxplot_fill <- 
    gDRutils::get_settings_from_json("BOXPLOT_FILL",
                                     system.file(package = "gDRplots", "settings.json"))
  hline_color <- 
    gDRutils::get_settings_from_json("HLINE_COLOR",
                                     system.file(package = "gDRplots", "settings.json"))
  jitter_poinst_color <- 
    gDRutils::get_settings_from_json("JITTER_POINST_COLOR",
                                     system.file(package = "gDRplots", "settings.json"))
  edge_color <- 
    gDRutils::get_settings_from_json("EDGE_COLOR",
                                     system.file(package = "gDRplots", "settings.json"))
  
  selected_metric <- setdiff(names(dt_response), 
                             c(cellline_name, "rId", "cId"))
  stopifnot("Provide `dt_response` for one metric." = NROW(selected_metric) == 1)
  
  # prep table with data to plot
  available_meta <- setdiff(names(dt_depmap), c("ModelID", "CCLEName"))
  X_dt <- dt_depmap[, c("CCLEName", available_meta), with = FALSE]
  Y_dt <- dt_response[, c(cellline_name, selected_metric), with = FALSE]
  tab_plot <- Y_dt[X_dt, on = .(CellLineName = CCLEName), nomatch = NULL]
  
  check_lbl <- rowSums(tab_plot[, .SD, .SDcols = available_meta], na.rm = TRUE)
  if (any(check_lbl > 1) ||
      typeof(unlist(tab_plot[, -c(cellline_name, selected_metric), with = FALSE])) != "integer") {
    warning(
      "The data does not appear to be categorical because there is no one-to-one relationship between ids and features."
    )
  }
  
  tab_plot <- data.table::melt(tab_plot,
                               id.vars = c(cellline_name, selected_metric),
                               variable.name = selected_feat_meta_col,
                               variable.factor = FALSE
  )[value == 1, !"value"]
  
  # handle -Inf (NA will be not shown on boxplots when with_inf = FALSE)
  if (!with_inf) {
    tab_plot[is.infinite(get(selected_metric)), (selected_metric) := NA] 
  }
  
  if (NROW(tab_plot) == 0 || all(is.na(tab_plot[[selected_feat_meta_col]]))) {
    plt <- 
      ggplot2::ggplot() + 
      ggplot2::labs(title = paste(selected_feat_meta_col, ": all NAs"),
                    x = "", 
                    y = selected_metric) +
      ggplot2::theme_bw()
  } else {
    # remove NA
    tab_plot <- stats::na.omit(tab_plot)
    
    # plot without group with only on item
    if (!with_1_item_grp) {
      multi_item_grp <- tab_plot[, .N, by = selected_feat_meta_col][N > 1, ][[selected_feat_meta_col]]
      tab_plot <- tab_plot[get(selected_feat_meta_col) %chin% multi_item_grp, ]
    }
    
    # prep value ranges for y-axis
    range_y <- range(tab_plot[!is.infinite(get(selected_metric)), ][[selected_metric]])
    min_val <- min(c(range_y[1], 0), na.rm = TRUE) - 0.05 * (range_y[2] - range_y[1]) 
    
    # prep the number of items in each category
    tab_count <- tab_plot[, .N, by = selected_feat_meta_col]
    
    # final plt
    plt <-        
      ggplot2::ggplot(
        data = tab_plot,
        mapping =  ggplot2::aes(x = get(selected_feat_meta_col), 
                                y = get(selected_metric))) +
      ggplot2::geom_hline(yintercept = 0, color = hline_color, linetype = "solid") +
      ggplot2::geom_boxplot(fill = boxplot_fill, color = edge_color, alpha = 0.25, staplewidth = 0.5, 
                            outliers = FALSE, na.rm = TRUE) +
      ggplot2::geom_jitter(width = 0.2, height = 0, color = jitter_poinst_color, na.rm = TRUE) + 
      ggplot2::geom_text(data = tab_count,
                         mapping = ggplot2::aes(x = get(selected_feat_meta_col), 
                                                y = min_val, 
                                                label = N), 
                         vjust = 0, size.unit = "pt", size = 8) +
      ggplot2::scale_y_continuous(limits = c(min_val, NA)) + 
      ggplot2::labs(title = selected_feat_meta_col,
                    x = "",
                    y = selected_metric, 
                    caption = unique(dt_response$rId)) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none",
                     axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
    
    if (!all(is.na(tab_plot[[selected_metric]])) && 
        data.table::between(max(tab_plot[[selected_metric]], na.rm = TRUE), 0.5, 1.5)) {
      plt <- plt +
        ggplot2::geom_hline(yintercept = 1, color = hline_color, linetype = "dashed")
    }
    
    # some labels may be too long to see the boxes 
    if (is.character(tab_plot[[selected_feat_meta_col]]) && 
        any(nchar(unique(tab_plot[[selected_feat_meta_col]])) > max_x_lbl_length)) {
      too_long_lbl <- which(nchar(tab_plot[[selected_feat_meta_col]]) > max_x_lbl_length)
      
      vec_lbl <- tab_plot[[selected_feat_meta_col]]
      names(vec_lbl) <- tab_plot[[selected_feat_meta_col]]
      vec_lbl[too_long_lbl] <- 
        paste0(substr(tab_plot[too_long_lbl, ][[selected_feat_meta_col]], 1, max_x_lbl_length - 3), "...")
      
      plt <- plt + 
        ggplot2::scale_x_discrete(labels = vec_lbl)
    } 
  }
  return(plt)
}

#' Plot panel with volcano plot and according to the data type - scatter plots or box plots
#'
#' @param dt_response \code{data.table} with the experimental response data (rows are samples) 
#'  for one metric outputted by one of functions: \code{\link[gDRplots]{prep_dt_response_metric_sa}},
#'  \code{\link[gDRplots]{prep_dt_response_dose_sa}}, \code{\link[gDRplots]{prep_dt_response_scores}}
#'  or \code{\link[gDRplots]{prep_dt_response_metric_diff}}, 
#'  must have at least a column with \code{CellLineName} and a numeric column with metric values.
#' @param dt_depmap \code{data.table} with dependent variables data loaded from DepMap
#'  where rows are samples, columns are features/metadata levels;
#'  one of: data for one feature outputted by \code{\link[gDRplots]{prep_dt_depmap_feat}} or data 
#'  or data for one metadata outputted by \code{\link[gDRplots]{prep_dt_depmap_meta}}
#' @param selected_metric string name of the metric in \code{dt_response}
#' @param selected_feat_meta_col string with name of selected feature from \code{dt_depmap} or 
#'  the name of the selected metadata from \code{dt_depmap} - respectively
#'
#' @return \code{ggplot} object containing a panel with volcano plot and  depending on data type:
#'  a scatter plots with correlation for top 4 variables or boxplots for variable levels
#' 
#' @keywords prism_plots
#' 
#' @export
plot_volcano_assoc_panel <- function(dt_response,
                                     dt_depmap,
                                     selected_metric,  
                                     selected_feat_meta_col) {
  
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  
  checkmate::assert_data_table(dt_response)
  checkmate::assert_data_table(dt_depmap)
  checkmate::assert_string(selected_metric)
  checkmate::assert_string(selected_feat_meta_col)
  checkmate::assert_names(names(dt_response), must.include = c(cellline_name, selected_metric))
  checkmate::assert_names(names(dt_depmap), must.include = "CCLEName")
  
  # plot data
  ls_cols <- intersect(names(dt_response), c("rId", "cId", cellline_name, selected_metric))
  dt_response_ <- dt_response[, ls_cols, with = FALSE]
  
  obj_assoc <- prep_dt_assoc(dt_response = dt_response_,
                             dt_depmap = dt_depmap,
                             selected_feat_meta_col = selected_feat_meta_col)
  top_4 <- .get_n_top_asssoc(obj_assoc[["dt_assoc"]])
  
  # volcano plot
  plt_vol <- plot_volcano_assoc(dt_assoc = obj_assoc[["dt_assoc"]],
                                selected_feat_meta_col = obj_assoc[["selected_feat_meta_col"]],
                                selected_metric = obj_assoc[["selected_metric"]]) +
    ggplot2::labs(title = "")
  
  # checking type of data: numeric or categorical
  data_type <- .get_data_type(dt_depmap, desc_col = c("ModelID", "CCLEName"))
  
  if (data_type == "numeric") {
    # scatter plot with corr
    plt_side <- plot_scatter_with_corr_panel(dt_response = dt_response_,
                                             dt_depmap = dt_depmap,
                                             selected_feats = top_4,
                                             selected_feat_meta_col = selected_feat_meta_col) + 
      ggplot2::labs(title = "", caption = "")
  } else {
    # boxplot for categorical &  boxplot for numeric as categorical
    plt_side <- plot_boxplot_num_panel(dt_response = dt_response_,
                                       dt_depmap = dt_depmap,
                                       selected_feats = top_4,
                                       selected_feat_meta_col = selected_feat_meta_col,
                                       ncol = 2) + 
      ggplot2::labs(title = "", caption = "")
  }
  
  # final panel
  panel_title <- 
    ifelse(is.null(obj_assoc[["condition_info"]]),
           sprintf("%s__%s", selected_metric, selected_feat_meta_col),
           sprintf("%s__%s\n%s", selected_metric, selected_feat_meta_col, obj_assoc[["condition_info"]]))
  
  panel <- ggpubr::annotate_figure(
    ggpubr::ggarrange(plotlist = list(plt_vol, plt_side), widths = c(1, 1)),
    top = panel_title)
  # final
  panel
}

#' Check data type
#' 
#' @param dt_ \code{data.table} with dependent variables data in the wide format,
#'    where rows are samples, columns are feature levels
#' @param desc_col a character vector with column names describing the data and which 
#'    do not contain data itself
#' 
#' @return a string describing type of data - "numeric" or "categorical"
#' 
#' @examples 
#' \dontrun{
#' tab_cat <- data.table::data.table(
#'   ID = sprintf("ID_%s", seq_len(5)),
#'   brown = c(0, 1, 1, 0, 0),
#'   blue = c(1, 0, NA, 0, 1),
#'   green = c(0, 0, 0, 1, 0)
#' )
#' .get_data_type(dt_ = tab_cat, desc_col = "ID")
#' 
#' tab_feat <- data.table::data.table(
#'   ID = sprintf("ID_%s", seq_len(5)),
#'   grp = LETTERS[seq_len(5)],
#'   low = c(0, 1, 1, NA, 0),
#'   med = c(1, 1, NA, 0, 1),
#'   high = c(0, 1, 0, 1, 0)
#' )
#' .get_data_type(dt_ = tab_feat, desc_col = c("ID", "grp"))
#' }
#' @keywords prism_plots
.get_data_type <- function(dt_,
                           desc_col = NULL) {
  
  checkmate::assert_data_table(dt_)
  checkmate::assert_character(desc_col, null.ok = TRUE)
  if (!is.null(desc_col)) checkmate::assert_subset(desc_col, names(dt_))
  
  # column names with features
  ls_col <- names(dt_)[!names(dt_) %chin% desc_col]
  
  if (all(vapply(dt_[, c(ls_col), with = FALSE], is.numeric, logical(1)))) {
    # checking whether relation in one-to-one or one-to-many
    one_to_one <- !any(rowSums(dt_[, .SD, .SDcols = ls_col], na.rm = TRUE) > 1)
    # unique values
    unique_val <- unique(unlist(lapply(ls_col, function(nm) unique(dt_[[nm]]))))
    
    data_type <- if (is.numeric(unique_val) && all(unique_val %in% c(0, 1, NA))) {
      # assumption: the presence of a feature is described by 0-1; NA means lack of information
      #     categorical - when one id has only one cat
      #     num_as_cat - when one id has many cat
      ifelse(one_to_one, "categorical", "num_as_cat")
    } else {
      "numeric"
    } 
  } else {
    data_type <- "unknown"
  }
  return(data_type)  
}

#' Get n-top linear associations
#' 
#' Currently, n-top linear associations are selected based on values of \code{q_value} 
#' (column sorted increasing) and values of \code{rho} (column sorted decreasing)
#' 
#' @param dt_assoc \code{data.table} with the calculated linear association between DepMap and metrics
#'     outputted by \code{calc_assoc}
#' @param n_top number of requested top linear associations
#' 
#' @return a vector with name of the n-top linear associations; 
#'     when \code{n_top} will be higher than number of available features - only available will be returned.
#' 
#' @keywords prism_plots
.get_n_top_asssoc <- function(dt_assoc,
                              n_top = 4) {
  
  checkmate::assert_data_table(dt_assoc)
  checkmate::assert_names(names(dt_assoc), must.include = c("feature", "q_value", "rho"))
  checkmate::assert_number(n_top, lower = 1) 
  
  dt_ <- data.table::copy(dt_assoc)
  dt_$abs_rho <- abs(dt_$rho)
  dt_ <- data.table::setorderv(dt_, cols = c("q_value", "abs_rho"), order = c(1L, -1L))[, abs_rho := NULL]
  top_n <- dt_[["feature"]][seq_len(n_top)]
  
  return(top_n[!is.na(top_n)])  
}
