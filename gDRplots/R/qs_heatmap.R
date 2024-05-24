mae <- gDRutils::get_synthetic_data("combo_matrix")

se <- mae[[gDRutils::get_supported_experiments("sa")]]


response_metrics <- gDRutils::convert_se_assay_to_dt(se = se, assay_name = "Metrics")
response_data <- gDRutils::convert_se_assay_to_dt(se = se, assay_name = "Averaged")

# groups <- c("normalization_type", "fit_source")
# wide_cols <- gDRutils::get_header("response_metrics")
# # dt_QCS_sa_metrics_0 <- gDRutils::flatten(dt_QCS_sa_metrics, groups = groups, wide_cols = wide_cols)
# response_metrics_flat <- gDRutils::flatten(tbl = response_metrics, 
#                                            groups = groups, 
#                                            wide_cols = wide_cols)


heatmap_QS <- function(tab_response,
                       normalization_type = "GR",
                       metric = "xc50",
                       fit_source = "gDR",
                       dataset_name = NULL,
                       mapcolor = c("firebrick2", "white"),
                       no_breaks = 50) {
  
  checkmate::assert_data_table(tab_response)
  checkmate::assert_choice(normalization_type, choices = c("GR", "RV"))
  checkmate::assert_choice(metric, choices = c("xc50", "x_max", "x_mean"))
  checkmate::assert_string(fit_source, null.ok = TRUE)
  checkmate::assert_string(dataset_name, null.ok = TRUE)
  checkmate::assert_character(mapcolor) # check valid name
  checkmate::assert_int(no_breaks, lower = 2) # check valid name
  
  data.table::setkeyv(tab_response, "normalization_type")
  tab_response <- tab_response[normalization_type]
  data.table::setkey(tab_response, NULL)
  
  qmfun <- switch(metric,
                  "xc50" = log10, 
                  "x_max" = identity, 
                  "x_mean" = identity)
  
  # select data for normalization type
  tab_plot <- data.table::dcast(
    tab_response, 
    factor(CellLineName, levels = unique(CellLineName)) ~  factor(DrugName), 
    value.var = metric)
  
  # prep matrix
  mat_cvd <- as.matrix(tab_plot[,-c("CellLineName")])
  rownames(mat_cvd) <- tab_plot$CellLineName
  rm_col <- vapply(colnames(mat_cvd), function(i) !all(is.na(mat_cvd[,i])), logical(1))
  rm_row <- vapply(seq_along(rownames(mat_cvd)), function(i) !all(is.na(mat_cvd[i,])), logical(1))
  if (!all(rm_col)) mat_cvd <- mat_cvd[,rm_col]
  if (!all(rm_row)) mat_cvd <- mat_cvd[rm_row,]
  mat_cvd[] <- vapply(mat_cvd, function(x) qmfun(x), numeric(1))
  
  # flip 
  t_mat_cvd <- t(mat_cvd)
  
  # prep hmcolors
  breaks <- seq(from = min(na.omit(mat_cvd)), to = 1.0, length.out = no_breaks)
  hmcol <- grDevices::colorRampPalette(mapcolor)(no_breaks + 1)

  hm <- pheatmap::pheatmap(t_mat_cvd,
                           scale = "none",
                           display_numbers = TRUE, 
                           fontsize_number = 6,
                           number_color = "black", 
                           color = hmcol,
                           breaks = breaks, 
                           angle_col = 90, 
                           fontsize = 6,
                           treeheight_row = 30, 
                           treeheight_col = 30,
                           # TODO manual annotation
                           # annotation_col = CL_ann_sa_ready,
                           # annotation_colors = col_ann_sa, 
                           main = paste(metric, dataset_name),
                           cluster_rows = FALSE,
                           cluster_cols = FALSE
  )
}
