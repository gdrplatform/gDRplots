context("Test qc_heatmap")

test_that("pheatmap_qc works as expected", {
  mae <- gDRutils::get_synthetic_data("combo_matrix")
  se <- mae[[gDRutils::get_supported_experiments("sa")]][2:5, ]
  dt_average <- gDRutils::convert_se_assay_to_dt(se = se,
                                                 assay_name = "Averaged")
  
  hm_1 <- pheatmap_qc(dt_average = dt_average) # default
  expect_is(hm_1, "pheatmap")
  expect_equal(sort(hm_1$gtable$grobs[[3]]$label), sort(unique(dt_average$clid))) # row names
  expect_equal(sort(hm_1$gtable$grobs[[5]]$label), 
               sort(unique(c(dt_average$Gnumber, dt_average$Gnumber_2)))) # annotation
  expect_is(hm_1[["tree_row"]], "hclust") # dendrogram
  expect_equal(hm_1[["tree_col"]], NA) # no dendrogram
  
  hm_2 <- pheatmap_qc(dt_average = dt_average,
                      normalization_type = "RV",
                      colors_vec = c("darkblue", "grey90"),
                      lbl_by_CellLineName = TRUE,
                      lbl_by_DrugName = TRUE)
  expect_is(hm_2, "pheatmap")
  expect_equal(sort(hm_2$gtable$grobs[[3]]$label), sort(unique(dt_average$CellLineName))) # row names
  expect_equal(sort(hm_2$gtable$grobs[[5]]$label), 
               sort(unique(c(dt_average$DrugName, dt_average$DrugName_2)))) # annotation
  expect_is(hm_2[["tree_row"]], "hclust") # dendrogram
  expect_equal(hm_2[["tree_col"]], NA) # no dendrogram
  
  se <- mae[[gDRutils::get_supported_experiments("combo")]]
  dt_average <- gDRutils::convert_se_assay_to_dt(se = se,
                                                 assay_name = "Averaged")
  hm_3 <- pheatmap_qc(dt_average = dt_average,
                      hm_title = "Combo Data",
                      cluster_rows = FALSE)
  expect_is(hm_3, "pheatmap")
  expect_equal(hm_3$gtable$grobs[[1]]$label, "Combo Data")
  expect_equal(hm_3[["tree_row"]], NA) # no dendrogram
  expect_equal(hm_3[["tree_col"]], NA) # no dendrogram
  
  mae <- gDRutils::get_synthetic_data("combo_codilution")
  se <- mae[[gDRutils::get_supported_experiments("cd")]]
  dt_average <- gDRutils::convert_se_assay_to_dt(se = se, assay_name = "Averaged")
  
  hm_4 <- pheatmap_qc(dt_average = dt_average) # default
  expect_is(hm_4, "pheatmap")
  expect_equal(sort(hm_4$gtable$grobs[[3]]$label), sort(unique(dt_average$clid))) # row names
  expect_equal(sort(hm_4$gtable$grobs[[5]]$label), 
               sort(unique(c(dt_average$Gnumber, dt_average$Gnumber_2)))) # annotation
  expect_is(hm_4[["tree_row"]], "hclust") # dendrogram
  expect_equal(hm_4[["tree_col"]], NA) # no dendrogram
  
  ls_col <- c("#000000", "#F0F0F0")
  hm_5 <- pheatmap_qc(dt_average = dt_average,
                      metric = "x_std",
                      colors_vec = ls_col) 
  expect_is(hm_5, "pheatmap")
  hm_5_data <- hm_5[["gtable"]][["grobs"]][[2]][["children"]][[1]][["gp"]][["fill"]]
  min_val <- data.table::setorderv(data.table::copy(dt_average)[normalization_type == "GR", ], "x_std")[1, ]
  expect_equal(hm_5_data[min_val$clid, 
                         grepl(sprintf("%s_%s_%s_", min_val$Gnumber, min_val$Gnumber_2, min_val$Concentration), 
                               colnames(hm_5_data))],
               ls_col[2]) # check rev
  
  # testing assertions
  expect_error(pheatmap_qc(dt_average = unlist(dt_average)),
               "Assertion on 'dt_average' failed: Must be a data.table")
  expect_error(pheatmap_qc(dt_average = dt_average,
                           normalization_type = "XX"),
               "Assertion on 'normalization_type' failed: Must be element of set")
  expect_error(pheatmap_qc(dt_average = dt_average,
                           metric = "xxx"),
               "Assertion on 'metric' failed: Must be element of set")
  expect_error(pheatmap_qc(dt_average = dt_average,
                           fit_source = 1),
               "Assertion on 'fit_source' failed: Must be of type 'string'")
  expect_error(pheatmap_qc(dt_average = dt_average,
                           hm_title = NULL),
               "Assertion on 'hm_title' failed: Must be of type 'string'")
  expect_error(pheatmap_qc(dt_average = dt_average,
                           colors_vec = 1:3),
               "Assertion on 'colors_vec' failed: Must be of type 'character'")
  expect_error(pheatmap_qc(dt_average = dt_average,
                           colors_vec = c("pinky", "blackish")),
               "Must be a valid color name")
  expect_error(pheatmap_qc(dt_average = dt_average,
                           no_breaks = "str"),
               "Assertion on 'no_breaks' failed: Must be of type 'single integerish value'")
  expect_error(pheatmap_qc(dt_average = dt_average,
                           cluster_rows = 1),
               "Assertion on 'cluster_rows' failed: Must be of type 'logical flag'")
  expect_error(pheatmap_qc(dt_average = dt_average,
                           lbl_by_CellLineName = "str"),
               "Assertion on 'lbl_by_CellLineName' failed: Must be of type 'logical flag'")
  expect_error(pheatmap_qc(dt_average = dt_average,
                           lbl_by_DrugName = 1),
               "Assertion on 'lbl_by_DrugName' failed: Must be of type 'logical flag'")
  
})

test_that("pheatmap_with_anno_sa works as expected", {
  mae <- gDRutils::get_synthetic_data("combo_matrix_small")
  se <- mae[[gDRutils::get_supported_experiments("sa")]]
  dt_metrics <- gDRutils::convert_se_assay_to_dt(se = se, assay_name = "Metrics")
  
  cdata <- SummarizedExperiment::colData(se)
  rdata <- SummarizedExperiment::rowData(se)
  
  res_1 <- data.table::dcast(data = dt_metrics[normalization_type == "GR", ],
                             formula = CellLineName ~ DrugName, value.var = "xc50")
  data.table::setkey(res_1, NULL)
  
  # scenario 1: default
  out_1 <- pheatmap_with_anno_sa(dt_metrics = dt_metrics)
  expect_length(out_1, 2)
  expect_equal(names(out_1), c("data", "heatmap"))
  data_1 <- out_1[["data"]]
  expect_is(data_1, "list")
  expect_equal(names(data_1), c("matrix", "annotation_col", "annotation_row"))
  expect_is(data_1[["matrix"]], "data.table")
  expect_equal(data_1[["matrix"]], res_1)
  expect_equal(data_1[["annotation_col"]], NULL)
  expect_equal(data_1[["annotation_row"]], NULL)
  plt_1 <- out_1[["heatmap"]]
  expect_is(plt_1, "pheatmap")
  expect_equal(plt_1$gtable$grobs[[2]]$label, cdata[["CellLineName"]])
  expect_equal(plt_1$gtable$grobs[[3]]$label, rdata[["DrugName"]])
  expect_true(all(vapply(plt_1$gtable$grobs[[1]]$children[[1]]$gp$fill, is_valid_color, logical(1))))
  
  # scenario 2: selected metric & normalization_type and NA
  dt_metrics_na <- data.table::copy(dt_metrics)
  dt_metrics_na[DrugName %in% c("drug_021", "drug_026")]$x_max <- NA
  res_2 <- data.table::dcast(data = dt_metrics_na[normalization_type == "RV", ],
                             formula = CellLineName ~ DrugName, value.var = "x_max")
  res_2 <- res_2[, .SD, .SDcols = !anyNA]
  data.table::setkey(res_2, NULL)
  
  out_2 <- pheatmap_with_anno_sa(dt_metrics = dt_metrics_na, 
                                 normalization_type = "RV",
                                 metric = "x_max",
                                 hm_title = "X MAX")
  expect_length(out_2, 2)
  expect_equal(names(out_2), c("data", "heatmap"))
  data_2 <- out_2[["data"]]
  expect_is(data_2, "list")
  expect_equal(names(data_2), c("matrix", "annotation_col", "annotation_row"))
  expect_is(data_2[["matrix"]], "data.table")
  expect_equal(data_2[["matrix"]], res_2)
  plt_2 <- out_2[["heatmap"]]
  expect_is(plt_2, "pheatmap")
  expect_equal(plt_2$gtable$grobs[[1]]$label, "X MAX")
  expect_equal(sort(plt_2$gtable$grobs[[6]]$label), 
               sort(unique(dt_metrics_na[!is.na(x_max)]$DrugName))) # no rows with NA
  expect_is(plt_2[["tree_row"]], "hclust") # dendrogram
  expect_is(plt_2[["tree_col"]], "hclust") # dendrogram
  
  
  # scenario 3: annotations for row and col
  annotation_manual_col <- data.table::data.table(
    CellLineName = c("cellline_GB", "cellline_HB"),
    mut_A = c(1, 0),
    mut_B = c("yes", "no"),
    mut_C = c("CC", "BB")
  )
  annotation_manual_row <- data.table::data.table(
    DrugName = c("drug_002", "drug_011", "drug_021", "drug_026"),
    group = c(1, 1, 2, 3)
  )
  annotation_manual_na <- annotation_manual_row[unique(dt_metrics[, "DrugName"]), on = "DrugName"]
  data.table::setorderv(annotation_manual_na, cols = "group", na.last = TRUE)
  annotation_manual_na$group[is.na(annotation_manual_na$group)] <- "NA"
  
  out_3 <- pheatmap_with_anno_sa(dt_metrics = dt_metrics, 
                                 annotation_row = annotation_manual_row,
                                 annotation_col = annotation_manual_col)
  expect_length(out_3, 2)
  expect_equal(names(out_3), c("data", "heatmap"))
  data_3 <- out_3[["data"]]
  expect_is(data_3, "list")
  expect_equal(names(data_3), c("matrix", "annotation_col", "annotation_row"))
  anno_col_3 <- data_3[["annotation_col"]]
  expect_is(anno_col_3, "data.table")
  expect_equal(anno_col_3, annotation_manual_col)
  anno_row_3 <- data_3[["annotation_row"]]
  expect_is(anno_row_3, "data.table")
  expect_equal(anno_row_3, annotation_manual_na)
  expect_is(data_3[["matrix"]], "data.table")
  expect_equal(data_3[["matrix"]], 
               res_1[order(match(CellLineName, anno_col_3$CellLineName)), 
                     c("CellLineName", anno_row_3$DrugName), with = FALSE])
  plt_3 <- out_3[["heatmap"]]
  expect_is(plt_3, "pheatmap")
  expect_equal(plt_3$gtable$grobs[[5]]$label, c("mut_A", "mut_B", "mut_C"))
  expect_equal(plt_3$gtable$grobs[[7]]$label, c("group"))
  expect_true(is.na(plt_3[["tree_row"]])) # no clustering due Inf
  expect_true(is.na(plt_3[["tree_col"]])) # no clustering due Inf
  
  # scenario 4: incomplete annotations for col and color maps
  annotation_map <- list(
    mut_A = c("1" = "coral", "0" = "cadetblue"),
    mut_B = c("yes" = "black", "no" = "grey90"),
    mut_C = c("AA" = "yellow", "BB" = "green")
  )
  annotation_manual_col_na <- data.table::copy(annotation_manual_col)
  annotation_manual_col_na[2, c("mut_A", "mut_B", "mut_C") := "NA"]
  
  out_4 <- pheatmap_with_anno_sa(dt_metrics = dt_metrics, 
                                 annotation_col = annotation_manual_col[1, ],
                                 annotation_colors = annotation_map)
  expect_length(out_4, 2)
  expect_equal(names(out_4), c("data", "heatmap"))
  data_4 <- out_4[["data"]]
  expect_is(data_4, "list")
  expect_equal(names(data_3), c("matrix", "annotation_col", "annotation_row"))
  anno_4 <- out_4[["data"]][["annotation_col"]]
  expect_is(anno_4, "data.table")
  expect_equal(anno_4, annotation_manual_col_na)
  expect_equal(data_4[["matrix"]], 
               res_1[order(match(CellLineName, anno_4$CellLineName))])
  plt_4 <- out_4[["heatmap"]]
  expect_is(plt_4, "pheatmap")
  expect_equal(plt_4$gtable$grobs[[5]]$label, c("mut_A", "mut_B", "mut_C"))
  
  # scenario 5: incomplete annotations for row
  annotation_manual_row <-
    unique(dt_metrics[, c("DrugName", "drug_moa"), with = FALSE])
  anno_test <- data.table::data.table(
    DrugName = c("drug_004", "drug_005", "drug_006"),
    tested_AB = c("yes", "yes", "no")
  )
  annotation_manual_row <- anno_test[annotation_manual_row, on = "DrugName"]
  
  out_5 <- pheatmap_with_anno_sa(dt_metrics = dt_metrics, 
                                 annotation_row = annotation_manual_row)
  expect_length(out_5, 2)
  expect_equal(names(out_5), c("data", "heatmap"))
  data_5 <- out_5[["data"]]
  expect_is(data_5, "list")
  expect_equal(names(data_5), c("matrix", "annotation_col", "annotation_row"))
  expect_is(data_5[["matrix"]], "data.table")
  anno_5 <- out_5[["data"]][["annotation_row"]]
  expect_is(anno_5, "data.table")
  expect_equal(anno_5, annotation_manual_row)
  expect_equal(data_5[["matrix"]], 
               res_1[, c("CellLineName", anno_5$DrugName), with = FALSE])
  plt_5 <- out_5[["heatmap"]]
  expect_is(plt_5, "pheatmap")
  expect_equal(plt_5$gtable$grobs[[5]]$label, c("tested_AB", "drug_moa"))
  expect_true(is.na(plt_5[["tree_row"]])) # no clustering due Inf
  expect_true(is.na(plt_5[["tree_col"]])) # no clustering due Inf
  
  out_6 <- pheatmap_with_anno_sa(dt_metrics = dt_metrics, 
                                 metric = "x_AOC_range",
                                 cluster_cols = FALSE)
  expect_length(out_6, 2)
  expect_equal(names(out_6), c("data", "heatmap"))
  data_6 <- out_6[["data"]]
  expect_is(data_6, "list")
  expect_true(all(vapply(c("annotation_col", "annotation_row"), 
                         function(i) is.null(data_6[[i]]), logical(1))))
  plt_6 <- out_6[["heatmap"]]
  expect_is(plt_6, "pheatmap")
  expect_is(plt_6[["tree_row"]], "hclust") # rows are clustered
  expect_true(is.na(plt_6[["tree_col"]])) # cols aren't clustered
  
  # testing assertions
  expect_error(pheatmap_with_anno_sa(dt_metrics = unlist(dt_metrics)),
               "Assertion on 'dt_metrics' failed: Must be a data.table")
  expect_error(pheatmap_with_anno_sa(dt_metrics = dt_metrics,
                                     normalization_type = "XX"),
               "Assertion on 'normalization_type' failed: Must be element of set")
  expect_error(pheatmap_with_anno_sa(dt_metrics = dt_metrics,
                                     metric = "xxx"),
               "Assertion on 'metric' failed: Must be element of set")
  expect_error(pheatmap_with_anno_sa(dt_metrics = dt_metrics,
                                     fit_source = 1),
               "Assertion on 'fit_source' failed: Must be of type 'string'")
  expect_error(pheatmap_with_anno_sa(dt_metrics = dt_metrics,
                                     hm_title = NULL),
               "Assertion on 'hm_title' failed: Must be of type 'string'")
  expect_error(pheatmap_with_anno_sa(dt_metrics = dt_metrics,
                                     colors_vec = 1:3),
               "Assertion on 'colors_vec' failed: Must be of type 'character'")
  expect_error(pheatmap_with_anno_sa(dt_metrics = dt_metrics,
                                     no_breaks = "str"),
               "Assertion on 'no_breaks' failed: Must be of type 'single integerish value'")
  expect_error(pheatmap_with_anno_sa(dt_metrics = dt_metrics,
                                     cluster_rows = 1),
               "Assertion on 'cluster_rows' failed: Must be of type 'logical flag'")
  expect_error(pheatmap_with_anno_sa(dt_metrics = dt_metrics,
                                     cluster_cols = "yes"),
               "Assertion on 'cluster_cols' failed: Must be of type 'logical flag'")
  expect_error(pheatmap_with_anno_sa(dt_metrics = dt_metrics,
                                     annotation_row = unlist(annotation_manual_row)),
               "Assertion on 'annotation_row' failed: Must be a data.table")
  expect_error(pheatmap_with_anno_sa(dt_metrics = dt_metrics,
                                     annotation_col = unlist(annotation_manual_col)),
               "Assertion on 'annotation_col' failed: Must be a data.table")
  expect_error(pheatmap_with_anno_sa(dt_metrics = dt_metrics,
                                     annotation_colors = unlist(annotation_map)),
               "Assertion on 'annotation_colors' failed: Must be of type 'list'")
})

test_that("pheatmap_with_anno_cd works as expected", {
  mae <- gDRutils::get_synthetic_data("combo_codilution")
  se <- mae[[gDRutils::get_supported_experiments("cd")]]
  dt_metrics <- gDRutils::convert_se_assay_to_dt(se = se, assay_name = "Metrics")
  
  cdata <- SummarizedExperiment::colData(se)
  rdata <- SummarizedExperiment::rowData(se)
  
  res_1 <- data.table::dcast(
    data = dt_metrics[normalization_type == "GR", ],
    formula = CellLineName ~ paste(DrugName, "x", paste0(DrugName_2, "__", Concentration_2)),
    value.var = "xc50")
  data.table::setkey(res_1, NULL)
  
  # scenario 1: default
  out_1 <- pheatmap_with_anno_cd(dt_metrics = dt_metrics)
  expect_length(out_1, 2)
  expect_equal(names(out_1), c("data", "heatmap"))
  data_1 <- out_1[["data"]]
  expect_is(data_1, "list")
  expect_equal(names(data_1), c("matrix", "annotation_col", "annotation_row"))
  expect_is(data_1[["matrix"]], "data.table")
  expect_equal(data_1[["matrix"]], res_1)
  expect_equal(data_1[["annotation_col"]], NULL)
  expect_equal(data_1[["annotation_row"]], NULL)
  plt_1 <- out_1[["heatmap"]]
  expect_is(plt_1, "pheatmap")
  expect_equal(plt_1$gtable$grobs[[2]]$label, cdata[["CellLineName"]])
  expect_equal(
    sort(plt_1$gtable$grobs[[3]]$label), 
    sort(paste(rdata[["DrugName"]], "x", paste0(rdata[["DrugName_2"]], "__", rdata[["Concentration_2"]]))))
  expect_true(all(vapply(plt_1$gtable$grobs[[1]]$children[[1]]$gp$fill, is_valid_color, logical(1))))
  
  annotation_manual_col <-
    unique(dt_metrics[, c("CellLineName", "Tissue"), with = FALSE])
  annotation_manual_row <-
    unique(dt_metrics[, c("DrugName", "DrugName_2", "Concentration_2", "drug_moa", "drug_moa_2"),
                      with = FALSE])
  annotation_map <-
    get_ann_color_map(unique(dt_metrics[, c("Tissue", "drug_moa", "drug_moa_2"), with = FALSE]))
  
  res_2 <- data.table::dcast(
    data = dt_metrics[normalization_type == "RV", ],
    formula = CellLineName ~ paste(DrugName, "x", paste0(DrugName_2, "__", Concentration_2)),
    value.var = "x_max")
  data.table::setkey(res_2, NULL)
  
  out_2 <- pheatmap_with_anno_cd(dt_metrics = dt_metrics, 
                                 metric = "x_max",
                                 normalization_type = "RV",
                                 cluster_cols = FALSE,
                                 annotation_row = annotation_manual_row,
                                 annotation_col = annotation_manual_col,
                                 annotation_colors = annotation_map)
  expect_length(out_2, 2)
  expect_equal(names(out_2), c("data", "heatmap"))
  data_2 <- out_2[["data"]]
  expect_is(data_2, "list")
  expect_equal(names(data_2), c("matrix", "annotation_col", "annotation_row"))
  expect_is(data_2[["matrix"]], "data.table")
  expect_equal(data_2[["matrix"]], 
               res_2[, names(data_2[["matrix"]]), with = FALSE])
  anno_2 <- out_2[["data"]][["annotation_row"]]
  expect_is(anno_2, "data.table")
  expect_equal(anno_2, annotation_manual_row)
  plt_2 <- out_2[["heatmap"]]
  expect_is(plt_2, "pheatmap")
  expect_is(plt_2[["tree_row"]], "hclust") # rows are clustered
  expect_true(is.na(plt_2[["tree_col"]])) # cols aren't clustered
  
  
  annotation_manual_row_res <- data.table::copy(annotation_manual_row)
  annotation_manual_row_res[1:20, c("drug_moa", "drug_moa_2")] <- "NA"
  annotation_manual_row_na <- annotation_manual_row_res[21:NROW(annotation_manual_row_res), ]
  
  out_3 <- pheatmap_with_anno_cd(dt_metrics = dt_metrics, 
                                 metric = "x_max",
                                 normalization_type = "RV",
                                 cluster_cols = FALSE,
                                 cluster_rows = FALSE,
                                 annotation_row = annotation_manual_row_na,
                                 annotation_colors = annotation_map)
  expect_length(out_3, 2)
  expect_equal(names(out_3), c("data", "heatmap"))
  data_3 <- out_3[["data"]]
  expect_is(data_3, "list")
  expect_equal(names(data_3), c("matrix", "annotation_col", "annotation_row"))
  expect_true(is.null(data_3[["annotation_col"]]))
  anno_3 <- out_3[["data"]][["annotation_row"]]
  expect_is(anno_3, "data.table")
  expect_equal(data.table::setorder(anno_3), data.table::setorder(annotation_manual_row_res))
  plt_3 <- out_3[["heatmap"]]
  expect_is(plt_3, "pheatmap")
  expect_true(is.na(plt_3[["tree_row"]])) # rows aren't clustered
  expect_true(is.na(plt_3[["tree_col"]])) # cols aren't clustered
  
  # testing assertions
  expect_error(pheatmap_with_anno_cd(dt_metrics = unlist(dt_metrics)),
               "Assertion on 'dt_metrics' failed: Must be a data.table")
  expect_error(pheatmap_with_anno_cd(dt_metrics = dt_metrics,
                                     normalization_type = "XX"),
               "Assertion on 'normalization_type' failed: Must be element of set")
  expect_error(pheatmap_with_anno_cd(dt_metrics = dt_metrics,
                                     metric = "xxx"),
               "Assertion on 'metric' failed: Must be element of set")
  expect_error(pheatmap_with_anno_cd(dt_metrics = dt_metrics,
                                     fit_source = 1),
               "Assertion on 'fit_source' failed: Must be of type 'string'")
  expect_error(pheatmap_with_anno_cd(dt_metrics = dt_metrics,
                                     hm_title = NULL),
               "Assertion on 'hm_title' failed: Must be of type 'string'")
  expect_error(pheatmap_with_anno_cd(dt_metrics = dt_metrics,
                                     colors_vec = 1:3),
               "Assertion on 'colors_vec' failed: Must be of type 'character'")
  expect_error(pheatmap_with_anno_cd(dt_metrics = dt_metrics,
                                     no_breaks = "str"),
               "Assertion on 'no_breaks' failed: Must be of type 'single integerish value'")
  expect_error(pheatmap_with_anno_cd(dt_metrics = dt_metrics,
                                     cluster_rows = 1),
               "Assertion on 'cluster_rows' failed: Must be of type 'logical flag'")
  expect_error(pheatmap_with_anno_cd(dt_metrics = dt_metrics,
                                     cluster_cols = "yes"),
               "Assertion on 'cluster_cols' failed: Must be of type 'logical flag'")
  expect_error(pheatmap_with_anno_cd(dt_metrics = dt_metrics,
                                     annotation_row = unlist(annotation_manual_row)),
               "Assertion on 'annotation_row' failed: Must be a data.table")
  expect_error(pheatmap_with_anno_cd(dt_metrics = dt_metrics,
                                     annotation_col = unlist(annotation_manual_col)),
               "Assertion on 'annotation_col' failed: Must be a data.table")
  expect_error(pheatmap_with_anno_cd(dt_metrics = dt_metrics,
                                     annotation_colors = unlist(annotation_map)),
               "Assertion on 'annotation_colors' failed: Must be of type 'list'")
})

test_that("pheatmap_with_anno_combo works as expected", {
  mae <- gDRutils::get_synthetic_data("combo_matrix_small")
  se <- mae[[gDRutils::get_supported_experiments("combo")]]
  dt_scores <- gDRutils::convert_se_assay_to_dt(se = se, assay_name = "scores")
  
  cdata <- SummarizedExperiment::colData(se)
  rdata <- SummarizedExperiment::rowData(se)
  
  res_1 <- data.table::dcast(data = dt_scores[normalization_type == "GR", ],
                             formula = CellLineName ~ paste(DrugName, "x", DrugName_2), 
                             value.var = "hsa_score")
  data.table::setkey(res_1, NULL)
  
  # scenario 1: default
  out_1 <- pheatmap_with_anno_combo(dt_scores = dt_scores)
  expect_length(out_1, 2)
  expect_equal(names(out_1), c("data", "heatmap"))
  data_1 <- out_1[["data"]]
  expect_is(data_1, "list")
  expect_equal(names(data_1), c("matrix", "annotation_col", "annotation_row"))
  expect_equal(data_1[["matrix"]], res_1)
  plt_1 <- out_1[["heatmap"]]
  expect_is(plt_1, "pheatmap")
  expect_equal(plt_1$gtable$grobs[[4]]$label, cdata[["CellLineName"]])
  expect_equal(sort(plt_1$gtable$grobs[[5]]$label), 
               sort(sprintf("%s x %s", rdata$DrugName, rdata$DrugName_2)))
  expect_true(all(vapply(plt_1$gtable$grobs[[1]]$children[[1]]$gp$fill, is_valid_color, logical(1))))
  expect_is(plt_1[["tree_row"]], "hclust") # rows are clustered
  expect_is(plt_1[["tree_col"]], "hclust") # cols are clustered
  
  # scenario 2: selected metric & normalization_type and NA
  dt_scores_na <- data.table::copy(dt_scores)
  dt_scores_na[DrugName == "drug_004"]$bliss_score <- NA
  res_2 <- data.table::dcast(data = dt_scores_na[normalization_type == "RV", ],
                             formula = CellLineName ~ paste(DrugName, "x", DrugName_2), 
                             value.var = "bliss_score")
  res_2 <- res_2[, .SD, .SDcols = !anyNA]
  data.table::setkey(res_2, NULL)
  drug_combo_names <- unique(sprintf("%s x %s",
                                     dt_scores_na[!is.na(bliss_score)]$DrugName, 
                                     dt_scores_na[!is.na(bliss_score)]$DrugName_2))
  
  out_2 <- pheatmap_with_anno_combo(dt_scores = dt_scores_na, 
                                    normalization_type = "RV",
                                    metric = "bliss_score",
                                    hm_title = "RV Bliss Score",
                                    cluster_cols = FALSE)
  expect_length(out_2, 2)
  expect_equal(names(out_2), c("data", "heatmap"))
  data_2 <- out_2[["data"]]
  expect_is(data_2, "list")
  expect_equal(names(data_2), c("matrix", "annotation_col", "annotation_row"))
  expect_equal(data_2[["matrix"]], res_2)
  expect_equal(data_2[["annotation_col"]], NULL)
  expect_equal(data_2[["annotation_row"]], NULL)
  plt_2 <- out_2[["heatmap"]]
  expect_is(plt_2, "pheatmap")
  expect_equal(plt_2$gtable$grobs[[1]]$label, "RV Bliss Score")
  expect_equal(sort(plt_2$gtable$grobs[[5]]$label), sort(drug_combo_names))
  expect_is(plt_2[["tree_row"]], "hclust") # rows are clustered
  expect_true(is.na(plt_2[["tree_col"]])) # cols aren't clustered
  
  # scenario 3: annotations for col and color maps
  annotation_manual_col <- data.table::data.table(
    CellLineName = c("cellline_GB", "cellline_HB"),
    mut_A = c(1, 0),
    mut_B = c("yes", "no")
  )
  annotation_map <- list(
    mut_A = c("1" = "coral", "0" = "cadetblue"),
    mut_B = c("yes" = "black", "no" = "grey90"),
    mut_C = c("AA" = "yellow", "BB" = "green")
  )
  
  out_3 <- pheatmap_with_anno_combo(dt_scores = dt_scores, 
                                    cluster_rows = FALSE,
                                    annotation_col = annotation_manual_col,
                                    annotation_colors = annotation_map)
  expect_length(out_3, 2)
  expect_equal(names(out_3), c("data", "heatmap"))
  data_3 <- out_3[["data"]]
  expect_is(data_3, "list")
  expect_equal(names(data_3), c("matrix", "annotation_col", "annotation_row"))
  expect_equal(data_3[["annotation_col"]], annotation_manual_col)
  expect_equal(data_3[["annotation_row"]], NULL)
  plt_3 <- out_3[["heatmap"]]
  expect_is(plt_3, "pheatmap")
  expect_equal(plt_3$gtable$grobs[[6]]$label, c("mut_A", "mut_B"))
  expect_true(is.na(plt_3[["tree_row"]])) # rows aren't clustered
  expect_is(plt_3[["tree_col"]], "hclust") # cols are clustered
  
  # scenario 4: incomplete annotations for row and incomplete color maps
  annotation_map <- list(
    grp_B = c("yes" = "black", "no" = "grey90"),
    grp_C = c("AA" = "yellow", "BB" = "blue")
  )
  ls_combo_col <- c("DrugName", "DrugName_2", "drug_moa", "drug_moa_2")
  annotation_manual_row <- unique(dt_scores[, ls_combo_col, with = FALSE])
  annotation_manual_row[["grp_B"]] <- c("yes", "yes", "yes", "no", "no", NA)
  annotation_manual_row[["grp_C"]] <- c("AA", "AA", "no_check", "BB", "no_check", "BB")
  
  annotation_manual_row_res <- data.table::setorderv(
    data.table::copy(annotation_manual_row), cols = c("drug_moa", "drug_moa_2", "grp_B", "grp_C"), na.last = TRUE)
  annotation_manual_row_res <- annotation_manual_row_res[, lapply(.SD, change_NA_into_char)]
  
  out_4 <- pheatmap_with_anno_combo(dt_scores = dt_scores,
                                    annotation_row = annotation_manual_row,
                                    annotation_colors = annotation_map)
  expect_length(out_4, 2)
  expect_equal(names(out_4), c("data", "heatmap"))
  data_4 <- out_4[["data"]]
  expect_is(data_4, "list")
  expect_equal(names(data_3), c("matrix", "annotation_col", "annotation_row"))
  anno_4 <- out_4[["data"]][["annotation_row"]]
  expect_is(anno_4, "data.table")
  expect_equal(anno_4, annotation_manual_row_res)
  expect_equal(data_4[["matrix"]], 
               res_1[, c("CellLineName", paste(anno_4$DrugName, "x", anno_4$DrugName_2)), with = FALSE])
  plt_4 <- out_4[["heatmap"]]
  expect_is(plt_4, "pheatmap")
  expect_equal(plt_4$gtable$grobs[[7]]$label, c("drug_moa", "drug_moa_2", "grp_B", "grp_C"))
  
  # scenario 5: incomplete annotations for row and incomplete color maps
  annotation_map <- list(
    drug_moa = c(moa_A = "lightblue", moa_B = "steelblue"),
    drug_moa_2 = c(moa_D = "black", moa_E = "grey")
  )
  annotation_manual_row_na <- annotation_manual_row[2:5, ls_combo_col, with = FALSE]
  annotation_manual_row_res <- merge(
    annotation_manual_row_na,
    unique(dt_scores[, c("DrugName", "DrugName_2"), with = FALSE]), 
    by = c("DrugName", "DrugName_2"), all = TRUE)
  data.table::setorderv(annotation_manual_row_res, cols = c("drug_moa", "drug_moa_2"), na.last = TRUE)
  annotation_manual_row_res <- annotation_manual_row_res[, lapply(.SD, change_NA_into_char)]
  
  out_5 <- pheatmap_with_anno_combo(dt_scores = dt_scores, 
                                    annotation_row = annotation_manual_row_na,
                                    annotation_colors = annotation_map)
  expect_length(out_5, 2)
  expect_equal(names(out_5), c("data", "heatmap"))
  data_5 <- out_5[["data"]]
  expect_is(data_5, "list")
  expect_equal(names(data_5), c("matrix", "annotation_col", "annotation_row"))
  expect_equal(data_5[["matrix"]], res_1[, names(data_5[["matrix"]]), with = FALSE])
  anno_5 <- out_5[["data"]][["annotation_row"]]
  expect_is(anno_5, "data.table")
  expect_equal(anno_5, annotation_manual_row_res) 
  plt_5 <- out_5[["heatmap"]]
  expect_is(plt_5, "pheatmap")
  expect_equal(plt_5$gtable$grobs[[7]]$label, c("drug_moa", "drug_moa_2"))
  
  # testing assertions
  expect_error(pheatmap_with_anno_combo(dt_scores = unlist(dt_scores)),
               "Assertion on 'dt_scores' failed: Must be a data.table")
  expect_error(pheatmap_with_anno_combo(dt_scores = dt_scores,
                                        normalization_type = "XX"),
               "Assertion on 'normalization_type' failed: Must be element of set")
  expect_error(pheatmap_with_anno_combo(dt_scores = dt_scores,
                                        metric = "xxx"),
               "Assertion on 'metric' failed: Must be element of set")
  expect_error(pheatmap_with_anno_combo(dt_scores = dt_scores,
                                        fit_source = 1),
               "Assertion on 'fit_source' failed: Must be of type 'string'")
  expect_error(pheatmap_with_anno_combo(dt_scores = dt_scores,
                                        hm_title = NULL),
               "Assertion on 'hm_title' failed: Must be of type 'string'")
  expect_error(pheatmap_with_anno_combo(dt_scores = dt_scores,
                                        colors_vec = 1:3),
               "Assertion on 'colors_vec' failed: Must be of type 'character'")
  expect_error(pheatmap_with_anno_combo(dt_scores = dt_scores,
                                        no_breaks = "str"),
               "Assertion on 'no_breaks' failed: Must be of type 'single integerish value'")
  expect_error(pheatmap_with_anno_combo(dt_scores = dt_scores,
                                        cluster_rows = 1),
               "Assertion on 'cluster_rows' failed: Must be of type 'logical flag'")
  expect_error(pheatmap_with_anno_combo(dt_scores = dt_scores,
                                        cluster_cols = "yes"),
               "Assertion on 'cluster_cols' failed: Must be of type 'logical flag'")
  expect_error(pheatmap_with_anno_combo(dt_scores = dt_scores,
                                        annotation_row = unlist(annotation_manual_row)),
               "Assertion on 'annotation_row' failed: Must be a data.table")
  expect_error(pheatmap_with_anno_combo(dt_scores = dt_scores,
                                        annotation_col = unlist(annotation_manual_col)),
               "Assertion on 'annotation_col' failed: Must be a data.table")
  expect_error(pheatmap_with_anno_combo(dt_scores = dt_scores,
                                        annotation_colors = unlist(annotation_map)),
               "Assertion on 'annotation_colors' failed: Must be of type 'list'")
})

test_that("get_hm_title works as expected", {
  metric <- "xc50"
  normalization_type <- "GR"
  expect_equal(get_hm_title(metric, normalization_type), "log10(GR50)")
  
  metric <- "bliss_score"
  normalization_type <- "RV"
  expect_equal(get_hm_title(metric, normalization_type), "Bliss Score RV")
  
  dataset_name <- "Dataset AB123"
  metric <- "xc50"
  normalization_type <- "RV"
  expect_equal(get_hm_title(metric, normalization_type, dataset_name), "Dataset AB123 (log10(IC50))")
  
  metric <- "hsa_score"
  normalization_type <- "GR"
  expect_equal(get_hm_title(metric, normalization_type, dataset_name), "Dataset AB123 (HSA Score GR)")
  
  metric <- "x_mean_sd_custom_metric"
  normalization_type <- "RV"
  expect_equal(get_hm_title(metric, normalization_type), "RV Mean Sd Custom Metric")
  
  metric <- "aggregated xc50"
  normalization_type <- "RV"
  expect_equal(get_hm_title(metric, normalization_type), "Aggregated IC50")
})

test_that("change_NA_into_char works", {
  test_vec <- list(11, " ", NA, "A", NULL)
  
  expect_equal(change_NA_into_char(test_vec),
               c("11", " ", "NA", "A", "NULL"))
  expect_equal(change_NA_into_char(test_vec, "Unavailable"),
               c("11", " ", "Unavailable", "A", "NULL"))
  
  expect_error(change_NA_into_char(test_vec, lbl_NA = 1),
               "Assertion on 'lbl_NA' failed: Must be of type 'string', not 'double'")
  expect_error(change_NA_into_char(test_vec, lbl_NA = NA),
               "Assertion on 'lbl_NA' failed: May not be NA")
  expect_error(change_NA_into_char(test_vec, lbl_NA = NULL),
               "Assertion on 'lbl_NA' failed: Must be of type 'string', not 'NULL'")
  expect_error(change_NA_into_char(test_vec, lbl_NA = c("test", "Unavailable")),
               "Assertion on 'lbl_NA' failed: Must have length 1.")
})

test_that("get_qual_colors works", {
  
  max_len <- sum(RColorBrewer::brewer.pal.info[
    RColorBrewer::brewer.pal.info$category == "qual" &
      RColorBrewer::brewer.pal.info$colorblind == TRUE, ]$maxcolors)
  expect_equal(NROW(get_qual_colors()), max_len)
  expect_equal(NROW(unique(get_qual_colors())), max_len)
  expect_equal(get_qual_colors(0), "#000000")
  expect_equal(get_qual_colors(5), c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E"))
  expect_equal(NROW(get_qual_colors(42)), 42)
  expect_equal(NROW(unique(get_qual_colors(42))), 42)
  
  expect_error(get_qual_colors("one"),
               "Assertion on 'n' failed: Must be of type 'single integerish value'")
  expect_error(get_qual_colors(c(2, 3)),
               "Assertion on 'n' failed: Must have length 1.")
  expect_error(get_qual_colors(2.3),
               "Assertion on 'n' failed: Must be of type 'single integerish value'")
  expect_error(get_qual_colors(-1),
               "Assertion on 'n' failed: Element 1 is not >= 0.")
})

test_that("get_ann_color_map works", {
  mae <- gDRutils::get_synthetic_data("small")
  se <- mae[[gDRutils::get_supported_experiments("sa")]][2:3, ]
  dt_averaged <- gDRutils::convert_se_assay_to_dt(se = se, assay_name = "Averaged")
  dt_ann <- unique(dt_averaged[, .SD, .SDcols = c("Tissue", "ReferenceDivisionTime")])
  
  result <- get_ann_color_map(dt_ann)
  expect_equal(NROW(result), NCOL(dt_ann))
  expect_equal(names(result), names(dt_ann))
  expect_equal(NROW(result[["Tissue"]]), NROW(unique(dt_ann[["Tissue"]])))
  expect_equal(NROW(result[["ReferenceDivisionTime"]]), NROW(unique(dt_ann[["ReferenceDivisionTime"]])))
  expect_equal(NROW(unlist(result)), 
               sum(dt_ann[, lapply(.SD, data.table::uniqueN), .SDcols = names(dt_ann)]))
  
  # the same lbl for different annotation column
  dt_ann_2 <- data.table::data.table(
    Tissue = c("tissue_x", "tissue_x", "unknown"),
    drug_moa = c("moa_A", "unknown", "unknown"),
    drug_moa_2 = c("unknown", "moa_AB", "moa_CD")
  )
  
  result_2 <- get_ann_color_map(dt_ann_2)
  expect_equal(NROW(result_2), NCOL(dt_ann_2))
  expect_equal(names(result_2), names(dt_ann_2))
  expect_equal(NROW(result_2[["Tissue"]]), NROW(unique(dt_ann_2[["Tissue"]])))
  expect_equal(NROW(result_2[["drug_moa"]]), NROW(unique(dt_ann_2[["drug_moa"]])))
  expect_equal(NROW(result_2[["drug_moa_2"]]), NROW(unique(dt_ann_2[["drug_moa_2"]])))
  expect_equal(NROW(unlist(result_2)), 
               sum(dt_ann_2[, lapply(.SD, data.table::uniqueN), .SDcols = names(dt_ann_2)]))
  
  expect_error(get_ann_color_map(list()),
               "Assertion on 'dt_ann' failed: Must be a data.table")
})

test_that("fill_ann_color_map works", {
  annotation_manual <- data.table::data.table(
    CellLineName = c("cellline_AA", "cellline_EA", "cellline_IB", "cellline_MC", "cellline_BC"),
    mut_A = c(0, 0, 1, 2, 3),
    mut_B = c("yes", "yes", "no", NA, NA),
    mut_C = c("AA", "AA", "AB", NA, "B")
  )
  
  annotation_map <- list(
    mut_A = c("1" = "coral", "0" = "cadetblue"), # should be removed as not enough long
    mut_B = c("yes" = "black", "no" = "grey90", "not_checked" = "lightblue"), 
    # should be filled with NA = "darkred";"not_checked" should be removed
    mut_C = c("AA" = "red", "AB" = "orange", "B" = "yellow") # should be filled with B = "black"
  )
  
  res <- fill_ann_color_map(dt_ann = annotation_manual, 
                            map_ann = annotation_map)
  
  expect_length(res, 2)
  expect_equal(names(res), c("mut_B", "mut_C"))
  expect_equal(res[["mut_B"]], 
               c(annotation_map[["mut_B"]][na.omit(unique(annotation_manual[["mut_B"]]))], "NA" = "darkred"))
  expect_equal(res[["mut_C"]], c(annotation_map[["mut_C"]], "NA" = "black"))
  
  expect_error(fill_ann_color_map(dt_ann = unlist(annotation_manual)),
               "Assertion on 'dt_ann' failed: Must be a data.table")
  expect_error(fill_ann_color_map(dt_ann = annotation_manual, map_ann = NULL),
               "Assertion on 'map_ann' failed: Must be of type 'list'")
})

