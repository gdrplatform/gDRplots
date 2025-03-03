context("Test pheatmap with anno")

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
  expect_equal(hm_1[["tree_col"]], NA) # no dendrogram for cols
  
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
  expect_equal(hm_2[["tree_col"]], NA) # no dendrogram for cols
  
  se <- mae[[gDRutils::get_supported_experiments("combo")]]
  dt_average <- gDRutils::convert_se_assay_to_dt(se = se,
                                                 assay_name = "Averaged")
  hm_3 <- pheatmap_qc(dt_average = dt_average,
                      hm_title = "Combo Data",
                      cluster_rows = FALSE)
  expect_is(hm_3, "pheatmap")
  expect_equal(hm_3$gtable$grobs[[1]]$label, "Combo Data")
  expect_equal(hm_3[["tree_row"]], NA) # no dendrogram due cluster_rows = FALSE
  expect_equal(hm_3[["tree_col"]], NA) # no dendrogram for cols
  
  mae <- gDRutils::get_synthetic_data("combo_codilution")
  se <- mae[[gDRutils::get_supported_experiments("cd")]]
  dt_average <- gDRutils::convert_se_assay_to_dt(se = se, assay_name = "Averaged")
  
  hm_4 <- pheatmap_qc(dt_average = dt_average) # default
  expect_is(hm_4, "pheatmap")
  expect_equal(sort(hm_4$gtable$grobs[[3]]$label), sort(unique(dt_average$clid))) # row names
  expect_equal(sort(hm_4$gtable$grobs[[5]]$label), 
               sort(unique(c(dt_average$Gnumber, dt_average$Gnumber_2)))) # annotation
  expect_is(hm_4[["tree_row"]], "hclust") # dendrogram
  expect_equal(hm_4[["tree_col"]], NA) # no dendrogram for cols
  
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

  # scenario: error is thrown on duplicates
  dt_average_dup <- data.table::rbindlist(list(dt_average, dt_average))
  expect_error(pheatmap_qc(dt_average = dt_average_dup),
               "Unexpected data aggregation")
  
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
                           distfun = "distfun"),
               "Assertion on 'distfun' failed: Must be a function")
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
  dt_averaged <- gDRutils::convert_se_assay_to_dt(se = se, assay_name = "Averaged")
  dt_metrics_capped <- 
    gDRutils::cap_assay_infinities(conc_assay_dt = dt_averaged,
                                   assay_dt = dt_metrics,
                                   experiment_name = gDRutils::get_supported_experiments("sa"),
                                   capping_fold = 5)
  
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
  expect_equal(plt_1$gtable$grobs[[4]]$label, cdata[["CellLineName"]])
  expect_equal(sort(plt_1$gtable$grobs[[5]]$label), sort(rdata[["DrugName"]]))
  expect_true(all(vapply(plt_1$gtable$grobs[[1]]$children[[1]]$gp$fill, is_valid_color, logical(1))))
  expect_is(plt_1[["tree_row"]], "hclust") # clustering despite Inf
  expect_is(plt_1[["tree_col"]], "hclust") # clustering despite Inf
  
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
                                 hm_title = "X MAX",
                                 colors_vec = c("blue", "grey90"))
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
  
  # scenario 3: annotations for row and col and capped metrics
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
                                 dt_metrics_capped = dt_metrics_capped,
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
                     c("CellLineName", anno_row_3$DrugName), with = FALSE]) # origin data
  plt_3 <- out_3[["heatmap"]]
  expect_is(plt_3, "pheatmap")
  expect_equal(plt_3$gtable$grobs[[7]]$label, c("mut_A", "mut_B", "mut_C"))
  expect_equal(plt_3$gtable$grobs[[9]]$label, c("group"))
  expect_is(plt_3[["tree_row"]], "hclust") # clustering despite Inf
  expect_is(plt_3[["tree_col"]], "hclust") # clustering despite Inf
  expect_false(any(is.infinite(as.numeric(
    plt_3[["gtable"]][["grobs"]][[3]][["children"]][[2]][["label"]])))) # capped data
  
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
  expect_equal(plt_4$gtable$grobs[[7]]$label, c("mut_A", "mut_B", "mut_C"))
  expect_is(plt_4[["tree_row"]], "hclust") # clustering despite Inf
  expect_is(plt_4[["tree_col"]], "hclust") # clustering despite Inf
  
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
  expect_equal(plt_5$gtable$grobs[[7]]$label, c("tested_AB", "drug_moa"))
  expect_is(plt_5[["tree_row"]], "hclust") # clustering despite Inf
  expect_is(plt_5[["tree_col"]], "hclust") # clustering despite Inf
  
  # scenario 6: rows are clustered, cols are not clustered
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
  expect_true(is.na(plt_6[["tree_col"]])) # cols aren't clustered due cluster_cols = FALSE
  
  # scenario 7: selected metric & normalization_type and -Inf in the data
  dt_metrics_inf <- data.table::copy(dt_metrics)
  dt_metrics_inf[DrugName %in% c("drug_021", "drug_026")]$x_max <- -Inf
  out_7 <- pheatmap_with_anno_sa(dt_metrics = dt_metrics_inf, 
                                 normalization_type = "RV",
                                 metric = "x_max",
                                 hm_title = "X MAX",
                                 colors_vec = c("orange", "grey90"))
  expect_length(out_7, 2)
  expect_true(any(out_7$data$matrix == -Inf))
  expect_equal(names(out_7), c("data", "heatmap"))
  data_7 <- out_7[["data"]]
  expect_is(data_7, "list")
  expect_equal(names(data_7), c("matrix", "annotation_col", "annotation_row"))
  expect_is(data_7[["matrix"]], "data.table")
  plt_7 <- out_7[["heatmap"]]
  expect_is(plt_7, "pheatmap")
  expect_equal(plt_7$gtable$grobs[[1]]$label, "X MAX")
  expect_is(plt_7[["tree_row"]], "hclust") # clustering despite Inf
  expect_is(plt_7[["tree_col"]], "hclust") # clustering despite Inf
  expect_false(any(as.numeric(plt_7[["gtable"]][["grobs"]][[7]][["children"]][["GRID.text.385"]][["label"]]))) # colbar
  
  # testing assertions
  expect_error(pheatmap_with_anno_sa(dt_metrics = unlist(dt_metrics)),
               "Assertion on 'dt_metrics' failed: Must be a data.table")
  expect_error(pheatmap_with_anno_sa(dt_metrics = dt_metrics,
                                     dt_metrics_capped = unlist(dt_metrics_capped)),
               "Assertion on 'dt_metrics_capped' failed: Must be a data.table")
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
                                     distfun = "distfun"),
               "Assertion on 'distfun' failed: Must be a function")
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
  expect_equal(sort(plt_1$gtable$grobs[[4]]$label), sort(cdata[["CellLineName"]]))
  expect_equal(
    sort(plt_1$gtable$grobs[[5]]$label),
    sort(paste(rdata[["DrugName"]], "x", paste0(rdata[["DrugName_2"]], "__", rdata[["Concentration_2"]]))))
  expect_true(all(vapply(plt_1$gtable$grobs[[3]]$children[[1]]$gp$fill, is_valid_color, logical(1))))
  
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
  expect_true(is.na(plt_2[["tree_col"]])) # cols aren't clustered due cluster_cols = FALSE
  
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
  expect_true(is.na(plt_3[["tree_row"]])) # rows aren't clustered due cluster_cols = FALSE
  expect_true(is.na(plt_3[["tree_col"]])) # cols aren't clustered due cluster_rows = FALSE
  
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
                                     distfun = "distfun"),
               "Assertion on 'distfun' failed: Must be a function")
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
  expect_true(is.na(plt_2[["tree_col"]])) # cols aren't clustered due cluster_cols = FALSE
  
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
  expect_true(is.na(plt_3[["tree_row"]])) # rows aren't clustered due cluster_rows = FALSE
  expect_is(plt_3[["tree_col"]], "hclust") # cols are clustered
  
  # scenario 4: incomplete annotations for row and incomplete color maps
  annotation_map_4 <- list(
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
                                    annotation_colors = annotation_map_4)
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
  expect_is(plt_4[["tree_row"]], "hclust") # rows are clustered
  expect_is(plt_4[["tree_col"]], "hclust") # cols are clustered
  
  # scenario 5: incomplete annotations for row and incomplete color maps
  annotation_map_5 <- list(
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
                                    annotation_colors = annotation_map_5)
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
  expect_is(plt_4[["tree_row"]], "hclust") # rows are clustered
  expect_is(plt_4[["tree_col"]], "hclust") # cols are clustered
  
  # scenario 6: incomplete annotations for col and color maps
  annotation_manual_col <- data.table::data.table(
    CellLineName = c("cellline_XX", "cellline_HB"),
    mut_A = c(1, 0),
    mut_B = c("yes", "no")
  )
  annotation_map <- list(
    mut_A = c("1" = "coral", "0" = "cadetblue")
  )
  
  annotation_manual_row_res <- 
    merge(unique(dt_scores[, c("CellLineName"), with = FALSE]),
          annotation_manual_col, all.x = TRUE)[, lapply(.SD, change_NA_into_char, "NA")]
  data.table::setorderv(annotation_manual_row_res, order = c(-1L))
  
  out_6 <- pheatmap_with_anno_combo(dt_scores = dt_scores, 
                                    cluster_rows = FALSE,
                                    annotation_col = annotation_manual_col,
                                    annotation_colors = annotation_map)
  expect_length(out_6, 2)
  expect_equal(names(out_6), c("data", "heatmap"))
  data_6 <- out_6[["data"]]
  expect_is(data_6, "list")
  expect_equal(names(data_6), c("matrix", "annotation_col", "annotation_row"))
  expect_equal(data_6[["annotation_col"]], annotation_manual_row_res)
  expect_equal(data_6[["annotation_row"]], NULL)
  plt_6 <- out_6[["heatmap"]]
  expect_is(plt_6, "pheatmap")
  expect_equal(plt_6$gtable$grobs[[6]]$label, c("mut_A", "mut_B"))
  expect_true(is.na(plt_6[["tree_row"]])) # rows aren't clustered
  expect_is(plt_6[["tree_col"]], "hclust") # cols are clustered
  
  # scenario 7: NA in matrix
  dt_scores_NA <- data.table::copy(dt_scores)
  dt_scores_NA[3:5, ][["bliss_score"]] <- NA
  
  out_7 <- pheatmap_with_anno_combo(dt_scores = dt_scores_NA, 
                                    metric = "bliss_score",
                                    normalization_type = "RV",
                                    annotation_col = annotation_manual_col,
                                    annotation_colors = annotation_map)
  expect_length(out_7, 2)
  expect_equal(names(out_7), c("data", "heatmap"))
  data_7 <- out_7[["data"]]
  expect_is(data_7, "list")
  expect_equal(names(data_7), c("matrix", "annotation_col", "annotation_row"))
  plt_7 <- out_7[["heatmap"]]
  expect_is(plt_7, "pheatmap")
  expect_is(plt_7[["tree_row"]], "hclust") # rows are clustered
  expect_is(plt_7[["tree_col"]], "hclust") # cols are clustered
  
  # scenario 8: NA in matrix
  dt_scores_Inf <- data.table::copy(dt_scores)
  dt_scores_Inf[c(1:3, 5:7), ][["bliss_score"]] <- Inf
  
  out_8 <- pheatmap_with_anno_combo(dt_scores = dt_scores_Inf, 
                                    metric = "bliss_score",
                                    normalization_type = "RV",
                                    annotation_col = annotation_manual_col,
                                    annotation_colors = annotation_map)
  expect_length(out_8, 2)
  expect_equal(names(out_8), c("data", "heatmap"))
  data_8 <- out_8[["data"]]
  expect_is(data_8, "list")
  expect_equal(names(data_8), c("matrix", "annotation_col", "annotation_row"))
  plt_8 <- out_8[["heatmap"]]
  expect_is(plt_8, "pheatmap")
  expect_is(plt_8[["tree_row"]], "hclust") # rows are clustered
  expect_is(plt_8[["tree_col"]], "hclust") # cols are clustered
  
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
                                        distfun = "distfun"),
               "Assertion on 'distfun' failed: Must be a function")
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


test_that("prep_pheatmap_matrix works as expected", {
  mae <- gDRutils::get_synthetic_data("combo_matrix")
  se_sa <- mae[[gDRutils::get_supported_experiments("sa")]]
  se_combo <- mae[[gDRutils::get_supported_experiments("combo")]]
  dt_metrics <- gDRutils::convert_se_assay_to_dt(se = se_sa,
                                                 assay_name = "Metrics")
  dt_scores <- gDRutils::convert_se_assay_to_dt(se = se_combo,
                                                assay_name = "scores")
  
  mat_1 <- prep_pheatmap_matrix(dt_response = dt_metrics) # default
  expect_is(mat_1, "matrix")
  expect_equal(sort(rownames(mat_1)), sort(unique(dt_metrics[["CellLineName"]])))
  expect_equal(sort(colnames(mat_1)), sort(unique(dt_metrics[["DrugName"]])))
  expect_true(all(dt_metrics[normalization_type == "GR"][["xc50"]] %in% mat_1))
  
  mat_2 <- prep_pheatmap_matrix(dt_response = dt_metrics,
                                metric = "p_value")
  expect_is(mat_2, "matrix")
  expect_true(all(dt_metrics[normalization_type == "GR"][["p_value"]] %in% mat_2))
  
  # scenario: combo data
  mat_3 <- prep_pheatmap_matrix(dt_response = dt_scores,
                                metric = "hsa_score",
                                normalization_type = "RV",
                                experiment_type = gDRutils::get_supported_experiments("combo"))
  expect_is(mat_3, "matrix")
  expect_equal(sort(rownames(mat_3)), sort(unique(dt_scores[["CellLineName"]])))
  expect_equal(sort(colnames(mat_3)), 
               sort(unique(unique(dt_scores[, paste(DrugName, "x", DrugName_2)]))))
  expect_true(all(dt_scores[normalization_type == "RV"][["hsa_score"]] %in% mat_3))
  
  # scenario: codilution data
  mae <- gDRutils::get_synthetic_data("combo_codilution_small")
  se_cd <- mae[[gDRutils::get_supported_experiments("cd")]]
  dt_metrics_cd <- gDRutils::convert_se_assay_to_dt(se = se_cd,
                                                    assay_name = "Metrics")
  
  mat_4 <- prep_pheatmap_matrix(dt_response = dt_metrics_cd,
                                metric = "x_max",
                                normalization_type = "RV",
                                experiment_type = gDRutils::get_supported_experiments("cd"))
  expect_is(mat_4, "matrix")
  expect_equal(sort(rownames(mat_4)), sort(unique(dt_metrics_cd[["CellLineName"]])))
  expect_equal(
    sort(colnames(mat_4)), 
    sort(unique(unique(dt_metrics_cd[, paste0(DrugName, " x ", DrugName_2, "__", Concentration_2)]))))
  expect_true(all(dt_metrics_cd[normalization_type == "RV"][["x_max"]] %in% mat_4))
  
  # scenario: NA-row
  dt_scores_NA <- data.table::copy(dt_scores)
  dt_scores_NA[CellLineName == "cellline_EA"][["bliss_score"]] <- NA
  
  mat_5 <- prep_pheatmap_matrix(dt_response = dt_scores_NA,
                                metric = "bliss_score",
                                normalization_type = "RV",
                                experiment_type = gDRutils::get_supported_experiments("combo"))
  expect_is(mat_5, "matrix")
  expect_false("cellline_EA" %in% rownames(mat_5))
  expect_equal(sort(colnames(mat_5)), 
               sort(unique(unique(dt_scores[, paste(DrugName, "x", DrugName_2)]))))
  expect_true(all(dt_scores[normalization_type == "RV" & is.na(hsa_score), ][["hsa_score"]] %in% mat_5))
  
  # scenario: matrix 0x0 (all values are NAs)
  dt_metrics_NA <- data.table::copy(dt_metrics)
  dt_metrics_NA[normalization_type == "GR"][["xc50"]] <- NA
  
  mat_6 <- prep_pheatmap_matrix(dt_response = dt_metrics_NA) # default
  expect_is(mat_6, "matrix")
  expect_equal(dim(mat_6), c(0, 0))
  
  # scenario: duplicates
  mat_7 <- purrr::quietly(prep_pheatmap_matrix)(dt_response = dt_scores,
                                                metric = "bliss_score")$result
  expect_is(mat_7, "matrix")
  expect_equal(sort(rownames(mat_7)), sort(unique(dt_scores[["CellLineName"]])))
  expect_equal(sort(colnames(mat_7)), sort(unique(dt_scores[["DrugName"]])))
  expect_false(all(dt_scores[normalization_type == "GR"][["bliss_score"]] %in% mat_7))
  expect_true(all(
    dt_scores[normalization_type == "GR", .N, by = c("CellLineName", "DrugName")][["N"]] %in% mat_7))
  
  # testing assertions
  expect_error(prep_pheatmap_matrix(dt_response = as.list(dt_metrics)),
               "Assertion on 'dt_response' failed: Must be a data.table")
  expect_error(prep_pheatmap_matrix(dt_response = dt_metrics,
                                    normalization_type = "XX"),
               "Assertion on 'normalization_type' failed: Must be element of set")
  expect_error(prep_pheatmap_matrix(dt_response = dt_metrics,
                                    metric = "IC50"),
               "Assertion on 'metric' failed: Must be element of set")
  expect_error(prep_pheatmap_matrix(dt_response = dt_metrics,
                                    fit_source = 1),
               "Assertion on 'fit_source' failed: Must be of type 'string'")
  expect_error(prep_pheatmap_matrix(dt_response = dt_metrics,
                                    experiment_type = 1),
               "Assertion on 'experiment_type' failed: Must be element of set")
  expect_error(prep_pheatmap_matrix(dt_response = dt_metrics,
                                    experiment_type = "sa"),
               "Assertion on 'experiment_type' failed: Must be element of set")
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

test_that(".get_pheatmap_cluster_param works as expected", {
  mat <- matrix(1:24, nrow = 4)
  rownames(mat) <- sprintf("row_%s", 1:4)
  colnames(mat) <- sprintf("col_%s", 1:6)
  
  out_1 <- .get_pheatmap_cluster_param(mat)
  expect_is(out_1, "hclust")
  expect_equal(out_1$labels, rownames(mat))
  expect_equal(out_1$dist.method, "euclidean") # default for stats::dist
  expect_equal(out_1$method, "complete")
  
  out_2 <- .get_pheatmap_cluster_param(t(mat))
  expect_is(out_2, "hclust")
  expect_equal(out_2$labels, colnames(mat))
  
  out_3 <- .get_pheatmap_cluster_param(t(mat), distfun = compute_distances)
  expect_is(out_3, "hclust")
  expect_equal(out_3$labels, colnames(mat))
  expect_equal(out_3$dist.method, NULL) # hidden in compute_distances
  expect_equal(out_3$method, "complete")
  
  add_cond <- any(dim(mat) > 10)
  out_4 <- .get_pheatmap_cluster_param(t(mat), additional_condition = add_cond)
  expect_is(out_4, "logical")
  expect_equal(out_4, add_cond)
  
  out_5 <- .get_pheatmap_cluster_param(mat[1, , drop = FALSE]) # only one row
  expect_is(out_5, "logical")
  expect_false(out_5)
  
  out_6 <- .get_pheatmap_cluster_param(mat[, 1, drop = FALSE]) # only one row
  expect_is(out_6, "hclust")
  expect_equal(out_6$labels, rownames(mat))
  
  # scenario: matrix with NA
  mat_NA <- mat
  mat_NA[2, ] <- NA
  
  out_7 <- .get_pheatmap_cluster_param(mat_NA)
  expect_is(out_7, "logical")
  expect_false(out_7) # default stats::dist does not handle NAs
  
  out_8 <- .get_pheatmap_cluster_param(t(mat_NA), distfun = compute_distances)
  expect_is(out_8, "hclust")
  expect_equal(out_8$labels, colnames(mat))
  
  out_9 <- .get_pheatmap_cluster_param(mat_NA, distfun = compute_distances, additional_condition = add_cond)
  expect_is(out_9, "logical")
  expect_equal(out_9, add_cond)
  
  # scenario: matrix with -Inf/Inf
  mat_inf <- mat
  mat_inf[2, 2:3] <- Inf
  mat_inf[4, 3:4] <- -Inf
  
  out_10 <- .get_pheatmap_cluster_param(mat_inf)
  expect_is(out_10, "logical")
  expect_false(out_10) # default stats::dist does not handle Inf
  
  out_11 <- .get_pheatmap_cluster_param(t(mat_inf), distfun = compute_distances)
  expect_is(out_11, "hclust")
  expect_equal(out_11$labels, colnames(mat))
  
  # scenario: matrix with -Inf/Inf and NAs
  mat_inf_na <- mat_inf
  mat_inf_na[3, 1:4] <- NA
  
  out_12 <- .get_pheatmap_cluster_param(mat_inf_na)
  expect_is(out_12, "logical")
  expect_false(out_12) # default stats::dist does not handle Inf and NA
  
  out_13 <- .get_pheatmap_cluster_param(t(mat_inf_na), distfun = compute_distances)
  expect_is(out_13, "hclust")
  expect_equal(out_13$labels, colnames(mat))
  
  # testing assertions
  expect_error(.get_pheatmap_cluster_param(mat_to_cluster = list(mat)),
               "Assertion on 'mat_to_cluster' failed: Must be of type 'matrix'")
  expect_error(.get_pheatmap_cluster_param(mat_to_cluster = matrix(1:6, nrow = 3)),
               "Assertion on 'mat_to_cluster' failed: Must have rownames")
  expect_error(.get_pheatmap_cluster_param(mat_to_cluster = matrix(LETTERS[1:6], nrow = 3)),
               "Assertion on 'mat_to_cluster' failed: Must store numerics")
  expect_error(.get_pheatmap_cluster_param(mat_to_cluster = mat, distfun = "dist"),
               "Assertion on 'distfun' failed: Must be a function")
  expect_error(.get_pheatmap_cluster_param(mat_to_cluster = mat, additional_condition = "yes"),
               "Assertion on 'additional_condition' failed: Must be of type 'logical flag'")
})

test_that(".fill_pheatmap_annotation works as expected", {
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  drug_name_2 <- gDRutils::get_env_identifiers("drug_name2")
  
  # single-agent ----
  mae <- gDRutils::get_synthetic_data("combo_matrix")
  se_sa <- mae[[gDRutils::get_supported_experiments("sa")]]
  dt_metrics <- gDRutils::convert_se_assay_to_dt(se = se_sa, assay_name = "Metrics")
  
  cl_names <- unique(dt_metrics[[cellline_name]])
  d_names <- unique(dt_metrics[[drug_name]])
  
  mat_sa <- prep_pheatmap_matrix(dt_response = dt_metrics,
                                 normalization_type = "GR",
                                 metric = "x_mean",
                                 experiment_type = gDRutils::get_supported_experiments("sa"))
  annotation_manual_row <-
    unique(dt_metrics[, c(drug_name, "drug_moa"), with = FALSE])
  
  annotation_1 <- .fill_pheatmap_annotation(dt_anno = annotation_manual_row,
                                            mat_with_metric = mat_sa,
                                            anno_var = drug_name) # default
  expect_equal(annotation_1, annotation_manual_row)
  
  annotation_manual_col <-
    unique(dt_metrics[, c("CellLineName", "Tissue"), with = FALSE])
  
  annotation_2 <- .fill_pheatmap_annotation(dt_anno = annotation_manual_col,
                                            mat_with_metric = mat_sa,
                                            anno_var = cellline_name) # default
  expect_equal(annotation_2, annotation_manual_col)
  
  # scenario: incomplete annotations for columns
  annotation_manual_col_na <- data.table::copy(annotation_manual_col)[c(4:8), ]
  
  annotation_3 <- .fill_pheatmap_annotation(dt_anno = annotation_manual_col_na,
                                            mat_with_metric = mat_sa,
                                            anno_var = cellline_name)
  res_3 <- data.table::copy(annotation_manual_col)
  res_3[!get(cellline_name) %in% annotation_manual_col_na[[cellline_name]], ][["Tissue"]] <- NA
  res_3[["Tissue"]] <- change_NA_into_char(res_3[["Tissue"]])
  expect_equal(annotation_3[order(get(cellline_name)), ], res_3[order(get(cellline_name)), ])
  
  # scenario: annotation with extra items
  annotation_manual_col_lng <- data.table::data.table(
    CellLineName = c("cellline_XX", "cellline_YY", cl_names),
    mut_A = rep(c(1, 0, 1), length.out = NROW(cl_names) + 2),
    mut_B = rep(c("yes", "no"), length.out = NROW(cl_names) + 2)
  )
  annotation_4 <- .fill_pheatmap_annotation(dt_anno = annotation_manual_col_lng,
                                            mat_with_metric = mat_sa,
                                            anno_var = cellline_name)
  expect_equal(annotation_4[order(get(cellline_name)), ], 
               annotation_manual_col_lng[get(cellline_name) %in% cl_names, ][order(get(cellline_name)), ])
  
  # codilution ----
  mae <- gDRutils::get_synthetic_data("combo_codilution_small")
  se_cd <- mae[[gDRutils::get_supported_experiments("cd")]]
  dt_metrics <- gDRutils::convert_se_assay_to_dt(se = se_cd, assay_name = "Metrics")
  
  mat_cd <- prep_pheatmap_matrix(dt_response = dt_metrics,
                                 normalization_type = "GR",
                                 metric = "x_max",
                                 experiment_type = gDRutils::get_supported_experiments("cd"))
  annotation_manual_cd_col <-
    unique(dt_metrics[, c("CellLineName", "Tissue"), with = FALSE])
  
  annotation_5 <- .fill_pheatmap_annotation(dt_anno = annotation_manual_cd_col,
                                            mat_with_metric = mat_cd,
                                            anno_var = cellline_name) # default
  expect_equal(annotation_5, annotation_manual_cd_col)
  
  annotation_manual_cd_row <- unique(dt_metrics[, c("DrugName", "DrugName_2", "Concentration_2",
                                                    "drug_moa", "drug_moa_2"),
                                                with = FALSE])
  # TODO GDR-2791 annotation_6
  
  annotation_manual_cd_row_na <- data.table::copy(annotation_manual_cd_row)[1:20, ]
  
  # TODO GDR-2791 annotation_7
  
  # combo ----
  mae <- gDRutils::get_synthetic_data("combo_matrix")
  se_combo <- mae[[gDRutils::get_supported_experiments("combo")]]
  dt_scores <- gDRutils::convert_se_assay_to_dt(se = se_combo, assay_name = "scores")
  
  cl_names_combo <- unique(dt_scores[[cellline_name]])
  
  mat_combo <- prep_pheatmap_matrix(dt_response = dt_scores,
                                    normalization_type = "GR",
                                    metric = "bliss_score",
                                    experiment_type = gDRutils::get_supported_experiments("combo"))
  drug_combo_names <- sprintf("%s x %s",
                              unique(dt_scores[, c(drug_name, drug_name_2), with = FALSE])[[drug_name]],
                              unique(dt_scores[, c(drug_name, drug_name_2), with = FALSE])[[drug_name_2]])
  annotation_manual_combo_row <-
    unique(dt_scores[, c("DrugName", "DrugName_2", "drug_moa", "drug_moa_2"), with = FALSE])
  
  annotation_8 <- .fill_pheatmap_annotation(dt_anno = annotation_manual_combo_row,
                                            mat_with_metric = mat_combo,
                                            anno_var = drug_name) 
  
  # scenario: missing anno
  annotation_manual_combo_col_na <- data.table::data.table(
    CellLineName =
      c("cellline_AA", "cellline_EA", "cellline_IB", "cellline_MC", "cellline_BC"),
    mut_A = c(1, 1, 1, 0, 0),
    mut_B = c("yes", "yes", "no", "no", "no")
  )
  
  annotation_9 <- .fill_pheatmap_annotation(dt_anno = annotation_manual_combo_col_na,
                                            mat_with_metric = mat_combo,
                                            anno_var = cellline_name)
  res_9 <- rbind(
    data.table::copy(annotation_manual_combo_col_na), 
    data.table::data.table(
      CellLineName = cl_names_combo[!cl_names_combo %in% annotation_manual_combo_col_na[[cellline_name]]]),
    fill = TRUE
  )
  res_9[, (c("mut_A", "mut_B")) := lapply(.SD, change_NA_into_char, "NA"), .SDcols = c("mut_A", "mut_B")]
  expect_equal(annotation_9[order(get(cellline_name)), ], res_9[order(get(cellline_name)), ])
  
  # scenario: incomplete annotations for columns
  annotation_manual_combo_row_na <- data.table::copy(annotation_manual_combo_row)[c(4:8), ]
  
  annotation_10 <- .fill_pheatmap_annotation(dt_anno = annotation_manual_combo_row_na,
                                            mat_with_metric = mat_combo,
                                            anno_var = drug_name)
  # TODO GDR-2791 
  

  # scenario: error is thrown on duplicates
  ## single-agent data
  expect_error(
    prep_pheatmap_matrix(dt_response = dt_scores, metric = "bliss_score"),
    "Unexpected data aggregation"
  )
  
  ## combination data
  dt_scores_dup <- data.table::rbindlist(list(dt_scores, dt_scores))
  expect_error(
    prep_pheatmap_matrix(
      dt_response = dt_scores_dup,
      metric = "bliss_score",
      experiment_type = "combination"
    ),
    "Unexpected data aggregation"
  )
  
  ## co-dilution data
  dt_average_dup <- data.table::rbindlist(list(dt_average, dt_average))
  expect_error(
    prep_pheatmap_matrix(
      dt_response = dt_average_dup,
      metric = "x",
      experiment_type = "co-dilution"
    ),
    "Unexpected data aggregation"
  )

  # testing assertions
  expect_error(.fill_pheatmap_annotation(dt_anno = as.list(annotation_manual_row),
                                         mat_with_metric = mat_sa,
                                         anno_var = drug_name),
               "Assertion on 'dt_anno' failed: Must be a data.table")
  expect_error(.fill_pheatmap_annotation(dt_anno = annotation_manual_row,
                                         mat_with_metric = data.table::data.table(mat_sa),
                                         anno_var = drug_name),
               "Assertion on 'mat_with_metric' failed: Must be of type 'matrix'")
  expect_error(.fill_pheatmap_annotation(dt_anno = annotation_manual_row,
                                         mat_with_metric = matrix(1:48, nrow = 8),
                                         anno_var = drug_name),
               "Assertion on 'mat_with_metric' failed: Must have rownames")
  expect_error(.fill_pheatmap_annotation(dt_anno = annotation_manual_row,
                                         mat_with_metric = matrix(1:48, nrow = 8,
                                                                  dimnames = list(letters[1:8], NULL)),
                                         anno_var = drug_name),
               "Assertion on 'mat_with_metric' failed: Must have colnames")
  expect_error(.fill_pheatmap_annotation(dt_anno = annotation_manual_row,
                                         mat_with_metric = mat_sa,
                                         anno_var = "str"),
               "Assertion on 'anno_var' failed: Must be element of set")
  expect_error(.fill_pheatmap_annotation(dt_anno = annotation_manual_row,
                                         mat_with_metric = mat_sa,
                                         anno_var = cellline_name),
               "Assertion on 'anno_var' failed: Must be a subset of")
})

test_that(".get_pheatmap_number_color works as expected", {
  mat <- matrix(-14:30, ncol = 5,
                dimnames = list(letters[1:9], LETTERS[1:5]))
  no_breaks <- 15
  breaks <- seq(from = min(mat), to = max(mat), length.out = no_breaks + 1)
  colors_vec <- c("limegreen", "darkblue", "orange")
  hm_colors <- grDevices::colorRampPalette(colors_vec)(no_breaks)
  
  number_color <- .get_pheatmap_number_color(mat_with_metric = mat, 
                                             colors_vec = hm_colors, 
                                             breaks = breaks) # default
  
  expect_equal(.get_pheatmap_number_color(mat_with_metric = mat, 
                                          colors_vec = hm_colors, 
                                          breaks = breaks,
                                          dark_color_font = "pinki"), number_color)
  expect_equal(.get_pheatmap_number_color(mat_with_metric = mat, 
                                          colors_vec = hm_colors, 
                                          breaks = breaks,
                                          light_color_font = "darkly"), number_color)
  expect_equal(unique(c(number_color)), c("black", "white"))
  ls_dark <- which(vapply(hm_colors, is_color_dark, logical(1)))
  dark_range <- breaks[c(min(ls_dark), max(ls_dark) + 1)]
  res <- vapply(c(mat), function(i) {
    dark_range[1] < i & i <= dark_range[2]
  }, FUN.VALUE = logical(1))
  expect_equal(c(number_color), ifelse(res, "white", "black"))
  
  number_color_1 <- .get_pheatmap_number_color(mat_with_metric = mat, 
                                               colors_vec = hm_colors, 
                                               breaks = breaks,
                                               dark_color_font = "yellow",
                                               light_color_font = "blue4")
  expect_equal(c(number_color_1), ifelse(res, "yellow", "blue4"))
  
  # scenario: invalid colors in colors_vec
  hm_col_invalid <- gsub("#", "", hm_colors)
  number_color_2 <- .get_pheatmap_number_color(mat_with_metric = mat, 
                                               colors_vec = hm_col_invalid, 
                                               breaks = breaks)
  expect_equal(number_color_2, "black")
  
  number_color_3 <- .get_pheatmap_number_color(mat_with_metric = mat, 
                                               colors_vec = hm_col_invalid, 
                                               breaks = breaks,
                                               light_color_font = "blue4")
  expect_equal(number_color_3, "blue4")
  
  # scenario: only not dark colors in colors_vec
  hm_col_light <- grDevices::colorRampPalette(c("deeppink", "lightblue"))(no_breaks)
  number_color_4 <- .get_pheatmap_number_color(mat_with_metric = mat, 
                                               colors_vec = hm_col_light, 
                                               breaks = breaks)
  expect_equal(number_color_4, "black")
  
  # scenario: short breaks list
  min_no_breaks <- 2
  min_breaks <- seq(from = min(mat), to = max(mat), length.out = min_no_breaks + 1)
  hm_col_short <- grDevices::colorRampPalette(c("deeppink", "navyblue"))(min_no_breaks)
  number_color_5 <- .get_pheatmap_number_color(mat_with_metric = mat, 
                                               colors_vec = hm_col_short,
                                               breaks = min_breaks)
  dark_range <- list(c(-Inf, min_breaks[2]), 
                     c(min_breaks[2], Inf))[[which(vapply(hm_col_short, is_color_dark, logical(1)))]]
  res_5 <- vapply(c(mat), function(i) {
    dark_range[1] < i & i <= dark_range[2]
  }, FUN.VALUE = logical(1))
  expect_equal(c(number_color_5), ifelse(res_5, "white", "black"))
  
  # scenario: NA in metric matrix
  mat_NA <- mat
  mat_NA[1, 1] <- NA
  mat_NA[1, 5] <- NA
  number_color_6 <- .get_pheatmap_number_color(mat_with_metric = mat_NA, 
                                               colors_vec = hm_col_short,
                                               breaks = min_breaks)
  dark_range <- list(c(-Inf, min_breaks[2]), 
                     c(min_breaks[2], Inf))[[which(vapply(hm_col_short, is_color_dark, logical(1)))]]
  res_6 <- vapply(c(mat_NA), function(i) {
    if (is.na(i)) {
      FALSE
    } else {
      dark_range[1] < i & i <= dark_range[2]
    }
  }, FUN.VALUE = logical(1))
  expect_equal(c(number_color_6), ifelse(res_6, "white", "black"))
  
  # scenario: Inf in metric matrix
  mat_Inf <- mat
  mat_Inf[3, 2:4] <- Inf # color for max value
  mat_Inf[8:9, 5] <- -Inf # color for min value
  number_color_7 <- .get_pheatmap_number_color(mat_with_metric = mat_Inf, 
                                               colors_vec = hm_col_short,
                                               breaks = min_breaks)
  res_7 <- vapply(c(mat_Inf), function(i) {
    dark_range[1] < i & i <= dark_range[2]
  }, FUN.VALUE = logical(1))
  expect_equal(c(number_color_7), ifelse(res_7, "white", "black"))
  
  expect_error(.get_pheatmap_number_color(data.table::data.table(),
                                          colors_vec = hm_colors,
                                          breaks = breaks),
               "Assertion on 'mat_with_metric' failed: Must be of type 'matrix'")
  expect_error(.get_pheatmap_number_color(mat_with_metric = matrix(-14:30, ncol = 5),
                                          colors_vec = hm_colors,
                                          breaks = breaks),
               "Assertion on 'mat_with_metric' failed: Must have rownames.")
  expect_error(.get_pheatmap_number_color(mat_with_metric = matrix(-14:30, ncol = 5,
                                                                   dimnames = list(letters[1:9], NULL)),
                                          colors_vec = hm_colors,
                                          breaks = breaks),
               "Assertion on 'mat_with_metric' failed: Must have colnames")
  expect_error(.get_pheatmap_number_color(mat_with_metric = mat,
                                          colors_vec = LETTERS[seq_along(breaks)],
                                          breaks = breaks),
               "Assertion on 'colors_vec' failed: Must have length")
  expect_error(.get_pheatmap_number_color(mat_with_metric = mat,
                                          colors_vec = hm_colors,
                                          breaks = 2),
               "Assertion on 'breaks' failed: Must have length")
  expect_error(.get_pheatmap_number_color(mat_with_metric = mat,
                                          colors_vec = hm_colors,
                                          breaks = "ten"),
               "Assertion on 'breaks' failed: Must be of type 'numeric'")
  expect_error(.get_pheatmap_number_color(mat_with_metric = mat,
                                          colors_vec = hm_colors,
                                          breaks = breaks,
                                          dark_color_font = 123),
               "Assertion on 'dark_color_font' failed: Must be of type 'string'")
  expect_error(.get_pheatmap_number_color(mat_with_metric = mat,
                                          colors_vec = hm_colors,
                                          breaks = breaks,
                                          light_color_font = 123),
               "Assertion on 'light_color_font' failed: Must be of type 'string'")
})
