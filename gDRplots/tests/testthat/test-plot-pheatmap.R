context("Test qc_heatmap")

test_that("pheatmap_with_anno_sa works as expected", {
  mae <- gDRutils::get_synthetic_data("combo_matrix_small")
  se <- mae[[gDRutils::get_supported_experiments("sa")]]
  response_metrics <- gDRutils::convert_se_assay_to_dt(se = se, assay_name = "Metrics")
  
  cdata <- SummarizedExperiment::colData(se)
  rdata <- SummarizedExperiment::rowData(se)
  
  res_1 <- data.table::dcast(data = response_metrics[normalization_type == "GR", ],
                             formula = CellLineName ~ DrugName, value.var = "xc50")
  data.table::setkey(res_1, NULL)
  
  out_1 <- pheatmap_with_anno_sa(tab_response = response_metrics) # default
  expect_length(out_1, 2)
  expect_equal(names(out_1), c("matrix", "heatmap"))
  data_1 <- out_1[["matrix"]]
  expect_is(data_1, "data.table")
  expect_equal(data_1, res_1)
  plt_1 <- out_1[["heatmap"]]
  expect_is(plt_1, "pheatmap")
  expect_equal(plt_1$gtable$grobs[[2]]$label, cdata[["CellLineName"]])
  expect_equal(plt_1$gtable$grobs[[3]]$label, rdata[["DrugName"]])
  expect_true(all(vapply(plt_1$gtable$grobs[[1]]$children[[1]]$gp$fill, is_valid_color, logical(1))))
  
  response_metrics_na <- data.table::copy(response_metrics)
  response_metrics_na[DrugName %in% c("drug_021", "drug_026")]$x_max <- NA
  res_2 <- data.table::dcast(data = response_metrics_na[normalization_type == "RV", ],
                             formula = CellLineName ~ DrugName, value.var = "x_max")
  res_2 <- res_2[, .SD, .SDcols = !anyNA]
  data.table::setkey(res_2, NULL)
  
  out_2 <- pheatmap_with_anno_sa(tab_response = response_metrics_na, 
                                 normalization_type = "RV",
                                 metric = "x_max")
  expect_length(out_2, 2)
  expect_equal(names(out_2), c("matrix", "heatmap"))
  data_2 <- out_2[["matrix"]]
  expect_is(data_2, "data.table")
  expect_equal(data_2, res_2)
  plt_2 <- out_2[["heatmap"]]
  expect_is(plt_2, "pheatmap")
  expect_equal(plt_2$gtable$grobs[[3]]$label, unique(response_metrics_na[!is.na(x_max)]$DrugName))
  
  annotation_manual_col <- data.table::data.table(
    CellLineName = c("cellline_GB", "cellline_HB"),
    mut_A = c(1, 0),
    mut_B = c("yes", "no"),
    mut_C = c("CC", "CC")
  )
  annotation_manual_row <- data.table::data.table(
    DrugName = c("drug_002", "drug_011", "drug_021", "drug_026"),
    group = c(1, 1, 2, 3)
  )
  annotation_manual_na <- annotation_manual_row[unique(response_metrics[, "DrugName"]), on = "DrugName"]
  data.table::setorderv(annotation_manual_na, cols = "group", na.last = TRUE)
  annotation_manual_na$group[is.na(annotation_manual_na$group)] <- "NA"
  
  out_3 <- pheatmap_with_anno_sa(tab_response = response_metrics, 
                                 annotation_row = annotation_manual_row,
                                 annotation_col = annotation_manual_col)
  expect_length(out_3, 4)
  expect_equal(sort(names(out_3)), sort(c("matrix", "annotation_col", "annotation_row", "heatmap")))
  anno_col_3 <- out_3[["annotation_col"]]
  expect_is(anno_col_3, "data.table")
  expect_equal(anno_col_3, annotation_manual_col)
  anno_row_3 <- out_3[["annotation_row"]]
  expect_is(anno_row_3, "data.table")
  expect_equal(anno_row_3, annotation_manual_na)
  plt_3 <- out_3[["heatmap"]]
  expect_is(plt_3, "pheatmap")
  expect_equal(plt_3$gtable$grobs[[5]]$label, c("mut_A", "mut_B", "mut_C"))
  expect_equal(plt_3$gtable$grobs[[7]]$label, c("group"))
  
  annotation_map <- list(
    mut_A = c("1" = "coral", "0" = "cadetblue"),
    mut_B = c("yes" = "black", "no" = "grey90"),
    mut_C = c("AA" = "yellow", "BB" = "green")
  )
  annotation_manual_col_na <- data.table::copy(annotation_manual_col)
  annotation_manual_col_na[2, c("mut_A", "mut_B", "mut_C") := "NA"]
  
  out_4 <- pheatmap_with_anno_sa(tab_response = response_metrics, 
                                 annotation_col = annotation_manual_col[1, ],
                                 annotation_colors = annotation_map)
  expect_length(out_4, 3)
  expect_equal(sort(names(out_4)), sort(c("matrix", "annotation_col", "heatmap")))
  anno_4 <- out_4[["annotation_col"]]
  expect_is(anno_4, "data.table")
  expect_equal(anno_4, annotation_manual_col_na)
  
  annotation_manual_row <-
    unique(response_metrics[, .SD, .SDcols = c("DrugName", "drug_moa")])
  anno_test <- data.table::data.table(
    DrugName = c("drug_004", "drug_005", "drug_006"),
    tested_AB = c("yes", "yes", "no")
  )
  annotation_manual_row <- anno_test[annotation_manual_row, on = "DrugName"]
  
  out_5 <- pheatmap_with_anno_sa(tab_response = response_metrics, 
                                 annotation_row = annotation_manual_row)
  expect_length(out_5, 3)
  expect_equal(sort(names(out_5)), sort(c("matrix", "annotation_row", "heatmap")))
  data_5 <- out_5[["matrix"]]
  expect_is(data_5, "data.table")
  expect_equal(data_5, res_1)
  anno_5 <- out_5[["annotation_row"]]
  expect_is(anno_5, "data.table")
  expect_equal(anno_5, annotation_manual_row)
  plt_5 <- out_5[["heatmap"]]
  expect_is(plt_5, "pheatmap")
  expect_equal(plt_5$gtable$grobs[[5]]$label, c("tested_AB", "drug_moa"))
  
  expect_error(pheatmap_with_anno_sa(tab_response = unlist(response_metrics)),
               "Assertion on 'tab_response' failed: Must be a data.table")
  expect_error(pheatmap_with_anno_sa(tab_response = response_metrics,
                                     normalization_type = "XX"),
               "Assertion on 'normalization_type' failed: Must be element of set")
  expect_error(pheatmap_with_anno_sa(tab_response = response_metrics,
                                     metric = "xxx"),
               "Assertion on 'metric' failed: Must be element of set")
  expect_error(pheatmap_with_anno_sa(tab_response = response_metrics,
                                     fit_source = 1),
               "Assertion on 'fit_source' failed: Must be of type 'string'")
  expect_error(pheatmap_with_anno_sa(tab_response = response_metrics,
                                     hm_title = NULL),
               "Assertion on 'hm_title' failed: Must be of type 'string'")
  expect_error(pheatmap_with_anno_sa(tab_response = response_metrics,
                                     colors_vec = 1:3),
               "Assertion on 'colors_vec' failed: Must be of type 'character'")
  expect_error(pheatmap_with_anno_sa(tab_response = response_metrics,
                                     colors_vec = c("pinky", "blackish")),
               "Must be a valid color name")
  expect_error(pheatmap_with_anno_sa(tab_response = response_metrics,
                                     no_breaks = "str"),
               "Assertion on 'no_breaks' failed: Must be of type 'single integerish value'")
  expect_error(pheatmap_with_anno_sa(tab_response = response_metrics,
                                     annotation_row = unlist(annotation_manual_row)),
               "Assertion on 'annotation_row' failed: Must be a data.table")
  expect_error(pheatmap_with_anno_sa(tab_response = response_metrics,
                                     annotation_col = unlist(annotation_manual_col)),
               "Assertion on 'annotation_col' failed: Must be a data.table")
  expect_error(pheatmap_with_anno_sa(tab_response = response_metrics,
                                     annotation_colors = unlist(annotation_map)),
               "Assertion on 'annotation_colors' failed: Must be of type 'list'")
})

test_that("pheatmap_with_anno_combo works as expected", {
  mae <- gDRutils::get_synthetic_data("combo_matrix_small")
  se <- mae[[gDRutils::get_supported_experiments("combo")]]
  response_metrics <- gDRutils::convert_se_assay_to_dt(se = se, assay_name = "scores")
  
  cdata <- SummarizedExperiment::colData(se)
  rdata <- SummarizedExperiment::rowData(se)
  
  res_1 <- data.table::dcast(data = response_metrics[normalization_type == "GR", ],
                             formula = CellLineName ~ paste(DrugName, "x", DrugName_2), 
                             value.var = "hsa_score")
  data.table::setkey(res_1, NULL)
  
  out_1 <- pheatmap_with_anno_combo(tab_response = response_metrics)
  expect_length(out_1, 2)
  expect_equal(names(out_1), c("matrix", "heatmap"))
  data_1 <- out_1[["matrix"]]
  expect_is(data_1, "data.table")
  expect_equal(data_1, res_1)
  plt_1 <- out_1[["heatmap"]]
  expect_is(plt_1, "pheatmap")
  expect_equal(plt_1$gtable$grobs[[2]]$label, cdata[["CellLineName"]])
  expect_equal(plt_1$gtable$grobs[[3]]$label, sprintf("%s x %s", rdata$DrugName, rdata$DrugName_2))
  expect_true(all(vapply(plt_1$gtable$grobs[[1]]$children[[1]]$gp$fill, is_valid_color, logical(1))))
  
  response_metrics_na <- data.table::copy(response_metrics)
  response_metrics_na[DrugName == "drug_004"]$bliss_score <- NA
  res_2 <- data.table::dcast(data = response_metrics_na[normalization_type == "RV", ],
                             formula = CellLineName ~ paste(DrugName, "x", DrugName_2), 
                             value.var = "bliss_score")
  res_2 <- res_2[, .SD, .SDcols = !anyNA]
  data.table::setkey(res_2, NULL)
  drug_combo <- unique(sprintf("%s x %s",
                               response_metrics_na[!is.na(bliss_score)]$DrugName, 
                               response_metrics_na[!is.na(bliss_score)]$DrugName_2))
  
  out_2 <- pheatmap_with_anno_combo(tab_response = response_metrics_na, 
                                    normalization_type = "RV",
                                    metric = "bliss_score")
  expect_length(out_2, 2)
  expect_equal(names(out_2), c("matrix", "heatmap"))
  data_2 <- out_2[["matrix"]]
  expect_is(data_2, "data.table")
  expect_equal(data_2, res_2)
  plt_2 <- out_2[["heatmap"]]
  expect_is(plt_2, "pheatmap")
  expect_equal(plt_2$gtable$grobs[[3]]$label, drug_combo)
  
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
  
  out_3 <- pheatmap_with_anno_combo(tab_response = response_metrics, 
                                    annotation_col = annotation_manual_col,
                                    annotation_colors = annotation_map)
  expect_length(out_3, 3)
  expect_equal(sort(names(out_3)), sort(c("matrix", "annotation_col", "heatmap")))
  expect_equal(out_3[["annotation_col"]], annotation_manual_col)
  plt_3 <- out_3[["heatmap"]]
  expect_is(plt_3, "pheatmap")
  expect_equal(plt_3$gtable$grobs[[5]]$label, c("mut_A", "mut_B"))
  
  annotation_manual_row <-
    unique(response_metrics[, .SD, .SDcols = c("DrugName", "DrugName_2", "drug_moa", "drug_moa_2")])
  
  out_5 <- pheatmap_with_anno_combo(tab_response = response_metrics, 
                                    annotation_row = annotation_manual_row)
  expect_length(out_5, 3)
  expect_equal(sort(names(out_5)), sort(c("matrix", "annotation_row", "heatmap")))
  data_5 <- out_5[["matrix"]]
  expect_is(data_5, "data.table")
  expect_equal(data_5, res_1)
  anno_5 <- out_5[["annotation_row"]]
  expect_is(anno_5, "data.table")
  expect_equal(anno_5, annotation_manual_row)
  plt_5 <- out_5[["heatmap"]]
  expect_is(plt_5, "pheatmap")
  expect_equal(plt_5$gtable$grobs[[5]]$label, c("drug_moa", "drug_moa_2"))
  
  expect_error(pheatmap_with_anno_combo(tab_response = unlist(response_metrics)),
               "Assertion on 'tab_response' failed: Must be a data.table")
  expect_error(pheatmap_with_anno_combo(tab_response = response_metrics,
                                        normalization_type = "XX"),
               "Assertion on 'normalization_type' failed: Must be element of set")
  expect_error(pheatmap_with_anno_combo(tab_response = response_metrics,
                                        metric = "xxx"),
               "Assertion on 'metric' failed: Must be element of set")
  expect_error(pheatmap_with_anno_combo(tab_response = response_metrics,
                                        fit_source = 1),
               "Assertion on 'fit_source' failed: Must be of type 'string'")
  expect_error(pheatmap_with_anno_combo(tab_response = response_metrics,
                                        hm_title = NULL),
               "Assertion on 'hm_title' failed: Must be of type 'string'")
  expect_error(pheatmap_with_anno_combo(tab_response = response_metrics,
                                        colors_vec = 1:3),
               "Assertion on 'colors_vec' failed: Must be of type 'character'")
  expect_error(pheatmap_with_anno_combo(tab_response = response_metrics,
                                        colors_vec = c("pinky", "blackish")),
               "Must be a valid color name")
  expect_error(pheatmap_with_anno_combo(tab_response = response_metrics,
                                        no_breaks = "str"),
               "Assertion on 'no_breaks' failed: Must be of type 'single integerish value'")
  # expect_error(pheatmap_with_anno_combo(tab_response = response_metrics,
  #                                       annotation_row = unlist(annotation_manual_row)),
  #              "Assertion on 'annotation_row' failed: Must be a data.table")
  expect_error(pheatmap_with_anno_combo(tab_response = response_metrics,
                                        annotation_col = unlist(annotation_manual_col)),
               "Assertion on 'annotation_col' failed: Must be a data.table")
  expect_error(pheatmap_with_anno_combo(tab_response = response_metrics,
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
  response_metrics <- gDRutils::convert_se_assay_to_dt(se = se, assay_name = "Averaged")
  dt_ann <- unique(response_metrics[, .SD, .SDcols = c("Tissue", "ReferenceDivisionTime")])
  
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