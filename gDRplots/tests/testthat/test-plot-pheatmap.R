context("Test qc_heatmap")

test_that("pheatmap_with_anno_sa works as expected", {
  mae <- gDRutils::get_synthetic_data("combo_matrix_small")
  se <- mae[[gDRutils::get_supported_experiments("sa")]]
  response_metrics <- gDRutils::convert_se_assay_to_dt(se = se, assay_name = "Metrics")
  
  cdata <- SummarizedExperiment::colData(se)
  rdata <- SummarizedExperiment::rowData(se)
  
  plt_1 <- pheatmap_with_anno_sa(tab_response = response_metrics)
  expect_is(plt_1, "pheatmap")
  expect_equal(plt_1$gtable$grobs[[2]]$label, cdata[["CellLineName"]])
  expect_equal(plt_1$gtable$grobs[[3]]$label, rdata[["DrugName"]])
  expect_true(all(vapply(plt_1$gtable$grobs[[1]]$children[[1]]$gp$fill, is_valid_color, logical(1))))
  
  response_metrics_na <- data.table::copy(response_metrics)
  response_metrics_na[DrugName %in% c("drug_021", "drug_026")]$x_max <- NA
  
  plt_2 <- pheatmap_with_anno_sa(tab_response = response_metrics_na, metric = "x_max")
  expect_is(plt_2, "pheatmap")
  expect_equal(plt_2$gtable$grobs[[3]]$label, unique(response_metrics_na[!is.na(x_max)]$DrugName))
  
  annotation_manual <- data.table::data.table(
    CellLineName = c("cellline_GB", "cellline_HB"),
    mut_A = c(1, 0),
    mut_B = c("yes", "no")
  )
  plt_3 <- pheatmap_with_anno_sa(tab_response = response_metrics, annotation_col = annotation_manual)
  expect_is(plt_3, "pheatmap")
  expect_equal(plt_3$gtable$grobs[[5]]$label, c("mut_A", "mut_B"))
})

test_that("pheatmap_with_anno_combo works as expected", {
  mae <- gDRutils::get_synthetic_data("combo_matrix_small")
  se <- mae[[gDRutils::get_supported_experiments("combo")]]
  response_metrics <- gDRutils::convert_se_assay_to_dt(se = se, assay_name = "scores")
  
  cdata <- SummarizedExperiment::colData(se)
  rdata <- SummarizedExperiment::rowData(se)
  
  plt_1 <- pheatmap_with_anno_combo(tab_response = response_metrics)
  expect_is(plt_1, "pheatmap")
  expect_equal(plt_1$gtable$grobs[[2]]$label, cdata[["CellLineName"]])
  expect_equal(plt_1$gtable$grobs[[3]]$label, sprintf("%s x %s", rdata$DrugName, rdata$DrugName_2))
  expect_true(all(vapply(plt_1$gtable$grobs[[1]]$children[[1]]$gp$fill, is_valid_color, logical(1))))
  
  response_metrics_na <- data.table::copy(response_metrics)
  response_metrics_na[DrugName == "drug_004"]$bliss_score <- NA
  drug_combo <- unique(sprintf("%s x %s", 
                               response_metrics_na[!is.na(bliss_score)]$DrugName, 
                               response_metrics_na[!is.na(bliss_score)]$DrugName_2))
  
  plt_2 <- pheatmap_with_anno_combo(tab_response = response_metrics_na, metric = "bliss_score")
  expect_is(plt_2, "pheatmap")
  expect_equal(plt_2$gtable$grobs[[3]]$label, drug_combo)
  
  annotation_manual <- data.table::data.table(
    CellLineName = c("cellline_GB", "cellline_HB"),
    mut_A = c(1, 0),
    mut_B = c("yes", "no")
  )
  plt_3 <- pheatmap_with_anno_combo(tab_response = response_metrics, annotation_col = annotation_manual)
  expect_is(plt_3, "pheatmap")
  expect_equal(plt_3$gtable$grobs[[5]]$label, c("mut_A", "mut_B"))
})

test_that("get_hm_title works as expected", {
  metric <- "xc50"
  metric_growth <- "GR"
  expect_equal(get_hm_title(metric, metric_growth), "log10(GR50)")
  
  metric <- "bliss_score"
  metric_growth <- "RV"
  expect_equal(get_hm_title(metric, metric_growth), "Bliss Score RV")
  
  dataset_name <- "Dataset AB123"
  metric <- "xc50"
  metric_growth <- "RV"
  expect_equal(get_hm_title(metric, metric_growth, dataset_name), "Dataset AB123 (log10(IC50))")
  
  metric <- "hsa_score"
  metric_growth <- "GR"
  expect_equal(get_hm_title(metric, metric_growth, dataset_name), "Dataset AB123 (HSA Score GR)")
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
})
