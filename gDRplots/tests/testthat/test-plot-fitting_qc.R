context("Test fitting_qc")

test_that("heatmap_control_mapping_qc works as expected", {
  mae <- gDRutils::get_synthetic_data("combo_matrix")
  se <- mae[[gDRutils::get_supported_experiments("sa")]]
  cdata <- SummarizedExperiment::colData(se)
  rdata <- SummarizedExperiment::rowData(se)
  
  treat <- gDRutils::convert_se_assay_to_dt(se, "RawTreated")
  controls <- gDRutils::convert_se_assay_to_dt(se, "Controls")
  
  plt_1 <- heatmap_control_mapping_qc(dt_treat = treat,
                                      dt_controls = controls)
  expect_is(plt_1, "pheatmap")
  expect_equal(plt_1$gtable$grobs[[1]]$label, "Counts of mapped controls") # title
  expect_equal(plt_1$gtable$grobs[[3]]$label, cdata[["CellLineName"]])
  expect_equal(plt_1$gtable$grobs[[4]]$label, rdata[["DrugName"]])
  expect_true(all(vapply(plt_1$gtable$grobs[[2]]$children[[1]]$gp$fill, is_valid_color, logical(1))))
  expect_length(unique(as.vector(plt_1$gtable$grobs[[2]]$children[[1]]$gp$fill)), 1)
  
  # lack of controls
  controls_2 <- controls[DrugName %in% c("drug_002", "drug_011")]
  plt_2 <- heatmap_control_mapping_qc(dt_treat = treat,
                                      dt_controls = controls_2)
  expect_length(unique(as.vector(plt_2$gtable$grobs[[2]]$children[[1]]$gp$fill)), 2)
  
  # random lack of controls
  controls_3 <- controls[CorrectedReadout %in% 98:99]
  frequency <- controls_3[, .N, by = .(rId, cId)]
  plt_3 <- heatmap_control_mapping_qc(dt_treat = treat,
                                      dt_controls = controls_3)
  expect_length(unique(as.vector(plt_3$gtable$grobs[[2]]$children[[1]]$gp$fill)), 
                NROW(unique(frequency$N)) + 1) # number of unique val + "red" for NA
  
  se <- mae[[gDRutils::get_supported_experiments("combo")]]
  cdata <- SummarizedExperiment::colData(se)
  rdata <- SummarizedExperiment::rowData(se)
  
  treat <- gDRutils::convert_se_assay_to_dt(se, "RawTreated")
  controls <- gDRutils::convert_se_assay_to_dt(se, "Controls")
  
  plt_4 <- heatmap_control_mapping_qc(dt_treat = treat,
                                      dt_controls = controls)
  expect_is(plt_4, "pheatmap")
  expect_equal(plt_4$gtable$grobs[[1]]$label, "Counts of mapped controls") # title
  expect_equal(plt_4$gtable$grobs[[3]]$label, cdata[["CellLineName"]])
  expect_equal(plt_4$gtable$grobs[[4]]$label,
               paste(rdata[["DrugName"]], rdata[["DrugName_2"]], sep = " x "))
  expect_true(all(vapply(plt_4$gtable$grobs[[2]]$children[[1]]$gp$fill, is_valid_color, logical(1))))
  expect_length(NROW(unique(plt_4$gtable$grobs[[2]]$children[[1]]$gp$fill)), 1)
  
  
  expect_error(heatmap_control_mapping_qc(dt_treat = unlist(treat),
                                          dt_controls = controls),
               "Check on 'dt_treat' failed: Must be a data.table")
  expect_error(heatmap_control_mapping_qc(dt_treat = treat,
                                          dt_controls = as.list(controls)),
               "Check on 'dt_controls' failed: Must be a data.table")
})
