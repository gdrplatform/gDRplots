context("Test fitting_qc")

test_that("plot_var_stat_qc works as expected", {
  mae <- gDRutils::get_synthetic_data("small")
  se <- mae[[gDRutils::get_supported_experiments("sa")]]
  
  dt_metrics <- gDRutils::convert_se_assay_to_dt(se, "Metrics")
  cl_name <- dt_metrics[["CellLineName"]][1]
  
  plt_1 <- plot_var_stat_qc(dt_assay = dt_metrics,
                            cl_name = cl_name)
  expect_is(plt_1, "gg")
  expect_length(plt_1[["layers"]], 4)
  expect_true(grepl(cl_name, plt_1[["labels"]][["title"]]))
  expect_true(grepl("r2", plt_1[["labels"]][["y"]]))
  expect_true(grepl("GR", plt_1[["labels"]][["y"]]))
  expect_equal(NROW(ggplot2::ggplot_build(plt_1)[["data"]][[3]]),
               NROW(unique(dt_metrics[["DrugName"]])))
  
  plt_2 <- plot_var_stat_qc(dt_assay = dt_metrics,
                            cl_name = cl_name,
                            normalization_type = "RV")
  expect_is(plt_2, "gg")
  expect_true(grepl("RV", plt_2[["labels"]][["y"]]))
  
  plt_3 <- plot_var_stat_qc(dt_assay = dt_metrics,
                            cl_name = cl_name,
                            metric = "x_AOC",
                            normalization_type = "RV",
                            with_table = TRUE)
  expect_is(plt_3, "gg")
  expect_length(plt_3[["layers"]], 2)
  
  expect_error(plot_var_stat_qc(dt_assay = unlist(dt_metrics),
                                cl_name = cl_name),
               "Assertion on 'dt_assay' failed: Must be a data.table")
  expect_error(plot_var_stat_qc(dt_assay = dt_metrics,
                                cl_name = c("cl_1", "cl_2")),
               "Assertion on 'cl_name' failed: Must have length 1.")
  expect_error(plot_var_stat_qc(dt_assay = dt_metrics,
                                cl_name = "unknown_cl"),
               "Assertion on 'cl_name' failed: Must be element of set")
  expect_error(plot_var_stat_qc(dt_assay = dt_metrics,
                                cl_name = cl_name, 
                                metric = "strange_metric"),
               "Assertion on 'metric' failed: Must be element of set")
  expect_error(plot_var_stat_qc(dt_assay = dt_metrics,
                                cl_name = cl_name, 
                                normalization_type = "XX"),
               "Assertion on 'normalization_type' failed: Must be element of set")
  expect_error(plot_var_stat_qc(dt_assay = dt_metrics,
                                cl_name = cl_name, 
                                with_table = "yes"),
               "Assertion on 'with_table' failed: Must be of type 'logical flag'")
})

test_that("plot_fitting_acc works as expected", {
  mae <- gDRutils::get_synthetic_data("small")
  se <- mae[[gDRutils::get_supported_experiments("sa")]]
  
  dt_metrics <- gDRutils::convert_se_assay_to_dt(se, "Metrics")
  cl_name <- dt_metrics[["CellLineName"]][1]
  
  plt_1 <- plot_fitting_acc(dt_assay = dt_metrics,
                            cl_name = cl_name)
  expect_is(plt_1, "gg")
  
  plt_2 <- plot_fitting_acc(dt_assay = dt_metrics,
                            cl_name = cl_name,
                            normalization_type = "RV")
  expect_is(plt_2, "gg")
  
  dt_metric_NA <- data.table::copy(dt_metrics)
  dt_metric_NA[CellLineName == cl_name]$r2 <- NA
  
  expect_warning({
    plt_3 <- plot_fitting_acc(dt_assay = dt_metric_NA,
                              cl_name = cl_name)
    } , sprintf("Missing data for %s in GR normalization type.", cl_name))
  expect_length(ggplot2::ggplot_build(plt_3)[["data"]][[1]], 0)
  
  expect_error(plot_fitting_acc(dt_assay = unlist(dt_metrics),
                                cl_name = cl_name),
               "Assertion on 'dt_assay' failed: Must be a data.table")
  expect_error(plot_fitting_acc(dt_assay = dt_metrics,
                                cl_name = 123),
               "Assertion on 'cl_name' failed: Must be of type 'string'")
  expect_error(plot_fitting_acc(dt_assay = dt_metrics,
                                cl_name = "unknown_cl"),
               "Assertion on 'cl_name' failed: Must be element of set")
  expect_error(plot_fitting_acc(dt_assay = dt_metrics,
                                cl_name = cl_name,
                                normalization_type = "XX"),
               "Assertion on 'normalization_type' failed: Must be element of set")
})

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
               "Assertion on 'dt_treat' failed: Must be a data.table")
  expect_error(heatmap_control_mapping_qc(dt_treat = treat,
                                          dt_controls = as.list(controls)),
               "Assertion on 'dt_controls' failed: Must be a data.table")
})
