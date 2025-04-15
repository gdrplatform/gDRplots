context("Test var_distribution")

test_that("plot_var_distribution_qc works as expected", {
  mae <- gDRutils::get_synthetic_data("combo_matrix")
  se_sa <- mae[[gDRutils::get_supported_experiments("sa")]]
  se_combo <- mae[[gDRutils::get_supported_experiments("combo")]]
  
  dt_norm <- gDRutils::convert_se_assay_to_dt(se_sa, "Normalized")
  cl_name <- dt_norm[["CellLineName"]][1]
  
  plt_1 <- plot_var_distribution_qc(dt_assay = dt_norm,
                                    cl_name = cl_name)
  expect_is(plt_1, "gg")
  expect_true(grepl("GR", plt_1[["labels"]][["y"]]))
  expect_length(plt_1[["layers"]], 4)
  
  dt_average <- gDRutils::convert_se_assay_to_dt(se_combo, "Averaged")
  cl_name <- dt_norm[["CellLineName"]][1]
  
  plt_2 <- plot_var_distribution_qc(dt_assay = dt_average,
                                    cl_name = cl_name,
                                    metric = "x_std",
                                    normalization_type = "RV")
  expect_is(plt_2, "gg")
  expect_true(grepl("RV", plt_2[["labels"]][["y"]]))
  expect_true(grepl("x_std", plt_2[["labels"]][["y"]]))
  expect_length(plt_2[["layers"]], 3)
  
  dt_norm_small <- 
    rbind(dt_norm[DrugName %in% c("drug_001", "drug_002", "drug_011")], 
          dt_norm[DrugName == "drug_021" & normalization_type == "GR", ][1, ])
  
  plt_3 <- plot_var_distribution_qc(dt_assay = dt_norm_small,
                                    cl_name = cl_name)
  expect_is(plt_3, "gg")
  expect_true(grepl("GR", plt_3[["labels"]][["y"]]))
  expect_length(plt_3[["plot_env"]][["color_palette"]], NROW(unique(dt_norm_small[["DrugName"]])))
})
