context("Test sa_plots")

test_that("plot_dose_response_sa works as expected", {
  mae <- gDRutils::get_synthetic_data("small")
  se <- mae[[1]]
  
  grouping <- gDRutils::get_env_identifiers("cellline_name")
  iR <- rownames(se)[1]
  dt_metrics <- gDRutils::convert_se_assay_to_dt(se[iR], "Metrics")
  dt_average <- gDRutils::convert_se_assay_to_dt(se[iR], "Averaged")
  
  plt_1 <- plot_dose_response_sa(dt_metrics = dt_metrics,
                                 dt_average = dt_average,
                                 grouping = grouping)
  expect_is(plt_1, "gg")
  expect_true(grepl("GR", plt_1[["labels"]][["y"]]))
  expect_length(plt_1[["layers"]], 3)
  

  iC <- colnames(se)[1:2]
  dt_metrics <- gDRutils::convert_se_assay_to_dt(se[iR, iC], "Metrics")
  dt_average <- gDRutils::convert_se_assay_to_dt(se[iR, iC], "Averaged")
  # lack of data for curve and small number of conc point
  dt_metrics[, c("ec50", "x_inf", "x_0", "h", "r2")] <- NA
  dt_average <- dt_average[Concentration %in% unique(dt_average$Concentration)[1:2], ]
  
  plt_2 <- plot_dose_response_sa(dt_metrics = dt_metrics,
                                 dt_average = dt_average,
                                 grouping = grouping)
  expect_is(plt_2, "gg")
  expect_true(grepl("GR", plt_2[["labels"]][["y"]]))
  expect_length(plt_2[["layers"]], 3)
  expect_equal(names(plt_2[["guides"]][["guides"]]), "colour")

  grouping <- gDRutils::get_env_identifiers("drug_name")
  iC <- colnames(se)[1]
  dt_metrics <- gDRutils::convert_se_assay_to_dt(se[, iC], "Metrics")
  dt_average <- gDRutils::convert_se_assay_to_dt(se[, iC], "Averaged")
  normalization_type <- "RV"
  
  plt_3 <- plot_dose_response_sa(dt_metrics = dt_metrics,
                                 dt_average = dt_average,
                                 grouping = grouping,
                                 normalization_type = normalization_type,
                                 colormap = c("cadetblue", "orange", "darkblue"),
                                 plot_fit_flag = FALSE)
  expect_is(plt_3, "gg")
  expect_true(grepl(normalization_type, plt_3[["labels"]][["y"]]))
  expect_length(plt_3[["layers"]], 2)
  
  expect_error(plot_dose_response_sa(dt_metrics = as.list(dt_metrics),
                                     dt_average = dt_average,
                                     grouping = grouping),
               "Assertion on 'dt_metrics' failed: Must be a data.table")
  expect_error(plot_dose_response_sa(dt_metrics = dt_metrics,
                                     dt_average = dt_average,
                                     grouping = "str"),
               "Assertion on 'grouping' failed: Must be element of set")
})

test_that("plot_dose_response_sa_by_CLs works as expected", {
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  clid <- gDRutils::get_env_identifiers("cellline")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  gnumber <- gDRutils::get_env_identifiers("drug")
  
  mae <- gDRutils::get_synthetic_data("small")
  se <- mae[[1]]
  dt_metrics <- gDRutils::convert_se_assay_to_dt(se, "Metrics")
  dt_average <- gDRutils::convert_se_assay_to_dt(se, "Averaged")
  
  cellline_name_vec <- unique(dt_metrics[[cellline_name]])[2:5]
  drug_name_vec <- unique(dt_metrics[[drug_name]])[5:7]
  
  plts <- plot_dose_response_sa_by_CLs(dt_metrics = dt_metrics, 
                                       dt_average = dt_average)
  expect_is(plts, "list")
  expect_equal(names(plts), 
               sprintf("%s (%s)", unique(dt_metrics[[drug_name]]), unique(dt_metrics[[gnumber]])))
  
  normalization_type <- "RV"
  plts <- plot_dose_response_sa_by_CLs(dt_metrics = dt_metrics, 
                                       dt_average = dt_average,
                                       cellline_name_vec = cellline_name_vec,
                                       drug_name_vec = drug_name_vec,
                                       normalization_type = normalization_type,
                                       colormap = c("#00008B", "#FF6347", "#4CBB17"))
  expect_is(plts, "list")
  plotted_ <- intersect(drug_name_vec, unique(dt_metrics[[drug_name]]))
  expect_equal(names(plts), 
               sprintf("%s (%s)", plotted_, unique(dt_metrics[get(drug_name) %in% plotted_][[gnumber]])))
  expect_true(all(vapply(seq_along(plts), 
                         function(i) grepl(normalization_type, plts[[i]]$labels$y), logical(1))))
  
  cellline_name_vec_2 <- c(cellline_name_vec, "cellline_XX")
  drug_name_vec_2 <- c(drug_name_vec, "drug_100")
  
  plts <- plot_dose_response_sa_by_CLs(dt_metrics = dt_metrics, 
                                       dt_average = dt_average,
                                       cellline_name_vec = cellline_name_vec_2,
                                       drug_name_vec = drug_name_vec_2)
  expect_is(plts, "list")
  plotted_ <- intersect(drug_name_vec_2, unique(dt_metrics[[drug_name]]))
  expect_equal(names(plts), 
               sprintf("%s (%s)", plotted_, unique(dt_metrics[get(drug_name) %in% plotted_][[gnumber]])))
  plotted <- intersect(cellline_name_vec_2, unique(dt_metrics[[cellline_name]]))
  expect_true(all(vapply(seq_along(plts), 
                         function(i) all(plts[[i]]$plot_env$group_names == plotted), logical(1))))
})

test_that("plot_dose_response_sa_qc works as expected", {
  mae <- gDRutils::get_synthetic_data("small")
  se <- mae[[gDRutils::get_supported_experiments("sa")]]
  
  dt_metrics <- gDRutils::convert_se_assay_to_dt(se, "Metrics")
  dt_average <- gDRutils::convert_se_assay_to_dt(se, "Averaged")
  cl_name <- dt_metrics[["CellLineName"]][1]
  d_name <- dt_metrics[["DrugName"]][1]
  
  plt_1 <- plot_dose_response_sa_qc(dt_metrics = dt_metrics,
                                    dt_average = dt_average,
                                    cl_name = cl_name,
                                    d_name = d_name)
  expect_is(plt_1, "gg")
  expect_true(grepl("GR", plt_1[["labels"]][["y"]]))
  expect_true(grepl(d_name, plt_1[["labels"]][["title"]]))
  expect_length(plt_1[["layers"]], 4)
  
  normalization_type <- "RV"
  plt_2 <- plot_dose_response_sa_qc(dt_metrics = dt_metrics,
                                    dt_average = dt_average,
                                    cl_name = cl_name,
                                    d_name = d_name, 
                                    normalization_type = normalization_type)
  expect_is(plt_2, "gg")
  expect_true(grepl(normalization_type, plt_2[["labels"]][["y"]]))
})

test_that("plot_dose_response_sa_qc works as expected", {
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  
  mae <- gDRutils::get_synthetic_data("small")
  se <- mae[[1]]
  
  dt_metrics <- gDRutils::convert_se_assay_to_dt(se, "Metrics")
  dt_average <- gDRutils::convert_se_assay_to_dt(se, "Averaged")
  cl_name <- dt_metrics[[cellline_name]][1]
  d_names <- unique(dt_metrics[[drug_name]])[1:3]
  
  plt_1 <- plot_dose_response_sa_qc_panel(dt_metrics = dt_metrics,
                                          dt_average = dt_average,
                                          cl_name = cl_name)
  expect_is(plt_1, "gg")
  expect_length(plt_1[["layers"]], 4)
  expect_true(grepl(cl_name, plt_1[["labels"]][["title"]]))
  expect_true(grepl("GR", plt_1[["labels"]][["y"]]))
  
  normalization_type <- "RV"
  plt_2 <- plot_dose_response_sa_qc_panel(dt_metrics = dt_metrics,
                                          dt_average = dt_average,
                                          cl_name = cl_name,
                                          d_names = d_names,
                                          normalization_type = normalization_type)
  expect_is(plt_2, "gg")
  expect_true(grepl(normalization_type, plt_2[["labels"]][["y"]]))
})
