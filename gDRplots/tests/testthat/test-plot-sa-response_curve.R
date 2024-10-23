context("Test sa_plots")

test_that("plot_dose_response_sa works as expected", {
  mae <- gDRutils::get_synthetic_data("small")
  se <- mae[[1]]
  
  group_var <- gDRutils::get_env_identifiers("cellline_name")
  selected_drug <- "drug_002"
  dt_metrics <- gDRutils::convert_se_assay_to_dt(se, "Metrics")
  dt_average <- gDRutils::convert_se_assay_to_dt(se, "Averaged")
  
  plt_1 <- plot_dose_response_sa(dt_metrics = dt_metrics,
                                 dt_average = dt_average,
                                 selection_name = selected_drug,
                                 group_var = group_var)
  expect_is(plt_1, "gg")
  expect_equal(plt_1[["labels"]][["y"]], "GR")
  expect_length(plt_1[["layers"]], 3)
  
  plt_2 <- plot_dose_response_sa(dt_metrics = dt_metrics,
                                 dt_average = dt_average,
                                 group_var = group_var,
                                 selection_name = selected_drug,
                                 colors_vec = rainbow(NROW(unique(dt_metrics[[group_var]]))))
  expect_is(plt_2, "gg")
  expect_length(unique(ggplot2::ggplot_build(plt_2)$data[[2]][["colour"]]),
                NROW(unique(dt_metrics[[group_var]])))
  
  subset_celline <- c("cellline_BA", "cellline_CA")
  # lack of data for curve and small number of conc point
  dt_metrics_lack <- data.table::copy(dt_metrics)[CellLineName %in% subset_celline, ]
  dt_metrics_lack[, c("ec50", "x_inf", "x_0", "h", "r2")] <- NA
  dt_average_lack <- 
    data.table::copy(dt_average)[Concentration %in% unique(dt_average$Concentration)[1:2], ]
  
  plt_3 <- plot_dose_response_sa(dt_metrics = dt_metrics_lack,
                                 dt_average = dt_average_lack,
                                 selection_name = selected_drug,
                                 group_var = group_var)
  expect_is(plt_3, "gg")
  expect_equal(plt_3[["labels"]][["y"]], "GR")
  expect_length(plt_3[["layers"]], 3)
  expect_equal(names(plt_3[["guides"]][["guides"]]), "colour")
  
  group_var <- gDRutils::get_env_identifiers("drug_name")
  selected_celline <- "cellline_BA"
  normalization_type <- "RV"
  
  plt_4 <- plot_dose_response_sa(dt_metrics = dt_metrics,
                                 dt_average = dt_average,
                                 selection_name = selected_celline,
                                 group_var = group_var,
                                 normalization_type = normalization_type,
                                 colors_vec = c("cadetblue", "orange", "darkblue"),
                                 plot_fit_flag = FALSE)
  expect_is(plt_4, "gg")
  expect_equal(plt_4[["labels"]][["y"]], normalization_type)
  expect_true(grepl(selected_celline, plt_4[["labels"]][["title"]]))
  expect_length(plt_4[["layers"]], 2)
  
  expect_error(plot_dose_response_sa(dt_metrics = as.list(dt_metrics),
                                     dt_average = dt_average,
                                     selection_name = selected_celline,
                                     group_var = group_var),
               "Assertion on 'dt_metrics' failed: Must be a data.table")
  expect_error(plot_dose_response_sa(dt_metrics = dt_metrics,
                                     dt_average = dt_average,
                                     selection_name = c("cl_1", "cl_2")),
               "Assertion on 'selection_name' failed: Must have length 1.")
  expect_error(plot_dose_response_sa(dt_metrics = dt_metrics,
                                     dt_average = dt_average,
                                     selection_name = selected_celline,
                                     group_var = "str"),
               "Assertion on 'group_var' failed: Must be element of set")
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
                                       colors_vec = c("#00008B", "#FF6347", "#4CBB17"))
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
  cl_name <- "cellline_BA"
  d_name <- "drug_002"
  
  plt_1 <- plot_dose_response_sa_qc(dt_metrics = dt_metrics,
                                    dt_average = dt_average,
                                    cl_name = cl_name,
                                    d_name = d_name)
  expect_is(plt_1, "gg")
  expect_equal(plt_1[["labels"]][["y"]], "GR")
  expect_true(grepl(cl_name, plt_1[["labels"]][["title"]]))
  expect_true(grepl(d_name, plt_1[["labels"]][["colour"]]))
  expect_length(plt_1[["layers"]], 5)
  
  normalization_type <- "RV"
  plt_2 <- plot_dose_response_sa_qc(dt_metrics = dt_metrics,
                                    dt_average = dt_average,
                                    cl_name = cl_name,
                                    d_name = d_name, 
                                    normalization_type = normalization_type)
  expect_is(plt_2, "gg")
  expect_equal(plt_2[["labels"]][["y"]], normalization_type)
  
  # lack of data in dt_metrics
  dt_metrics_na <- data.table::copy(dt_metrics)
  dt_metrics_na$x_inf <- NA
  plt_3 <- plot_dose_response_sa_qc(dt_metrics = dt_metrics_na,
                                    dt_average = dt_average,
                                    cl_name = cl_name,
                                    d_name = d_name)
  expect_is(plt_3, "gg")
  expect_length(ggplot2::ggplot_build(plt_3)$data[[2]], 0) # no "Fitted Curve"
  
  # one concentration
  plt_4 <- plot_dose_response_sa_qc(dt_metrics = dt_metrics,
                                    dt_average = dt_average[Concentration == 10, ],
                                    cl_name = cl_name,
                                    d_name = d_name)
  expect_is(plt_4, "gg")
  expect_length(ggplot2::ggplot_build(plt_4)$data[[2]], 0) # no "Fitted Curve"
})

test_that("plot_dose_response_sa_qc works as expected", {
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  
  mae <- gDRutils::get_synthetic_data("small")
  se <- mae[[1]]
  
  dt_metrics <- gDRutils::convert_se_assay_to_dt(se, "Metrics")
  dt_average <- gDRutils::convert_se_assay_to_dt(se, "Averaged")
  cl_name <- "cellline_BA"
  d_names <- c("drug_002", "drug_003", "drug_004")
  
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
  
  ls_drug <- c(d_names, "drug_YY")
  no_comb_err <- NROW(unique(
    dt_average[DrugName %in% ls_drug & CellLineName == cl_name, c("CellLineName", "DrugName"), with = FALSE]
  ))
  plt_3 <- plot_dose_response_sa_qc_panel(dt_metrics = dt_metrics,
                                          dt_average = dt_average,
                                          cl_name = cl_name,
                                          d_names = ls_drug)
  expect_is(plt_3, "gg")
  expect_equal(plt_3[["labels"]][["y"]], "GR")
  expect_equal(NROW(ggplot2::ggplot_build(plt_3)$data[[1]]), 2 * no_comb_err)
  
  dt_metrics_na <- data.table::copy(dt_metrics)
  dt_metrics_na[DrugName == "drug_002"]$x_0 <- NA
  no_curve <- NROW(unique(stats::na.omit(dt_metrics_na[DrugName %in% d_names])[["DrugName"]]))
  plt_4 <- plot_dose_response_sa_qc_panel(dt_metrics = dt_metrics_na,
                                          dt_average = dt_average,
                                          cl_name = cl_name,
                                          d_names = d_names)
  expect_is(plt_4, "gg")
  expect_equal(plt_4[["labels"]][["y"]], "GR")
  expect_equal(NROW(unique(ggplot2::ggplot_build(plt_4)$data[[2]]["PANEL"])), no_curve)
  
  plt_5 <- plot_dose_response_sa_qc_panel(dt_metrics = dt_metrics,
                                          dt_average = dt_average[Concentration == 10, ],
                                          cl_name = cl_name,
                                          d_names = d_names)
  expect_is(plt_5, "gg")
  expect_length(ggplot2::ggplot_build(plt_5)$data[[2]], 0) # no "Fitted Curve"
  
  dt_average_na <- data.table::copy(dt_average)
  dt_average_na <- dt_average_na[DrugName %in% d_names[1:2] | 
                                   (DrugName == d_names[3] & Concentration == 10)]
  no_curve <- NROW(dt_average_na[, data.table::uniqueN(Concentration), by = DrugName][V1 > 1])
  plt_6 <- plot_dose_response_sa_qc_panel(dt_metrics = dt_metrics,
                                          dt_average = dt_average_na,
                                          cl_name = cl_name,
                                          d_names = d_names)
  expect_is(plt_6, "gg")
  expect_equal(NROW(unique(ggplot2::ggplot_build(plt_6)$data[[2]]["PANEL"])), no_curve)
})
