context("Test sa_plots")

test_that("plot_dose_response_sa works as expected", {
  mae <- gDRutils::get_synthetic_data("small")
  se <- mae[[1]]
  
  grouping <- gDRutils::get_env_identifiers("cellline_name")
  iR <- rownames(se)[1]
  dt_metrics <- gDRutils::convert_se_assay_to_dt(se[iR], "Metrics")
  dt_average <- gDRutils::convert_se_assay_to_dt(se[iR], "Averaged")
  
  plt <- plot_dose_response_sa(dt_metrics = dt_metrics,
                               dt_average = dt_average,
                               grouping = grouping)
  expect_is(plt, "gg")
  expect_true(grepl("GR", plt[["labels"]][["y"]]))
  expect_length(plt[["layers"]], 4)
  
  grouping <- gDRutils::get_env_identifiers("drug_name")
  iC <- colnames(se)[1]
  dt_metrics <- gDRutils::convert_se_assay_to_dt(se[, iC], "Metrics")
  dt_average <- gDRutils::convert_se_assay_to_dt(se[, iC], "Averaged")
  metric_growth <- "RV"
  
  plt <- plot_dose_response_sa(dt_metrics = dt_metrics,
                               dt_average = dt_average,
                               grouping = grouping,
                               metric_growth = metric_growth,
                               colormap = c("cadetblue", "orange", "darkblue"),
                               plot_fit_flag = FALSE)
  expect_is(plt, "gg")
  expect_true(grepl(metric_growth, plt[["labels"]][["y"]]))
  expect_length(plt[["layers"]], 3)
  
  expect_error(plot_dose_response_sa(dt_metrics = as.list(dt_metrics),
                                     dt_average = dt_average,
                                     grouping = grouping),
               "Check on 'dt_metrics' failed: Must be a data.table")
  expect_error(plot_dose_response_sa(dt_metrics = dt_metrics,
                                     dt_average = dt_average,
                                     grouping = "str"),
               "Check on 'grouping' failed: Must be element of set")
})

test_that("plot_sa_by_CLs works as expected", {
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  
  mae <- gDRutils::get_synthetic_data("small")
  se <- mae[[1]]
  dt_metrics <- gDRutils::convert_se_assay_to_dt(se, "Metrics")
  dt_average <- gDRutils::convert_se_assay_to_dt(se, "Averaged")
  
  cellline_name_vec <- unique(dt_metrics[[cellline_name]])[2:5]
  drug_name_vec <- unique(dt_metrics[[drug_name]])[5:7]
  
  plts <- plot_sa_by_CLs(dt_metrics = dt_metrics, 
                         dt_average = dt_average)
  expect_is(plts, "list")
  expect_equal(names(plts), paste("GR", unique(dt_metrics[[drug_name]])))
  
  metric_growth <- "RV"
  plts <- plot_sa_by_CLs(dt_metrics = dt_metrics, 
                         dt_average = dt_average,
                         cellline_name_vec = cellline_name_vec,
                         drug_name_vec = drug_name_vec,
                         metric_growth = metric_growth,
                         colormap = c("#B9D3EE", "#FF6347", "#C2F970"))
  expect_is(plts, "list")
  expect_equal(names(plts), paste("RV", drug_name_vec))
  expect_true(all(vapply(seq_along(plts), 
                         function(i) grepl(metric_growth, plts[[i]]$labels$y), logical(1))))
  
  cellline_name_vec_2 <- c(cellline_name_vec, "cellline_XX")
  drug_name_vec_2 <- c(drug_name_vec, "drug_100")
  
  plts <- plot_sa_by_CLs(dt_metrics = dt_metrics, 
                         dt_average = dt_average,
                         cellline_name_vec = cellline_name_vec_2,
                         drug_name_vec = drug_name_vec_2)
  expect_is(plts, "list")
  expect_equal(names(plts), 
               paste("GR", intersect(drug_name_vec_2, unique(dt_metrics[[drug_name]]))))
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
  
  metric_growth <- "RV"
  plt_2 <- plot_dose_response_sa_qc(dt_metrics = dt_metrics,
                                    dt_average = dt_average,
                                    cl_name = cl_name,
                                    d_name = d_name, 
                                    metric_growth = metric_growth)
  expect_is(plt_2, "gg")
  expect_true(grepl(metric_growth, plt_2[["labels"]][["y"]]))
})

test_that("plot_dose_response_sa_qc works as expected", {
  mae <- gDRutils::get_synthetic_data("small")
  se <- mae[[1]]
  
  dt_metrics <- gDRutils::convert_se_assay_to_dt(se, "Metrics")
  dt_average <- gDRutils::convert_se_assay_to_dt(se, "Averaged")
  cl_name <- dt_metrics[["CellLineName"]][1]
  d_names <- unique(dt_metrics[["DrugName"]])[1:3]
  
  plt_1 <- plot_dose_response_sa_qc_panel(dt_metrics = dt_metrics,
                                          dt_average = dt_average,
                                          cl_name = cl_name)
  expect_is(plt_1, "gg")
  expect_equal(as.character(plt_1[["layers"]][[1]][["constructor"]][[1]]), "draw_grob")
  
  metric_growth <- "RV"
  plt_2 <- plot_dose_response_sa_qc_panel(dt_metrics = dt_metrics,
                                          dt_average = dt_average,
                                          cl_name = cl_name,
                                          d_names = d_names,
                                          metric_growth = metric_growth)
  expect_is(plt_2, "gg")
})
