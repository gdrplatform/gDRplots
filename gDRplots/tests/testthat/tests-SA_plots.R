context("Test SA_plots")

test_that("grob_SA works as expected", {
  mae <- gDRutils::get_synthetic_data("small")
  SE <- mae[[1]]
  
  grouping <- "cId"
  iR <- rownames(SE)[1]
  dt_metrics <- gDRutils::convert_se_assay_to_dt(SE[iR], "Metrics")
  dt_average <- gDRutils::convert_se_assay_to_dt(SE[iR], "Averaged")
  
  plt <- grob_SA(dt_metrics = dt_metrics,
                 dt_average = dt_average,
                 grouping = grouping)
  expect_is(plt, "gg")
  expect_true(grepl("GR", plt[["labels"]][["y"]]))
  expect_length(plt[["layers"]], 4)
  
  grouping <- "rId"
  iC <- colnames(SE)[1]
  dt_metrics <- gDRutils::convert_se_assay_to_dt(SE[, iC], "Metrics")
  dt_average <- gDRutils::convert_se_assay_to_dt(SE[, iC], "Averaged")
  normalization_type <- "RV"
  
  plt <- grob_SA(dt_metrics = dt_metrics,
                 dt_average = dt_average,
                 grouping = grouping,
                 normalization_type = normalization_type,
                 colormap = c("cadetblue", "orange", "darkblue"),
                 plot_fit_flag = FALSE)
  expect_is(plt, "gg")
  expect_true(grepl(normalization_type, plt[["labels"]][["y"]]))
  expect_length(plt[["layers"]], 3)
  
  expect_error(grob_SA(dt_metrics = as.list(dt_metrics),
                       dt_average = dt_average,
                       grouping = grouping),
               "Check on 'dt_metrics' failed: Must be a data.table")
  expect_error(grob_SA(dt_metrics = dt_metrics,
                       dt_average = dt_average,
                       grouping = "str"),
               "Check on 'grouping' failed: Must be element of set")
})

test_that("plot_SA_byCLs works as expected", {
  mae <- gDRutils::get_synthetic_data("small")
  SE <- mae[[1]]
  cellline_name <- colnames(SE)[2:5]
  drug_name <- rownames(SE)[5:7]
  
  plts <- plot_SA_byCLs(SE = SE)
  expect_is(plts, "list")
  expect_equal(names(plts), rownames(SE))
  
  normalization_type <- "RV"
  
  plts <- plot_SA_byCLs(SE = SE,
                       cellline_name = cellline_name,
                       drug_name = drug_name,
                       normalization_type = normalization_type,
                       colormap = c("#B9D3EE", "#FF6347", "#C2F970"))

  expect_is(plts, "list")
  expect_equal(names(plts), drug_name)
  expect_true(all(vapply(seq_along(plts), 
                         function(i) grepl(normalization_type, plts[[i]]$labels$y), logical(1))))
})

test_that("plot_SA_1CL works as expected", {
  mae <- gDRutils::get_synthetic_data("small")
  SE <- mae[[1]]
  iC <- colnames(SE)[1]
  
  plt <- plot_SA_1CL(SE = SE[, iC], 
                     colormap = c("cadetblue", "orange", "darkblue"))
  expect_is(plt, "gg")
  expect_true(grepl("GR", plt[["labels"]][["y"]]))
  expect_length(plt[["layers"]], 4)
  
  normalization_type <- "RV"
  
  plt <- plot_SA_1CL(SE = SE[, iC], 
                     normalization_type = normalization_type,
                     plot_averaged_flag = FALSE)
  expect_is(plt, "gg")
  expect_true(grepl(normalization_type, plt[["labels"]][["y"]]))
  expect_length(plt[["layers"]], 3)
})