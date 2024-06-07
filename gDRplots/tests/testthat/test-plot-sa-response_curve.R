context("Test sa_plots")

test_that("grob_sa works as expected", {
  mae <- gDRutils::get_synthetic_data("small")
  se <- mae[[1]]
  
  grouping <- "cId"
  iR <- rownames(se)[1]
  dt_metrics <- gDRutils::convert_se_assay_to_dt(se[iR], "Metrics")
  dt_average <- gDRutils::convert_se_assay_to_dt(se[iR], "Averaged")
  
  plt <- grob_sa(dt_metrics = dt_metrics,
                 dt_average = dt_average,
                 grouping = grouping)
  expect_is(plt, "gg")
  expect_true(grepl("GR", plt[["labels"]][["y"]]))
  expect_length(plt[["layers"]], 4)
  
  grouping <- "rId"
  iC <- colnames(se)[1]
  dt_metrics <- gDRutils::convert_se_assay_to_dt(se[, iC], "Metrics")
  dt_average <- gDRutils::convert_se_assay_to_dt(se[, iC], "Averaged")
  normalization_type <- "RV"
  
  plt <- grob_sa(dt_metrics = dt_metrics,
                 dt_average = dt_average,
                 grouping = grouping,
                 normalization_type = normalization_type,
                 colormap = c("cadetblue", "orange", "darkblue"),
                 plot_fit_flag = FALSE)
  expect_is(plt, "gg")
  expect_true(grepl(normalization_type, plt[["labels"]][["y"]]))
  expect_length(plt[["layers"]], 3)
  
  expect_error(grob_sa(dt_metrics = as.list(dt_metrics),
                       dt_average = dt_average,
                       grouping = grouping),
               "Check on 'dt_metrics' failed: Must be a data.table")
  expect_error(grob_sa(dt_metrics = dt_metrics,
                       dt_average = dt_average,
                       grouping = "str"),
               "Check on 'grouping' failed: Must be element of set")
})

test_that("plot_sa_byCLs works as expected", {
  mae <- gDRutils::get_synthetic_data("small")
  se <- mae[[1]]
  cellline_name <- colnames(se)[2:5]
  drug_name <- rownames(se)[5:7]
  
  plts <- plot_sa_byCLs(se = se)
  expect_is(plts, "list")
  expect_equal(names(plts), paste("GR", rownames(se)))
  
  normalization_type <- "RV"
  
  plts <- plot_sa_byCLs(se = se,
                        cellline_name = cellline_name,
                        drug_name = drug_name,
                        normalization_type = normalization_type,
                        colormap = c("#B9D3EE", "#FF6347", "#C2F970"))
  expect_is(plts, "list")
  expect_equal(names(plts), paste("RV", drug_name))
  expect_true(all(vapply(seq_along(plts), 
                         function(i) grepl(normalization_type, plts[[i]]$labels$y), logical(1))))
  
  cellline_name_2 <- c(colnames(se)[2:3], "CL00014_cellline_XX_tissue_x_38")
  drug_name_2 <- c(rownames(se)[5:6], "G00008_drug_100_moa_A_72")
  
  plts <- plot_sa_byCLs(se = se,
                        cellline_name = cellline_name_2,
                        drug_name = drug_name_2)
  expect_is(plts, "list")
  expect_equal(names(plts), paste("GR", intersect(drug_name_2, rownames(se))))
  plotted <- intersect(cellline_name_2, colnames(se))
  expect_true(all(vapply(seq_along(plts), 
                         function(i) all(plts[[i]]$plot_env$group_names == plotted), logical(1))))
})

test_that("plot_sa_1CL works as expected", {
  mae <- gDRutils::get_synthetic_data("small")
  se <- mae[[1]]
  iC <- colnames(se)[1]
  
  plt <- plot_sa_1CL(se = se[, iC], 
                     colormap = c("cadetblue", "orange", "darkblue"))
  expect_is(plt, "gg")
  expect_true(grepl("GR", plt[["labels"]][["y"]]))
  expect_length(plt[["layers"]], 4)
  
  normalization_type <- "RV"
  
  plt <- plot_sa_1CL(se = se[, iC], 
                     normalization_type = normalization_type,
                     plot_averaged_flag = FALSE)
  expect_is(plt, "gg")
  expect_true(grepl(normalization_type, plt[["labels"]][["y"]]))
  expect_length(plt[["layers"]], 3)
})