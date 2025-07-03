context("Test dose-response sa")

test_that("plot_dose_response_sa works as expected", {
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  
  mae <- gDRutils::get_synthetic_data("small")
  se <- mae[[gDRutils::get_supported_experiments("sa")]]
  
  group_var <- cellline_name
  selected_drug <- "drug_002"
  dt_metrics <- gDRutils::convert_se_assay_to_dt(se, "Metrics")
  dt_average <- gDRutils::convert_se_assay_to_dt(se, "Averaged")
  
  plt_1 <- plot_dose_response_sa(dt_metrics = dt_metrics,
                                 dt_average = dt_average,
                                 selection_name = selected_drug,
                                 group_var = group_var) # default
  expect_is(plt_1, "gg")
  expect_equal(plt_1[["labels"]][["y"]], "GR")
  expect_true(grepl(selected_drug,  plt_1[["labels"]][["title"]]))
  expect_length(plt_1[["layers"]], 3)
  expect_equal(ggplot2::get_guide_data(plt_1, "colour")[[".label"]],
               unique(dt_average[[group_var]]))
  expect_equal(plt_1[["labels"]][["colour"]], "group_var")
  
  normalization_type <- "RV"
  plt_2 <- plot_dose_response_sa(dt_metrics = dt_metrics,
                                 dt_average = dt_average,
                                 selection_name = selected_drug,
                                 group_var = group_var,
                                 normalization_type = normalization_type,
                                 colors_vec = rainbow(NROW(unique(dt_metrics[[group_var]]))))
  expect_is(plt_2, "gg")
  expect_equal(plt_2[["labels"]][["y"]], normalization_type)
  expect_length(unique(ggplot2::ggplot_build(plt_2)$data[[2]][["colour"]]),
                NROW(unique(dt_metrics[[group_var]])))
  
  subset_celline <- c("cellline_BA", "cellline_CA")
  # lack of data for curve and small number of conc point
  dt_metrics_lack <- data.table::copy(dt_metrics)
  dt_metrics_lack[CellLineName %in% subset_celline, c("ec50", "x_inf", "x_0", "h", "r2")] <- NA
  dt_average_lack <- 
    data.table::copy(dt_average)[Concentration %in% unique(dt_average$Concentration)[1:2], ]
  sel_grp_names <- c(subset_celline, "cellline_KB")
  
  plt_3 <- plot_dose_response_sa(dt_metrics = dt_metrics_lack,
                                 dt_average = dt_average,
                                 selection_name = selected_drug,
                                 group_names = sel_grp_names,
                                 group_var = group_var)
  expect_is(plt_3, "gg")
  expect_equal(plt_3[["labels"]][["y"]], "GR")
  expect_length(plt_3[["layers"]], 3)
  expect_equal(names(plt_3[["guides"]][["guides"]]), "colour")
  expect_equal(sort(ggplot2::get_guide_data(plt_3, "colour")[[".label"]]), sort(sel_grp_names))
  expect_length(unique(ggplot2::ggplot_build(plt_3)$data[[3]]$colour), 1) # curve data only for 1 cell line
  
  # scenario: lack of metric data at all -> plot only observation
  plt_4 <- plot_dose_response_sa(dt_metrics = dt_metrics[normalization_type == "RV", ],
                                 dt_average = dt_average,
                                 selection_name = selected_drug,
                                 group_var = group_var)
  expect_is(plt_4, "gg")
  expect_equal(plt_4[["labels"]][["y"]], "GR")
  expect_true(any(grepl("Concentration", plt_4[["labels"]][["x"]])))
  expect_length(plt_4[["layers"]], 2)
  expect_length(ggplot2::ggplot_build(plt_4)[["data"]], 2) # lack of fit lines
  
  # scenario: plot by drugs 
  group_var <- drug_name
  selected_celline <- "cellline_BA"
  normalization_type <- "RV"
  ls_col <- c("#5F9EA0", "#FFA500", "#00008B")
  
  plt_5 <- plot_dose_response_sa(dt_metrics = dt_metrics,
                                 dt_average = dt_average,
                                 selection_name = selected_celline,
                                 group_var = group_var,
                                 normalization_type = normalization_type,
                                 colors_vec = ls_col,
                                 plot_fit_flag = FALSE)
  expect_is(plt_5, "gg")
  expect_equal(plt_5[["labels"]][["y"]], normalization_type)
  expect_true(grepl(selected_celline, plt_5[["labels"]][["title"]]))
  expect_length(plt_5[["layers"]], 2) # lack of fit lines
  expect_length(ggplot2::ggplot_build(plt_5)[["data"]], 2) # lack of fit lines
  expect_true(all(ls_col[c(1, 3)] %in% unique(ggplot2::ggplot_build(plt_5)$data[[2]][["colour"]])))
  
  # scenario: lack of metric data for selected `group_names`
  drug_name_subset <- c("drug_002", "drug_003", "drug_004", "drug_005", "drug_006")
  dt_metrics_lack <- data.table::copy(dt_metrics)[DrugName %in% drug_name_subset]
  drug_name_vec <- c("drug_005", "drug_006", "drug_007", "drug_008", "drug_009", "drug_010")
  
  plt_6 <- plot_dose_response_sa(dt_metrics = dt_metrics_lack,
                                 dt_average = dt_average,
                                 selection_name = selected_celline,
                                 group_var = group_var,
                                 group_names = drug_name_vec,
                                 normalization_type = normalization_type)
  expect_is(plt_6, "gg")
  expect_length(plt_6[["layers"]], 3)
  expect_equal(NROW(unique(ggplot2::ggplot_build(plt_6)[["data"]][[2]][["colour"]])), 
               NROW(drug_name_vec))  # avg 
  expect_equal(NROW(unique(ggplot2::ggplot_build(plt_6)[["data"]][[3]][["colour"]])), 
               NROW(intersect(drug_name_subset, drug_name_vec))) # metric 
  
  # scenario: lack of averager data for selected `group_names`
  drug_name_subset <- c("drug_002", "drug_003", "drug_004", "drug_005", "drug_006")
  dt_average_lack <- data.table::copy(dt_average)[get(drug_name) %in% drug_name_subset]
  drug_name_vec <- c("drug_005", "drug_006", "drug_007", "drug_008", "drug_009", "drug_010")
  
  plt_7 <- plot_dose_response_sa(dt_metrics = dt_metrics,
                                 dt_average = dt_average_lack,
                                 selection_name = selected_celline,
                                 group_var = group_var,
                                 group_names = drug_name_vec,
                                 normalization_type = normalization_type)
  expect_is(plt_7, "gg")
  expect_length(plt_7[["layers"]], 3)
  expect_equal(NROW(unique(ggplot2::ggplot_build(plt_7)[["data"]][[2]][["colour"]])), 
               NROW(intersect(drug_name_subset, drug_name_vec))) # avg 
  expect_equal(NROW(unique(ggplot2::ggplot_build(plt_7)[["data"]][[3]][["colour"]])), 
               NROW(drug_name_vec)) # metric
  
  # scenario: one drug is not available in data - plot for drugs
  plt_8 <- plot_dose_response_sa(dt_metrics = dt_metrics,
                                 dt_average = dt_average,
                                 selection_name = selected_celline,
                                 group_var = group_var,
                                 group_names = c(drug_name_vec, "drug_xx"))
  expect_is(plt_8, "gg")
  expect_equal(ggplot2::get_guide_data(plt_8, "colour")[[".label"]],
               drug_name_vec) # legend is present without "drug_xx"
  
  # scenario: NA in concentration column -> only fitted curve
  dt_average_NA <- data.table::copy(dt_average)
  dt_average_NA[["Concentration"]] <- NA
  
  plt_9 <- plot_dose_response_sa(dt_metrics = dt_metrics,
                                 dt_average = dt_average_NA,
                                 selection_name = selected_celline,
                                 group_var = group_var)
  expect_is(plt_9, "gg")
  expect_equal(plt_9[["labels"]][["y"]], "GR")
  expect_true(any(grepl("Concentration", plt_9[["labels"]][["x"]])))
  expect_true(grepl(selected_celline, plt_9[["labels"]][["title"]]))
  expect_length(ggplot2::ggplot_build(plt_9)$data[[1]], 7)
  expect_equal(ggplot2::get_guide_data(plt_9, "colour")[[".label"]],
               sort(unique(dt_metrics[["DrugName"]])))
  
  # scenario: values for drug are NAs for selected cell line
  dt_metrics_NA <- data.table::copy(dt_metrics)
  ls_col_met <- intersect(names(dt_metrics_NA), gDRutils::get_header("response_metrics"))
  dt_metrics_NA[get(cellline_name) == selected_celline, (ls_col_met) := NA]
  
  dt_average_NA <- data.table::copy(dt_average)
  ls_col_avg <- intersect(names(dt_average_NA), gDRutils::get_header("averaged_results"))
  dt_average_NA[get(cellline_name) == selected_celline, (ls_col_avg) := NA]
  
  plt_10 <- plot_dose_response_sa(dt_metrics = dt_metrics_NA,
                                  dt_average = dt_average_NA,
                                  selection_name = selected_celline,
                                  group_var = group_var,
                                  group_names = drug_name_vec)
  expect_is(plt_10, "gg")
  expect_equal(plt_10[["labels"]][["y"]], "GR")
  expect_true(any(grepl("Concentration", plt_10[["labels"]][["x"]])))
  expect_true(grepl(selected_celline, plt_10[["labels"]][["title"]]))
  expect_length(ggplot2::ggplot_build(plt_10)$data[[1]], 0) # no data at al
  
  # scenario: there is no data in metrics for selected cell line -> plot only observations
  dt_metrics_miss <- data.table::copy(dt_metrics)
  dt_metrics_miss <- dt_metrics_miss[get(cellline_name) != selected_celline, ]
  
  plt_11 <- plot_dose_response_sa(dt_metrics = dt_metrics_miss,
                                  dt_average = dt_average,
                                  selection_name = selected_celline,
                                  group_var = group_var,
                                  group_names = drug_name_vec)
  
  expect_is(plt_11, "gg")
  expect_true(grepl(selected_celline, plt_11[["labels"]][["title"]]))
  expect_length(ggplot2::ggplot_build(plt_11)$data[[1]], 7) 
  
  # scenario: there is no data in average for selected cell line -> only fiited curve
  dt_average_miss <- data.table::copy(dt_average)
  dt_average_miss <- dt_average_miss[get(cellline_name) != selected_celline, ]
  
  plt_12 <- plot_dose_response_sa(dt_metrics = dt_metrics,
                                  dt_average = dt_average_miss,
                                  selection_name = selected_celline,
                                  group_var = group_var,
                                  group_names = drug_name_vec)
  expect_is(plt_12, "gg")
  expect_true(grepl(selected_celline, plt_12[["labels"]][["title"]]))
  expect_length(ggplot2::ggplot_build(plt_12)$data[[1]], 7) # no data at all
  expect_equal(ggplot2::get_guide_data(plt_12, "colour")[[".label"]],
               drug_name_vec)
  
  # scenario: there is no data in metrics and average for selected cell line
  plt_13 <- plot_dose_response_sa(dt_metrics = dt_metrics_miss,
                                  dt_average = dt_average_miss,
                                  selection_name = selected_celline,
                                  group_var = group_var,
                                  group_names = drug_name_vec)
  
  expect_is(plt_13, "gg")
  expect_equal(plt_13[["labels"]][["title"]], selected_celline)
  expect_length(ggplot2::ggplot_build(plt_13)$data[[1]], 0) # no data at all
  
  # scenario: fitted curve has bigger y-range than avg points
  drug_nm <- "drug_0A"
  cl_nm <- "cl_14"
  tab_met <- data.table::data.table(
    xc50 = c(2.688, 3.755),
    x_inf = c(0, -1),
    x_0 = c(1, 1),
    ec50 = c(2.688, 5.016),
    h = c(5.000, 3.794),
    normalization_type = c("RV", "GR")
  )
  tab_met$DrugName <- drug_nm
  tab_met$Gnumber <- "G0012"
  tab_met$CellLineName <- cl_nm
  tab_met$clid <- "CLXX"
  
  tab_avg <- data.table::data.table(
    x = c(0.971, 0.989, 0.961, 0.985, 0.834, 0.933, 0.846, 0.938, 1,
          1, 1.012, 1.004, 0.884, 0.954, 0.926, 0.971, 0.026, 0.006),
    normalization_type = rep(c("RV", "GR"), 9),
    Concentration = rep(c(0.0008, 0.0023, 0.0069, 0.0206, 0.0617, 
                          0.1852, 0.5556, 1.6667, 5.0000), each = 2)
  )
  tab_avg$DrugName <- drug_nm
  tab_avg$CellLineName <- cl_nm
  
  conc_range <- range(tab_avg[normalization_type == "GR", ]$Concentration)
  sel_conc <- 10 ^ (seq(conc_range[1], conc_range[2], 0.05))
  sel_metrics <- tab_met[normalization_type == "GR", ]
  fitted_range <- range(gDRutils::predict_efficacy_from_conc(sel_conc,
                                                             sel_metrics$x_inf,
                                                             sel_metrics$x_0,
                                                             sel_metrics$ec50,
                                                             sel_metrics$h))
  
  plt_14 <- plot_dose_response_sa(dt_metrics = tab_met, 
                                  dt_average = tab_avg,
                                  selection_name = "drug_0A",
                                  group_var = cellline_name,
                                  group_names = "cl_14",
                                  normalization_type = "GR")
  expect_is(plt_14, "gg")
  expect_true(grepl(drug_nm, plt_14[["labels"]][["title"]]))
  plot_range <- range(as.numeric(ggplot2::layer_scales(plt_14)$y$get_labels()))
  expect_true(all(data.table::between(fitted_range, plot_range[1], plot_range[2])))
  
  # testing assertion
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
  expect_error(plot_dose_response_sa(dt_metrics = dt_metrics,
                                     dt_average = dt_average,
                                     selection_name = selected_celline,
                                     group_var = group_var,
                                     group_names = 1:3),
               "Assertion on 'group_names' failed: Must be of type 'character'")
  expect_error(plot_dose_response_sa(dt_metrics = dt_metrics,
                                     dt_average = dt_average,
                                     selection_name = selected_celline,
                                     group_var = group_var,
                                     normalization_type = "str"),
               "Assertion on 'normalization_type' failed: Must be element of set")
  expect_error(plot_dose_response_sa(dt_metrics = dt_metrics,
                                     dt_average = dt_average,
                                     selection_name = selected_celline,
                                     group_var = group_var,
                                     colors_vec = 1:10),
               "Assertion on 'colors_vec' failed: Must be of type 'character'")
  expect_error(plot_dose_response_sa(dt_metrics = dt_metrics,
                                     dt_average = dt_average,
                                     selection_name = selected_celline,
                                     group_var = group_var,
                                     plot_averaged_flag = "str"),
               "Assertion on 'plot_averaged_flag' failed: Must be of type 'logical flag'")
  expect_error(plot_dose_response_sa(dt_metrics = dt_metrics,
                                     dt_average = dt_average,
                                     selection_name = selected_celline,
                                     group_var = group_var,
                                     plot_fit_flag = NULL),
               "Assertion on 'plot_fit_flag' failed: Must be of type 'logical flag'")
  expect_error(plot_dose_response_sa(dt_metrics = dt_metrics,
                                     dt_average = dt_average,
                                     selection_name = selected_celline,
                                     group_var = group_var,
                                     fit_source = 1),
               "Assertion on 'fit_source' failed: Must be of type 'string'")
})

test_that("plot_dose_response_sa_by_CLs works as expected", {
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  
  mae <- gDRutils::get_synthetic_data("small")
  se <- mae[[gDRutils::get_supported_experiments("sa")]]
  dt_metrics <- gDRutils::convert_se_assay_to_dt(se, "Metrics")
  dt_average <- gDRutils::convert_se_assay_to_dt(se, "Averaged")
  
  cellline_name_vec <- unique(dt_metrics[[cellline_name]])[2:5]
  drug_name_vec <- unique(dt_metrics[[drug_name]])[5:7]
  
  plts_1 <- plot_dose_response_sa_by_CLs(dt_metrics = dt_metrics, 
                                         dt_average = dt_average)
  expect_is(plts_1, "list")
  expect_equal(names(plts_1), unique(dt_metrics[[drug_name]]))
  expect_true(all(vapply(names(plts_1), 
                         function(nm) grepl(nm, plts_1[[nm]][["labels"]][["title"]]), logical(1))))
  
  normalization_type <- "RV"
  plts_2 <- plot_dose_response_sa_by_CLs(dt_metrics = dt_metrics, 
                                         dt_average = dt_average,
                                         cellline_name_vec = cellline_name_vec,
                                         drug_name_vec = drug_name_vec,
                                         normalization_type = normalization_type,
                                         colors_vec = c("#00008B", "#FF6347", "#4CBB17"))
  expect_is(plts_2, "list")
  plotted_ <- intersect(drug_name_vec, unique(dt_metrics[[drug_name]]))
  expect_equal(names(plts_2), plotted_)
  expect_true(all(vapply(seq_along(plts_2), 
                         function(i) grepl(normalization_type, plts_2[[i]]$labels$y), logical(1))))
  expect_true(all(vapply(names(plts_2), 
                         function(nm) grepl(nm, plts_2[[nm]][["labels"]][["title"]]), logical(1))))
  
  
  # scenario: selected cell line is not available in data
  cellline_name_vec_2 <- c(cellline_name_vec, "cellline_XX")
  drug_name_vec_2 <- c(drug_name_vec, "drug_100")
  
  plts_3 <- plot_dose_response_sa_by_CLs(dt_metrics = dt_metrics, 
                                         dt_average = dt_average,
                                         cellline_name_vec = cellline_name_vec_2,
                                         drug_name_vec = drug_name_vec_2)
  expect_is(plts_3, "list")
  plotted_ <- intersect(drug_name_vec_2, unique(dt_metrics[[drug_name]]))
  expect_equal(names(plts_3), plotted_)
  plotted <- intersect(cellline_name_vec_2, unique(dt_metrics[[cellline_name]])) # only available
  expect_true(all(vapply(seq_along(plts_3), 
                         function(i) all(plts_3[[i]]$plot_env$group_names == plotted), logical(1))))
  
  # scenario: values for selected drug are NAs
  dt_metrics_NA <- data.table::copy(dt_metrics)
  ls_col_met <- intersect(names(dt_metrics_NA), gDRutils::get_header("response_metrics"))
  dt_metrics_NA[get(drug_name) == drug_name_vec[1], (ls_col_met) := NA]
  
  dt_average_NA <- data.table::copy(dt_average)
  ls_col_avg <- intersect(names(dt_average_NA), gDRutils::get_header("averaged_results"))
  dt_average_NA[get(drug_name) == drug_name_vec[1], (ls_col_avg) := NA]
  
  plts_4 <- plot_dose_response_sa_by_CLs(dt_metrics = dt_metrics_NA, 
                                         dt_average = dt_average_NA,
                                         cellline_name_vec = cellline_name_vec,
                                         drug_name_vec = drug_name_vec)
  expect_is(plts_4, "list")
  expect_equal(names(plts_4), drug_name_vec)
  expect_length(unique(ggplot2::ggplot_build(plts_4[[drug_name_vec[2]]])[["data"]][[2]][["colour"]]),
                NROW(cellline_name_vec))
  expect_length(ggplot2::ggplot_build(plts_4[[drug_name_vec[1]]])[["data"]][[1]], 0) # no data for 1st drug
  expect_equal(ggplot2::get_guide_data(plts_4[[drug_name_vec[2]]], "colour")[[".label"]],
               cellline_name_vec) # legend is present
})

test_that("plot_dose_response_sa_by_drugs works as expected", {
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  clid <- gDRutils::get_env_identifiers("cellline")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  
  mae <- gDRutils::get_synthetic_data("small")
  se <- mae[[gDRutils::get_supported_experiments("sa")]]
  dt_metrics <- gDRutils::convert_se_assay_to_dt(se, "Metrics")
  dt_average <- gDRutils::convert_se_assay_to_dt(se, "Averaged")
  
  cellline_name_vec <- unique(dt_metrics[[cellline_name]])[2:5]
  drug_name_vec <- unique(dt_metrics[[drug_name]])[5:7]
  
  plts_1 <- plot_dose_response_sa_by_drugs(dt_metrics = dt_metrics, 
                                           dt_average = dt_average)
  expect_is(plts_1, "list")
  expect_equal(names(plts_1), unique(dt_metrics[[cellline_name]]))
  expect_true(all(vapply(names(plts_1), 
                         function(nm) grepl(nm, plts_1[[nm]][["labels"]][["title"]]), logical(1))))
  
  normalization_type <- "RV"
  plts_2 <- plot_dose_response_sa_by_drugs(dt_metrics = dt_metrics, 
                                           dt_average = dt_average,
                                           cellline_name_vec = cellline_name_vec,
                                           drug_name_vec = drug_name_vec,
                                           normalization_type = normalization_type,
                                           colors_vec = c("#00008B", "#FF6347", "#4CBB17"))
  expect_is(plts_2, "list")
  plotted_ <- intersect(cellline_name_vec, unique(dt_metrics[[cellline_name]]))
  expect_equal(names(plts_2), plotted_)
  expect_true(all(vapply(seq_along(plts_2), 
                         function(i) grepl(normalization_type, plts_2[[i]]$labels$y), logical(1))))
  expect_true(all(vapply(names(plts_2), 
                         function(nm) grepl(nm, plts_2[[nm]][["labels"]][["title"]]), logical(1))))
  
  cellline_name_vec_2 <- c(cellline_name_vec, "cellline_XX")
  drug_name_vec_2 <- c(drug_name_vec, "drug_100")
  
  plts_3 <- plot_dose_response_sa_by_drugs(dt_metrics = dt_metrics, 
                                           dt_average = dt_average,
                                           cellline_name_vec = cellline_name_vec_2,
                                           drug_name_vec = drug_name_vec_2)
  expect_is(plts_3, "list")
  plotted_ <- intersect(cellline_name_vec_2, unique(dt_metrics[[cellline_name]]))
  expect_equal(names(plts_3),  plotted_)
  plotted <- intersect(drug_name_vec_2, unique(dt_metrics[[drug_name]]))
  expect_true(all(vapply(seq_along(plts_3), 
                         function(i) all(plts_3[[i]]$plot_env$group_names == plotted), logical(1))))
  
  # scenario: values for selected cell line are NAs
  dt_metrics_NA <- data.table::copy(dt_metrics)
  ls_col_met <- intersect(names(dt_metrics_NA), gDRutils::get_header("response_metrics"))
  dt_metrics_NA[get(cellline_name) == cellline_name_vec[1], (ls_col_met) := NA]
  
  dt_average_NA <- data.table::copy(dt_average)
  ls_col_avg <- intersect(names(dt_average_NA), gDRutils::get_header("averaged_results"))
  dt_average_NA[get(cellline_name) == cellline_name_vec[1], (ls_col_avg) := NA]
  
  plts_4 <- plot_dose_response_sa_by_drugs(dt_metrics = dt_metrics_NA, 
                                           dt_average = dt_average_NA,
                                           cellline_name_vec = cellline_name_vec,
                                           drug_name_vec = drug_name_vec)
  expect_is(plts_4, "list")
  expect_equal(names(plts_4), cellline_name_vec)
  expect_length(unique(ggplot2::ggplot_build(plts_4[[cellline_name_vec[2]]])[["data"]][[2]][["colour"]]),
                NROW(drug_name_vec))
  expect_length(ggplot2::ggplot_build(plts_4[[cellline_name_vec[1]]])[["data"]][[1]], 0) # no data for 1st cell line
  expect_equal(ggplot2::get_guide_data(plts_4[[cellline_name_vec[2]]], "colour")[[".label"]],
               drug_name_vec) # legend is present
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
  
  # empty plot
  plt_5 <- plot_dose_response_sa_qc(dt_metrics = dt_metrics[normalization_type == "RV", ],
                                    dt_average = dt_average[normalization_type == "RV", ],
                                    cl_name = cl_name,
                                    d_name = d_name)
  expect_is(plt_5, "gg")
  expect_equal(plt_5[["labels"]][["y"]], "GR")
  expect_true(grepl(cl_name, plt_5[["labels"]][["title"]]))
  expect_true(any(grepl("Concentration", plt_5[["labels"]][["x"]])))
  expect_length(ggplot2::ggplot_build(plt_5)[["data"]][[1]], 0)
})

test_that("plot_dose_response_sa_qc_panel works as expected", {
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  
  mae <- gDRutils::get_synthetic_data("small")
  se <- mae[[gDRutils::get_supported_experiments("sa")]]
  
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
    dt_average[get(drug_name) %in% ls_drug & get(cellline_name) == cl_name, c(cellline_name, drug_name), with = FALSE]
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
