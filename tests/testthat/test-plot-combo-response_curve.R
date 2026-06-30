context("Test dose-response combo")

test_that("plot_dose_response_combo works as expected", {
  mae <- gDRutils::get_synthetic_data("combo_matrix")
  se <- mae[[gDRutils::get_supported_experiments("combo")]]
  dt_average <- gDRutils::convert_se_assay_to_dt(se, "Averaged")

  cl_name <- "cellline_BC"
  drug1_name <- "drug_011"
  drug2_name <- "drug_021"
  cl_clid <- unique(dt_average[CellLineName == cl_name, ]$clid)

  plt_1 <- plot_dose_response_combo(dt_average = dt_average,
                                    drug1_name = drug1_name,
                                    drug2_name = drug2_name,
                                    cl_name = cl_name)
  expect_is(plt_1, "gg")
  expect_equal(plt_1[["labels"]][["y"]], "GR")
  expect_length(plt_1[["layers"]], 3)
  expect_equal(plt_1[["labels"]][["title"]], sprintf("%s (%s)", cl_name, cl_clid))

  normalization_type <- "RV"
  plt_2 <- plot_dose_response_combo(dt_average = dt_average,
                                    drug1_name = drug1_name,
                                    drug2_name = drug2_name,
                                    cl_name = cl_name,
                                    normalization_type = normalization_type)
  expect_is(plt_2, "gg")
  expect_equal(plt_2[["labels"]][["y"]], normalization_type)

  plt_3 <- plot_dose_response_combo(dt_average = dt_average,
                                    drug1_name = drug1_name,
                                    drug2_name = drug2_name,
                                    cl_name = cl_name,
                                    colors_vec = c("#FFA500", "#8B0000"))
  expect_is(plt_3, "gg")
  expect_equal(plt_3[["labels"]][["y"]], "GR")
  expect_true(any(grepl("#FFA500", plt_3[["plot_env"]][["colormap"]])))
  expect_true(any(grepl("#8B0000", plt_3[["plot_env"]][["colormap"]])))

  plt_4 <- plot_dose_response_combo(dt_average = dt_average,
                                    drug1_name = drug1_name,
                                    drug2_name = drug2_name,
                                    cl_name = cl_name,
                                    colors_vec = c("pinky", "blackish"))
  expect_is(plt_4, "gg")
  expect_equal(plt_4[["labels"]][["y"]], "GR")
  expect_true(all(plt_4[["plot_env"]][["colormap"]] ==
                    .get_combo_curves_colors(as.factor(unique(dt_average[["Concentration_2"]])))
  )) # default colors when invalid `colors_vec`

  plt_5 <- plot_dose_response_combo(dt_average = dt_average,
                                    drug1_name = drug1_name,
                                    drug2_name = drug2_name,
                                    cl_name = cl_name,
                                    split_by_conc = TRUE)
  expect_is(plt_5, "gg")
  expect_equal(plt_5[["labels"]][["y"]], "GR")
  expect_equal(plt_5[["labels"]][["title"]], sprintf("%s (%s)", cl_name, cl_clid))
  expect_true(NROW(plt_5[["facet"]][["params"]]) > 0) # plot is faceted
  expect_true(grepl("conc_2", names(plt_5[["facet"]][["params"]][["facets"]])))

  sel_conc <- c(0, 0.1, 1.0)
  sel_colors <- c("#00008B", "#FF8C00", "#008B8B")
  dt_average_conc <- dt_average[Concentration_2 %in% sel_conc]
  plt_6 <- plot_dose_response_combo(dt_average = dt_average_conc,
                                    drug1_name = drug1_name,
                                    drug2_name = drug2_name,
                                    cl_name = cl_name,
                                    normalization_type = normalization_type,
                                    colors_vec = sel_colors)
  expect_is(plt_6, "gg")
  expect_equal(plt_6[["labels"]][["y"]], "RV")
  expect_equal(as.numeric(names(plt_6[["plot_env"]][["colormap"]])), sel_conc)
  expect_equal(unname(plt_6[["plot_env"]][["colormap"]]), sel_colors)

  expect_error(plot_dose_response_combo(dt_average = unlist(dt_average),
                                        drug1_name = drug1_name,
                                        drug2_name = drug2_name,
                                        cl_name = cl_name),
               "Assertion on 'dt_average' failed: Must be a data.table")

  expect_error(plot_dose_response_combo(dt_average = dt_average,
                                        drug1_name = 1,
                                        drug2_name = drug2_name,
                                        cl_name = cl_name),
               "Assertion on 'drug1_name' failed: Must be element of set")

  expect_error(plot_dose_response_combo(dt_average = dt_average,
                                        drug1_name = drug1_name,
                                        drug2_name = "drug_XX",
                                        cl_name = cl_name),
               "Assertion on 'drug2_name' failed: Must be element of set")

  expect_error(plot_dose_response_combo(dt_average = dt_average,
                                        drug1_name = drug1_name,
                                        drug2_name = drug2_name,
                                        cl_name = 1),
               "Assertion on 'cl_name' failed: Must be of type 'string'")

  expect_error(plot_dose_response_combo(dt_average = dt_average,
                                        drug1_name = drug1_name,
                                        drug2_name = drug2_name,
                                        cl_name = "cellline_XX"),
               "Assertion on 'cl_name' failed: Must be element of set")

  expect_error(plot_dose_response_combo(dt_average = dt_average,
                                        drug1_name = drug1_name,
                                        drug2_name = drug2_name,
                                        cl_name = cl_name,
                                        normalization_type = "AA"),
               "Assertion on 'normalization_type' failed: Must be element of set")

  expect_error(plot_dose_response_combo(dt_average = dt_average,
                                        drug1_name = drug1_name,
                                        drug2_name = drug2_name,
                                        cl_name = cl_name,
                                        colors_vec = 1:5),
               "Assertion on 'colors_vec' failed: Must be of type 'character'")

  expect_error(plot_dose_response_combo(dt_average = dt_average,
                                        drug1_name = drug1_name,
                                        drug2_name = drug2_name,
                                        cl_name = cl_name,
                                        split_by_conc = 1),
               "Assertion on 'split_by_conc' failed: Must be of type 'logical flag'")

  # auto-orient: reduce drug1 to 2 dose levels so drug2 has more
  concs_d1 <- sort(unique(
    dt_average[DrugName == drug1_name & Concentration > 0]$Concentration
  ))
  dt_asym <- dt_average[
    !(DrugName == drug1_name & Concentration > 0 &
        !Concentration %in% concs_d1[1:2])
  ]
  plt_7 <- plot_dose_response_combo(dt_average = dt_asym,
                                    drug1_name = drug1_name,
                                    drug2_name = drug2_name,
                                    cl_name = cl_name)
  expect_is(plt_7, "gg")
  expect_true(grepl(drug2_name, deparse(plt_7[["labels"]][["x"]])))
})

test_that("plot_dose_response_combo_panel works as expected", {
  mae <- gDRutils::get_synthetic_data("combo_matrix")
  se <- mae[[gDRutils::get_supported_experiments("combo")]]
  dt_average <- gDRutils::convert_se_assay_to_dt(se, "Averaged")

  cl_name <- "cellline_IB"
  d_names <- unique(dt_average$DrugName)[2:3]

  no_comb_all <- NROW(unique(
    dt_average[CellLineName == cl_name,  c("CellLineName", "DrugName", "DrugName_2"), with = FALSE]
  ))
  no_comb <- NROW(unique(
    dt_average[DrugName %in% d_names & CellLineName == cl_name,
               c("CellLineName", "DrugName", "DrugName_2"), with = FALSE]
  ))

  plt_1 <- plot_dose_response_combo_panel(dt_average = dt_average,
                                          cl_name = cl_name)
  expect_is(plt_1, "gg")
  expect_equal(plt_1[["labels"]][["y"]], "GR")
  expect_length(plt_1[["layers"]], 3)
  expect_true(grepl(cl_name, plt_1[["labels"]][["title"]]))

  normalization_type <- "RV"
  plt_2 <- plot_dose_response_combo_panel(dt_average = dt_average,
                                          cl_name = cl_name,
                                          d_names = d_names,
                                          normalization_type = normalization_type)
  expect_is(plt_2, "gg")
  expect_equal(plt_2[["labels"]][["y"]], normalization_type)
  expect_equal(NROW(ggplot2::ggplot_build(plt_2)$data[[1]]), 2 * no_comb) # 0 & 1 line

  plt_3 <- plot_dose_response_combo_panel(dt_average = dt_average,
                                          cl_name = cl_name,
                                          colors_vec = c("#FFA500", "#8B0000"))
  expect_is(plt_3, "gg")
  expect_equal(plt_1[["labels"]][["y"]], "GR")
  expect_true(any(grepl("#FFA500", plt_3[["plot_env"]][["colormap"]])))
  expect_true(any(grepl("#8B0000", plt_3[["plot_env"]][["colormap"]])))

  plt_4 <- plot_dose_response_combo_panel(dt_average = dt_average,
                                          cl_name = cl_name,
                                          colors_vec = c("pinky", "blackish"))
  expect_is(plt_4, "gg")
  expect_equal(plt_4[["labels"]][["y"]], "GR")
  expect_true(all(plt_4[["plot_env"]][["colormap"]] ==
                    .get_combo_curves_colors(as.factor(unique(dt_average[["Concentration_2"]])))
  )) # default colors when invalid `colors_vec`

  plt_5 <- plot_dose_response_combo_panel(dt_average = dt_average,
                                          d_names = c("drug_XX", "drug_YY"),
                                          cl_name = cl_name)
  expect_is(plt_5, "gg")
  expect_equal(plt_5[["labels"]][["y"]], "GR")
  expect_equal(NROW(ggplot2::ggplot_build(plt_5)$data[[1]]), 2 * no_comb_all) # default all drugs when invalid `d_names`

  ls_drug <- c(d_names, "drug_YY")
  no_comb_err <- NROW(unique(
    dt_average[DrugName %in% ls_drug & CellLineName == cl_name,
               c("CellLineName", "DrugName", "DrugName_2"), with = FALSE]
  ))

  plt_6 <- plot_dose_response_combo_panel(dt_average = dt_average,
                                          d_names = ls_drug,
                                          cl_name = cl_name)
  expect_is(plt_6, "gg")
  expect_equal(plt_6[["labels"]][["y"]], "GR")
  expect_equal(NROW(ggplot2::ggplot_build(plt_6)$data[[1]]), 2 * no_comb_err)

  expect_error(plot_dose_response_combo_panel(dt_average = unlist(dt_average),
                                              cl_name = cl_name),
               "Assertion on 'dt_average' failed: Must be a data.table")

  expect_error(plot_dose_response_combo_panel(dt_average = dt_average,
                                              d_names = 1,
                                              cl_name = cl_name),
               "Assertion on 'd_names' failed: Must be of type 'character'")

  expect_error(plot_dose_response_combo_panel(dt_average = dt_average,
                                              cl_name = 1),
               "Assertion on 'cl_name' failed: Must be of type 'string'")

  expect_error(plot_dose_response_combo_panel(dt_average = dt_average,
                                              cl_name = "cellline_XX"),
               "Assertion on 'cl_name' failed: Must be element of set")

  expect_error(plot_dose_response_combo_panel(dt_average = dt_average,
                                              cl_name = cl_name,
                                              normalization_type = "AA"),
               "Assertion on 'normalization_type' failed: Must be element of set")

  expect_error(plot_dose_response_combo_panel(dt_average = dt_average,
                                              cl_name = cl_name,
                                              colors_vec = 1:5),
               "Assertion on 'colors_vec' failed: Must be of type 'character'")

  # auto-orient: asymmetric concentrations across multiple drug pairs
  cl_name_2 <- "cellline_BC"
  drug1_name_2 <- "drug_011"
  drug2_name_2 <- "drug_026"
  concs_d1 <- sort(unique(
    dt_average[DrugName == drug1_name_2 & Concentration > 0]$Concentration
  ))
  concs_d2 <- sort(unique(
    dt_average[DrugName_2 == drug2_name_2 & Concentration_2 > 0]$Concentration_2
  ))
  dt_asym <- dt_average[
    !(DrugName == drug1_name_2 & Concentration > 0 &
        !Concentration %in% concs_d1[1:2]) &
      !(DrugName_2 == drug2_name_2 & Concentration_2 > 0 &
          !Concentration_2 %in% concs_d2[3:4])
  ]
  plt_7 <- plot_dose_response_combo_panel(dt_average = dt_asym,
                                          cl_name = cl_name_2)
  expect_is(plt_7, "gg")
})

test_that("plot_dose_response_combo_panel handles RV-only cell line with GR request", {
  mae <- gDRutils::get_synthetic_data("combo_matrix")
  se <- mae[[gDRutils::get_supported_experiments("combo")]]
  dt_average <- gDRutils::convert_se_assay_to_dt(se, "Averaged")

  cl_name <- "cellline_BC"
  dt_rv_only <- dt_average[
    !(CellLineName == cl_name & normalization_type == "GR")
  ]

  result <- plot_dose_response_combo_panel(
    dt_average = dt_rv_only,
    cl_name = cl_name,
    normalization_type = "GR"
  )
  expect_null(result)
})

test_that("plot_dose_response_combo handles RV-only cell line with GR request", {
  mae <- gDRutils::get_synthetic_data("combo_matrix")
  se <- mae[[gDRutils::get_supported_experiments("combo")]]
  dt_average <- gDRutils::convert_se_assay_to_dt(se, "Averaged")

  cl_name <- "cellline_BC"
  drug1_name <- "drug_011"
  drug2_name <- "drug_021"
  dt_rv_only <- dt_average[
    !(CellLineName == cl_name & normalization_type == "GR")
  ]

  result <- plot_dose_response_combo(
    dt_average = dt_rv_only,
    drug1_name = drug1_name,
    drug2_name = drug2_name,
    cl_name = cl_name,
    normalization_type = "GR"
  )
  expect_null(result)
})

test_that(".get_combo_curves_colors works as expected", {
  json_path <- system.file(package = "gDRplots", "settings.json")
  s <- gDRutils::get_settings_from_json(json_path = json_path)
  ls_cur_p <- s$COMBO_CURVES_PALETTE

  ls_conc <- factor(c("0.001", "0.005", "0.01", "0.05", "1", "5"))
  res <- .get_combo_curves_colors(ls_conc)
  ls_col <- grDevices::colorRampPalette(ls_cur_p)(2 * NROW(ls_conc))
  expect_length(res, NROW(ls_conc))
  expect_named(res, levels(ls_conc))
  expect_true(all(res %in% ls_col))

  expect_error(.get_combo_curves_colors(c("0.001", "0.01", "1")),
               "Assertion on 'ls_conc_2' failed: Must be of type 'factor'")
})
