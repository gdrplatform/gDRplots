context("Test combo_plots")

test_that("heatmap_combo_metrics works as expected", {
  cl_name <- "cellline_FD"
  drug1_name <- "drug_011"
  drug2_name <- "drug_026"
  
  mae <- gDRutils::get_synthetic_data("combo_matrix")
  se <- mae[[gDRutils::get_supported_experiments("combo")]]
  dt_excess <- gDRutils::convert_se_assay_to_dt(se, "excess")
  dt_isobolograms <- gDRutils::convert_se_assay_to_dt(se, "isobolograms")
  
  plts_1 <- heatmap_combo_metrics(dt_excess, dt_isobolograms,
                                  drug1_name, drug2_name, cl_name)
  expect_is(plts_1, "gg")
  
  normalization_type <- "RV"
  plts_2 <- heatmap_combo_metrics(dt_excess, dt_isobolograms,
                                  drug1_name, drug2_name, cl_name,
                                  normalization_type, as_panel = FALSE)
  expect_is(plts_2, "list")
  expect_equal(names(plts_2), c(names(gDRutils::get_combo_excess_field_names()), "iso_compare"))
  expect_true(all(vapply(seq_along(plts_2), 
                         function(i) grepl(normalization_type, plts_2[[i]]$labels$title), logical(1))))
  
  expect_error(heatmap_combo_metrics(dt_excess = unlist(dt_excess),
                                     dt_isobolograms = dt_isobolograms,
                                     drug1_name = drug1_name,
                                     drug2_name = drug2_name,
                                     cl_name = cl_name),
               "Assertion on 'dt_excess' failed: Must be a data.table")
  expect_error(heatmap_combo_metrics(dt_excess = dt_excess,
                                     dt_isobolograms = dt_isobolograms,
                                     drug1_name = "unknown_drug",
                                     drug2_name = drug2_name,
                                     cl_name = cl_name),
               "Assertion on 'drug1_name' failed: Must be element of set")
  expect_error(heatmap_combo_metrics(dt_excess = dt_excess,
                                     dt_isobolograms = dt_isobolograms,
                                     drug1_name = drug1_name,
                                     drug2_name = drug2_name,
                                     cl_name = cl_name,
                                     normalization_type = "AB"),
               "Assertion on 'normalization_type' failed: Must be element of set")
})


test_that("heatmap_combo_with_isoref works as expected", {
})

test_that("heatmap_combo_with_isoref_qc_panel works as expected", {
  drug1_name <- "drug_001"
  drug2_name <- "drug_026"
  
  mae <- gDRutils::get_synthetic_data("combo_matrix")
  se <- mae[[gDRutils::get_supported_experiments("combo")]]
  dt_excess <- gDRutils::convert_se_assay_to_dt(se, "excess")
  dt_isobolograms <- gDRutils::convert_se_assay_to_dt(se, "isobolograms")
  
  cl_names <- unique(dt_excess[["CellLineName"]])[seq_len(6)]
  
  plt_1 <- heatmap_combo_with_isoref_qc_panel(dt_excess,
                                              dt_isobolograms,
                                              drug1_name, drug2_name,
                                              cl_names)
  expect_is(plt_1, "gg")
  expect_true(grepl("GR", plt_1[["labels"]][["fill"]])) 
  expect_true(grepl(drug1_name, plt_1[["labels"]][["y"]][2]))
  expect_true(grepl(drug2_name, plt_1[["labels"]][["x"]][2]))
  expect_length(plt_1[["layers"]], 2)
  expect_true(grepl(drug1_name, plt_1[["labels"]][["title"]]))
  expect_true(grepl(drug2_name, plt_1[["labels"]][["title"]]))
  expect_equal(plt_1[["labels"]][["fill"]], "Smooth GR")
  expect_equal(plt_1[["labels"]][["colour"]], "iso_level")
  expect_equal(plt_1[["labels"]][["linetype"]], "iso_source")
  
  plt_2 <- heatmap_combo_with_isoref_qc_panel(dt_excess,
                                              dt_isobolograms,
                                              drug1_name, drug2_name,
                                              cl_names = "cellline_XX")
  
  expect_is(plt_2, "gg")
  expect_length(unique(ggplot2::ggplot_build(plt_2)$data[[1]]$PANEL),
                NROW(unique(dt_excess[["CellLineName"]])))
  
  cl_names_NA <- c("cellline_AA", "cellline_XX")
  plt_3 <- heatmap_combo_with_isoref_qc_panel(dt_excess,
                                              dt_isobolograms,
                                              drug1_name, drug2_name,
                                              cl_names = cl_names_NA)
  
  expect_is(plt_3, "gg")
  expect_length(unique(ggplot2::ggplot_build(plt_3)$data[[1]]$PANEL), 
                NROW(intersect(cl_names_NA, unique(dt_excess[["CellLineName"]]))))
  
  normalization_type <- "RV"
  plt_4 <- heatmap_combo_with_isoref_qc_panel(dt_excess,
                                              dt_isobolograms,
                                              drug1_name, drug2_name,
                                              cl_names,
                                              normalization_type = normalization_type,
                                              iso_levels = NULL,
                                              colors_vec = c("#003366", "#FFFAFA", "#FF8800"))
  expect_is(plt_4, "gg")
  expect_true(grepl(normalization_type, plt_4[["labels"]][["fill"]]))
  expect_true(grepl(drug1_name, plt_4[["labels"]][["y"]][2]))
  expect_true(grepl(drug2_name, plt_4[["labels"]][["x"]][2]))
  expect_length(plt_4[["layers"]], 1) # without isoline
  expect_true(grepl(drug1_name, plt_4[["labels"]][["title"]]))
  expect_true(grepl(drug2_name, plt_4[["labels"]][["title"]]))
  expect_equal(plt_4[["labels"]][["fill"]], "Smooth RV")
  expect_equal(plt_4[["labels"]][["colour"]], NULL)
  expect_equal(plt_4[["labels"]][["linetype"]], NULL)
  expect_true(any(grepl("#003366", plt_4[["plot_env"]][["hm_color_palette"]])))
  expect_true(any(grepl("#FF8800", plt_4[["plot_env"]][["hm_color_palette"]])))
  
  expect_error(heatmap_combo_with_isoref_qc_panel(dt_excess = unlist(dt_excess),
                                                  dt_isobolograms = dt_isobolograms,
                                                  drug1_name = drug1_name,
                                                  drug2_name = drug2_name,
                                                  cl_names = cl_names),
               "Assertion on 'dt_excess' failed: Must be a data.table")
  
  expect_error(heatmap_combo_with_isoref_qc_panel(dt_excess = dt_excess,
                                                  dt_isobolograms = unlist(dt_isobolograms),
                                                  drug1_name = drug1_name,
                                                  drug2_name = drug2_name,
                                                  cl_names = cl_names),
               "Assertion on 'dt_isobolograms' failed: Must be a data.table")
  
  expect_error(heatmap_combo_with_isoref_qc_panel(dt_excess = dt_excess,
                                                  dt_isobolograms = dt_isobolograms,
                                                  drug1_name = 1,
                                                  drug2_name = drug2_name,
                                                  cl_names = cl_names),
               "Assertion on 'drug1_name' failed: Must be of type 'string'")
  
  expect_error(heatmap_combo_with_isoref_qc_panel(dt_excess = dt_excess,
                                                  dt_isobolograms = dt_isobolograms,
                                                  drug1_name = drug1_name,
                                                  drug2_name = "drug_XX",
                                                  cl_names = cl_names),
               "Assertion on 'drug2_name' failed: Must be element of set")
  
  expect_error(heatmap_combo_with_isoref_qc_panel(dt_excess = dt_excess,
                                                  dt_isobolograms = dt_isobolograms,
                                                  drug1_name = drug1_name,
                                                  drug2_name = drug2_name,
                                                  cl_names = 1),
               "Assertion on 'cl_names' failed: Must be of type 'character'")

  expect_error(heatmap_combo_with_isoref_qc_panel(dt_excess = dt_excess,
                                                  dt_isobolograms = dt_isobolograms,
                                                  drug1_name = drug1_name,
                                                  drug2_name = drug2_name,
                                                  cl_names = cl_names,
                                                  normalization_type = "XX"),
               "Assertion on 'normalization_type' failed: Must be element of set")
  
  expect_error(heatmap_combo_with_isoref_qc_panel(dt_excess = dt_excess,
                                                  dt_isobolograms = dt_isobolograms,
                                                  drug1_name = drug1_name,
                                                  drug2_name = drug2_name,
                                                  cl_names = cl_names,
                                                  iso_levels = LETTERS[seq_len(6)]),
               "`iso_levels` must be a valid numeric value")
  
  expect_error(heatmap_combo_with_isoref_qc_panel(dt_excess = dt_excess,
                                                  dt_isobolograms = dt_isobolograms,
                                                  drug1_name = drug1_name,
                                                  drug2_name = drug2_name,
                                                  cl_names = cl_names,
                                                  colors_vec = c("pinky", "blackish")),
               "`colors_vec` must be a valid color name")
  
  expect_error(heatmap_combo_with_isoref_qc_panel(dt_excess = dt_excess,
                                                  dt_isobolograms = dt_isobolograms,
                                                  drug1_name = drug1_name,
                                                  drug2_name = drug2_name,
                                                  cl_names = cl_names,
                                                  no_breaks = "10"),
               "Assertion on 'no_breaks' failed: Must be of type 'single integerish value'")
})

test_that("prep_hm_limits works as expected", {
  vec <- c(1.0089, 0.9806, 0.1174, -0.1657, -0.2826)
  # GR smooth
  expect_equal(prep_hm_limits(vec), c(-0.2826, 1.0089))
  # RV smooth
  expect_equal(prep_hm_limits(vec, normalization_type = "RV"), c(0, 1.0089))
  
  # hsa_excess
  vec <- c(-0.1016, -0.0647, 0.0021, 0.0328, 0.6824)
  expect_equal(prep_hm_limits(vec, metric = "hsa_excess"), c(-0.25, 0.6824))
  # bliss_excess
  vec <- c(-0.2651, -0.1289, 0.0051, 0.0202, 0.0394)
  expect_equal(prep_hm_limits(vec, metric = "bliss_excess"), c(-0.2651, 0.25))
  
  expect_error(prep_hm_limits(LETTERS[1:5]),
               "Assertion on 'num_vec' failed: Must be of type 'numeric'")
  expect_error(prep_hm_limits(vec, metric = "str"),
               "Assertion on 'metric' failed: Must be element of set")
  expect_error(prep_hm_limits(vec, normalization_type = "str"),
               "Assertion on 'normalization_type' failed: Must be element of set")
})

test_that("transform_log_conc works as expected", {
  vec <- c(0, 0.003, 0.01, 0.03)
  result <- log10(vec)
  result[1] <- result[2] + (result[2] - result[3])
  expect_equal(transform_log_conc(vec), result)
  expect_equal(transform_log_conc(vec), result)
  
  expect_equal(transform_log_conc(c(0, 0.001)), c(-3.5, -3.0))
  
  expect_error(transform_log_conc(LETTERS[seq_len(5)]), 
               "Assertion on 'conc_vec' failed: Must be of type 'numeric'")
  expect_error(transform_log_conc(-1:2),
               "Assertion on 'conc_vec' failed: Element 1 is not >= 0")
})
