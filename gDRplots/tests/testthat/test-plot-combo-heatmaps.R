context("Test combo_plots")

test_that("heatmap_combo_metrics works as expected", {
  cl_name <- "cellline_FD"
  drug1_name <- "drug_011"
  drug2_name <- "drug_026"
  
  mae <- gDRutils::get_synthetic_data("combo_matrix")
  se <- mae[[gDRutils::get_supported_experiments("combo")]]
  dt_excess <- gDRutils::convert_se_assay_to_dt(se, "excess")
  dt_isobolograms <- gDRutils::convert_se_assay_to_dt(se, "isobolograms")
  
  dt_excess_sub <- 
    dt_excess[DrugName == drug1_name & DrugName_2 == drug2_name & CellLineName == cl_name, ]
  
  plt_1 <- heatmap_combo_metrics(dt_excess, 
                                 dt_isobolograms,
                                 drug1_name, drug2_name, 
                                 cl_name) # default
  expect_is(plt_1, "gg")
  expect_true(grepl("GR", plt_1[["labels"]][["title"]]))
  expect_true(grepl("GR", plt_1[["labels"]][["fill"]]))
  expect_true(grepl("Smooth", plt_1[["labels"]][["title"]]))
  expect_true(grepl("Smooth", plt_1[["labels"]][["fill"]]))
  expect_true(any(grepl(drug1_name, plt_1[["labels"]][["y"]])))
  expect_true(any(grepl(drug2_name, plt_1[["labels"]][["x"]])))
  expect_length(plt_1[["layers"]], 2) # heatmap + isobolograms # nolint
  expect_length(unique(ggplot2::ggplot_build(plt_1)[["data"]][[2]][["colour"]]), 3) # iso_levels
  default_lim <- round(prep_hm_limits(dt_excess_sub[normalization_type == "GR", ][["smooth"]]), 0)
  expect_equal(range(as.numeric(ggplot2::get_guide_data(plt_1, "fill")[[".label"]])), default_lim)
  
  norm_type <- "RV"
  plt_2 <- heatmap_combo_metrics(dt_excess, 
                                 dt_isobolograms = NULL,
                                 drug1_name, drug2_name, 
                                 cl_name,
                                 metric = "bliss_excess",
                                 normalization_type = norm_type,
                                 show_values = TRUE)
  expect_is(plt_2, "gg")
  expect_true(grepl(norm_type, plt_2[["labels"]][["title"]]))
  expect_true(grepl(norm_type, plt_2[["labels"]][["fill"]]))
  expect_true(grepl("Bliss Excess", plt_2[["labels"]][["title"]]))
  expect_true(grepl("Bliss Excess", plt_2[["labels"]][["fill"]]))
  expect_length(plt_2[["layers"]], 2) # heatmap + numbers # nolint
  expect_length(ggplot2::ggplot_build(plt_2)[["data"]], 2)
  expect_equal(
    NROW(ggplot2::ggplot_build(plt_2)[["data"]][[2]][["label"]]),
    NROW(dt_excess[DrugName == drug1_name & DrugName_2 == drug2_name & CellLineName == cl_name &
                     normalization_type == norm_type]$bliss_excess))
  
  ls_col_1 <- c("#008B8B", "#FFA500", "#FFFFFF")
  plt_3 <- heatmap_combo_metrics(dt_excess, 
                                 dt_isobolograms,
                                 drug1_name, drug2_name, 
                                 cl_name,
                                 metric = "hsa_excess",
                                 colors_vec_excess = ls_col_1,
                                 limit = c(NA, NA),
                                 iso_levels = NULL)
  expect_is(plt_3, "gg")
  expect_true(grepl("HSA Excess", plt_3[["labels"]][["title"]]))
  expect_true(grepl("HSA Excess", plt_3[["labels"]][["fill"]]))
  expect_length(plt_3[["layers"]], 1) # heatmap # nolint
  expect_length(ggplot2::ggplot_build(plt_3)[["data"]], 1)
  expect_true(any(ls_col_1 %in% unique(ggplot2::ggplot_build(plt_3)[["data"]][[1]][["fill"]])))
  lim_fill <- round(range(dt_excess_sub[normalization_type == "GR", ] [["hsa_excess"]], na.rm = TRUE), 1)
  expect_equal(range(as.numeric(ggplot2::get_guide_data(plt_3, "fill")[[".label"]])), lim_fill)
  
  plt_4 <- heatmap_combo_metrics(dt_excess, 
                                 dt_isobolograms,
                                 drug1_name, drug2_name, 
                                 cl_name,
                                 iso_levels = c("0.2", "0.4"),
                                 colors_vec_excess = ls_col_1)
  expect_false(any(ls_col_1 %in% unique(ggplot2::ggplot_build(plt_4)[["data"]][[1]][["fill"]])))
  expect_true(any(.get_smooth_palette(50) %in% unique(ggplot2::ggplot_build(plt_4)[["data"]][[1]][["fill"]])))
  expect_equal(NROW(unique(ggplot2::ggplot_build(plt_4)[["data"]][[2]][["colour"]])),
               NROW(c("0.2", "0.4")))
  expect_equal(plt_4[["labels"]][["colour"]], "Iso Levels")
  
  plt_5 <- heatmap_combo_metrics(dt_excess, 
                                 dt_isobolograms,
                                 drug1_name, drug2_name, 
                                 cl_name,
                                 metric = "bliss_excess",
                                 colors_vec_smooth = ls_col_1,
                                 limit = c(-0.1, 0.1))
  expect_false(any(ls_col_1 %in% unique(ggplot2::ggplot_build(plt_5)[["data"]][[1]][["fill"]])))
  expect_true(any(.get_excess_palette(50) %in% unique(ggplot2::ggplot_build(plt_5)[["data"]][[1]][["fill"]])))
  expect_equal(range(as.numeric(ggplot2::get_guide_data(plt_5, "fill")[[".label"]])), c(-0.1, 0.1))
  
  plt_1_swap_axes <- heatmap_combo_metrics(dt_excess, 
                                           dt_isobolograms,
                                           drug1_name, drug2_name, 
                                           cl_name,
                                           swap_axes = TRUE)
  expect_true(any(grepl(drug1_name, plt_1_swap_axes[["labels"]][["x"]])))
  expect_true(any(grepl(drug2_name, plt_1_swap_axes[["labels"]][["y"]])))
  expect_equal(plt_1_swap_axes[["data"]][["pos_y"]], plt_1[["data"]][["pos_x"]])
  expect_equal(plt_1_swap_axes[["data"]][["pos_x"]], plt_1[["data"]][["pos_y"]])
  expect_equal(plt_1_swap_axes[["data"]][["smooth"]], plt_1[["data"]][["smooth"]])
  
  # lack of data
  plt_6 <- heatmap_combo_metrics(dt_excess[normalization_type == "RV"], 
                                 dt_isobolograms,
                                 drug1_name, drug2_name, 
                                 cl_name)
  expect_is(plt_6, "gg")
  expect_true(any(grepl(drug1_name, plt_6[["labels"]][["y"]])))
  expect_true(any(grepl(drug2_name, plt_6[["labels"]][["x"]])))
  expect_length(ggplot2::ggplot_build(plt_6)[["data"]][[1]], 0)
})

test_that("heatmap_combo_metrics_panel works as expected", {
  cl_name <- "cellline_FD"
  drug1_name <- "drug_011"
  drug2_name <- "drug_026"
  
  mae <- gDRutils::get_synthetic_data("combo_matrix")
  se <- mae[[gDRutils::get_supported_experiments("combo")]]
  dt_excess <- gDRutils::convert_se_assay_to_dt(se, "excess")
  dt_isobolograms <- gDRutils::convert_se_assay_to_dt(se, "isobolograms")
  
  plts_1 <- heatmap_combo_metrics_panel(dt_excess, 
                                        dt_isobolograms,
                                        drug1_name, 
                                        drug2_name, 
                                        cl_name) # default
  expect_is(plts_1, "gg")
  expect_length(ggplot2::ggplot_build(plts_1)$data, 1)
  
  # scenario: list with default isolines
  normalization_type <- "RV"
  plts_2 <- heatmap_combo_metrics_panel(dt_excess, 
                                        dt_isobolograms,
                                        drug1_name, 
                                        drug2_name, 
                                        cl_name,
                                        normalization_type,
                                        show_values = TRUE,
                                        as_list = TRUE)
  expect_is(plts_2, "list")
  expect_equal(names(plts_2), c(names(gDRutils::get_combo_excess_field_names()), "iso_compare"))
  expect_true(all(vapply(
    seq_along(plts_2), 
    function(i) grepl(normalization_type, plts_2[[i]]$labels$title), logical(1))))
  expect_true(all(vapply(
    seq_along(plts_2), 
    function(i) grepl("Iso Levels", plts_2[[i]]$labels$linetype), logical(1))))
  expect_true(all(vapply(
    names(gDRutils::get_combo_excess_field_names()),
    function(i) "label" %in% names(ggplot2::ggplot_build(plts_2[[i]])[["data"]][[2]]), logical(1))))
  
  # scenario: list with defined isolines
  iso_lvl <- c("0.2", "0.4")
  plts_3 <- heatmap_combo_metrics_panel(dt_excess, 
                                        dt_isobolograms,
                                        drug1_name, 
                                        drug2_name, 
                                        cl_name,
                                        iso_levels = iso_lvl, 
                                        as_list = TRUE)
  expect_is(plts_3, "list")
  expect_equal(names(plts_3), c(names(gDRutils::get_combo_excess_field_names()), "iso_compare"))
  expect_equal(NROW(unique(ggplot2::ggplot_build(plts_3[["smooth"]])[["data"]][[2]][["colour"]])),
               NROW(iso_lvl))
  expect_equal(NROW(unique(ggplot2::ggplot_build(plts_3[["iso_compare"]])[["data"]][[3]][["colour"]])),
               NROW(iso_lvl))
  expect_equal(NROW(unique(ggplot2::ggplot_build(plts_3[["hsa_excess"]])[["data"]][[2]][["colour"]])),
               NROW(iso_lvl))
  expect_equal(NROW(unique(ggplot2::ggplot_build(plts_3[["bliss_excess"]])[["data"]][[2]][["colour"]])),
               NROW(iso_lvl))
  
  # scenario: list without isolines (iso_levels = NULL) and defined color palette for smooth heatmap
  ls_col_1 <- c("#008B8B", "#FFA500", "#FFFFFF")
  plts_4 <- heatmap_combo_metrics_panel(dt_excess, 
                                        dt_isobolograms,
                                        drug1_name, 
                                        drug2_name, 
                                        cl_name,
                                        iso_levels = NULL,
                                        colors_vec_smooth = ls_col_1,
                                        as_list = TRUE)
  expect_is(plts_4, "list")
  expect_true(all(vapply(seq_along(plts_4), function(x) is(plts_4[[x]], "gg"), logical(1))))
  expect_length(plts_4, 3) # smooth, hsa_excess, bliss_excess
  expect_equal(names(plts_4), c(names(gDRutils::get_combo_excess_field_names())))
  expect_true(all(ls_col_1 %in% plts_4[["smooth"]][["plot_env"]][["hm_color_palette"]]))
  
  # scenario: list without isolines (t_isobolograms = NULL) and defined color palette for excess heatmaps
  ls_col_2 <- c("#FF1493", "#350040")
  plts_5 <- heatmap_combo_metrics_panel(dt_excess, 
                                        dt_isobolograms = NULL,
                                        drug1_name, 
                                        drug2_name, 
                                        cl_name,
                                        colors_vec_excess = ls_col_2,
                                        as_list = TRUE)
  expect_is(plts_5, "list")
  expect_true(all(vapply(seq_along(plts_5), function(x) is(plts_5[[x]], "gg"), logical(1))))
  expect_length(plts_5, 3) # smooth, hsa_excess, bliss_excess
  expect_true(all(ls_col_2 %in% plts_5[["hsa_excess"]][["plot_env"]][["hm_color_palette"]]))
  
  # scenario: lack data for selection
  dt_excess_ <- 
    data.table::copy(dt_excess)[DrugName == drug1_name & DrugName_2 == drug2_name &
                                  CellLineName == cl_name][1:2, ]
  iso_cols <- gDRutils::get_header("iso_position")
  dt_isobolograms_ <- 
    data.table::copy(dt_isobolograms)[DrugName == drug1_name & DrugName_2 == drug2_name &
                                        CellLineName == cl_name][, (iso_cols) := NA]
  plts_6 <- heatmap_combo_metrics_panel(dt_excess = dt_excess_, 
                                        dt_isobolograms = dt_isobolograms_,
                                        drug1_name, 
                                        drug2_name, 
                                        cl_name,
                                        as_list = TRUE)
  expect_is(plts_6, "list")
  expect_true(all(vapply(seq_along(plts_6), function(x) is(plts_6[[x]], "gg"), logical(1))))
  expect_length(plts_6, 3) # smooth, hsa_excess, bliss_excess
  expect_true(all(vapply(names(plts_6), 
                         function(nm) NROW(ggplot2::ggplot_build(plts_6[[nm]])[["data"]][[1]]) == 0,
                         logical(1)))) # no excess data
  
  # scenario: unlikely excess data to be missing data, but isodata exists
  plts_7 <- heatmap_combo_metrics_panel(dt_excess = dt_excess_, 
                                        dt_isobolograms,
                                        drug1_name, 
                                        drug2_name, 
                                        cl_name,
                                        as_list = TRUE)
  expect_is(plts_7, "list")
  expect_true(all(vapply(seq_along(plts_7), function(x) is(plts_7[[x]], "gg"), logical(1))))
  expect_true(all(vapply(names(gDRutils::get_combo_excess_field_names()), 
                         function(nm) NROW(ggplot2::ggplot_build(plts_7[[nm]])[["data"]][[1]]) == 0,
                         logical(1)))) # no excess data
  
  # scenario: as panel but 2x2
  plts_8 <- heatmap_combo_metrics_panel(dt_excess = dt_excess, 
                                        dt_isobolograms,
                                        drug1_name, 
                                        drug2_name, 
                                        cl_name,
                                        one_row_panel = FALSE)
  expect_is(plts_8, "gg")
  expect_length(ggplot2::ggplot_build(plts_8)$data, 1)
  
  # scenario: as panel but 3x1
  plts_9 <- heatmap_combo_metrics_panel(dt_excess = dt_excess, 
                                        dt_isobolograms,
                                        drug1_name, 
                                        drug2_name, 
                                        cl_name,
                                        one_row_panel = TRUE)
  expect_is(plts_9, "gg")
  expect_length(ggplot2::ggplot_build(plts_9)$data, 3)
  
  # scenario: as list even if one_row_panel = TRUE
  plts_10 <- heatmap_combo_metrics_panel(dt_excess = dt_excess, 
                                         dt_isobolograms,
                                         drug1_name, 
                                         drug2_name, 
                                         cl_name,
                                         iso_levels = c("0.4", "0.6", "0.8"),
                                         as_list = TRUE,
                                         one_row_panel = TRUE)
  expect_is(plts_10, "list")
  expect_equal(names(plts_10), c(names(gDRutils::get_combo_excess_field_names()), "iso_compare"))
  
  # checking assertions
  expect_error(heatmap_combo_metrics_panel(dt_excess = unlist(dt_excess),
                                           dt_isobolograms = dt_isobolograms,
                                           drug1_name = drug1_name,
                                           drug2_name = drug2_name,
                                           cl_name = cl_name),
               "Assertion on 'dt_excess' failed: Must be a data.table")
  expect_error(heatmap_combo_metrics_panel(dt_excess = dt_excess,
                                           dt_isobolograms = dt_isobolograms,
                                           drug1_name = "unknown_drug",
                                           drug2_name = drug2_name,
                                           cl_name = cl_name),
               "Assertion on 'drug1_name' failed: Must be element of set")
  expect_error(heatmap_combo_metrics_panel(dt_excess = dt_excess,
                                           dt_isobolograms = dt_isobolograms,
                                           drug1_name = drug1_name,
                                           drug2_name = drug2_name,
                                           cl_name = cl_name,
                                           normalization_type = "AB"),
               "Assertion on 'normalization_type' failed: Must be element of set")
  expect_error(heatmap_combo_metrics_panel(dt_excess = dt_excess,
                                           dt_isobolograms = dt_isobolograms,
                                           drug1_name = drug1_name,
                                           drug2_name = drug2_name,
                                           cl_name = cl_name,
                                           as_list = "TRUE"),
               "Assertion on 'as_list' failed: Must be of type 'logical flag'")
  expect_error(heatmap_combo_metrics_panel(dt_excess = dt_excess,
                                           dt_isobolograms = dt_isobolograms,
                                           drug1_name = drug1_name,
                                           drug2_name = drug2_name,
                                           cl_name = cl_name,
                                           one_row_panel = 0),
               "Assertion on 'one_row_panel' failed: Must be of type 'logical flag'")
  expect_error(heatmap_combo_metrics_panel(dt_excess = dt_excess,
                                           dt_isobolograms = dt_isobolograms,
                                           drug1_name = drug1_name,
                                           drug2_name = drug2_name,
                                           cl_name = cl_name,
                                           swap_axes = "yes"),
               "Assertion on 'swap_axes' failed: Must be of type 'logical flag'")
  expect_error(heatmap_combo_metrics_panel(dt_excess = dt_excess,
                                           dt_isobolograms = dt_isobolograms,
                                           drug1_name = drug1_name,
                                           drug2_name = drug2_name,
                                           cl_name = cl_name,
                                           show_values = "yes"),
               "Assertion on 'show_values' failed: Must be of type 'logical flag'")
  
  # scenarios: heatmap_combo_metrics_panel works as expected when dt_isobolograms is NULL
  # test with dt_isobolograms as NULL
  plts_iso_1 <- heatmap_combo_metrics_panel(dt_excess, 
                                            iso_levels = NULL,
                                            dt_isobolograms = NULL,
                                            drug1_name, 
                                            drug2_name, 
                                            cl_name)
  
  # check if the output is a ggplot object
  expect_is(plts_iso_1, "gg")
  
  normalization_type <- "RV"
  plts_iso_2 <- heatmap_combo_metrics_panel(dt_excess, 
                                            iso_levels = NULL,
                                            dt_isobolograms = NULL,
                                            drug1_name, 
                                            drug2_name, 
                                            cl_name, 
                                            as_list = TRUE)
  
  # check if the output is a list
  expect_is(plts_iso_2, "list")
  expect_length(plts_iso_2, 3)
  # check if names of the list are as expected
  expect_equal(names(plts_iso_2), names(gDRutils::get_combo_excess_field_names()))
  
  # check if switch_axes works as expected
  plts_iso_2_swap_axes <- heatmap_combo_metrics_panel(dt_excess,
                                                      iso_levels = NULL,
                                                      dt_isobolograms = NULL,
                                                      drug1_name,
                                                      drug2_name,
                                                      cl_name,
                                                      as_list = TRUE,
                                                      swap_axes = TRUE)
  expect_equal(plts_iso_2_swap_axes[["smooth"]][["labels"]][["x"]], plts_iso_2[["smooth"]][["labels"]][["y"]])
  expect_equal(plts_iso_2_swap_axes[["smooth"]][["labels"]][["y"]], plts_iso_2[["smooth"]][["labels"]][["x"]])
  expect_equal(plts_iso_2_swap_axes[["smooth"]][["data"]][["pos_y"]], plts_iso_2[["smooth"]][["data"]][["pos_x"]])
  expect_equal(plts_iso_2_swap_axes[["smooth"]][["data"]][["pos_x"]], plts_iso_2[["smooth"]][["data"]][["pos_y"]])
  expect_equal(plts_iso_2_swap_axes[["smooth"]][["data"]][["smooth"]], plts_iso_2[["smooth"]][["data"]][["smooth"]])
})

test_that("plot_combination_index works as expected", {
  cl_name <- "cellline_BC"
  drug1_name <- "drug_001"
  drug2_name <- "drug_026"
  
  mae <- gDRutils::get_synthetic_data("combo_matrix")
  se <- mae[[gDRutils::get_supported_experiments("combo")]]
  dt_excess <- gDRutils::convert_se_assay_to_dt(se, "excess")
  dt_isobolograms <- gDRutils::convert_se_assay_to_dt(se, "isobolograms")
  
  plt_1 <- plot_combination_index(dt_excess,
                                  dt_isobolograms,
                                  drug1_name, drug2_name,
                                  cl_name) # default
  expect_is(plt_1, "gg")
  expect_length(plt_1[["layers"]], 3)
  expect_true(grepl(drug1_name, plt_1[["labels"]][["x"]]))
  expect_true(grepl(drug2_name, plt_1[["labels"]][["x"]]))
  expect_true(grepl("T=", plt_1[["labels"]][["title"]]))
  expect_true(grepl("GR", plt_1[["labels"]][["title"]]))
  
  plt_2 <- plot_combination_index(dt_excess = NULL,
                                  dt_isobolograms,
                                  drug1_name, drug2_name,
                                  cl_name,
                                  normalization_type = "RV", 
                                  iso_levels = c("0.25", "0.75"),
                                  colors_vec_iso = c("#00008B", "#008B8B"))
  expect_is(plt_2, "gg")
  expect_length(plt_2[["layers"]], 3)
  expect_false(grepl("T=", plt_2[["labels"]][["title"]]))
  expect_true(grepl("RV", plt_2[["labels"]][["title"]]))
  expect_true(
    all(c("#00008B", "#008B8B") %in% unique(ggplot2::ggplot_build(plt_2)[["data"]][[3]][["colour"]])))
  
  plt_3 <- plot_combination_index(dt_excess,
                                  dt_isobolograms,
                                  drug1_name, drug2_name,
                                  cl_name,
                                  iso_levels = c("0.85", "0.9")) # not avialable
  expect_is(plt_3, "gg")
  expect_length(plt_3[["layers"]], 2)
  expect_true(grepl("T=", plt_3[["labels"]][["title"]]))
  expect_true(grepl("GR", plt_3[["labels"]][["title"]]))
  expect_false("colour" %in% names(plt_3[["labels"]])) # lack of `iso_levels`
  
  plt_4 <- plot_combination_index(dt_excess,
                                  dt_isobolograms,
                                  drug1_name, drug2_name,
                                  cl_name,
                                  normalization_type = "RV",
                                  iso_levels = c("0.2", "0.4"))
  expect_is(plt_4, "gg")
  expect_true(grepl("RV", plt_4[["labels"]][["title"]]))
  expect_equal(NROW(unique(ggplot2::ggplot_build(plt_4)[["data"]][[3]][["colour"]])),
               NROW(c("0.2", "0.4")))
  
  expect_error(plot_combination_index(dt_excess = unlist(dt_excess),
                                      dt_isobolograms = dt_isobolograms,
                                      drug1_name = drug1_name,
                                      drug2_name = drug2_name,
                                      cl_name = cl_name),
               "Assertion on 'dt_excess' failed: Must be a data.table")
  expect_error(plot_combination_index(dt_excess = dt_excess,
                                      dt_isobolograms = unlist(dt_isobolograms),
                                      drug1_name = drug1_name,
                                      drug2_name = drug2_name,
                                      cl_name = cl_name),
               "Assertion on 'dt_isobolograms' failed: Must be a data.table")
  expect_error(plot_combination_index(dt_excess = dt_excess,
                                      dt_isobolograms = dt_isobolograms,
                                      drug1_name = "unknown_drug",
                                      drug2_name = drug2_name,
                                      cl_name = cl_name),
               "Assertion on 'drug1_name' failed: Must be element of set")
  expect_error(plot_combination_index(dt_excess = dt_excess,
                                      dt_isobolograms = dt_isobolograms,
                                      drug1_name = drug1_name,
                                      drug2_name = drug2_name,
                                      cl_name = 1),
               "Assertion on 'cl_name' failed: Must be of type 'string'")
  expect_error(plot_combination_index(dt_excess = dt_excess,
                                      dt_isobolograms = dt_isobolograms,
                                      drug1_name = drug1_name,
                                      drug2_name = drug2_name,
                                      cl_name = cl_name,
                                      normalization_type = "AB"),
               "Assertion on 'normalization_type' failed: Must be element of set")
  expect_error(plot_combination_index(dt_excess = dt_excess,
                                      dt_isobolograms = dt_isobolograms,
                                      drug1_name = drug1_name,
                                      drug2_name = drug2_name,
                                      cl_name = cl_name,
                                      iso_levels = NULL),
               "Assertion on 'iso_levels' failed: Must be of type 'character', not 'NULL'")
})

test_that("heatmap_combo_with_isoref works as expected", {
  cl_name <- "cellline_BC"
  drug1_name <- "drug_001"
  drug2_name <- "drug_026"
  
  mae <- gDRutils::get_synthetic_data("combo_matrix")
  se <- mae[[gDRutils::get_supported_experiments("combo")]]
  dt_excess <- gDRutils::convert_se_assay_to_dt(se, "excess")
  dt_isobolograms <- gDRutils::convert_se_assay_to_dt(se, "isobolograms")
  
  plt_1 <- heatmap_combo_with_isoref(dt_excess,
                                     dt_isobolograms,
                                     drug1_name, 
                                     drug2_name,
                                     cl_name)
  expect_is(plt_1, "gg")
  expect_true(grepl("GR", plt_1[["labels"]][["fill"]])) 
  expect_true(grepl(cl_name, plt_1[["labels"]][["title"]])) 
  expect_true(any(grepl(drug1_name, plt_1[["labels"]][["y"]])))
  expect_true(any(grepl(drug2_name, plt_1[["labels"]][["x"]])))
  
  normalization_type <- "RV"
  iso_lvl <- c("0.2", "0.4")
  ls_col <- c("#003366", "#FFFAFA", "#FF8800")
  plt_2 <- heatmap_combo_with_isoref(dt_excess,
                                     dt_isobolograms,
                                     drug1_name, 
                                     drug2_name,
                                     cl_name,
                                     normalization_type = normalization_type,
                                     iso_levels = iso_lvl,
                                     colors_vec = ls_col)
  expect_is(plt_2, "gg")
  expect_true(grepl(normalization_type, plt_2[["labels"]][["fill"]]))
  expect_true(all(ls_col %in% plt_2[["plot_env"]][["hm_color_palette"]]))
  expect_equal(NROW(unique(ggplot2::ggplot_build(plt_2)[["data"]][[2]][["colour"]])),
               NROW(iso_lvl))
  
  plt_3 <- heatmap_combo_with_isoref(dt_excess,
                                     dt_isobolograms,
                                     drug1_name, 
                                     drug2_name,
                                     cl_name,
                                     normalization_type = normalization_type,
                                     iso_levels = NULL,
                                     colors_vec = ls_col)
  expect_is(plt_3, "gg")
  expect_true(grepl(normalization_type, plt_3[["labels"]][["fill"]]))
  expect_true(all(ls_col %in% plt_3[["plot_env"]][["hm_color_palette"]]))
  expect_length(ggplot2::ggplot_build(plt_3)[["data"]], 1) # no isoline data
  
  dt_excess_ <- dt_excess[DrugName == drug1_name & DrugName_2 == drug2_name & 
                            CellLineName == cl_name][1:2, ]
  plt_4 <- heatmap_combo_with_isoref(dt_excess_,
                                     dt_isobolograms,
                                     drug1_name, 
                                     drug2_name,
                                     cl_name)
  expect_is(plt_4, "gg")
  expect_true(grepl(cl_name, plt_4[["labels"]][["title"]])) 
  expect_true(grepl(drug1_name, plt_4[["labels"]][["y"]]))
  expect_true(grepl(drug2_name, plt_4[["labels"]][["x"]]))
  expect_length(ggplot2::ggplot_build(plt_4)[["data"]][[1]], 0) # no excess data
  
  plt_5 <- heatmap_combo_with_isoref(dt_excess,
                                     dt_isobolograms,
                                     drug1_name, 
                                     drug2_name,
                                     cl_name, 
                                     metric = "hsa_excess")
  expect_is(plt_5, "gg")
  expect_true(grepl(drug1_name, plt_5[["labels"]][["y"]]))
  expect_true(grepl(drug2_name, plt_5[["labels"]][["x"]]))
  expect_true(grepl("HSA Excess GR", plt_5[["labels"]][["fill"]]))
  expect_true(grepl("GR", plt_5[["labels"]][["linetype"]]))
  expect_length(names(plt_5[["guides"]][["guides"]]), NROW(c("fill", "linetype", "colour")))
  
  plt_6 <- heatmap_combo_with_isoref(dt_excess,
                                     dt_isobolograms,
                                     drug1_name, 
                                     drug2_name,
                                     cl_name, 
                                     metric = "hsa_excess",
                                     iso_levels = c("0.25", "0.75"),
                                     swap_axes = TRUE)
  expect_is(plt_6, "gg")
  expect_true(grepl(drug1_name, plt_6[["labels"]][["x"]]))
  expect_true(grepl(drug2_name, plt_6[["labels"]][["y"]]))
  expect_true(grepl("HSA Excess GR", plt_6[["labels"]][["fill"]]))
  expect_true(grepl("Iso Levels", plt_6[["labels"]][["colour"]]))
  expect_length(names(plt_6[["guides"]][["guides"]]), NROW(c("fill", "linetype", "colour")))
  
  plt_1_swap_axes <- heatmap_combo_with_isoref(dt_excess,
                                               dt_isobolograms,
                                               drug1_name, 
                                               drug2_name,
                                               cl_name,
                                               swap_axes = TRUE)
  
  expect_equal(plt_1_swap_axes[["labels"]][["x"]], plt_1[["labels"]][["y"]])
  expect_equal(plt_1_swap_axes[["labels"]][["y"]], plt_1[["labels"]][["x"]])
  expect_equal(plt_1_swap_axes[["data"]][["pos_y"]], plt_1[["data"]][["pos_x"]])
  expect_equal(plt_1_swap_axes[["data"]][["pos_x"]], plt_1[["data"]][["pos_y"]])
  expect_equal(plt_1_swap_axes[["data"]][["smooth"]], plt_1[["data"]][["smooth"]])
})

test_that("heatmap_combo_with_isoref_panel works as expected", {
  cl_names <- 
    c("cellline_AA", "cellline_EA", "cellline_IB", "cellline_MC", "cellline_BC", "cellline_FD")
  drug1_name <- "drug_001"
  drug2_name <- "drug_026"
  
  mae <- gDRutils::get_synthetic_data("combo_matrix")
  se <- mae[[gDRutils::get_supported_experiments("combo")]]
  dt_excess <- gDRutils::convert_se_assay_to_dt(se, "excess")
  dt_isobolograms <- gDRutils::convert_se_assay_to_dt(se, "isobolograms")
  
  plt_1 <- heatmap_combo_with_isoref_panel(dt_excess,
                                           dt_isobolograms,
                                           drug1_name, 
                                           drug2_name,
                                           cl_names)
  expect_is(plt_1, "gg")
  expect_true(grepl(drug1_name, plt_1[["labels"]][["y"]]))
  expect_true(grepl(drug2_name, plt_1[["labels"]][["x"]]))
  expect_length(plt_1[["layers"]], 2)
  expect_true(grepl(drug1_name, plt_1[["labels"]][["title"]]))
  expect_true(grepl(drug2_name, plt_1[["labels"]][["title"]]))
  expect_equal(plt_1[["labels"]][["fill"]], "Smooth GR")
  expect_equal(plt_1[["labels"]][["colour"]], "Iso Levels")
  expect_equal(plt_1[["labels"]][["linetype"]], "GR")
  
  plt_2 <- heatmap_combo_with_isoref_panel(dt_excess,
                                           dt_isobolograms,
                                           drug1_name, 
                                           drug2_name,
                                           cl_names = "cellline_XX")
  
  expect_is(plt_2, "gg")
  expect_length(unique(ggplot2::ggplot_build(plt_2)$data[[1]]$PANEL),
                NROW(unique(dt_excess[["CellLineName"]])))
  
  cl_names_NA <- c("cellline_AA", "cellline_XX")
  plt_3 <- heatmap_combo_with_isoref_panel(dt_excess,
                                           dt_isobolograms,
                                           drug1_name, 
                                           drug2_name,
                                           cl_names = cl_names_NA)
  
  expect_is(plt_3, "gg")
  expect_length(unique(ggplot2::ggplot_build(plt_3)$data[[1]]$PANEL), 
                NROW(intersect(cl_names_NA, unique(dt_excess[["CellLineName"]]))))
  
  normalization_type <- "RV"
  plt_4 <- heatmap_combo_with_isoref_panel(dt_excess,
                                           dt_isobolograms,
                                           drug1_name, 
                                           drug2_name,
                                           cl_names,
                                           normalization_type = normalization_type,
                                           iso_levels = NULL,
                                           colors_vec = c("#003366", "#FFFAFA", "#FF8800"))
  expect_is(plt_4, "gg")
  expect_true(grepl(drug1_name, plt_4[["labels"]][["y"]]))
  expect_true(grepl(drug2_name, plt_4[["labels"]][["x"]]))
  expect_length(plt_4[["layers"]], 1) # without isoline
  expect_true(grepl(drug1_name, plt_4[["labels"]][["title"]]))
  expect_true(grepl(drug2_name, plt_4[["labels"]][["title"]]))
  expect_equal(plt_4[["labels"]][["fill"]], "Smooth RV")
  expect_equal(plt_4[["labels"]][["colour"]], NULL)
  expect_equal(plt_4[["labels"]][["linetype"]], NULL)
  expect_true(any(grepl("#003366", plt_4[["plot_env"]][["hm_color_palette"]])))
  expect_true(any(grepl("#FF8800", plt_4[["plot_env"]][["hm_color_palette"]])))
  
  plt_5 <- heatmap_combo_with_isoref_panel(dt_excess,
                                           dt_isobolograms,
                                           drug1_name, 
                                           drug2_name,
                                           cl_names,
                                           metric = "hsa_excess",
                                           swap_axes = TRUE)
  expect_is(plt_5, "gg")
  expect_true(grepl(drug2_name, plt_5[["labels"]][["y"]]))
  expect_true(grepl(drug1_name, plt_5[["labels"]][["x"]]))
  expect_length(plt_5[["layers"]], 2)
  expect_true(grepl(drug1_name, plt_5[["labels"]][["title"]]))
  expect_true(grepl(drug2_name, plt_5[["labels"]][["title"]]))
  expect_equal(plt_5[["labels"]][["fill"]], "HSA Excess GR")
  expect_equal(plt_5[["labels"]][["colour"]], "Iso Levels")
  expect_equal(plt_5[["labels"]][["linetype"]], "GR")
  
  plt_1_swap_axes <- heatmap_combo_with_isoref_panel(dt_excess,
                                                     dt_isobolograms,
                                                     drug1_name, 
                                                     drug2_name,
                                                     cl_names,
                                                     swap_axes = TRUE)
  
  expect_equal(plt_1_swap_axes[["labels"]][["x"]], plt_1[["labels"]][["y"]])
  expect_equal(plt_1_swap_axes[["labels"]][["y"]], plt_1[["labels"]][["x"]])
  expect_equal(plt_1_swap_axes[["data"]][["pos_y"]], plt_1[["data"]][["pos_x"]])
  expect_equal(plt_1_swap_axes[["data"]][["pos_x"]], plt_1[["data"]][["pos_y"]])
  expect_equal(plt_1_swap_axes[["data"]][["smooth"]], plt_1[["data"]][["smooth"]])
  
  expect_error(heatmap_combo_with_isoref_panel(dt_excess = unlist(dt_excess),
                                               dt_isobolograms = dt_isobolograms,
                                               drug1_name = drug1_name,
                                               drug2_name = drug2_name,
                                               cl_names = cl_names),
               "Assertion on 'dt_excess' failed: Must be a data.table")
  
  expect_error(heatmap_combo_with_isoref_panel(dt_excess = dt_excess,
                                               dt_isobolograms = unlist(dt_isobolograms),
                                               drug1_name = drug1_name,
                                               drug2_name = drug2_name,
                                               cl_names = cl_names),
               "Assertion on 'dt_isobolograms' failed: Must be a data.table")
  
  expect_error(heatmap_combo_with_isoref_panel(dt_excess = dt_excess,
                                               dt_isobolograms = dt_isobolograms,
                                               drug1_name = 1,
                                               drug2_name = drug2_name,
                                               cl_names = cl_names),
               "Assertion on 'drug1_name' failed: Must be of type 'string'")
  
  expect_error(heatmap_combo_with_isoref_panel(dt_excess = dt_excess,
                                               dt_isobolograms = dt_isobolograms,
                                               drug1_name = drug1_name,
                                               drug2_name = "drug_XX",
                                               cl_names = cl_names),
               "Assertion on 'drug2_name' failed: Must be element of set")
  
  expect_error(heatmap_combo_with_isoref_panel(dt_excess = dt_excess,
                                               dt_isobolograms = dt_isobolograms,
                                               drug1_name = drug1_name,
                                               drug2_name = drug2_name,
                                               cl_names = 1),
               "Assertion on 'cl_names' failed: Must be of type 'character'")
  
  expect_error(heatmap_combo_with_isoref_panel(dt_excess = dt_excess,
                                               dt_isobolograms = dt_isobolograms,
                                               drug1_name = drug1_name,
                                               drug2_name = drug2_name,
                                               cl_names = cl_names,
                                               normalization_type = "XX"),
               "Assertion on 'normalization_type' failed: Must be element of set")
  
  expect_error(heatmap_combo_with_isoref_panel(dt_excess = dt_excess,
                                               dt_isobolograms = dt_isobolograms,
                                               drug1_name = drug1_name,
                                               drug2_name = drug2_name,
                                               cl_names = cl_names,
                                               iso_levels = LETTERS[seq_len(6)]),
               "`iso_levels` must be a valid numeric value")
  
  expect_error(heatmap_combo_with_isoref_panel(dt_excess = dt_excess,
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
  
  # GR smooth symmetric
  expect_equal(prep_hm_limits(vec, symmetric = TRUE), c(-1.0089, 1.0089))
  # RV smooth symmetric
  expect_equal(prep_hm_limits(vec, normalization_type = "RV", symmetric = TRUE), c(-1.0089, 1.0089))
  
  # hsa_excess
  vec <- c(-0.1016, -0.0647, 0.0021, 0.0328, 0.6824)
  expect_equal(prep_hm_limits(vec, metric = "hsa_excess"), c(-0.25, 0.6824))
  # hsa_excess symmetric
  expect_equal(prep_hm_limits(vec, metric = "hsa_excess", symmetric = TRUE), c(-0.6824, 0.6824))
  
  # bliss_excess
  vec <- c(-0.2651, -0.1289, 0.0051, 0.0202, 0.0394)
  expect_equal(prep_hm_limits(vec, metric = "bliss_excess"), c(-0.2651, 0.25))
  # bliss_excess symmetric
  expect_equal(prep_hm_limits(vec, metric = "bliss_excess", symmetric = TRUE), c(-0.2651, 0.2651))
  
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

test_that(".get_tile_size works as expected", {
  res_1 <- .get_tile_size(c(-0.35, 0, 0.349, 0.7))
  expect_equal(res_1, 0.351)
  
  res_2 <- .get_tile_size(c(0.349, 0.7, -0.35, 0))
  expect_equal(res_1, res_2)
  
  res_3 <- .get_tile_size(c(0.349, 0.7))
  expect_equal(res_3, 0.351)
  
  res_4 <- .get_tile_size(c(0.349))
  expect_equal(res_4, 0.5)
  
  res_5 <- .get_tile_size(c(0.349, 0.349))
  expect_equal(res_5, 0.5)
  
  res_6 <- .get_tile_size(c(0.349, 0.7, -0.35, 0, -0.35, 0))
  expect_equal(res_6, 0.351)
  
  res_7 <- .get_tile_size(c(0.5, 0.75, NA))
  expect_equal(res_7, 0.25)
  expect_error(.get_tile_size(c("0.5", "0.75")),
               "Assertion on 'pos_vec' failed: Must be of type 'numeric'")
})

test_that(".get_iso_colors works as expected", {
  json_path <- system.file(package = "gDRplots", "settings.json")
  s <- gDRutils::get_settings_from_json(json_path = json_path)
  ls_iso_p <- s$ISOLINE_PALETTE
  
  ls_iso_lvl <- c("0.25", "0.5", "0.75")
  res <- .get_iso_colors(ls_iso_lvl)
  ls_col <- grDevices::colorRampPalette(ls_iso_p)(2 * NROW(ls_iso_lvl))
  expect_length(res, NROW(ls_iso_lvl))
  expect_named(res, ls_iso_lvl)
  expect_true(all(res %in% ls_col))
  
  expect_error(.get_iso_colors(c(0.5, 0.75)),
               "Assertion on 'iso_levels' failed: Must be of type 'character'")
})

test_that(".get_smooth_palette works as expected", {
  json_path <- system.file(package = "gDRplots", "settings.json")
  s <- gDRutils::get_settings_from_json(json_path = json_path)
  ls_smooth <- s$SMOOTH_PALETTE
  
  no_br <- 25
  res <- .get_smooth_palette(no_br)
  ls_col <- grDevices::colorRampPalette(ls_smooth)(no_br)
  expect_length(res, NROW(ls_col))
  expect_true(all(res %in% ls_col))
  
  expect_error(.get_smooth_palette("all"),
               "Assertion on 'no_breaks' failed: Must be of type 'single integerish value'")
  expect_error(.get_smooth_palette(0),
               "Assertion on 'no_breaks' failed: Element 1 is not >= 2.")
})

test_that(".get_excess_palette works as expected", {
  json_path <- system.file(package = "gDRplots", "settings.json")
  s <- gDRutils::get_settings_from_json(json_path = json_path)
  ls_excess <- s$EXCESS_PALETTE
  
  no_br <- 2
  res <- .get_excess_palette(no_br)
  ls_col <- grDevices::colorRampPalette(ls_excess)(no_br)
  expect_length(res, NROW(ls_col))
  expect_true(all(res %in% ls_col))
  
  expect_error(.get_excess_palette("all"),
               "Assertion on 'no_breaks' failed: Must be of type 'single integerish value'")
  expect_error(.get_excess_palette(1),
               "Assertion on 'no_breaks' failed: Element 1 is not >= 2.")
})
