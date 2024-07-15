context("Test drugCombo plot")

test_mae <- gDRutils::get_synthetic_data("finalMAE_combo_matrix_small")
test_exp <- gDRutils::convert_combo_data_to_dt(test_mae[[gDRutils::get_supported_experiments("combo")]][1:2, 1:2])

test_that("get_combo_score_matrix", {
  expect_error(
    get_combo_score_matrix(data = data.table::data.table()),
    "Assertion on 'data' failed: Must be of type 'list'"
  )
  err_msg_1 <- "Assertion on 'normalization_type' failed: Must be element of set"
  err_msg_2 <- "{'RV','GR'}"
  expect_error(
    get_combo_score_matrix(test_exp, normalization_type = "dummy"),
    sprintf("%s %s", err_msg_1, err_msg_2), fixed = TRUE
  )
})

test_that("plolty_drug_combo_heatmap", {
  
  # data
  score_field <- "hsa_score"
  p_name <- gDRutils::get_combo_score_field_names()[score_field]
  metric_growth <- "RV"
  p_metric_growth <- switch(metric_growth,
                            "GR" = "GR value",
                            "RV" = "Relative Viability")
  
  # plot
  pl1 <- plolty_drug_combo_heatmap(test_exp, score_field, metric_growth)
  checkmate::expect_class(pl1, "plotly")
  
  ref_data <- data.table::data.table(get_combo_score_matrix(test_exp, score_field, metric_growth))
  ob_data <- data.table::data.table(pl1$x$data[[pl1$x$data_index]]$z)
  expect_equal(dim(ob_data), dim(ref_data))
  expect_equal(ref_data[c(2, 1), ], ob_data)
  
  expect_false(is.null(pl1$x$data[[3]]$text)) # hovers exist
  expect_true(is.matrix(pl1$x$data[[3]]$text)) # matrix with hovers for heatmam fields
  expect_true(all(grepl(p_name, pl1$x$data[[3]]$text))) # valid hovers for selected metric
  
  exp_title <- sprintf("%s with %s", p_name, p_metric_growth)
  expect_equal(pl1$x$layout$title, exp_title)
  
  # errors
  expect_error(
    plolty_drug_combo_heatmap(list(a = 1)),
    "Assertion on 'data[[1]]' failed: Must be a data.table, not double.",
    fixed = TRUE
  )
  err_msg <- "Assertion on 'metric_combo' failed: Must be element of set "
  err_msg2 <- "{'hsa_score','bliss_score','CIScore_50','CIScore_80'}, but is 'dummy'."
  expect_error(
    plolty_drug_combo_heatmap(test_exp, metric_combo = "dummy", metric_growth = "RV"),
    sprintf("%s%s", err_msg, err_msg2),
    fixed = TRUE
  )
  expect_error(
    plolty_drug_combo_heatmap(test_exp, 
                              metric_combo = names(gDRutils::get_combo_assay_names())[1], 
                              metric_growth = "relval"),
    "Assertion on 'metric_growth' failed: Must be element of set {'GR','RV'}, but is 'relval'.",
    fixed = TRUE
  )
})

test_that("plolty_drug_combo_conc_heatmap", {
  
  # data
  field_type <- "hsa_excess"
  drug1 <- "drug_004"
  drug2 <- "drug_021"
  cl_name <- "cellline_HB"
  i_levels <- c("0.25", "0.5", "0.75")
  
  assay_type <- "excess"
  a_name <- gDRutils::get_combo_excess_field_names()[field_type]
  metric_growth <- "GR"
  
  # plot
  pl1 <- plolty_drug_combo_conc_heatmap(data = test_exp,
                                        drug1,
                                        drug2,
                                        cl_name,
                                        normalization_type = "GR",
                                        c_assay = field_type,
                                        iso_levels = i_levels
  )
  checkmate::expect_class(pl1, "plotly")
  
  ref_data <- get_combo_base_data(test_exp,
                                  drug1,
                                  drug2,
                                  cl_name,
                                  field_type,
                                  metric_growth)
  pidfs <- gDRutils::get_prettified_identifiers()
  ref_data <- get_combo_base_data(test_exp, drug1, drug2, cl_name, field_type, metric_growth)
  ref_mat <- ref_data$matrix
  expect_equal(nrow(ref_mat) * ncol(ref_mat), length(pl1$x$attrs[[1]]$text))
  
  exp_title <- sprintf("%s (%s for %s, T = %sh)",
                       cl_name,
                       a_name,
                       metric_growth,
                       ref_data$condition[[pidfs[["duration"]]]])
  expect_equal(pl1$x$layoutAttrs[[1]]$title, exp_title)
  
  # errors
  expect_error(
    plolty_drug_combo_conc_heatmap(data = data.table::data.table(a = 1), drug1, drug2, cl_name),
    "Assertion on 'data' failed: Must be of type 'list', not 'data.table/data.frame'.",
    fixed = TRUE
  )
  expect_error(
    plolty_drug_combo_conc_heatmap(data = list(a = 1), drug1, drug2, cl_name),
    "Assertion on 'data[[1]]' failed: Must be a data.table, not double.",
    fixed = TRUE
  )
  expect_error(
    plolty_drug_combo_conc_heatmap(data = test_exp, drug1 = 1, drug2, cl_name),
    "Assertion on 'drug1_name' failed: Must be of type 'string', not 'double'.",
    fixed = TRUE
  )
  expect_error(
    plolty_drug_combo_conc_heatmap(data = test_exp, drug1, drug2 = c("01", "02"), cl_name),
    "Assertion on 'drug2_name' failed: Must have length 1.",
    fixed = TRUE
  )
  expect_error(
    plolty_drug_combo_conc_heatmap(data = test_exp, drug1, drug2, cell_line = list(c = "cl2")),
    "Assertion on 'cell_line' failed: Must be of type 'string', not 'list'.",
    fixed = TRUE
  )
  
  expect_error(
    plolty_drug_combo_conc_heatmap(data = test_exp, drug1, drug2, cl_name, c_assay = "dummy"),
    "Assertion on 'c_assay' failed: Must be element of set {'smooth','hsa_excess','bliss_excess'}, but is 'dummy'.",
    fixed = TRUE
  )
  
  expect_error(
    plolty_drug_combo_conc_heatmap(data = test_exp, drug1, drug2, cl_name, normalization_type = "GRvalue"),
    "Assertion on 'normalization_type' failed: Must be element of set {'GR','RV'}, but is 'GRvalue'.",
    fixed = TRUE
  )
  
  expect_error(
    plolty_drug_combo_conc_heatmap(data = test_exp, drug1, drug2, cl_name, iso_levels = c(50, 75)),
    "Assertion on 'iso_levels' failed: Must be of type 'character', not 'double'.",
    fixed = TRUE
  )
  
  expect_error(
    plolty_drug_combo_conc_heatmap(data = test_exp, drug1, drug2, cl_name, iso_levels = c("0.223", "0.33")),
    "Assertion on '0.223' failed. "
  )
  
  test_exp_NA <- test_exp
  test_exp$isobolograms$`Pos x` <- NA
  test_exp$isobolograms$`Pos y` <- NA
  
  
  pl1_NA <- plolty_drug_combo_conc_heatmap(
    data = test_exp_NA,
    drug1,
    drug2,
    cl_name,
    normalization_type = "GR",
    c_assay = "smooth",
    iso_levels = i_levels
  )
  checkmate::expect_class(pl1_NA, "plotly")
  
})


test_that("calc_up_limes", {
  expect_equal(calc_up_limes(1:10),
               data.table::data.table(ux = 1:10, ulx = c(1:9 + 0.5, 10)))
  
  expect_equal(calc_up_limes(c(1, 3, 22, 88), cnames = c("a", "b")),
               data.table::data.table(a = c(1, 3, 22, 88), b = c(2, 12.5, 55, 88)))
  
  expect_error(
    calc_up_limes(x = c(letters[1:10])),
    "Assertion on 'x' failed: Must be of type 'numeric', not 'character'."
  )
  expect_error(
    calc_up_limes(x = 1:10, cnames = 1:2),
    "Assertion on 'cnames' failed: Must be of type 'character', not 'integer'."
  )
})

test_that(
  "map_coords", {
    set.seed(1234)
    ref <- data.table::data.table(Pos_x = rep(0:5, each = 5), Pos_y = rep(0:4, 6))
    data <- data.table::data.table(Pos_x = 4:1 + 0.1, Pos_y = 0:3 + 0.1)
    ref_rand <- ref[sample(nrow(ref)), ]
    expect_equal(map_coords(ref, data), list(data_idx = 1:4, ref_idx = c(21, 17, 13, 9)))
    
    data2 <- data.table::data.table(Pos_x = 4:1 + 0.1, Pos_y = 0:3 + 0.6)
    expect_equal(map_coords(data.table::copy(ref_rand), data2), list(data_idx = 1:4, ref_idx = c(4, 30, 27, 16)))
    
    # expected records with x,y coords in data (outside x,y pos in ref) filetered out
    data3 <- data.table::data.table(Pos_x = 5:0 + 0.1, Pos_y = -1:4 + 0.6)
    expect_equal(map_coords(data.table::copy(ref_rand), data3),
                 list(data_idx = 2:5, ref_idx = c(4, 30, 27, 16)))
    
    expect_error(
      map_coords(list()),
      "Assertion on 'ref' failed: Must be a data.table, not list."
    )
    expect_error(
      map_coords(
        ref = data.table::data.table(a = 1, b = 2),
        "Assertion on 'ref' failed: Must have at least 2 rows, but has 1 rows. "
      )
    )
    expect_error(
      map_coords(ref, 1:2),
      "Assertion on 'data' failed: Must be a data.table, not integer."
    )
    expect_error(
      map_coords(
        ref,
        data = data.table::data.table(a = 1, b = 2),
        "Assertion on 'data' failed: Must have at least 2 rows, but has 1 rows. "
      )
    )
    expect_error(
      map_coords(ref, data, r_cols = 1:2),
      "Assertion on 'r_cols' failed: Must be of type 'character', not 'integer'."
    )
    expect_error(
      map_coords(ref, data, d_cols = 1:2),
      "Assertion on 'd_cols' failed: Must be of type 'character', not 'integer'."
    )
  })

test_that("get_isobologram_data", {
  expect_error(
    get_isobologram_data(data = data.table::data.table()),
    "Assertion on 'data' failed: Must be of type 'list'"
  )
  drug1 <- "drug_004"
  drug2 <- "drug_021"
  cl_name <- "cellline_HB"
  assay_type <- "smooth"
  metric_growth <- "GR"
  sdata <- get_combo_base_data(test_exp,
                               drug1,
                               drug2,
                               cl_name,
                               assay_type,
                               metric_growth)
  i_levels <- c("0.25", "0.5", "0.75")
  
  err_msg_1 <- "Assertion on 'normalization_type' failed: Must be element of set"
  err_msg_2 <- "{'GR','RV'}, but is 'dummy'."
  expect_error(
    get_isobologram_data(
      test_exp,
      drug1,
      drug2,
      cl_name,
      c_assay = names(gDRutils::get_combo_assay_names(group = "combo_iso"))[1],
      normalization_type = "dummy",
      drug1_axis = sdata$drug1_axis,
      drug2_axis = sdata$drug2_axis,
      iso_levels = i_levels
    ),
    sprintf("%s %s", err_msg_1, err_msg_2), fixed = TRUE
  )
})

test_that("plolty_drug_combo_ratio", {
  assay_type <- "excess"
  a_name <- gDRutils::get_assay_names(type = assay_type, prettify = TRUE)
  metric_growth <- "GR"
  drug1 <- "drug_004"
  drug2 <- "drug_021"
  cl_name <- "cellline_HB"
  i_levels <- c("0.25", "0.5", "0.75")
  normalization_type <- "GR"
  assay_name <- "smooth"
  
  pl1 <- plolty_drug_combo_ratio(
    data = test_exp,
    drug1_name = drug1,
    drug2_name = drug2,
    cell_line = cl_name,
    normalization_type = normalization_type,
    c_assay = assay_name,
    iso_levels = i_levels
  )
  
  checkmate::expect_class(pl1, "plotly")
  
  sdata <- get_combo_base_data(test_exp,
                               drug1,
                               drug2,
                               cl_name,
                               assay_name,
                               metric_growth)
  
  ref_data <- get_isobologram_data(
    test_exp,
    drug1,
    drug2,
    cl_name,
    c_assay = names(gDRutils::get_combo_assay_names(group = "combo_iso"))[1],
    normalization_type = normalization_type,
    drug1_axis = sdata$drug1_axis,
    drug2_axis = sdata$drug2_axis,
    iso_levels = i_levels
  )
  
  expect_equal(length(pl1$x$attrs[[1]]$text), NROW(ref_data))
  
  # errors
  expect_error(
    plolty_drug_combo_ratio(data = data.table::data.table(a = 1), drug1, drug2, cl_name),
    "Assertion on 'data' failed: Must be of type 'list', not 'data.table/data.frame'.",
    fixed = TRUE
  )
  expect_error(
    plolty_drug_combo_ratio(data = list(a = 1), drug1, drug2, cl_name),
    "Assertion on 'data[[1]]' failed: Must be a data.table, not double.",
    fixed = TRUE
  )
  expect_error(
    plolty_drug_combo_ratio(data = test_exp, drug1 = 1, drug2, cl_name),
    "Assertion on 'drug1_name' failed: Must be of type 'string', not 'double'.",
    fixed = TRUE
  )
  expect_error(
    plolty_drug_combo_ratio(data = test_exp, drug1, drug2 = c("01", "02"), cl_name),
    "Assertion on 'drug2_name' failed: Must have length 1.",
    fixed = TRUE
  )
  expect_error(
    plolty_drug_combo_ratio(data = test_exp, drug1, drug2, cell_line = list(c = "cl2")),
    "Assertion on 'cell_line' failed: Must be of type 'string', not 'list'.",
    fixed = TRUE
  )
  
  expect_error(
    plolty_drug_combo_ratio(data = test_exp, drug1, drug2, cl_name, c_assay = "dummy"),
    "Assertion on 'c_assay' failed: Must be element of set {'smooth','hsa_excess','bliss_excess'}, but is 'dummy'.",
    fixed = TRUE
  )
  
  expect_error(
    plolty_drug_combo_ratio(data = test_exp, drug1, drug2, cl_name, normalization_type = "GRvalue"),
    "Assertion on 'normalization_type' failed: Must be element of set {'GR','RV'}, but is 'GRvalue'.",
    fixed = TRUE
  )
  
  expect_error(
    plolty_drug_combo_ratio(data = test_exp, drug1, drug2, cl_name, iso_levels = c(50, 75)),
    "Assertion on 'iso_levels' failed: Must be of type 'character', not 'double'.",
    fixed = TRUE
  )
  
  expect_error(
    plolty_drug_combo_ratio(data = test_exp, drug1, drug2, cl_name, iso_levels = c("0.223", "0.33")),
    "Assertion on '0.223' failed. "
  )
  
})

test_that(".round.conc works as expected", {
  conc <- c("0", "0.00000179", "0.00000541", "0.0000163", "0.0000493",
            "0.000149", "0.00045", "0.00136", "0.0041")
  expect_equal(length(unique(conc)), length(unique(.round_conc(conc))))
  
  # old solution
  expect_false(length(unique(conc)) == length(unique(round(as.numeric(conc), 5))))
})

test_that("get_iso_colors", {
  ### expected values
  gic <- get_iso_colors()
  expect_true(length(gic) > 2)
  expect_identical("character", class(gic))
  gic2 <- get_iso_colors(formals(get_iso_colors)[[1]][[3]])
  expect_true(any(gic != gic2))
  expect_identical(length(gic), length(gic2))
  
  ### errors
  expect_error(get_iso_colors("inv_param"), "'arg' should be one of ")
})

test_that("assert_RGB_format", {
  color_vector <- c(25, 56, 189)
  expect_equal(assert_RGB_format(color_vector), NULL)
  color_vector <- c(201, 128, 352)
  expect_error(assert_RGB_format(color_vector), 
               "Some value is greater than 255. Not valid RGB format.")
})

test_that("get_combo_col_settings",  {
  ### expected values
  gcan <- names(gDRutils::get_combo_assay_names()[1])
  gcc <-
    get_combo_col_settings(g_metric = "GR", assay_type = gcan)
  expect_true(inherits(gcc, "list"))
  expect_identical(sort(names(gcc)), c("breaks", "colors", "limits"))
  
  ### errors
  err_msg <- "Assertion on 'assay_type' failed: "
  expect_error(get_combo_col_settings("GR", 8), err_msg)
  err_msg <- "Assertion on 'g_metric' failed: "
  expect_error(get_combo_col_settings("grvalue", 8), err_msg)
})

