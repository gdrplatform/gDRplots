context("Test metricClustering plot")

test_that("build_label_clustering works as expected", {
})

test_that("prepareDataMH works as expected", {
  
  test_mae <- gDRutils::get_synthetic_data("finalMAE_combo_matrix_small")
  test_exp <- gDRutils::convert_combo_data_to_dt(test_mae[[gDRutils::get_supported_experiments("combo")]]
                                                 [1:2, 1:2])[[1]]
  
  expected_length <-
    length(unique(test_exp$Concentration)) * length(unique(test_exp$`r Id`)) * length(unique(test_exp$`c Id`))
  
  expect_error(prepareDataMH(NULL), "Assertion on 'data' failed: Must be a data.table, not 'NULL'.")
  expect_error(prepareDataMH(test_exp), "argument \"variable\" is missing, with no default")
  expect_error(prepareDataMH(test_exp, NULL), "Assertion on 'variable' failed: Must be of type 'string'")
  res <- prepareDataMH(test_exp, "Smooth")
  expect_length(res, 2)
  expect_equal(names(res), c("data_matrix", "annotations"))
  expect_length(res$data_matrix, expected_length)
  expect_length(res$annotations, 2)
  expect_equal(names(res$annotations), c("row_annotation", "col_annotation"))
  expect_equal(unique(unlist(res$annotations$col_annotation[, "Tissue"])), unique(test_exp$Tissue))
  
  ic50 <- c(0.1, 10, 1, 100)
  dt <- data.table::data.table(
    `Drug Name` = c(1, 1, 2, 2),
    `Cell Line Name` = c("a", "b", "a", "b"),
    `Drug MOA` = ("moa1"),
    IC50 = ic50
  )
  res <- prepareDataMH(dt, "IC50")
  expect_equal(length(res), 2)
  expect_equal(c(t(res$data_matrix)), log10(ic50))
  
})

test_that("plotly_metric_clustering works as expected", {
})

test_that("map_annotations works as expected", {
})

test_that(".tidy_cell_metadata works as expected", {
})

test_that("create_formula works as expected", {
})

test_that("calc_duplicate_freq works as expected", {
})


test_that("convert_factor_to_character works as expected", {
  dt <- data.table::data.table(a = LETTERS, b = as.factor(LETTERS))
  
  expect_equal(unname(unlist(lapply(dt, class))), c("character", "factor"))
  obs <- convert_factor_to_character(dt)
  expect_equal(dim(obs), dim(dt))
  expect_equal(names(obs), names(dt))
  expect_equal(class(obs), class(dt))
  expect_equal(unname(unlist(lapply(obs, class))), c("character", "character"))
  
  dt_2 <- data.table::data.table(
    a = LETTERS[1:5], b = factor(LETTERS[1:5], levels = LETTERS)
  )
  obs_2 <- convert_factor_to_character(dt_2)
  expect_equal(dim(obs_2), dim(dt_2))
  expect_equal(names(obs_2), names(dt_2))
  expect_equal(unname(unlist(lapply(obs_2, class))), c("character", "character"))
  
  expect_error(convert_factor_to_character(as.list(dt)))
})
