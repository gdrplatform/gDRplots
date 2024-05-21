context("Test fun internals")

test_that("coerce_cotreatment_data works as expected", {
  cotrt_data <- data.table::data.table(
    `Drug Name` = rep(sprintf("drug_00%s", 1:5), 2),
    `Drug Name 2` = c(rep(NA, 5), sprintf("drug_00%s", 5:9)),
    Concentration = 1e-3 * 1:10
  )
  
  conc_2_vals <-  1e-3 * 11:20
  res_cotrt_data <- cotrt_data
  res_cotrt_data$`Concentration 2` <- conc_2_vals
  
  expect_equal(coerce_cotreatment_data(cotrt_data), cotrt_data)
  
  cotrt_data$`Concentration 2` <- factor(conc_2_vals)
  expect_equal(coerce_cotreatment_data(cotrt_data), res_cotrt_data)
  
  cotrt_data$`Concentration 2` <- as.character(conc_2_vals)
  expect_equal(coerce_cotreatment_data(cotrt_data), res_cotrt_data)
  
  cotrt_data$`Concentration 2` <- factor(as.character(conc_2_vals))
  expect_equal(coerce_cotreatment_data(cotrt_data), res_cotrt_data)
  
  cotrt_data$`Concentration 2` <- factor(letters[1:10])
  res_cotrt_data$`Concentration 2` <- as.numeric(NA)
  res <- expect_warning(coerce_cotreatment_data(cotrt_data))
  expect_equal(res, res_cotrt_data)
  
  expect_error(coerce_cotreatment_data(as.list(cotrt_data)),
                         "Assertion on 'data' failed: Must be a data.table")
  expect_error(coerce_cotreatment_data(NULL),
                         "Assertion on 'data' failed: Must be a data.table")
})


test_that("reformat_untreated_cases works as expected", {
  expect_equal(
    reformat_untreated_cases("Cell Line Name: d\nDrug Name: a\n(untreated at 0 &mu;M)\nTitle: 1.00"),
    "Cell Line Name: d\nDrug Name: a\nTitle: 1.00"
  )
  expect_equal(
    reformat_untreated_cases(c(
      "Cell Line Name: d\nDrug Name: a\n(untreated at 0 &mu;M)\nTitle: 1.00",
      "Cell Line Name: d\nDrug Name: a\n(untreated at 4 &mu;M)\nTitle: 1.00"
    )),
    c("Cell Line Name: d\nDrug Name: a\nTitle: 1.00",
      "Cell Line Name: d\nDrug Name: a\n(untreated at 4 &mu;M)\nTitle: 1.00")
  )
  
  expect_error(
    reformat_untreated_cases(1L), 
    "Assertion on 'label' failed: Must inherit from class 'character'/'numeric'")
  expect_error(
    reformat_untreated_cases(NULL), 
    "Assertion on 'label' failed: Must inherit from class 'character'/'numeric'")
})


test_that("replaceValues works as expected", {
  vals <- data.table::data.table(value = runif(10, 1, 2),
                                 group = rep(c("A", "B"), each = 10))
  nts <- data.table::data.table(value = runif(2, 0, 0.1),
                                group = rep("nt", 2))
  dt <- rbind(vals, nts)
  dtr <- replaceValues(x = dt, column = "group", replacee = "nt", replacement = c("A", "B"))
  
  # replacee is absent
  expect_false("nt" %in% unique(dtr$group))
  # rows are not missing
  expect_true(all(dt$value %in% dtr$value))
  # number of observations per group checks out
  expect_identical(tapply(dtr$value, dtr$group, length),
                   tapply(vals$value, vals$group, length) + nrow(nts))
  expect_identical(length(which(dtr$value == nts$value)),
                   nrow(nts) * length(unique(vals$group)))
  for (i in nts$value) {
    expect_identical(sum(i == dtr), 2L)
  }
  
  dt$group <- factor(dt$group)
  dtr_2 <- replaceValues(x = dt, column = "group", replacee = "nt", replacement = c("A", "B"))
  
  # replacee is absent
  expect_false("nt" %in% unique(dtr_2$group))
  # rows are not missing
  expect_true(all(dt$value %in% dtr_2$value))
  # replaced level is absent in data
  expect_false("nt" %in% levels(dtr_2$group))
  expect_equal(levels(dtr_2$group), setdiff(levels(dt$group), "nt"))
  
  
  expect_error(
    replaceValues(x = as.list(dt), column = "group", replacee = "nt", replacement = c("A", "B")),
    "Assertion on 'x' failed: Must be a data.table")
  expect_error(
    replaceValues(x = dt, column = 1L, replacee = "nt", replacement = c("A", "B")),
    "Assertion on 'column' failed: Must be of type 'string'")
  expect_error(
    replaceValues(x = dt, column = "str", replacee = "nt", replacement = c("A", "B")),
    "Assertion on 'column' failed: Must be element of set")
  expect_error(
    replaceValues(x = dt, column = "group", replacee = 1L, replacement = c("A", "B")),
    "Assertion on 'replacee' failed: Must be of type 'string'")
  expect_error(
    replaceValues(x = dt, column = "group", replacee = "str", replacement = c("A", "B")),
    "Assertion on 'replacee' failed: Must be element of set")
  expect_error(
    replaceValues(x = dt, column = "group", replacee = "nt", replacement = 1L),
    "Assertion on 'replacement' failed: Must be of type 'character'")
})
