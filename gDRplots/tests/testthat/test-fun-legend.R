context("Test fun legend")

test_that("get_legend_title works as expected", {
  expect_equal(get_legend_title("Drug 2"), list(text = "<b>Drug 2</b>"))
  expect_equal(get_legend_title("Drug", has_codrug_data = TRUE),
               list(text = "<b>Concentration 2</b> && <b>Drug</b>"))
  expect_equal(get_legend_title("Concentration 2", has_codrug_data = TRUE),
               list(text = "<b>Concentration 2</b>"))
  expect_equal(get_legend_title("Concentration 2"),
               list(text = "<b>Concentration 2</b>"))
  expect_equal(get_legend_title("Conc_2", default_var = "Conc_2"),
               list(text = "<b>Conc_2</b>"))
  expect_equal(get_legend_title("Conc 2", has_codrug_data = TRUE, default_var = "Conc 2"),
               list(text = "<b>Conc 2</b>"))
  expect_equal(get_legend_title("Drug 2", has_codrug_data = TRUE, default_var = "Conc_2"),
               list(text = "<b>Conc_2</b> && <b>Drug 2</b>"))
  expect_equal(get_legend_title("Drug 2", has_codrug_data = TRUE, default_var = "Drug 2"),
               list(text = "<b>Drug 2</b>"))
  expect_equal(get_legend_title("Drug 2", default_var = "Conc_2"),
               list(text = "<b>Drug 2</b>"))
  expect_equal(get_legend_title("None", has_codrug_data = FALSE, default_var = "Conc_2"),
               NULL)
  expect_equal(get_legend_title("None", has_codrug_data = TRUE, default_var = "Conc_2"),
               list(text = "<b>Conc_2</b>"))
  expect_equal(get_legend_title(var = NULL), NULL)
  expect_equal(get_legend_title(var = NULL, has_codrug_data = TRUE),
               list(text = "<b>Concentration 2</b>"))
  
  expect_error(
    get_legend_title(1),
    "Assertion on 'var' failed: Must be of type 'string' \\(or 'NULL'\\), not 'double'."
  )
  expect_error(
    get_legend_title("test", has_codrug_data = 1),
    "Assertion on 'has_codrug_data' failed: Must be of type 'logical flag', not 'double'."
  )
  expect_error(
    get_legend_title("test", default_var = c("d", "b")),
    "Assertion on 'default_var' failed: Must have length 1."
  )
  expect_error(
    get_legend_title("test", default_var = NULL),
    "Assertion on 'default_var' failed: Must be of type 'string', not 'NULL'."
  )
})


test_that("do_show_legend works as expected", {
  data <-
    data.table::data.table(
      x = seq_len(10),
      y = seq_len(10) + 0.5,
      `Tissue` = c("lung", "brain"),
      `Drug MOA` = c("DHFR", "DNA"),
      `Concentration 2` = 0:1,
      check.names = FALSE
    )
  
  expect_true(do_show_legend("none", data))
  expect_true(do_show_legend("Primary Tissue", data))
  
  data[["Concentration 2"]] <- 2
  expect_false(do_show_legend("none", data))
  expect_true(do_show_legend("Tissue", data))
  
  expect_true(do_show_legend("Drug MOA", data))
  
  data[["Concentration 2"]] <- NULL
  expect_false(do_show_legend("none", data))
  
  data[["Drug MOA"]] <- "single value"
  expect_true(do_show_legend("Drug MOA", data))
  
  expect_error(
    do_show_legend(NULL, data),
    "Assertion on 'var_col' failed: Must be of type 'string', not 'NULL'."
  )
  expect_error(
    do_show_legend(c("none", "str"), data),
    "Assertion on 'var_col' failed: Must have length 1."
  )
  expect_error(
    do_show_legend("none", list(a = 2)),
    "Assertion on 'data' failed: Must be a data.table, not list."
  )
})
