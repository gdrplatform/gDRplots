context("Test fun utils")

test_that("isColDark works as expected", {
  color_names <- c("#33cc33", "#d6f5d6", "#0a290a")
  expect_false(isColDark(color_names[1]))
  expect_false(isColDark(color_names[2]))
  expect_true(isColDark(color_names[3]))
  
  color_names <- c("darkblue", "yellow", "tomato")
  expect_true(isColDark(color_names[1]))
  expect_false(isColDark(color_names[2]))
  expect_false(isColDark(color_names[3]))
  
  expect_error(isColDark(1), "Must be of type 'string'")
  expect_error(isColDark("nice pink"), "Must be valid color name")
})


test_that("getColLuminance works as expected", {
  color_names <- c("#33cc33", "#d6f5d6", "#0a290a")
  expect_equal(getColLuminance(color_names[1]), 0.429871, tolerance = 1e-5)
  expect_equal(getColLuminance(color_names[2]), 0.839746, tolerance = 1e-5)
  expect_equal(getColLuminance(color_names[3]), 0.016340, tolerance = 1e-5)
  
  color_names <- c("darkblue", "yellow", "tomato")
  expect_equal(getColLuminance(color_names[1]), 0.0186408, tolerance = 1e-5)
  expect_equal(getColLuminance(color_names[2]), 0.9278000, tolerance = 1e-5)
  expect_equal(getColLuminance(color_names[3]), 0.3238907, tolerance = 1e-5)
  
  expect_error(getColLuminance(1), "Must be of type 'string'")
  expect_error(getColLuminance("nice pink"), "Must be valid color name")
})


test_that("isValidColor works as expected", {
  color_names <- c("#33cc33", "#d6f5d6", "#0a290a", "#F9B42DFF", "#714D6932", "#C2F970DC")
  expect_true(all(vapply(color_names, isValidColor, logical(1))))
  
  color_names <- c("darkblue", "yellow", "tomato")
  expect_true(all(vapply(color_names, isValidColor, logical(1))))
  
  color_names <- c("nice pink", "RED", "#0a290", "#C2F970D")
  expect_false(all(vapply(color_names, isValidColor, logical(1))))
  
  expect_error(isValidColor(1), "Must be of type 'string'")
  expect_error(isValidColor(NULL), "Must be of type 'string'")
  expect_error(isValidColor(NA), "Assertion on 'col_name' failed: May not be NA.")
})

test_that("colorToHex works as expected", {
  color_names <- c("orange", "darkblue", "green", "lavenderblush", "gray66", "slategray2", "tomato")
  expect_identical(
    vapply(color_names, colorToHex, character(1), USE.NAMES = FALSE),
    c("#FFA500", "#00008B", "#00FF00", "#FFF0F5", "#A8A8A8", "#B9D3EE", "#FF6347"))
  
  expect_error(colorToHex(1), "Must be of type 'string'")
  expect_error(colorToHex(NULL), "Must be of type 'string'")
  expect_error(colorToHex("pinki"), "Must be valid color name")
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


test_that("buildLabel works as expected", {
  dt1 <-
    data.table::data.table(
      "Drug Name" = letters[seq_len(3)],
      "Concentration" = seq_len(3),
      "Drug Name 2" = "untreated",
      "Concentration 2" = 4:6,
      "Drug Name 3" = "untreated",
      "Concentration 3" = 7:9,
      "GR AOC" = seq_len(3),
      "Cell Line Name" = letters[4:6]
    )
  dt2 <-
    data.table::data.table(
      "Drug Name" = letters[seq_len(3)],
      "Concentration" = seq_len(3),
      "Drug Name 2" = "untreated",
      "GR AOC" = seq_len(3),
      "Cell Line Name" = letters[4:6]
    )
  
  var_x <- "Cell Line Name"
  var_y <- "GR AOC"
  var_col <- "none"
  title_x <- "Title"
  
  output <- c("Cell")
  
  expect_equal(buildLabel(dt1, "distribution")[1],
               "Cell Line Name: d\nDrug Name: a\n(untreated at 4 &mu;M)\nTitle: 1.00")
  expect_equal(buildLabel(dt2, "distribution")[1], "Cell Line Name: d\nDrug Name: a\nTitle: 1.00")
  
  var_col <- "Cell Line Name"
  var_not_col <- "Drug Name"
  
  expect_equal(buildLabel(dt1, "distribution")[1],
               "Cell Line Name: d\nDrug Name: a\nCell Line Name: d\n(untreated at 4 &mu;M)\nTitle: 1.00")
  
  dt3 <-
    data.table::data.table(
      "var_not_col" = letters[seq_len(3)],
      "Concentration" = seq_len(3),
      "Drug Name 2" = "untreated",
      "Concentration 2" = 4:6,
      "Drug Name 3" = "untreated",
      "Concentration 3" = 7:9,
      "var_y" = seq_len(3),
      "var_col" = letters[4:6]
    )
  dt4 <-
    data.table::data.table(
      "var_not_col" = letters[seq_len(3)],
      "Concentration" = seq_len(3),
      "Drug Name 2" = "untreated",
      "var_y" = seq_len(3),
      "var_col" = letters[4:6]
    )
  
  expect_equal(
    buildLabel(dt3, "curve")[1],
    "Cell Line Name: d\nDrug Name: a\nConcentration: 1 &mu;M\n(untreated at 4 &mu;M)\nGR AOC: 1.00")
  expect_equal(buildLabel(dt4, "curve")[1], 
               "Cell Line Name: d\nDrug Name: a\nConcentration: 1 &mu;M\nGR AOC: 1.00")
  
  expect_equal(buildLabel(dt1, "grid")[1], "Cell Line Name: d\nDrug Name: a")
  
  dt4 <-
    data.table::data.table(
      "Cell Line" = "Cell Line A",
      "Drug A" = 2,
      "Drug B" = 3,
      "Drug Name 2" = "untreated",
      "Concentration 2" = 4:6
    )
  dt5 <-
    data.table::data.table(
      "Cell Line" = "Cell Line A",
      "Drug A" = 2,
      "Drug B" = 3
    )
  var_x <- "Drug A"
  var_y <- "Drug B"
  var_txt <- "Cell Line"
  
  expect_equal(buildLabel(dt4, "contrast")[1],
               "Cell Line: Cell Line A\n(untreated at 4 &mu;M)\nDrug A: 2\nDrug B: 3")
  expect_equal(buildLabel(dt5, "contrast"), "Cell Line: Cell Line A\nDrug A: 2\nDrug B: 3")
  
  dt6 <-
    data.table::data.table(
      "Drug Name" = letters[seq_len(3)],
      "Concentration" = seq_len(3),
      "Drug Name 2" = "untreated",
      "Concentration 2" = 4:6,
      "Drug Name 3" = "untreated",
      "Concentration 3" = 7:9,
      "GR AOC" = seq_len(3),
      "Cell Line Name" = letters[4:6],
      "Drug MOA" = paste0("moa_", LETTERS[seq_len(3)]),
      "Tissue" = paste0("tissue_", LETTERS[seq_len(3)])
    )
  dt7 <-
    data.table::data.table(
      "Drug Name" = letters[seq_len(3)],
      "Concentration" = seq_len(3),
      "Drug Name 2" = "untreated",
      "GR AOC" = seq_len(3),
      "Cell Line Name" = letters[4:6],
      "Drug MOA" = paste0("moa_", LETTERS[seq_len(3)]),
      "Tissue" = paste0("tissue_", LETTERS[seq_len(3)])
    )
  var_col <- "none"
  var_grp <- "none"
  var_x <- "Drug Name"
  var_y <- "GR AOC"
  expect_equal(buildLabel(dt6, "ranking")[1],
               "Drug Name: a\nDrug MOA: moa_A\nTissue: tissue_A\n(untreated at 4 &mu;M)\nTitle: 1.00")
  expect_equal(buildLabel(dt7, "ranking")[1],
               "Drug Name: a\nDrug MOA: moa_A\nTissue: tissue_A\nTitle: 1.00")
  var_x <- "Cell Line Name"
  expect_equal(buildLabel(dt6, "ranking")[1],
               "Cell Line Name: d\nTissue: tissue_A\nDrug MOA: moa_A\n(untreated at 4 &mu;M)\nTitle: 1.00")
  var_x <- "Concentration"
  expect_error(buildLabel(dt6, "ranking"), "bad value provided for 'var_col'")
  var_x <- "Drug Name"
  var_col <- "Tissue"
  var_grp <- "Tissue"
  expect_equal(buildLabel(dt6, "ranking")[1],
               "Drug Name: a\nTissue: tissue_A\nDrug MOA: moa_A\n(untreated at 4 &mu;M)\nTitle: 1.00")
  var_col <- "Concentration" 
  var_grp <- "Concentration" 
  expect_error(buildLabel(dt6, "ranking"), "bad value provided for 'var_col'")
  var_col <- "Drug MOA"
  var_grp <- "none"
  expect_equal(buildLabel(dt6, "ranking")[1],
               "Drug Name: a\nDrug MOA: moa_A\nTissue: tissue_A\n(untreated at 4 &mu;M)\nTitle: 1.00")
  var_col <- "Tissue"
  expect_equal(buildLabel(dt6, "ranking")[1],
               "Drug Name: a\nTissue: tissue_A\nDrug MOA: moa_A\n(untreated at 4 &mu;M)\nTitle: 1.00")
  var_col <- "none"
  var_grp <- "Drug MOA"
  expect_equal(buildLabel(dt6, "ranking")[1],
               "Drug Name: a\nTissue: tissue_A\nDrug MOA: moa_A\n(untreated at 4 &mu;M)\nTitle: 1.00")
  var_grp <- "Tissue"
  expect_equal(buildLabel(dt6, "ranking")[1],
               "Drug Name: a\nDrug MOA: moa_A\nTissue: tissue_A\n(untreated at 4 &mu;M)\nTitle: 1.00")
  var_grp <- "Concentration"
  expect_error(buildLabel(dt6, "ranking"), "bad value provided for 'var_grp'")
  
  var_grp <- "none"
  var_col <- "Drug MOA"
  expect_equal(buildLabel(dt6, "ranking")[1],
               "Drug Name: a\nDrug MOA: moa_A\nTissue: tissue_A\n(untreated at 4 &mu;M)\nTitle: 1.00")
  var_col <- "Tissue"
  expect_equal(buildLabel(dt6, "ranking")[1],
               "Drug Name: a\nTissue: tissue_A\nDrug MOA: moa_A\n(untreated at 4 &mu;M)\nTitle: 1.00")
  var_col <- "Concentration"
  expect_error(buildLabel(dt6, "ranking"), "bad value provided for 'var_grp'")
  
  dt8 <-
    data.table::data.table(
      "value" = 1,
      "Pos_x" = -3.2,
      "Pos_y" = -2,
      "Pos_x_Ref" = -3,
      "Pos_y_Ref" = -1.8,
      "x2_off" = 1.2,
      "name" = "GR75",
      "Log10_Ratio_Conc" = -1.5,
      "Log2_CI" = 0,
      "Iso_Level" = 0.75,
      "level" = 0.25
    )
  matrix_pretty <- "MX full"
  condition <- list(
    "Drug Name" = "a",
    "Drug Name 2" = "b",
    "Gnumber" = "G1",
    "Gnumber 2" = "G2",
    "Duration" = "72",
    "Cell Line Name" = "cellline_A"
  )
  
  var_col <- "none"
  var_grp <- "none"
  var_x <- "Drug Name"
  expect_equal(buildLabel(dt8, "combo1-heatmap"), 
               "Cell Line: cellline_A\nb: 0.00063\na: 0.01\nMX full: 1")
  expect_equal(buildLabel(dt8, "combo1-lines_ref"), 
               "Cell Line: cellline_A\nb: 0.001\na: 0.016\nreference GR75")
  expect_equal(buildLabel(dt8, "combo1-lines"), 
               "Cell Line: cellline_A\nb: 0.00063\na: 0.01\nIsobol: GR75")
  expect_equal(buildLabel(dt8, "combo1-points"), 
               "Cell Line: cellline_A\nb: 0.00063\na: 0.01\nMX full: 1.2")
  expect_equal(buildLabel(dt8, "combo-ratios"),
               "Cell Line Name: cellline_A\nlog10_ratio_conc: -1.5\nlog2_CI: 0\niso: 0.75")
  expect_equal(buildLabel(dt8, "combo3"), "level: 0.25\nlog2_CI: 0")
  
  expect_error(
    buildLabel(dt1, "dist"),
    paste0("Assertion on 'view' failed: Must be element of set {'clustering','distribution','ranking',",
           "'contrast','curve','grid','combo1-heatmap','combo1-points','combo1-lines_ref',",
           "'combo1-lines','combo-ratios','combo3'}, but is 'dist'."),
    fixed = TRUE
  )
  expect_error(
    buildLabel(dt1, 2),
    "Assertion on 'view' failed: Must be of type 'string', not 'double'.",
    fixed = TRUE
  )
  
})