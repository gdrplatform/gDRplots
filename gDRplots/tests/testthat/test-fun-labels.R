context("Test fun labels")

test_that("adjustLabel works as expected", {
  ls_lbl <- c(
    "G03405395 (G03642866 at 0.0007621 &mu;M)", 
    "G03580756 (G03642866 at 0.0007621 &mu;M)",
    "G03405395 (G03642866 at 0.002286 &mu;M)")
  ls_lbl_res <- c(
    "G03405395 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(G03642866 at 0.0007621 &mu;M)",
    "G03580756 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(G03642866 at 0.0007621 &mu;M)",
    "G03405395 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(G03642866 at 0.002286 &mu;M)"
  )
  expect_equal(unname(adjustLabel(ls_lbl)), ls_lbl_res) 
  
  ls_lbl_2 <- c(
    "G03405395\n(G03642866 at 0.0007621 &mu;M)", 
    "G03580756\n(G03642866 at 0.0007621 &mu;M)",
    "G03405395\n(G03642866 at 0.002286 &mu;M)")
  expect_equal(adjustLabel(ls_lbl_2), ls_lbl_2) 
  
  expect_error(adjustLabel(1:5),
               "Assertion on 'x' failed: Must be of type 'character'")
  expect_error(adjustLabel(ls_lbl, pattern = "[0-9]{2}"),
               "Assertion on 'pattern' failed: Must comply to pattern")
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
