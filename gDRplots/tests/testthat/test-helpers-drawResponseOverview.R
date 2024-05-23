context("Test drawResponseOverview helpers")

# single-agent test data ----
synthetic_data <- gDRutils::get_testdata()
drug_names <- synthetic_data$drug_names
cell_names <- synthetic_data$cell_line_names
dt <- synthetic_data$dt

prepared_curves <- prepareCurves(dt)
var_y <- "GR value"
key_cells_drugs <- data.table::data.table(
  "Cell Line Name" = cell_names[1],
  "Drug Name" = drug_names[1]
)

prepared_extras <- prepareExtras(dt)
subset_data <- dt[key_cells_drugs, on = intersect(names(dt), names(key_cells_drugs))]
subset_curves <- prepared_curves[key_cells_drugs, on = intersect(names(prepared_curves),
                                                                 names(key_cells_drugs))]
# prepare extras for selected
subset_extras_points <- prepared_extras$points[key_cells_drugs, 
                                               on = intersect(names(prepared_curves),
                                                              names(key_cells_drugs))]
subset_extras_lines <- prepared_extras$lines[grepl(paste0(key_cells_drugs, collapse = "."), 
                                                   names(prepared_extras$lines))]
subset_extras <- list(points = subset_extras_points,
                      lines = subset_extras_lines)

# combo test data ----
combo_synthetic_data <- gDRutils::get_testdata_combo()
combo_drug_names <- combo_synthetic_data$drug_names
combo_cell_names <- combo_synthetic_data$cell_line_names
combo_dt <- combo_synthetic_data$dt

combo_prepared_curves <- prepareCurves(combo_dt)
combo_key_cells_drugs <- data.table::data.table(
  "Cell Line Name" = combo_cell_names[2:4],
  "Drug Name" = combo_drug_names[2:4]
)

combo_prepared_extras <- prepareExtras(combo_dt)
combo_subset_data <- 
  combo_dt[combo_key_cells_drugs, on = intersect(names(combo_dt), names(combo_key_cells_drugs))]
combo_subset_curves <- 
  combo_prepared_curves[combo_key_cells_drugs, 
                        on = intersect(names(combo_prepared_curves), names(combo_key_cells_drugs))]
# prepare extras for selected
combo_subset_extras_points <- 
  combo_prepared_extras$points[
    combo_key_cells_drugs, 
    on = intersect(names(combo_prepared_curves), names(combo_key_cells_drugs))]
combo_subset_extras_lines <- 
  combo_prepared_extras$lines[grepl(paste0(combo_key_cells_drugs, collapse = "."), 
                                    names(combo_prepared_extras$lines))]
combo_subset_extras <- list(points = combo_subset_extras_points,
                            lines = combo_subset_extras_lines)

# prepareCurves tests ----
test_that("prepareCurves works as expected", {
  prepared_curves <- prepareCurves(dt)
  expect_is(prepared_curves, "data.table")
  expect_equal(NROW(prepared_curves), 100 * NROW(dt))
  expect_equal(min(unique(prepared_curves$Concentration)), 1e-3)
  expect_equal(max(unique(prepared_curves$Concentration)), 50e+0)
  
  n_density <- 200
  expect_equal(NROW(prepareCurves(dt, density = n_density)), n_density * NROW(dt))
  
  con_range <- c(1e-2, 5)
  expect_equal(range(prepareCurves(dt, range_x = con_range)$Concentration), con_range)
  expect_equal(unique(prepareCurves(dt, range_x = c(2, 2), density = 1)$Concentration), 2)
  expect_equal(NROW(prepareCurves(dt, range_x = c(2, 2), density = 1)), NROW(dt))
  
  dt2 <- dt
  dt2[2, ]$`Drug Name` <- dt2[1, `Drug Name`] # duplicated Drug Name
  expect_error(prepareCurves(dt2), "sth wrong with the data model")
  expect_error(prepareCurves(as.list(dt)))
  expect_error(prepareCurves(dt, range_x = 1))
  expect_error(prepareCurves(dt, density = "str"))
})

# prepareExtras tests ----
test_that("check output type, data, values for prepareExtras", {
  ext_1 <- prepareExtras(dt)
  
  expect_type(ext_1, "list")
  expect_identical(names(ext_1), c("points", "lines"))
  expect_identical(ext_1$points$`Cell Line`, rep(dt$`Cell Line Name`, 2))
  expect_identical(ext_1$points$Drug, rep(dt$`Drug Name`, 2))
  expect_identical(ext_1$points$Concentration, c(dt$EC50, dt$GR50))
  expect_identical(ext_1$lines[[1]]$type, "line")
  
  r_x <- c(3.8e-2, 4e+0)
  ext_2 <- prepareExtras(dt, range_x = r_x)
  
  expect_type(ext_2, "list")
  expect_identical(names(ext_2), c("points", "lines"))
  expect_identical(ext_2$points, ext_1$points)
  expect_identical(ext_2$lines[[1]], ext_1$lines[[1]])
  horizontal_1st <- NROW(ext_2$lines) / 2 + 1
  expect_identical(ext_2$lines[[horizontal_1st]]$x0, r_x[1])
  expect_identical(ext_2$lines[[horizontal_1st]]$x1, r_x[2])
})

test_that("returns error with wrong argument for prepareExtras", {
  expect_error(
    prepareExtras(),
    "argument \"metrics\" is missing, with no default"
  )
  expect_error(
    prepareExtras(5),
    "Assertion on 'metrics' failed: Must be a data.table, not double."
  )
  expect_error(
    prepareExtras(dt, 1),
    "Assertion on 'range_x' failed: Must have length 2, but has length 1."
  )
  expect_error(
    prepareExtras(dt, "str"),
    "Assertion on 'range_x' failed: Must be of type 'numeric', not 'character'."
  )
})

# plotlyRCAll tests ----
test_that("check output type for plotlyRCAll", {
  plt_all <- plotlyRCAll(prepared_curves, var_y) # default
  
  expect_type(plt_all, "list")
  expect_is(plt_all, "plotly")
  expect_equal(plt_all$x$source, "curvePlot")
  expect_equal(plt_all$width, 400)
  expect_equal(plt_all$height, 300)
  expect_equal(plt_all$x$layoutAttrs[[1]]$xaxis$range, log10(c(1e-3, 50e+0)))
  expect_equal(plt_all$x$layoutAttrs[[1]]$title$text,
               paste0("Dose response curves for Drug Name: ", 
                      paste(drug_names, collapse = ", ")))
  
  plt_w <- 200
  plt_h <- 600
  plt <- plotlyRCAll(prepared_curves, var_y, plot_width = plt_w, plot_height = plt_h)
  expect_type(plt, "list")
  expect_is(plt, "plotly")
  expect_equal(plt$width, plt_w)
  expect_equal(plt$height, plt_h)
  
  r_x <- c(3.8e-2, 4e+0)
  plt <- plotlyRCAll(prepared_curves, var_y, range_x = r_x)
  expect_type(plt, "list")
  expect_is(plt, "plotly")
  expect_equal(plt$x$layoutAttrs[[1]]$xaxis$range, log10(r_x))
  
  selected_cl <- cell_names[1:5]
  dt2 <- dt[`Cell Line Name` %in% selected_cl]
  prepared_curves_2 <- prepareCurves(dt2)
  
  plt_all_2 <- plotlyRCAll(prepared_curves_2, var_y)
  expect_is(plt_all_2, "plotly")
  expect_equal(plt_all_2$x$layoutAttrs[[1]]$title$text,
               paste0("Dose response curves for Cell Line Name: ", 
                      paste(selected_cl, collapse = ", ")))
})

test_that("check output for wrong data for plotlyRCAll", {
  prepared_curves_out <- prepared_curves
  prepared_curves_out[[var_y]] <- prepared_curves_out[[var_y]] + 100
  
  plt <- plotlyRCAll(prepared_curves_out, var_y = var_y)
  plt_msg <- plt$x$layoutAttrs[[1]]$annotations$text
  
  pidfs <- gDRutils::get_prettified_identifiers(simplify = TRUE)
  drug_name <- pidfs[["drug_name"]]
  cell_name <- pidfs[["cellline_name"]]
  comb_name <- paste0(prepared_curves_out[[cell_name]][1], " x ", prepared_curves_out[[drug_name]][1])
  
  expect_equal(plt$x$attrs[[1]]$mode, "text")
  expect_true(grepl("Invalid averaged data", plt_msg))
  expect_true(grepl(comb_name, plt_msg))
  expect_true(grepl("contact gdrplatform team", plt_msg))
})

test_that("returns error with missing and wrong argument for plotlyRCAll", {
  expect_error(
    plotlyRCAll(dt, var_y),
    paste0("Assertion on 'names\\(curve_data\\)' failed:")
  )
  expect_error(
    plotlyRCAll(5, var_y),
    paste0("Assertion on 'curve_data' failed: Must be a data.table, not double.")
  )
  expect_error(
    plotlyRCAll(prepared_curves, 5),
    "Assertion on 'var_y' failed: Must be of type 'string', not 'double'."
  )
  expect_error(
    plotlyRCAll(prepared_curves, var_y, range_x = 1),
    "Assertion on 'range_x' failed: Must have length 2, but has length 1."
  )
  expect_error(
    plotlyRCAll(prepared_curves, var_y, range_x = "str"),
    "Assertion on 'range_x' failed: Must be of type 'numeric', not 'character'."
  )
  expect_error(
    plotlyRCAll(prepared_curves, var_y, plot_width = "str"),
    "Assertion on 'plot_width' failed: Must be of type 'numeric', not 'character'."
  )
  expect_error(
    plotlyRCAll(prepared_curves, var_y, plot_height = "str"),
    "Assertion on 'plot_height' failed: Must be of type 'numeric', not 'character'."
  )
  
  prepared_curves_2 <- data.table::copy(prepared_curves)
  data.table::setnames(prepared_curves_2, "Drug Name", "Bad Name")
  expect_error(
    plotlyRCAll(prepared_curves_2, var_y),
    "failed: Names must include the elements"
  )
})

# plotlyRCSelected tests ----
test_that("check output type for plotlyRCSelected", {
  plt_sel_1 <- plotlyRCSelected(data = subset_data, 
                                var_y = var_y, 
                                layers = c("curve", "average", "error"),
                                curves = subset_curves, 
                                extras = subset_extras)
  
  expect_is(plt_sel_1, "plotly")
  expect_equal(plt_sel_1$x$attrs[[2]]$type, "scatter")
  expect_equal(plt_sel_1$x$attrs[[2]]$mode, "lines")
  expect_equal(plt_sel_1$width, 400)
  expect_equal(plt_sel_1$height, 300)
  expect_equal(plt_sel_1$x$layoutAttrs[[1]]$xaxis$range, log10(c(1e-3, 50e+0)))
  expect_equal(plt_sel_1$x$attrs[[2]]$color, rep(key_cells_drugs$`Cell Line Name`, 100))
  
  plt_w <- 200
  plt_h <- 600
  plt_sel_2 <- plotlyRCSelected(data = subset_data, 
                                var_y = var_y, 
                                layers = c("curve", "extras"),
                                curves = subset_curves, 
                                extras = subset_extras,
                                plot_width = plt_w,
                                plot_height = plt_h)
  
  plt_id <- names(plt_sel_2$x$visdat)
  expect_is(plt_sel_2, "plotly")
  expect_equal(plt_sel_2$width, plt_w)
  expect_equal(plt_sel_2$height, plt_h)
  
  # curve data
  expect_equal(plt_sel_2$x$attrs[[2]]$x, subset_curves[["Concentration"]])
  expect_equal(plt_sel_2$x$attrs[[2]]$y, subset_curves[[var_y]])
  # points for half response
  expect_equal(plt_sel_2$x$attrs[[3]]$x, subset_extras$points[metric == var_y]$Concentration)
  expect_equal(plt_sel_2$x$attrs[[3]]$y, subset_extras$points[metric == var_y]$Response)
  # points ticks for half response
  expect_equal(plt_sel_2$x$attrs[[4]]$y, subset_extras$points[metric == var_y]$MaxEffectiveness)
  # points ticks for max effectiveness
  expect_equal(plt_sel_2$x$attrs[[5]]$x, subset_extras$points[metric == var_y]$Concentration)
  # horizontal and vertical line for Extras
  sub_ex <- subset_extras$lines[grepl(var_y, names(subset_extras$lines))]
  names(sub_ex) <- NULL
  expect_equal(plt_sel_2$x$layoutAttrs[[plt_id]]$shapes[[5]], sub_ex[[1]])
  expect_equal(plt_sel_2$x$layoutAttrs[[plt_id]]$shapes[[6]], sub_ex[[2]])
  expect_equal(plt_sel_2$x$config$edits, get_plotly_edits())
  
  plt_sel_3 <- plotlyRCSelected(data = subset_data, 
                                var_y = var_y, 
                                layers = c("observations", "error"),
                                curves = subset_curves, 
                                extras = subset_extras)
  
  expect_is(plt_sel_3, "plotly")
  expect_equal(plt_sel_3$width, 400)
  expect_equal(plt_sel_3$height, 300)
  # observations data
  expect_equal(plt_sel_3$x$attrs[[2]]$type, "scatter")
  expect_equal(plt_sel_3$x$attrs[[2]]$mode, "markers")
  expect_equal(plt_sel_3$x$attrs[[2]]$x, subset_data[["Concentration"]])
  expect_equal(plt_sel_3$x$attrs[[2]]$y, subset_data[[var_y]])
  
  r_x <- c(3.8e-2, 4e+0)
  plt_sel_4 <- plotlyRCSelected(data = subset_data, 
                                var_y = var_y, 
                                layers = c("observations", "error"),
                                curves = subset_curves, 
                                range_x = r_x,
                                extras = subset_extras)
  
  expect_is(plt_sel_4, "plotly")
  expect_equal(plt_sel_4$x$layoutAttrs[[1]]$xaxis$range, log10(r_x))
  
  selected_cl <- cell_names[1:3]
  plt_sel_5 <- plotlyRCSelected(data = dt[`Cell Line Name` %in% selected_cl, ], 
                                var_y = var_y, 
                                layers = c("observations", "error"),
                                curves = prepared_curves[`Cell Line Name` %in% selected_cl, ],
                                extras = subset_extras)
  
  expect_is(plt_sel_5, "plotly")
  expect_equal(plt_sel_5$x$layoutAttrs[[1]]$title$text,
               paste0("Drug dose response for Cell Line Name: ", 
                      paste(selected_cl, collapse = ", ")))
  
  # fast end - plotly_empty
  plt_sel_4 <- plotlyRCSelected(data = subset_data, 
                                var_y = var_y, 
                                layers = NULL,
                                curves = subset_curves, 
                                extras = subset_extras)
  expect_is(plt_sel_4, "plotly")
  expect_equal(plt_sel_4$x$attrs[[1]]$mode, "markers")
  expect_equal(plt_sel_4$x$attrs[[1]]$x, NULL)
  
  subset_data_na <- data.table::copy(subset_data)
  subset_data_na[[var_y]] <- NA
  plt_sel_5 <- plotlyRCSelected(data = subset_data_na, 
                                var_y = var_y, 
                                layers = c("observations", "error"),
                                curves = subset_curves, 
                                extras = subset_extras)
  expect_is(plt_sel_5, "plotly")
  expect_equal(plt_sel_5$x$attrs[[1]]$mode, "markers")
  expect_equal(plt_sel_5$x$attrs[[1]]$x, NULL)
})

test_that("check output type for plotlyRCSelected for combo", {
  plt_sel_1 <- plotlyRCSelected(data = combo_subset_data, 
                                var_y = var_y, 
                                layers = c("curve", "average", "error"),
                                curves = combo_subset_curves, 
                                extras = combo_subset_extras)
  
  expect_is(plt_sel_1, "plotly")
  expect_equal(plt_sel_1$x$attrs[[2]]$type, "scatter")
  expect_equal(plt_sel_1$x$attrs[[2]]$mode, "lines")
  expect_equal(plt_sel_1$width, 400)
  expect_equal(plt_sel_1$height, 300)
  expect_equal(plt_sel_1$x$layoutAttrs[[1]]$xaxis$range, log10(c(1e-3, 50e+0)))
  expect_equal(sort(plt_sel_1$x$attrs[[2]]$color), 
               sort(rep(combo_key_cells_drugs$`Cell Line Name`, 100)))
  
})

test_that("check output for wrong data for plotlyRCSelected", {
  subset_data_out <- subset_data
  subset_data_out[[var_y]] <- subset_data_out[[var_y]] + 100
  
  
  plt <- plotlyRCSelected(data = subset_data_out, 
                          var_y = "GR value", 
                          layers = c("curve", "average", "error"),
                          curves = subset_curves, 
                          extras = subset_extras)
  plt_msg <- plt$x$layoutAttrs[[1]]$annotations$text
  
  pidfs <- gDRutils::get_prettified_identifiers(simplify = TRUE)
  drug_name <- pidfs[["drug_name"]]
  cell_name <- pidfs[["cellline_name"]]
  comb_name <- paste0(subset_data_out[[cell_name]][1], " x ", subset_data_out[[drug_name]][1])

  expect_equal(plt$x$attrs[[1]]$mode, "text")
  expect_true(grepl("Invalid averaged data", plt_msg))
  expect_true(grepl(comb_name, plt_msg))
  expect_true(grepl("contact gdrplatform team", plt_msg))
})

test_that("expect error with wrong input", {
  expect_error(
    plotlyRCSelected(5, var_y, layers = c("curve", "average", "error"),
                     subset_curves, extras = subset_extras),
    "Assertion on 'data' failed: Must be a data.table, not double."
  )
  expect_error(
    plotlyRCSelected(subset_data, 5, layers = c("curve", "average", "error"),
                     subset_curves, extras = subset_extras),
    "Assertion on 'var_y' failed: Must be of type 'string', not 'double'."
  )
  expect_error(
    plotlyRCSelected(subset_data, var_y, layers = 5,
                     subset_curves, extras = subset_extras),
    "Assertion on 'layers' failed: Must be of type 'character' (or 'NULL'), not 'double'.",
    fixed = TRUE
  )
  expect_error(
    plotlyRCSelected(subset_data, var_y, layers = "curves",
                     subset_curves, extras = subset_extras),
    paste0("Assertion on 'X[[i]]' failed: Must be element of set {'curve',",
           "'average','error','observations','extras'}, but is 'curves'."),
    fixed = TRUE
  )
  expect_error(
    plotlyRCSelected(subset_data, var_y, layers = c("curve", "average", "error"),
                     5, extras = subset_extras),
    "Assertion on 'curves' failed: Must be a data.table, not double."
  )
  expect_error(
    plotlyRCSelected(subset_data, var_y, layers = c("curve", "average", "error"),
                     subset_curves, extras = 5),
    "Assertion on 'extras' failed: Must be of type 'list', not 'double'."
  )
  expect_error(
    plotlyRCSelected(subset_data, var_y, layers = c("curve", "average", "error"),
                     subset_curves, extras = subset_extras, range = 1),
    "Assertion on 'range_x' failed: Must have length 2, but has length 1."
  )
  expect_error(
    plotlyRCSelected(subset_data, var_y, layers = c("curve", "average", "error"),
                     subset_curves, extras = subset_extras, range = "str"),
    "Assertion on 'range_x' failed: Must be of type 'numeric', not 'character'."
  )
  expect_error(
    plotlyRCSelected(subset_data, var_y, layers = c("curve", "average", "error"),
                     subset_curves, extras = subset_extras, plot_width = "str"),
    "Assertion on 'plot_width' failed: Must be of type 'numeric', not 'character'."
  )
  expect_error(
    plotlyRCSelected(subset_data, var_y, layers = c("curve", "average", "error"),
                     subset_curves, extras = subset_extras, plot_height = "str"),
    "Assertion on 'plot_height' failed: Must be of type 'numeric', not 'character'."
  )
})


test_that("logSeq works as expected", {
  sequence <- c(1, 1.18920711500272, 1.4142135623731, 1.68179283050743, 2)
  
  # sequence is generated
  expect_equal(logSeq(1, 2, 5), sequence)
  # sequence has proper length
  expect_equal(length(logSeq(1, 2, 5)), 5)
  # sequence has proper limits
  expect_equal(logSeq(1, 2, 5)[1], 1)
  expect_equal(logSeq(1, 2, 5)[length(logSeq(1, 2, 5))], 2)
  # sequence grows linearly in log domain
  s <- logSeq(1, 2, 5)
  sLog <- log10(s)
  differences <- diff(sLog)
  differences <- signif(differences, 10)
  expect_true(length(unique(differences)) == 1L)
  
  expect_error(logSeq(NULL, 2, 5))
  expect_error(logSeq(-2, 2, 5))
  expect_error(logSeq(1, -5, 5))
  expect_error(logSeq(1, "str", 5))
  expect_error(logSeq(1, 2, TRUE))
  expect_error(logSeq(1, 2, 0))
})


test_that("getLongest works as expected", {
  factorList <- list("letters" = letters, "numbers" = seq_len(10))
  expect_identical(getLongest(factorList), "letters")
  expect_identical(getLongest(factorList, "letters"), "letters")
  
  factorList <- c(factorList, list("LETTERS" = LETTERS))
  expect_identical(getLongest(factorList, "letters"), "letters")
  expect_identical(getLongest(factorList, 3), names(factorList)[3])
  expect_identical(getLongest(factorList, "LETTERS"), "LETTERS")
  
  factorDatatable <- data.table::data.table(
    "letters" = letters,
    "numbers" = rep(seq_len(10), length.out = NROW(letters))
  )
  expect_identical(getLongest(factorDatatable), "letters")
  
  factorDatatable <- data.table::data.table(
    "letters" = letters,
    "numbers" = seq_len(NROW(letters))
  )
  expect_identical(getLongest(factorDatatable), "letters")
  expect_identical(getLongest(factorDatatable, "letters"), "letters")
  expect_identical(getLongest(factorDatatable, 2L), names(factorDatatable)[2L])
  
  expect_error(getLongest(letters), "Must inherit from class")
  expect_error(getLongest(unname(factorList)), "Must have Object")
  expect_error(getLongest(list(f = factorList, i = iris)), 
               "Must be of type 'atomic vector'")
  expect_error(getLongest(factorDatatable, TRUE),
               "Must inherit from class 'character'/'integer'/'numeric'")
  expect_error(getLongest(factorDatatable, "str"), "Must be element of set")
  expect_error(getLongest(factorDatatable, 5), "Assertion on 'default' failed")
})

