context("Test metricContrast plot")

m <- 9
n <- 5

synthetic_data <- gDRutils::gen_synthetic_data(m, n)
drug_names <- synthetic_data$drug_names
cell_names <- synthetic_data$cell_names
dt <- synthetic_data$dt

choices <- list(
  primary = "Drug Name",
  secondary = "Cell Line Name",
  drug = drug_names[(n - 3):(n - 2)],
  cell_line = cell_names
)
names(choices) <- c("primary", "secondary", "Drug Name", "Cell Line Name")

variable <- "GR_AOC"

dt_prep <- data.table::data.table(
  "Cell Line Name" = cell_names,
  d1 = dt$GR_AOC[dt$`Drug Name` == drug_names[n - 3]],
  d2 = dt$GR_AOC[dt$`Drug Name` == drug_names[n - 2]],
  check.names = FALSE
)
data.table::setnames(dt_prep, c("d1", "d2"), c(drug_names[n - 3], drug_names[n - 2]))

var_txt <- names(dt_prep)[1]
var_x <- drug_names[2]
var_y <- drug_names[3]
var_col <- "none"
metric <- "GR AOC within range"

dtMissingDrug <- data.table::copy(dt)
dtMissingDrug$`Drug Name` <- NULL
dtMissingDrugMOA <- dt
dtMissingDrugMOA$`Drug MOA` <- NULL
dtMissingCellLine <- dt
dtMissingCellLine$`Cell Line Name` <- NULL
dtMissingPrimaryTissue <- data.table::copy(dt)
dtMissingPrimaryTissue$`Tissue` <- NULL
dtMissingGR_AOC <- dt
dtMissingGR_AOC$GR_AOC <- NULL

choicesMissingPrimary <- choices
choicesMissingPrimary$primary <- NULL
choicesMissingSecondary <- choices
choicesMissingSecondary$secondary <- NULL
choicesMissingDrug <- choices
choicesMissingDrug$`Drug Name` <- NULL
choicesMissingCellLine <- choices
choicesMissingCellLine$`Cell Line Name` <- NULL

# prepare_data_metric_contrast tests
test_that("testing output", {
  prepared_dt <- prepare_data_metric_contrast(dt, choices, variable, var_col)
  data.table::setkey(prepared_dt, NULL)
  expect_equal(prepared_dt, dt_prep)
})

test_that("prepare_data_metric_contrast returns error with wrong argument", {
  expect_error(
    prepare_data_metric_contrast(n, choices, variable, var_col),
    "Assertion on 'data' failed: Must be a data.table, not double."
  )
  expect_error(
    prepare_data_metric_contrast(dt, n, variable, var_col),
    "Assertion on 'choices' failed: Must be of type 'list', not 'double'."
  )
  expect_error(
    prepare_data_metric_contrast(dt, choices, n, var_col),
    "Assertion on 'variable' failed: Must be of type 'string', not 'double'."
  )
  expect_error(
    prepare_data_metric_contrast(dt, choices, variable, n),
    paste0("Assertion on 'var_col' failed: Must be element of set {'Drug Name',",
           "'Drug MOA','Cell Line Name','Tissue','GR_AOC','GR Inf','GR 0','GEC50',",
           "'h GR','E Inf','E0','EC50','h RV','GR50','IC50','GR Max','E Max','GR value',",
           "'Concentration','none'}, but types do not match (numeric != character)."),
    fixed = TRUE
  )
})

test_that("test missing data", {
  expect_error(
    prepare_data_metric_contrast(dtMissingDrug, choices, variable, var_col),
    paste0("Assertion on 'choices[[\"primary\"]]' failed: Must be element of ",
           "set {'Drug MOA','Cell Line Name','Tissue','GR_AOC','GR Inf','GR 0',",
           "'GEC50','h GR','E Inf','E0','EC50','h RV','GR50','IC50','GR Max','E Max',",
           "'GR value','Concentration'}, but is 'Drug Name'."),
    fixed = TRUE
  )
  prepared_dt_1 <- prepare_data_metric_contrast(dtMissingDrugMOA, choices, variable, var_col)
  data.table::setkey(prepared_dt_1, NULL)
  expect_identical(prepared_dt_1, dt_prep)
  expect_error(
    prepare_data_metric_contrast(dtMissingCellLine, choices, variable, var_col),
    paste0("Assertion on 'choices[[\"secondary\"]]' failed: Must be element of ",
           "set {'Drug Name','Drug MOA','Tissue','GR_AOC','GR Inf','GR 0','GEC50',",
           "'h GR','E Inf','E0','EC50','h RV','GR50','IC50','GR Max','E Max','GR value',",
           "'Concentration'}, but is 'Cell Line Name'."),
    fixed = TRUE
  )
  prepared_dt_2 <- prepare_data_metric_contrast(dtMissingPrimaryTissue, choices, variable, var_col)
  data.table::setkey(prepared_dt_2, NULL)
  expect_identical(prepared_dt_2, dt_prep)
  expect_error(
    prepare_data_metric_contrast(dtMissingGR_AOC, choices, variable, var_col),
    paste0("Assertion on 'variable' failed: Must be element of set {'Drug Name',",
           "'Drug MOA','Cell Line Name','Tissue','GR Inf','GR 0','GEC50','h GR',",
           "'E Inf','E0','EC50','h RV','GR50','IC50','GR Max','E Max','GR value',",
           "'Concentration'}, but is 'GR_AOC'."),
    fixed = TRUE
  )
  expect_error(
    prepare_data_metric_contrast(dt, choicesMissingPrimary, variable, var_col),
    "Assertion on 'length(choices) == 4' failed: Must be TRUE.",
    fixed = TRUE
  )
  expect_error(
    prepare_data_metric_contrast(dt, choicesMissingSecondary, variable, var_col),
    "Assertion on 'length(choices) == 4' failed: Must be TRUE.",
    fixed = TRUE
  )
  expect_error(
    prepare_data_metric_contrast(dt, choicesMissingDrug, variable, var_col),
    "Assertion on 'length(choices) == 4' failed: Must be TRUE.",
    fixed = TRUE
  )
  expect_error(
    prepare_data_metric_contrast(dt, choicesMissingCellLine, variable, var_col),
    "Assertion on 'length(choices) == 4' failed: Must be TRUE.",
    fixed = TRUE
  )
})

# plotly_metric_contrast tests
test_that("check output type and data",  {
  plt_0 <- plotly_metric_contrast(dt_prep, var_x, var_y, var_col, var_txt, metric)
  checkmate::expect_class(plt_0, "plotly")
  expect_identical(plt_0$x$attrs[[1]]$type, "scatter")
  expect_identical(plt_0$x$attrs[[1]]$x, dt_prep[[var_x]])
  expect_identical(plt_0$x$attrs[[1]]$y, dt_prep[[var_y]])
  expect_identical(plt_0$x$layoutAttrs[[1]]$xaxis$type, "linear")
  expect_identical(plt_0$x$layoutAttrs[[1]]$yaxis$type, "linear")
  
  plt_1 <- plotly_metric_contrast(dt_prep, var_x, var_y, var_col, var_txt, "IC50")
  checkmate::expect_class(plt_1, "plotly")
  expect_identical(plt_1$x$attrs[[1]]$type, "scatter")
  expect_identical(plt_1$x$attrs[[1]]$x, dt_prep[[var_x]])
  expect_identical(plt_1$x$attrs[[1]]$y, dt_prep[[var_y]])
  expect_identical(plt_1$x$layoutAttrs[[1]]$xaxis$type, "log")
  expect_identical(plt_1$x$layoutAttrs[[1]]$yaxis$type, "log")
  expect_true(grepl("IC", plt_1$x$layoutAttrs[[1]]$title$text))
  
  plt_2 <- plotly_metric_contrast(dt_prep, var_x, var_y, var_col, var_txt, metric, identity = TRUE)
  checkmate::expect_class(plt_2, "plotly")
  expect_identical(plt_2$x$attrs[[1]]$mode, "markers")
  checkmate::expect_class(plt_2$x$attrs[[2]], "plotly_segment") # identity line
  expect_identical(plt_2$x$attrs[[2]]$mode, "lines")
  
  plt_3 <- plotly_metric_contrast(dt_prep, var_x, var_y, var_col, var_txt, metric, correlation = TRUE)
  checkmate::expect_class(plt_3, "plotly")
  expect_identical(plt_3$x$attrs[[1]]$mode, "markers")
  checkmate::expect_class(plt_3$x$attrs[[2]], "plotly_line") # correlation line
  expect_identical(plt_3$x$attrs[[2]]$mode, "lines")
  expect_true(grepl("correlation", plt_3$x$attrs[[3]]$text)) # correlation text
  
  dt_tissue <- data.table::copy(dt_prep)[, Tissue := c("tissue_x", "tissue_x", "tissue_y", "tissue_y", "tissue_z")]
  plt_4 <- plotly_metric_contrast(dt_tissue, var_x, var_y, var_col = "Tissue", var_txt,
                    metric, identity = TRUE, correlation = TRUE)
  checkmate::expect_class(plt_4, "plotly")
  expect_identical(plt_4$x$attrs[[1]]$mode, "markers")
  expect_identical(plt_4$x$attrs[[1]]$color, dt_tissue$Tissue)
  expect_true(plt_4$x$attrs[[1]]$showlegend)
  expect_true(all(grepl("Tissue", plt_4$x$attrs[[1]]$text))) # color info in hover
  checkmate::expect_class(plt_4$x$attrs[[2]], "plotly_line") # correlation line
  expect_identical(plt_4$x$attrs[[2]]$mode, "lines")
  expect_true(grepl("correlation", plt_4$x$attrs[[3]]$text)) # correlation text
  checkmate::expect_class(plt_4$x$attrs[[4]], "plotly_segment") # identity line
  
  dt_conc_2 <- data.table::CJ(
    `Drug Name 2` = "drug_AB",
    `Cell Line Name` = c("cellline_XY", "cellline_XZ", "cellline_YX"),
    `Concentration 2` = c(0.0022, 0.0458, 0.0007, 0.6173, 1.6122, 4.1234))
  dt_conc_2[[var_x]] <- dt[["GR_AOC"]][seq_len(NROW(dt_conc_2))] * 0.50
  dt_conc_2[[var_y]] <- dt[["GR_AOC"]][seq_len(NROW(dt_conc_2))] * 1.25
  plt_5 <- plotly_metric_contrast(dt_conc_2, var_x, var_y, var_col, var_txt,
                    metric, identity = TRUE, correlation = TRUE)
  checkmate::expect_class(plt_5, "plotly")
  expect_identical(plt_5$x$attrs[[1]]$mode, "markers")
  expect_true(plt_5$x$attrs[[1]]$showlegend)
  expect_true(all(grepl("drug_AB at ", plt_5$x$attrs[[1]]$text))) # add concentration info
  expect_true(all(unique(plt_5$x$attrs[[1]]$symbol) %in% unique(dt_conc_2$`Concentration 2`)))
  checkmate::expect_class(plt_5$x$attrs[[2]], "plotly_line") # correlation line
  expect_identical(plt_5$x$attrs[[2]]$mode, "lines")
  expect_true(grepl("correlation", plt_5$x$attrs[[3]]$text)) # correlation text
  checkmate::expect_class(plt_5$x$attrs[[4]], "plotly_segment") # identity line
  
  dt_conc_2_tiss <- data.table::CJ(
    `Drug Name 2` = "drug_001",
    `Cell Line Name` = c("cellline_XY", "cellline_XZ", "cellline_YX", "cellline_YZ"),
    `Tissue` = c("tissue_x", "tissue_x", "tissue_y", "tissue_z"),
    `Concentration 2` = c(0.0007, 0.1851, 1.6666))
  dt_conc_2_tiss$drug_002 <- dt[["GR_AOC"]][seq_len(NROW(dt_conc_2_tiss))] * 0.50
  dt_conc_2_tiss$drug_003 <- dt[["GR_AOC"]][seq_len(NROW(dt_conc_2_tiss))] * 2.25
  dt_conc_2_tiss <- na.omit(dt_conc_2_tiss)
  plt_6 <- plotly_metric_contrast(dt_conc_2_tiss,
                    var_x, var_y, var_col = "Tissue", var_txt,
                    metric, identity = TRUE, correlation = TRUE)
  checkmate::expect_class(plt_6, "plotly")
  expect_identical(plt_6$x$attrs[[2]]$mode, "markers")
  expect_true(all(unique(plt_6$x$attrs[[2]]$symbol) %in% unique(dt_conc_2_tiss$`Concentration 2`)))
  expect_true(plt_6$x$attrs[[2]]$showlegend)
  expect_identical(plt_6$x$attrs[[3]]$color, dt_conc_2_tiss$Tissue) # legend
  expect_true(plt_6$x$attrs[[3]]$showlegend)
  expect_identical(plt_6$x$attrs[[4]]$color, dt_conc_2_tiss$Tissue) # plot
  expect_true(all(unique(plt_6$x$attrs[[4]]$symbol) %in% unique(dt_conc_2_tiss$`Concentration 2`)))
  expect_false(plt_6$x$attrs[[4]]$showlegend)
  expect_true(all(grepl("Tissue", plt_6$x$attrs[[4]]$text))) # color info in hover
  checkmate::expect_class(plt_6$x$attrs[[5]], "plotly_line") # correlation line
  expect_identical(plt_6$x$attrs[[5]]$mode, "lines")
  expect_true(grepl("correlation", plt_6$x$attrs[[6]]$text)) # correlation text
  checkmate::expect_class(plt_6$x$attrs[[7]], "plotly_segment") # identity line
})

test_that("check input arguments ", {
  expect_error(
    plotly_metric_contrast(n, var_x, var_y, var_col, var_txt),
    "Assertion on 'data' failed: Must be a data.table, not double."
  )
  expect_error(
    plotly_metric_contrast(dt_prep, n, var_y, var_col, var_txt),
    "Assertion on 'var_x' failed: Must be of type 'string', not 'double'."
  )
  expect_error(
    plotly_metric_contrast(dt_prep, var_x, n, var_col, var_txt),
    "Assertion on 'var_y' failed: Must be of type 'string', not 'double'."
  )
  expect_error(
    plotly_metric_contrast(dt_prep, var_x, var_y, n, var_txt),
    paste0("Assertion on 'var_col' failed: Must be element of set ",
           "{'Cell Line Name','drug_002','drug_003','none'}, but types do not match ",
           "(numeric != character)."),
    fixed = TRUE
  )
  expect_error(
    plotly_metric_contrast(dt_prep, var_x, var_y, var_col, n),
    "Assertion on 'var_txt' failed: Must be of type 'string', not 'double'."
  )
})
