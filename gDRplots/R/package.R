#' @importFrom data.table :=
NULL


# Make sure data.table knows we know we're using it
.datatable.aware <- TRUE

# Prevent R CMD check from complaining about the standard data.table variables
utils::globalVariables(
  c(
    ".",
    ".N",
    ".SD",
    "E Inf",
    "E Max",
    "E0",
    "EC50",
    "GEC50",
    "GR 0",
    "GR Inf",
    "GR Max",
    "GR50",
    "h GR",
    "h RV",
    "IC50",
    "iso_level",
    "log10_ratio_conc",
    "log2_CI",
    "MaxEffectiveness",
    "metric",
    "pos_x",
    "pos_y",
    "Response",
    "x"
  ),
  utils::packageName())
