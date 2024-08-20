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
    "iso_level",
    "log10_ratio_conc",
    "log2_CI",
    "pos_x",
    "pos_y",
    "x",
    "x_std",
    "cId",
    "rId",
    "p_value"
  ),
  utils::packageName())

