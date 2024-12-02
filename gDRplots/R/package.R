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
    "p_value",
    "CCLEName",
    "cotrt_value",
    "cotrt_value_zero",
    "DrugCombination",
    "feature",
    "feat_val",
    "iso_source",
    "label",
    "legend_lbl_iso",
    "ModelID",
    "N",
    "stat_sig",
    "q_value",
    "value",
    "ReadoutValue",
    "WellColumn",
    "WellRow"
  ),
  utils::packageName())

