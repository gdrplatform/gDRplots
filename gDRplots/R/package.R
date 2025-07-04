#' @importFrom data.table :=
#' @importFrom data.table %chin%
NULL


# Make sure data.table knows we know we're using it
.datatable.aware <- TRUE

# Prevent R CMD check from complaining about the standard data.table variables
utils::globalVariables(
  c(".",
    ".N",
    ".SD" ,
    "abs_rho",
    "CCLEName",
    "cId",
    "cotrt_value",
    "cotrt_value_zero",
    "DrugCombination" ,
    "ec50",
    "feat_val",
    "feature",
    "h",
    "iso_level" ,
    "iso_source",
    "label",
    "legend_lbl_iso",
    "log10_ratio_conc",
    "log2_CI",
    "ModelID",
    "N",
    "NES" ,
    "p_value",
    "padj",
    "pos_x",
    "pos_y",
    "pval",
    "q_value",
    "ReadoutValue",
    "rId",
    "stat_sig",
    "value",
    "WellColumn",
    "WellRow",
    "x",
    "x_0" ,
    "x_inf",
    "x_pos",
    "x_std",
    "y_pos"
  ),
  utils::packageName())
