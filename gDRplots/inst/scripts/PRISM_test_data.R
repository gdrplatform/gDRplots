# Model.csv ----
cell_lines <- gDRtestData::create_synthetic_cell_lines()[["CellLineName"]]
no_cell_lines <- NROW(cell_lines)

tab_model <- data.table::data.table(
  ModelID = sprintf("ACH-%06d", seq_along(cell_lines)),
  CCLEName = cell_lines,
  OncotreeLineage = 
    withr::with_seed(42, 
                     sample(c("Soft Tissue", "Skin", "Lung", "Liver", "Breast", "Kidney", "Other"),
                            no_cell_lines, replace = TRUE)),
  Age = withr::with_seed(42, sample(c(18:98, NA), no_cell_lines, replace = TRUE)),
  GrowthPattern = 
    rep(c("Adherent", "Adherent",  NA, "Suspension", "Mixed", "Unknown", "Neurosphere", "Organoid", ""),
        length.out = no_cell_lines, replace = TRUE),
  PatientRace = rep(c("asian", "african", "caucasian", "caucasian", "unknown", "hispanic_or_latino"), 
                    length.out = no_cell_lines, replace = TRUE),
  Sex = as.factor(withr::with_seed(42, sample(c("Female",  "Male", "Female", "Male", "Unknown"), 
                                              no_cell_lines, replace = TRUE))),
  SourceDetail = withr::with_seed(42, sample(c(TRUE, FALSE, NA), no_cell_lines, replace = TRUE)),
  TreatmentStatus = 
    withr::with_seed(42, 
                     sample(c( "", "Unknown", "Post-treatment", "Pre-treatment", "Active treatment", NA), 
                            no_cell_lines, replace = TRUE))
)

data.table::fwrite(tab_model, "./gDRplots/inst/testdata/Model.csv", row.names = FALSE)


# CRISPRGeneEffect.csv ----
tab_cko <- data.table::data.table(
  V1 = tab_model$ModelID,
  "XZ_A1QW (123)" = withr::with_seed(42, rnorm(n = NROW(tab_model), mean = -0.05, sd = 0.11)),
  "XZ_A2GH (456987)" = withr::with_seed(42, rnorm(n = NROW(tab_model), mean = -0.03, sd = 0.13)),
  "XZ_A3OP (Unknown)" = withr::with_seed(42, rnorm(n = NROW(tab_model), mean = 0.045, sd = 0.10)),
  "XZ_A4RT ()" = withr::with_seed(42, rnorm(n = NROW(tab_model), mean = 0.05, sd = 0.10)),
  "XZ_A5BN " = withr::with_seed(42, rnorm(n = NROW(tab_model), mean = 0.11, sd = 0.13))
)
tab_cko[[2]][4:5] <- NA

data.table::fwrite(tab_cko, "./gDRplots/inst/testdata/CRISPRGeneEffect.csv", row.names = FALSE)

# OmicsSomaticMutationsMatrixHotspot.csv ----
tab_hot <- data.table::data.table(
  V1 = tab_model$ModelID,
  "NU_X1QW (456)" = withr::with_seed(42, sample(c(0, 1), size = NROW(tab_model), replace = TRUE)),
  "NU_X2GH (unknown)" = withr::with_seed(314, sample(c(0, 1), size = NROW(tab_model), replace = TRUE)),
  "NU_X3OP (523)" = withr::with_seed(271, sample(c(0, 1), size = NROW(tab_model), replace = TRUE)),
  "NU_X4RT (4789)" = withr::with_seed(981, sample(c(0, 1, 2), size = NROW(tab_model), replace = TRUE)),
  "NU_X5BN" = rep(1, size = NROW(cell_lines))
)
tab_hot[[3]][11:12] <- NA

data.table::fwrite(tab_hot, "./gDRplots/inst/testdata/OmicsSomaticMutationsMatrixHotspot.csv", row.names = FALSE)
