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
  Age = withr::with_seed(42, sample(c(18, 25, 45, 65:68, 89, 90, NA), no_cell_lines, replace = TRUE)),
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

# OmicsSignaturesProfile ----
chars <- c(letters, LETTERS, 0:9)
no_profils <- 100
tab_profil <- data.table::data.table(
  V1 = sprintf("PR-%s", 
               vapply(seq_len(no_profils), function(i) {
                 withr::with_seed(i, paste(sample(chars, 8), collapse = "")) 
                 }, character(1))),
  MSIScore = round(withr::with_seed(42, rnorm(n = no_profils, mean = 2.56, sd = 1.25)), 2), 
  Ploidy = withr::with_seed(42, rnorm(n = no_profils, mean = 6.5, sd = 2.11)),
  CIN = withr::with_seed(42, rnorm(n = no_profils, mean = 0.5, sd = 0.02))
)

data.table::fwrite(tab_profil, "./gDRplots/inst/testdata/OmicsSignaturesProfile.csv", row.names = FALSE)


# assoc table ----
tab_drug_001_met1 <- data.table::data.table(
  feature = sprintf("NU_%03d_X1%s", 1:25, LETTERS[1:25]),
  response = rep("RV_gDR_x_max", 25),
  rho = withr::with_seed(42, sample(seq(-1, 1, 0.35), 25, replace = TRUE)),
  q_value = withr::with_seed(314, sample(seq(0.001, 0.1, 0.0025), 25, replace = TRUE))
)
tab_drug_001_met1$neglog_q_value <- -log10(tab_drug_001_met1$q_value)

writexl::write_xlsx(tab_drug_001_met1, 
                    "./gDRplots/inst/testdata/tab_assoc_RV__FEAT_NU_drug_001_RV_gDR_x_max.xlsx")

tab_drug_001_met2 <- data.table::data.table(
  feature = sprintf("NU_%03d_X1%s", 1:25, LETTERS[1:25]),
  response = rep("RV_gDR_x_mean", 25),
  rho = withr::with_seed(42, sample(seq(-0.85, 1.55, 0.35), 25, replace = TRUE)),
  q_value = withr::with_seed(42, sample(seq(0.001, 0.55, 0.0015), 25, replace = TRUE))
)
tab_drug_001_met2$neglog_q_value <- -log10(tab_drug_001_met2$q_value)

writexl::write_xlsx(tab_drug_001_met2, 
                    "./gDRplots/inst/testdata/tab_assoc_RV__META_GRP_drug_001_RV_gDR_x_mean.xlsx")


tab_drug_002 <- data.table::data.table(
  feature = sprintf("GRP_%03d_XC", 1:25),
  response = rep("RV_gDR_x_mean", 25),
  rho = withr::with_seed(314, sample(seq(-0.85, 1.55, 0.35), 25, replace = TRUE)),
  q_value = withr::with_seed(42, sample(seq(0.001, 0.55, 0.0025), 25, replace = TRUE))
)
tab_drug_002$neglog_q_value <- -log10(tab_drug_002$q_value)

writexl::write_xlsx(tab_drug_002, 
                    "./gDRplots/inst/testdata/tab_assoc_RV__META_GRP_drug_002_RV_gDR_x_mean.xlsx")

tab_drug_001_drug_021_met3 <- data.table::data.table(
  feature = sprintf("GRP_%03d_XC", 1:25),
  response = rep("GR_gDR_hsa_score", 25),
  rho = withr::with_seed(42, sample(seq(-1, 1, 0.35), 25, replace = TRUE)),
  q_value = withr::with_seed(314, sample(seq(0.001, 0.1, 0.0025), 25, replace = TRUE))
)
tab_drug_001_drug_021_met3$neglog_q_value <- -log10(tab_drug_001_drug_021_met3$q_value)

writexl::write_xlsx(tab_drug_001_drug_021_met3, 
                    "./gDRplots/inst/testdata/tab_assoc_GR__META_GRP_drug_001_x_drug_021_GR_gDR_hsa_score.xlsx")

tab_drug_001_drug_021_met4 <- data.table::data.table(
  feature = sprintf("GRP_%03d_XC", 1:25),
  response = rep("GR_gDR_bliss_score", 25),
  rho = withr::with_seed(314, sample(seq(-1, 1, 0.35), 25, replace = TRUE)),
  q_value = withr::with_seed(314, sample(seq(0.001, 0.55, 0.0015), 25, replace = TRUE))
)
tab_drug_001_drug_021_met4$neglog_q_value <- -log10(tab_drug_001_drug_021_met4$q_value)

writexl::write_xlsx(tab_drug_001_drug_021_met4, 
                    "./gDRplots/inst/testdata/tab_assoc_GR__META_GRP_drug_001_x_drug_021_GR_GR_gDR_bliss_score.xlsx")


