testdata_dir <- feat_data_path <- system.file("testdata", package = "gDRplots")

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

data.table::fwrite(tab_model, 
                   file = file.path(testdata_dir, "Model.csv"),
                   row.names = FALSE)


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

data.table::fwrite(tab_cko, 
                   file = file.path(testdata_dir, "CRISPRGeneEffect.csv"), 
                   row.names = FALSE)

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

data.table::fwrite(tab_hot, 
                   file = file.path(testdata_dir, "OmicsSomaticMutationsMatrixHotspot.csv"), 
                   row.names = FALSE)

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

data.table::fwrite(tab_profil, 
                   file = file.path(testdata_dir, "OmicsSignaturesProfile.csv"), 
                   row.names = FALSE)


# OmicsArmLevelCNA ----
tab_arm <- data.table::data.table(
  V1 = tab_model$ModelID,
  `1p` = withr::with_seed(42, sample(c(-1, 0, 1), size = NROW(tab_model), prob = c(0.20, 0.75, 0.05), replace = TRUE)),
  `1q` = withr::with_seed(42, sample(c(-1, 0, 1), size = NROW(tab_model), prob = c(0.08, 0.8, 0.12), replace = TRUE)),
  `8p` = withr::with_seed(314, sample(c(-1, 0, 1), size = NROW(tab_model), prob = c(0.20, 0.75, 0.05), replace = TRUE)),
  `8q` = withr::with_seed(314, sample(c(-1, 0, 1), size = NROW(tab_model), prob = c(0.08, 0.8, 0.12), replace = TRUE)),
  `16p` = withr::with_seed(42, sample(c(-1, 0, 1), size = NROW(tab_model), prob = c(0.15, 0.70, 0.15), replace = TRUE)),
  `16q` = withr::with_seed(42, sample(c(-1, 0, 1), size = NROW(tab_model), prob = c(0.09, 0.7, 0.21), replace = TRUE)),
  `22p` = withr::with_seed(314, sample(c(-1, 0, 1), size = NROW(tab_model), prob = c(0.15, 0.70, 0.15), replace = TRUE)),
  `22q` = withr::with_seed(314, sample(c(-1, 0, 1), size = NROW(tab_model), prob = c(0.09, 0.7, 0.21), replace = TRUE))
)

data.table::fwrite(tab_arm, 
                   file = file.path(testdata_dir, "OmicsArmLevelCNA"), 
                   row.names = FALSE)

# assoc table ----
tab_drug_001_met1_RV <- data.table::data.table(
  feature = sprintf("NU_%03d_X1%s", 1:25, LETTERS[1:25]),
  response = rep("RV_gDR_x_max", 25),
  rho = withr::with_seed(42, sample(seq(-1, 1, 0.35), 25, replace = TRUE)),
  q_value = withr::with_seed(314, sample(seq(0.001, 0.1, 0.0025), 25, replace = TRUE))
)
tab_drug_001_met1_RV$neglog_q_value <- -log10(tab_drug_001_met1_RV$q_value)

writexl::write_xlsx(
  tab_drug_001_met1_RV, 
  path = file.path(testdata_dir, "tab_assoc_RV__featNUX_drug_001_RV_gDR_x_max.xlsx")
)

tab_drug_001_met2_RV <- data.table::data.table(
  feature = sprintf("NU_%03d_X1%s", 1:25, LETTERS[1:25]),
  response = rep("RV_gDR_x_mean", 25),
  rho = withr::with_seed(42, sample(seq(-0.85, 1.55, 0.35), 25, replace = TRUE)),
  q_value = withr::with_seed(42, sample(seq(0.001, 0.55, 0.0015), 25, replace = TRUE))
)
tab_drug_001_met2_RV$neglog_q_value <- -log10(tab_drug_001_met2_RV$q_value)

writexl::write_xlsx(
  tab_drug_001_met2_RV, 
  path = file.path(testdata_dir, "tab_assoc_RV__featNUX_drug_001_RV_gDR_x_mean.xlsx")
)


tab_drug_002_met2_RV <- data.table::data.table(
  feature = sprintf("GRP_%03d_XC", 1:25),
  response = rep("RV_gDR_x_mean", 25),
  rho = withr::with_seed(314, sample(seq(-0.85, 1.55, 0.35), 25, replace = TRUE)),
  q_value = withr::with_seed(42, sample(seq(0.001, 0.55, 0.0025), 25, replace = TRUE))
)
tab_drug_002_met2_RV$neglog_q_value <- -log10(tab_drug_002_met2_RV$q_value)

writexl::write_xlsx(
  tab_drug_002_met2_RV, 
  "./gDRplots/inst/testdata/tab_assoc_RV__metaGRP_drug_002_RV_gDR_x_mean.xlsx")

tab_drug_001_met1_GR <- data.table::data.table(
  feature = sprintf("GRP_%03d_XC", 1:25),
  response = rep("GR_gDR_x_max", 25),
  rho = withr::with_seed(42, sample(seq(-1, 1, 0.35), 25, replace = TRUE)),
  q_value = withr::with_seed(314, sample(seq(0.001, 0.1, 0.0025), 25, replace = TRUE))
)
tab_drug_001_met1_GR$neglog_q_value <- -log10(tab_drug_001_met1_GR$q_value)

writexl::write_xlsx(
  tab_drug_001_met1_GR, 
  path = file.path(testdata_dir, "tab_assoc_GR__metaGRP_drug_001_GR_gDR_x_max.xlsx")
)

tab_drug_001_met4_GR <- data.table::data.table(
  feature = sprintf("GRP_%03d_XC", 1:25),
  response = rep("GR_gDR_log10_xc50", 25),
  rho = withr::with_seed(314, sample(seq(-1, 1, 0.35), 25, replace = TRUE)),
  q_value = withr::with_seed(314, sample(seq(0.001, 0.55, 0.0015), 25, replace = TRUE))
)
tab_drug_001_met4_GR$neglog_q_value <- -log10(tab_drug_001_met4_GR$q_value)

writexl::write_xlsx(
  tab_drug_001_met4_GR, 
  path = file.path(testdata_dir, "tab_assoc_GR__metaGRP_drug_001_GR_gDR_log10_xc50.xlsx")
)

tab_drug_003_met4_GR <- data.table::data.table(
  feature = sprintf("NU_%03d_X1%s", 1:25, LETTERS[1:25]),
  response = rep("GR_gDR_log10_xc50", 25),
  rho = withr::with_seed(42, sample(seq(-0.85, 1.55, 0.35), 25, replace = TRUE)),
  q_value = withr::with_seed(42, sample(seq(0.001, 0.55, 0.0015), 25, replace = TRUE))
)
tab_drug_003_met4_GR$neglog_q_value <- -log10(tab_drug_003_met4_GR$q_value)

writexl::write_xlsx(
  tab_drug_003_met4_GR, 
  path = file.path(testdata_dir, "tab_assoc_GR__featNUX_drug_003_GR_gDR_log10_xc50.xlsx")
)
