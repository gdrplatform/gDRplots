#' Help function to change column suffix
#'
#' @param dt_ data.table with data
#'
#' @keywords internal
#' @examples
#' \dontrun{
#' mae <- gDRutils::get_synthetic_data("small")
#' se <- mae[[1]]
#' dt_ <- gDRutils::convert_se_assay_to_dt(se, assay_name = "Metrics", wide_structure = TRUE)
#' 
#' swapped_dt_ <- swap_drugs_1_2(dt_)
#' }
#' 
#' @return swapped data.table
#' 
swap_drugs_1_2 <- function(dt_) {
  
  checkmate::assert_data_table(dt_)
  
  conc <- gDRutils::get_env_identifiers("concentration")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  drug <- gDRutils::get_env_identifiers("drug")
  drug_moa <- gDRutils::get_env_identifiers("drug_moa")

  for (dr_var in intersect(c(drug_name, drug, drug_moa, conc), colnames(dt_))) {
    dt_[[paste0("temp_", dr_var)]] <- dt_[[paste0(dr_var, "_2")]]
    dt_[[paste0(dr_var, "_2")]] <- dt_[[dr_var]]
    dt_[[dr_var]] <- dt_[[paste0("temp_", dr_var)]]
    dt_[[paste0("temp_", dr_var)]] <- NULL
  }
  dt_
}
