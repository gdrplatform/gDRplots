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

  for (dr_var in intersect(c("DrugName", "Gnumber", "drug_moa", "Concentration"), colnames(dt_))) {
    dt_[[paste0("temp_", dr_var)]] <- dt_[[paste0(dr_var, "_2")]]
    dt_[[paste0(dr_var, "_2")]] <- dt_[[dr_var]] 
    dt_[[dr_var]] <- dt_[[paste0("temp_", dr_var)]]
    dt_[[paste0("temp_", dr_var)]] <- NULL
  }
  dt_
}
