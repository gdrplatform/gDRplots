#' Prep table with metric values for single-agent experiment
#' 
#' @param dt_metrics \code{data.table} representing data from the \code{Metrics} assay,
#'  outputted by \code{gDRutils::convert_se_assay_to_dt(se, "Metrics")}
#'  and single-agent \code{SummarizedExperiment}
#' @param d_name string with drug name to be plotted (identifiers \code{DrugName})
#' @param normalization_type string with normalization types to be selected
#'                           one of: "GR" ("GRvalue") or "RV" ("RelativeViability")
#' @param metric string name of metric;
#'  one of: "xc50" ("GR50" or "IC50" - respectively depending on \code{normalization_type}), 
#'  "x_max" ("GR Max" or "E Max") or "x_mean" ("GR Mean" or "RV Mean")
#' @param fit_source string source name for metrics
#' 
#' @return \code{data.table} with selected metric, input to \code{\link[gDRplots]{prep_dt_assoc}}
#' @keywords prism_plots
#' 
#' @examples
#' mae <- gDRutils::get_synthetic_data("combo_matrix_small")
#' se <- mae[[gDRutils::get_supported_experiments("sa")]]
#' dt_metrics <- gDRutils::convert_se_assay_to_dt(se = se,
#'                                                assay_name = "Metrics")
#' d_name <- "drug_004"
#' dt_response <- .prep_dt_response_metric_sa(dt_metrics, d_name)
#' dt_response <- .prep_dt_response_metric_sa(dt_metrics, d_name,
#'                                            metric = c("xc50", "x_mean", "x_max"))
#' 
#' @export
.prep_dt_response_metric_sa <- function(dt_metrics,
                                        d_name,
                                        normalization_type = "RV",
                                        metric = "xc50",
                                        fit_source = "gDR"
                                        
) {
  
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  clid <- gDRutils::get_env_identifiers("cellline")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  gnumber <- gDRutils::get_env_identifiers("drug")
  
  checkmate::assert_data_table(dt_metrics)
  checkmate::assert_string(d_name)
  checkmate::assert_choice(d_name, choices = dt_metrics[[drug_name]])
  checkmate::assert_choice(normalization_type, choices = c("GR", "RV"))
  checkmate::assert_character(metric, any.missing = FALSE)
  checkmate::assert_subset(metric, choices = c("xc50", "x_mean", "x_max"), empty.ok = FALSE)
  checkmate::assert_string(fit_source, null.ok = TRUE)
  
  # select data for normalization type
  filter_expr <- substitute(normalization_type == norm_type & fit_source == fit_src,
                            list(norm_type = normalization_type, fit_src = fit_source))
  dt_response_metric <- dt_metrics[eval(filter_expr)]
  
  # select required columns
  dt_response_metric <- dt_response_metric[get(drug_name) == d_name, ]
  
  # take care of Inf and NaN values in IC50 metrics
  if (any(metric == "xc50")) {
    inf_xc50 <- is.infinite(dt_response_metric[["xc50"]]) # TODO check: Inf & -Inf
    if (any(inf_xc50)) {
      dt_response_metric[inf_xc50, ][["xc50"]] <- 10^dt_response_metric[inf_xc50, ][["maxlog10Concentration"]]
      # check whether all metric are below 10 ^ maxlog10Concentration
      over_xc50 <- dt_response_metric[["xc50"]] > 10^dt_response_metric[["maxlog10Concentration"]]
      if (any(over_xc50)) {
        dt_response_metric[over_xc50, ][["xc50"]] <- 10^dt_response_metric[over_xc50, ][["maxlog10Concentration"]]
      }
    }
  }
  
  # final
  meta_col <- c(cellline_name, clid, drug_name, gnumber)
  dt_response_metric <- dt_response_metric[, c(meta_col, metric), with = FALSE]
  data.table::setnames(dt_response_metric, metric, sprintf("%s_%s_%s", normalization_type, fit_source, metric))
}

#' Prep table with metric values by doses for single-agent experiment
#' 
#' @param dt_average  \code{data.table} representing data from the \code{Averaged} assay,
#'  outputted by \code{gDRutils::convert_se_assay_to_dt(se, "Averaged")}
#'  and \code{SummarizedExperiment} with chosen data type: single-agent or combo
#' @param d_name string with drug name to be plotted (identifiers \code{DrugName})
#' @param normalization_type string with normalization types to be selected
#'                           one of: "GR" ("GRvalue") or "RV" ("RelativeViability")
#' @param metric string name of metric;
#'    one of: "x" (value of "GR" or "RV" itself - respectively depending on \code{normalization_type}),
#'    or "x_std" (standard deviation)
#' @param fit_source string source name for metrics
#' 
#' @return \code{data.table} with selected metric, input to \code{\link[gDRplots]{prep_dt_assoc}}
#' @keywords prism_plots
#' 
#' @examples
#' mae <- gDRutils::get_synthetic_data("combo_matrix_small")
#' se <- mae[[gDRutils::get_supported_experiments("sa")]]
#' dt_average <- gDRutils::convert_se_assay_to_dt(se = se,
#'                                                assay_name = "Averaged")
#' d_name <- "drug_004"
#' dt_response <- .prep_dt_response_dose_sa(dt_average, d_name)
#' 
#' @export
.prep_dt_response_dose_sa <- function(dt_average,
                                      d_name,
                                      normalization_type = "RV",
                                      metric = "x",
                                      fit_source = "gDR"
                                      
) {
  # TODO add ls_conc -> user can selec conc
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  clid <- gDRutils::get_env_identifiers("cellline")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  gnumber <- gDRutils::get_env_identifiers("drug")
  conc <- gDRutils::get_env_identifiers("concentration")
  
  checkmate::assert_data_table(dt_average)
  checkmate::assert_string(d_name)
  checkmate::assert_choice(d_name, choices = dt_average[[drug_name]])
  checkmate::assert_choice(normalization_type, choices = c("GR", "RV"))
  checkmate::assert_choice(metric, choices = c("x", "x_std"))
  checkmate::assert_string(fit_source, null.ok = TRUE)
  
  # select data for normalization type
  filter_expr <- substitute(normalization_type == norm_type & fit_source == fit_src,
                            list(norm_type = normalization_type, fit_src = fit_source))
  dt_response_dose <- dt_average[eval(filter_expr)]
  
  # select required columns
  dt_response_dose <- dt_response_dose[get(drug_name) == d_name, ]
  
  dt_response_dose_fin <- data.table::dcast(
    data = dt_response_dose,
    formula = get(cellline_name) ~ get(conc),
    value.var = metric,
    fill = NA
  )
  data.table::setkey(dt_response_dose_fin, NULL)
  ls_con <- names(dt_response_dose_fin)[names(dt_response_dose_fin) != "cellline_name"]
  data.table::setnames(dt_response_dose_fin, names(dt_response_dose_fin), 
                       c(cellline_name, sprintf("%s_%s_%s", normalization_type, fit_source, ls_con)))
  
  # final
  meta_col <- c(cellline_name, clid, drug_name, gnumber)
  unique(dt_response_dose[, meta_col, with = FALSE])[dt_response_dose_fin, on = cellline_name]
}


#' Prep table with metric values for combination experiment
#' 
#' @param dt_scores \code{data.table} representing data from the \code{scores} assay,
#'  outputted by \code{gDRutils::convert_se_assay_to_dt(se, "scores")}
#'  and combo \code{SummarizedExperiment}
#' @param d_name string with drug name to be plotted (identifiers \code{DrugName})
#' @param normalization_type string with normalization types to be selected
#'                           one of: "GR" ("GRvalue") or "RV" ("RelativeViability")
#' @param metric string name of combo metric;
#'   one of: "hsa_score"("Bliss Excess GR" or "Bliss Excess RV" - respectively 
#'   depending on \code{normalization_type}), "bliss_score" ("Bliss Score GR" or "Bliss Score RV")
#' @param fit_source string source name for metrics
#' 
#' @return \code{data.table} with selected metric, input to \code{\link[gDRplots]{prep_dt_assoc}}
#' @keywords prism_plots
#' 
#' @examples
#' mae <- gDRutils::get_synthetic_data("combo_matrix_small")
#' se <- mae[[gDRutils::get_supported_experiments("combo")]]
#' dt_scores <- gDRutils::convert_se_assay_to_dt(se = se,
#'                                               assay_name = "scores")
#' d_name <- "drug_004"
#' dt_response <- .prep_dt_response_scores(dt_scores, d_name)
#' dt_response <- .prep_dt_response_scores(dt_scores, d_name,
#'                                         metric = c("hsa_score", "bliss_score"))
#' 
#' @export
.prep_dt_response_scores <- function(dt_scores,
                                     d_name,
                                     normalization_type = "RV",
                                     metric = "hsa_score",
                                     fit_source = "gDR"
                                     
) {
  
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  clid <- gDRutils::get_env_identifiers("cellline")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  gnumber <- gDRutils::get_env_identifiers("drug")
  drug_name_2 <- gDRutils::get_env_identifiers("drug_name2")
  gnumber_2 <- gDRutils::get_env_identifiers("drug2")
  
  checkmate::assert_data_table(dt_scores)
  checkmate::assert_string(d_name)
  checkmate::assert_choice(d_name, choices = dt_scores[[drug_name]])
  checkmate::assert_choice(normalization_type, choices = c("GR", "RV"))
  checkmate::assert_character(metric, any.missing = FALSE)
  checkmate::assert_subset(metric, choices = c("hsa_score", "bliss_score"), empty.ok = FALSE)
  checkmate::assert_string(fit_source, null.ok = TRUE)
  
  # select data for normalization type
  filter_expr <- substitute(normalization_type == norm_type & fit_source == fit_src,
                            list(norm_type = normalization_type, fit_src = fit_source))
  dt_response_scores <- dt_scores[eval(filter_expr)]
  
  # select required columns
  dt_response_scores <- dt_response_scores[get(drug_name) == d_name, ]
  
  # final
  meta_col <- c(cellline_name, clid, drug_name, gnumber, drug_name_2, gnumber_2)
  dt_response_scores <- dt_response_scores[, c(meta_col, metric), with = FALSE]
  data.table::setnames(dt_response_scores, metric, sprintf("%s_%s_%s", normalization_type, fit_source, metric))
}


#' Prep table with calculated linear associations
#'
#' @param dt_response \code{data.table} with experimental response data (rows are samples).
#' @param dt_depmap \code{data.table} with dependent variables data matrix 
#'    (rows are samples, columns are features or meta).
#'
#' @return \code{data.table} with selected metric, input to \code{\link[gDRplots]{plot_volcano_assoc}}
#' @keywords prism_plots
#' 
prep_dt_assoc <- function(dt_response,
                          dt_depmap) {
  
  checkmate::assert_data_table(dt_response)
  checkmate::assert_data_table(dt_depmap)
  checkmate::assert_names(names(dt_depmap), must.include = "CCLEName")
  
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  clid <- gDRutils::get_env_identifiers("cellline")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  gnumber <- gDRutils::get_env_identifiers("drug")
  drug_name_2 <- gDRutils::get_env_identifiers("drug_name2")
  gnumber_2 <- gDRutils::get_env_identifiers("drug2")
  CCLEName <- NULL # due to NSE notes in R CMD check
  
  # shared cell line
  depmap_lines <- dt_depmap[CCLEName != "", unique(CCLEName)]
  response_lines <- dt_response[[cellline_name]]
  shared_lines <- intersect(depmap_lines, response_lines)
  
  # subset the data.table and order it
  X_dt <- dt_depmap[CCLEName %in% shared_lines]
  data.table::setorder(X_dt, "CCLEName")
  Y_dt <- dt_response[get(cellline_name) %in% shared_lines]
  data.table::setorder(Y_dt, cellline_name)
  
  # convert to a matrix
  x_col <- setdiff(names(X_dt), c("ModelID", "CCLEName"))
  X <- as.matrix(
    X_dt[, .SD, .SDcols = c("CCLEName", x_col)]
    , rownames = "CCLEName"
  )
  y_col <- setdiff(names(Y_dt), c(cellline_name, clid, drug_name, gnumber, drug_name_2, gnumber_2))
  Y <- as.matrix(
    Y_dt[, .SD, .SDcols = c("CellLineName", y_col)]
    , rownames = "CellLineName"
  )
  
  # # create dt_assoc # nolint start WIP
  # dt_assoc <- kaleidoscope::calc_assoc(X, Y)
  # return(dt_assoc) # nolint end
}