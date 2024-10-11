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
#' dt_response <- prep_dt_response_metric_sa(dt_metrics, d_name)
#' dt_response <-
#'   prep_dt_response_metric_sa(dt_metrics, d_name,
#'                              metric = c("xc50", "x_mean", "x_max"))
#' 
#' @export
prep_dt_response_metric_sa <- function(dt_metrics,
                                       d_name,
                                       normalization_type = "RV",
                                       metric = "xc50",
                                       fit_source = "gDR") {
  
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  
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
  
  # select required drug
  dt_response_metric <- dt_response_metric[get(drug_name) == d_name, ]
  
  # take care of Inf and NaN values in IC50 metrics
  if (any(metric == "xc50")) {
    inf_xc50 <- is.infinite(dt_response_metric[["xc50"]]) # TODO check: Inf & -Inf
    if (any(inf_xc50, na.rm = TRUE)) {
      dt_response_metric[inf_xc50, ][["xc50"]] <- 10^dt_response_metric[inf_xc50, ][["maxlog10Concentration"]]
      # check whether all metric are below 10 ^ maxlog10Concentration
      over_xc50 <- dt_response_metric[["xc50"]] > 10^dt_response_metric[["maxlog10Concentration"]]
      if (any(over_xc50, na.rm = TRUE)) {
        dt_response_metric[over_xc50, ][["xc50"]] <- 10^dt_response_metric[over_xc50, ][["maxlog10Concentration"]]
      }
    }
  }
  
  # final
  meta_col <- c("rId", "cId", cellline_name)
  dt_response_metric <- dt_response_metric[, c(meta_col, metric), with = FALSE]
  data.table::setnames(dt_response_metric, 
                       old = metric, 
                       new = sprintf("%s_%s_%s", normalization_type, fit_source, metric))
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
#' dt_response <- prep_dt_response_dose_sa(dt_average, d_name)
#' 
#' @export
prep_dt_response_dose_sa <- function(dt_average,
                                     d_name,
                                     normalization_type = "RV",
                                     metric = "x",
                                     fit_source = "gDR") {
  # TODO add ls_conc -> user can select conc
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
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
  
  # select required drug
  dt_response_dose <- dt_response_dose[get(drug_name) == d_name, ]
  
  dt_response_dose_fin <- data.table::dcast(
    data = dt_response_dose,
    formula = get(cellline_name) ~ get(conc),
    value.var = metric,
    fill = NA
  )
  data.table::setkey(dt_response_dose_fin, NULL)
  ls_conc <- names(dt_response_dose_fin)[names(dt_response_dose_fin) != "cellline_name"]
  data.table::setnames(
    dt_response_dose_fin, 
    old = names(dt_response_dose_fin), 
    new = c(cellline_name, sprintf("%s_%s_%s_%s", normalization_type, fit_source, metric, ls_conc)))
  
  # final
  meta_col <- c("rId", "cId", cellline_name)
  dt_response_dose_fin <- 
    unique(dt_response_dose[, meta_col, with = FALSE])[dt_response_dose_fin, on = cellline_name]
  dt_response_dose_fin
}

#' Prep table with metric values for single-agent experiment
#' 
#' @param dt_average  \code{data.table} representing data from the \code{Averaged} assay,
#'  outputted by \code{gDRutils::convert_se_assay_to_dt(se, "Averaged")}
#'  and \code{SummarizedExperiment} with chosen data type: single-agent or combo
#' @param dt_metrics \code{data.table} representing data from the \code{Metrics} assay,
#'  outputted by \code{gDRutils::convert_se_assay_to_dt(se, "Metrics")}
#'  and single-agent \code{SummarizedExperiment}
#' @param d_name string with drug name to be plotted (identifiers \code{DrugName})
#' @param normalization_type string with normalization types to be selected
#'                           one of: "GR" ("GRvalue") or "RV" ("RelativeViability")
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
#' dt_metrics <- gDRutils::convert_se_assay_to_dt(se = se,
#'                                                assay_name = "Metrics")
#' d_name <- "drug_004"
#' dt_response_sa <- prep_dt_response_sa(dt_average, dt_metrics, d_name)
#' 
#' @export
prep_dt_response_sa <- function(dt_average,
                                dt_metrics,
                                d_name,
                                normalization_type = "RV",
                                fit_source = "gDR") {
  
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  
  checkmate::assert_data_table(dt_average)
  checkmate::assert_data_table(dt_metrics)
  checkmate::assert_string(d_name)
  checkmate::assert_choice(d_name, choices = dt_average[[drug_name]])
  checkmate::assert_choice(d_name, choices = dt_metrics[[drug_name]])
  checkmate::assert_choice(normalization_type, choices = c("GR", "RV"))
  checkmate::assert_string(fit_source, null.ok = TRUE)
  
  dt_response_met <-
    prep_dt_response_metric_sa(dt_metrics = dt_metrics,
                               d_name = d_name,
                               normalization_type = normalization_type,
                               metric = c("xc50", "x_mean", "x_max"))
  dt_response_dose <- 
    prep_dt_response_dose_sa(dt_average = dt_average, 
                             d_name = d_name,
                             normalization_type = normalization_type)
  
  id_col <- c("rId", "cId", cellline_name)
  dt_response_sa <- merge(dt_response_met, dt_response_dose, by = id_col)
  return(dt_response_sa)
}

#' Prep table with metric values for combination experiment
#' 
#' @param dt_scores \code{data.table} representing data from the \code{scores} assay,
#'  outputted by \code{gDRutils::convert_se_assay_to_dt(se, "scores")}
#'  and combo \code{SummarizedExperiment}
#' @param d_name string with drug name to be plotted (identifiers \code{DrugName})
#' @param d_name2 string with drug name to be plotted (identifiers \code{DrugName_2})
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
#' d_name2 <- "drug_026"
#' dt_response <- prep_dt_response_scores(dt_scores, d_name, d_name2)
#' dt_response <- 
#'   prep_dt_response_scores(dt_scores, d_name, d_name2,
#'                           metric = c("hsa_score", "bliss_score"))
#' 
#' @export
prep_dt_response_scores <- function(dt_scores,
                                    d_name,
                                    d_name2,
                                    normalization_type = "RV",
                                    metric = "hsa_score",
                                    fit_source = "gDR") {
  
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  drug_name_2 <- gDRutils::get_env_identifiers("drug_name2")
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  
  checkmate::assert_data_table(dt_scores)
  checkmate::assert_string(d_name)
  checkmate::assert_choice(d_name, choices = dt_scores[[drug_name]])
  checkmate::assert_string(d_name2)
  checkmate::assert_choice(d_name2, choices = dt_scores[[drug_name_2]])
  checkmate::assert_choice(normalization_type, choices = c("GR", "RV"))
  checkmate::assert_character(metric, any.missing = FALSE)
  checkmate::assert_subset(metric, choices = c("hsa_score", "bliss_score"), empty.ok = FALSE)
  checkmate::assert_string(fit_source, null.ok = TRUE)
  
  # select data for normalization type
  filter_expr <- substitute(normalization_type == norm_type & fit_source == fit_src,
                            list(norm_type = normalization_type, fit_src = fit_source))
  dt_response_scores <- dt_scores[eval(filter_expr)]
  
  # select required drugs combination
  dt_response_scores <- dt_response_scores[get(drug_name) == d_name & get(drug_name_2) == d_name2, ]
  
  # final
  meta_col <- c("rId", "cId", cellline_name)
  dt_response_scores <- dt_response_scores[, c(meta_col, metric), with = FALSE]
  data.table::setnames(dt_response_scores, old = metric, 
                       new = sprintf("%s_%s_%s", normalization_type, fit_source, metric))
}


#' Prep table with metric values for combination experiment
#' 
#' @param dt_metrics \code{data.table} representing data from the \code{Metrics} assay,
#'  outputted by \code{gDRutils::convert_se_assay_to_dt(se, "Metrics")}
#'  and combo \code{SummarizedExperiment}
#' @param d_name string with drug name to be plotted (identifiers \code{DrugName})
#' @param d_name2 string with drug name to be plotted (identifiers \code{DrugName_2})
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
#' dt_metrics <- gDRutils::convert_se_assay_to_dt(se = se,
#'                                                assay_name = "Metrics")
#' d_name <- "drug_004"
#' d_name2 <- "drug_026"
#' dt_response <- prep_dt_response_metric_diff(dt_metrics, d_name, d_name2)
#' dt_response <- 
#'   prep_dt_response_metric_diff(dt_metrics, d_name, d_name2,
#'                                metric = c("xc50", "x_mean", "x_max"))
#' 
#' @export
prep_dt_response_metric_diff <- function(dt_metrics,
                                         d_name,
                                         d_name2,
                                         normalization_type = "RV",
                                         metric = "xc50",
                                         fit_source = "gDR") {
  
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  drug_name_2 <- gDRutils::get_env_identifiers("drug_name2")
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  
  checkmate::assert_data_table(dt_metrics)
  checkmate::assert_string(d_name)
  checkmate::assert_choice(d_name, choices = dt_metrics[[drug_name]])
  checkmate::assert_string(d_name2)
  checkmate::assert_choice(d_name2, choices = dt_metrics[[drug_name_2]])
  checkmate::assert_choice(normalization_type, choices = c("GR", "RV"))
  checkmate::assert_character(metric, any.missing = FALSE)
  checkmate::assert_subset(metric, choices = c("xc50", "x_mean", "x_max"), empty.ok = FALSE)
  checkmate::assert_string(fit_source, null.ok = TRUE)
  
  # select data for normalization type
  filter_expr <- substitute(normalization_type == norm_type & fit_source == fit_src,
                            list(norm_type = normalization_type, fit_src = fit_source))
  dt_response_metric <- dt_metrics[eval(filter_expr)]
  
  # select required drugs combination
  dt_response_metric <- dt_response_metric[get(drug_name) == d_name & get(drug_name_2) == d_name2, ]
  
  # take care of Inf and NaN values in IC50 metrics
  if (any(metric == "xc50")) {
    inf_xc50 <- is.infinite(dt_response_metric[["xc50"]]) # TODO check: Inf & -Inf
    if (any(inf_xc50, na.rm = TRUE)) {
      dt_response_metric[inf_xc50, ][["xc50"]] <- 10^dt_response_metric[inf_xc50, ][["maxlog10Concentration"]]
      # check whether all metric are below 10 ^ maxlog10Concentration
      over_xc50 <- dt_response_metric[["xc50"]] > 10^dt_response_metric[["maxlog10Concentration"]]
      if (any(over_xc50, na.rm = TRUE)) {
        dt_response_metric[over_xc50, ][["xc50"]] <- 10^dt_response_metric[over_xc50, ][["maxlog10Concentration"]]
      }
    }
  }
  
  # create entries of non-zero co-trt
  meta_col <- c("rId", "cId", cellline_name)
  ls_cols <- c(meta_col, "cotrt_value", "source", metric)
  dt_non_zero <- data.table::copy(dt_response_metric)[cotrt_value != 0, .SD, .SDcols = ls_cols]
  data.table::setnames(dt_non_zero, metric, paste0(metric, "_cotrt"))
  
  # create entries of zero co-trt (single agent)
  dt_zero <- data.table::copy(dt_response_metric)[cotrt_value == 0, .SD, .SDcols = ls_cols]
  data.table::setnames(dt_zero, 
                       old = c("cotrt_value", metric), 
                       new = c("cotrt_value_zero", paste0(metric, "_cotrt_zero")))
  
  # merge zero and non zero
  dt_combo_merged <- dt_zero[dt_non_zero, on = c(meta_col, "source"), nomatch = NULL]
  dt_combo_merged[, cotrt_value_zero := NULL]
  
  # calculate differences
  dt_combo_diff <- 
    dt_combo_merged[, (paste0(metric, "_cotrt_diff")) := Map("-", 
                                                             mget(paste0(metric, "_cotrt")), 
                                                             mget(paste0(metric, "_cotrt_zero")))]
  ls_col_met <- 
    colnames(dt_combo_diff)[!colnames(dt_combo_diff) %in% c(meta_col, "cotrt_value", "source")]
  ls_col_met_fin <- sprintf("%s_%s_%s", normalization_type, fit_source, ls_col_met)
  data.table::setnames(dt_combo_diff, ls_col_met, ls_col_met_fin)
  
  # final
  dt_combo_diff <- data.table::dcast(
    data = dt_combo_diff, 
    formula = rId + cId + get(cellline_name) ~ cotrt_value + source, 
    value.var = ls_col_met_fin)
  data.table::setnames(dt_combo_diff, "cellline_name", cellline_name)
  data.table::setkey(dt_combo_diff, NULL)
  (dt_combo_diff)
}

#' Prep table with metric values for combination experiment
#' 
#' @param dt_metrics \code{data.table} representing data from the \code{Metrics} assay,
#'  outputted by \code{gDRutils::convert_se_assay_to_dt(se, "Metrics")}
#'  and combo \code{SummarizedExperiment}
#' @param dt_scores \code{data.table} representing data from the \code{scores} assay,
#'  outputted by \code{gDRutils::convert_se_assay_to_dt(se, "scores")}
#'  and combo \code{SummarizedExperiment}
#' @param d_name string with drug name to be plotted (identifiers \code{DrugName})
#' @param d_name2 string with drug name to be plotted (identifiers \code{DrugName_2})
#' @param normalization_type string with normalization types to be selected
#'                           one of: "GR" ("GRvalue") or "RV" ("RelativeViability")
#' @param fit_source string source name for metrics
#' 
#' @return \code{data.table} with selected metric, input to \code{\link[gDRplots]{prep_dt_assoc}}
#' @keywords prism_plots
#' 
#' @examples
#' mae <- gDRutils::get_synthetic_data("combo_matrix_small")
#' se <- mae[[gDRutils::get_supported_experiments("combo")]]
#' dt_metrics <- gDRutils::convert_se_assay_to_dt(se = se,
#'                                                assay_name = "Metrics")
#' dt_scores <- gDRutils::convert_se_assay_to_dt(se = se,
#'                                               assay_name = "scores")
#' d_name <- "drug_004"
#' d_name2 <- "drug_026"
#' dt_response_combo <-
#'   prep_dt_response_combo(dt_metrics, dt_scores, d_name, d_name2)
#' 
#' @export
prep_dt_response_combo <- function(dt_metrics,
                                   dt_scores,
                                   d_name,
                                   d_name2,
                                   normalization_type = "RV",
                                   fit_source = "gDR") {
  
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  drug_name_2 <- gDRutils::get_env_identifiers("drug_name2")
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  
  checkmate::assert_data_table(dt_metrics)
  checkmate::assert_data_table(dt_scores)
  checkmate::assert_string(d_name)
  checkmate::assert_choice(d_name, choices = dt_metrics[[drug_name]])
  checkmate::assert_choice(d_name, choices = dt_scores[[drug_name]])
  checkmate::assert_string(d_name2)
  checkmate::assert_choice(d_name2, choices = dt_metrics[[drug_name_2]])
  checkmate::assert_choice(d_name2, choices = dt_scores[[drug_name_2]])
  checkmate::assert_choice(normalization_type, choices = c("GR", "RV"))
  checkmate::assert_string(fit_source, null.ok = TRUE)
  
  dt_response_scores <-
    prep_dt_response_scores(dt_scores = dt_scores,
                            d_name = d_name,
                            d_name2 = d_name2,
                            normalization_type = normalization_type,
                            metric = c("hsa_score", "bliss_score"))
  
  dt_response_met_diff <- 
    prep_dt_response_metric_diff(dt_metrics = dt_metrics, 
                                 d_name = d_name,
                                 d_name2 = d_name2,
                                 normalization_type = normalization_type,
                                 metric = c("xc50", "x_mean", "x_max"))
  
  id_col <- c("rId", "cId", cellline_name)
  dt_response_combo <- merge(dt_response_scores, dt_response_met_diff, by = id_col)
  return(dt_response_combo)
}

#' Load DepMap merged data for one selected feature
#'
#' @param feature_set string name of the molecular feature set to load from DepMap.
#' @param prefix string prefixes to use for the each feature set in \code{feature_set};
#'    has to be the same length as \code{feature_sets}
#'
#' @return A named list with elements, that may be input to \code{\link[gDRplots]{prep_dt_assoc}}
#' \itemize{
#'   \item \code{dt_depmap} \code{data.table} with feature data from DepMap (wide format),.
#'   \item \code{selected_feat_meta_col} string name of feature.
#' }
#' 
#' @keywords internal
#'
#' @seealso \code{kaleidoscope::load_depmap_merged}
#'
#' @examples
#' \dontrun{
#' dt_depmap_feat <- prep_dt_depmap_feat() 
#' }
#' 
#' @export
prep_dt_depmap_feat <- function(
    feature_set = "CRISPRGeneEffect",
    prefix = "KO_") {
  
  checkmate::assert_string(feature_set)
  checkmate::assert_string(prefix)
  
  stopifnot("`prefix` has to be the same length as `feature_sets`" = NROW(feature_set) == NROW(prefix))
  
  # TODO in GDR-2710 # nolint start
  # dt_depmap <- kaleidoscope::load_depmap_merged(
  #   feature_sets = feature_set,
  #   prefix = prefix,
  #   metadata_columns = "CCLEName") 
  # 
  # data.table::setkey(dt_depmap, NULL)
  # dt_depmap["CCLEName" != ""]
  # 
  # return(list(dt_depmap = dt_depmap, selected_feat_meta_col = feature_set)) # nolint end
}

#' Load DepMap merged data for one selected metadata
#'
#' @param metadata_col character vector with the metadata columns to load for DepMap cell lines
#'
#' @return A named list with elements, that may be input to \code{\link[gDRplots]{prep_dt_assoc}}
#' \itemize{
#'   \item \code{dt_depmap} \code{data.table} with feature data from DepMap (wide format),
#'   \item \code{selected_feat_meta_col} string name of metadata column..
#' }
#' 
#' @keywords internal
#'
#' @seealso \code{kaleidoscope::load_depmap_merged}
#'
#' @examples
#' \dontrun{
#' dt_depmap_meta <- prep_dt_depmap_meta() 
#' }
#' @export
prep_dt_depmap_meta <- function(metadata_col = "OncotreeLineage") {
  
  checkmate::assert_string(metadata_col)
  
  # TODO in GDR-2710 # nolint start
  # ls_depmap <- kaleidoscope::load_depmap_list(
  #   feature_sets = "OmicsCNGene",
  #   prefix = "CN_",
  #   metadata_columns = unique(c(metadata_col, "CCLEName"))) # nolint end
  ls_depmap <- ls_depmap[unique(c(metadata_col, "CCLEName"))]
  
  dt_depmap <- data.table::data.table(
    merge(ls_depmap[["CCLEName"]], ls_depmap[[metadata_col]], by = "row.names", all = "TRUE")
  )
  data.table::setnames(dt_depmap, c("V1", "Row.names"), c("CCLEName", "ModelID"))
  
  data.table::setkey(dt_depmap, NULL)
  dt_depmap["CCLEName" != ""]
  
  return(list(dt_depmap = dt_depmap, selected_feat_meta_col = metadata_col))
}

#' Prep table with calculated linear associations
#'
#' @param dt_response \code{data.table} with experimental response data (rows are samples) for one metric
#' @param dt_depmap \code{data.table} with dependent variables data load from DepMap.
#'   (rows are samples, columns are features or meta);  
#'   outputted by one of \code{\link[gDRplots]{prep_dt_depmap_feat}} or
#'   \code{\link[gDRplots]{prep_dt_depmap_meta}}
#' @param selected_feat_meta_col string name of feature/meta column in DepMap
#'   
#' @return A named list with elements, that may be input to \code{\link[gDRplots]{plot_volcano_assoc}}
#' \itemize{
#'   \item \code{dt_assoc} \code{data.table} with calculated association values between 
#'      feature/meta of DepMap and selected metric,
#'   \item \code{condition_info} string describing experiment condition (drugs),
#'   \item \code{selected_feat_meta_col} string name of feature/meta.
#' }
#' 
#' @keywords prism_plots
#' 
#' @export
prep_dt_assoc <- function(dt_response,
                          dt_depmap,
                          selected_feat_meta_col = NULL) {
  
  checkmate::assert_data_table(dt_response)
  checkmate::assert_data_table(dt_depmap)
  checkmate::assert_names(names(dt_depmap), must.include = "CCLEName")
  checkmate::assert_string(selected_feat_meta_col, null.ok = TRUE)
  
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  
  # checking input format
  selected_metric <- setdiff(names(dt_response), c("rId", "cId", cellline_name))
  stopifnot("Provide `dt_response` with for one metric." = NROW(selected_metric) == 1)
  selected_feat_meta <- setdiff(names(dt_depmap), c("ModelID", "CCLEName"))
  stopifnot("Provide `dt_depmap` for one feature or one meta only." = 
              all(vapply(dt_depmap[, selected_feat_meta, with = FALSE], is.numeric, logical(1))))
  
  # result
  # will be returned in this format when: 
  # 1) length (shared lines) < 6 
  # 2) all/one of the association calculation inputs will have only NA values
  # 3) selected_metric values have no variance (due to cdsrmodels::lin_associations)
  obj_assoc <- list(dt_assoc = data.table::data.table(feature = character(0),
                                                      response = character(0),
                                                      rho = numeric(0),
                                                      q_value = numeric(0)),
                    condition_info = unique(dt_response[["rId"]]),
                    selected_metric = selected_metric,
                    selected_feat_meta_col = selected_feat_meta_col)
  
  # shared cell line
  depmap_lines <- dt_depmap[CCLEName != "", unique(CCLEName)]
  response_lines <- dt_response[[cellline_name]]
  shared_lines <- intersect(depmap_lines, response_lines)
  
  if (NROW(shared_lines) >= 6) { # (the minimum degrees of freedom = 4) + 2
    # subset the data.table and order it
    X_dt <- dt_depmap[CCLEName %in% shared_lines, ]
    data.table::setorder(X_dt, "CCLEName")
    Y_dt <- dt_response[get(cellline_name) %in% shared_lines, ]
    data.table::setorderv(Y_dt, cellline_name)
    
    # association can only be calculated for a value other than NA
    Y_condition <- 
      all(NROW(stats::na.omit(Y_dt)) >= 6 && stats::sd(Y_dt[[selected_metric]], na.rm = TRUE) > 0)
    X_condition <-
      any(vapply(shared_lines, function(nm) {
        stats::sd(X_dt[CCLEName == nm, .SD, 
                       .SDcols = vapply(X_dt, is.numeric, logical(1))], na.rm = TRUE) > 0 && 
          all(!is.na(X_dt[CCLEName == nm]))
      }, logical(1)))
    
    if (Y_condition && X_condition) {
      # convert to a matrix
      X <- as.matrix(
        X_dt[, .SD, .SDcols = c("CCLEName", selected_feat_meta)], rownames = "CCLEName"
      )
      Y <- as.matrix(
        Y_dt[, .SD, .SDcols = c("CellLineName", selected_metric)], rownames = "CellLineName"
      )
      
      # create dt_assoc
      # TODO in GDR-2710
      # dt_assoc <- kaleidoscope::calc_assoc(X, Y)  # nolint start
      # 
      # # final
      # obj_assoc[["dt_assoc"]] <- dt_assoc[, c("feature", "response", "rho", "q_value"), with = FALSE] # nolint end
    }
  }
  # return
  return(obj_assoc)
}
