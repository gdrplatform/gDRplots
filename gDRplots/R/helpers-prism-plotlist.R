#' Create a nested list of plots for PRISM data with single-agent metrics
#' 
#' @param drug_name_vec character vector with drug names to be plotted (identifiers \code{DrugName})
#' @param dt_metrics \code{data.table} representing data from the \code{Metrics} assay,
#'  outputted by \code{gDRutils::convert_se_assay_to_dt(se, "Metrics")}
#'  and single-agent \code{SummarizedExperiment}
#' @param dt_average  \code{data.table} representing data from the \code{Averaged} assay,
#'  outputted by \code{gDRutils::convert_se_assay_to_dt(se, "Averaged")}
#'  and \code{SummarizedExperiment} with chosen data type: single-agent or combo 
#' @param normalization_type_vec character vector with normalization types to be selected
#'                               one of: "GR" ("GRvalue") or "RV" ("RelativeViability") or both
#' @param metric character vector with names of metric;
#'   chosen from: "xc50" ("GR50" or "IC50" - respectively depending on \code{normalization_type}), 
#'  "x_max" ("GR Max" or "E Max") or "x_mean" ("GR Mean" or "RV Mean")
#' @param fit_source string source name for metrics
#' @param meta_data_path string path to metadata file describing all cancer models/cell lines
#'  which are referenced by a dataset contained within the DepMap portal. 
#'  It is usually a file named \code{Model.csv}.
#' @param feature_sets character vector names of the molecular feature set to load from DepMap.
#' @param prefixes character vector prefixes to use for the each feature set in \code{feature_sets};
#'    has to be the same length as \code{feature_sets}
#' @param metadata_columns character vector with the metadata columns to load for DepMap cell lines
#' 
#' @return nested list of plots
#' @keywords prism_plots
#'   
#' @export
create_PRISM_plot_list_sa <- function(drug_name_vec,
                                      dt_metrics,
                                      dt_average,
                                      normalization_type_vec = "RV",
                                      metric = c("xc50", "x_mean", "x_max"),
                                      fit_source = "gDR",
                                      meta_data_path,
                                      feature_sets,
                                      prefixes,
                                      metadata_columns = NULL) {
  
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  id_col <- c("rId", "cId", cellline_name)
  
  checkmate::assert_character(drug_name_vec)
  checkmate::assert_data_table(dt_metrics, min.rows = 1, null.ok = TRUE)
  checkmate::assert_data_table(dt_average, min.rows = 1, null.ok = TRUE)
  checkmate::assert_subset(normalization_type_vec, choices = c("GR", "RV"))
  checkmate::assert_character(metric, any.missing = FALSE, null.ok = TRUE)
  stopifnot("Provide `metric` for `dt_metrics`." = !is.null(dt_metrics) && !is.null(metric))
  if (!is.null(metric)) checkmate::assert_subset(metric, choices = c("xc50", "x_mean", "x_max"))
  if (!is.null(dt_metrics)) checkmate::assert_subset(metric, choices = names(dt_metrics))
  checkmate::assert_string(fit_source, null.ok = TRUE)
  checkmate::assert_string(meta_data_path, pattern = ".*Model.csv$")
  checkmate::assert_file_exists(meta_data_path)
  checkmate::assert_character(feature_sets, null.ok = TRUE)
  checkmate::assert_character(prefixes, null.ok = TRUE)
  stopifnot("`prefixes` has to be the same length as `feature_sets`" = NROW(feature_sets) == NROW(prefixes))
  checkmate::assert_character(metadata_columns, null.ok = TRUE)
  stopifnot("Provide `feature_sets` or `metadata_columns` for DepMam subset." =
              !is.null(feature_sets) || !is.null(metadata_columns))
  
  # adjust drug list
  available_drugs <- unique(dt_metrics[[drug_name]])
  if (is.null(drug_name_vec) || all(!drug_name_vec %in% available_drugs)) {
    drug_name_vec <- available_drugs
  } else if (!all(drug_name_vec %in% available_drugs)) {
    drug_name_vec <- drug_name_vec[drug_name_vec %in% available_drugs]
  }
  
  ls_plot <- list()
  
  if (!is.null(feature_sets)) {
    # 1st level
    for (feat in feature_sets) {
      obj_depmap <- prep_dt_depmap_feat(feature_set = feat,
                                        prefix = prefixes[which(feat == feature_sets)])
      # 2nd level
      for (d_name in drug_name_vec) {
        # 3rd level
        for (norm in normalization_type_vec) {
          # prep data for drug
          if (!is.null(dt_metrics)) {
            dt_response_met <-
              prep_dt_response_metric_sa(dt_metrics = dt_metrics,
                                         d_name = d_name,
                                         normalization_type = norm,
                                         metric = metric)
          }
          if (!is.null(dt_average)) {
            dt_response_dose <- 
              prep_dt_response_dose_sa(dt_average = dt_average, 
                                       d_name = d_name,
                                       normalization_type = norm)
          }
          
          dt_response_sa <- if (!is.null(dt_metrics) && !is.null(dt_average)) {
            merge(dt_response_met, dt_response_dose, by = id_col)
          } else if (!is.null(dt_metrics)) {
            dt_response_met
          } else {
            dt_response_dose
          }
          ls_selected_met <- list(selected_metric = setdiff(names(dt_response_sa), id_col))
          
          # 4th level - prep vis
          ls_vol <- purrr::pmap(ls_selected_met,
                                plot_volcano_assoc_panel,
                                dt_response = dt_response_sa,
                                dt_depmap = obj_depmap[["dt_depmap"]],
                                selected_feat_meta_col = feat)
          names(ls_vol) <- ls_selected_met$selected_metric
          
          ls_plot[[feat]][[d_name]][[norm]] <- ls_vol
        }
        
      }
      rm(obj_depmap)
    }
  }
  if (!is.null(metadata_columns)) {
    # 1st level
    for (meta in metadata_columns) {
      obj_depmap <- prep_dt_depmap_meta(meta_data_path = meta_data_path,
                                        metadata_col = meta)
      # 2nd level
      for (d_name in drug_name_vec) {
        # 3rd level
        for (norm in normalization_type_vec) {
          # prep data for drug
          if (!is.null(dt_metrics)) {
            dt_response_met <-
              prep_dt_response_metric_sa(dt_metrics = dt_metrics,
                                         d_name = d_name,
                                         normalization_type = norm,
                                         metric = metric)
          }
          if (!is.null(dt_average)) {
            dt_response_dose <- 
              prep_dt_response_dose_sa(dt_average = dt_average, 
                                       d_name = d_name,
                                       normalization_type = norm)
          }
          
          dt_response_sa <- if (!is.null(dt_metrics) && !is.null(dt_average)) {
            merge(dt_response_met, dt_response_dose, by = id_col)
          } else if (!is.null(dt_metrics)) {
            dt_response_met
          } else {
            dt_response_dose
          }
          ls_selected_met <- list(selected_metric = setdiff(names(dt_response_sa), id_col))
          
          # 4th level - prep vis
          ls_vol <- purrr::pmap(ls_selected_met,
                                plot_volcano_assoc_panel,
                                dt_response = dt_response_sa,
                                dt_depmap = obj_depmap[["dt_depmap"]],
                                selected_feat_meta_col = meta)
          names(ls_vol) <- ls_selected_met$selected_metric
          
          ls_plot[[meta]][[d_name]][[norm]] <- ls_vol
        }
      }
      rm(obj_depmap)
    }
  }
  return(ls_plot)
}


#' Create a nested list of plots for PRISM data with combo metrics
#' 
#' @param drug1_name_vec character vector with drug names to be plotted (identifiers \code{DrugName})
#' @param drug2_name_vec character vector with co-drug names to be plotted (identifiers \code{DrugName_2})
#' @param dt_metrics \code{data.table} representing data from the \code{Metrics} assay,
#'  outputted by \code{gDRutils::convert_se_assay_to_dt(se, "Metrics")}
#'  and combo \code{SummarizedExperiment}
#' @param dt_scores \code{data.table} representing data from the \code{scores} assay,
#'  outputted by \code{gDRutils::convert_se_assay_to_dt(se, "scores")}
#'  and combo \code{SummarizedExperiment}
#' @param normalization_type_vec character vector with normalization types to be selected
#'                               one of: "GR" ("GRvalue") or "RV" ("RelativeViability") or both
#' @param metric character vector with names of metric for difference:
#'   chosen from: "xc50" ("GR50" or "IC50" - respectively depending on \code{normalization_type}), 
#'  "x_max" ("GR Max" or "E Max") or "x_mean" ("GR Mean" or "RV Mean")
#' @param metric_scores character vector with names of combo metric;
#'   chosen from: "hsa_score"("Bliss Excess GR" or "Bliss Excess RV" - respectively 
#'   depending on \code{normalization_type}), "bliss_score" ("Bliss Score GR" or "Bliss Score RV")
#' @param fit_source string source name for metrics
#' @param meta_data_path string path to metadata file describing all cancer models/cell lines
#'  which are referenced by a dataset contained within the DepMap portal. 
#'  It is usually a file named \code{Model.csv}.
#' @param feature_sets character vector names of the molecular feature set to load from DepMap.
#' @param prefixes character vector prefixes to use for the each feature set in \code{feature_sets};
#'    has to be the same length as \code{feature_sets}
#' @param metadata_columns character vector with the metadata columns to load for DepMap cell lines
#' 
#' @return nested list of plots
#' @keywords prism_plots
#'   
#' @export
create_PRISM_plot_list_combo <- function(drug1_name_vec,
                                         drug2_name_vec,
                                         dt_metrics,
                                         dt_scores,
                                         normalization_type_vec = "RV",
                                         metric =  c("xc50", "x_mean", "x_max"),
                                         metric_scores = c("hsa_score", "bliss_score"),
                                         fit_source = "gDR",
                                         meta_data_path,
                                         feature_sets,
                                         prefixes,
                                         metadata_columns = NULL) {
  
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  drug_name_2 <- gDRutils::get_env_identifiers("drug_name2")
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  id_col <- c("rId", "cId", cellline_name)
  
  checkmate::assert_character(drug1_name_vec, all.missing = FALSE)
  checkmate::assert_character(drug2_name_vec, all.missing = FALSE)
  checkmate::assert_data_table(dt_metrics, min.rows = 1, null.ok = TRUE)
  checkmate::assert_data_table(dt_scores, min.rows = 1, null.ok = TRUE)
  checkmate::assert_subset(normalization_type_vec, choices = c("GR", "RV"))
  checkmate::assert_character(metric, any.missing = FALSE, null.ok = TRUE)
  stopifnot("Provide `metric` for `dt_metrics`." = !is.null(dt_metrics) && !is.null(metric))
  if (!is.null(metric)) checkmate::assert_subset(metric, choices = c("xc50", "x_mean", "x_max"))
  if (!is.null(dt_metrics)) checkmate::assert_subset(metric, choices = names(dt_metrics))
  checkmate::assert_character(metric_scores, any.missing = FALSE, null.ok = TRUE)
  stopifnot("Provide `metric_scores` for `dt_scores`." = !is.null(dt_scores) && !is.null(metric_scores))
  if (!is.null(metric_scores)) checkmate::assert_subset(metric_scores, choices = c("hsa_score", "bliss_score"))
  if (!is.null(dt_scores)) checkmate::assert_subset(metric_scores, choices = names(dt_scores))
  checkmate::assert_string(meta_data_path, pattern = ".*Model.csv$")
  checkmate::assert_file_exists(meta_data_path)
  checkmate::assert_character(feature_sets, null.ok = TRUE)
  checkmate::assert_character(prefixes, null.ok = TRUE)
  stopifnot("`prefixes` has to be the same length as `feature_sets`" = NROW(feature_sets) == NROW(prefixes))
  checkmate::assert_character(metadata_columns, null.ok = TRUE)
  stopifnot("Provide `feature_sets` or `metadata_columns` for DepMam subset." =
              !is.null(feature_sets) || !is.null(metadata_columns))
  
  # prep drugs combinations
  drug_name_grid <- if (!is.null(dt_metrics)) {
    unique(dt_metrics[get(drug_name) %in% drug1_name_vec & get(drug_name_2) %in% drug2_name_vec, 
                      c(drug_name, drug_name_2), with = FALSE])
  } else {
    unique(dt_scores[get(drug_name) %in% drug1_name_vec & get(drug_name_2) %in% drug2_name_vec,
                     c(drug_name, drug_name_2), with = FALSE])
  }
  data.table::setorder(drug_name_grid)
  drug_name_grid[, DrugCombination := paste(get(drug_name), get(drug_name_2), sep = " x ")]
  
  ls_plot <- list()
  
  if (!is.null(feature_sets)) {
    # 1st level
    for (feat in feature_sets) {
      obj_depmap <- prep_dt_depmap_feat(feature_set = feat,
                                        prefix = prefixes[which(feat == feature_sets)])
      # 2nd level
      for (d_combo in drug_name_grid$DrugCombination) {
        d_name <- drug_name_grid[DrugCombination == d_combo, ][[drug_name]]
        d_name2 <- drug_name_grid[DrugCombination == d_combo, ][[drug_name_2]]
        # 3rd level
        for (norm in normalization_type_vec) {
          # prep data for drugs combinations
          if (!is.null(dt_metrics)) {
            dt_response_met_diff <-
              prep_dt_response_metric_diff(dt_metrics = dt_metrics,
                                           d_name = d_name,
                                           d_name2 = d_name2,
                                           normalization_type = norm,
                                           metric = metric)
          }
          if (!is.null(dt_scores)) {
            dt_response_scores <- 
              prep_dt_response_scores(dt_scores = dt_scores,
                                      d_name = d_name,
                                      d_name2 = d_name2,
                                      normalization_type = norm,
                                      metric = metric_scores)
          }
          
          dt_response_combo <- if (!is.null(dt_metrics) && !is.null(dt_scores)) {
            merge(dt_response_scores, dt_response_met_diff, by = id_col)
          } else if (!is.null(dt_metrics)) {
            dt_response_met_diff
          } else {
            dt_response_scores
          }
          ls_selected_met <- list(selected_metric = setdiff(names(dt_response_combo),
                                                            c(id_col, drug_name, drug_name_2)))
          
          # 4th level - prep vis
          ls_vol <- purrr::pmap(ls_selected_met,
                                plot_volcano_assoc_panel,
                                dt_response = dt_response_combo,
                                dt_depmap = obj_depmap[["dt_depmap"]],
                                selected_feat_meta_col = feat)
          names(ls_vol) <- ls_selected_met$selected_metric
          
          ls_plot[[feat]][[d_combo]][[norm]] <- ls_vol
        }
        
      }
      rm(obj_depmap)
    }
  }
  if (!is.null(metadata_columns)) {
    # 1st level
    for (meta in metadata_columns) {
      obj_depmap <- prep_dt_depmap_meta(meta_data_path = meta_data_path,
                                        metadata_col = meta)
      # 2nd level
      for (d_combo in drug_name_grid$DrugCombination) {
        d_name <- drug_name_grid[DrugCombination == d_combo, ][[drug_name]]
        d_name2 <- drug_name_grid[DrugCombination == d_combo, ][[drug_name_2]]
        # 3rd level
        for (norm in normalization_type_vec) {
          # prep data for drugs combinations
          if (!is.null(dt_metrics)) {
            dt_response_met_diff <-
              prep_dt_response_metric_diff(dt_metrics = dt_metrics,
                                           d_name = d_name,
                                           d_name2 = d_name2,
                                           normalization_type = norm,
                                           metric = metric)
          }
          if (!is.null(dt_scores)) {
            dt_response_scores <- 
              prep_dt_response_scores(dt_scores = dt_scores,
                                      d_name = d_name,
                                      d_name2 = d_name2,
                                      normalization_type = norm,
                                      metric = metric_scores)
          }
          
          dt_response_combo <- if (!is.null(dt_metrics) && !is.null(dt_scores)) {
            merge(dt_response_scores, dt_response_met_diff, by = id_col)
          } else if (!is.null(dt_metrics)) {
            dt_response_met_diff
          } else {
            dt_response_scores
          }
          ls_selected_met <- list(selected_metric = setdiff(names(dt_response_combo),
                                                            c(id_col, drug_name, drug_name_2)))
          
          # 4th level - prep vis
          ls_vol <- purrr::pmap(ls_selected_met,
                                plot_volcano_assoc_panel,
                                dt_response = dt_response_combo,
                                dt_depmap = obj_depmap[["dt_depmap"]],
                                selected_feat_meta_col = meta)
          names(ls_vol) <- ls_selected_met$selected_metric
          
          ls_plot[[meta]][[d_combo]][[norm]] <- ls_vol
        }
      }
      rm(obj_depmap)
    }
  }
  return(ls_plot)
}
