#' Create a list of PRISM plots for a selected type of experiment
#' 
#' This internal function was created to avoid code duplication and maintain consistency.
#' For user, it is contained in a single-purpose function, depending on the experiment type.
#' 
#' @param experiment_type string with the type of experiment
#'  one of: "sa" for the single-agent experiment "combo" for the combination experiment
#' @param dt_metrics_sa \code{data.table} representing data from the \code{Metrics} assay,
#'  outputted by \code{gDRutils::convert_se_assay_to_dt(se, "Metrics")}
#'  and single-agent \code{SummarizedExperiment}
#' @param dt_metrics_combo \code{data.table} representing data from the \code{Metrics} assay,
#'  outputted by \code{gDRutils::convert_se_assay_to_dt(se, "Metrics")}
#'  and combo \code{SummarizedExperiment}
#' @param dt_average  \code{data.table} representing data from the \code{Averaged} assay,
#'  outputted by \code{gDRutils::convert_se_assay_to_dt(se, "Averaged")}
#'  and single-agent \code{SummarizedExperiment}
#' @param dt_scores \code{data.table} representing data from the \code{scores} assay,
#'  outputted by \code{gDRutils::convert_se_assay_to_dt(se, "scores")}
#'  and combo \code{SummarizedExperiment}
#' @param metric character vector with names of single-agent metric;
#'   chosen from: "xc50" ("GR50" or "IC50" - respectively depending on \code{normalization_type}),
#'  "x_max" ("GR Max" or "E Max") or "x_mean" ("GR Mean" or "RV Mean")
#' @param metric_scores character vector with names of combo metric;
#'   chosen from: "hsa_score"("Bliss Excess GR" or "Bliss Excess RV" - respectively
#'   depending on \code{normalization_type}), "bliss_score" ("Bliss Score GR" or "Bliss Score RV")
#' @param metric character vector with names of single-agent metric;
#'   chosen from: "xc50" ("GR50" or "IC50" - respectively depending on \code{normalization_type}),
#'  "x_max" ("GR Max" or "E Max") or "x_mean" ("GR Mean" or "RV Mean")
#' @param drug1_name_vec character vector with drug names to be plotted (identifiers \code{DrugName}) 
#' @param drug2_name_vec character vector with co-drug names to be plotted (identifiers \code{DrugName_2})
#' @param normalization_type_vec character vector with normalization types to be selected
#'                               one of: "GR" ("GRvalue") or "RV" ("RelativeViability") or both
#' @param metric character vector with names of metric;
#'   chosen from: "xc50" ("GR50" or "IC50" - respectively depending on \code{normalization_type}),
#'  "x_max" ("GR Max" or "E Max") or "x_mean" ("GR Mean" or "RV Mean")
#' @param metric_scores character vector with names of combo metric;
#'   chosen from: "hsa_score"("Bliss Excess GR" or "Bliss Excess RV" - respectively
#'   depending on \code{normalization_type}), "bliss_score" ("Bliss Score GR" or "Bliss Score RV")
#' @param fit_source string source name for metrics
#' @param meta_data_path string path to metadata file describing all cancer models/cell lines
#'  which are referenced by a dataset contained within the DepMap portal.
#'  It is usually a file named \code{Model.csv}.
#' @param feat_data_path string path to the directory containing the molecular feature set file to load from DepMap.
#' @param feature_sets character vector containing the names of the molecular feature sets to load from DepMap.
#'  These names should also correspond to the file names containing the feature data
#'  (without the extension, which is assumed to be \code{csv})
#' @param metadata_columns character vector with the metadata columns to load for DepMap cell lines
#' @param clear_taxonomy_info logical flag whether to remove taxonomy information for gene names in table
#'  with the molecular feature sets from DepMap.
#' 
#' @return A named list with elements:
#' \itemize{
#'   \item \code{ls_plot} nested list of plots for selected type of experiment
#'   \item \code{ls_assoc_data} nested list of table with association data
#' }
#'
#' @author Janina Smoła \email{janina.smola@@contractors.roche.com}
#'
#' @keywords prism_plots
.create_PRISM_plot_list <- function(experiment_type,
                                    dt_metrics_sa,
                                    dt_metrics_combo = NULL,
                                    dt_average = NULL,
                                    dt_scores = NULL,
                                    drug1_name_vec = NULL,
                                    drug2_name_vec = NULL,
                                    normalization_type_vec = "RV",
                                    metric = c("xc50", "x_mean", "x_max"),
                                    metric_scores = c("hsa_score", "bliss_score"),
                                    fit_source = "gDR",
                                    meta_data_path,
                                    feat_data_path,
                                    feature_sets,
                                    metadata_columns = NULL,
                                    clear_taxonomy_info = TRUE) {
  
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  drug_name_2 <- gDRutils::get_env_identifiers("drug_name2")
  cellline_name <- gDRutils::get_env_identifiers("cellline_name")
  id_col <- c("rId", "cId", cellline_name)
  
  checkmate::assert_choice(experiment_type, choices = c("sa", "combo"))
  checkmate::assert_subset(normalization_type_vec, choices = c("GR", "RV"))
  checkmate::assert_string(fit_source, null.ok = TRUE)
  checkmate::assert_string(meta_data_path)
  checkmate::assert_true(tools::file_ext(meta_data_path) == "csv", .var.name = "File ext must be csv")
  checkmate::assert_file_exists(meta_data_path)
  checkmate::assert_character(metadata_columns, null.ok = TRUE)
  checkmate::assert_string(feat_data_path, null.ok = TRUE)
  checkmate::assert_character(feature_sets, null.ok = TRUE)
  checkmate::assert_flag(clear_taxonomy_info)
  
  if (experiment_type == "sa") {
    checkmate::assert_character(drug1_name_vec, null.ok = TRUE)
    checkmate::assert_data_table(dt_metrics_sa, min.rows = 1, null.ok = TRUE)
    checkmate::assert_data_table(dt_average, min.rows = 1, null.ok = TRUE)
    stopifnot("Provide response data - at least one of `dt_metrics` or `dt_average`." =
                !is.null(dt_metrics_sa) || !is.null(dt_average))
    if (!is.null(metric)) checkmate::assert_subset(metric, 
                                                   choices = c("xc50", "x_mean", "x_max"))
    if (!is.null(dt_metrics_sa)) checkmate::assert_subset(metric, 
                                                          choices = names(dt_metrics_sa), 
                                                          empty.ok = FALSE)
  } else { # combo
    checkmate::assert_character(drug1_name_vec, all.missing = FALSE)
    checkmate::assert_character(drug2_name_vec, all.missing = FALSE)
    checkmate::assert_data_table(dt_metrics_combo, min.rows = 1, null.ok = TRUE)
    checkmate::assert_data_table(dt_scores, min.rows = 1, null.ok = TRUE)
    stopifnot("Provide response data - at least one of `dt_metrics` or `dt_scores`." =
                !is.null(dt_metrics_combo) || !is.null(dt_scores))
    if (!is.null(metric)) checkmate::assert_subset(metric, 
                                                   choices = c("xc50", "x_mean", "x_max"))
    if (!is.null(dt_metrics_combo)) checkmate::assert_subset(metric, 
                                                             choices = names(dt_metrics_combo), 
                                                             empty.ok = FALSE)
    if (!is.null(metric_scores)) checkmate::assert_subset(metric_scores, 
                                                          choices = c("hsa_score", "bliss_score"))
    if (!is.null(dt_scores)) checkmate::assert_subset(metric_scores, 
                                                      choices = names(dt_scores), 
                                                      empty.ok = FALSE)
  }
  
  # depmap specyfic
  if (!is.null(metadata_columns)) {
    meta_data <- data.table::fread(meta_data_path, nrow = 1)
    checkmate::assert_subset(c("ModelID", "CCLEName"), names(meta_data))
    metadata_columns <- intersect(metadata_columns, names(meta_data))
    if (NROW(metadata_columns) == 0) metadata_columns <- NULL
    rm(meta_data)
  }
  if (!is.null(feat_data_path) && !is.null(feature_sets)) {
    checkmate::assert_directory_exists(feat_data_path)
    feature_sets <- feature_sets[vapply(feature_sets, function(f) {
      file.exists(file.path(feat_data_path, paste0(f, ".csv")))
    }, logical(1))]
    if (NROW(feature_sets) == 0) feature_sets <- NULL
  }
  if (xor(is.null(feature_sets), is.null(feat_data_path))) {
    stopifnot("Provide consistent values for `feature_sets` and `feat_data_path` for DepMap subset." =
                !is.null(metadata_columns))
    if (!is.null(metadata_columns)) {
      feature_sets <- feat_data_path <- NULL
    }
  }
  stopifnot("Provide `feature_sets` or `metadata_columns` for DepMap subset." =
              !is.null(feature_sets) || !is.null(metadata_columns))
  
  # select data for normalization type
  filter_expr <- substitute(normalization_type %in% norm_type & fit_source == fit_src,
                            list(norm_type = normalization_type_vec, fit_src = fit_source))
  
  if (experiment_type == "sa") {
    if (!is.null(dt_metrics_sa)) dt_metrics_sa <- dt_metrics_sa[eval(filter_expr)]
    if (!is.null(dt_average)) dt_average <- dt_average[eval(filter_expr)]
  } else { # combo
    if (!is.null(dt_metrics_combo)) dt_metrics_combo <- dt_metrics_combo[eval(filter_expr)]
    if (!is.null(dt_scores)) dt_scores <- dt_scores[eval(filter_expr)]
  }
  
  # adjust drug list
  iter_id <- NULL
  if (experiment_type == "sa") {
    dt_drug_source <- if (!is.null(dt_metrics_sa)) dt_metrics_sa else dt_average
    available_drugs <- unique(dt_drug_source[[drug_name]])
    
    if (is.null(drug1_name_vec) || all(!drug1_name_vec %in% available_drugs)) {
      drug1_name_vec <- available_drugs
    } else {
      drug1_name_vec <- drug1_name_vec[drug1_name_vec %in% available_drugs]
    }
    drug_name_grid <- data.table::data.table(iter_id = drug1_name_vec, 
                                             drug_name = drug1_name_vec)
    data.table::setnames(drug_name_grid, "drug_name", drug_name)
  } else { # combo
    dt_drug_source <- if (!is.null(dt_metrics_combo)) dt_metrics_combo else dt_scores
    drug_name_grid <- 
      unique(dt_drug_source[get(drug_name) %in% drug1_name_vec & get(drug_name_2) %in% drug2_name_vec,
                            c(drug_name, drug_name_2), with = FALSE])
    data.table::setorder(drug_name_grid)
    drug_name_grid[, iter_id := paste(get(drug_name), get(drug_name_2), sep = " x ")]
    drug_name_grid
  }
  
  # fast end when there are no available drugs
  if (NROW(drug_name_grid) == 0) {
    if (experiment_type == "sa") {
      message("There was no data for selected drugs.")
    } else {
      message("There was no data for selected drugs combination.")
    }
    return(list())
  }
  
  # features and meta for looping
  depmap_items <- data.table::data.table()
  if (!is.null(feature_sets)) {
    depmap_items <- rbind(depmap_items,
                          data.table::data.table(feat_meta_name = feature_sets, 
                                                 feat_meta_type = "feature"))
  }
  if (!is.null(metadata_columns)) {
    depmap_items <- rbind(depmap_items,
                          data.table::data.table(feat_meta_name = metadata_columns, 
                                                 feat_meta_type = "meta"))
  }
  
  # final object
  ls_plot <- list()
  ls_assoc_data <- list()
  
  # 1st level
  for (i_feat_meta in seq_len(NROW(depmap_items))) {
    obj_depmap <- if (depmap_items$feat_meta_type[[i_feat_meta]] == "feature") {
      prep_dt_depmap_feat(feat_data_path = feat_data_path, 
                          meta_data_path = meta_data_path, 
                          feature_set = depmap_items$feat_meta_name[[i_feat_meta]])
    } else {
      prep_dt_depmap_meta(meta_data_path = meta_data_path,
                          metadata_col = depmap_items$feat_meta_name[[i_feat_meta]])
    }
    
    dt_depmap <- obj_depmap[["dt_depmap"]]
    if (clear_taxonomy_info && depmap_items$feat_meta_type[[i_feat_meta]] == "feature") {
      names(dt_depmap) <- gsub(" \\((.*))", "", names(dt_depmap))
    }
    
    # 2nd level
    for (i_drug in seq_len(NROW(drug_name_grid))) {
      # 3rd level
      for (norm in normalization_type_vec) {
        if (experiment_type == "sa") {
          d_name <- drug_name_grid[[drug_name]][i_drug]
          
          dt_met <- 
            if (!is.null(dt_metrics_sa)) prep_dt_response_metric_sa(dt_metrics = dt_metrics_sa, 
                                                                    d_name = d_name, 
                                                                    normalization_type = norm, 
                                                                    metric = metric,
                                                                    fit_source = fit_source)
          dt_dose <- 
            if (!is.null(dt_average)) prep_dt_response_dose_sa(dt_average = dt_average, 
                                                               d_name = d_name, 
                                                               normalization_type = norm,
                                                               fit_source = fit_source)
          
          dt_response <- if (!is.null(dt_met) && !is.null(dt_dose)) {
            merge(dt_met, dt_dose, by = id_col)
          } else if (!is.null(dt_met)) {
            dt_met 
          } else {
            dt_dose
          }
          
          id_cols_to_remove <- id_col
          
        } else { # combo
          d_name <- drug_name_grid[[drug_name]][i_drug]
          d_name2 <- drug_name_grid[[drug_name_2]][i_drug]
          
          dt_met_diff <- 
            if (!is.null(dt_metrics_combo)) prep_dt_response_metric_diff(dt_metrics = dt_metrics_combo, 
                                                                         d_name = d_name, 
                                                                         d_name2 = d_name2, 
                                                                         normalization_type = norm, 
                                                                         metric = metric,
                                                                         fit_source = fit_source)
          dt_scores_res <- 
            if (!is.null(dt_scores)) prep_dt_response_scores(dt_scores = dt_scores, 
                                                             d_name = d_name, 
                                                             d_name2 = d_name2, 
                                                             normalization_type = norm, 
                                                             metric = metric_scores,
                                                             fit_source = fit_source)
          
          dt_response <- if (!is.null(dt_scores_res) && !is.null(dt_met_diff)) {
            merge(dt_scores_res, dt_met_diff, by = id_col)
          } else if (!is.null(dt_met_diff)) {
            dt_met_diff 
          } else {
            dt_scores_res
          }
          
          id_cols_to_remove <- c(id_col, drug_name, drug_name_2)
        }
        # 4th level - prep vis for each metric
        selected_metrics <- list(selected_metric = setdiff(names(dt_response), id_cols_to_remove))
        
        obj_vol <- purrr::pmap(selected_metrics,
                               plot_volcano_assoc_panel,
                               dt_response = dt_response,
                               dt_depmap = dt_depmap,
                               selected_feat_meta_col = obj_depmap[["selected_feat_meta_col"]])
        names(obj_vol) <- selected_metrics$selected_metric
        ls_vol <- purrr::map(obj_vol, "panel")
        ls_tab <- purrr::map(obj_vol, "assoc_data")

        id_feat_meta <- depmap_items$feat_meta_name[[i_feat_meta]]
        id_drug <- drug_name_grid[["iter_id"]][i_drug]
        
        ls_plot[[id_feat_meta]][[id_drug]][[norm]] <- ls_vol
        ls_assoc_data[[id_feat_meta]][[id_drug]][[norm]] <- ls_tab
      }
    }
    rm(obj_depmap)
  }
  return(list(ls_plot = ls_plot,
              ls_assoc_data = ls_assoc_data))
}

#' Create a nested list of plots for PRISM data with single-agent metrics
#'
#' @inheritParams .create_PRISM_plot_list
#' @param drug_name_vec character vector with drug names to be plotted (identifiers \code{DrugName}) 
#' @param dt_metrics \code{data.table} representing data from the \code{Metrics} assay,
#'  outputted by \code{gDRutils::convert_se_assay_to_dt(se, "Metrics")}
#'  and single-agent \code{SummarizedExperiment}
#'
#' @inherit .create_PRISM_plot_list return
#'
#' @author Janina Smoła \email{janina.smola@@contractors.roche.com}
#'
#' @keywords prism_plots
#'
#' @export
create_PRISM_plot_list_sa <- function(drug_name_vec,
                                      dt_metrics,
                                      dt_average = NULL,
                                      normalization_type_vec = "RV",
                                      metric = c("xc50", "x_mean", "x_max"),
                                      fit_source = "gDR",
                                      meta_data_path,
                                      feat_data_path,
                                      feature_sets,
                                      metadata_columns = NULL,
                                      clear_taxonomy_info = TRUE) {
  .create_PRISM_plot_list(
    experiment_type = "sa",
    drug1_name_vec = drug_name_vec,
    dt_metrics_sa = dt_metrics,
    dt_metrics_combo = NULL,
    dt_average = dt_average,
    normalization_type_vec = normalization_type_vec,
    metric = metric,
    fit_source = fit_source,
    meta_data_path = meta_data_path,
    feat_data_path = feat_data_path,
    feature_sets = feature_sets,
    metadata_columns = metadata_columns,
    clear_taxonomy_info = clear_taxonomy_info
  )
}

#' Create a nested list of plots for PRISM data with combo metrics
#'
#' @inheritParams .create_PRISM_plot_list
#' @param drug1_name_vec character vector with drug names to be plotted (identifiers \code{DrugName})
#' @param dt_metrics \code{data.table} representing data from the \code{Metrics} assay,
#'  outputted by \code{gDRutils::convert_se_assay_to_dt(se, "Metrics")}
#'  and combo \code{SummarizedExperiment}
#'
#' @inherit .create_PRISM_plot_list return
#'
#' @author Janina Smoła \email{janina.smola@@contractors.roche.com}
#'
#' @keywords prism_plots
#'
#' @export
create_PRISM_plot_list_combo <- function(drug1_name_vec,
                                         drug2_name_vec,
                                         dt_metrics,
                                         dt_scores = NULL,
                                         normalization_type_vec = "RV",
                                         metric = c("xc50", "x_mean", "x_max"),
                                         metric_scores = c("hsa_score", "bliss_score"),
                                         fit_source = "gDR",
                                         meta_data_path,
                                         feat_data_path,
                                         feature_sets,
                                         metadata_columns = NULL,
                                         clear_taxonomy_info = TRUE) {
  .create_PRISM_plot_list(
    experiment_type = "combo",
    drug1_name_vec = drug1_name_vec,
    drug2_name_vec = drug2_name_vec,
    dt_metrics_sa = NULL,
    dt_metrics_combo = dt_metrics,
    dt_scores = dt_scores,
    normalization_type_vec = normalization_type_vec,
    metric = metric,
    metric_scores = metric_scores,
    fit_source = fit_source,
    meta_data_path = meta_data_path,
    feat_data_path = feat_data_path,
    feature_sets = feature_sets,
    metadata_columns = metadata_columns,
    clear_taxonomy_info = clear_taxonomy_info
  )
}
