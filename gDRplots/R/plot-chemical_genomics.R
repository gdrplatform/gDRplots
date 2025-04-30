#' Analyze Chemical Genomics (CGS) Data and Perform GSEA
#'
#' This function analyzes chemical genomics (CGS) data, filters by compound mechanism of action (MOA),
#' prepares the data for Gene Set Enrichment Analysis (GSEA), and performs GSEA for specified cell lines and metrics.
#'
#' @param dt_metrics A data.table containing screening data. Requires columns: `drug_moa`, `DrugName`, `CellLineName`,
#' and columns specified in the `metrics` argument.
#' @param metrics A character vector specifying the response metrics to analyze:
#' "x_mean", "x_AOC_range", "xc50", "ec50", "x_max".
#' @param cl_name An optional string specifying a single cell line to analyze. If NULL (default),
#' all cell lines in the data are analyzed. Should be NULL if `resistant_cl` and `sensitive_cl` are provided.
#' @param resistant_cl An optional string representing the resistant cell line name. 
#' Should be specified alongside `sensitive_cl`.
#' @param sensitive_cl An optional string representing the sensitive cell line name. 
#' Should be specified alongside `resistant_cl`.
#' @param normalization_type A string with normalization types to be selected, one of:
#' "GR" ("GRvalue") or "RV" ("RelativeViability").
#' Passed to `gDRplots::prep_dt_response_metric_diff`.
#'
#' @return A list of results, where each element corresponds to a cell line or cell line difference.
#' Each result contains:
#' \itemize{
#'   \item \code{fgsea}: GSEA results for the specified metrics,
#'   \item \code{metrics_diff}: The prepped data.table used for the GSEA analysis,
#'   \item \code{moa_list}: A list of DrugNames grouped by drug_moa.
#' }
#' @keywords cgs_plots
#' 
#' @examples
#' dt_metrics <- qs::qread(system.file("testdata/cgs_data.qs", package = "gDRplots"))
#' analyze_cgs(dt_metrics, metrics = "xc50", cl_name = "CellLineName_1")
#' @export
#'
analyze_cgs <- function(dt_metrics,
                        metrics,
                        cl_name = NULL,
                        resistant_cl = NULL,
                        sensitive_cl = NULL,
                        normalization_type = "RV") {
  
  # identifiers
  drug_moa <- gDRutils::get_env_identifiers("drug_moa")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  cellline <- gDRutils::get_env_identifiers("cellline_name")
  
  # asserts
  checkmate::assert_data_table(dt_metrics)
  checkmate::assert_character(metrics, any.missing = FALSE)
  checkmate::assert_subset(metrics, choices = c("x_mean", "x_AOC_range", "xc50", "ec50", "x_max"), empty.ok = FALSE)
  checkmate::assert_true(all(metrics %in% names(dt_metrics)))
  checkmate::assert_choice(normalization_type, choices = c("GR", "RV"))
  
  if (!is.null(resistant_cl) || !is.null(sensitive_cl)) {
    checkmate::assert_string(resistant_cl, null.ok = FALSE)
    checkmate::assert_string(sensitive_cl, null.ok = FALSE)
    checkmate::assert_subset(resistant_cl, choices = unique(dt_metrics[[cellline]]), empty.ok = FALSE)
    checkmate::assert_subset(sensitive_cl, choices = unique(dt_metrics[[cellline]]), empty.ok = FALSE)
  } else {
    checkmate::assert_string(cl_name, null.ok = TRUE)
    checkmate::assert_subset(cl_name, choices = unique(dt_metrics[[cellline]]), empty.ok = TRUE)
  }
  
  # filter out unwanted drug moa
  dt_metrics <- dt_metrics[!eval(drug_moa) %in% c("unknown", "Unknown"), ]
  
  # take care of Inf and NaN values in IC50 metrics
  if (any(is.infinite(dt_metrics$xc50))) {
    message("In `dt_metrics` some xc50 values are infinite.")
    dt_metrics <- dt_metrics[!is.infinite(get("xc50")), ]
  }
  
  # prepare the data with specified metric differences
  metrics_diff <- prep_dt_response_metric_diff(dt_metrics,
                                               metric = metrics,
                                               d_name = NULL,
                                               d_name2 = NULL,
                                               resistant_cl = resistant_cl,
                                               sensitive_cl = sensitive_cl,
                                               normalization_type = normalization_type,
                                               additional_cols = drug_moa)
  
  to_remove <- names(metrics_diff)[grepl("_fittings$", names(metrics_diff))]
  to_remove <- to_remove[!grepl("cotrt_diff", to_remove)]
  metrics_diff <- metrics_diff[, -c(to_remove), with = FALSE]
  
  original <- grep("cotrt_diff|cellline_diff", names(metrics_diff), value = TRUE)
  new <- gsub(".*gDR_(log10_)?(.*)_cotrt_diff.*", "\\2", original)
  new <- gsub("_cellline_diff", "", new)
  
  data.table::setnames(metrics_diff, original, new)
  
  # create a list of DrugNames grouped by drug_moa, filtering out those with less than 4 unique drugs
  moa_list <- lapply(split(metrics_diff[[drug_name]], metrics_diff[[drug_moa]]), unique)
  moa_list <- moa_list[vapply(moa_list, length, FUN.VALUE = numeric(1)) > 3]
  
  metrics_diff <- metrics_diff[eval(drug_moa) %chin% names(moa_list)]
  
  # Determine FGSEA analysis path
  if (!is.null(resistant_cl) && !is.null(sensitive_cl)) {
    # Perform FGSEA on the difference between the specified cell lines
    list_results <- lapply(metrics, function(metric) {
      metric_values <- metrics_diff[[metric]]
      names(metric_values) <- metrics_diff[[drug_name]]
      fgsea_result <- purrr::quietly(fgsea::fgsea)(pathways = moa_list,
                                                   stats = metric_values,
                                                   maxSize = 500,
                                                   minSize = 4,
                                                   nPermSimple = 1e5)$result
      
      median_values <- metrics_diff[, stats::median(get(metric), na.rm = TRUE), by = drug_moa]$V1
      names(median_values) <- metrics_diff[, unique(drug_moa)]
      
      fgsea_result$median <- median_values[fgsea_result$pathway]
      data.table::setorder(fgsea_result, -NES)
      fgsea_result
    })
    names(list_results) <- metrics
    results <- list(fgsea = list_results,
                    metrics_diff = metrics_diff,
                    moa_list = moa_list)
    results <- list(results)
    names(results) <- unique(paste0(metrics_diff$CellLineName_c1, "_", metrics_diff$CellLineName_c2))
  } else {
    # Analyze all or a single specified cell line
    cell_lines <- if (!is.null(cl_name)) {
      cl_name
    } else {
      unique(metrics_diff[[cellline]])
    }
    
    results <- lapply(cell_lines, function(cl) {
      data_subset <- metrics_diff[get(cellline) == cl]
      data_subset$normalization_type <- normalization_type
      list_results <- lapply(metrics, function(metric) {
        metric_values <- data_subset[[metric]]
        names(metric_values) <- data_subset[[drug_name]]
        fgsea_result <- purrr::quietly(fgsea::fgsea)(pathways = moa_list,
                                                     stats = metric_values,
                                                     maxSize = 500,
                                                     minSize = 4,
                                                     nPermSimple = 1e5)$result
        
        median_values <- data_subset[, stats::median(get(metric), na.rm = TRUE), by = drug_moa]$V1
        names(median_values) <- data_subset[, unique(drug_moa)]
        
        fgsea_result$median <- median_values[fgsea_result$pathway]
        data.table::setorder(fgsea_result, -NES)
        return(fgsea_result)
      })
      names(list_results) <- metrics
      list(fgsea = list_results,
           metrics_diff = data_subset,
           moa_list = moa_list)
    })
    names(results) <- cell_lines
  }
  results
}


#' Plot Chemical Genomics Screen GSEA Results
#'
#' Generates a ggplot2 visualization of chemical genomics screening data, highlighting GSEA results.
#'
#' @param results A list object returned from `analyze_cgs`.
#' @param cl_name A string specifying the cell line included in the \code{results} to prepare a visualization.
#' @param metric A string specifying the metric included in the \code{results} to prepare a visualization.
#' @param padj_threshold A numeric value specifying the threshold for filtering significant GSEA results
#' based on adjusted p-value.
#' @param top_results_no_sig A numeric value specifying the number of top results to plot
#' if there are no significant values.
#' @param max_results_with_sig A numeric value specifying the maximum number of results
#' to plot when there are more than this number of significant values.
#'
#' @return A ggplot2 object with cgs results
#' @keywords cgs_plots
#' @examples
#' dt_metrics <- qs::qread(system.file("testdata/cgs_data.qs", package = "gDRplots"))
#' results <- analyze_cgs(dt_metrics, metrics = c("xc50"), cl_name = "CellLineName_1")
#' plot_cgs_ranking(results,
#'   cl_name = "CellLineName_1",
#'   metric = "xc50",
#'   padj_threshold = 0.1,
#'   top_results_no_sig = 5,
#'   max_results_with_sig = 15)
#' @export
plot_cgs_ranking <- function(results,
                             cl_name,
                             metric,
                             padj_threshold = 0.1,
                             top_results_no_sig = 5,
                             max_results_with_sig = 15) {
  
  # identifiers
  drug_moa <- gDRutils::get_env_identifiers("drug_moa")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  norm_type <- gDRutils::get_env_identifiers("normalization_type")
  
  # asserts
  checkmate::assert_list(results)
  checkmate::assert_string(cl_name, null.ok = FALSE)
  checkmate::assert_subset(cl_name, choices = names(results))
  checkmate::assert_string(metric, null.ok = FALSE)
  checkmate::assert_subset(metric, choices = names(results[[cl_name]]$fgsea), empty.ok = FALSE)
  
  
  # extract relevant data
  metrics_diff <- results[[cl_name]]$metrics_diff
  fgsea_results <- results[[cl_name]]$fgsea[[metric]]
  moa_groups_drugs <- results[[cl_name]]$moa_list
  
  # prepare data for plotting
  plot_data <- data.table::copy(metrics_diff)
  plot_data$x_pos <- NROW(plot_data) - rank(plot_data[[metric]]) + 1
  stats <- plot_data[[metric]]
  stats <- pmin(2, pmax(-2, stats))
  names(stats) <- plot_data[[drug_name]]
  stats <- stats[!is.na(plot_data[[drug_moa]])]
  
  norm_type <- unique(metrics_diff[[norm_type]])
  
  # filter significant GSEA results
  gsea_sign <- fgsea_results[padj < padj_threshold & !pathway %in% c("", "unknown")]
  if (NROW(gsea_sign) == 0) {
    # if there are no significant values, plot only top_results_no_sig top results
    gsea_sign <- fgsea_results[pval < sort(pval)[top_results_no_sig]]
  } else if (NROW(gsea_sign) > max_results_with_sig) {
    # if there are more than max_results_with_sig significant values, plot only max_results_with_sig top results
    gsea_sign <- utils::head(gsea_sign[order(padj)], max_results_with_sig)
  }
  gsea_sign[, y_pos := seq_len(.N)]
  gsea_sign[NES < 0, y_pos := -(seq_len(.N))]
  
  # pre-calculate statistics
  stats <- plot_data[[metric]]
  threshold_count <- sum(stats > mean(stats))
  yrange <- diff(range(stats))
  mean_effect <- mean(stats)
  
  # create the ggplot object
  plt <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x_pos, y = !!rlang::sym(metric))) +
    ggplot2::geom_col(color = gDRutils::get_settings_from_json("EDGE_COLOR",
                                                               system.file(package = "gDRplots", "settings.json"))) +
    ggplot2::labs(title = cl_name,
                  y = bquote(~ Delta ~ .(metric) ~ "for" ~ .(norm_type)),
                  x = "Ranked drugs",
                  caption = sprintf("Top results with FDR < %.2f are shown. If no results meet this threshold,
                                    the top %d results by p-value are displayed.",
                                    padj_threshold, top_results_no_sig)
                  ) +
    ggplot2::theme_bw() +
    ggplot2::geom_hline(yintercept = 0,
                        color = gDRutils::get_settings_from_json("HLINE_COLOR",
                                                                 system.file(package = "gDRplots", "settings.json"))) +
    ggplot2::geom_hline(yintercept = mean_effect, color = "black") +
    ggplot2::geom_segment(x = threshold_count, xend = threshold_count, 
                          y = 0, yend = mean_effect + 0.2 * yrange,
                          color = "black") +
    ggplot2::annotate(geom = "text", x = threshold_count, y = mean_effect + 0.25 * yrange,
                      label = sprintf("Mean effect = %.2f", mean_effect),
                      hjust = 0, color = "black") +
    ggplot2::coord_cartesian(xlim = c(-2, NROW(plot_data) + 3),
                             ylim = c(-1.01 * yrange - 0.15 * yrange * NROW(gsea_sign), yrange + 0.01),
                             expand = FALSE, clip = "off") +
    ggplot2::theme(plot.margin = ggplot2::unit(c(1, 16, 1, 1), "lines"))
  
  # define color palettes for the loop
  loop_colors <- gDRplots::get_qual_colors(NROW(gsea_sign))
  
  for (i in seq_len(NROW(gsea_sign))) {
    pathway <- gsea_sign$pathway[i]
    x <- plot_data[get(drug_name) %in% moa_groups_drugs[[pathway]]]$x_pos
    stats_moa <- plot_data[get(drug_name) %in% moa_groups_drugs[[pathway]]][[metric]]
    
    median_moa <- stats::median(stats_moa)
    count_above_median <- sum(stats > median_moa)
    
    current_color <- loop_colors[i]
    
    plt <- plt +
      ggplot2::geom_segment(
        data = data.frame(x = x),
        ggplot2::aes(x = x, xend = x),
        y = -yrange - (0.15 * yrange * (i - 1)),
        yend = -yrange - (0.15 * yrange * i),
        linewidth = 0.8, inherit.aes = FALSE, color = current_color) +
      ggplot2::annotate(
        geom = "text",
        x = max(plot_data$x_pos) + 10,
        y = -yrange - (0.15 * yrange * (i - 0.5)),
        label = sprintf("%s, NES=%.2f, FDR=%.2g", gsub("_", " ", pathway), gsea_sign$NES[i], gsea_sign$padj[i]),
        hjust = 0, color = current_color) +
      ggplot2::geom_segment(
        x = count_above_median - 0.5,
        xend = count_above_median - 0.5,
        y = median_moa, yend = -(gsea_sign$y_pos[i] + 0.5 * sign(gsea_sign$y_pos[i])) * 0.15 * yrange,
        color = current_color) +
      ggplot2::annotate(
        geom = "label",
        x = count_above_median,
        y = -(gsea_sign$y_pos[i] + 0.5 * sign(gsea_sign$y_pos[i])) * 0.185 * yrange,
        label = sprintf(" %s median = %.2f ", pathway, median_moa),
        hjust = 1 * (gsea_sign$NES[i] < 0),
        color = current_color,
        fill = "white", alpha = 0.65)
  }
  return(plt)
}
