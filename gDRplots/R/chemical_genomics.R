#' Analyze Chemical Genomics (CGS) Data and Perform GSEA
#'
#' This function analyzes chemical genomics (CGS) data, filters by compound mechanism of action (MOA),
#' prepares the data for Gene Set Enrichment Analysis (GSEA), and performs GSEA for specified cell lines and metrics.
#'
#' @param metrics_data A data.table containing screening data. Requires columns: `drug_moa`, `DrugName`, `CellLineName`,
#' and columns specified in the `metric` argument.
#' @param metrics A character vector specifying the response metrics to analyze (e.g., "xc50", "x_max").
#' @param cell_line An optional character string specifying a single cell line to analyze. If NULL (default),
#' all cell lines in the data are analyzed.
#' @param normalization_type A character string specifying the normalization type. Default is "RV". 
#' Passed to `gDRplots::prep_dt_response_metric_diff`.
#'
#' @return A list of results, where each element corresponds to a cell line. Each cell line's results contain:
#'   - `fgsea`: A list of GSEA results for each metric.
#'   - `metrics_diff`: The prepped data.table used for the GSEA analysis.
#'   - `moa_list`: A list of DrugNames grouped by drug_moa.
#' @examples
#' metrics_data <- qs::qread(system.file("testdata/cgs_data.qs", package = "gDRplots"))
#' analyze_cgs(metrics_data, metrics = c("xc50"), cell_line = "CellLineName_1")
#' @export
#'
analyze_cgs <- function(metrics_data, metrics, cell_line = NULL, normalization_type = "RV") {
  
  # identifiers
  drug_moa <- gDRutils::get_env_identifiers("drug_moa")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  cl <- gDRutils::get_env_identifiers("cellline_name")
  
  # asserts
  checkmate::assert_data_table(metrics_data)
  checkmate::assert_character(metrics, any.missing = FALSE)
  checkmate::assert_subset(metrics, choices = c("x_mean", "x_AOC_range", "xc50", "ec50", "x_max"), empty.ok = FALSE)
  checkmate::assert_subset(cell_line, choices = unique(metrics_data[[cl]]), empty.ok = TRUE)
  checkmate::assert_choice(normalization_type, choices = c("GR", "RV"))
  
  # filter out unwanted drug moa
  metrics_data <- metrics_data[!eval(drug_moa) %in% c("unknown", "Unknown"), ]
  
  # prepare the data with specified metric differences
  metrics_diff <- prep_dt_response_metric_diff(metrics_data,
                                               metric = metrics,
                                               d_name = NULL,
                                               d_name2 = NULL,
                                               normalization_type = normalization_type,
                                               additional_cols = drug_moa)
  
  original <- grep("diff", names(metrics_diff), value = TRUE)
  new <- gsub(".*gDR_(.*)_cotrt_diff.*", "\\1", original)
  
  data.table::setnames(metrics_diff, original, new)
  
  # create a list of DrugNames grouped by drug_moa, filtering out those with less than 4 unique drugs
  moa_list <- lapply(split(metrics_diff[[drug_name]], metrics_diff[[drug_moa]]), unique)
  moa_list <- moa_list[vapply(moa_list, length, FUN.VALUE = numeric(1)) > 3]
  
  # determine which cell lines to analyze
  if (!is.null(cell_line)) {
    cell_lines <- cell_line
  } else {
    cell_lines <- unique(metrics_diff[[cl]])
  }
  
  # Run fgsea analysis for each specified cell line
  results <- lapply(cell_lines, function(cl) {
    data_subset <- metrics_diff[eval(cl) == cl & eval(drug_moa) %in% names(moa_list), ]
    data_subset$normalization_type <- normalization_type
    list_results <- lapply(metrics, function(metric) {
      metric_values <- data_subset[[metric]]
      names(metric_values) <- data_subset[[drug_name]]
      fgsea_result <- suppressWarnings(fgsea::fgsea(moa_list,
                                                    metric_values,
                                                    500,
                                                    minSize = 4,
                                                    nPermSimple = 1e5))
      
      median_values <- data_subset[, median(get(metric), na.rm = TRUE), by = drug_moa]$V1
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
  return(results)
}


#' Plot Chemical Genomics Screen GSEA Results
#'
#' Generates a ggplot2 visualization of chemical genomics screening data, highlighting GSEA results.
#'
#' @param results A list object returned from `analyze_cgs`.
#' @param cell_line A character string specifying the cell line to prepare data for.
#' @param metric A character string specifying the metric to prepare data for.#'
#'
#' @return A ggplot2 object. The plot is also printed to the console.
#' @examples
#' metrics_data <- qs::qread(system.file("testdata/cgs_data.qs", package = "gDRplots"))
#' results <- analyze_cgs(metrics_data, metrics = c("xc50"), cell_line = "CellLineName_1")
#' plot_cgs_ranking(results, cell_line = "CellLineName_1", metric = "xc50")
#' @export
plot_cgs_ranking <- function(results, cell_line, metric) {
  
  # identifiers
  drug_moa <- gDRutils::get_env_identifiers("drug_moa")
  drug_name <- gDRutils::get_env_identifiers("drug_name")
  cl <- gDRutils::get_env_identifiers("cellline_name")
  norm_type <- gDRutils::get_env_identifiers("normalization_type")
  
  # asserts
  checkmate::assert_list(results)
  checkmate::assert_subset(cell_line, choices = names(results))
  checkmate::assert_subset(metric, choices = names(results[[cell_line]]$fgsea), empty.ok = FALSE)
  
  
  # extract relevant data
  metrics_diff <- results[[cell_line]]$metrics_diff
  fgsea_results <- results[[cell_line]]$fgsea[[metric]]
  moa_groups_drugs <- results[[cell_line]]$moa_list
  
  # prepare data for plotting
  plot_data <- data.table::copy(metrics_diff)
  plot_data$x_pos <- NROW(plot_data) - rank(plot_data[[metric]]) + 1
  stats <- plot_data[[metric]]
  stats <- pmin(2, pmax(-2, stats))
  names(stats) <- plot_data[[drug_name]]
  stats <- stats[!is.na(plot_data[[drug_moa]])]
  
  norm_type <- unique(metrics_diff[[norm_type]])
  
  # filter significant GSEA results
  gsea_sign <- fgsea_results[padj < 0.1 & !pathway %in% c("", "unknown")]
  if (NROW(gsea_sign) == 0) {
    gsea_sign <- fgsea_results[pval < sort(pval)[5]]
  } else if (NROW(gsea_sign) > 15) {
    gsea_sign <- head(gsea_sign[order(padj)], 15)
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
    ggplot2::geom_col(color = "#777777") +
    ggplot2::labs(title = cell_line,
         y = paste0("\u0394 ", metric, " for ", norm_type),
         x = "Ranked drugs",
         caption = "Top 15 results with FDR < 0.1 are shown. If no results meet this threshold,
         the top 4 results by p-value are displayed."
    ) +
    ggplot2::theme_bw() +
    ggplot2::geom_hline(yintercept = 0, color = "#555555") +
    ggplot2::geom_hline(yintercept = mean_effect, color = "black") +
    ggplot2::geom_segment(x = threshold_count, xend = threshold_count, 
                 y = 0, yend = mean_effect + 0.2 * yrange,
                 color = "black") +
    ggplot2::annotate("text", x = threshold_count, y = mean_effect + 0.25 * yrange,
             label = sprintf("Mean effect = %.2f", mean_effect),
             hjust = 0, color = "black") +
    ggplot2::coord_cartesian(xlim = c(-2, NROW(plot_data) + 3),
                             ylim = c(-1.01 * yrange - 0.15 * yrange * NROW(gsea_sign), yrange + 0.01),
                             expand = FALSE, clip = "off") +
    ggplot2::theme(plot.margin = ggplot2::unit(c(1, 16, 1, 1), "lines"))
  
  # define color palettes for the loop (using both Set1 and Set2 if needed)
  n_colors_needed <- NROW(gsea_sign)
  loop_colors <- c(RColorBrewer::brewer.pal(9, "Set1"),
                   RColorBrewer::brewer.pal(8, "Set2"))
  loop_colors <- loop_colors[-6][seq_len(n_colors_needed)]
  
  for (i in seq_len(NROW(gsea_sign))) {
    pathway <- gsea_sign$pathway[i]
    x <- plot_data[get(drug_name) %in% moa_groups_drugs[[pathway]]]$x_pos
    stats_moa <- plot_data[get(drug_name) %in% moa_groups_drugs[[pathway]]][[metric]]
    
    median_moa <- median(stats_moa)
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
        "text",
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
        "label",
        x = count_above_median,
        y = -(gsea_sign$y_pos[i] + 0.5 * sign(gsea_sign$y_pos[i])) * 0.185 * yrange,
        label = sprintf(" %s median = %.2f ", pathway, median_moa),
        hjust = 1 * (gsea_sign$NES[i] < 0),
        color = current_color,
        fill = "white", alpha = 0.65)
  }
  return(plt)
}
