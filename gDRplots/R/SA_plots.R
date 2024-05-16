#' Plot drug response curves for single-agent data
#'
#' @param dt_metrics data.table representation of the data in \code{Metrics} assay
#'    output from \code{gDRutils::convert_se_assay_to_dt(SE, "Metrics")}
#' @param dt_average data.table representation of the data in \code{Averaged} assay
#'    output from \code{gDRutils::convert_se_assay_to_dt(SE, "Averaged")}
#' @param grouping string name of dim of SE to plot; it has to be opppsoite to this used to filter SE
#'    (that is when rownames(SE) == 1 it has to be "cId", otherwise - "rId")
#' @param group_names character vector with names to subset from SE (the same dim as \code{grouping});
#'    if \code{NULL} then all values will be plotted
#' @param normalization_type string with normalization_types to be selected
#'                           one of: "GR" ("GRvalue") or "RV" ("RelativeViability")
#' @param colormap character vector with colors for \code{group_names} - name or hex value
#' @param plot_averaged_flag logical flag whether plot points with average values 
#' @param plot_fit_flag logical flag whether plot points with fitted values 
#'
#' @return plot with dose-response curves
#'
#' @keywords single-agent_plots
#' @examples
#' mae <- gDRutils::get_synthetic_data("small")
#' SE <- mae[[1]]
#' iR <-  rownames(SE)[1]
#' grouping <- "cId"
#' dt_metrics <- gDRutils::convert_se_assay_to_dt(SE[iR], "Metrics")
#' dt_average <- gDRutils::convert_se_assay_to_dt(SE[iR], "Averaged")
#' group_names <- colnames(SE)[2:5]
#' 
#' grob_SA(dt_metrics = dt_metrics, 
#'         dt_average = dt_average, 
#'         grouping = grouping,
#'         group_names = group_names)
#' 
#' @export
grob_SA <- function(dt_metrics, 
                    dt_average, 
                    grouping, 
                    group_names = NULL, 
                    normalization_type = "GR", 
                    colormap = NULL, 
                    plot_averaged_flag = TRUE, 
                    plot_fit_flag = TRUE) {
  
  checkmate::expect_data_table(dt_metrics)
  checkmate::expect_data_table(dt_average)
  checkmate::expect_choice(grouping, choices = c("cId", "rId"))
  checkmate::expect_character(group_names, null.ok = TRUE)
  checkmate::expect_choice(normalization_type, choices = c("GR", "RV"))
  checkmate::expect_character(colormap, null.ok = TRUE)
  checkmate::expect_flag(plot_averaged_flag)
  checkmate::expect_flag(plot_fit_flag)
  
  # check input data
  if (grouping == "cId") {
    stopifnot("grouping` does not fit to `dt_metrics` and `dt_average`" =
                (NROW(unique(dt_metrics[["rId"]])) == 1 && NROW(unique(dt_average[["rId"]])) == 1))
  } else if (grouping == "rId") {
    stopifnot("grouping` does not fit to `dt_metrics` and `dt_average`" =
                (NROW(unique(dt_metrics[["cId"]])) == 1 && NROW(unique(dt_average[["cId"]])) == 1))
  }
  stopifnot("empty plot was selected" = any(plot_averaged_flag, plot_fit_flag))
  
  # filter data for normalization type
  data.table::setkeyv(dt_metrics, "normalization_type")
  dt_met_norm <- dt_metrics[normalization_type]
  data.table::setkeyv(dt_average, "normalization_type")
  dt_avg_norm <- dt_average[normalization_type]
  
  # prep value ranges for plot
  data_range <- c(min(min(dt_avg_norm$x), 0) - 0.05, max(max(dt_avg_norm$x), 1) + 0.05)
  min_conc <- min(dt_avg_norm[dt_avg_norm$Concentration > 0, ]$Concentration)
  max_conc <- max(dt_avg_norm$Concentration)
  conc_range <- 0.5 * c(floor(2 * log10(min_conc) - 0.5), ceiling(2 * log10(max(max_conc)) + 0.3)) 
  # remove conc = 0
  dt_avg_norm$Concentration[dt_avg_norm$Concentration == 0] <- min_conc / 10
  
  # group
  if (is.null(group_names)) group_names <- unique(dt_met_norm[[grouping]])
  
  # prep fitted data
  sel_conc <- 10 ^ (seq(conc_range[1], conc_range[2], 0.05))
  dt_fit <- data.table::data.table()
  
  for (icol in group_names) {
    sel_metrics <- dt_met_norm[dt_met_norm[[grouping]] == icol, ]
    dt_fit <- rbind(dt_fit, 
                    cbind(sel_metrics[, grouping, with = FALSE],
                          data.table::data.table(
                            Concentration = sel_conc,
                            x = gDRutils::predict_efficacy_from_conc(sel_conc, 
                                                                     sel_metrics$x_inf, 
                                                                     sel_metrics$x_0, 
                                                                     sel_metrics$ec50, 
                                                                     sel_metrics$h)
                          )
                    )
    )
  }
  
  # colors
  if (is.null(colormap) || !all(vapply(colormap, gDRcomponents::isValidColor, logical(1)))) {
    colormap <- grDevices::colorRampPalette(c("#dd0000", "#bbdd33", "#0000dd"))(NROW(group_names))
    
  } else if (NROW(colormap) != NROW(group_names)) {
    colormap <- grDevices::colorRampPalette(colormap)(NROW(group_names))
  }
  names(colormap) <- group_names
  color_values <- colormap[names(colormap) %in% unique(dt_fit[[grouping]])]
  
  # prep averaged data
  dt_avg <- dt_avg_norm[get(grouping) %in% group_names, ]
  
  # levels
  dt_avg$grouping <- factor(dt_avg[[grouping]], levels = group_names)
  dt_fit$grouping <- factor(dt_fit[[grouping]], levels = group_names)
  
  plt_title <- 
    sprintf("Drug dose response for %s: %s",
            ifelse(grouping == "cId", "Drug", "Cell Line"),
            ifelse(grouping == "cId", unique(dt_metrics[["rId"]]), unique(dt_metrics[["cId"]]))
    )
  
  # final plot
  plt <- 
    ggplot2::ggplot(mapping = ggplot2::aes(x = log10(Concentration), y = x, color = grouping, group = grouping)) +
    ggplot2::geom_hline(yintercept = c(0, 1), color = "#555555") +
    ggplot2::geom_vline(xintercept = 0, color = "#555555") +
    ggplot2::scale_color_manual(values = color_values,
                                name = ifelse(grouping == "cId", "Cell line", "Drug")) +
    ggplot2::coord_cartesian(xlim = conc_range, ylim = data_range) +
    ggplot2::scale_x_continuous(breaks = -5:2, labels = c("1e-5", "1e-4", 10 ^ (-3:2))) +
    ggplot2::xlab(expression(paste("Concentration [", mu, "M]"))) +
    ggplot2::ylab(paste(normalization_type, "values")) +
    ggplot2::ggtitle(plt_title) +
    ggplot2::theme_bw()
  
  if (plot_averaged_flag) {
    plt <- plt + ggplot2::geom_point(data = dt_avg)
  }
  
  if (plot_fit_flag) {
    plt <- plt + ggplot2::geom_line(data = dt_fit)
  }
  
  return(plt)
}


#' Plot drug response curves for single-agent data for selected call lines and drugs
#' 
#' @inheritParams grob_SA
#' @param SE single-agent \code{SummarizedExperiment} object holding raw and/or processed 
#'    dose-response data in its assays for one cell line
#' @param cellline_name character vector with cell line to be plotted (colnames of SE)
#' @param drug_name character vector with cell line to be plotted (rownames of SE)
#'    
#' @return list of plots with dose-response curves
#' 
#' @examples
#' \dontrun{
#' mae <- gDRutils::get_synthetic_data("small")
#' SE <- mae[[1]]
#' cellline_name <- colnames(SE)[2:5]
#' drug_name <- rownames(SE)[5:7]
#' 
#' plot_SA_byCLs(SE = SE, 
#'               cellline_name = cellline_name, 
#'               drug_name = drug_name, 
#'               normalization_type = "RV", 
#'               colormap = c("#B9D3EE", "#FF6347", "#C2F970"))
#' }
#' 
#' @keywords internal
plot_SA_byCLs <-  function(SE, 
                           cellline_name = NULL, 
                           drug_name = NULL,
                           normalization_type = "GR", 
                           colormap = NULL, 
                           plot_averaged_flag = TRUE, 
                           plot_fit_flag = TRUE) {
  
  checkmate::assert_class(SE, "SummarizedExperiment")
  checkmate::expect_character(cellline_name, null.ok = TRUE)
  checkmate::expect_character(drug_name, null.ok = TRUE)
  checkmate::expect_choice(normalization_type, choices = c("GR", "RV"))
  checkmate::expect_character(colormap, null.ok = TRUE)
  checkmate::expect_flag(plot_averaged_flag)
  checkmate::expect_flag(plot_fit_flag)
  
  if (is.null(drug_name) || all(!drug_name %in% rownames(SE))) {
    drug_name  <- rownames(SE)
  } else if (!all(drug_name  %in% rownames(SE))) {
    drug_name <- drug_name[drug_name  %in% rownames(SE)]
  }  
  
  plt_list <- list()
  for (iR in drug_name) {
    
    if (is.null(cellline_name) || all(!cellline_name %in% colnames(SE))) {
      cellline_name <- colnames(SE)
    } else if (!all(cellline_name %in% colnames(SE))) {
      cellline_name <- cellline_name[cellline_name  %in% colnames(SE)]
    }  
    
    subset_SE <- SE[iR, cellline_name]
    
    plt <- 
      grob_SA(dt_metrics = gDRutils::convert_se_assay_to_dt(subset_SE, "Metrics"), 
              dt_average = gDRutils::convert_se_assay_to_dt(subset_SE, "Averaged"), 
              grouping = "cId",
              group_names = cellline_name,
              normalization_type = normalization_type,
              colormap = colormap,
              plot_averaged_flag = plot_averaged_flag,
              plot_fit_flag = plot_fit_flag) +
      ggplot2::ggtitle(iR) 
    
    plt_list[[iR]] <- plt
    
  }
  
  return(plt_list)
}

#' Plot drug response curves for single-agent data for one selected cell line
#' 
#' @inheritParams grob_SA
#' @param SE single-agent \code{SummarizedExperiment} object holding raw and/or processed 
#'    dose-response data in its assays for one cell line
#'    
#' @return plot with dose-response curves
#' 
#' @examples
#' \dontrun{
#' mae <- gDRutils::get_synthetic_data("small")
#' SE <- mae[[1]]
#' 
#' plot_SA_1CL(SE = SE[,colnames(SE)[1]], colormap = c("cadetblue", "orange", "darkblue"))
#' }
#' 
#' @keywords internal
plot_SA_1CL <- function(SE, 
                        normalization_type = "GR", 
                        colormap = NULL, 
                        plot_averaged_flag = TRUE, 
                        plot_fit_flag = TRUE) {
  
  stopifnot(NCOL(SE) == 1) # plot for 1 cell line
  
  checkmate::assert_class(SE, "SummarizedExperiment")
  checkmate::expect_choice(normalization_type, choices = c("GR", "RV"))
  checkmate::expect_character(colormap, null.ok = TRUE)
  checkmate::expect_flag(plot_averaged_flag)
  checkmate::expect_flag(plot_fit_flag)
  
  grob_SA(
    dt_metrics = gDRutils::convert_se_assay_to_dt(SE, "Metrics"), 
    dt_average = gDRutils::convert_se_assay_to_dt(SE, "Averaged"), 
    grouping = "rId",
    normalization_type = normalization_type,
    colormap = colormap,
    plot_averaged_flag = plot_averaged_flag,
    plot_fit_flag = plot_fit_flag
  ) 
  
}
