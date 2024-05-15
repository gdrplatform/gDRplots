# selecting data subset ----
cell_line <- "cellline_AA"
drug1_name <- "drug_001"
drug2_name <- "drug_026"
norm_type <- "RV"

mae <- gDRutils::get_synthetic_data("combo_matrix")
SE <- mae[["combination"]]

selected_col <- 
  SummarizedExperiment::colData(SE)[SummarizedExperiment::colData(SE)$CellLineName == cell_line, ]
selected_row <- 
  SummarizedExperiment::rowData(SE)[SummarizedExperiment::rowData(SE)$DrugName == drug1_name & 
                                      SummarizedExperiment::rowData(SE)$DrugName_2 == drug2_name, ]
se1 <- SE[rownames(selected_row), rownames(selected_col)]

# prep data source ----
dt_mx <- data.table::as.data.table(
  BumpyMatrix::unsplitAsDataFrame(SummarizedExperiment::assay(se1, "excess")))
dt_mx <- dt_mx[normalization_type == norm_type, ]

dt_iso <- gDRutils::convert_se_assay_to_dt(se1, "isobolograms")
dt_iso <- dt_iso[normalization_type == norm_type, ][, .(iso_level, pos_x, pos_y)]
all_iso <- unique(dt_iso$iso_level)
iso_colors <- gDRutils::get_iso_colors()[all_iso]

# nolint start
# # note: this is now wrapped in gDRutils::get_iso_colors 
# all_iso <- unique(dt_iso$iso_level)
# iso_cutoff <- as.numeric(all_iso)
# 
# if (normalization_type == "GR") {
#   iso_colors = sapply(iso_cutoff, function(x)
#     sprintf("#%s",paste(as.hexmode(c(70, round((0.85 - x * 0.7) * 170), round((1.1 - x * 0.7)*200))),
#                         collapse = "")))
# } else {
#   iso_colors = sapply(iso_cutoff, function(x)
#     sprintf("#%s", paste(as.hexmode(c(70, round((1 - x * 0.85) * 170), round((1.1 - x * 0.85) * 240))),
#                         collapse = "")))
# }
# names(iso_colors) <- iso_cutoff
# nolint end

mx_names <- names(gDRutils::get_combo_excess_field_names())
ls_mx <- lapply(mx_names, function(mx_name) {
  dt_ <- dt_mx[, c("Concentration", "Concentration_2", mx_name), with = FALSE]
  dt_[, mx_name] <- pmin(1.1, dt_[[mx_name]])
  dt_
})
names(ls_mx) <- mx_names

# heatmaps with average value of combo metrics - conc drug vs conc codrug  ----
plts <- lapply(mx_names, function(mx_name) {
  df_ <- ls_mx[[mx_name]]
  
  x_conc <- sort(unique(df_$Concentration_2[df_$Concentration_2 > 0]))
  x_gap <- max(diff(log10(x_conc)))
  df_$pos_x <- pmax(log10(df_$Concentration_2), 
                    log10(min(df_$Concentration_2[df_$Concentration_2 > 0])) - x_gap) # - 0.1)
  df_$marks_x <- sprintf("%.2g", df_$Concentration_2)
  
  y_conc  <- sort(unique(df_$Concentration[df_$Concentration > 0]))
  y_gap <- max(diff(log10(y_conc)))
  df_$pos_y <- pmax(log10(df_$Concentration), 
                    log10(min(df_$Concentration[df_$Concentration > 0])) - y_gap) # - 0.1)
  df_$marks_y <- sprintf("%.2g", df_$Concentration)
  
  plt_title <- sprintf("%s (%s) : %s for %s, T=%sh",
                       cell_line,
                       selected_col$clid, 
                       gDRutils::get_combo_excess_field_names()[[mx_name]],
                       norm_type,
                       selected_row$Duration)
  # base plot
  plt <- 
    ggplot2::ggplot(df_, ggplot2::aes(x = pos_x, y = pos_y)) +
    ggplot2::geom_tile(ggplot2::aes(fill = get(mx_name)), height = y_gap, width = x_gap) +
    ggplot2::labs(x = bquote(.(drug2_name) ~ "[" ~ mu * M ~ "]"),
                  y = bquote(.(drug1_name) ~ "[" ~ mu * M ~ "]"),
                  fill = gDRutils::get_combo_excess_field_names()[[mx_name]],
                  title = plt_title) +
    ggplot2::theme_bw() + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 9, angle = 45, vjust = 1, hjust = 1),
                   axis.text.y = ggplot2::element_text(size = 9),
                   plot.title = ggplot2::element_text(size = 11)) +
    ggplot2::scale_x_continuous(breaks = df_$pos_x, labels = df_$marks_x,
                                expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = df_$pos_y, labels = df_$marks_y,
                                expand = c(0, 0)) +
    ggplot2::scale_shape_discrete(name = paste0(ifelse(norm_type == "GR", "GR", "IC"), "50"))
  
  # add color scale
  if (!(mx_name %in% c("hsa_excess", "bliss_excess"))) { # heatmaps with readout values
    if (norm_type == "GR") {
      plt <- plt +
        ggplot2::scale_fill_gradientn(
          colors = c("black", "#b06000", "#c07700", "white"),
          values = c(0, 0.59 / 1.7, 0.61 / 1.7, 1), 
          limits = c(-0.6, 1.1), 
          name = "GR val",
          oob = scales::squish)
    } else {
      plt <- plt +
        ggplot2::scale_fill_gradientn(
          colors = c("#440000", "#ff5500", "white"), 
          values = c(0, 0.4, 1),
          name = "RV",
          limits = c(0, 1.1), )
    }
  } else { # bliss/hsa excess matrix
    plt <- plt +
      ggplot2::scale_fill_gradientn(
        colors = c("black", "#ffffaa", "white", "white", "#aaffff", "blue"),
        values = c(0, 0.35, 0.48, 0.51, 0.65, 1), 
        limits = c(-0.6, 0.6), 
        name = gDRutils::get_combo_excess_field_names()[[mx_name]],
        oob = scales::squish)
  }

    plt <- plt +
      ggplot2::geom_path(data = dt_iso, linewidth = 0.8,
                         ggplot2::aes(x = pos_x, y = pos_y, color = iso_level)) + 
      ggplot2::scale_color_manual(values = iso_colors[all_iso],
                                  name = "Iso levels")

  
  return(plt)
})


# compare iso levels ----
# isobolograms across range of concentration ratios
dt_iso <- gDRutils::convert_se_assay_to_dt(se1, "isobolograms")
dt_iso <- dt_iso[normalization_type == norm_type, ][, .(iso_level, log10_ratio_conc, log2_CI)]

all_iso <- unique(dt_iso$iso_level)
iso_colors <- gDRutils::get_iso_colors()[all_iso]

# base plot
plt_iso_compare <- 
  ggplot2::ggplot(mapping = ggplot2::aes(x = log10_ratio_conc, y = log2_CI)) +
  ggplot2::geom_line(data = data.frame(log10_ratio_conc = c(-2, 2), log2_CI = c(0, 0))) +
  ggplot2::geom_hline(yintercept = 0)

for (iso in all_iso) {
  plt_iso_compare <- plt_iso_compare +
    ggplot2::geom_path(data = dt_iso[iso_level == iso, ], color = iso_colors[iso])
}

plt_iso_compare <- plt_iso_compare  +
  ggplot2::scale_y_continuous(breaks = -5:4, labels = c(paste0("1/", 2 ** (5:1)), 2 ** (0:4))) +
  ggplot2::scale_x_continuous(breaks = -3:3, labels = c(paste0("1/", 10 ** (3:1)), 10 ** (0:3))) +
  ggplot2::coord_cartesian(ylim = c(-5, 4)) +
  ggplot2::ylab("CI") +
  ggplot2::xlab(paste(drug2_name, "/", drug1_name, "ratio")) +
  ggplot2::theme_bw()

if (length(all_iso) > 1) {
  plt_iso_compare <- plt_iso_compare +
    ggplot2::scale_color_gradientn(
      colors = iso_colors,
      values = scales::rescale(as.numeric(names(iso_colors))),
      breaks = quantile(as.numeric(names(iso_colors)), seq(0, 1, 0.25)),
      labels = 100 - round(100 * (quantile(as.numeric(names(iso_colors)), seq(0, 1, 0.25)))),
      name = ifelse(norm_type == "GR", "GR", "IC"))
}
