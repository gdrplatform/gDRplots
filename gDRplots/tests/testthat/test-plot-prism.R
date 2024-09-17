context("Test plot-prism")

# data ----
mae <- gDRutils::get_synthetic_data("combo_matrix")
se_sa <- mae[[gDRutils::get_supported_experiments("sa")]]
dt_metrics <- gDRutils::convert_se_assay_to_dt(se = se_sa,
                                               assay_name = "Metrics")
dt_average <- gDRutils::convert_se_assay_to_dt(se = se_sa,
                                               assay_name = "Averaged")

d_name <- "drug_002"
dt_response_met <- prep_dt_response_metric_sa(dt_metrics, d_name,
                                              metric = c("xc50", "x_mean", "x_max"))
dt_response_dose <- prep_dt_response_dose_sa(dt_average, d_name)

# fake depmap data
cell_lines <- gDRtestData::create_synthetic_cell_lines()[["CellLineName"]]
drugs <- gDRtestData::create_synthetic_drugs()[["DrugName"]]

dt_depmap_meta_lng <- data.table::data.table(
  CCLEName = cell_lines,
  meta_xx = withr::with_seed(42, sample(sprintf("meta_%s", c("AA", "BB", "CC")), size = NROW(cell_lines), replace = TRUE))
)
dt_depmap_meta_lng[CCLEName %in% c("cellline_OO", "cellline_AA"), ][["meta_xx"]] <- NA
dt_depmap_meta_lng[CCLEName == "cellline_FD", ][["meta_xx"]] <- "meta_DD"
dt_depmap_meta_lng[CCLEName == "cellline_NE", ][["meta_xx"]] <- "longer_than_other_meta_EE"

dt_depmap_meta <- data.table::dcast(data = dt_depmap_meta_lng, 
                                    formula = CCLEName ~ meta_xx, 
                                    fun.aggregate = length)
dt_depmap_meta$ModelID = sprintf("ACH-%06d", seq_along(cell_lines))
obj_depmap_meta <- list(
  dt_depmap = dt_depmap_meta,
  selected_feat_meta_col = "fake_mata_data"
)

# tests ----
test_that("plot_volcano_assoc works as expected", {
})

test_that("plot_scatter_with_corr works as expected", {
}) 

test_that("plot_boxplot_meta works as expected", {
  selected_metrics <- "RV_gDR_xc50"
  selected_meta <- "meta_xx"
  dt_response <- dt_response_met[, c("rId", "cId", "CellLineName", selected_metrics), with = FALSE]
  
  grp_stat <- dt_depmap_meta_lng[CCLEName %in% dt_response$CellLineName, .N, by = meta_xx]
  
  plt_1 <- plot_boxplot_meta(dt_response = dt_response,
                             dt_depmap_lng = dt_depmap_meta_lng, 
                             selected_meta = selected_meta)
  expect_is(plt_1, "gg")
  expect_equal(plt_1[["labels"]][["y"]], selected_metrics)
  expect_equal(plt_1[["labels"]][["title"]], selected_meta)
  expect_length(plt_1[["layers"]], 4)
  expect_length(ggplot2::ggplot_build(plt_1)$data[[3]]$xid,
                NROW(grp_stat[!is.na(meta_xx)]))
  expect_equal(sort(ggplot2::layer_scales(plt_1)$x$range$range), 
               sort(grp_stat[!is.na(meta_xx)]$meta_xx))
  
  plt_2 <- plot_boxplot_meta(dt_response = dt_response,
                             dt_depmap_lng = dt_depmap_meta_lng, 
                             selected_meta = selected_meta,
                             with_1_item_grp = FALSE)
  expect_is(plt_2, "gg")
  expect_length(plt_2[["layers"]], 4)
  expect_length(ggplot2::ggplot_build(plt_2)$data[[3]]$xid,
                NROW(grp_stat[!is.na(meta_xx) & N > 1]))
  expect_equal(sort(ggplot2::layer_scales(plt_2)$x$range$range), 
               sort(grp_stat[!is.na(meta_xx) & N > 1]$meta_xx))
  
  plt_3 <- plot_boxplot_meta(dt_response = dt_response,
                             dt_depmap_lng = dt_depmap_meta_lng, 
                             selected_meta = selected_meta,
                             max_x_lbl_length = 8)
  expect_is(plt_3, "gg")
  expect_length(plt_3[["layers"]], 4)
  expect_length(ggplot2::ggplot_build(plt_3)$data[[3]]$xid,
                NROW(grp_stat[!is.na(meta_xx)]))
  ls_x_lbl <- ggplot2::layer_scales(plt_3)$x$labels
  short_lbl <- paste0(substr(grp_stat[!is.na(meta_xx) & nchar(meta_xx) > 8]$meta_xx, 1, 8 - 3), "...")
  expect_true(all(grp_stat[!is.na(meta_xx) & nchar(meta_xx) < 8]$meta_xx %in% ls_x_lbl))
  expect_true(short_lbl %in% ls_x_lbl)
  expect_equal(sort(ggplot2::layer_scales(plt_3)$x$range$range), 
               sort(grp_stat[!is.na(meta_xx)]$meta_xx))
  
}) 
