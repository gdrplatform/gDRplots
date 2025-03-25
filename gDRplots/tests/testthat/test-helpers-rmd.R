context("Test helpers-rmd")

test_that("prep_plot_chunk works as expected", {
  plotlist <- lapply(unique(iris$Species), function(iris_name) {
    ggplot2::ggplot(iris[iris$Species == iris_name, c("Sepal.Length", "Sepal.Width")]) +
      ggplot2::geom_point(ggplot2::aes(x = Sepal.Length, y = Sepal.Width))
  })
  names(plotlist) <- unique(iris$Species)
  linklist <- lapply(unique(iris$Species), function(iris_name) {
    file.path("plot", paste0(iris_name, ".png"))
  })
  names(linklist) <- unique(iris$Species)
  
  res_1 <- prep_plot_chunk(plt_list = plotlist, chunk_name = "iris")
  expect_is(res_1, "list")
  expect_length(res_1, NROW(plotlist))
  expect_true(all(vapply(seq_along(res_1), function(i) grepl("###", res_1[[i]]), logical(1))))
  expect_true(all(vapply(seq_along(res_1), 
                         function(i) grepl(unique(iris$Species)[i], res_1[[i]]), logical(1))))
  
  res_2 <- prep_plot_chunk(plt_list = plotlist, chunk_name = "iris", 
                           link_list = linklist, header_level = 2)
  expect_true(all(vapply(seq_along(res_2), function(i) grepl("##", res_2[[i]]), logical(1))))
  expect_true(all(vapply(seq_along(res_2), function(i) grepl("a href", res_2[[i]]), logical(1))))
  
  # scenario: incomplete list of links
  res_3 <- prep_plot_chunk(plt_list = plotlist, chunk_name = "iris", 
                           link_list = linklist[1:2], header_level = 2)
  expect_false(all(vapply(seq_along(res_3), function(i) grepl("a href", res_3[[i]]), logical(1))))
  
  # scenario: plotlist without names
  plotlist_noname <- plotlist
  names(plotlist_noname) <- NULL
  
  res_4 <- prep_plot_chunk(plt_list = plotlist_noname, chunk_name = "iris")
  expect_is(res_4, "list")
  expect_length(res_4, NROW(plotlist))
  expect_true(all(vapply(seq_along(res_4), 
                         function(i) grepl(sprintf("### %s", i), res_4[[i]]), logical(1))))
  
  # scenario: both - plotlist and linklist - without names
  linklist_noname <- linklist
  names(linklist_noname) <- NULL
  
  res_5 <- prep_plot_chunk(plt_list = plotlist_noname, 
                           link_list = linklist_noname, 
                           chunk_name = "iris")
  expect_is(res_5, "list")
  expect_length(res_5, NROW(plotlist))
  expect_true(all(vapply(seq_along(res_5), 
                         function(i) grepl(sprintf("### %s", i), res_5[[i]]), logical(1))))
  expect_true(all(vapply(seq_along(res_5), function(i) grepl("a href", res_5[[i]]), logical(1))))
  
  expect_error(prep_plot_chunk(plt_list = iris, chunk_name = "iris"), 
               "Assertion on 'plt_list' failed: Must be of type 'list'")
  expect_error(prep_plot_chunk(plt_list = plotlist, chunk_name = 123), 
               "Assertion on 'chunk_name' failed: Must be of type 'string'")
  expect_error(prep_plot_chunk(plt_list = plotlist, 
                               chunk_name = "iris", 
                               link_list = unique(iris$Species)), 
               "Assertion on 'link_list' failed: Must be of type 'list'")
  expect_error(prep_plot_chunk(plt_list = plotlist, 
                               chunk_name = "iris", 
                               header_level = "1"), 
               "Assertion on 'header_level' failed: Must be of type 'single integerish value'")
  
  # nested plotlist
  plotlist_nest <- list(someCategory = c(plotlist), anotherCategory = c(plotlist))
  linklist_nest <- list(someCategory = c(linklist), anotherCategory = c(linklist))
  
  res_6 <- prep_plot_chunk(plt_list = plotlist_nest, chunk_name = "iris")
  expect_is(res_6, "list")
  expect_length(res_6, NROW(plotlist_nest))
  expect_true(all(vapply(res_6, function(i) is.character(i), logical(1))))
  expect_length(res_6, NROW(plotlist_nest))
  expect_equal(unlist(lapply(seq_along(res_6), function(i) sum(grepl("####", res_6[[i]])))), 
               c(NROW(plotlist_nest[[1]]), NROW(plotlist_nest[[2]])))
  expect_equal(sum(unlist(lapply(seq_along(res_6), 
                                 function(i) grepl("\\{.tabset .tabset-dropdown\\}", res_6[[i]])))),
               NROW(plotlist_nest))
  
  res_7 <- prep_plot_chunk(plt_list = plotlist_nest, 
                           link_list = linklist_nest,
                           chunk_name = "iris")
  expect_equal(unlist(lapply(seq_along(res_7), function(i) sum(grepl("####", res_7[[i]])))), 
               c(NROW(plotlist_nest[[1]]), NROW(plotlist_nest[[2]])))
  expect_equal(unlist(lapply(seq_along(res_7), function(i) sum(grepl("a href", res_7[[i]])))), 
               c(NROW(linklist_nest[[1]]), NROW(linklist_nest[[2]])))
  
  # scenario: incomplete list of links
  linklist_nest_incom <- list(someCategory = c(linklist[2:3]), anotherCategory = c(linklist))
  res_8 <- prep_plot_chunk(plt_list = plotlist_nest, 
                           link_list = linklist_nest_incom,
                           chunk_name = "iris")
  expect_equal(unlist(lapply(seq_along(res_8), function(i) sum(grepl("####", res_8[[i]])))), 
               c(NROW(plotlist_nest[[1]]), NROW(plotlist_nest[[2]])))
  expect_equal(unlist(lapply(seq_along(res_8), function(i) sum(grepl("a href", res_8[[i]])))), 
               c(0, 0))
  
  # scenario: plotlist without names (partially)
  plotlist_nest_noname <- plotlist_nest
  names(plotlist_nest_noname) <- NULL
  res_9 <- prep_plot_chunk(plt_list = plotlist_nest, 
                           link_list = linklist_nest,
                           chunk_name = "iris",
                           tabset_options = "unnumbered")
  expect_equal(unlist(lapply(seq_along(res_9), function(i) sum(grepl("####", res_9[[i]])))), 
               c(NROW(plotlist_nest_noname[[1]]), NROW(plotlist_nest_noname[[2]])))
  expect_equal(unlist(lapply(seq_along(res_9), function(i) sum(grepl("a href", res_9[[i]])))), 
               c(NROW(linklist_nest[[1]]), NROW(linklist_nest[[2]])))
  expect_equal(sum(unlist(lapply(seq_along(res_9), 
                                 function(i) grepl("\\{.unnumbered\\}", res_9[[i]])))),
               NROW(plotlist_nest_noname))
})

test_that("prep_nested_plot_chunk works as expected", {
  mae <- gDRutils::get_synthetic_data("small")
  se <- mae[[gDRutils::get_supported_experiments("sa")]]
  dt_metrics <- gDRutils::convert_se_assay_to_dt(se, "Metrics")
  
  # help function
  plot_col <- function(tab_plt, norm_type, col = "red") {
    tab_plt <- data.table::melt(
      data = tab_plt[normalization_type == norm_type][, c("rId", "xc50", "x_mean", "x_max")],
      id = "rId")
    plt <- ggplot2::ggplot(tab_plt, ggplot2::aes(x = variable, y = value)) +
      ggplot2::geom_col(fill = col)
    return(plt)
  }
  
  # creating nested list with plots
  plotlist <- list()
  ls_color <- c("darkred", "orange", "darkcyan")
  for (drug in unique(dt_metrics$DrugName)) {
    for (cl in unique(dt_metrics$CellLineName)) {
      tab_plot <- dt_metrics[DrugName == drug & CellLineName == cl]
      
      plt_GR <- lapply(ls_color, function(col) plot_col(tab_plot, "RV", col))
      names(plt_GR) <- sprintf("%s_%s", "GR", ls_color)
      plt_RV <- lapply(ls_color, function(col) plot_col(tab_plot, "RV", col))
      names(plt_RV) <- sprintf("%s_%s", "RV", ls_color)
      
      plotlist[[drug]][[cl]][["RV"]] <- plt_RV
      plotlist[[drug]][[cl]][["GR"]] <- plt_GR
    }
  }
  
  res_1 <- prep_nested_plot_chunk(plt_list = plotlist, 
                                  chunk_name = "metric_col")
  expect_is(res_1, "list")
  expect_length(res_1, NROW(plotlist))
  expect_length(res_1, NROW(unique(dt_metrics$DrugName)))
  expect_equal(sum(grepl("#####", res_1[[1]])), # the lowest lvl with plots - default
               NROW(unique(dt_metrics$CellLineName)) * NROW(c("GR", "RV")) * NROW(ls_color))
  
  res_2 <- prep_nested_plot_chunk(plt_list = plotlist, 
                                  chunk_name = "metric_col", 
                                  header_level = 1)
  expect_is(res_2, "list")
  expect_equal(sum(grepl("####", res_2[[1]])),
               NROW(unique(dt_metrics$CellLineName)) * NROW(c("GR", "RV")) * NROW(ls_color))
  
  expect_error(prep_nested_plot_chunk(plt_list = dt_metrics, chunk_name = "metric_col"), 
               "Assertion on 'plt_list' failed: Must be of type 'list'")
  expect_error(prep_nested_plot_chunk(plt_list = plotlist, chunk_name = 123), 
               "Assertion on 'chunk_name' failed: Must be of type 'string'")
  expect_error(prep_nested_plot_chunk(plt_list = plotlist, chunk_name = "metric_col", header_level = "1"), 
               "Assertion on 'header_level' failed: Must be of type 'single integerish value'")
  expect_error(prep_nested_plot_chunk(plt_list = plotlist, chunk_name = "metric_col", header_level = 0), 
               "Assertion on 'header_level' failed: Element 1 is not >= 1")
})

test_that("escape_special_characters works as expected", {
  expect_equal(escape_special_characters("ABC:123"), "ABC[colon]123")
  expect_equal(escape_special_characters("AD_12"), "AD_12")
  expect_equal(escape_special_characters("AD#12"), "AD[hash]12")
  expect_equal(escape_special_characters("AD/12"), "AD[slash]12")
  
  expect_error(escape_special_characters(123), "Assertion on 'x' failed: Must be of type 'string'")
})

test_that("neutralize_spaces works as expected", {
  expect_equal(neutralize_spaces("GDC-123|Abc x G01234"), "GDC-123|Abc_x_G01234")
  expect_equal(neutralize_spaces("MNO-321P 789R YY#1 "), "MNO-321P_789R_YY#1")
  expect_equal(neutralize_spaces("drug_001 x drug_002", replacement = "."), "drug_001.x.drug_002")
  
  expect_error(neutralize_spaces(123), "Assertion on 'x' failed: Must be of type 'string'")
  expect_error(neutralize_spaces("ABC EF", replacement = 1), 
               "Assertion on 'replacement' failed: Must be of type 'string'")
})

# Example ggplot object

p1 <- ggplot2::ggplot(ggplot2::mpg, ggplot2::aes(displ, hwy, color = class)) + ggplot2::geom_point()

# Example pheatmap object
mat <- matrix(rnorm(100), 10, 10)
p2 <- pheatmap::pheatmap(mat)

test_that("estimate_plot_size works for ggplot object", {
  size <- estimate_plot_size(p1)
  expect_named(size, c("width", "height"))
  expect_equal(size[["width"]], 13.5)
  expect_equal(size[["height"]], 9.5)
})

test_that("estimate_plot_size works for pheatmap object", {
  size <- estimate_plot_size(p2)
  expect_named(size, c("width", "height"))
  expect_equal(size[["width"]], 15)
  expect_equal(size[["height"]], 11)
})

test_that("estimate_plot_size handles invalid inputs", {
  invalid_plot <- list()
  expect_error(estimate_plot_size(invalid_plot),
               "Assertion on 'plt' failed: Must inherit from class 'ggplot'/'pheatmap', but has class 'list'.")
  
  invalid_plot_2 <- plot(1:10, 1:10, cex = 2, pch = 20, col = "pink")
  expect_error(estimate_plot_size(invalid_plot_2),
               "Assertion on 'plt' failed: Must inherit from class 'ggplot'/'pheatmap'")
  
  invalid_base_width <- -5
  expect_error(estimate_plot_size(p1, base_width = invalid_base_width),
               "Assertion on 'base_width' failed: Element 1 is not >= 0.")
  
  invalid_base_height <- -5
  expect_error(estimate_plot_size(p1, base_height = invalid_base_height),
               "Assertion on 'base_height' failed: Element 1 is not >= 0.")
  
  invalid_scale_factor <- -0.5
  expect_error(estimate_plot_size(p1, scale_factor = invalid_scale_factor),
               "Assertion on 'scale_factor' failed: Element 1 is not >= 0.")
})


test_that("save_plot saves ggplot2 object in a correct format", {
  temp_dir <- tempdir()
  file_path <- file.path(temp_dir, "test_plot")
  
  save_plot(p1, file_path)
  expect_true(file.exists(paste0(file_path, ".svg")))
  save_plot(p1, file_path, "png")
  expect_true(file.exists(paste0(file_path, ".png")))
  save_plot(p1, file_path, "pdf")
  expect_true(file.exists(paste0(file_path, ".pdf")))
})

test_that("save_plot throws error for unsupported plot type", {
  p <- list()  # Not a ggplot2 or pheatmap object
  temp_dir <- tempdir()
  file_path <- file.path(temp_dir, "test_plot")
  
  expect_error(save_plot(p, file_path, "svg"),
               "Assertion on 'plt' failed: Must inherit from class 'ggplot'/'pheatmap', but has class 'list'.")
})

test_that("save_plot throws error for non-existent directory", {
  p <- ggplot2::ggplot(datasets::mtcars, ggplot2::aes(mpg, wt)) + ggplot2::geom_point()
  file_path <- "non_existent_directory/test_plot"
  
  expect_error(save_plot(p, file_path, "svg"), "The specified directory does not exist.")
})

test_that("get_r_file_path works as expected", {
  
  r_path <- "test-helpers-rmd.R" 
  ca1 <- c("a",
           "b=3",
           paste0("--file=", r_path),
           "--file-path-ext=1234")
  
  a_path <- system.file(package = "gDRplots", "DESCRIPTION")
  ca2 <- c("a", "b=3", paste0("--file=", a_path), "--file-path-ext=1234")
  
  mockery::stub(
    where = get_r_file_path,
    what = "commandArgs",
    how = function() {
      ca1
    },
    depth = 1
  )
  fr_path <- get_r_file_path(test_mode = TRUE)
  checkmate::test_file_exists(fr_path)
  expect_identical(fr_path, tools::file_path_as_absolute(fr_path))
  expect_false(fr_path == r_path)
  
  mockery::stub(
    where = get_r_file_path,
    what = "commandArgs",
    how = function() {
      ca2
    },
    depth = 1
  )
  fa_path <- get_r_file_path(test_mode = TRUE)
  checkmate::test_file_exists(fa_path)
  expect_identical(fa_path, tools::file_path_as_absolute(fa_path))
  expect_identical(fa_path, a_path)
  
  expect_error(get_r_file_path(test_mode = 1),
               "Assertion on 'test_mode' failed")
})

test_that("create_zoom_link works as expected", {
  i_path <- "./folder/file.png"
  zoom_txt <- "Click to see bigger picture"
  
  res_1 <- create_zoom_link(img_path = i_path) # default
  expect_true(grepl(i_path, res_1))
  expect_true(grepl("a href", res_1))
  expect_true(grepl("Zoom In for Details", res_1))
  
  res_2 <- create_zoom_link(img_path = i_path, link_txt = zoom_txt)
  expect_true(grepl(i_path, res_2))
  expect_true(grepl("a href", res_2))
  expect_true(grepl(zoom_txt, res_2))

  expect_error(create_zoom_link(img_path = 1),
               "Assertion on 'img_path' failed: Must be of type 'string'")
  expect_error(create_zoom_link(img_path = c("A", "B")),
               "Assertion on 'img_path' failed: Must have length 1")
  expect_error(create_zoom_link(img_path = i_path,
                                link_txt = 123),
               "Assertion on 'link_txt' failed: Must be of type 'string'")
})

test_that("create_download_link works as expected", {
  file_path <- "./folder/file.xlsx"
  dwn_txt <- "Click to download"
  
  res_1 <- create_download_link(dwn_path = file_path) # default
  expect_true(grepl(file_path, res_1))
  expect_true(grepl("download", res_1))
  expect_true(grepl("Download Table", res_1))
  
  res_2 <- create_download_link(dwn_path = file_path, link_txt = dwn_txt)
  expect_true(grepl(file_path, res_2))
  expect_true(grepl("download", res_2))
  expect_true(grepl(dwn_txt, res_2))
    
  expect_error(create_download_link(dwn_path = 1),
               "Assertion on 'dwn_path' failed: Must be of type 'string'")
  expect_error(create_download_link(dwn_path = c("A", "B")),
               "Assertion on 'dwn_path' failed: Must have length 1")
  expect_error(create_download_link(dwn_path = file_path,
                                    link_txt = 123),
               "Assertion on 'link_txt' failed: Must be of type 'string'")
})
