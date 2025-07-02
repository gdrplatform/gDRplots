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
  dwnlist <- lapply(unique(iris$Species), function(iris_name) {
    file.path("tables", paste0(iris_name, ".xlsx"))
  })
  names(dwnlist) <- unique(iris$Species)
  
  res_1 <- prep_plot_chunk(plt_list = plotlist,
                           chunk_name = "iris")
  expect_is(res_1, "list")
  expect_length(res_1, NROW(plotlist))
  expect_true(all(vapply(seq_along(res_1), function(i) grepl("###", res_1[[i]]), logical(1))))
  expect_true(all(vapply(seq_along(res_1), function(i) {
    grepl(unique(iris$Species)[i], res_1[[i]]) }, logical(1))))
  
  res_2 <- prep_plot_chunk(plt_list = plotlist, 
                           chunk_name = "iris", 
                           link_list = linklist, 
                           dwn_list = dwnlist,
                           header_level = 2)
  expect_true(all(vapply(seq_along(res_2), function(i) grepl("##", res_2[[i]]), logical(1))))
  expect_true(all(vapply(seq_along(res_2), function(i) grepl("_blank", res_2[[i]]), logical(1)))) # link_list
  expect_true(all(vapply(seq_along(res_2), function(i) grepl("download>", res_2[[i]]), logical(1)))) # dwn_list
  
  # scenario: incomplete list of links
  res_3 <- prep_plot_chunk(plt_list = plotlist, 
                           chunk_name = "iris", 
                           link_list = linklist[1:2], 
                           dwn_list = dwnlist,
                           header_level = 2)
  expect_false(all(vapply(seq_along(res_3), function(i) grepl("_blank", res_3[[i]]), logical(1))))
  expect_true(all(vapply(seq_along(res_3), function(i) grepl("download>", res_3[[i]]), logical(1))))
  
  # scenario: plotlist without names
  plotlist_noname <- plotlist
  names(plotlist_noname) <- NULL
  
  res_4 <- prep_plot_chunk(plt_list = plotlist_noname, 
                           chunk_name = "CHUNK")
  expect_is(res_4, "list")
  expect_length(res_4, NROW(plotlist))
  expect_true(all(vapply(seq_along(res_4), function(i) {
    grepl(sprintf("### %s", i), res_4[[i]]) }, logical(1))))
  expect_true(all(vapply(seq_along(res_4), function(i) {
    grepl(sprintf("CHUNK_%s", sprintf("%02d", i)), res_4[[i]], fixed = TRUE) }, logical(1))))
  
  # scenario: all - plotlist, linklist and dwn_list - without names
  linklist_noname <- linklist
  names(linklist_noname) <- NULL
  dwnlist_noname <- dwnlist
  names(dwnlist_noname) <- NULL
  
  res_5 <- prep_plot_chunk(plt_list = plotlist_noname, 
                           link_list = linklist_noname, 
                           dwn_list = dwnlist_noname, 
                           chunk_name = "iris")
  expect_is(res_5, "list")
  expect_length(res_5, NROW(plotlist))
  expect_true(all(vapply(seq_along(res_5), 
                         function(i) grepl(sprintf("### %s", i), res_5[[i]]), logical(1))))
  expect_true(all(vapply(seq_along(res_5), function(i) grepl("_blank", res_5[[i]]), logical(1))))
  expect_true(all(vapply(seq_along(res_5), function(i) grepl("download>", res_5[[i]]), logical(1))))
  
  # scenario: linklist with NA
  dwnlist_NA <- dwnlist
  dwnlist_NA[2] <- NA
  res_6 <- prep_plot_chunk(plt_list = plotlist, 
                           chunk_name = "iris", 
                           link_list = linklist, 
                           dwn_list = dwnlist_NA,
                           header_level = 2)
  expect_is(res_6, "list")
  expect_length(res_6, NROW(plotlist))
  expect_true(all(vapply(seq_along(res_6), function(i) {
    grepl(sprintf("## %s", names(plotlist)[i]), res_6[[i]]) }, logical(1))))
  expect_true(all(vapply(seq_along(res_6), function(i) {
    grepl(sprintf("iris_%s", sprintf("%02d", i)), res_6[[i]], fixed = TRUE) }, logical(1))))
  expect_true(all(vapply(seq_along(res_6), function(i) grepl("_blank", res_6[[i]]), logical(1))))
  expect_equal(sum(vapply(seq_along(res_6), function(i) {
    grepl("download>", res_6[[i]]) }, logical(1))), 2)
  
  # scenario: plot names has spaces
  plotlist_lng <- lapply(sort(unique(iris$Sepal.Length)), function(sl) {
    ggplot2::ggplot(iris[iris$Sepal.Length == sl, c("Petal.Length", "Petal.Width", "Species")]) +
      ggplot2::geom_point(ggplot2::aes(x = Petal.Length, y = Petal.Width, color = Species))
  })
  names(plotlist_lng) <- sprintf("Sepal.Length %s (%s of %s)", 
                                 sort(unique(iris$Sepal.Length)), 
                                 seq_along(names(plotlist_lng)), 
                                 rep(NROW(plotlist_lng), NROW(plotlist_lng)))
  
  res_7 <- prep_plot_chunk(plt_list = plotlist_lng, 
                           chunk_name = "iris",
                           header_level = 2)
  expect_is(res_7, "list")
  expect_length(res_7, NROW(plotlist_lng))
  expect_true(all(vapply(seq_along(res_7), function(i) { 
    grepl(sprintf("## %s", names(plotlist_lng)[i]), res_7[[i]], fixed = TRUE) }, logical(1))))
  expect_true(all(vapply(seq_along(res_7), function(i) {
    grepl(sprintf("iris_%s", sprintf("%03d", i)), res_7[[i]], fixed = TRUE) }, logical(1))))
  
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
  expect_error(prep_plot_chunk(plt_list = plotlist, 
                               chunk_name = "iris", 
                               tabset_options = 2), 
               "Assertion on 'tabset_options' failed: Must be of type 'character'")
  
  # nested plotlist
  plotlist_nest <- list(someCategory = c(plotlist), anotherCategory = c(plotlist))
  linklist_nest <- list(someCategory = c(linklist), anotherCategory = c(linklist))
  
  res_1n <- prep_plot_chunk(plt_list = plotlist_nest, 
                            chunk_name = "iris") # default
  expect_is(res_1n, "list")
  expect_length(res_1n, NROW(plotlist_nest))
  expect_true(all(vapply(res_1n, function(i) is.character(i), logical(1))))
  expect_length(res_1n, NROW(plotlist_nest))
  expect_equal(unlist(lapply(seq_along(res_1n), function(i) sum(grepl("####", res_1n[[i]])))), 
               c(NROW(plotlist_nest[[1]]), NROW(plotlist_nest[[2]])))
  expect_equal(sum(unlist(lapply(seq_along(res_1n), 
                                 function(i) grepl("\\{.tabset .tabset-dropdown\\}", res_1n[[i]])))),
               NROW(plotlist_nest))
  
  res_2n <- prep_plot_chunk(plt_list = plotlist_nest, 
                            link_list = linklist_nest,
                            chunk_name = "iris",
                            tabset_options = NULL)
  expect_equal(unlist(lapply(seq_along(res_2n), function(i) sum(grepl("####", res_2n[[i]])))), 
               c(NROW(plotlist_nest[[1]]), NROW(plotlist_nest[[2]])))
  expect_equal(unlist(lapply(seq_along(res_2n), function(i) sum(grepl("a href", res_2n[[i]])))), 
               c(NROW(linklist_nest[[1]]), NROW(linklist_nest[[2]])))
  expect_equal(
    vapply(seq_along(res_2n), function(i) {
      sum(grepl(sprintf("### %s\n\n", names(plotlist_nest)[i]), res_2n[[i]])) }, numeric(1)),
    c(1, 1)) # headers
  
  # scenario: incomplete list of links
  linklist_nest_incom <- list(someCategory = c(linklist[2:3]), anotherCategory = c(linklist))
  res_3n <- prep_plot_chunk(plt_list = plotlist_nest, 
                            link_list = linklist_nest_incom,
                            chunk_name = "iris")
  expect_equal(unlist(lapply(seq_along(res_3n), function(i) sum(grepl("####", res_3n[[i]])))), 
               c(NROW(plotlist_nest[[1]]), NROW(plotlist_nest[[2]])))
  expect_equal(unlist(lapply(seq_along(res_3n), function(i) sum(grepl("a href", res_3n[[i]])))), 
               c(0, 0))
  
  # scenario: plotlist without names (partially)
  plotlist_nest_noname <- plotlist_nest
  names(plotlist_nest_noname) <- NULL
  res_4n <- prep_plot_chunk(plt_list = plotlist_nest, 
                            link_list = linklist_nest,
                            chunk_name = "iris",
                            tabset_options = "unnumbered")
  expect_equal(unlist(lapply(seq_along(res_4n), function(i) sum(grepl("####", res_4n[[i]])))), 
               c(NROW(plotlist_nest_noname[[1]]), NROW(plotlist_nest_noname[[2]])))
  expect_equal(unlist(lapply(seq_along(res_4n), function(i) sum(grepl("a href", res_4n[[i]])))), 
               c(NROW(linklist_nest[[1]]), NROW(linklist_nest[[2]])))
  expect_equal(sum(unlist(lapply(seq_along(res_4n), 
                                 function(i) grepl("\\{.unnumbered\\}", res_4n[[i]])))),
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
  linklist <- list()
  ls_color <- c("darkred", "orange", "darkcyan")
  for (drug in unique(dt_metrics$DrugName)) {
    for (cl in unique(dt_metrics$CellLineName)) {
      tab_plot <- dt_metrics[DrugName == drug & CellLineName == cl]
      
      plt_GR <- lapply(ls_color, function(col) plot_col(tab_plot, "GRV", col))
      names(plt_GR) <- sprintf("%s_%s", "GR", ls_color)
      plt_RV <- lapply(ls_color, function(col) plot_col(tab_plot, "RV", col))
      names(plt_RV) <- sprintf("%s_%s", "RV", ls_color)
      
      plotlist[[drug]][[cl]][["GR"]] <- plt_GR
      plotlist[[drug]][[cl]][["RV"]] <- plt_RV
    }
  }
  
  for (drug in unique(dt_metrics$DrugName)) {
    for (cl in unique(dt_metrics$CellLineName)) {
      link_GR <- lapply(ls_color, function(col) {
        name_GR <- sprintf("%s_%s", "GR", col)
        file.path("plot", paste0(name_GR, ".png"))
      })
      names(link_GR) <- sprintf("%s_%s", "GR", ls_color)
      link_RV <- lapply(ls_color, function(col) {
        name_RV <- sprintf("%s_%s", "RV", col)
        file.path("plot", paste0(name_RV, ".png"))
      })
      names(link_RV) <- sprintf("%s_%s", "RV", ls_color)
      
      linklist[[drug]][[cl]][["GR"]] <- link_GR
      linklist[[drug]][[cl]][["RV"]] <- link_RV
    }
  }  
  no_plots_in_section <-  NROW(unique(dt_metrics$CellLineName)) * NROW(c("GR", "RV")) * NROW(ls_color)
  
  res_1 <- prep_nested_plot_chunk(plt_list = plotlist, 
                                  chunk_name = "metric_col") # default
  expect_is(res_1, "list")
  expect_length(res_1, NROW(plotlist))
  expect_length(res_1, NROW(unique(dt_metrics$DrugName)))
  expect_equal(sum(grepl("#####", res_1[[1]])), no_plots_in_section) # the lowest lvl with plots - default
  
  chunk_name_2 <- "CHUNK"
  res_2 <- prep_nested_plot_chunk(plt_list = plotlist, 
                                  chunk_name = chunk_name_2, 
                                  link_list = linklist,
                                  header_level = 1)
  expect_is(res_2, "list")
  expect_equal(sum(vapply(seq_along(res_2), function(i) grepl("# drug_", res_2[[i]][1]), logical(1))),
               NROW(res_2))
  expect_equal(sum(grepl(chunk_name_2, res_2[[1]])), no_plots_in_section)
  expect_equal(sum(grepl("####", res_2[[1]])), no_plots_in_section)
  expect_true(all(vapply(seq_along(res_2), 
                         function(i) sum(grepl("_blank", res_2[[i]])) == no_plots_in_section, logical(1)))) # link_list
  
  linklist_noname <- linklist[1:5]
  names(linklist_noname) <- NULL
  res_3 <- prep_nested_plot_chunk(plt_list = plotlist[1:5], 
                                  chunk_name = "metric_col", 
                                  link_list = linklist_noname,
                                  header_level = 2)
  expect_is(res_3, "list")
  expect_equal(sum(vapply(seq_along(res_3), function(i) grepl("## drug_", res_3[[i]][1]), logical(1))),
               NROW(res_3))
  expect_equal(sum(grepl("#####", res_3[[1]])), no_plots_in_section)
  expect_false(all(vapply(seq_along(res_3), 
                          function(i) sum(grepl("_blank", res_3[[i]])) == no_plots_in_section, logical(1)))) # link_list
  
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

test_that("prep_double_table_chunk works as expected", {
  nested_tables <- list(
    CellLine1 = list(MetricA = mtcars[1:5, ], 
                     MetricB = mtcars[6:10, ]),
    CellLine2 = list(MetricC = iris[1:5, ], 
                     MetricD = iris[6:10, ])
  )
  
  download_link <- lapply(names(nested_tables), function(nm) {
    file.path("tables", paste0("cgs_tables_RV__", nm, ".xlsx"))
  })
  names(download_link) <- names(nested_tables)
  
  res_1 <- prep_double_table_chunk(tbl_list = nested_tables, 
                                   chunk_name = "nested_tables")
  expect_is(res_1, "list")
  expect_length(res_1, NROW(nested_tables))
  expect_true(all(vapply(seq_along(res_1), function(i) all(grepl("###", res_1[[i]])), logical(1))))
  expect_equal(
    vapply(seq_along(res_1), 
           function(i) sum(grepl(sprintf("### %s", names(nested_tables)[i]), res_1[[i]])), numeric(1)),
    c(1, 1)) # headers
  expect_equal(
    vapply(seq_along(res_1),
           function(i) sum(grepl("DT::formatRound", res_1[[i]])), numeric(1)),
    vapply(seq_along(nested_tables), function(i) NROW(nested_tables[[i]]), numeric(1)))
  
  res_2 <- prep_double_table_chunk(tbl_list = nested_tables, 
                                   chunk_name = "nested_tables", 
                                   dwn_list = download_link,
                                   header_level = 4, 
                                   tabset_options = "tabset")
  expect_is(res_2, "list")
  expect_true(all(vapply(seq_along(res_1), function(i) all(grepl("####", res_2[[i]])), logical(1))))
  expect_equal(
    vapply(seq_along(res_2), 
           function(i) sum(grepl("\\{\\.tabset\\}", res_2[[i]])), numeric(1)), c(1, 1)) # headers
  expect_equal(
    vapply(seq_along(res_2), 
           function(i) sum(grepl("download>", res_2[[i]])), numeric(1)), c(1, 1)) # dwn_list
  
  download_link_noname <- download_link
  names(download_link_noname) <- NULL
  res_3 <- prep_double_table_chunk(nested_tables, 
                                   chunk_name = "dt_tables", 
                                   dwn_list = download_link_noname,
                                   header_level = 2, 
                                   tabset_options = NULL)
  expect_is(res_3, "list")
  expect_equal(
    vapply(seq_along(res_3), 
           function(i) sum(grepl(sprintf("## %s \n\n", names(nested_tables)[i]), res_3[[i]])), numeric(1)),
    c(1, 1)) # headers
  expect_equal(
    vapply(seq_along(res_3), 
           function(i) sum(grepl("download>", res_3[[i]])), numeric(1)), c(0, 0)) # lack of dwn_list
  
  # sorting opts
  res_4 <- prep_double_table_chunk(tbl_list = nested_tables, 
                                   chunk_name = "nested_tables",
                                   sorting_opts = c("Sepal.Length", "-gear"))
  expect_is(res_4, "list")
  expect_true(all(grepl("order", res_4)))
  expect_true(all(grepl(sprintf('(%sL, "desc")', which(colnames(mtcars) == "gear")), res_4[[1]][2:3])))
  expect_true(all(grepl(sprintf('(%sL, "asc")', which(colnames(iris) == "Sepal.Length")), res_4[[2]][2:3])))
  
  res_5 <- prep_double_table_chunk(tbl_list = nested_tables, 
                                   chunk_name = "nested_tables",
                                   sorting_opts = c("-non_existen_col"))
  expect_is(res_5, "list")
  expect_false(all(grepl("order", res_5)))
  
  expect_error(prep_double_table_chunk(tbl_list = data.table::data.table(iris), chunk_name = "iris"), 
               "Assertion on 'tbl_list' failed: Must be of type 'list'")
  expect_error(prep_double_table_chunk(tbl_list = nested_tables, chunk_name = 123), 
               "Assertion on 'chunk_name' failed: Must be of type 'string'")
  expect_error(prep_double_table_chunk(tbl_list = nested_tables, 
                                       chunk_name = "nested_tables",
                                       dwn_list = unique(iris$Species)), 
               "Assertion on 'dwn_list' failed: Must be of type 'list'")
  expect_error(prep_double_table_chunk(tbl_list = nested_tables, 
                                       chunk_name = "nested_tables",
                                       header_level = "1"), 
               "Assertion on 'header_level' failed: Must be of type 'single integerish value'")
  expect_error(prep_double_table_chunk(tbl_list = nested_tables, 
                                       chunk_name = "nested_tables",
                                       tabset_options = 2), 
               "Assertion on 'tabset_options' failed: Must be of type 'character'")
  expect_error(prep_double_table_chunk(tbl_list = nested_tables, 
                                       chunk_name = "nested_tables",
                                       sorting_opts = 2), 
               "Assertion on 'sorting_opts' failed: Must be of type 'character'")
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
  expect_true(grepl("_blank", res_2))
  expect_true(grepl(zoom_txt, res_2))
  
  res_3 <- create_zoom_link(img_path = NA, link_txt = zoom_txt)
  expect_equal(res_3, "")
  expect_false(grepl("_blank", res_3))
  expect_false(grepl(zoom_txt, res_3))
  
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
  
  res_3 <- create_download_link(dwn_path = NA, link_txt = dwn_txt)
  expect_equal(res_3, "")
  expect_false(grepl("download", res_3))
  expect_false(grepl(dwn_txt, res_3))
  
  expect_error(create_download_link(dwn_path = 1),
               "Assertion on 'dwn_path' failed: Must be of type 'string'")
  expect_error(create_download_link(dwn_path = c("A", "B")),
               "Assertion on 'dwn_path' failed: Must have length 1")
  expect_error(create_download_link(dwn_path = file_path,
                                    link_txt = 123),
               "Assertion on 'link_txt' failed: Must be of type 'string'")
})

test_that("prep_filename_path works as expected", {
  # simple list
  plotlist <- lapply(unique(iris$Species), function(iris_name) {
    ggplot2::ggplot(iris[iris$Species == iris_name, c("Sepal.Length", "Sepal.Width")]) +
      ggplot2::geom_point(ggplot2::aes(x = Sepal.Length, y = Sepal.Width))
  })
  names(plotlist) <- unique(iris$Species)
  
  
  res_1 <- prep_filename_path(plt_list = plotlist) # default 
  expect_equal(names(res_1), names(plotlist))
  expect_equal(unlist(res_1, use.names = FALSE), names(plotlist))
  
  prefix_i <- "iris__"
  path_i <- file.path(".", "plots")
  format_i <- "png"
  res_2 <- prep_filename_path(plt_list = plotlist,
                              prefix = prefix_i,
                              path_file = path_i)
  expect_equal(names(res_2), names(plotlist))
  expect_true(all(vapply(seq_along(res_2), function(i) grepl(prefix_i, res_2[[i]]), logical(1))))
  expect_true(all(vapply(seq_along(res_2), function(i) grepl(path_i, res_2[[i]]), logical(1))))
  
  # scenario: list without name
  noname_plotlist <- plotlist
  names(noname_plotlist) <- NULL
  res_3 <- prep_filename_path(plt_list = noname_plotlist,
                              prefix = prefix_i,
                              file_format = format_i)
  expect_equal(names(res_3), as.character(seq_along(plotlist)))
  expect_true(all(vapply(seq_along(res_3), function(i) grepl(prefix_i, res_3[[i]]), logical(1))))
  expect_true(all(vapply(seq_along(res_3), function(i) grepl(format_i, res_3[[i]]), logical(1))))
  
  # nested list
  nested_plotlist <- list()
  for (species in unique(iris$Species)) {
    nested_plotlist[[species]] <- list()
    nested_plotlist[[species]][["Sepal"]] <-
      ggplot2::ggplot(iris[iris$Species == species, ],
                      ggplot2::aes(x = Sepal.Length, y = Sepal.Width)) + ggplot2::geom_point()
    nested_plotlist[[species]][["Petal"]] <-
      ggplot2::ggplot(iris[iris$Species == species, ],
                      ggplot2::aes(x = Petal.Length, y = Petal.Width)) + ggplot2::geom_point()
  }
  
  res_4 <- prep_filename_path(plt_list = nested_plotlist)
  expect_equal(names(res_4), names(nested_plotlist))
  expect_true(all(
    vapply(seq_along(res_4), 
           function(i) all(names(res_4[[i]]) == names(nested_plotlist[[i]])), logical(1))))
  expect_true(all(
    vapply(seq_along(res_4), 
           function(i) all(names(res_4[[i]]) == unlist(res_4[[i]], use.names = FALSE)), logical(1))))
  
  path_t <- file.path(".", "tables")
  format_t <- "xlsx"
  res_5 <- prep_filename_path(plt_list = nested_plotlist,
                              path_file = path_t,
                              file_format = format_t)
  expect_equal(names(res_5), names(nested_plotlist))
  expect_true(all(vapply(seq_along(unlist(res_5)), function(i) grepl(path_t, unlist(res_5)[[i]]), logical(1))))
  expect_true(all(vapply(seq_along(unlist(res_5)), function(i) grepl(format_t, unlist(res_5)[[i]]), logical(1))))
  
  # scenario: list without name
  noname_nested_plotlist <- nested_plotlist
  names(noname_nested_plotlist) <- NULL
  res_6 <- prep_filename_path(plt_list = noname_nested_plotlist,
                              prefix = prefix_i)
  expect_equal(names(res_6), as.character(seq_along(nested_plotlist)))
  expect_true(all(vapply(seq_along(unlist(res_6)), function(i) grepl(prefix_i, unlist(res_6)[[i]]), logical(1))))
  
  # scenario: nested level without names
  noname_nested_plotlist_2 <- nested_plotlist
  names(noname_nested_plotlist_2[[3]]) <- NULL
  names(noname_nested_plotlist_2) <- LETTERS[1:3]
  res_7 <- prep_filename_path(plt_list = noname_nested_plotlist_2)
  expect_equal(names(res_7), names(noname_nested_plotlist_2))
  expect_true(all(
    vapply(seq_along(res_7), 
           function(i) all(names(res_7[[i]]) == as.character(names(noname_nested_plotlist_2[[i]]))), logical(1))))
  
  expect_error(prep_filename_path(plt_list = data.table::data.table()),
               "Assertion on 'plt_list' failed: Must be of type 'list'")
  expect_error(prep_filename_path(plt_list = plotlist,
                                  prefix = 1),
               "Assertion on 'prefix' failed: Must be of type 'string'")
  expect_error(prep_filename_path(plt_list = plotlist,
                                  path_file = file.path()),
               "Assertion on 'path_file' failed: Must have length 1.")
  expect_error(prep_filename_path(plt_list = plotlist,
                                  path_file = 123),
               "Assertion on 'path_file' failed: Must be of type 'string'")
  expect_error(prep_filename_path(plt_list = plotlist,
                                  file_format = TRUE),
               "Assertion on 'file_format' failed: Must be of type 'string'")
})


test_that("generate_datatable works as expected", {
  result_df <- generate_datatable(iris)
  expect_s3_class(result_df, "datatables")
  
  expect_error(
    generate_datatable(matrix(1:10, ncol = 2)), 
    "Assertion failed"
  )
  expect_error(
    generate_datatable(iris, options = "invalid_options"),
    "Must be of type 'list'"
  )
  expect_error(
    generate_datatable(iris, width = 100), 
    "Must be of type 'string'"
  )
  result_custom_options <- generate_datatable(iris, options = list(scrollX = TRUE, pageLength = 5))
  expect_s3_class(result_custom_options, "datatables")
  result_with_caption <- generate_datatable(iris, caption = "Iris Dataset")
  expect_s3_class(result_with_caption, "datatables")
})
