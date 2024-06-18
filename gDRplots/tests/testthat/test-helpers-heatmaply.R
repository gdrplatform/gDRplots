context("Test heatmaply helpers")


test_that("compute_distances works as expected", {
  set.seed(1234)
  
  # matrix without row names
  y <- matrix(LETTERS[1:9], nrow = 3, ncol = 3, dimnames = list(letters[1:3]))
  expect_error(compute_distances(y),
               "Assertion on 'x' failed: Must be of type 'numeric'")
  y <- matrix(seq_len(9), nrow = 3, ncol = 3)
  rownames(y) <- letters[seq(nrow(y))]
  expect_error(compute_distances(y, method = "str"),
               "Assertion on 'method' failed: Must be element of set")
  expect_error(compute_distances(y, use = "str"),
               "Assertion on 'use' failed: Must be element of set")
  expect_error(compute_distances(y, stand = "str"),
               "Assertion on 'stand' failed: Must be of type 'logical flag', not 'character'.")
  expect_error(compute_distances(y, dummy = "str"),
               "Assertion on 'dummy' failed: Must be of type 'number', not 'character'.")
  
  # all rows with variance
  x <- matrix(seq_len(9), nrow = 3, ncol = 3)
  rownames(x) <- letters[seq(nrow(x))]
  expect_equal(as.numeric(compute_distances(x)), rep(0, 3))
  expect_equal(as.numeric(compute_distances(x, stand = TRUE)), rep(1, 3))
  diag(x) <- NA
  expect_equal(as.numeric(compute_distances(x)), rep(0, 3))
  expect_equal(as.numeric(compute_distances(x, dummy = NA)), as.numeric(rep(NA, 3)))
  expect_equal(as.numeric(compute_distances(x, dummy = -Inf)), rep(-Inf, 3))
  
  xx <- matrix(-4:4, nrow = 3, ncol = 3)
  rownames(xx) <- letters[seq(nrow(xx))]
  expect_equal(as.numeric(compute_distances(xx)), rep(0, 3)) 
  expect_equal(as.numeric(compute_distances(xx, stand = TRUE)), rep(1, 3)) 
  expect_equal(as.numeric(compute_distances(xx, method = "manhattan")), c(3, 6, 3)) 
  
  # single row with no variance
  a <- t(matrix(c(rep(1, 3), sample(seq_len(100), 15)), nrow = 3, ncol = 6))
  rownames(a) <- letters[seq(nrow(a))]
  expect_equal(sort(unique(as.numeric(compute_distances(a)))), seq(0, 2, 0.5))
  
  # single row with variance
  b <- t(matrix(c(rep(1, 15), sample(seq_len(100), 3)), nrow = 3, ncol = 6))
  rownames(b) <- letters[seq(nrow(b))]
  expect_equal(as.numeric(compute_distances(b)), rep(1, 15))
  
  # single row with variance - small
  y <- matrix(c(1, 1, 2, 1), nrow = 2, ncol = 2)
  rownames(y) <- letters[seq(nrow(y))]
  expect_equal(as.numeric(compute_distances(y)), 1)
  
  # all rows with variance
  z <- matrix(sample(seq_len(100), 9), nrow = 3, ncol = 3)
  rownames(z) <- letters[seq(nrow(z))]
  expect_equal(as.numeric(compute_distances(z)), c(1.5, 1.5, 0))
  
  # all rows with no variance
  y <- matrix(1, nrow = 3, ncol = 3)
  rownames(y) <- letters[seq(nrow(y))]
  expect_equal(as.numeric(compute_distances(y)), rep(1, 3))
  
})


test_that("heatmaply force scale limits works correctly", {
  heatmaply <- heatmaply::heatmaply(
    mtcars,
    k_row = 3,
    k_col = 2,
    plot_method = "plotly",
    show_dendrogram = rep(all(dim(mtcars) >= 2), 2),
    colors = create_color_palette(c("blue", "green", "red"), range(mtcars), 5)
  )
  heatmaply$x$data_index <- 4
  
  expect_equal(heatmaply$x$data[[4]]$zmin, NULL)
  expect_equal(heatmaply$x$data[[4]]$zmax, NULL)
  
  limited_range <- c(100, 300)
  hm <- force_heatmaply_limits(heatmaply, limited_range)
  
  expect_equal(hm$x$data[[4]]$zmin, limited_range[1])
  expect_equal(hm$x$data[[4]]$zmax, limited_range[2])
  
  expect_error(force_heatmaply_limits(1, limited_range))
  expect_error(force_heatmaply_limits(heatmaply, 100))
})

test_that("create color palette works correctly", {
  test_colors <- c("blue", "green", "red")
  test_color_hex <- c("#0000FF", "#00FF00", "#FF0000")
  
  expect_equal(
    create_color_palette(colors_vec = test_colors, limits = c(-0.05, 0.1)), test_color_hex)
  
  pos_palette <- create_color_palette(colors_vec = test_colors,
                                      limits = c(1, 1.5),
                                      breaks = 0.1)
  expect_equal(pos_palette, c("#00FF00", "#3FBF00", "#7F7F00", "#BF3F00", "#FF0000"))
  expect_equal(pos_palette[1], test_color_hex[2])
  
  pos_0_palette <- create_color_palette(colors_vec = test_colors,
                                        limits = c(0, 0.5),
                                        breaks = 0.1)
  expect_equal(pos_0_palette, c("#00FF00", "#3FBF00", "#7F7F00", "#BF3F00", "#FF0000"))
  expect_equal(pos_0_palette[1], test_color_hex[2])
  
  neg_palette <- create_color_palette(colors_vec = test_colors,
                                      limits = c(-3, -1),
                                      breaks = 0.5)
  expect_equal(neg_palette, c("#0000FF", "#0055AA", "#00AA55", "#00FF00"))
  expect_equal(neg_palette[NROW(neg_palette)], test_color_hex[2])
  
  neg_0_palette <- create_color_palette(colors_vec = test_colors,
                                        limits = c(-2, 0),
                                        breaks = 0.5)
  expect_equal(neg_0_palette, c("#0000FF", "#0055AA", "#00AA55", "#00FF00"))
  expect_equal(neg_0_palette[NROW(neg_0_palette)], test_color_hex[2])
  
  full_palette <- create_color_palette(colors_vec = test_colors,
                                       limits = c(-2, 3),
                                       breaks = 0.75)
  expect_equal(
    full_palette,
    c("#0000FF", "#0055AA", "#00AA55", "#00FF00", "#55AA00", "#AA5500", "#FF0000"))
  expect_equal(full_palette[ceiling(length(full_palette) / 2)], test_color_hex[2])
  
  full_palette_x <- create_color_palette(colors_vec = test_colors,
                                         limits = c(3, -2),
                                         breaks = 0.75)
  expect_identical(full_palette_x, full_palette)
  
  full_palette_y <- create_color_palette(colors_vec = test_colors,
                                         limits = c(-2, 3),
                                         breaks = -0.75)
  expect_identical(full_palette_y, full_palette)
  
  full_palette_z <- create_color_palette(colors_vec = test_colors,
                                         limits = c(-2, 3),
                                         breaks = -0.75)
  expect_identical(full_palette_z, full_palette)
  
  zero_palette <- create_color_palette(colors_vec = test_colors,
                                       limits = c(0, 0),
                                       breaks = -0.75)
  expect_identical(zero_palette, test_color_hex[2])
  
  expect_error(create_color_palette(1, c(0, 1)), 
                         "Assertion on 'colors_vec' failed")
  expect_error(create_color_palette("blue", c(0, 1)), 
                         "Assertion on 'colors_vec' failed")
  expect_error(create_color_palette(c("red", "blue"), 1), 
                         "Assertion on 'limits' failed")
  expect_error(create_color_palette(c("red", "blue"), c(0, 1), "str"), 
                         "Assertion on 'breaks' failed: Must be of type 'number', not 'character'.")
  expect_error(create_color_palette(c("red", "blue"), c(0, 1), c(0, 1)), 
                         "Assertion on 'breaks' failed: Must have length 1.")
})


test_that("get_visualization_range works fine", {
  json_path <- system.file(package = "gDRplots", "settings.json")
  s <- gDRutils::get_settings_from_json(json_path = json_path)
  expect_identical(get_visualization_range(), s$VIS_RANGE)
})

