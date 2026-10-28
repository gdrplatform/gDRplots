context("Test common fun colors")

test_that("is_valid_color works as expected", {
  color_names <- c("#33cc33", "#d6f5d6", "#0a290a", "#F9B42DFF", "#714D6932", "#C2F970DC")
  expect_true(all(vapply(color_names, is_valid_color, logical(1))))

  color_names <- c("darkblue", "yellow", "tomato")
  expect_true(all(vapply(color_names, is_valid_color, logical(1))))

  color_names <- c("nice pink", "RED", "#0a290", "#C2F970D")
  expect_false(all(vapply(color_names, is_valid_color, logical(1))))

  expect_error(is_valid_color(1), "Must be of type 'string'")
  expect_error(is_valid_color(NULL), "Must be of type 'string'")
  expect_error(is_valid_color(NA), "Assertion on 'col_name' failed: May not be NA.")
})

test_that("get_iso_colors works as expected", {
  ### expected values
  gic <- get_iso_colors()
  expect_true(length(gic) > 2)
  expect_identical(class(gic), "character")
  gic2 <- get_iso_colors(formals(get_iso_colors)[[1]][[3]])
  expect_true(any(gic != gic2))
  expect_identical(length(gic), length(gic2))

  ### errors
  expect_error(get_iso_colors("inv_param"), "'arg' should be one of ")
})

test_that("assert_RGB_format as expected", {
  color_vector <- c(25, 56, 189)
  expect_equal(assert_RGB_format(color_vector), NULL)
  color_vector <- c(201, 128, 352)
  expect_error(assert_RGB_format(color_vector),
               "Some value is greater than 255. Not valid RGB format.")
})

test_that("get_qual_colors works as expected", {
  max_len <- sum(RColorBrewer::brewer.pal.info[
    RColorBrewer::brewer.pal.info$category == "qual" &
      RColorBrewer::brewer.pal.info$colorblind == TRUE, ]$maxcolors)

  expect_equal(NROW(get_qual_colors()), max_len)
  expect_equal(NROW(unique(get_qual_colors())), max_len)
  expect_equal(get_qual_colors(0), "#000000")
  expect_equal(get_qual_colors(5), c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E"))
  expect_equal(NROW(get_qual_colors(42)), 42)
  expect_equal(NROW(unique(get_qual_colors(42))), 42)

  expect_error(get_qual_colors("one"),
               "Assertion on 'n' failed: Must be of type 'single integerish value'")
  expect_error(get_qual_colors(c(2, 3)),
               "Assertion on 'n' failed: Must have length 1.")
  expect_error(get_qual_colors(2.3),
               "Assertion on 'n' failed: Must be of type 'single integerish value'")
  expect_error(get_qual_colors(-1),
               "Assertion on 'n' failed: Element 1 is not >= 0.")
})


test_that("is_color_dark works as expected", {
  color_names <- c("#33cc33", "#d6f5d6", "#0a290a")
  expect_false(is_color_dark(color_names[1]))
  expect_false(is_color_dark(color_names[2]))
  expect_true(is_color_dark(color_names[3]))

  color_names <- c("darkblue", "yellow", "tomato")
  expect_true(is_color_dark(color_names[1]))
  expect_false(is_color_dark(color_names[2]))
  expect_false(is_color_dark(color_names[3]))

  expect_error(is_color_dark(1), "Must be of type 'string'")
  expect_error(is_color_dark("nice pink"), "Must be a valid color name")
})


test_that("get_col_luminance works as expected", {
  color_names <- c("#33cc33", "#d6f5d6", "#0a290a")
  expect_equal(get_col_luminance(color_names[1]), 0.429871, tolerance = 1e-5)
  expect_equal(get_col_luminance(color_names[2]), 0.839746, tolerance = 1e-5)
  expect_equal(get_col_luminance(color_names[3]), 0.016340, tolerance = 1e-5)

  color_names <- c("darkblue", "yellow", "tomato")
  expect_equal(get_col_luminance(color_names[1]), 0.0186408, tolerance = 1e-5)
  expect_equal(get_col_luminance(color_names[2]), 0.9278000, tolerance = 1e-5)
  expect_equal(get_col_luminance(color_names[3]), 0.3238907, tolerance = 1e-5)

  expect_error(get_col_luminance(1), "Must be of type 'string'")
  expect_error(get_col_luminance("nice pink"), "Must be a valid color name")
})

test_that(".resolve_group_colors works as expected", {
  req_grp_name <- c("cell_DC", "cell_XY", "cell_AB")
  data_order <- c("cell_AB", "cell_DC", "cell_XY")

  res_1 <- .resolve_group_colors(
    group_names_requested = req_grp_name,
    group_names = data_order)
  expect_equal(names(res_1), data_order)
  expect_true(all(vapply(res_1, is_valid_color, logical(1))))

  col_2 <- c("blue", "red", "grey")
  res_2 <- .resolve_group_colors(
    group_names_requested = req_grp_name,
    group_names = data_order,
    colors_vec = col_2)

  expect_equal(res_2, setNames(col_2, req_grp_name)[data_order])
  expect_equal(names(res_2), data_order)
  expect_true(all(vapply(res_2, is_valid_color, logical(1))))

  col_3 <- c("blue", "pinky", "grey")
  res_3 <- .resolve_group_colors(
    group_names_requested = req_grp_name,
    group_names = data_order,
    colors_vec = col_3)

  expect_equal(res_3[1:2], res_2[1:2])
  expect_true(is_valid_color(res_3[3]))

  col_4 <- c(col_2, "orange")
  res_4 <- .resolve_group_colors(
    group_names_requested = req_grp_name,
    group_names = data_order,
    colors_vec = col_4)

  expect_equal(res_4, res_2)

  res_5 <- .resolve_group_colors(
    group_names_requested = req_grp_name,
    group_names = data_order[1:2],
    colors_vec = col_2[1:2])

  col_vec_5 <- setNames(col_2[1:2], req_grp_name[1:2])
  expect_true(col_vec_5[names(col_vec_5) %in% data_order[1:2]] %in% res_5)

  col_6 <- c("blue", "lightblue")
  res_6 <- .resolve_group_colors(
    group_names_requested = req_grp_name,
    group_names = data_order,
    colors_vec = col_6)

  col_vec_6 <- setNames(col_6, req_grp_name[1:2])
  expect_equal(res_6[names(res_6) %in% names(col_vec_6)], col_vec_6)

  col_7 <- c("pinky", "blackish", "dark")
  res_7 <- .resolve_group_colors(
    group_names_requested = req_grp_name,
    group_names = data_order,
    colors_vec = col_7)

  expect_equal(names(res_7), data_order)
  expect_true(all(vapply(res_7, is_valid_color, logical(1))))
  expect_false(all(col_7 %in% res_7))

  req_grp_name_8 <- sprintf("%s_123", req_grp_name)
  res_8 <- .resolve_group_colors(
    group_names_requested = req_grp_name_8,
    group_names = data_order,
    colors_vec = col_2)

  expect_equal(names(res_8), data_order)
  expect_true(all(vapply(res_8, is_valid_color, logical(1))))
  expect_false(all(req_grp_name_8 %in% names(res_8)))
  expect_false(all(col_2 %in% res_8))

  col_9 <- c(col_6, NA)
  res_9 <- .resolve_group_colors(
    group_names_requested = req_grp_name,
    group_names = data_order,
    colors_vec = col_9)

  expect_equal(names(res_9), data_order)
  expect_equal(res_9[names(res_9) %in% names(col_vec_6)], col_vec_6)
  expect_true(is_valid_color(res_9[req_grp_name[which(is.na(col_9))]]))

  col_10 <- col_2
  names(col_10) <- data_order
  res_10 <- .resolve_group_colors(
    group_names_requested = req_grp_name,
    group_names = data_order,
    colors_vec = col_10)

  expect_equal(names(res_10), data_order)
  expect_equal(res_10, col_10)

  col_11 <- col_2
  names(col_11) <- req_grp_name
  res_11 <- .resolve_group_colors(
    group_names_requested = req_grp_name,
    group_names = data_order,
    colors_vec = col_11)

  expect_equal(names(res_11), data_order)
  expect_equal(res_11, col_11[data_order])

  col_12 <- c(col_10[1:2], "orange")
  expect_warning({
    res_12 <- .resolve_group_colors(
      group_names_requested = req_grp_name,
      group_names = data_order,
      colors_vec = col_12)
  }, "'colors_vec' has some, but not all, elements named")

  expect_equal(names(res_12), data_order)
  expect_equal(res_12, setNames(col_12, req_grp_name)[data_order])

  # named colors_vec ignores group_names_requested entirely - names alone decide the mapping
  col_13 <- col_2
  names(col_13) <- data_order
  res_13 <- .resolve_group_colors(
    group_names_requested = req_grp_name,
    group_names = data_order,
    colors_vec = col_13)
  res_13_no_request <- .resolve_group_colors(
    group_names_requested = NULL,
    group_names = data_order,
    colors_vec = col_13)

  expect_equal(res_13, col_13[data_order])
  expect_equal(res_13, res_13_no_request)

  col_14 <- col_13
  names(col_14)[3] <- NA
  expect_warning({
    res_14 <- .resolve_group_colors(
      group_names_requested = req_grp_name,
      group_names = data_order,
      colors_vec = col_14)
  }, "'colors_vec' has some, but not all, elements named")

  expect_equal(res_14, setNames(col_2, req_grp_name)[data_order])

  col_15 <- col_13
  names(col_15)[3] <- ""
  expect_warning({
    res_15 <- .resolve_group_colors(
      group_names_requested = req_grp_name,
      group_names = data_order,
      colors_vec = col_15)
  }, "'colors_vec' has some, but not all, elements named")

  expect_equal(res_15, setNames(col_2, req_grp_name)[data_order])

  # scenario: fully named `colors_vec` covering only part of group_names
  col_16 <- c(cell_AB = "grey")
  res_16 <- .resolve_group_colors(
    group_names_requested = req_grp_name,
    group_names = data_order,
    colors_vec = col_16)

  expect_equal(names(res_16), data_order)
  expect_equal(res_16[names(col_16)], col_16)
  expect_true(all(
    vapply(res_16[!names(res_16) %in% names(col_16)], is_valid_color, logical(1))))

  # scenario: named `colors_vec` with a name not present in group_names
  col_17 <- col_2
  names(col_17) <- req_grp_name
  names(col_17)[3] <- "cell_EF"
  res_17 <- .resolve_group_colors(
    group_names_requested = req_grp_name,
    group_names = data_order,
    colors_vec = col_17)

  def_col <- intersect(req_grp_name, names(col_17))
  expect_equal(names(res_17), data_order)
  expect_equal(res_17[def_col], col_17[def_col])
  expect_true(is_valid_color(res_17[setdiff(names(res_17), def_col)]))
  expect_false("grey" %in% res_17)

  # scenario: group_names / group_names_requested arrive as factors
  res_18 <- .resolve_group_colors(
    group_names_requested = factor(req_grp_name),
    group_names = data_order,
    colors_vec = col_2)
  expect_equal(res_18, res_2)

  res_19 <- .resolve_group_colors(
    group_names_requested = req_grp_name,
    group_names = factor(data_order),
    colors_vec = col_2)
  expect_equal(res_19, res_2)

  res_20 <- .resolve_group_colors(
    group_names_requested = factor(req_grp_name),
    group_names = factor(data_order),
    colors_vec = col_2)
  expect_equal(res_20, res_2)

  res_21 <- .resolve_group_colors(
    group_names_requested = NULL,
    group_names = factor(data_order),
    colors_vec = col_2)
  expect_equal(names(res_21), data_order)
  expect_true(all(vapply(res_21, is_valid_color, logical(1))))

  expect_error(.resolve_group_colors(group_names_requested = 1:3,
                                     group_names = data_order),
               "Assertion on 'group_names_requested' failed: Must be of type 'character'")
  expect_error(.resolve_group_colors(group_names_requested = rep("A", 3),
                                     group_names = data_order),
               "Assertion on 'group_names_requested' failed: Contains duplicated values")
  expect_error(.resolve_group_colors(group_names_requested = req_grp_name,
                                     group_names = 1:3),
               "Assertion on 'group_names' failed: Must be of type 'character'")
  expect_error(.resolve_group_colors(group_names_requested = req_grp_name,
                                     group_names = rep("A", 3)),
               "Assertion on 'group_names' failed: Contains duplicated values")
  expect_error(.resolve_group_colors(group_names_requested = req_grp_name,
                                     group_names = NULL),
               "Assertion on 'group_names' failed: Must be of type 'character'")
  expect_error(.resolve_group_colors(group_names_requested = req_grp_name,
                                     group_names = data_order,
                                     colors_vec = 1:5),
               "Assertion on 'colors_vec' failed: Must be of type 'character'")

})
