test_that("sep() is the bias-corrected standard deviation of errors", {
  df <- data.frame(obs = c(10, 20, 30, 40, 50), pred = c(13, 22, 31, 45, 52))
  e <- df$pred - df$obs
  expect_equal(sep_vec(df$obs, df$pred), sd(e))
  res <- sep(df, obs, pred)
  expect_equal(res$.metric, "sep")
  expect_equal(res$.estimate, sd(e))
  shifted <- data.frame(obs = df$obs, pred = df$obs + 5)
  expect_equal(sep_vec(shifted$obs, shifted$pred), 0)
  expect_equal(sep_vec(c(1, NA, 3), c(2, 2, 5)), sd(c(1, 2)))
  expect_true(is.na(sep_vec(c(1, NA, 3), c(2, 2, 5), na_rm = FALSE)))
  expect_equal(sep_vec(df$obs, df$pred, case_weights = rep(2, 5)), sd(e))
  expect_equal(attr(sep, "direction"), "minimize")
})

test_that("sep() works in a metric set with yardstick metrics", {
  df <- data.frame(obs = c(10, 20, 30, 40, 50), pred = c(13, 22, 31, 45, 52))
  ms <- yardstick::metric_set(yardstick::rmse, sep, yardstick::rpiq, yardstick::ccc)
  res <- ms(df, obs, pred)
  expect_equal(res$.metric, c("rmse", "sep", "rpiq", "ccc"))
  grouped <- dplyr::group_by(data.frame(df, g = c(1, 1, 1, 2, 2)), g)
  expect_equal(nrow(sep(grouped, obs, pred)), 2)
})

test_that("texture_adjacent_accuracy() credits neighbouring classes", {
  lv <- levels(soil_texture_polygons("usda")$class)
  df <- data.frame(
    truth = factor(c("clay", "silty clay loam", "loam", "sand"), levels = lv),
    pred = factor(c("clay", "clay loam", "sandy loam", "silt"), levels = lv)
  )
  res <- texture_adjacent_accuracy(df, truth, pred)
  expect_equal(res$.estimate, 0.75)
  expect_equal(res$.estimator, "multiclass")
  expect_equal(texture_adjacent_accuracy_vec(df$truth, df$pred, case_weights = c(1, 1, 1, 5)),
    3 / 8)
  fl <- levels(soil_texture_polygons("fao")$class)
  f <- data.frame(t = factor(c("coarse", "fine"), levels = fl),
    p = factor(c("medium", "coarse"), levels = fl))
  expect_equal(texture_adjacent_accuracy(f, t, p, system = "fao")$.estimate, 0.5)
  bad <- data.frame(t = factor(c("a", "b")), p = factor(c("a", "b")))
  expect_error(texture_adjacent_accuracy(bad, t, p), "texture class")
  expect_equal(attr(texture_adjacent_accuracy, "direction"), "maximize")
})

test_that("adjacency is symmetric and excludes point contacts", {
  adj <- tidysoilr:::get_adjacency("usda")
  expect_true("clay loam|silty clay loam" %in% adj)
  expect_false("clay|silty clay loam" %in% adj)
  expect_false("sand|silt" %in% adj)
})
