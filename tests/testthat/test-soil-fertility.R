test_that("base saturation, ESP and ECEC", {
  df <- data.frame(ca = 8, mg = 3, k = 0.6, na = 0.4, cec = 20, acid = 1)
  res <- soil_base_saturation(df, ca, mg, k, na, cec, acidity = acid)
  expect_equal(res$sum_bases, 12)
  expect_equal(res$base_saturation, 60)
  expect_equal(res$esp, 2)
  expect_equal(res$ecec, 13)
  expect_equal(res$acid_saturation, 100 / 13)
  expect_equal(res$ca_mg_ratio, 8 / 3)
  high <- data.frame(ca = 30, mg = 3, k = 0.6, na = 0.4, cec = 20)
  expect_message(soil_base_saturation(high, ca, mg, k, na, cec), "exceed 100")
})

test_that("pH reaction classes", {
  df <- data.frame(ph = c(3, 4, 4.8, 5.3, 5.8, 6.3, 7, 7.6, 8.2, 8.7, 9.5))
  res <- soil_ph_class(df, ph)
  expect_s3_class(res$ph_class, "ordered")
  expect_equal(as.integer(res$ph_class), 1:11)
  expect_equal(as.character(soil_ph_class(data.frame(ph = 7.4), ph)$ph_class),
    "slightly alkaline")
  expect_error(soil_ph_class(data.frame(ph = 15), ph), "between")
})

test_that("fertility ratings use the threshold table", {
  df <- data.frame(cec = c(3, 8, 15, 30, 50), exch_k = c(0.1, 0.25, 0.5, 1, 3))
  res <- soil_fertility_rating(df, cec, k = exch_k)
  expect_equal(as.character(res$cec_rating), c("very low", "low", "medium", "high", "very high"))
  expect_equal(as.character(res$k_rating), c("very low", "low", "medium", "high", "very high"))
  expect_error(soil_fertility_rating(df, phosphorus = cec), "No thresholds")
  expect_error(soil_fertility_rating(df), "at least one")
  expect_error(soil_fertility_rating(df, cec, thresholds = data.frame(property = "cec")), "missing")
  thr <- soil_fertility_thresholds()
  expect_true(all(thr$low < thr$medium & thr$medium < thr$high & thr$high < thr$very_high))
})

test_that("Parker nutrient index", {
  df <- data.frame(g = rep(c("a", "b"), each = 4),
    r = c("low", "low", "medium", "high", "high", "very high", "high", "medium"))
  res <- soil_nutrient_index(df, r, .by = g)
  expect_equal(res$nutrient_index, c(7 / 4, 11 / 4))
  expect_equal(as.character(res$status), c("medium", "high"))
  all_rows <- soil_nutrient_index(df, r)
  expect_equal(nrow(all_rows), 1)
  expect_error(soil_nutrient_index(data.frame(r = "great"), r), "ratings")
})
