test_that("scoring functions", {
  expect_equal(sqi_more(0, 10)$fn(c(-1, 5, 20)), c(0, 0.5, 1))
  expect_equal(sqi_less(1, 2)$fn(c(0.5, 1.5, 3)), c(1, 0.5, 0))
  opt <- sqi_optimum(4, 5, 7, 8)
  expect_equal(opt$fn(c(3, 4.5, 6, 7.5, 9)), c(0, 0.5, 1, 0.5, 0))
  expect_error(sqi_more(10, 0), "increasing")
  expect_error(sqi_optimum(4, 6, 5, 8), "increasing")
  expect_snapshot(print(opt))
})

test_that("soil_quality_index() combines scores", {
  df <- data.frame(oc = c(0, 10, 20), bd = c(1, 1.5, 2))
  res <- soil_quality_index(df, oc = sqi_more(0, 20), bd = sqi_less(1, 2))
  expect_equal(res$oc_score, c(0, 0.5, 1))
  expect_equal(res$sqi, c(0.5, 0.5, 0.5))
  w <- soil_quality_index(df, oc = sqi_more(0, 20), bd = sqi_less(1, 2),
    weights = c(oc = 3, bd = 1))
  expect_equal(w$sqi, c(0.25, 0.5, 0.75))
  df$oc[2] <- NA
  expect_true(is.na(soil_quality_index(df, oc = sqi_more(0, 20), bd = sqi_less(1, 2))$sqi[2]))
  expect_equal(soil_quality_index(df, oc = sqi_more(0, 20), bd = sqi_less(1, 2),
    na_rm = TRUE)$sqi[2], 0.5)
  expect_error(soil_quality_index(df, sqi_more(0, 1)), "named")
  expect_error(soil_quality_index(df, oc = 1), "scoring")
  expect_error(soil_quality_index(df, zz = sqi_more(0, 1)), "not found")
  expect_error(soil_quality_index(df, oc = sqi_more(0, 20), weights = c(x = 1)), "weights")
})

test_that("land uses differ in the built-in trial", {
  res <- soil_quality_index(soil_quality_trial, oc = sqi_more(10, 40),
    bulk_density = sqi_less(1, 1.6), agg_stability = sqi_more(20, 80))
  by_lu <- tapply(res$sqi, res$land_use, mean)
  expect_gt(by_lu[["forest"]], by_lu[["food garden"]])
})

test_that("soil_mds() selects a minimum data set", {
  mds <- soil_mds(soil_quality_trial, oc:agg_stability)
  expect_named(mds, c("component", "eigenvalue", "variance", "indicator", "loading",
    "selected", "weight"))
  expect_true(all(mds$eigenvalue >= 1))
  expect_equal(sum(mds$weight, na.rm = TRUE), 1)
  expect_true(any(mds$selected))
  expect_error(soil_mds(soil_quality_trial, oc), "at least two")
  expect_error(soil_mds(soil_quality_trial, c(oc, land_use)), "numeric")
})
