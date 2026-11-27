test_that("equal-area spline preserves horizon means with small lambda", {
  df <- data.frame(id = "A", top = c(0, 10, 30, 60), bottom = c(10, 30, 60, 100),
    oc = c(30, 20, 10, 5))
  prof <- soil_profile(df, id, top, bottom)
  res <- soil_spline(prof, oc, depths = c(0, 10, 30, 60, 100), lambda = 1e-6)
  expect_s3_class(res, "soil_profile")
  expect_equal(res$oc, c(30, 20, 10, 5), tolerance = 1e-4)
  smooth <- soil_spline(prof, oc, depths = c(0, 10, 30, 60, 100))
  expect_equal(smooth$oc, c(30, 20, 10, 5), tolerance = 0.05)
})

test_that("spline handles gaps, single horizons and depth limits", {
  df <- data.frame(id = c("A", "A", "A", "B"), top = c(0, 15, 40, 0),
    bottom = c(10, 30, 60, 20), clay = c(30, 20, 10, 25))
  res <- soil_spline(df, clay, depths = c(0, 5, 15, 30, 100), id = id, top = top, bottom = bottom)
  a <- res[res$id == "A", ]
  expect_true(all(!is.na(a$clay[1:3])))
  expect_true(is.na(a$clay[4]))
  b <- res[res$id == "B", ]
  expect_equal(b$clay[1:2], c(25, 25))
  expect_true(all(is.na(b$clay[3:4])))
  capped <- soil_spline(df, clay, depths = c(0, 5), upper = 26, id = id, top = top, bottom = bottom)
  expect_true(all(capped$clay <= 26))
  expect_error(soil_spline(df, clay, depths = 5, id = id, top = top, bottom = bottom), "increasing")
  expect_error(soil_spline(df, id = id, top = top, bottom = bottom), "at least one")
})

test_that("spline runs on the PNG profiles", {
  prof <- soil_profile(png_profiles, profile, top, bottom)
  res <- soil_spline(prof, oc, clay)
  expect_equal(nrow(res), 32 * 6)
  expect_true(all(res$oc >= 0, na.rm = TRUE))
})

test_that("carbon stock per layer", {
  df <- data.frame(top = 0, bottom = 30, oc = 10, bd = 1.2, cf = 10)
  res <- soil_carbon_stock(df, oc, bd, top = top, bottom = bottom)
  expect_equal(res$soc_stock, 36)
  expect_equal(soil_carbon_stock(df, oc, bd, coarse = cf, top = top, bottom = bottom)$soc_stock,
    36 * 0.9)
  pct <- data.frame(top = 0, bottom = 30, oc = 1, bd = 1.2)
  expect_equal(soil_carbon_stock(pct, oc, bd, oc_units = "percent", top = top,
    bottom = bottom)$soc_stock, 36)
  expect_error(soil_carbon_stock(df, oc, bd, coarse = cf, coarse_units = "fraction",
    top = top, bottom = bottom), "exceed")
})

test_that("carbon totals prorate layers to a fixed depth", {
  df <- data.frame(id = c("A", "A"), top = c(0, 20), bottom = c(20, 40), s = c(40, 20))
  res <- soil_carbon_total(df, s, from = 0, to = 30, id = id, top = top, bottom = bottom)
  expect_equal(res$s, 40 + 10)
  expect_equal(res$coverage, 1)
  short <- soil_carbon_total(df, s, to = 60, id = id, top = top, bottom = bottom)
  expect_equal(short$coverage, 40 / 60)
  prof <- soil_profile(png_profiles, profile, top, bottom) |> soil_carbon_stock(oc, bulk_density)
  tot <- soil_carbon_total(prof, soc_stock)
  expect_equal(nrow(tot), 32)
  expect_error(soil_carbon_total(df, s, from = 10, to = 5, id = id, top = top, bottom = bottom),
    "greater")
})
