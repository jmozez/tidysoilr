test_that("soil_profile_check() finds gaps, overlaps and inversions", {
  bad <- data.frame(
    id = c("A", "A", "A", "B", "B", "C"),
    top = c(0, 10, 25, 0, 8, NA),
    bottom = c(10, 20, 22, 10, 30, 5)
  )
  res <- soil_profile_check(bad, id = id, top = top, bottom = bottom)
  expect_setequal(res$issue, c("missing depth", "inverted depths", "gap", "overlap"))
  expect_equal(res$size[res$issue == "overlap"], 2)
  ok <- data.frame(id = "A", top = c(0, 10), bottom = c(10, 20))
  expect_equal(nrow(soil_profile_check(ok, id = id, top = top, bottom = bottom)), 0)
})

test_that("soil_profile_check() reads metadata from a soil_profile", {
  prof <- soil_profile(png_profiles, profile, top, bottom)
  res <- soil_profile_check(prof)
  expect_s3_class(res, "tbl_df")
  expect_named(res, c("id", "top", "bottom", "issue", "size"))
  expect_error(soil_profile_check(png_profiles), "id")
})

test_that("soil_slice() slices into regular increments", {
  df <- data.frame(id = "A", top = c(0, 10), bottom = c(10, 25), clay = c(20, 40))
  prof <- soil_profile(df, id, top, bottom)
  s <- soil_slice(prof, clay)
  expect_s3_class(s, "soil_profile")
  expect_equal(nrow(s), 25)
  expect_equal(s$clay[s$top == 9], 20)
  expect_equal(s$clay[s$top == 10], 40)
  s5 <- soil_slice(prof, interval = 5, max_depth = 20)
  expect_equal(s5$top, c(0, 5, 10, 15))
  expect_error(soil_slice(prof, interval = -1), "positive")
})

test_that("soil_depth_mean() weights by thickness and prorates", {
  df <- data.frame(id = c("A", "A", "B"), top = c(0, 10, 0), bottom = c(10, 40, 20),
    clay = c(20, 50, 30))
  res <- soil_depth_mean(df, clay, from = 0, to = 30, id = id, top = top, bottom = bottom)
  expect_equal(res$clay[res$id == "A"], (20 * 10 + 50 * 20) / 30)
  expect_equal(res$clay[res$id == "B"], 30)
  expect_equal(res$coverage, c(1, 2 / 3))
  expect_error(soil_depth_mean(df, clay, from = 30, to = 10, id = id, top = top, bottom = bottom),
    "greater")
  expect_error(soil_depth_mean(df, id = id, top = top, bottom = bottom), "at least one")
  expect_error(soil_depth_mean(df, id, id = id, top = top, bottom = bottom), "numeric")
})
