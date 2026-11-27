test_that("soil_profile_cv() keeps profiles together", {
  set.seed(1)
  prof <- soil_profile(png_profiles, profile, top, bottom)
  folds <- soil_profile_cv(prof, v = 4)
  expect_s3_class(folds, c("soil_profile_cv", "rset"))
  expect_equal(nrow(folds), 4)
  for (s in folds$splits) {
    a <- rsample::analysis(s)
    b <- rsample::assessment(s)
    expect_length(intersect(a$profile, b$profile), 0)
    expect_equal(nrow(a) + nrow(b), nrow(png_profiles))
  }
  expect_output(print(folds), "Profile-grouped")
  expect_error(soil_profile_cv(png_profiles), "id")
  expect_error(soil_profile_cv(prof, v = 100), "profiles")
  ok <- soil_profile_cv(png_profiles, id = profile, v = 3)
  expect_equal(nrow(ok), 3)
})

test_that("soil_spatial_block_cv() makes spatially separate folds", {
  set.seed(2)
  folds <- soil_spatial_block_cv(soil_nir, lon = lon, lat = lat, v = 4)
  expect_s3_class(folds, c("soil_spatial_block_cv", "rset"))
  ids <- lapply(folds$splits, function(s) as.integer(s, data = "assessment"))
  expect_setequal(unlist(ids), seq_len(nrow(soil_nir)))
  expect_output(print(folds), "Spatial block")
  expect_error(soil_spatial_block_cv(soil_nir, lon = lon, lat = lat, v = 4, block_size = 100),
    "block")
})

test_that("soil_spatial_block_cv() accepts sf points", {
  skip_if_not_installed("sf")
  pts <- sf::st_as_sf(soil_nir[c("sample", "lon", "lat", "clay")], coords = c("lon", "lat"),
    crs = 4326)
  set.seed(3)
  folds <- soil_spatial_block_cv(pts, v = 3)
  expect_equal(nrow(folds), 3)
})

test_that("soil_kennard_stone() spans the predictor space", {
  spec <- soil_spectra(soil_nir, dplyr::starts_with("nm_"))
  split <- soil_kennard_stone(spec, prop = 0.7)
  expect_s3_class(split, "rsplit")
  expect_equal(nrow(rsample::training(split)), 84)
  expect_equal(nrow(rsample::testing(split)), 36)
  df <- data.frame(x = c(5, 2, 9, 1, 7, 10, 3, 6, 4, 8))
  ks <- soil_kennard_stone(df, x, prop = 0.3)
  expect_setequal(df$x[as.integer(ks, "analysis")], c(1, 10, 5))
  expect_error(soil_kennard_stone(df, x, prop = 0.01), "adjust")
  expect_error(soil_kennard_stone(df), "Select")
})
