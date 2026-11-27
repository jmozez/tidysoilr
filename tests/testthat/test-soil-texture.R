test_that("USDA classes match reference points", {
  df <- data.frame(
    sand = c(92, 82, 65, 40, 20, 5, 60, 30, 10, 50, 5, 20),
    silt = c(4, 12, 25, 40, 65, 90, 12, 35, 55, 5, 50, 20),
    clay = c(4, 6, 10, 20, 15, 5, 28, 35, 35, 45, 45, 60)
  )
  res <- soil_texture_class(df, sand, silt, clay)
  expect_s3_class(res$texture, "factor")
  expect_equal(as.character(res$texture), c(
    "sand", "loamy sand", "sandy loam", "loam", "silt loam", "silt",
    "sandy clay loam", "clay loam", "silty clay loam", "sandy clay", "silty clay", "clay"
  ))
})

test_that("every point of the triangle gets exactly one class", {
  g <- expand.grid(sand = seq(0, 100, 1), clay = seq(0, 100, 1))
  g <- g[g$sand + g$clay <= 100, ]
  g$silt <- 100 - g$sand - g$clay
  usda <- soil_texture_class(g, sand, silt, clay)$texture
  fao <- soil_texture_class(g, sand, silt, clay, system = "fao")$texture
  expect_false(anyNA(usda))
  expect_false(anyNA(fao))
  expect_equal(nlevels(usda), 12)
  expect_equal(nlevels(fao), 5)
})

test_that("FAO classes follow the HWSD definitions", {
  df <- data.frame(sand = c(80, 40, 10, 30, 10), silt = c(10, 40, 70, 25, 20),
    clay = c(10, 20, 20, 45, 70))
  res <- soil_texture_class(df, sand, silt, clay, system = "fao")
  expect_equal(as.character(res$texture),
    c("coarse", "medium", "medium fine", "fine", "very fine"))
})

test_that("texture inputs are checked and normalised", {
  df <- data.frame(sand = 0.4, silt = 0.4, clay = 0.2)
  expect_error(soil_texture_class(df, sand, silt, clay), "fractions")
  off <- data.frame(sand = 50, silt = 30, clay = 30)
  expect_warning(res <- soil_texture_class(off, sand, silt, clay), "rescaled")
  expect_s3_class(res$texture, "factor")
  neg <- data.frame(sand = -1, silt = 60, clay = 41)
  expect_error(soil_texture_class(neg, sand, silt, clay), "negative")
  expect_error(soil_texture_class(off, sand, silt, clay, system = "uk"))
})

test_that("polygons are consistent with the classification", {
  for (sys in c("usda", "fao")) {
    p <- soil_texture_polygons(sys)
    expect_named(p, c("class", "vertex", "sand", "silt", "clay"))
    expect_true(all(abs(p$sand + p$silt + p$clay - 100) < 1e-9))
    cent <- dplyr::summarise(dplyr::group_by(p, class), sand = mean(sand), clay = mean(clay))
    cent$silt <- 100 - cent$sand - cent$clay
    res <- soil_texture_class(cent, sand, silt, clay, system = sys)
    expect_equal(as.character(res$texture), as.character(cent$class))
  }
})

test_that("ilr coordinates round trip", {
  df <- data.frame(sand = c(20, 60, 33), silt = c(30, 25, 33), clay = c(50, 15, 34))
  z <- soil_texture_ilr(df, sand, silt, clay)
  expect_true(all(c("ilr_1", "ilr_2") %in% names(z)))
  back <- soil_texture_ilr_inv(z, ilr_1, ilr_2)
  expect_equal(back$sand_ilr, df$sand / rowSums(df) * 100, tolerance = 1e-8)
  expect_equal(back$clay_ilr, df$clay / rowSums(df) * 100, tolerance = 1e-8)
  zero <- soil_texture_ilr(data.frame(sand = 100, silt = 0, clay = 0), sand, silt, clay)
  expect_true(all(is.finite(c(zero$ilr_1, zero$ilr_2))))
  expect_error(soil_texture_ilr(df, sand, silt, clay, names = "a"), "length 2")
  expect_error(soil_texture_ilr_inv(z, ilr_1, ilr_2, names = "a"), "length 3")
})
