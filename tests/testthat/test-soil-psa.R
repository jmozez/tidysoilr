test_that("hydrometer calculation recovers known fractions", {
  df <- data.frame(mass = 50, r1 = 30 + 5, r2 = 10 + 4, b1 = 5, b2 = 4)
  res <- soil_hydrometer(df, mass = mass, reading_silt_clay = r1, reading_clay = r2,
    blank_silt_clay = b1, blank_clay = b2)
  expect_equal(res$sand, 40)
  expect_equal(res$clay, 20)
  expect_equal(res$silt, 40)
})

test_that("hydrometer temperature correction and messages", {
  df <- data.frame(mass = 50, r1 = 30, r2 = 10, t1 = 25, t2 = 25)
  res <- soil_hydrometer(df, mass, r1, r2, temp_silt_clay = t1, temp_clay = t2)
  expect_equal(res$clay, 100 * (10 + 0.36 * 5) / 50)
  expect_message(soil_hydrometer(df, mass, r1, r2), "uncorrected")
  swapped <- data.frame(mass = 50, r1 = 10, r2 = 30, t1 = 20, t2 = 20)
  expect_warning(soil_hydrometer(swapped, mass, r1, r2, temp_silt_clay = t1, temp_clay = t2),
    "exceed")
  expect_error(soil_hydrometer(df, mass, r1, r2, volume = 0), "positive")
})

test_that("the built-in hydrometer data give sensible fractions", {
  res <- soil_hydrometer(soil_hydrometer_lab, mass, r_40s, r_2h, blank_40s, blank_2h)
  expect_equal(res$sand + res$silt + res$clay, rep(100, nrow(res)))
  expect_true(all(res$clay >= 0 & res$silt >= 0))
})

test_that("pipette method computes fractions and recovery", {
  df <- data.frame(mass = 10, sand = 4, sc = 0.15 + 0.0125, cl = 0.05 + 0.0125)
  res <- soil_pipette(df, mass, sand, sc, cl, aliquot = 25, volume = 1000, dispersant = 0.0125)
  expect_equal(res$recovery, 100)
  expect_equal(res$sand, 40)
  expect_equal(res$clay, 20)
  expect_equal(res$silt, 40)
  raw <- soil_pipette(df, mass, sand, sc, cl, dispersant = 0.0125, normalise = FALSE)
  expect_equal(raw$clay, 20)
  expect_error(soil_pipette(df, mass, sand, sc, cl, aliquot = -1), "positive")
})
