test_that("bulk density and fine-earth bulk density", {
  df <- data.frame(m = c(130, 142), v = 100, cm = c(0, 26.5))
  res <- soil_bulk_density(df, m, v, coarse_mass = cm)
  expect_equal(res$bulk_density, c(1.30, 1.42))
  expect_equal(res$coarse_vol, c(0, 0.1))
  expect_equal(res$fine_earth_bd[2], (142 - 26.5) / 90)
  expect_false("fine_earth_bd" %in% names(soil_bulk_density(df, m, v)))
  expect_error(soil_bulk_density(data.frame(m = 1, v = 1, cm = 2), m, v, cm), "exceed")
})

test_that("porosity, air-filled porosity and WFPS", {
  df <- data.frame(bd = 1.325, w = 0.25)
  res <- soil_porosity(df, bd, water = w)
  expect_equal(res$porosity, 0.5)
  expect_equal(res$air_filled, 0.25)
  expect_equal(res$wfps, 0.5)
  expect_warning(soil_porosity(data.frame(bd = 3), bd), "negative")
})

test_that("Saxton and Rawls (2006) reproduces published values", {
  # Table 3 of Saxton and Rawls (2006), organic matter 2.5 %:
  # loamy sand 5/12/46 %, Ksat 96.7 mm/h; loam 12/27/46 % (approximately).
  df <- data.frame(sand = c(80, 40), clay = c(5, 20), om = 2.5)
  res <- soil_ptf_saxton(df, sand, clay, om)
  expect_equal(round(res$theta_1500 * 100), c(5, 14))
  expect_equal(round(res$theta_33 * 100), c(12, 28))
  expect_equal(round(res$theta_sat * 100), c(46, 46))
  expect_equal(round(res$ksat, 1)[1], 96.7)
  expect_equal(res$paw, res$theta_33 - res$theta_1500)
  expect_true(all(res$bulk_density > 1 & res$bulk_density < 1.8))
  expect_warning(soil_ptf_saxton(data.frame(sand = 40, clay = 20, om = 12), sand, clay, om),
    "outside")
})

test_that("bulk density PTFs", {
  df <- data.frame(oc = c(0, 1, 4))
  expect_equal(soil_ptf_bd(df, oc)$bd_ptf, 1.660 - 0.318 * sqrt(c(0, 1, 4)))
  expect_equal(soil_ptf_bd(df, oc, method = "alexander")$bd_ptf, 1.72 - 0.294 * sqrt(c(0, 1, 4)))
  expect_equal(soil_ptf_bd(df, oc, method = "honeysett_ratkowsky")$bd_ptf,
    1 / (0.548 + 0.0588 * c(0, 1, 4)))
  ad <- soil_ptf_bd(df, oc, method = "adams", .name = "bd")
  expect_equal(ad$bd[1], 1.64)
  expect_true(all(diff(ad$bd) < 0))
  expect_warning(soil_ptf_bd(data.frame(oc = 120), oc), "g/kg")
  expect_error(soil_ptf_bd(df, oc, method = "nope"))
})
