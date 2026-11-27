esm_layers <- function() {
  tibble::tibble(
    id = rep(c("a", "b"), each = 3),
    top = rep(c(0, 10, 30), 2),
    bottom = rep(c(10, 30, 60), 2),
    oc = c(30, 20, 10, 30, 20, 10),
    bd = c(1.0, 1.2, 1.4, 1.3, 1.4, 1.5)
  )
}

test_that("soil_carbon_esm() matches a hand calculation", {
  d <- esm_layers()
  res <- soil_carbon_esm(d, oc, bd, depths = 30, method = "linear", id = id, top = top,
    bottom = bottom)
  # Masses to 30 cm: a = 1000 + 2400 = 3400, b = 1300 + 2800 = 4100 Mg/ha
  expect_equal(res$ref_mass, c(3400, 3400))
  expect_equal(res$soc_esm[1], 30 + 48)
  # b: 39 Mg C in the first 1300 Mg, then 2100 Mg at 20 g/kg
  expect_equal(res$soc_esm[2], 39 + 2100 * 0.02)
  expect_equal(res$esm_depth[2], 10 + 2100 / 140)
  expect_equal(res$soc_fixed, c(78, 39 + 56))
  expect_false(any(res$extrapolated))
  # The equivalent mass removes the bulk density effect on stocks
  expect_lt(res$soc_esm[2], res$soc_fixed[2])
})

test_that("soil_carbon_esm() agrees with soil_carbon_stock() at layer limits", {
  d <- esm_layers()
  st <- soil_carbon_stock(d, oc, bd, top = top, bottom = bottom)
  for (method in c("linear", "cubic")) {
    res <- soil_carbon_esm(d, oc, bd, ref_mass = c(1000, 3400, 7600), method = method,
      id = id, top = top, bottom = bottom)
    a <- res[res$id == "a", ]
    expect_equal(a$soc_esm, cumsum(st$soc_stock[1:3]))
    expect_equal(a$esm_depth, c(10, 30, 60))
  }
})

test_that("soil_carbon_esm() handles references, extrapolation and bad profiles", {
  d <- esm_layers()
  res <- soil_carbon_esm(d, oc, bd, depths = 30, reference = "b", id = id, top = top,
    bottom = bottom, method = "linear")
  expect_equal(res$ref_mass[1], 4100)
  expect_equal(res$soc_esm[2], res$soc_fixed[2])
  big <- soil_carbon_esm(d, oc, bd, ref_mass = 8000, id = id, top = top, bottom = bottom)
  expect_equal(big$extrapolated, c(TRUE, FALSE))
  expect_equal(big$soc_esm[1], 30 + 48 + 42 + 400 * 0.01)
  d$bd[2] <- NA
  expect_warning(out <- soil_carbon_esm(d, oc, bd, id = id, top = top, bottom = bottom), "gaps")
  expect_true(is.na(out$soc_esm[1]))
  expect_error(soil_carbon_esm(d, oc, bd, ref_mass = -1, id = id, top = top, bottom = bottom))
  expect_error(soil_carbon_esm(esm_layers(), oc, bd, reference = "z", id = id, top = top,
    bottom = bottom))
})

test_that("soil_carbon_esm() runs on the real PNG profiles", {
  prof <- png_profiles[png_profiles$profile %in% c("82P0793", "82P0794", "82P0800"), ]
  prof$coarse_frag[is.na(prof$coarse_frag)] <- 0
  res <- soil_carbon_esm(prof, oc, bulk_density, coarse = coarse_frag, depths = c(30, 100),
    id = profile, top = top, bottom = bottom)
  expect_equal(nrow(res), 6L)
  expect_true(all(res$esm_depth <= res$depth + 1e-8))
  expect_true(all(res$soc_esm > 0))
})

test_that("soil_vg_conductivity() matches the closed form in h", {
  h <- c(0, 1, 10, 100, 330, 1000, 15000)
  ks <- 24.96
  alpha <- 0.036
  n <- 1.56
  m <- 1 - 1 / n
  ah <- alpha * h
  ref <- ks * (1 - ah^(n - 1) * (1 + ah^n)^(-m))^2 / (1 + ah^n)^(m / 2)
  out <- soil_vg_conductivity(data.frame(h = h), h, ks = ks, alpha = alpha, n = n)
  expect_equal(out$k_mvg, ref)
  expect_equal(out$k_mvg[1], ks)
  expect_true(all(diff(out$k_mvg) < 0))
  expect_error(soil_vg_conductivity(data.frame(h = 1), h, ks = 1, alpha = 0.1, n = 0.9), "n")
  expect_error(soil_vg_conductivity(data.frame(h = 1), h, ks = -1, alpha = 0.1, n = 1.5), "ks")
})

test_that("soil_infiltration_fit() recovers known parameters", {
  t <- c(1, 2, 5, 10, 20, 30, 45, 60, 90, 120)
  d <- rbind(
    data.frame(test = "philip", t = t, i = 2 * sqrt(t) + 0.3 * t),
    data.frame(test = "horton", t = t, i = 0.5 * t + (3 - 0.5) / 0.1 * (1 - exp(-0.1 * t))),
    data.frame(test = "kostiakov", t = t, i = 1.5 * t^0.7)
  )
  fit <- soil_infiltration_fit(d, t, i, .by = test)
  expect_equal(nrow(fit), 9L)
  p <- fit[fit$test == "philip" & fit$model == "philip", ]
  expect_equal(c(p$sorptivity, p$transmissivity), c(2, 0.3), tolerance = 1e-8)
  h <- fit[fit$test == "horton" & fit$model == "horton", ]
  expect_equal(c(h$f0, h$fc, h$k), c(3, 0.5, 0.1), tolerance = 1e-4)
  k <- fit[fit$test == "kostiakov" & fit$model == "kostiakov", ]
  expect_equal(c(k$a, k$b), c(1.5, 0.7), tolerance = 1e-5)
  best <- vapply(split(fit, fit$test), function(x) x$model[which.min(x$aic)], character(1))
  expect_equal(unname(best), names(best))
})

test_that("soil_infiltration_fit() agrees with nls for Horton", {
  t <- c(1, 2, 5, 10, 20, 30, 45, 60, 90)
  i <- c(4.1, 6.2, 10.4, 15.6, 23.9, 31.0, 40.6, 49.7, 67.2)
  fit <- soil_infiltration_fit(data.frame(t = t, i = i), t, i, models = "horton")
  ref <- stats::nls(i ~ fc * t + (f0 - fc) / k * (1 - exp(-k * t)),
    start = list(f0 = 3, fc = 0.6, k = 0.2))
  expect_equal(c(fit$f0, fit$fc, fit$k), unname(stats::coef(ref)), tolerance = 1e-4)
  expect_error(soil_infiltration_fit(data.frame(t = 1:3, i = 1:3), t, i), "4 readings")
  expect_error(soil_infiltration_fit(data.frame(t = t, i = i), t, i, models = "green"))
})
