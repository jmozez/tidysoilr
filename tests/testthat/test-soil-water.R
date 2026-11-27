test_that("soil_vg_theta() evaluates the model", {
  df <- data.frame(h = c(0, 100, 15000))
  res <- soil_vg_theta(df, h, theta_r = 0.05, theta_s = 0.45, alpha = 0.02, n = 1.5)
  expect_equal(res$theta_vg[1], 0.45)
  expect_true(all(diff(res$theta_vg) < 0))
  expect_error(soil_vg_theta(df, h, 0.05, 0.45, 0.02, n = 1), "greater than 1")
  df$a <- c(0.01, 0.02, 0.03)
  expect_length(soil_vg_theta(df, h, 0.05, 0.45, a, 1.5)$theta_vg, 3)
})

test_that("soil_vg_fit() recovers known parameters", {
  h <- c(0, 10, 30, 60, 100, 330, 1000, 3000, 15000)
  df <- data.frame(h = h)
  df <- soil_vg_theta(df, h, 0.078, 0.43, 0.036, 1.56, .name = "theta")
  fit <- soil_vg_fit(df, h, theta)
  expect_s3_class(fit, "soil_retention")
  expect_equal(fit$theta_s, 0.43, tolerance = 1e-3)
  expect_equal(fit$alpha, 0.036, tolerance = 1e-2)
  expect_equal(fit$n, 1.56, tolerance = 1e-2)
  expect_true(fit$converged)
  expect_gt(fit$r_squared, 0.999)
})

test_that("soil_vg_fit() fits grouped curves and supports fixed theta_r", {
  fits <- soil_vg_fit(soil_retention_lab, h_cm, theta, .by = sample)
  expect_equal(nrow(fits), 12)
  expect_true(all(fits$converged))
  expect_true(all(fits$rmse < 0.02))
  grouped <- soil_vg_fit(dplyr::group_by(soil_retention_lab, sample), h_cm, theta)
  expect_equal(grouped$alpha, fits$alpha)
  fixed <- soil_vg_fit(soil_retention_lab, h_cm, theta, .by = sample, fix_theta_r = 0)
  expect_true(all(fixed$theta_r == 0))
  expect_output(print(fits[1:2, ]), "van Genuchten")
})

test_that("soil_vg_fit() validates input", {
  expect_error(soil_vg_fit(soil_retention_lab, h_cm, h_kpa), "volumetric")
  few <- data.frame(id = c("a", "a", "b", "b", "b", "b", "b", "b"),
    h = c(0, 100, 0, 10, 100, 330, 1000, 15000),
    theta = c(0.4, 0.3, 0.45, 0.44, 0.35, 0.3, 0.25, 0.15))
  expect_warning(res <- soil_vg_fit(few, h, theta, .by = id), "too few")
  expect_true(is.na(res$alpha[1]))
})

test_that("predict() and soil_awc() use the fitted curves", {
  fits <- soil_vg_fit(soil_retention_lab, h_cm, theta, .by = sample)
  pr <- predict(fits, h = c(330, 15000))
  expect_equal(nrow(pr), 24)
  expect_named(pr, c("sample", "h_cm", ".pred"))
  nd <- data.frame(sample = "S07", h_cm = 330)
  expect_equal(predict(fits, new_data = nd)$.pred, pr$.pred[pr$sample == "S07" & pr$h_cm == 330])
  expect_error(predict(fits, new_data = data.frame(h_cm = 1)), "missing")
  awc <- soil_awc(fits, thickness = 30)
  expect_true(all(awc$awc > 0))
  expect_equal(awc$awc_mm, awc$awc * 300)
  expect_s3_class(awc, "soil_retention")
  expect_error(soil_awc(soil_retention_lab), "soil_retention")
  expect_error(soil_awc(fits, fc = 20000), "higher")
})
