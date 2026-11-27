nir <- soil_nir[c("clay", grep("^nm_", names(soil_nir), value = TRUE))]
x <- nir[-1]

test_that("soil_pls_fit() matches least squares when components are complete", {
  set.seed(4)
  df <- data.frame(a = rnorm(30), b = rnorm(30), c = rnorm(30))
  y <- 1 + 2 * df$a - df$b + 0.5 * df$c + rnorm(30, 0, 0.1)
  fit <- soil_pls_fit(df, y, num_comp = 3)
  ols <- lm(y ~ ., data = df)
  expect_equal(predict(fit, df), unname(fitted(ols)), tolerance = 1e-8)
  expect_equal(tidy(fit)$estimate, unname(coef(ols)), tolerance = 1e-8)
  scaled <- soil_pls_fit(df, y, num_comp = 3, scale = TRUE)
  expect_equal(predict(scaled, df), unname(fitted(ols)), tolerance = 1e-8)
})

test_that("soil_pls_fit() calibrates the simulated spectra", {
  fit <- soil_pls_fit(x, nir$clay, num_comp = 6)
  expect_s3_class(fit, "soil_pls_fit")
  expect_gt(fit$r_squared, 0.9)
  expect_length(predict(fit, x, num_comp = 2), nrow(x))
  expect_true(all(diff(fit$explained_y) >= -1e-12))
  expect_error(predict(fit, x, num_comp = 10), "exceeds")
  expect_error(predict(fit, x[1:3]), "missing")
  expect_error(soil_pls_fit(x, nir$clay[1:5]), "one value per row")
  expect_error(soil_pls_fit(data.frame(a = c(1, NA), b = 1:2), 1:2), "missing")
  expect_output(print(fit), "PLS regression")
})

test_that("soil_pls() works through parsnip", {
  spec <- soil_pls(num_comp = 5) |> parsnip::set_engine("nipals")
  expect_s3_class(spec, "soil_pls")
  expect_output(print(spec), "Partial least squares")
  fit <- parsnip::fit(spec, clay ~ ., data = nir)
  pr <- predict(fit, nir[1:4, ])
  expect_named(pr, ".pred")
  expect_equal(pr$.pred, predict(parsnip::extract_fit_engine(fit), nir[1:4, ]))
  up <- update(spec, num_comp = 3)
  expect_equal(rlang::eval_tidy(up$args$num_comp), 3)
  expect_error(soil_pls(mode = "classification"), "regression")
})

test_that("soil_pls() works in a workflow with recipes steps", {
  skip_if_not_installed("workflows")
  rec <- recipes::recipe(clay ~ ., data = nir) |>
    step_soil_snv(starts_with("nm_")) |>
    step_soil_savgol(starts_with("nm_"), window = 7, deriv = 1)
  wf <- workflows::workflow(rec, soil_pls(num_comp = 4) |> parsnip::set_engine("nipals"))
  fit <- parsnip::fit(wf, data = nir)
  pr <- predict(fit, nir)
  expect_gt(cor(pr$.pred, nir$clay), 0.9)
})
