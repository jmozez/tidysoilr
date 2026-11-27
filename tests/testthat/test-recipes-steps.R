library(recipes)
nir <- soil_nir[c("clay", grep("^nm_", names(soil_nir), value = TRUE))]

test_that("step_soil_snv() matches soil_snv()", {
  rec <- recipe(clay ~ ., data = nir) |> step_soil_snv(starts_with("nm_"))
  expect_message(print(rec), "Standard normal variate")
  prepped <- prep(rec)
  baked <- bake(prepped, new_data = NULL)
  direct <- soil_snv(soil_spectra(nir, starts_with("nm_")))
  expect_equal(baked$nm_1000, direct$nm_1000)
  expect_equal(bake(prepped, new_data = nir[1:5, ])$nm_400, direct$nm_400[1:5])
  td <- tidy(prepped, number = 1)
  expect_equal(nrow(td), 105)
  expect_named(td, c("terms", "id"))
  expect_equal(tidy(rec, number = 1)$terms, "starts_with(\"nm_\")")
  expect_equal(required_pkgs(prepped$steps[[1]]), "tidysoilr")
  expect_error(prep(recipe(clay ~ ., data = nir) |> step_soil_snv(nm_400)), "three")
})

test_that("step_soil_savgol() trims edges and filters", {
  rec <- recipe(clay ~ ., data = nir) |> step_soil_savgol(starts_with("nm_"), window = 7,
    deriv = 1)
  prepped <- prep(rec)
  baked <- bake(prepped, new_data = NULL)
  expect_equal(ncol(baked), 1 + 105 - 6)
  direct <- soil_savgol(soil_spectra(nir, starts_with("nm_")), window = 7, deriv = 1)
  expect_equal(baked$nm_1000, direct$nm_1000)
  td <- tidy(prepped, number = 1)
  expect_equal(nrow(td), 99)
  expect_equal(unique(td$window), 7)
  expect_message(print(prepped), "Savitzky-Golay")
  expect_equal(required_pkgs(prepped$steps[[1]]), "tidysoilr")
  expect_error(step_soil_savgol(recipe(clay ~ ., data = nir), starts_with("nm_"), window = 4), "odd")
})

test_that("step_soil_ilr() replaces texture with log-ratios", {
  dat <- na.omit(png_profiles[c("oc", "sand", "silt", "clay")])
  prepped <- recipe(oc ~ sand + silt + clay, data = dat) |> step_soil_ilr(sand, silt, clay) |> prep()
  baked <- bake(prepped, new_data = NULL)
  expect_named(baked, c("oc", "ilr_1", "ilr_2"))
  direct <- soil_texture_ilr(dat, sand, silt, clay)
  expect_equal(baked$ilr_1, direct$ilr_1)
  kept <- recipe(oc ~ ., data = dat) |>
    step_soil_ilr(sand, silt, clay, keep_original_cols = TRUE) |>
    prep() |>
    bake(new_data = NULL)
  expect_true(all(c("sand", "ilr_1") %in% names(kept)))
  expect_equal(tidy(prepped, number = 1)$terms, c("sand", "silt", "clay"))
  expect_message(print(prepped), "Texture log-ratios")
  expect_equal(required_pkgs(prepped$steps[[1]]), "tidysoilr")
  expect_error(prep(recipe(oc ~ ., data = dat) |> step_soil_ilr(sand, silt)), "exactly three")
})

test_that("step_soil_texture() adds a fixed-level factor", {
  dat <- na.omit(png_profiles[c("oc", "sand", "silt", "clay")])
  prepped <- recipe(oc ~ ., data = dat) |> step_soil_texture(sand, silt, clay) |> prep()
  baked <- bake(prepped, new_data = dat[1:3, ])
  expect_s3_class(baked$texture_class, "factor")
  expect_equal(nlevels(baked$texture_class), 12)
  expect_equal(baked$texture_class, soil_texture_class(dat[1:3, ], sand, silt, clay)$texture)
  fao <- recipe(oc ~ ., data = dat) |>
    step_soil_texture(sand, silt, clay, system = "fao", keep_original_cols = FALSE) |>
    prep() |>
    bake(new_data = NULL)
  expect_false("sand" %in% names(fao))
  expect_equal(nlevels(fao$texture_class), 5)
  expect_equal(tidy(prepped, number = 1)$system, rep("usda", 3))
  expect_message(print(prepped), "USDA texture class")
  expect_equal(required_pkgs(prepped$steps[[1]]), "tidysoilr")
  expect_error(prep(recipe(oc ~ ., data = dat) |> step_soil_texture(sand)), "exactly three")
})
