test_that("soil_profile() builds a classed, sorted tibble", {
  prof <- soil_profile(png_profiles, id = profile, top = top, bottom = bottom, horizon = horizon)
  expect_s3_class(prof, "soil_profile")
  expect_s3_class(prof, "tbl_df")
  expect_true(is_soil_profile(prof))
  expect_false(is_soil_spectra(prof))
  expect_equal(nrow(prof), nrow(png_profiles))
  expect_false(is.unsorted(prof$profile))
  expect_snapshot(print(prof, n = 3, width = 80))
})

test_that("soil_profile() keeps its class through dplyr verbs", {
  prof <- soil_profile(png_profiles, profile, top, bottom)
  expect_s3_class(dplyr::filter(prof, top < 30), "soil_profile")
  expect_s3_class(dplyr::mutate(prof, oc_pct = oc / 10), "soil_profile")
  expect_s3_class(dplyr::arrange(prof, dplyr::desc(top)), "soil_profile")
  expect_s3_class(prof[1:5, ], "soil_profile")
  expect_s3_class(vctrs::vec_rbind(prof[1:3, ], prof[4:6, ]), "soil_profile")
  dropped <- dplyr::select(prof, profile, clay)
  expect_false(inherits(dropped, "soil_profile"))
  expect_s3_class(dropped, "tbl_df")
})

test_that("renaming keeps metadata in step", {
  prof <- soil_profile(png_profiles, profile, top, bottom)
  renamed <- dplyr::rename(prof, upper = top)
  expect_s3_class(renamed, "soil_profile")
  expect_equal(attr(renamed, "soil_meta")$top, "upper")
  expect_s3_class(soil_carbon_stock(renamed, oc, bulk_density), "soil_profile")
})

test_that("soil_profile() validates input", {
  expect_error(soil_profile(1:3, a, b, c), "data frame")
  bad <- data.frame(id = c("a", NA), top = c(0, 10), bottom = c(10, 20))
  expect_error(soil_profile(bad, id, top, bottom), "missing")
  chr <- data.frame(id = "a", top = "0", bottom = 10)
  expect_error(soil_profile(chr, id, top, bottom), "numeric")
  inv <- data.frame(id = "a", top = c(0, 20), bottom = c(10, 15))
  expect_warning(soil_profile(inv, id, top, bottom), "lower depth")
})

test_that("soil_spectra() reads wavelengths and validates", {
  spec <- soil_spectra(soil_nir, dplyr::starts_with("nm_"))
  expect_s3_class(spec, "soil_spectra")
  meta <- attr(spec, "soil_meta")
  expect_equal(meta$wavelengths, seq(400, 2480, by = 20))
  expect_equal(meta$type, "reflectance")
  expect_s3_class(dplyr::filter(spec, clay > 30), "soil_spectra")
  expect_false(inherits(dplyr::select(spec, sample, clay), "soil_spectra"))
  expect_error(soil_spectra(soil_nir, c(clay, oc)), "at least three")
  expect_error(soil_spectra(soil_nir, c(sample, nm_400, nm_420)), "numeric")
  expect_error(soil_spectra(soil_nir, dplyr::starts_with("nm_"), wavelengths = 1:3), "one value")
  expect_error(
    soil_spectra(soil_nir, c(nm_420, nm_400, nm_440)),
    "increasing"
  )
  named <- data.frame(a = 1, b = 2, c = 3)
  expect_error(soil_spectra(named, a:c), "wavelengths")
  expect_snapshot(print(spec, n = 2, width = 80))
})

test_that("soil_as_aqp() converts to a SoilProfileCollection", {
  skip_if_not_installed("aqp")
  prof <- soil_profile(png_profiles, profile, top, bottom, horizon = horizon)
  spc <- soil_as_aqp(prof)
  expect_s4_class(spc, "SoilProfileCollection")
  expect_equal(length(spc), length(unique(png_profiles$profile)))
  expect_error(soil_as_aqp(png_profiles), "soil_profile")
})
