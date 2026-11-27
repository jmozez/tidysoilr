spec <- soil_spectra(soil_nir, dplyr::starts_with("nm_"))
bands <- attr(spec, "soil_meta")$spectra

test_that("absorbance transform", {
  ab <- soil_absorbance(spec)
  expect_equal(attr(ab, "soil_meta")$type, "absorbance")
  expect_equal(ab$nm_400, log10(1 / soil_nir$nm_400))
  expect_error(soil_absorbance(ab), "reflectance")
  expect_error(soil_absorbance(soil_nir), "soil_spectra")
})

test_that("SNV gives rows with mean 0 and sd 1", {
  s <- soil_snv(spec)
  m <- as.matrix(as.data.frame(s)[bands])
  expect_equal(unname(rowMeans(m)), rep(0, nrow(m)), tolerance = 1e-10)
  expect_equal(unname(apply(m, 1, sd)), rep(1, nrow(m)), tolerance = 1e-10)
  expect_equal(attr(s, "soil_meta")$history, "snv")
  expect_equal(s$clay, soil_nir$clay)
  expect_equal(names(s), names(soil_nir))
})

test_that("MSC aligns spectra to the reference", {
  s <- soil_msc(spec)
  ref <- attr(s, "msc_reference")
  expect_length(ref, length(bands))
  m <- as.matrix(as.data.frame(s)[bands])
  expect_lt(mean(apply(m, 1, function(x) sd(x - ref))), mean(apply(
    as.matrix(as.data.frame(spec)[bands]), 1, function(x) sd(x - ref))))
  expect_error(soil_msc(spec, reference = 1:3), "one value per band")
})

test_that("Savitzky-Golay smoothing and derivatives", {
  x <- seq(400, 2480, by = 20)
  lin <- tibble::as_tibble(as.data.frame(matrix(rep(2 + 0.5 * seq_along(x), 2), nrow = 2,
    byrow = TRUE, dimnames = list(NULL, paste0("nm_", x)))))
  ls <- soil_spectra(lin, dplyr::everything())
  sm <- soil_savgol(ls, window = 5, poly = 2, deriv = 0)
  expect_equal(ncol(sm), length(x) - 4)
  expect_equal(unname(unlist(sm[1, ])), 2 + 0.5 * (3:(length(x) - 2)))
  d1 <- soil_savgol(ls, window = 5, poly = 2, deriv = 1)
  expect_equal(unname(unlist(d1[1, ])), rep(0.5, length(x) - 4), tolerance = 1e-10)
  d2 <- soil_savgol(ls, window = 7, poly = 3, deriv = 2)
  expect_equal(unname(unlist(d2[1, ])), rep(0, length(x) - 6), tolerance = 1e-10)
  expect_error(soil_savgol(spec, window = 10), "odd")
  expect_error(soil_savgol(spec, window = 5, poly = 5), "less than")
  expect_error(soil_savgol(spec, window = 5, poly = 1, deriv = 2), "exceed")
  uneven <- soil_spectra(soil_nir, c(nm_400, nm_420, nm_480, nm_500, nm_520))
  expect_error(soil_savgol(uneven, window = 3), "evenly")
})

test_that("resampling and continuum removal", {
  r <- soil_resample(spec, seq(500, 2400, by = 100))
  expect_equal(attr(r, "soil_meta")$wavelengths, seq(500, 2400, by = 100))
  expect_equal(r$nm_500, soil_nir$nm_500)
  expect_true("nm_2400" %in% names(r))
  expect_error(soil_resample(spec, c(300, 500, 600)), "within")
  cr <- soil_continuum_removal(spec)
  m <- as.matrix(as.data.frame(cr)[bands])
  expect_true(all(m <= 1 + 1e-9))
  expect_equal(unname(m[, 1]), rep(1, nrow(m)))
  expect_error(soil_continuum_removal(soil_snv(spec)), "reflectance")
})

test_that("preprocessing chains with the pipe", {
  out <- spec |> soil_absorbance() |> soil_savgol(window = 11, deriv = 1) |> soil_snv()
  expect_s3_class(out, "soil_spectra")
  expect_equal(attr(out, "soil_meta")$history, c("absorbance", "savgol(w=11,p=2,d=1)", "snv"))
})
