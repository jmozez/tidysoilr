# Soil spectroscopy preprocessing ------------------------------------------

check_spectra <- function(x, arg = rlang::caller_arg(x), call = rlang::caller_env()) {
  if (!is_soil_spectra(x)) {
    cli::cli_abort(c(
      "{.arg {arg}} must be a {.cls soil_spectra} table.",
      "i" = "Create one with {.fn soil_spectra}."
    ), call = call)
  }
  invisible(x)
}

spectra_matrix <- function(x) {
  meta <- soil_meta(x)
  m <- as.matrix(as.data.frame(x)[meta$spectra])
  storage.mode(m) <- "double"
  m
}

spectra_prefix <- function(cols) {
  p <- sub("[0-9.]+$", "", cols[1])
  if (identical(p, cols[1])) "" else p
}

# Replace the spectral block of a soil_spectra table with a new matrix.
set_spectra <- function(x, m, wavelengths, step, type = NULL, cols = NULL) {
  meta <- soil_meta(x)
  data <- as_plain_tibble(x)
  old <- meta$spectra
  if (is.null(cols)) {
    keep <- match(wavelengths, meta$wavelengths)
    cols <- ifelse(is.na(keep), paste0(spectra_prefix(old), format(wavelengths, trim = TRUE)),
      old[keep])
  }
  colnames(m) <- cols
  first <- match(old[1], names(data))
  other <- setdiff(names(data), old)
  before <- other[seq_along(other) < first]
  after <- setdiff(other, before)
  out <- dplyr::bind_cols(data[before], tibble::as_tibble(m), data[after])
  meta$spectra <- cols
  meta$wavelengths <- wavelengths
  meta$required <- cols
  meta$history <- c(meta$history, step)
  if (!is.null(type)) meta$type <- type
  new_soil_tbl(out, meta, "soil_spectra")
}

snv_matrix <- function(m) {
  mu <- rowMeans(m, na.rm = TRUE)
  s <- apply(m, 1, stats::sd, na.rm = TRUE)
  s[!is.finite(s) | s == 0] <- 1
  (m - mu) / s
}

savgol_coefs <- function(window, poly, deriv) {
  h <- (window - 1) / 2
  z <- -h:h
  x <- outer(z, 0:poly, `^`)
  cmat <- solve(crossprod(x), t(x))
  factorial(deriv) * cmat[deriv + 1, ]
}

savgol_matrix <- function(m, window, poly, deriv) {
  cf <- savgol_coefs(window, poly, deriv)
  h <- (window - 1) / 2
  p <- ncol(m)
  idx <- (h + 1):(p - h)
  out <- matrix(0, nrow(m), length(idx))
  for (j in seq_len(window)) {
    out <- out + cf[j] * m[, idx - h - 1 + j, drop = FALSE]
  }
  list(m = out, idx = idx)
}

check_savgol <- function(window, poly, deriv, n_bands, call = rlang::caller_env()) {
  check_integerish(window, min = 3, call = call)
  check_integerish(poly, min = 0, call = call)
  check_integerish(deriv, min = 0, call = call)
  if (window %% 2 != 1) {
    cli::cli_abort("{.arg window} must be odd, not {window}.", call = call)
  }
  if (poly >= window) {
    cli::cli_abort("{.arg poly} must be less than {.arg window}.", call = call)
  }
  if (deriv > poly) {
    cli::cli_abort("{.arg deriv} must not exceed {.arg poly}.", call = call)
  }
  if (!is.null(n_bands) && window > n_bands) {
    cli::cli_abort("{.arg window} ({window}) is wider than the spectrum ({n_bands} bands).",
      call = call)
  }
}

#' Convert reflectance to absorbance
#'
#' Transforms reflectance (0 to 1) to apparent absorbance
#' \eqn{A = \log_{10}(1/R)}, which is closer to linear in the concentration
#' of absorbing constituents.
#'
#' @param .data A [soil_spectra()] table of reflectance.
#' @return A `<soil_spectra>` of absorbance.
#' @family soil spectroscopy
#' @export
#' @examples
#' spec <- soil_spectra(soil_nir, dplyr::starts_with("nm_"))
#' soil_absorbance(spec)
soil_absorbance <- function(.data) {
  check_spectra(.data)
  meta <- soil_meta(.data)
  if (!identical(meta$type, "reflectance")) {
    cli::cli_abort("Spectra must be reflectance, not {meta$type}.")
  }
  m <- spectra_matrix(.data)
  if (any(m <= 0, na.rm = TRUE) || any(m > 1.5, na.rm = TRUE)) {
    cli::cli_abort("Reflectance must be positive fractions (0 to 1); divide percentages by 100.")
  }
  set_spectra(.data, log10(1 / m), meta$wavelengths, "absorbance", type = "absorbance",
    cols = meta$spectra)
}

#' Standard normal variate
#'
#' Centres and scales each spectrum by its own mean and standard deviation
#' (Barnes et al. 1989), removing multiplicative scatter and baseline
#' effects from particle size and surface roughness.
#'
#' @param .data A [soil_spectra()] table.
#' @return A `<soil_spectra>` with corrected spectra.
#' @references Barnes, R. J., Dhanoa, M. S. and Lister, S. J. (1989).
#'   Standard normal variate transformation and de-trending of
#'   near-infrared diffuse reflectance spectra. *Applied Spectroscopy*, 43,
#'   772-777. \doi{10.1366/0003702894202201}
#' @family soil spectroscopy
#' @export
#' @examples
#' spec <- soil_spectra(soil_nir, dplyr::starts_with("nm_"))
#' soil_snv(spec)
soil_snv <- function(.data) {
  check_spectra(.data)
  meta <- soil_meta(.data)
  set_spectra(.data, snv_matrix(spectra_matrix(.data)), meta$wavelengths, "snv",
    type = "preprocessed", cols = meta$spectra)
}

#' Multiplicative scatter correction
#'
#' Regresses each spectrum on a reference spectrum (by default the mean
#' spectrum) and removes the fitted offset and slope (Geladi et al. 1985).
#' Pass the reference from calibration samples when correcting new samples.
#'
#' @param .data A [soil_spectra()] table.
#' @param reference Optional reference spectrum (numeric vector with one
#'   value per band). Defaults to the mean spectrum of `.data`.
#' @return A `<soil_spectra>` with corrected spectra. The reference used is
#'   stored in the attribute `"msc_reference"`.
#' @references Geladi, P., MacDougall, D. and Martens, H. (1985).
#'   Linearization and scatter-correction for near-infrared reflectance
#'   spectra of meat. *Applied Spectroscopy*, 39, 491-500.
#'   \doi{10.1366/0003702854248656}
#' @family soil spectroscopy
#' @export
#' @examples
#' spec <- soil_spectra(soil_nir, dplyr::starts_with("nm_"))
#' soil_msc(spec)
soil_msc <- function(.data, reference = NULL) {
  check_spectra(.data)
  meta <- soil_meta(.data)
  m <- spectra_matrix(.data)
  ref <- reference %||% colMeans(m, na.rm = TRUE)
  if (!is.numeric(ref) || length(ref) != ncol(m)) {
    cli::cli_abort("{.arg reference} must be numeric with one value per band ({ncol(m)}).")
  }
  out <- t(apply(m, 1, function(x) {
    fit <- stats::lm.fit(cbind(1, ref), x)$coefficients
    (x - fit[1]) / fit[2]
  }))
  res <- set_spectra(.data, out, meta$wavelengths, "msc", type = "preprocessed",
    cols = meta$spectra)
  attr(res, "msc_reference") <- ref
  res
}

#' Savitzky-Golay smoothing and derivatives
#'
#' Fits a local polynomial of degree `poly` in a moving window of `window`
#' bands and returns the smoothed spectrum (`deriv = 0`) or its first or
#' second derivative (Savitzky and Golay 1964). Derivatives remove baseline
#' offsets (first) and slopes (second) and sharpen overlapping absorption
#' features. Derivatives are per band, as in `prospectr::savitzkyGolay()`.
#' The `(window - 1) / 2` bands at each end cannot be computed and are
#' dropped.
#'
#' @param .data A [soil_spectra()] table with evenly spaced bands.
#' @param window Window width in bands (odd).
#' @param poly Polynomial degree.
#' @param deriv Derivative order (0, 1 or 2).
#' @return A `<soil_spectra>` with the filtered spectra.
#' @references Savitzky, A. and Golay, M. J. E. (1964). Smoothing and
#'   differentiation of data by simplified least squares procedures.
#'   *Analytical Chemistry*, 36, 1627-1639. \doi{10.1021/ac60214a047}
#' @family soil spectroscopy
#' @export
#' @examples
#' spec <- soil_spectra(soil_nir, dplyr::starts_with("nm_"))
#' soil_savgol(spec, window = 11, poly = 2, deriv = 1)
soil_savgol <- function(.data, window = 11, poly = 2, deriv = 0) {
  check_spectra(.data)
  meta <- soil_meta(.data)
  check_savgol(window, poly, deriv, length(meta$wavelengths))
  sp <- diff(meta$wavelengths)
  if (max(sp) - min(sp) > 1e-6 * max(sp)) {
    cli::cli_abort(c(
      "Savitzky-Golay filtering needs evenly spaced bands.",
      "i" = "Use {.fn soil_resample} first."
    ))
  }
  res <- savgol_matrix(spectra_matrix(.data), window, poly, deriv)
  step <- paste0("savgol(w=", window, ",p=", poly, ",d=", deriv, ")")
  set_spectra(.data, res$m, meta$wavelengths[res$idx], step, type = "preprocessed",
    cols = meta$spectra[res$idx])
}

#' Resample spectra to new wavelengths
#'
#' Linearly interpolates every spectrum to a new set of wavelengths, for
#' example to align instruments or to thin a 1 nm spectrum to 10 nm.
#'
#' @param .data A [soil_spectra()] table.
#' @param wavelengths New wavelengths, within the measured range.
#' @return A `<soil_spectra>` at the new wavelengths.
#' @family soil spectroscopy
#' @export
#' @examples
#' spec <- soil_spectra(soil_nir, dplyr::starts_with("nm_"))
#' soil_resample(spec, seq(500, 2400, by = 50))
soil_resample <- function(.data, wavelengths) {
  check_spectra(.data)
  meta <- soil_meta(.data)
  if (!is.numeric(wavelengths) || length(wavelengths) < 3 || any(diff(wavelengths) <= 0)) {
    cli::cli_abort("{.arg wavelengths} must be an increasing numeric vector of at least 3 values.")
  }
  rng <- range(meta$wavelengths)
  if (min(wavelengths) < rng[1] || max(wavelengths) > rng[2]) {
    cli::cli_abort("{.arg wavelengths} must lie within the measured range {rng[1]} to {rng[2]}.")
  }
  m <- spectra_matrix(.data)
  out <- t(apply(m, 1, function(y) stats::approx(meta$wavelengths, y, xout = wavelengths)$y))
  if (length(wavelengths) == 1) out <- t(out)
  cols <- paste0(spectra_prefix(meta$spectra), format(wavelengths, trim = TRUE))
  set_spectra(.data, out, wavelengths, "resample", cols = cols)
}

# Upper convex hull of (x, y), evaluated at x.
upper_hull <- function(x, y) {
  hull <- integer()
  for (i in seq_along(x)) {
    while (length(hull) >= 2) {
      a <- hull[length(hull) - 1]
      b <- hull[length(hull)]
      cross <- (x[b] - x[a]) * (y[i] - y[a]) - (y[b] - y[a]) * (x[i] - x[a])
      if (cross >= 0) hull <- hull[-length(hull)] else break
    }
    hull <- c(hull, i)
  }
  stats::approx(x[hull], y[hull], xout = x)$y
}

#' Continuum removal
#'
#' Divides each reflectance spectrum by its continuum, the upper convex
#' hull, so that absorption features (for example clay minerals near 2200
#' nm and water near 1400 and 1900 nm) are expressed as band depths
#' relative to 1 (Clark and Roush 1984).
#'
#' @param .data A [soil_spectra()] table of reflectance.
#' @return A `<soil_spectra>` with continuum-removed spectra (1 = no
#'   absorption).
#' @references Clark, R. N. and Roush, T. L. (1984). Reflectance
#'   spectroscopy: quantitative analysis techniques for remote sensing
#'   applications. *Journal of Geophysical Research*, 89, 6329-6340.
#'   \doi{10.1029/JB089iB07p06329}
#' @family soil spectroscopy
#' @export
#' @examples
#' spec <- soil_spectra(soil_nir, dplyr::starts_with("nm_"))
#' soil_continuum_removal(spec)
soil_continuum_removal <- function(.data) {
  check_spectra(.data)
  meta <- soil_meta(.data)
  if (!identical(meta$type, "reflectance")) {
    cli::cli_abort("Continuum removal needs reflectance spectra, not {meta$type}.")
  }
  m <- spectra_matrix(.data)
  wl <- meta$wavelengths
  out <- t(apply(m, 1, function(y) y / upper_hull(wl, y)))
  set_spectra(.data, out, wl, "continuum removal", type = "preprocessed", cols = meta$spectra)
}
