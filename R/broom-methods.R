# broom methods -----------------------------------------------------------------

#' Tidy, glance at and augment water retention fits
#'
#' * `tidy()` returns one row per curve and parameter.
#' * `glance()` returns a one-row summary of all curves.
#' * `augment()` returns the measured points with fitted values and
#'   residuals.
#'
#' @param x A `<soil_retention>` from [soil_vg_fit()].
#' @param ... Not used.
#' @return A tibble.
#' @family broom methods
#' @name soil_retention_broom
#' @examples
#' fits <- soil_vg_fit(soil_retention_lab, h_cm, theta, .by = sample)
#' tidy(fits)
#' glance(fits)
#' augment(fits)
NULL

#' @rdname soil_retention_broom
#' @export
tidy.soil_retention <- function(x, ...) {
  meta <- soil_meta(x)
  fits <- as_plain_tibble(x)
  terms <- c("theta_r", "theta_s", "alpha", "n")
  idx <- rep(seq_len(nrow(fits)), each = length(terms))
  out <- fits[idx, meta$groups, drop = FALSE]
  out$term <- rep(terms, nrow(fits))
  out$estimate <- as.vector(t(as.matrix(fits[terms])))
  out
}

#' @rdname soil_retention_broom
#' @export
glance.soil_retention <- function(x, ...) {
  tibble::tibble(
    n_curves = nrow(x), n_converged = sum(x$converged, na.rm = TRUE),
    n_obs = sum(x$n_obs), mean_rmse = mean(x$rmse, na.rm = TRUE),
    median_r_squared = stats::median(x$r_squared, na.rm = TRUE)
  )
}

#' @rdname soil_retention_broom
#' @export
augment.soil_retention <- function(x, ...) {
  meta <- soil_meta(x)
  obs <- meta$data
  pred <- stats::predict(x, new_data = obs)
  obs$.fitted <- pred$.pred
  obs$.resid <- obs[[meta$theta]] - obs$.fitted
  obs
}

#' Tidy, glance at and augment PLS fits
#'
#' * `tidy()` returns the regression coefficients (`type = "coefficients"`)
#'   or the variable importance in projection scores (`type = "vip"`; values
#'   above 1 are conventionally important).
#' * `glance()` returns a one-row summary.
#' * `augment()` adds predictions to `new_data`.
#'
#' @param x A `soil_pls_fit` (or a parsnip fit of [soil_pls()]; use
#'   `parsnip::extract_fit_engine()`).
#' @param type `"coefficients"` or `"vip"`.
#' @param num_comp Number of components for the coefficients.
#' @param new_data A data frame with the predictors.
#' @param ... Not used.
#' @return A tibble.
#' @family broom methods
#' @name soil_pls_broom
#' @examples
#' nir <- soil_nir[grep("^nm_", names(soil_nir))]
#' fit <- soil_pls_fit(nir, soil_nir$clay, num_comp = 5)
#' tidy(fit)
#' tidy(fit, type = "vip")
#' glance(fit)
#' augment(fit, soil_nir)[c("clay", ".pred")]
NULL

#' @rdname soil_pls_broom
#' @export
tidy.soil_pls_fit <- function(x, type = c("coefficients", "vip"), num_comp = NULL, ...) {
  type <- rlang::arg_match(type)
  if (type == "vip") {
    return(tibble::tibble(term = x$predictors, vip = unname(x$vip)))
  }
  nc <- num_comp %||% x$num_comp
  check_integerish(nc)
  if (nc > x$num_comp) {
    cli::cli_abort("{.arg num_comp} ({nc}) exceeds the {x$num_comp} component{?s} fitted.")
  }
  tibble::tibble(term = c("(Intercept)", x$predictors),
    estimate = c(x$y_mean - sum(x$center / x$scale * x$coefficients[, nc]),
      unname(x$coefficients[, nc] / x$scale)))
}

#' @rdname soil_pls_broom
#' @export
glance.soil_pls_fit <- function(x, ...) {
  tibble::tibble(num_comp = x$num_comp, n = x$n, n_predictors = length(x$predictors),
    r_squared = x$r_squared, rmse = x$rmse)
}

#' @rdname soil_pls_broom
#' @export
augment.soil_pls_fit <- function(x, new_data, num_comp = NULL, ...) {
  check_data_frame(new_data)
  out <- as_plain_tibble(new_data)
  out$.pred <- stats::predict(x, new_data, num_comp = num_comp)
  out
}

#' Glance at a soil profile table
#'
#' @param x A [soil_profile()].
#' @param ... Not used.
#' @return A one-row tibble with the number of profiles and layers, the
#'   median and maximum profile depth and the number of depth problems found
#'   by [soil_profile_check()].
#' @family broom methods
#' @export
#' @examples
#' glance(soil_profile(png_profiles, profile, top, bottom))
glance.soil_profile <- function(x, ...) {
  meta <- soil_meta(x)
  depth <- tapply(x[[meta$bottom]], x[[meta$id]], max, na.rm = TRUE)
  tibble::tibble(
    n_profiles = length(unique(x[[meta$id]])), n_layers = nrow(x),
    median_depth = stats::median(depth), max_depth = max(depth),
    n_issues = nrow(soil_profile_check(x))
  )
}

#' Tidy a soil spectra table into long format
#'
#' @param x A [soil_spectra()] table.
#' @param id <[`tidy-select`][tidyselect::language]> Optional sample
#'   identifier column to keep; otherwise the row number is used.
#' @param ... Not used.
#' @return A tibble with columns `.row` (or the identifier), `wavelength`
#'   and `value`.
#' @family broom methods
#' @export
#' @examples
#' spec <- soil_spectra(soil_nir, dplyr::starts_with("nm_"))
#' tidy(spec, id = sample)
tidy.soil_spectra <- function(x, id = NULL, ...) {
  meta <- soil_meta(x)
  id_col <- select_one(rlang::enquo(id), x, "id", required = FALSE)
  m <- spectra_matrix(x)
  key <- if (is.null(id_col)) seq_len(nrow(m)) else x[[id_col]]
  out <- tibble::tibble(.row = rep(key, times = ncol(m)),
    wavelength = rep(meta$wavelengths, each = nrow(m)), value = as.vector(m))
  if (!is.null(id_col)) names(out)[1] <- id_col
  out
}
