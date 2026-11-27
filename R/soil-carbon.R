# Equal-area depth splines and carbon stocks ------------------------------

gauss_nodes <- c(-sqrt(3 / 5), 0, sqrt(3 / 5))
gauss_weights <- c(5 / 9, 8 / 9, 5 / 9)

# Integral of every quadratic B-spline basis function over [lo, hi]. The
# basis is piecewise quadratic, so three-point Gauss-Legendre quadrature on
# each knot interval is exact.
basis_integral <- function(lo, hi, knots, breaks, derivs = 0) {
  cuts <- sort(unique(c(lo, hi, breaks[breaks > lo & breaks < hi])))
  total <- 0
  for (k in seq_len(length(cuts) - 1)) {
    a <- cuts[k]
    b <- cuts[k + 1]
    x <- (a + b) / 2 + (b - a) / 2 * gauss_nodes
    bm <- splines::splineDesign(knots, x, ord = 3, derivs = rep(derivs, 3), outer.ok = TRUE)
    total <- total + colSums(bm * gauss_weights) * (b - a) / 2
  }
  total
}

# Fit the equal-area quadratic spline of Bishop et al. (1999) to one
# profile and return its mean over each target interval.
ea_spline_one <- function(u, v, y, lambda, target_top, target_bottom) {
  ok <- !is.na(u) & !is.na(v) & !is.na(y) & v > u
  u <- u[ok]
  v <- v[ok]
  y <- y[ok]
  out <- rep(NA_real_, length(target_top))
  if (!length(y)) {
    return(out)
  }
  o <- order(u)
  u <- u[o]
  v <- v[o]
  y <- y[o]
  n <- length(y)
  breaks <- sort(unique(c(u, v)))
  a <- breaks[1]
  b <- breaks[length(breaks)]
  inside <- target_top >= a - 1e-9 & target_bottom <= b + 1e-9
  if (n == 1) {
    out[inside] <- y
    return(out)
  }
  interior <- breaks[-c(1, length(breaks))]
  knots <- c(rep(a, 3), interior, rep(b, 3))
  nb <- length(knots) - 3
  amat <- t(vapply(seq_len(n), function(i) {
    basis_integral(u[i], v[i], knots, breaks) / (v[i] - u[i])
  }, numeric(nb)))
  pen <- matrix(0, nb, nb)
  for (k in seq_len(length(breaks) - 1)) {
    lo <- breaks[k]
    hi <- breaks[k + 1]
    x <- (lo + hi) / 2 + (hi - lo) / 2 * gauss_nodes
    d1 <- splines::splineDesign(knots, x, ord = 3, derivs = rep(1, 3), outer.ok = TRUE)
    pen <- pen + crossprod(d1 * sqrt(gauss_weights * (hi - lo) / 2))
  }
  lhs <- crossprod(amat) / n + lambda * pen + diag(1e-10, nb)
  beta <- solve(lhs, crossprod(amat, y) / n)
  for (j in which(inside)) {
    out[j] <- sum(basis_integral(target_top[j], target_bottom[j], knots, breaks) * beta) /
      (target_bottom[j] - target_top[j])
  }
  out
}

#' Equal-area quadratic splines for depth harmonisation
#'
#' Harmonises soil profile data measured on irregular horizons to standard
#' depth intervals (for example the GlobalSoilMap depths 0-5, 5-15, 15-30,
#' 30-60, 60-100 and 100-200 cm) with the mass-preserving equal-area
#' quadratic spline of Ponce-Hernandez et al. (1986) and Bishop et al.
#' (1999), as in `ithir::ea_spline()` and `mpspline2::mpspline()`. The
#' spline minimises
#' \deqn{\frac{1}{n}\sum_{i=1}^{n}(y_i - \bar f_i)^2 + \lambda \int [f'(x)]^2\,dx,}
#' where \eqn{\bar f_i} is the mean of the spline over horizon \eqn{i}, so
#' with small `lambda` the horizon means are reproduced. It is solved here
#' as penalised least squares in a quadratic B-spline basis with knots at
#' the horizon boundaries ([splines::splineDesign()]).
#'
#' @param .data A [soil_profile()] or a data frame of layers.
#' @param ... <[`tidy-select`][tidyselect::language]> Numeric property
#'   columns to harmonise.
#' @param depths Boundaries of the target intervals (cm), increasing. Use
#'   `0:100` for 1 cm resolution.
#' @param lambda Smoothing parameter; 0.1 is the usual default.
#' @param lower,upper Limits for the predicted values (for example 0 for
#'   concentrations, 100 for percentages).
#' @inheritParams soil_profile_check
#'
#' @return A `<soil_profile>` with one row per profile and target interval,
#'   with the identifier, `top`, `bottom` and one column per property.
#'   Intervals that are not fully covered by the measured layers are `NA`.
#' @references Bishop, T. F. A., McBratney, A. B. and Laslett, G. M. (1999).
#'   Modelling soil attribute depth functions with equal-area quadratic
#'   smoothing splines. *Geoderma*, 91, 27-45.
#'   \doi{10.1016/S0016-7061(99)00003-8}
#'
#'   Malone, B. P., McBratney, A. B., Minasny, B. and Laslett, G. M. (2009).
#'   Mapping continuous depth functions of soil carbon storage and available
#'   water capacity. *Geoderma*, 154, 138-152.
#'   \doi{10.1016/j.geoderma.2009.10.007}
#' @family soil carbon
#' @export
#' @examples
#' prof <- soil_profile(png_profiles, profile, top, bottom)
#' soil_spline(prof, oc, clay, depths = c(0, 5, 15, 30, 60, 100))
soil_spline <- function(.data, ..., depths = c(0, 5, 15, 30, 60, 100, 200), lambda = 0.1,
                        lower = 0, upper = Inf, id = NULL, top = NULL, bottom = NULL) {
  check_data_frame(.data)
  check_number_positive(lambda, allow_zero = TRUE)
  if (!is.numeric(depths) || length(depths) < 2 || any(diff(depths) <= 0)) {
    cli::cli_abort("{.arg depths} must be an increasing numeric vector of at least two values.")
  }
  cols <- profile_cols(.data, rlang::enquo(id), rlang::enquo(top), rlang::enquo(bottom))
  data <- as_plain_tibble(.data)
  vars <- select_dots(data, ...)
  if (!length(vars)) {
    cli::cli_abort("Select at least one property column in {.arg ...}.")
  }
  num <- vapply(data[vars], is.numeric, logical(1))
  if (!all(num)) {
    cli::cli_abort("Property columns must be numeric; {.field {vars[!num]}} {?is/are} not.")
  }
  tt <- depths[-length(depths)]
  tb <- depths[-1]
  ids <- unique(data[[cols$id]])
  idx <- split(seq_len(nrow(data)), factor(data[[cols$id]], levels = ids))
  res <- lapply(seq_along(ids), function(k) {
    r <- idx[[k]]
    out <- tibble::tibble(!!cols$id := rep(ids[k], length(tt)), !!cols$top := tt,
      !!cols$bottom := tb)
    for (v in vars) {
      val <- ea_spline_one(data[[cols$top]][r], data[[cols$bottom]][r], data[[v]][r],
        lambda, tt, tb)
      out[[v]] <- pmin(pmax(val, lower), upper)
    }
    out
  })
  out <- dplyr::bind_rows(res)
  meta <- list(id = cols$id, top = cols$top, bottom = cols$bottom, horizon = NULL,
    required = unlist(cols, use.names = FALSE))
  new_soil_tbl(out, meta, "soil_profile")
}

#' Soil organic carbon stock of each layer
#'
#' Computes the organic carbon stock of each layer in Mg C/ha (t/ha) as
#' \deqn{SOC = OC \times \rho_b \times d \times (1 - CF) / 10,}
#' with OC in g/kg, fine-earth bulk density \eqn{\rho_b} in g/cm3, layer
#' thickness \eqn{d} in cm and volumetric coarse fragment fraction CF
#' (Poeplau et al. 2017). Sum the layers over a fixed depth with
#' [soil_carbon_total()].
#'
#' @param .data A [soil_profile()] or a data frame of layers.
#' @param oc <[`tidy-select`][tidyselect::language]> Organic carbon.
#' @param bulk_density <[`tidy-select`][tidyselect::language]> Fine-earth
#'   bulk density (g/cm3).
#' @param coarse <[`tidy-select`][tidyselect::language]> Optional volumetric
#'   coarse fragment content.
#' @param oc_units Units of `oc`: `"g_kg"` or `"percent"`.
#' @param coarse_units Units of `coarse`: `"percent"` or `"fraction"`.
#' @param top,bottom <[`tidy-select`][tidyselect::language]> Depth columns
#'   (cm). Taken from a `<soil_profile>` when not given.
#' @param .name Name of the new column.
#'
#' @return `.data` with a new column of carbon stock (Mg C/ha).
#' @inherit soil_bulk_density references
#' @family soil carbon
#' @export
#' @examples
#' prof <- soil_profile(png_profiles, profile, top, bottom)
#' prof |>
#'   soil_carbon_stock(oc, bulk_density, coarse = coarse_frag) |>
#'   dplyr::select(profile, top, bottom, oc, bulk_density, soc_stock)
soil_carbon_stock <- function(.data, oc, bulk_density, coarse = NULL,
                              oc_units = c("g_kg", "percent"),
                              coarse_units = c("percent", "fraction"),
                              top = NULL, bottom = NULL, .name = "soc_stock") {
  check_data_frame(.data)
  oc_units <- rlang::arg_match(oc_units)
  coarse_units <- rlang::arg_match(coarse_units)
  check_name(.name)
  meta <- if (is_soil_profile(.data)) soil_meta(.data) else list()
  top <- select_one(rlang::enquo(top), .data, "top", required = is.null(meta$top)) %||% meta$top
  bottom <- select_one(rlang::enquo(bottom), .data, "bottom", required = is.null(meta$bottom)) %||%
    meta$bottom
  x <- numeric_col(rlang::enquo(oc), .data, "oc", non_negative = TRUE)
  bd <- numeric_col(rlang::enquo(bulk_density), .data, "bulk_density", non_negative = TRUE)
  cf <- numeric_col(rlang::enquo(coarse), .data, "coarse", required = FALSE, non_negative = TRUE)
  if (oc_units == "percent") x <- x * 10
  cf <- if (is.null(cf)) 0 else if (coarse_units == "percent") cf / 100 else cf
  cf[is.na(cf)] <- 0
  if (any(cf > 1)) {
    cli::cli_abort("Coarse fragment fractions must not exceed 1; check {.arg coarse_units}.")
  }
  thick <- .data[[bottom]] - .data[[top]]
  .data[[.name]] <- x * bd * thick * (1 - cf) / 10
  .data
}

#' Total stock to a fixed depth
#'
#' Sums a per-layer stock (for example from [soil_carbon_stock()]) over a
#' fixed depth interval for every profile, prorating layers that cross the
#' interval limits by thickness. Comparing stocks to a fixed depth, such as
#' 0 to 30 cm for IPCC reporting, needs this proration.
#'
#' @param .data A [soil_profile()] or a data frame of layers.
#' @param stock <[`tidy-select`][tidyselect::language]> Per-layer stock
#'   column.
#' @param from,to Limits of the interval (cm).
#' @inheritParams soil_profile_check
#'
#' @return A tibble with one row per profile: the identifier, `from`, `to`,
#'   the total stock (named after `stock`) and `coverage`, the share of the
#'   interval with data (a total with coverage below 1 is an
#'   underestimate).
#' @family soil carbon
#' @export
#' @examples
#' prof <- soil_profile(png_profiles, profile, top, bottom) |>
#'   soil_carbon_stock(oc, bulk_density)
#' soil_carbon_total(prof, soc_stock, from = 0, to = 30)
soil_carbon_total <- function(.data, stock, from = 0, to = 30, id = NULL, top = NULL,
                              bottom = NULL) {
  check_data_frame(.data)
  check_number_positive(from, allow_zero = TRUE)
  check_number_positive(to)
  if (to <= from) {
    cli::cli_abort("{.arg to} ({to}) must be greater than {.arg from} ({from}).")
  }
  cols <- profile_cols(.data, rlang::enquo(id), rlang::enquo(top), rlang::enquo(bottom))
  stock_col <- select_one(rlang::enquo(stock), .data, "stock")
  s <- .data[[stock_col]]
  if (!is.numeric(s)) {
    cli::cli_abort("{.arg stock} must be a numeric column.")
  }
  t <- .data[[cols$top]]
  b <- .data[[cols$bottom]]
  ov <- layer_overlap(t, b, from, to)
  share <- ifelse(b > t, ov / (b - t), 0)
  ids <- unique(.data[[cols$id]])
  idx <- split(seq_len(nrow(.data)), factor(.data[[cols$id]], levels = ids))
  total <- vapply(idx, function(r) {
    ok <- !is.na(s[r]) & !is.na(share[r]) & share[r] > 0
    if (!any(ok)) NA_real_ else sum(s[r][ok] * share[r][ok])
  }, numeric(1))
  cover <- vapply(idx, function(r) {
    ok <- !is.na(s[r]) & !is.na(ov[r])
    sum(ov[r][ok]) / (to - from)
  }, numeric(1))
  out <- tibble::tibble(!!cols$id := ids, from = from, to = to)
  out[[stock_col]] <- unname(total)
  out$coverage <- unname(cover)
  out
}
