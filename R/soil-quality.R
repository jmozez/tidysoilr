# Soil quality indices ------------------------------------------------------

new_sqi_scoring <- function(fn, type, limits) {
  structure(list(fn = fn, type = type, limits = limits), class = "sqi_scoring")
}

#' Scoring functions for soil quality indicators
#'
#' Linear scoring functions that map an indicator onto a 0 to 1 score for
#' [soil_quality_index()] (Andrews et al. 2002):
#'
#' * `sqi_more()`: "more is better" (for example organic carbon or
#'   aggregate stability): 0 at or below `lower`, 1 at or above `upper`.
#' * `sqi_less()`: "less is better" (for example bulk density or
#'   exchangeable sodium): 1 at or below `lower`, 0 at or above `upper`.
#' * `sqi_optimum()`: "optimum range" (for example pH): 0 outside
#'   `lower`-`upper`, 1 between `low_opt` and `high_opt`, linear between.
#'
#' @param lower,upper Limits of the scoring function.
#' @param low_opt,high_opt Limits of the optimum range for `sqi_optimum()`.
#' @return An object of class `sqi_scoring`.
#' @references Andrews, S. S., Karlen, D. L. and Mitchell, J. P. (2002). A
#'   comparison of soil quality indexing methods for vegetable production
#'   systems in Northern California. *Agriculture, Ecosystems and
#'   Environment*, 90, 25-45. \doi{10.1016/S0167-8809(01)00174-8}
#' @family soil quality
#' @export
#' @examples
#' sqi_more(5, 30)
#' sqi_less(1.1, 1.6)
#' sqi_optimum(4.5, 5.5, 7, 8.5)
sqi_more <- function(lower, upper) {
  check_limits(c(lower, upper))
  new_sqi_scoring(function(x) pmin(pmax((x - lower) / (upper - lower), 0), 1), "more is better",
    c(lower = lower, upper = upper))
}

#' @rdname sqi_more
#' @export
sqi_less <- function(lower, upper) {
  check_limits(c(lower, upper))
  new_sqi_scoring(function(x) pmin(pmax((upper - x) / (upper - lower), 0), 1), "less is better",
    c(lower = lower, upper = upper))
}

#' @rdname sqi_more
#' @export
sqi_optimum <- function(lower, low_opt, high_opt, upper) {
  check_limits(c(lower, low_opt, high_opt, upper), strict = FALSE)
  fn <- function(x) {
    up <- if (low_opt > lower) (x - lower) / (low_opt - lower) else as.numeric(x >= lower)
    down <- if (upper > high_opt) (upper - x) / (upper - high_opt) else as.numeric(x <= upper)
    pmin(pmax(pmin(up, down, 1), 0), 1)
  }
  new_sqi_scoring(fn, "optimum",
    c(lower = lower, low_opt = low_opt, high_opt = high_opt, upper = upper))
}

check_limits <- function(x, strict = TRUE, call = rlang::caller_env()) {
  ok <- is.numeric(x) && !anyNA(x) && all(is.finite(x)) &&
    (if (strict) all(diff(x) > 0) else all(diff(x) >= 0) && x[length(x)] > x[1])
  if (!ok) {
    cli::cli_abort("Scoring limits must be finite numbers in increasing order.", call = call)
  }
}

#' @export
print.sqi_scoring <- function(x, ...) {
  lim <- paste(names(x$limits), format(x$limits), sep = " = ", collapse = ", ")
  cat("<sqi_scoring> ", x$type, ": ", lim, "\n", sep = "")
  invisible(x)
}

#' Additive soil quality index
#'
#' Scores each indicator with a scoring function ([sqi_more()],
#' [sqi_less()] or [sqi_optimum()]) and combines the scores into a soil
#' quality index (SQI) between 0 and 1 as their (weighted) mean, following
#' the additive indexing approach of Andrews et al. (2002). Weights can come
#' from [soil_mds()].
#'
#' @param .data A data frame.
#' @param ... Name-value pairs: the name is an indicator column and the
#'   value its scoring function, for example `oc = sqi_more(5, 30)`.
#' @param weights Optional named numeric vector of weights, one per
#'   indicator. Equal weights by default.
#' @param na_rm Compute the index from the available scores when some are
#'   missing?
#' @param .name Name of the index column.
#'
#' @return `.data` with one `<indicator>_score` column per indicator and the
#'   index column.
#' @inherit sqi_more references
#' @family soil quality
#' @export
#' @examples
#' soil_quality_index(soil_quality_trial,
#'   oc = sqi_more(10, 40),
#'   bulk_density = sqi_less(1.0, 1.6),
#'   ph = sqi_optimum(4.5, 5.5, 7, 8.5),
#'   agg_stability = sqi_more(20, 80)
#' ) |>
#'   dplyr::group_by(land_use) |>
#'   dplyr::summarise(sqi = mean(sqi))
soil_quality_index <- function(.data, ..., weights = NULL, na_rm = FALSE, .name = "sqi") {
  check_data_frame(.data)
  check_name(.name)
  scoring <- rlang::list2(...)
  if (!length(scoring) || !rlang::is_named(scoring)) {
    cli::cli_abort("Supply indicators as named scoring functions, for example {.code oc = sqi_more(5, 30)}.")
  }
  bad <- !vapply(scoring, inherits, logical(1), "sqi_scoring")
  if (any(bad)) {
    cli::cli_abort("{.arg {names(scoring)[bad]}} must be scoring functions such as {.fn sqi_more}.")
  }
  miss <- setdiff(names(scoring), names(.data))
  if (length(miss)) {
    cli::cli_abort("Column{?s} {.field {miss}} not found in {.arg .data}.")
  }
  w <- weights %||% stats::setNames(rep(1, length(scoring)), names(scoring))
  if (!is.numeric(w) || !rlang::is_named(w) || !setequal(names(w), names(scoring)) ||
    any(w < 0)) {
    cli::cli_abort("{.arg weights} must be non-negative and named after the indicators.")
  }
  w <- w[names(scoring)]
  scores <- vapply(names(scoring), function(nm) {
    x <- .data[[nm]]
    if (!is.numeric(x)) cli::cli_abort("Indicator {.field {nm}} must be numeric.")
    scoring[[nm]]$fn(x)
  }, numeric(nrow(.data)))
  scores <- matrix(scores, nrow = nrow(.data), dimnames = list(NULL, names(scoring)))
  for (nm in names(scoring)) .data[[paste0(nm, "_score")]] <- scores[, nm]
  wm <- matrix(w, nrow(scores), ncol(scores), byrow = TRUE)
  wm[is.na(scores)] <- if (na_rm) 0 else NA
  s0 <- scores
  s0[is.na(s0)] <- 0
  .data[[.name]] <- rowSums(s0 * wm) / rowSums(wm)
  .data
}

#' Minimum data set of soil quality indicators by principal components
#'
#' Selects a minimum data set (MDS) of indicators with the principal
#' component procedure of Andrews et al. (2002): indicators are
#' standardised, principal components with eigenvalues of at least
#' `eigen_min` are kept, and in each component the indicators with absolute
#' loadings within `within` (10 percent) of the largest are highly
#' weighted. When several highly weighted indicators in a component are
#' correlated (|r| > `max_cor`), only the one with the largest loading is
#' kept. Each selected indicator is weighted by the share of variance
#' explained by its component, ready for [soil_quality_index()].
#'
#' @param .data A data frame.
#' @param ... <[`tidy-select`][tidyselect::language]> Numeric indicator
#'   columns.
#' @param eigen_min Minimum eigenvalue of retained components.
#' @param within Relative tolerance for highly weighted loadings.
#' @param max_cor Correlation above which highly weighted indicators are
#'   treated as redundant.
#'
#' @return A tibble with one row per retained component and highly weighted
#'   indicator: `component`, `eigenvalue`, `variance` (proportion), `indicator`,
#'   `loading`, `selected` and `weight` (for selected indicators).
#' @inherit sqi_more references
#' @family soil quality
#' @export
#' @examples
#' mds <- soil_mds(soil_quality_trial, oc:agg_stability)
#' mds
#' w <- with(mds[mds$selected, ], tapply(weight, indicator, sum))
#' w
soil_mds <- function(.data, ..., eigen_min = 1, within = 0.1, max_cor = 0.6) {
  check_data_frame(.data)
  check_number_positive(eigen_min, allow_zero = TRUE)
  check_proportion(within)
  check_proportion(max_cor)
  vars <- select_dots(.data, ...)
  if (length(vars) < 2) {
    cli::cli_abort("Select at least two indicator columns in {.arg ...}.")
  }
  x <- as.data.frame(.data)[vars]
  num <- vapply(x, is.numeric, logical(1))
  if (!all(num)) {
    cli::cli_abort("Indicators must be numeric; {.field {vars[!num]}} {?is/are} not.")
  }
  x <- stats::na.omit(x)
  if (nrow(x) <= length(vars)) {
    cli::cli_abort("Need more complete rows ({nrow(x)}) than indicators ({length(vars)}).")
  }
  pca <- stats::prcomp(x, center = TRUE, scale. = TRUE)
  ev <- pca$sdev^2
  keep <- which(ev >= eigen_min)
  if (!length(keep)) keep <- 1L
  prop <- ev / sum(ev)
  r <- stats::cor(x)
  res <- lapply(keep, function(k) {
    ld <- pca$rotation[, k]
    hw <- names(ld)[abs(ld) >= (1 - within) * max(abs(ld))]
    hw <- hw[order(-abs(ld[hw]))]
    sel <- hw[1]
    for (v in hw[-1]) {
      if (all(abs(r[v, sel]) <= max_cor)) sel <- c(sel, v)
    }
    tibble::tibble(component = paste0("PC", k), eigenvalue = ev[k], variance = prop[k],
      indicator = hw, loading = unname(ld[hw]), selected = hw %in% sel)
  })
  out <- dplyr::bind_rows(res)
  tot <- sum(prop[keep])
  n_sel <- stats::ave(as.numeric(out$selected), out$component, FUN = sum)
  out$weight <- ifelse(out$selected, out$variance / tot / n_sel, NA_real_)
  out
}
