# yardstick metrics -----------------------------------------------------------

#' Standard error of prediction (bias-corrected)
#'
#' The standard error of prediction (SEP) is the standard deviation of the
#' prediction errors after removing their mean (the bias):
#' \deqn{SEP = \sqrt{\frac{\sum (e_i - \bar e)^2}{n - 1}}, \quad e_i = \hat y_i - y_i.}
#' It is reported alongside RMSE and bias in soil spectroscopy, because
#' \eqn{RMSE^2 \approx SEP^2 + bias^2} separates random error from
#' systematic error (Bellon-Maurel et al. 2010). RPIQ, RPD and Lin's
#' concordance correlation coefficient are already provided by yardstick
#' as [yardstick::rpiq()], [yardstick::rpd()] and [yardstick::ccc()], and
#' combine with `sep()` in a [yardstick::metric_set()].
#'
#' @param data A data frame with the `truth` and `estimate` columns.
#' @param truth The column of observed values.
#' @param estimate The column of predicted values.
#' @param na_rm Remove missing values?
#' @param case_weights Optional column of case weights.
#' @param ... Not used.
#'
#' @return For `sep()`, a tibble with columns `.metric`, `.estimator` and
#'   `.estimate`. For `sep_vec()`, a single number.
#' @references Bellon-Maurel, V., Fernandez-Ahumada, E., Palagos, B., Roger,
#'   J.-M. and McBratney, A. (2010). Critical review of chemometric
#'   indicators commonly used for assessing the quality of the prediction of
#'   soil attributes by NIR spectroscopy. *TrAC Trends in Analytical
#'   Chemistry*, 29, 1073-1081. \doi{10.1016/j.trac.2010.05.006}
#' @family soil metrics
#' @export
#' @examples
#' df <- data.frame(obs = c(10, 20, 30, 40, 50), pred = c(13, 22, 31, 45, 52))
#' sep(df, obs, pred)
#' sep_vec(df$obs, df$pred)
#'
#' # Combine with yardstick's RPIQ and Lin's CCC
#' soil_metrics <- yardstick::metric_set(yardstick::rmse, sep, yardstick::rpiq, yardstick::ccc)
#' soil_metrics(df, obs, pred)
sep <- function(data, ...) {
  UseMethod("sep")
}
sep <- yardstick::new_numeric_metric(sep, direction = "minimize")

#' @rdname sep
#' @export
sep.data.frame <- function(data, truth, estimate, na_rm = TRUE, case_weights = NULL, ...) {
  yardstick::numeric_metric_summarizer(
    name = "sep", fn = sep_vec, data = data,
    truth = !!rlang::enquo(truth), estimate = !!rlang::enquo(estimate),
    na_rm = na_rm, case_weights = !!rlang::enquo(case_weights)
  )
}

#' @rdname sep
#' @export
sep_vec <- function(truth, estimate, na_rm = TRUE, case_weights = NULL, ...) {
  yardstick::check_numeric_metric(truth, estimate, case_weights)
  if (na_rm) {
    res <- yardstick::yardstick_remove_missing(truth, estimate, case_weights)
    truth <- res$truth
    estimate <- res$estimate
    case_weights <- res$case_weights
  } else if (yardstick::yardstick_any_missing(truth, estimate, case_weights)) {
    return(NA_real_)
  }
  e <- estimate - truth
  if (length(e) < 2) {
    return(NA_real_)
  }
  if (is.null(case_weights)) {
    return(stats::sd(e))
  }
  w <- as.numeric(case_weights)
  mu <- sum(w * e) / sum(w)
  # Reliability-weighted variance with an n - 1 style correction.
  v <- sum(w * (e - mu)^2) / (sum(w) - sum(w^2) / sum(w))
  sqrt(v)
}

# Pairs of texture classes that share a boundary on the triangle, found by
# classifying a fine lattice and recording neighbouring points that differ.
texture_adjacency <- function(system = c("usda", "fao")) {
  system <- match.arg(system)
  step <- 0.25
  g <- expand.grid(sand = seq(0, 100, step), clay = seq(0, 100, step))
  g <- g[g$sand + g$clay <= 100 + 1e-9, ]
  fn <- if (system == "usda") classify_usda else classify_fao
  key <- paste(g$sand, g$clay)
  cls <- stats::setNames(as.character(fn(g$sand, 100 - g$sand - g$clay, g$clay)), key)
  pairs <- character()
  for (d in list(c(step, 0), c(0, step), c(step, -step))) {
    nb <- paste(g$sand + d[1], g$clay + d[2])
    ok <- nb %in% key
    a <- cls[key[ok]]
    b <- cls[nb[ok]]
    diff <- a != b
    pairs <- c(pairs, paste(pmin(a[diff], b[diff]), pmax(a[diff], b[diff]), sep = "|"))
  }
  counts <- table(pairs)
  # Classes touching at a single point give only a handful of pairs.
  names(counts)[counts >= 8]
}

adjacency_cache <- new.env(parent = emptyenv())

get_adjacency <- function(system) {
  if (is.null(adjacency_cache[[system]])) {
    adjacency_cache[[system]] <- texture_adjacency(system)
  }
  adjacency_cache[[system]]
}

#' Texture class accuracy allowing adjacent classes
#'
#' The proportion of samples whose predicted texture class is either the
#' true class or a class that shares a boundary with it on the texture
#' triangle. Texture classes are arbitrary divisions of a continuum, so a
#' prediction of "clay loam" for a "silty clay loam" sample is a smaller
#' error than a prediction of "sand"; this metric complements plain
#' [yardstick::accuracy()] for classified texture predictions.
#'
#' @param data A data frame with the `truth` and `estimate` columns.
#' @param truth,estimate Columns of true and predicted texture classes
#'   (factors with the same levels, for example from
#'   [soil_texture_class()]).
#' @param system `"usda"` or `"fao"`: which triangle defines adjacency.
#' @param estimator Not used; always `"multiclass"`.
#' @param na_rm Remove missing values?
#' @param case_weights Optional column of case weights.
#' @param event_level Not used.
#' @param ... Not used.
#'
#' @return For `texture_adjacent_accuracy()`, a tibble with `.metric`,
#'   `.estimator` and `.estimate`. For `texture_adjacent_accuracy_vec()`, a
#'   single number between 0 and 1.
#' @family soil metrics
#' @export
#' @examples
#' lv <- soil_texture_polygons("usda")$class |> levels()
#' df <- data.frame(
#'   truth = factor(c("clay", "silty clay loam", "loam", "sand"), levels = lv),
#'   pred = factor(c("clay", "clay loam", "sandy loam", "silt"), levels = lv)
#' )
#' texture_adjacent_accuracy(df, truth, pred)
#' yardstick::accuracy(df, truth, pred)
texture_adjacent_accuracy <- function(data, ...) {
  UseMethod("texture_adjacent_accuracy")
}
texture_adjacent_accuracy <- yardstick::new_class_metric(texture_adjacent_accuracy,
  direction = "maximize")

#' @rdname texture_adjacent_accuracy
#' @export
texture_adjacent_accuracy.data.frame <- function(data, truth, estimate, system = "usda",
                                                 estimator = NULL, na_rm = TRUE,
                                                 case_weights = NULL, event_level = NULL, ...) {
  yardstick::class_metric_summarizer(
    name = "texture_adjacent_accuracy", fn = texture_adjacent_accuracy_vec, data = data,
    truth = !!rlang::enquo(truth), estimate = !!rlang::enquo(estimate),
    estimator = estimator, na_rm = na_rm, case_weights = !!rlang::enquo(case_weights),
    fn_options = list(system = system)
  )
}

#' @rdname texture_adjacent_accuracy
#' @export
texture_adjacent_accuracy_vec <- function(truth, estimate, system = "usda", estimator = NULL,
                                          na_rm = TRUE, case_weights = NULL, ...) {
  system <- rlang::arg_match0(system, c("usda", "fao"))
  estimator <- yardstick::finalize_estimator(truth, estimator,
    metric_class = "texture_adjacent_accuracy")
  yardstick::check_class_metric(truth, estimate, case_weights, estimator)
  if (na_rm) {
    res <- yardstick::yardstick_remove_missing(truth, estimate, case_weights)
    truth <- res$truth
    estimate <- res$estimate
    case_weights <- res$case_weights
  } else if (yardstick::yardstick_any_missing(truth, estimate, case_weights)) {
    return(NA_real_)
  }
  valid <- if (system == "usda") usda_levels else fao_levels
  bad <- setdiff(levels(truth), valid)
  if (length(bad)) {
    cli::cli_abort("{.val {bad}} {?is not a/are not} {toupper(system)} texture class{?es}.")
  }
  a <- as.character(truth)
  b <- as.character(estimate)
  adj <- paste(pmin(a, b), pmax(a, b), sep = "|") %in% get_adjacency(system)
  hit <- as.numeric(a == b | adj)
  w <- if (is.null(case_weights)) rep(1, length(hit)) else as.numeric(case_weights)
  sum(w * hit) / sum(w)
}

#' @exportS3Method yardstick::finalize_estimator_internal
finalize_estimator_internal.texture_adjacent_accuracy <- function(metric_dispatcher, x, estimator,
                                                                  call = rlang::caller_env()) {
  "multiclass"
}
