# Unsaturated hydraulic conductivity and infiltration -----------------------

mvg_k <- function(h, ks, alpha, n, l = 0.5) {
  m <- 1 - 1 / n
  se <- (1 + (alpha * abs(h))^n)^(-m)
  ks * se^l * (1 - (1 - se^(1 / m))^m)^2
}

#' Mualem-van Genuchten unsaturated hydraulic conductivity
#'
#' Predicts the unsaturated hydraulic conductivity from the van Genuchten
#' retention parameters with the Mualem (1976) pore-size model, as in
#' van Genuchten (1980), HYDRUS and Rosetta:
#' \deqn{K(h) = K_s S_e^{L} \left[1 - \left(1 - S_e^{1/m}\right)^m\right]^2,
#' \quad S_e = \left[1 + (\alpha h)^n\right]^{-m}, \quad m = 1 - 1/n.}{K(h)
#' = Ks Se^L [1 - (1 - Se^(1/m))^m]^2, Se = [1 + (alpha h)^n]^(-m), m = 1 -
#' 1/n.}
#' Each argument can be a column of `.data` or a single number, so the
#' parameters from [soil_vg_fit()] can be joined and used directly.
#'
#' @param .data A data frame.
#' @param h <[`data-masking`][rlang::args_data_masking]> Matric suction
#'   (positive, in the inverse unit of `alpha`).
#' @param ks <[`data-masking`][rlang::args_data_masking]> Saturated
#'   hydraulic conductivity (any unit; the result has the same unit).
#' @param alpha,n <[`data-masking`][rlang::args_data_masking]> van Genuchten
#'   shape parameters (`n` > 1).
#' @param l <[`data-masking`][rlang::args_data_masking]> Pore connectivity
#'   (tortuosity) parameter; 0.5 by default (Mualem 1976).
#' @param .name Name of the new column.
#'
#' @return `.data` with a new column of hydraulic conductivity.
#' @references Mualem, Y. (1976). A new model for predicting the hydraulic
#'   conductivity of unsaturated porous media. Water Resources Research, 12,
#'   513-522. \doi{10.1029/WR012i003p00513}
#' @inherit soil_vg_theta references
#' @family water retention
#' @export
#' @examples
#' # Loam parameters from Carsel and Parrish (1988), Ks in cm/day
#' df <- data.frame(h = c(0, 10, 100, 330, 1000, 15000))
#' soil_vg_conductivity(df, h, ks = 24.96, alpha = 0.036, n = 1.56)
soil_vg_conductivity <- function(.data, h, ks, alpha, n, l = 0.5, .name = "k_mvg") {
  check_data_frame(.data)
  check_name(.name)
  ev <- function(q, arg) {
    x <- rlang::eval_tidy(q, .data)
    if (!is.numeric(x)) {
      cli::cli_abort("{.arg {arg}} must be numeric.", call = rlang::caller_env(2))
    }
    x
  }
  nn <- ev(rlang::enquo(n), "n")
  if (any(nn <= 1, na.rm = TRUE)) {
    cli::cli_abort("{.arg n} must be greater than 1.")
  }
  k <- ev(rlang::enquo(ks), "ks")
  if (any(k < 0, na.rm = TRUE)) {
    cli::cli_abort("{.arg ks} must not be negative.")
  }
  .data[[.name]] <- mvg_k(ev(rlang::enquo(h), "h"), k, ev(rlang::enquo(alpha), "alpha"), nn,
    ev(rlang::enquo(l), "l"))
  .data
}

infiltration_models <- c("philip", "horton", "kostiakov")

infil_value <- function(model, p, t) {
  switch(model,
    philip = p[["sorptivity"]] * sqrt(t) + p[["transmissivity"]] * t,
    horton = p[["fc"]] * t + (p[["f0"]] - p[["fc"]]) / p[["k"]] * (1 - exp(-p[["k"]] * t)),
    kostiakov = p[["a"]] * t^p[["b"]]
  )
}

infil_rate <- function(model, p, t) {
  switch(model,
    philip = 0.5 * p[["sorptivity"]] / sqrt(t) + p[["transmissivity"]],
    horton = p[["fc"]] + (p[["f0"]] - p[["fc"]]) * exp(-p[["k"]] * t),
    kostiakov = p[["a"]] * p[["b"]] * t^(p[["b"]] - 1)
  )
}

infil_fit_one <- function(t, y, model) {
  pos <- t > 0
  if (model == "philip") {
    fit <- stats::lm.fit(cbind(sqrt(t), t), y)
    return(c(sorptivity = fit$coefficients[[1]], transmissivity = fit$coefficients[[2]]))
  }
  if (model == "kostiakov") {
    ok <- pos & y > 0
    st <- stats::lm.fit(cbind(1, log(t[ok])), log(y[ok]))$coefficients
    rss <- function(q) sum((y - exp(q[[1]]) * t^q[[2]])^2)
    res <- stats::optim(st, rss, method = "BFGS", control = list(maxit = 1000, reltol = 1e-14))
    return(c(a = exp(res$par[[1]]), b = res$par[[2]]))
  }
  # Horton, with fc >= 0, f0 > fc and k > 0 through log parameters.
  ord <- order(t)
  slopes <- diff(y[ord]) / diff(t[ord])
  slopes <- slopes[is.finite(slopes)]
  fc0 <- max(utils::tail(slopes, 1), 1e-6)
  f00 <- max(slopes[[1]], fc0 * 2)
  to_p <- function(q) c(f0 = exp(q[[1]]) + exp(q[[2]]), fc = exp(q[[1]]), k = exp(q[[3]]))
  rss <- function(q) {
    s <- sum((y - infil_value("horton", to_p(q), t))^2)
    if (is.finite(s)) s else 1e100
  }
  best <- NULL
  for (kk in c(0.5, 2, 8) / max(t)) {
    q0 <- c(log(fc0), log(f00 - fc0), log(kk))
    res <- stats::optim(q0, rss, method = "Nelder-Mead", control = list(maxit = 4000, reltol = 1e-12))
    res <- stats::optim(res$par, rss, method = "BFGS", control = list(maxit = 1000, reltol = 1e-14))
    if (is.null(best) || res$value < best$value) best <- res
  }
  to_p(best$par)
}

#' Fit infiltration models to cumulative infiltration
#'
#' Fits the classical infiltration equations to cumulative infiltration
#' measured with a ring or disc infiltrometer, by least squares:
#'
#' * Philip (1957): \eqn{I = S t^{1/2} + A t}, with sorptivity \eqn{S} and
#'   the transmissivity (gravity) term \eqn{A}, related to the saturated
#'   conductivity;
#' * Horton (1940): \eqn{I = f_c t + (f_0 - f_c)(1 - e^{-k t}) / k}, with
#'   initial and final infiltration rates \eqn{f_0} and \eqn{f_c};
#' * Kostiakov (1932): \eqn{I = a t^b}.
#'
#' The models are compared with the root mean square error and AIC for
#' least squares. The steady infiltration rate is \eqn{A} for Philip and
#' \eqn{f_c} for Horton; the infiltration rate at the last time is given
#' for all models.
#'
#' @param .data A data frame with one row per reading.
#' @param time,infiltration <[`tidy-select`][tidyselect::language]> Elapsed
#'   time and cumulative infiltration (for example minutes and mm, or hours
#'   and cm).
#' @param .by <[`tidy-select`][tidyselect::language]> Columns identifying
#'   each test. Groups of a grouped data frame are used when `.by` is not
#'   given.
#' @param models Models to fit, any of `"philip"`, `"horton"` and
#'   `"kostiakov"`.
#'
#' @return A tibble with one row per test and model: the group columns,
#'   `model`, the parameters (`sorptivity` and `transmissivity` for
#'   Philip, `f0`, `fc` and `k` for Horton, `a` and `b` for Kostiakov;
#'   missing where they do not apply), `steady_rate`, `final_rate`,
#'   `n_obs`, `rmse`, `r_squared` and `aic`.
#' @references Philip, J. R. (1957). The theory of infiltration: 4.
#'   Sorptivity and algebraic infiltration equations. Soil Science, 84,
#'   257-264. \doi{10.1097/00010694-195709000-00010}
#'
#'   Horton, R. E. (1940). An approach toward a physical interpretation of
#'   infiltration capacity. Soil Science Society of America Proceedings, 5,
#'   399-417. \doi{10.2136/sssaj1941.036159950005000C0075x}
#'
#'   Kostiakov, A. N. (1932). On the dynamics of the coefficient of water
#'   percolation in soils. Transactions of the 6th Commission of the
#'   International Society of Soil Science, Part A, 17-21.
#' @family water retention
#' @export
#' @examples
#' # Double-ring readings (minutes, mm) at two sites
#' rings <- data.frame(
#'   site = rep(c("garden", "pasture"), each = 9),
#'   minutes = rep(c(1, 2, 5, 10, 20, 30, 45, 60, 90), 2),
#'   mm = c(4.1, 6.2, 10.4, 15.6, 23.9, 31.0, 40.6, 49.7, 67.2,
#'          1.9, 2.8, 4.7, 6.8, 9.9, 12.6, 16.1, 19.4, 25.6)
#' )
#' soil_infiltration_fit(rings, minutes, mm, .by = site)
soil_infiltration_fit <- function(.data, time, infiltration, .by = NULL,
                                  models = c("philip", "horton", "kostiakov")) {
  check_data_frame(.data)
  models <- rlang::arg_match(models, infiltration_models, multiple = TRUE)
  t_col <- select_one(rlang::enquo(time), .data, "time")
  i_col <- select_one(rlang::enquo(infiltration), .data, "infiltration")
  numeric_col(rlang::quo(!!rlang::sym(t_col)), .data, "time", non_negative = TRUE)
  numeric_col(rlang::quo(!!rlang::sym(i_col)), .data, "infiltration", non_negative = TRUE)
  by_cols <- resolve_by(rlang::enquo(.by), .data)
  by_apply(.data, by_cols, function(d) {
    ok <- !is.na(d[[t_col]]) & !is.na(d[[i_col]])
    t <- d[[t_col]][ok]
    y <- d[[i_col]][ok]
    if (length(unique(t[t > 0])) < 4L) {
      cli::cli_abort("Each test needs at least 4 readings after time zero.", call = NULL)
    }
    rows <- lapply(models, function(m) {
      p <- infil_fit_one(t, y, m)
      fitted <- infil_value(m, p, t)
      rss <- sum((y - fitted)^2)
      n <- length(y)
      get <- function(nm) if (nm %in% names(p)) unname(p[[nm]]) else NA_real_
      tibble::tibble(
        model = m, sorptivity = get("sorptivity"), transmissivity = get("transmissivity"),
        a = get("a"), b = get("b"),
        f0 = get("f0"), fc = get("fc"), k = get("k"),
        steady_rate = switch(m, philip = get("transmissivity"), horton = get("fc"), NA_real_),
        final_rate = infil_rate(m, p, max(t)),
        n_obs = n, rmse = sqrt(rss / n), r_squared = 1 - rss / sum((y - mean(y))^2),
        aic = n * log(rss / n) + 2 * (length(p) + 1)
      )
    })
    dplyr::bind_rows(rows)
  })
}
