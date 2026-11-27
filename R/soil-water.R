# Water retention: van Genuchten model ------------------------------------

vg_theta <- function(h, theta_r, theta_s, alpha, n) {
  m <- 1 - 1 / n
  theta_r + (theta_s - theta_r) / (1 + (alpha * abs(h))^n)^m
}

#' van Genuchten water retention function
#'
#' Evaluates the van Genuchten (1980) water retention function
#' \deqn{\theta(h) = \theta_r + \frac{\theta_s - \theta_r}{[1 + (\alpha h)^n]^{1 - 1/n}}}
#' for every row of `.data`. Each argument can be a column of `.data` or a
#' single number.
#'
#' @param .data A data frame.
#' @param h <[`data-masking`][rlang::args_data_masking]> Matric suction
#'   (positive, for example in cm of water; `alpha` must use the inverse
#'   unit).
#' @param theta_r,theta_s <[`data-masking`][rlang::args_data_masking]>
#'   Residual and saturated water content (m3/m3).
#' @param alpha,n <[`data-masking`][rlang::args_data_masking]> Shape
#'   parameters (`alpha` > 0, `n` > 1).
#' @param .name Name of the new column.
#'
#' @return `.data` with a new column of volumetric water content.
#' @references van Genuchten, M. Th. (1980). A closed-form equation for
#'   predicting the hydraulic conductivity of unsaturated soils. *Soil
#'   Science Society of America Journal*, 44, 892-898.
#'   \doi{10.2136/sssaj1980.03615995004400050002x}
#' @family water retention
#' @export
#' @examples
#' df <- data.frame(h = c(0, 10, 100, 330, 1000, 15000))
#' soil_vg_theta(df, h, theta_r = 0.078, theta_s = 0.43, alpha = 0.036, n = 1.56)
soil_vg_theta <- function(.data, h, theta_r, theta_s, alpha, n, .name = "theta_vg") {
  check_data_frame(.data)
  check_name(.name)
  ev <- function(q, arg) {
    x <- rlang::eval_tidy(q, .data)
    if (!is.numeric(x)) {
      cli::cli_abort("{.arg {arg}} must be numeric.", call = rlang::caller_env(2))
    }
    x
  }
  hh <- ev(rlang::enquo(h), "h")
  a <- ev(rlang::enquo(alpha), "alpha")
  nn <- ev(rlang::enquo(n), "n")
  if (any(nn <= 1, na.rm = TRUE)) {
    cli::cli_abort("{.arg n} must be greater than 1.")
  }
  .data[[.name]] <- vg_theta(hh, ev(rlang::enquo(theta_r), "theta_r"),
    ev(rlang::enquo(theta_s), "theta_s"), a, nn)
  .data
}

# Fit one curve by bounded least squares with several starting values.
vg_fit_one <- function(h, theta, fix_theta_r = NULL) {
  ok <- !is.na(h) & !is.na(theta)
  h <- abs(h[ok])
  theta <- theta[ok]
  n_par <- if (is.null(fix_theta_r)) 4 else 3
  if (length(theta) < n_par + 1) {
    return(NULL)
  }
  tmax <- max(theta)
  tmin <- min(theta)
  unpack <- function(p) {
    if (is.null(fix_theta_r)) {
      list(theta_r = p[1], theta_s = p[2], alpha = 10^p[3], n = 1 + exp(p[4]))
    } else {
      list(theta_r = fix_theta_r, theta_s = p[1], alpha = 10^p[2], n = 1 + exp(p[3]))
    }
  }
  sse <- function(p) {
    q <- unpack(p)
    sum((theta - vg_theta(h, q$theta_r, q$theta_s, q$alpha, q$n))^2)
  }
  pos <- h[h > 0]
  mid <- if (length(pos)) stats::median(pos) else 100
  alpha_starts <- unique(log10(c(1 / mid, 0.002, 0.02, 0.2)))
  best <- NULL
  for (a0 in alpha_starts) {
    for (n0 in log(c(0.2, 0.6, 1.5))) {
      if (is.null(fix_theta_r)) {
        start <- c(tmin * 0.5, tmax, a0, n0)
        lower <- c(0, tmax * 0.8, -6, log(1e-3))
        upper <- c(tmin, 1, 3, log(15))
      } else {
        start <- c(tmax, a0, n0)
        lower <- c(max(fix_theta_r, tmax * 0.8), -6, log(1e-3))
        upper <- c(1, 3, log(15))
      }
      res <- tryCatch(
        stats::optim(start, sse, method = "L-BFGS-B", lower = lower, upper = upper),
        error = function(e) NULL
      )
      if (!is.null(res) && (is.null(best) || res$value < best$value)) best <- res
    }
  }
  if (is.null(best)) {
    return(NULL)
  }
  q <- unpack(best$par)
  fitted <- vg_theta(h, q$theta_r, q$theta_s, q$alpha, q$n)
  res <- theta - fitted
  tibble::tibble(
    theta_r = q$theta_r, theta_s = q$theta_s, alpha = q$alpha, n = q$n, m = 1 - 1 / q$n,
    n_obs = length(theta), rmse = sqrt(mean(res^2)),
    r_squared = 1 - sum(res^2) / sum((theta - mean(theta))^2),
    converged = best$convergence == 0
  )
}

#' Fit van Genuchten water retention curves
#'
#' Fits the van Genuchten (1980) model with the Mualem restriction
#' \eqn{m = 1 - 1/n} to measured water retention data by bounded nonlinear
#' least squares ([stats::optim()], method `"L-BFGS-B"`) from several
#' starting values, one curve per group. This is the base R equivalent of
#' the RETC program, with no compiled dependencies.
#'
#' @param .data A data frame in long format: one row per measured point.
#' @param h <[`tidy-select`][tidyselect::language]> Matric suction (positive;
#'   cm of water is conventional and gives `alpha` in 1/cm).
#' @param theta <[`tidy-select`][tidyselect::language]> Volumetric water
#'   content (m3/m3).
#' @param .by <[`tidy-select`][tidyselect::language]> Columns that identify
#'   each curve (for example sample identifier). Groups of a grouped data
#'   frame are used when `.by` is not given.
#' @param fix_theta_r Optional fixed residual water content (for example 0)
#'   for curves with few dry-end points.
#'
#' @return A `<soil_retention>` tibble with one row per curve: the group
#'   columns, `theta_r`, `theta_s`, `alpha`, `n`, `m`, `n_obs`, `rmse`,
#'   `r_squared` and `converged`. The data are kept for [augment()],
#'   [predict()] and [autoplot()].
#' @inherit soil_vg_theta references
#' @family water retention
#' @export
#' @examples
#' fits <- soil_vg_fit(soil_retention_lab, h = h_cm, theta = theta, .by = sample)
#' fits
#' predict(fits, h = c(100, 330, 15000))
soil_vg_fit <- function(.data, h, theta, .by = NULL, fix_theta_r = NULL) {
  check_data_frame(.data)
  check_number_positive(fix_theta_r, allow_zero = TRUE, allow_null = TRUE)
  h_col <- select_one(rlang::enquo(h), .data, "h")
  theta_col <- select_one(rlang::enquo(theta), .data, "theta")
  numeric_col(rlang::quo(!!rlang::sym(h_col)), .data, "h")
  th <- numeric_col(rlang::quo(!!rlang::sym(theta_col)), .data, "theta", non_negative = TRUE)
  if (any(th > 1, na.rm = TRUE)) {
    cli::cli_abort("{.arg theta} must be a volumetric fraction (0 to 1); divide percentages by 100.")
  }
  by_cols <- resolve_by(rlang::enquo(.by), .data)
  data <- as_plain_tibble(dplyr::ungroup(.data))
  fits <- by_apply(data, by_cols, function(d) {
    f <- vg_fit_one(d[[h_col]], d[[theta_col]], fix_theta_r)
    if (is.null(f)) {
      tibble::tibble(theta_r = NA_real_, theta_s = NA_real_, alpha = NA_real_, n = NA_real_,
        m = NA_real_, n_obs = sum(!is.na(d[[theta_col]])), rmse = NA_real_,
        r_squared = NA_real_, converged = FALSE)
    } else {
      f
    }
  })
  failed <- sum(!fits$converged | is.na(fits$theta_r))
  if (failed) {
    cli::cli_warn("{failed} curve{?s} did not converge or had too few points.")
  }
  meta <- list(groups = by_cols, h = h_col, theta = theta_col,
    data = data[c(by_cols, h_col, theta_col)])
  new_soil_retention(fits, meta)
}

#' Predict water content from fitted retention curves
#'
#' @param object A `<soil_retention>` from [soil_vg_fit()].
#' @param new_data Optional data frame with the group columns and a suction
#'   column named as in the fit. If `NULL`, every curve is evaluated at `h`.
#' @param h Suctions used when `new_data` is `NULL`.
#' @param ... Not used.
#' @return A tibble with the group columns, the suction and `.pred`.
#' @family water retention
#' @export
#' @examples
#' fits <- soil_vg_fit(soil_retention_lab, h_cm, theta, .by = sample)
#' predict(fits, h = c(10, 330, 15000))
predict.soil_retention <- function(object, new_data = NULL, h = 10^seq(0, 4.2, by = 0.1), ...) {
  meta <- soil_meta(object)
  fits <- as_plain_tibble(object)
  groups <- meta$groups
  if (is.null(new_data)) {
    idx <- rep(seq_len(nrow(fits)), each = length(h))
    new_data <- fits[idx, groups, drop = FALSE]
    new_data[[meta$h]] <- rep(h, nrow(fits))
  } else {
    check_data_frame(new_data)
    miss <- setdiff(c(groups, meta$h), names(new_data))
    if (length(miss)) {
      cli::cli_abort("{.arg new_data} is missing column{?s} {.field {miss}}.")
    }
    new_data <- as_plain_tibble(new_data)
  }
  if (length(groups)) {
    joined <- dplyr::left_join(new_data, fits[c(groups, "theta_r", "theta_s", "alpha", "n")],
      by = groups)
  } else {
    joined <- new_data
    for (p in c("theta_r", "theta_s", "alpha", "n")) joined[[p]] <- fits[[p]][1]
  }
  new_data$.pred <- vg_theta(joined[[meta$h]], joined$theta_r, joined$theta_s, joined$alpha,
    joined$n)
  new_data
}

#' Available water capacity from fitted retention curves
#'
#' Evaluates each fitted curve at field capacity and permanent wilting
#' point and returns their difference, the available water capacity
#' (AWC). When a layer `thickness` (cm) is given, the plant-available water
#' held by the layer is also returned in millimetres.
#'
#' @param .data A `<soil_retention>` from [soil_vg_fit()].
#' @param fc Suction at field capacity, in the units of the fit (330 cm,
#'   about 33 kPa, by default; use 100 cm for sandy soils).
#' @param pwp Suction at permanent wilting point (15000 cm, about 1500 kPa).
#' @param thickness Optional layer thickness in cm.
#'
#' @return `.data` with `theta_fc`, `theta_pwp`, `awc` (m3/m3) and, with
#'   `thickness`, `awc_mm`.
#' @family water retention
#' @export
#' @examples
#' fits <- soil_vg_fit(soil_retention_lab, h_cm, theta, .by = sample)
#' soil_awc(fits, thickness = 30)
soil_awc <- function(.data, fc = 330, pwp = 15000, thickness = NULL) {
  if (!is_soil_retention(.data)) {
    cli::cli_abort("{.arg .data} must be a {.cls soil_retention} from {.fn soil_vg_fit}.")
  }
  check_number_positive(fc)
  check_number_positive(pwp)
  check_number_positive(thickness, allow_null = TRUE)
  if (pwp <= fc) {
    cli::cli_abort("{.arg pwp} must be a higher suction than {.arg fc}.")
  }
  tfc <- vg_theta(fc, .data$theta_r, .data$theta_s, .data$alpha, .data$n)
  tpwp <- vg_theta(pwp, .data$theta_r, .data$theta_s, .data$alpha, .data$n)
  out <- list(theta_fc = tfc, theta_pwp = tpwp, awc = tfc - tpwp)
  if (!is.null(thickness)) out$awc_mm <- (tfc - tpwp) * thickness * 10
  add_cols(.data, out)
}
