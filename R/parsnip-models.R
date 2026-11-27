# parsnip model specification: partial least squares -------------------------

#' Partial least squares regression for soil spectroscopy
#'
#' `soil_pls()` defines a partial least squares (PLS) regression model for
#' the tidymodels framework, the standard calibration model in soil
#' spectroscopy. The `"nipals"` engine fits a single-response PLS model
#' with the NIPALS algorithm (Wold et al. 2001) in base R, so unlike
#' `parsnip::pls()` it needs no Bioconductor packages. Predictors are
#' centred, and scaled when the engine argument `scale = TRUE`. Tune
#' `num_comp` with resampling, for example [soil_profile_cv()] or
#' [soil_kennard_stone()].
#'
#' @param mode Must be `"regression"`.
#' @param engine `"nipals"`.
#' @param num_comp Number of latent components.
#'
#' @return A model specification of class `soil_pls`.
#' @references Wold, S., Sjostrom, M. and Eriksson, L. (2001). PLS-regression:
#'   a basic tool of chemometrics. *Chemometrics and Intelligent Laboratory
#'   Systems*, 58, 109-130. \doi{10.1016/S0169-7439(01)00155-1}
#' @family parsnip models
#' @export
#' @examples
#' nir <- soil_nir[c("clay", grep("^nm_", names(soil_nir), value = TRUE))]
#' spec <- soil_pls(num_comp = 5) |> parsnip::set_engine("nipals")
#' fit <- parsnip::fit(spec, clay ~ ., data = nir)
#' predict(fit, nir[1:5, ])
soil_pls <- function(mode = "regression", engine = "nipals", num_comp = NULL) {
  if (mode != "regression") {
    cli::cli_abort("{.fn soil_pls} only supports {.val regression} mode.")
  }
  args <- list(num_comp = rlang::enquo(num_comp))
  parsnip::new_model_spec(
    "soil_pls", args = args, eng_args = NULL, mode = mode, user_specified_mode = !missing(mode),
    method = NULL, engine = engine, user_specified_engine = !missing(engine)
  )
}

#' @export
print.soil_pls <- function(x, ...) {
  cat("Partial least squares regression model specification (", x$mode, ")\n\n", sep = "")
  parsnip::model_printer(x, ...)
  if (!is.null(x$method$fit$args)) {
    cat("Model fit template:\n")
    print(parsnip::show_call(x))
  }
  invisible(x)
}

#' @export
update.soil_pls <- function(object, parameters = NULL, num_comp = NULL, fresh = FALSE, ...) {
  args <- list(num_comp = rlang::enquo(num_comp))
  parsnip::update_spec(
    object = object, parameters = parameters, args_enquo_list = args, fresh = fresh,
    cls = "soil_pls", ...
  )
}

#' Fit engine for `soil_pls()`
#'
#' The fitting function called by parsnip for [soil_pls()]. It is exported
#' so that parsnip can call it and for use without tidymodels; most users
#' should use [soil_pls()].
#'
#' @param x A data frame or matrix of numeric predictors.
#' @param y A numeric response.
#' @param num_comp Number of latent components (reduced to the number the
#'   data support).
#' @param scale Scale predictors to unit variance?
#' @return An object of class `soil_pls_fit`.
#' @family parsnip models
#' @export
#' @examples
#' nir <- soil_nir[grep("^nm_", names(soil_nir))]
#' fit <- soil_pls_fit(nir, soil_nir$clay, num_comp = 4)
#' fit
#' head(predict(fit, nir))
soil_pls_fit <- function(x, y, num_comp = 10, scale = FALSE) {
  x <- as.matrix(as.data.frame(x))
  if (!is.numeric(x)) {
    cli::cli_abort("All predictors must be numeric.")
  }
  if (!is.numeric(y) || length(y) != nrow(x)) {
    cli::cli_abort("{.arg y} must be numeric with one value per row of {.arg x}.")
  }
  if (anyNA(x) || anyNA(y)) {
    cli::cli_abort("Predictors and response must not contain missing values.")
  }
  num_comp <- num_comp %||% 10
  check_integerish(num_comp)
  num_comp <- min(num_comp, nrow(x) - 1, ncol(x))
  center <- colMeans(x)
  scl <- if (isTRUE(scale)) apply(x, 2, stats::sd) else rep(1, ncol(x))
  scl[!is.finite(scl) | scl == 0] <- 1
  xm <- sweep(sweep(x, 2, center), 2, scl, "/")
  ym <- mean(y)
  yr <- y - ym
  w <- p <- matrix(0, ncol(x), num_comp)
  q <- numeric(num_comp)
  scores <- matrix(0, nrow(x), num_comp)
  xr <- xm
  for (a in seq_len(num_comp)) {
    wa <- crossprod(xr, yr)
    nw <- sqrt(sum(wa^2))
    if (nw < 1e-12) {
      num_comp <- a - 1
      break
    }
    wa <- wa / nw
    ta <- xr %*% wa
    tt <- sum(ta^2)
    pa <- crossprod(xr, ta) / tt
    qa <- sum(yr * ta) / tt
    xr <- xr - ta %*% t(pa)
    yr <- yr - qa * ta
    w[, a] <- wa
    p[, a] <- pa
    q[a] <- qa
    scores[, a] <- ta
  }
  if (num_comp < 1) {
    cli::cli_abort("The predictors carry no information about the response.")
  }
  keep <- seq_len(num_comp)
  w <- w[, keep, drop = FALSE]
  p <- p[, keep, drop = FALSE]
  q <- q[keep]
  scores <- scores[, keep, drop = FALSE]
  # Coefficients for every number of components, one column each.
  coefs <- vapply(keep, function(k) {
    wk <- w[, seq_len(k), drop = FALSE]
    pk <- p[, seq_len(k), drop = FALSE]
    as.vector(wk %*% solve(crossprod(pk, wk), q[seq_len(k)]))
  }, numeric(ncol(x)))
  coefs <- matrix(coefs, ncol = num_comp)
  rownames(coefs) <- colnames(x)
  fitted <- ym + xm %*% coefs[, num_comp]
  ss_t <- sum((y - ym)^2)
  # Variable importance in projection (Wold et al. 2001).
  ssy <- q^2 * colSums(scores^2)
  wn <- sweep(w, 2, sqrt(colSums(w^2)), "/")
  vip <- sqrt(ncol(x) * as.vector(wn^2 %*% ssy) / sum(ssy))
  structure(
    list(
      coefficients = coefs, center = center, scale = scl, y_mean = ym, num_comp = num_comp,
      predictors = colnames(x), vip = stats::setNames(vip, colnames(x)),
      r_squared = 1 - sum((y - fitted)^2) / ss_t, rmse = sqrt(mean((y - fitted)^2)),
      explained_y = cumsum(ssy) / ss_t, n = nrow(x)
    ),
    class = "soil_pls_fit"
  )
}

#' @export
print.soil_pls_fit <- function(x, ...) {
  cat_line("{.cls soil_pls_fit} PLS regression with {x$num_comp} component{?s} on {length(x$predictors)} predictor{?s}")
  cat_line("Calibration R squared {signif(x$r_squared, 3)}, RMSE {signif(x$rmse, 3)} (n = {x$n})")
  invisible(x)
}

#' Predict from a fitted PLS model
#'
#' @param object A `soil_pls_fit` object.
#' @param new_data A data frame or matrix with the predictor columns.
#' @param num_comp Number of components to use (at most the number fitted).
#' @param ... Not used.
#' @return A numeric vector of predictions.
#' @family parsnip models
#' @export
#' @examples
#' nir <- soil_nir[grep("^nm_", names(soil_nir))]
#' fit <- soil_pls_fit(nir, soil_nir$clay, num_comp = 6)
#' head(predict(fit, nir, num_comp = 3))
predict.soil_pls_fit <- function(object, new_data, num_comp = NULL, ...) {
  nc <- num_comp %||% object$num_comp
  check_integerish(nc)
  if (nc > object$num_comp) {
    cli::cli_abort("{.arg num_comp} ({nc}) exceeds the {object$num_comp} component{?s} fitted.")
  }
  nd <- as.data.frame(new_data)
  miss <- setdiff(object$predictors, names(nd))
  if (length(miss)) {
    cli::cli_abort("{.arg new_data} is missing {length(miss)} predictor{?s}, for example {.field {miss[1]}}.")
  }
  x <- as.matrix(nd[object$predictors])
  xm <- sweep(sweep(x, 2, object$center), 2, object$scale, "/")
  as.vector(object$y_mean + xm %*% object$coefficients[, nc])
}

# Registration with parsnip ---------------------------------------------------

make_soil_pls <- function() {
  parsnip::set_new_model("soil_pls")
  parsnip::set_model_mode("soil_pls", "regression")
  parsnip::set_model_engine("soil_pls", "regression", "nipals")
  parsnip::set_dependency("soil_pls", "nipals", "tidysoilr", mode = "regression")
  parsnip::set_model_arg(
    model = "soil_pls", eng = "nipals", parsnip = "num_comp", original = "num_comp",
    func = list(pkg = "dials", fun = "num_comp"), has_submodel = FALSE
  )
  parsnip::set_fit(
    model = "soil_pls", eng = "nipals", mode = "regression",
    value = list(
      interface = "data.frame", protect = c("x", "y"),
      func = c(pkg = "tidysoilr", fun = "soil_pls_fit"),
      defaults = list()
    )
  )
  parsnip::set_encoding(
    model = "soil_pls", eng = "nipals", mode = "regression",
    options = list(
      predictor_indicators = "traditional", compute_intercept = FALSE,
      remove_intercept = TRUE, allow_sparse_x = FALSE
    )
  )
  parsnip::set_pred(
    model = "soil_pls", eng = "nipals", mode = "regression", type = "numeric",
    value = list(
      pre = NULL, post = NULL, func = c(fun = "predict"),
      args = list(object = quote(object$fit), new_data = quote(new_data))
    )
  )
}
