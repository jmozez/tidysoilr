# Custom recipes steps ------------------------------------------------------

#' Standard normal variate recipe step
#'
#' `step_soil_snv()` creates a *specification* of a recipe step that
#' applies the standard normal variate transform ([soil_snv()]) row by row
#' across the selected spectral columns, replacing them in place. Select
#' the columns in order of wavelength.
#'
#' @param recipe A recipe object.
#' @param ... One or more selector functions choosing the spectral columns.
#' @param role Not used by this step since no new variables are created.
#' @param trained A logical to indicate if the step has been trained.
#' @param columns A character vector of the selected columns, filled in by
#'   [prep()].
#' @param skip A logical. Should the step be skipped when the recipe is baked
#'   by [bake()]?
#' @param id A character string that is unique to this step.
#' @param x A step object.
#' @param training A data frame used to train the step.
#' @param info Variable information from the recipe.
#' @param object A trained step.
#' @param new_data A data frame to process.
#' @param width Print width.
#'
#' @return An updated recipe with the new step added. [tidy()] returns a
#'   tibble with columns `terms` and `id`.
#' @inherit soil_snv references
#' @family recipes steps
#' @export
#' @examples
#' library(recipes)
#' rec <- recipe(clay ~ ., data = soil_nir[c("clay", paste0("nm_", seq(400, 2480, 20)))]) |>
#'   step_soil_snv(starts_with("nm_"))
#' prepped <- prep(rec)
#' bake(prepped, new_data = NULL)[1:3, 1:5]
#' tidy(prepped, number = 1)
step_soil_snv <- function(recipe, ..., role = NA, trained = FALSE, columns = NULL,
                          skip = FALSE, id = recipes::rand_id("soil_snv")) {
  recipes::add_step(
    recipe,
    step_soil_snv_new(terms = rlang::enquos(...), role = role, trained = trained,
      columns = columns, skip = skip, id = id)
  )
}

step_soil_snv_new <- function(terms, role, trained, columns, skip, id) {
  recipes::step(subclass = "soil_snv", terms = terms, role = role, trained = trained,
    columns = columns, skip = skip, id = id)
}

#' @rdname step_soil_snv
#' @export
prep.step_soil_snv <- function(x, training, info = NULL, ...) {
  col_names <- recipes::recipes_eval_select(x$terms, training, info)
  recipes::check_type(training[, col_names], types = c("double", "integer"))
  if (length(col_names) < 3) {
    cli::cli_abort("{.fn step_soil_snv} needs at least three spectral columns.")
  }
  step_soil_snv_new(terms = x$terms, role = x$role, trained = TRUE, columns = unname(col_names),
    skip = x$skip, id = x$id)
}

#' @rdname step_soil_snv
#' @export
bake.step_soil_snv <- function(object, new_data, ...) {
  recipes::check_new_data(object$columns, object, new_data)
  m <- as.matrix(new_data[object$columns])
  out <- snv_matrix(m)
  for (j in seq_along(object$columns)) new_data[[object$columns[j]]] <- out[, j]
  new_data
}

#' @rdname step_soil_snv
#' @export
print.step_soil_snv <- function(x, width = max(20, options()$width - 30), ...) {
  recipes::print_step(x$columns, x$terms, x$trained, "Standard normal variate on ", width)
  invisible(x)
}

#' @rdname step_soil_snv
#' @export
tidy.step_soil_snv <- function(x, ...) {
  terms <- if (recipes::is_trained(x)) x$columns else recipes::sel2char(x$terms)
  tibble::tibble(terms = terms, id = x$id)
}

#' @rdname step_soil_snv
#' @export
required_pkgs.step_soil_snv <- function(x, ...) {
  c("tidysoilr")
}

# Savitzky-Golay ------------------------------------------------------------

#' Savitzky-Golay recipe step
#'
#' `step_soil_savgol()` creates a *specification* of a recipe step that
#' applies a Savitzky-Golay filter ([soil_savgol()]) across the selected
#' spectral columns. The `(window - 1) / 2` columns at each end of the
#' spectrum are removed, because the filter cannot be computed there; the
#' remaining columns are replaced in place.
#'
#' @inheritParams step_soil_snv
#' @inheritParams soil_savgol
#' @param removed Edge columns removed, filled in by [prep()].
#'
#' @return An updated recipe with the new step added. [tidy()] returns a
#'   tibble with columns `terms`, `window`, `poly`, `deriv` and `id`.
#' @inherit soil_savgol references
#' @family recipes steps
#' @export
#' @examples
#' library(recipes)
#' rec <- recipe(clay ~ ., data = soil_nir[c("clay", paste0("nm_", seq(400, 2480, 20)))]) |>
#'   step_soil_savgol(starts_with("nm_"), window = 7, deriv = 1)
#' prepped <- prep(rec)
#' ncol(bake(prepped, new_data = NULL))
#' tidy(prepped, number = 1)
step_soil_savgol <- function(recipe, ..., role = NA, trained = FALSE, window = 11, poly = 2,
                             deriv = 0, columns = NULL, removed = NULL, skip = FALSE,
                             id = recipes::rand_id("soil_savgol")) {
  check_savgol(window, poly, deriv, NULL)
  recipes::add_step(
    recipe,
    step_soil_savgol_new(terms = rlang::enquos(...), role = role, trained = trained,
      window = window, poly = poly, deriv = deriv, columns = columns, removed = removed,
      skip = skip, id = id)
  )
}

step_soil_savgol_new <- function(terms, role, trained, window, poly, deriv, columns, removed,
                                 skip, id) {
  recipes::step(subclass = "soil_savgol", terms = terms, role = role, trained = trained,
    window = window, poly = poly, deriv = deriv, columns = columns, removed = removed,
    skip = skip, id = id)
}

#' @rdname step_soil_savgol
#' @export
prep.step_soil_savgol <- function(x, training, info = NULL, ...) {
  col_names <- unname(recipes::recipes_eval_select(x$terms, training, info))
  recipes::check_type(training[, col_names], types = c("double", "integer"))
  check_savgol(x$window, x$poly, x$deriv, length(col_names))
  h <- (x$window - 1) / 2
  removed <- col_names[c(seq_len(h), length(col_names) - seq_len(h) + 1)]
  step_soil_savgol_new(terms = x$terms, role = x$role, trained = TRUE, window = x$window,
    poly = x$poly, deriv = x$deriv, columns = col_names, removed = removed, skip = x$skip,
    id = x$id)
}

#' @rdname step_soil_savgol
#' @export
bake.step_soil_savgol <- function(object, new_data, ...) {
  recipes::check_new_data(object$columns, object, new_data)
  res <- savgol_matrix(as.matrix(new_data[object$columns]), object$window, object$poly,
    object$deriv)
  kept <- object$columns[res$idx]
  for (j in seq_along(kept)) new_data[[kept[j]]] <- res$m[, j]
  new_data[setdiff(names(new_data), object$removed)]
}

#' @rdname step_soil_savgol
#' @export
print.step_soil_savgol <- function(x, width = max(20, options()$width - 30), ...) {
  title <- paste0("Savitzky-Golay (window ", x$window, ", poly ", x$poly, ", deriv ", x$deriv,
    ") on ")
  recipes::print_step(x$columns, x$terms, x$trained, title, width)
  invisible(x)
}

#' @rdname step_soil_savgol
#' @export
tidy.step_soil_savgol <- function(x, ...) {
  terms <- if (recipes::is_trained(x)) setdiff(x$columns, x$removed) else recipes::sel2char(x$terms)
  tibble::tibble(terms = terms, window = x$window, poly = x$poly, deriv = x$deriv, id = x$id)
}

#' @rdname step_soil_savgol
#' @export
required_pkgs.step_soil_savgol <- function(x, ...) {
  c("tidysoilr")
}

# Texture log-ratios -----------------------------------------------------------

#' Isometric log-ratio recipe step for soil texture
#'
#' `step_soil_ilr()` creates a *specification* of a recipe step that
#' replaces sand, silt and clay with their two isometric log-ratio
#' coordinates ([soil_texture_ilr()]), so that models see unconstrained
#' predictors instead of a closed composition.
#'
#' @inheritParams step_soil_snv
#' @param ... Selector functions choosing exactly three columns, in the
#'   order sand, silt, clay.
#' @param role Role of the new columns.
#' @param prefix Prefix of the new column names (`<prefix>1`, `<prefix>2`).
#' @param zero Replacement value for zeros (percent).
#' @param keep_original_cols Keep the sand, silt and clay columns?
#'
#' @return An updated recipe with the new step added. [tidy()] returns a
#'   tibble with columns `terms`, `value` (the new column names) and `id`.
#' @inherit soil_texture_ilr references
#' @family recipes steps
#' @export
#' @examples
#' library(recipes)
#' dat <- na.omit(png_profiles[c("oc", "sand", "silt", "clay")])
#' rec <- recipe(oc ~ sand + silt + clay, data = dat) |>
#'   step_soil_ilr(sand, silt, clay)
#' bake(prep(rec), new_data = NULL)
step_soil_ilr <- function(recipe, ..., role = "predictor", trained = FALSE, prefix = "ilr_",
                          zero = 0.5, keep_original_cols = FALSE, columns = NULL,
                          skip = FALSE, id = recipes::rand_id("soil_ilr")) {
  check_name(prefix)
  check_number_positive(zero)
  recipes::add_step(
    recipe,
    step_soil_ilr_new(terms = rlang::enquos(...), role = role, trained = trained,
      prefix = prefix, zero = zero, keep_original_cols = keep_original_cols,
      columns = columns, skip = skip, id = id)
  )
}

step_soil_ilr_new <- function(terms, role, trained, prefix, zero, keep_original_cols, columns,
                              skip, id) {
  recipes::step(subclass = "soil_ilr", terms = terms, role = role, trained = trained,
    prefix = prefix, zero = zero, keep_original_cols = keep_original_cols, columns = columns,
    skip = skip, id = id)
}

#' @rdname step_soil_ilr
#' @export
prep.step_soil_ilr <- function(x, training, info = NULL, ...) {
  col_names <- unname(recipes::recipes_eval_select(x$terms, training, info))
  if (length(col_names) != 3) {
    cli::cli_abort("{.fn step_soil_ilr} needs exactly three columns (sand, silt, clay), not {length(col_names)}.")
  }
  recipes::check_type(training[, col_names], types = c("double", "integer"))
  step_soil_ilr_new(terms = x$terms, role = x$role, trained = TRUE, prefix = x$prefix,
    zero = x$zero, keep_original_cols = x$keep_original_cols, columns = col_names,
    skip = x$skip, id = x$id)
}

#' @rdname step_soil_ilr
#' @export
bake.step_soil_ilr <- function(object, new_data, ...) {
  recipes::check_new_data(object$columns, object, new_data)
  cols <- object$columns
  z <- texture_ilr(new_data[[cols[1]]], new_data[[cols[2]]], new_data[[cols[3]]], object$zero)
  nm <- paste0(object$prefix, 1:2)
  new_data[[nm[1]]] <- unname(z[, 1])
  new_data[[nm[2]]] <- unname(z[, 2])
  if (!isTRUE(object$keep_original_cols)) {
    new_data <- new_data[setdiff(names(new_data), cols)]
  }
  new_data
}

#' @rdname step_soil_ilr
#' @export
print.step_soil_ilr <- function(x, width = max(20, options()$width - 30), ...) {
  recipes::print_step(x$columns, x$terms, x$trained, "Texture log-ratios from ", width)
  invisible(x)
}

#' @rdname step_soil_ilr
#' @export
tidy.step_soil_ilr <- function(x, ...) {
  terms <- if (recipes::is_trained(x)) x$columns else recipes::sel2char(x$terms)
  tibble::tibble(terms = terms, value = paste(paste0(x$prefix, 1:2), collapse = ", "), id = x$id)
}

#' @rdname step_soil_ilr
#' @export
required_pkgs.step_soil_ilr <- function(x, ...) {
  c("tidysoilr")
}

# Texture class -----------------------------------------------------------------

#' Texture class recipe step
#'
#' `step_soil_texture()` creates a *specification* of a recipe step that
#' adds a factor of USDA or FAO texture classes ([soil_texture_class()])
#' computed from sand, silt and clay. The factor always has every class as
#' a level, so dummy variables are consistent between data sets.
#'
#' @inheritParams step_soil_ilr
#' @param system `"usda"` or `"fao"`.
#' @param name Name of the new column.
#' @param keep_original_cols Keep the sand, silt and clay columns?
#'
#' @return An updated recipe with the new step added. [tidy()] returns a
#'   tibble with columns `terms`, `system` and `id`.
#' @inherit soil_texture_class references
#' @family recipes steps
#' @export
#' @examples
#' library(recipes)
#' dat <- na.omit(png_profiles[c("oc", "sand", "silt", "clay")])
#' rec <- recipe(oc ~ sand + silt + clay, data = dat) |>
#'   step_soil_texture(sand, silt, clay)
#' table(bake(prep(rec), new_data = NULL)$texture_class)
step_soil_texture <- function(recipe, ..., role = "predictor", trained = FALSE,
                              system = c("usda", "fao"), name = "texture_class",
                              keep_original_cols = TRUE, columns = NULL, skip = FALSE,
                              id = recipes::rand_id("soil_texture")) {
  system <- rlang::arg_match(system)
  check_name(name)
  recipes::add_step(
    recipe,
    step_soil_texture_new(terms = rlang::enquos(...), role = role, trained = trained,
      system = system, name = name, keep_original_cols = keep_original_cols,
      columns = columns, skip = skip, id = id)
  )
}

step_soil_texture_new <- function(terms, role, trained, system, name, keep_original_cols,
                                  columns, skip, id) {
  recipes::step(subclass = "soil_texture", terms = terms, role = role, trained = trained,
    system = system, name = name, keep_original_cols = keep_original_cols, columns = columns,
    skip = skip, id = id)
}

#' @rdname step_soil_texture
#' @export
prep.step_soil_texture <- function(x, training, info = NULL, ...) {
  col_names <- unname(recipes::recipes_eval_select(x$terms, training, info))
  if (length(col_names) != 3) {
    cli::cli_abort("{.fn step_soil_texture} needs exactly three columns (sand, silt, clay), not {length(col_names)}.")
  }
  recipes::check_type(training[, col_names], types = c("double", "integer"))
  step_soil_texture_new(terms = x$terms, role = x$role, trained = TRUE, system = x$system,
    name = x$name, keep_original_cols = x$keep_original_cols, columns = col_names,
    skip = x$skip, id = x$id)
}

#' @rdname step_soil_texture
#' @export
bake.step_soil_texture <- function(object, new_data, ...) {
  recipes::check_new_data(object$columns, object, new_data)
  cols <- object$columns
  s <- new_data[[cols[1]]]
  si <- new_data[[cols[2]]]
  cl <- new_data[[cols[3]]]
  tot <- s + si + cl
  k <- ifelse(tot > 0, 100 / tot, NA_real_)
  fn <- if (object$system == "usda") classify_usda else classify_fao
  new_data[[object$name]] <- fn(s * k, si * k, cl * k)
  if (!isTRUE(object$keep_original_cols)) {
    new_data <- new_data[setdiff(names(new_data), cols)]
  }
  new_data
}

#' @rdname step_soil_texture
#' @export
print.step_soil_texture <- function(x, width = max(20, options()$width - 30), ...) {
  title <- paste0(toupper(x$system), " texture class from ")
  recipes::print_step(x$columns, x$terms, x$trained, title, width)
  invisible(x)
}

#' @rdname step_soil_texture
#' @export
tidy.step_soil_texture <- function(x, ...) {
  terms <- if (recipes::is_trained(x)) x$columns else recipes::sel2char(x$terms)
  tibble::tibble(terms = terms, system = x$system, id = x$id)
}

#' @rdname step_soil_texture
#' @export
required_pkgs.step_soil_texture <- function(x, ...) {
  c("tidysoilr")
}
