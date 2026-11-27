# Profile utilities -------------------------------------------------------

#' Check soil profiles for depth errors
#'
#' Lists layers whose depths are missing or inverted, and gaps or overlaps
#' between consecutive layers of the same profile. Run it before
#' [soil_spline()], [soil_carbon_stock()] or [soil_depth_mean()], which
#' assume contiguous, non-overlapping layers.
#'
#' @param .data A [soil_profile()] or a data frame of layers.
#' @param id,top,bottom <[`tidy-select`][tidyselect::language]> Profile
#'   identifier and depth columns. Taken from a `<soil_profile>` when not
#'   given.
#'
#' @return A tibble with one row per problem and columns `id`, `top`,
#'   `bottom`, `issue` (`"missing depth"`, `"inverted depths"`, `"gap"` or
#'   `"overlap"`) and `size` (the size of the gap or overlap, in depth
#'   units). It has zero rows when no problems are found.
#' @family soil profiles
#' @export
#' @examples
#' prof <- soil_profile(png_profiles, profile, top, bottom)
#' soil_profile_check(prof)
#'
#' bad <- data.frame(id = "A", top = c(0, 10, 25), bottom = c(10, 20, 22))
#' soil_profile_check(bad, id = id, top = top, bottom = bottom)
soil_profile_check <- function(.data, id = NULL, top = NULL, bottom = NULL) {
  check_data_frame(.data)
  cols <- profile_cols(.data, rlang::enquo(id), rlang::enquo(top), rlang::enquo(bottom))
  ids <- as.character(.data[[cols$id]])
  t <- .data[[cols$top]]
  b <- .data[[cols$bottom]]
  out <- list()
  miss <- is.na(t) | is.na(b)
  if (any(miss)) {
    out[[length(out) + 1]] <- tibble::tibble(id = ids[miss], top = t[miss], bottom = b[miss],
      issue = "missing depth", size = NA_real_)
  }
  inv <- !miss & b <= t
  if (any(inv)) {
    out[[length(out) + 1]] <- tibble::tibble(id = ids[inv], top = t[inv], bottom = b[inv],
      issue = "inverted depths", size = t[inv] - b[inv])
  }
  ok <- !miss
  for (p in unique(ids[ok])) {
    sel <- which(ok & ids == p)
    sel <- sel[order(t[sel])]
    if (length(sel) < 2) next
    nxt <- sel[-1]
    prv <- sel[-length(sel)]
    d <- t[nxt] - b[prv]
    if (any(d > 0)) {
      g <- d > 0
      out[[length(out) + 1]] <- tibble::tibble(id = p, top = b[prv][g], bottom = t[nxt][g],
        issue = "gap", size = d[g])
    }
    if (any(d < 0)) {
      o <- d < 0
      out[[length(out) + 1]] <- tibble::tibble(id = p, top = t[nxt][o], bottom = b[prv][o],
        issue = "overlap", size = -d[o])
    }
  }
  if (!length(out)) {
    return(tibble::tibble(id = character(), top = numeric(), bottom = numeric(),
      issue = character(), size = numeric()))
  }
  dplyr::bind_rows(out)
}

#' Slice profiles into regular depth increments
#'
#' Splits each profile into regular slices (1 cm by default) and gives every
#' slice the values of the layer that contains its mid-point, like
#' `aqp::dice()` but for plain data frames. Slices below the deepest layer or
#' inside gaps are dropped.
#'
#' @param .data A [soil_profile()] or a data frame of layers.
#' @param ... <[`tidy-select`][tidyselect::language]> Columns to carry into
#'   the slices. Defaults to every column other than the identifier and
#'   depths.
#' @param interval Slice thickness, in depth units.
#' @param max_depth Optional maximum depth.
#' @inheritParams soil_profile_check
#'
#' @return A `<soil_profile>` with one row per slice and the identifier,
#'   `top`, `bottom` and selected columns.
#' @family soil profiles
#' @export
#' @examples
#' prof <- soil_profile(png_profiles, profile, top, bottom)
#' slices <- soil_slice(prof, clay, oc, interval = 5, max_depth = 100)
#' head(slices)
soil_slice <- function(.data, ..., interval = 1, max_depth = NULL, id = NULL, top = NULL,
                       bottom = NULL) {
  check_data_frame(.data)
  check_number_positive(interval)
  check_number_positive(max_depth, allow_null = TRUE)
  cols <- profile_cols(.data, rlang::enquo(id), rlang::enquo(top), rlang::enquo(bottom))
  data <- as_plain_tibble(.data)
  vars <- select_dots(data, ...)
  if (!length(vars)) {
    vars <- setdiff(names(data), unlist(cols))
  }
  vars <- setdiff(vars, unlist(cols))
  ok <- !is.na(data[[cols$top]]) & !is.na(data[[cols$bottom]]) &
    data[[cols$bottom]] > data[[cols$top]]
  data <- data[ok, , drop = FALSE]
  pieces <- split(seq_len(nrow(data)), factor(data[[cols$id]], levels = unique(data[[cols$id]])))
  res <- lapply(pieces, function(rows) {
    t <- data[[cols$top]][rows]
    b <- data[[cols$bottom]][rows]
    deepest <- max(b)
    if (!is.null(max_depth)) deepest <- min(deepest, max_depth)
    start <- floor(min(t) / interval) * interval
    if (deepest <= start) return(NULL)
    tops <- seq(start, deepest - interval / 2, by = interval)
    mid <- tops + interval / 2
    hit <- vapply(mid, function(m) {
      w <- which(t <= m & b > m)
      if (length(w)) w[1] else NA_integer_
    }, integer(1))
    keep <- !is.na(hit)
    out <- data[rows[hit[keep]], c(cols$id, vars), drop = FALSE]
    out[[cols$top]] <- tops[keep]
    out[[cols$bottom]] <- tops[keep] + interval
    out[c(cols$id, cols$top, cols$bottom, vars)]
  })
  out <- dplyr::bind_rows(res)
  meta <- list(id = cols$id, top = cols$top, bottom = cols$bottom, horizon = NULL,
    required = unlist(cols, use.names = FALSE))
  new_soil_tbl(out, meta, "soil_profile")
}

# Overlap (in depth units) of each layer with the interval [from, to].
layer_overlap <- function(top, bottom, from, to) {
  pmax(0, pmin(bottom, to) - pmax(top, from))
}

#' Depth-weighted mean of soil properties
#'
#' Computes the thickness-weighted mean of one or more properties over a
#' depth interval (for example 0 to 30 cm) for every profile, prorating
#' layers that cross the interval limits. Missing values are skipped, and
#' `coverage` reports the share of the interval that contributed.
#'
#' @param .data A [soil_profile()] or a data frame of layers.
#' @param ... <[`tidy-select`][tidyselect::language]> Numeric property
#'   columns.
#' @param from,to Upper and lower limits of the interval.
#' @inheritParams soil_profile_check
#'
#' @return A tibble with one row per profile: the identifier, `from`, `to`,
#'   the weighted mean of each property and `coverage` (0 to 1, for the
#'   first property).
#' @family soil profiles
#' @export
#' @examples
#' prof <- soil_profile(png_profiles, profile, top, bottom)
#' soil_depth_mean(prof, clay, oc, ph_water, from = 0, to = 30)
soil_depth_mean <- function(.data, ..., from = 0, to = 30, id = NULL, top = NULL, bottom = NULL) {
  check_data_frame(.data)
  check_number_positive(from, allow_zero = TRUE)
  check_number_positive(to)
  if (to <= from) {
    cli::cli_abort("{.arg to} ({to}) must be greater than {.arg from} ({from}).")
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
  w_all <- layer_overlap(data[[cols$top]], data[[cols$bottom]], from, to)
  w_all[is.na(w_all)] <- 0
  ids <- unique(data[[cols$id]])
  idx <- split(seq_len(nrow(data)), factor(data[[cols$id]], levels = ids))
  res <- lapply(vars, function(v) {
    vapply(idx, function(r) {
      w <- w_all[r]
      x <- data[[v]][r]
      ok <- !is.na(x) & w > 0
      if (!any(ok)) NA_real_ else sum(x[ok] * w[ok]) / sum(w[ok])
    }, numeric(1))
  })
  cover <- vapply(idx, function(r) {
    w <- w_all[r]
    ok <- !is.na(data[[vars[1]]][r]) & w > 0
    sum(w[ok]) / (to - from)
  }, numeric(1))
  out <- tibble::tibble(!!cols$id := ids, from = from, to = to)
  for (i in seq_along(vars)) out[[vars[i]]] <- unname(res[[i]])
  out$coverage <- unname(cover)
  out
}
