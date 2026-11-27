# Internal validation and helper functions --------------------------------

check_data_frame <- function(x, arg = rlang::caller_arg(x),
                             call = rlang::caller_env()) {
  if (!is.data.frame(x)) {
    cli::cli_abort(
      "{.arg {arg}} must be a data frame or tibble, not {.obj_type_friendly {x}}.",
      call = call
    )
  }
  invisible(x)
}

check_number_positive <- function(x, arg = rlang::caller_arg(x),
                                  call = rlang::caller_env(),
                                  allow_zero = FALSE, allow_null = FALSE) {
  if (allow_null && is.null(x)) {
    return(invisible(x))
  }
  ok <- is.numeric(x) && length(x) == 1 && !is.na(x) && is.finite(x) &&
    (if (allow_zero) x >= 0 else x > 0)
  if (!ok) {
    cli::cli_abort(
      paste0("{.arg {arg}} must be a single ", if (allow_zero) "non-negative" else "positive", " number."),
      call = call
    )
  }
  invisible(x)
}

check_integerish <- function(x, arg = rlang::caller_arg(x),
                             call = rlang::caller_env(), min = 1) {
  ok <- is.numeric(x) && length(x) == 1 && !is.na(x) && x == round(x) && x >= min
  if (!ok) {
    cli::cli_abort("{.arg {arg}} must be a single whole number of at least {min}.", call = call)
  }
  invisible(x)
}

check_proportion <- function(x, arg = rlang::caller_arg(x),
                             call = rlang::caller_env()) {
  ok <- is.numeric(x) && length(x) == 1 && !is.na(x) && x > 0 && x < 1
  if (!ok) {
    cli::cli_abort("{.arg {arg}} must be a single number between 0 and 1.", call = call)
  }
  invisible(x)
}

check_name <- function(x, arg = rlang::caller_arg(x), call = rlang::caller_env()) {
  if (!rlang::is_string(x) || !nzchar(x)) {
    cli::cli_abort("{.arg {arg}} must be a single non-empty string.", call = call)
  }
  invisible(x)
}

# Resolve a single column from a tidyselect quosure. Returns the column name,
# or NULL when the quosure is NULL and the column is optional.
select_one <- function(quo, data, arg, required = TRUE,
                       call = rlang::caller_env()) {
  if (rlang::quo_is_null(quo) || rlang::quo_is_missing(quo)) {
    if (required) {
      cli::cli_abort("{.arg {arg}} is required and must name a column.", call = call)
    }
    return(NULL)
  }
  loc <- tidyselect::eval_select(quo, data, allow_rename = FALSE, error_call = call)
  if (length(loc) != 1) {
    cli::cli_abort(
      "{.arg {arg}} must select exactly one column, not {length(loc)}.",
      call = call
    )
  }
  names(loc)
}

# Resolve any number of columns from a tidyselect quosure.
select_many <- function(quo, data, arg, call = rlang::caller_env()) {
  if (rlang::quo_is_null(quo) || rlang::quo_is_missing(quo)) {
    return(character())
  }
  names(tidyselect::eval_select(quo, data, allow_rename = FALSE, error_call = call))
}

# Resolve a numeric column and return its values.
numeric_col <- function(quo, data, arg, required = TRUE, non_negative = FALSE,
                        call = rlang::caller_env()) {
  col <- select_one(quo, data, arg, required = required, call = call)
  if (is.null(col)) {
    return(NULL)
  }
  x <- data[[col]]
  if (!is.numeric(x)) {
    cli::cli_abort("Column {.field {col}} ({.arg {arg}}) must be numeric, not {.obj_type_friendly {x}}.",
      call = call
    )
  }
  if (non_negative && any(x < 0, na.rm = TRUE)) {
    cli::cli_abort("Column {.field {col}} ({.arg {arg}}) must not contain negative values.",
      call = call
    )
  }
  x
}

# Add or replace columns on a data frame while keeping its class.
add_cols <- function(.data, cols) {
  for (nm in names(cols)) {
    .data[[nm]] <- cols[[nm]]
  }
  .data
}

# Drop tidysoilr metadata and return a plain tibble.
as_plain_tibble <- function(x) {
  if (inherits(x, "sf")) {
    rlang::check_installed("sf", reason = "to drop the geometry of an sf object.")
    x <- sf::st_drop_geometry(x)
  }
  attr(x, "soil_meta") <- NULL
  tibble::as_tibble(as.data.frame(x, stringsAsFactors = FALSE))
}

# Split-apply-combine over `.by` groups without data masking.
by_apply <- function(.data, by_cols, fn) {
  data <- as_plain_tibble(dplyr::ungroup(.data))
  if (!length(by_cols)) {
    return(tibble::as_tibble(fn(data)))
  }
  pieces <- vctrs::vec_split(data, data[by_cols])
  res <- lapply(seq_along(pieces$val), function(i) {
    r <- tibble::as_tibble(fn(pieces$val[[i]]))
    key <- pieces$key[rep(i, nrow(r)), , drop = FALSE]
    dplyr::bind_cols(key, r[setdiff(names(r), by_cols)])
  })
  dplyr::bind_rows(res)
}

# Resolve `.by` columns, falling back to dplyr groups.
resolve_by <- function(by_quo, .data, call = rlang::caller_env()) {
  by_cols <- select_many(by_quo, .data, ".by", call = call)
  if (!length(by_cols) && inherits(.data, "grouped_df")) {
    by_cols <- dplyr::group_vars(.data)
  }
  by_cols
}

# Planar or great-circle distance in kilometres.
dist_km <- function(x1, y1, x2, y2, geographic = TRUE) {
  if (!geographic) {
    return(sqrt((x2 - x1)^2 + (y2 - y1)^2) / 1000)
  }
  rad <- pi / 180
  dlat <- (y2 - y1) * rad
  dlon <- (x2 - x1) * rad
  a <- sin(dlat / 2)^2 + cos(y1 * rad) * cos(y2 * rad) * sin(dlon / 2)^2
  2 * 6371.0088 * asin(pmin(1, sqrt(a)))
}

# Print a line formatted with cli inline markup to standard output.
cat_line <- function(x, .envir = parent.frame()) {
  cat(cli::format_inline(x, .envir = .envir), "\n", sep = "")
}

# Resolve columns selected through `...`.
select_dots <- function(data, ..., call = rlang::caller_env()) {
  names(tidyselect::eval_select(rlang::expr(c(...)), data, allow_rename = FALSE,
    error_call = call))
}
