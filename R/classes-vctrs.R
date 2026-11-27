# S3 tibble subclasses built on vctrs ------------------------------------
#
# Every tidysoilr class stores its metadata in a single "soil_meta"
# attribute. `required` lists the columns that must survive a dplyr verb for
# the class to be kept; if any is dropped the result falls back to a tibble.

new_soil_tbl <- function(x, meta, class) {
  if (!is.list(x)) {
    cli::cli_abort("{.arg x} must be a list or data frame.")
  }
  n <- if (is.data.frame(x)) nrow(x) else vctrs::vec_size_common(!!!x)
  x <- as.list(x)
  attr(x, "soil_meta") <- NULL
  out <- tibble::new_tibble(x, nrow = n, class = c(class, "soil_tbl"))
  attr(out, "soil_meta") <- meta
  out
}

soil_meta <- function(x) {
  attr(x, "soil_meta", exact = TRUE) %||% list()
}

soil_class <- function(x) {
  cls <- class(x)
  cls[seq_len(match("soil_tbl", cls, nomatch = 1L) - 1L)][1]
}

soil_restore <- function(x, to) {
  meta <- soil_meta(to)
  cls <- soil_class(to)
  required <- meta$required %||% character()
  plain <- x
  attr(plain, "soil_meta") <- NULL
  if (is.na(cls) || !all(required %in% names(x))) {
    class(plain) <- c("tbl_df", "tbl", "data.frame")
    return(plain)
  }
  new_soil_tbl(plain, meta, cls)
}

same_meta <- function(x, y) {
  keys <- c("required", "id", "top", "bottom", "spectra", "wavelengths")
  identical(soil_meta(x)[keys], soil_meta(y)[keys])
}

soil_ptype2 <- function(x, y, ...) {
  out <- vctrs::tib_ptype2(x, y, ...)
  if (same_meta(x, y)) soil_restore(out, x) else out
}

soil_cast <- function(x, to, ...) {
  out <- vctrs::tib_cast(x, to, ...)
  soil_restore(out, to)
}

#' @export
vec_restore.soil_tbl <- function(x, to, ...) {
  soil_restore(x, to)
}

#' @export
dplyr_reconstruct.soil_tbl <- function(data, template) {
  soil_restore(data, template)
}

#' @export
`[.soil_tbl` <- function(x, i, j, drop = FALSE, ...) {
  out <- NextMethod()
  if (!is.data.frame(out)) {
    return(out)
  }
  soil_restore(out, x)
}

#' @export
`names<-.soil_tbl` <- function(x, value) {
  old <- names(x)
  meta <- soil_meta(x)
  out <- NextMethod()
  if (length(old) == length(value)) {
    map <- stats::setNames(value, old)
    for (f in c("required", "id", "top", "bottom", "horizon", "spectra", "groups", "h", "theta")) {
      v <- meta[[f]]
      if (is.character(v) && length(v)) {
        hit <- v %in% old
        v[hit] <- unname(map[v[hit]])
        meta[[f]] <- v
      }
    }
    attr(out, "soil_meta") <- meta
  }
  out
}

#' @export
tbl_sum.soil_tbl <- function(x, ...) {
  out <- NextMethod()
  names(out)[1] <- paste0("A <", soil_class(x), ">")
  summary_fn <- switch(soil_class(x),
    soil_profile = profile_summary,
    soil_spectra = spectra_summary,
    soil_retention = retention_summary,
    function(x) character()
  )
  c(out, summary_fn(x))
}

profile_summary <- function(x) {
  meta <- soil_meta(x)
  n <- length(unique(x[[meta$id]]))
  rng <- suppressWarnings(range(c(x[[meta$top]], x[[meta$bottom]]), na.rm = TRUE))
  c(
    Profiles = paste0(n, " (id: ", meta$id, ")"),
    Depths = paste0(format(rng[1]), " to ", format(rng[2]), " (", meta$top, ", ", meta$bottom, ")")
  )
}

spectra_summary <- function(x) {
  meta <- soil_meta(x)
  wl <- meta$wavelengths
  c(
    Spectra = paste0(length(wl), " bands, ", format(min(wl)), " to ", format(max(wl)), " (", meta$type, ")"),
    Preprocessing = if (length(meta$history)) paste(meta$history, collapse = " > ") else "none"
  )
}

retention_summary <- function(x) {
  c(Model = "van Genuchten (m = 1 - 1/n)")
}

# Double dispatch methods. vctrs dispatches on the first class only, so each
# subclass gets its own set of methods.

#' @export
vec_ptype2.soil_profile.soil_profile <- function(x, y, ...) soil_ptype2(x, y, ...)
#' @export
vec_ptype2.soil_profile.tbl_df <- function(x, y, ...) vctrs::tib_ptype2(x, y, ...)
#' @export
vec_ptype2.tbl_df.soil_profile <- function(x, y, ...) vctrs::tib_ptype2(x, y, ...)
#' @export
vec_ptype2.soil_profile.data.frame <- function(x, y, ...) vctrs::tib_ptype2(x, y, ...)
#' @export
vec_ptype2.data.frame.soil_profile <- function(x, y, ...) vctrs::tib_ptype2(x, y, ...)
#' @export
vec_cast.soil_profile.soil_profile <- function(x, to, ...) soil_cast(x, to, ...)
#' @export
vec_cast.soil_profile.tbl_df <- function(x, to, ...) soil_cast(x, to, ...)
#' @export
vec_cast.tbl_df.soil_profile <- function(x, to, ...) vctrs::tib_cast(x, to, ...)
#' @export
vec_cast.soil_profile.data.frame <- function(x, to, ...) soil_cast(x, to, ...)
#' @export
vec_cast.data.frame.soil_profile <- function(x, to, ...) vctrs::df_cast(x, to, ...)

#' @export
vec_ptype2.soil_spectra.soil_spectra <- function(x, y, ...) soil_ptype2(x, y, ...)
#' @export
vec_ptype2.soil_spectra.tbl_df <- function(x, y, ...) vctrs::tib_ptype2(x, y, ...)
#' @export
vec_ptype2.tbl_df.soil_spectra <- function(x, y, ...) vctrs::tib_ptype2(x, y, ...)
#' @export
vec_ptype2.soil_spectra.data.frame <- function(x, y, ...) vctrs::tib_ptype2(x, y, ...)
#' @export
vec_ptype2.data.frame.soil_spectra <- function(x, y, ...) vctrs::tib_ptype2(x, y, ...)
#' @export
vec_cast.soil_spectra.soil_spectra <- function(x, to, ...) soil_cast(x, to, ...)
#' @export
vec_cast.soil_spectra.tbl_df <- function(x, to, ...) soil_cast(x, to, ...)
#' @export
vec_cast.tbl_df.soil_spectra <- function(x, to, ...) vctrs::tib_cast(x, to, ...)
#' @export
vec_cast.soil_spectra.data.frame <- function(x, to, ...) soil_cast(x, to, ...)
#' @export
vec_cast.data.frame.soil_spectra <- function(x, to, ...) vctrs::df_cast(x, to, ...)

#' @export
vec_ptype2.soil_retention.soil_retention <- function(x, y, ...) soil_ptype2(x, y, ...)
#' @export
vec_ptype2.soil_retention.tbl_df <- function(x, y, ...) vctrs::tib_ptype2(x, y, ...)
#' @export
vec_ptype2.tbl_df.soil_retention <- function(x, y, ...) vctrs::tib_ptype2(x, y, ...)
#' @export
vec_ptype2.soil_retention.data.frame <- function(x, y, ...) vctrs::tib_ptype2(x, y, ...)
#' @export
vec_ptype2.data.frame.soil_retention <- function(x, y, ...) vctrs::tib_ptype2(x, y, ...)
#' @export
vec_cast.soil_retention.soil_retention <- function(x, to, ...) soil_cast(x, to, ...)
#' @export
vec_cast.soil_retention.tbl_df <- function(x, to, ...) soil_cast(x, to, ...)
#' @export
vec_cast.tbl_df.soil_retention <- function(x, to, ...) vctrs::tib_cast(x, to, ...)
#' @export
vec_cast.soil_retention.data.frame <- function(x, to, ...) soil_cast(x, to, ...)
#' @export
vec_cast.data.frame.soil_retention <- function(x, to, ...) vctrs::df_cast(x, to, ...)

# Constructors ---------------------------------------------------------------

#' Soil profile (horizon) table
#'
#' `soil_profile()` marks a data frame of horizons or depth layers as a
#' `<soil_profile>`: a tibble that remembers which columns hold the profile
#' identifier and the upper and lower depths of each layer. Rows are sorted
#' by profile and upper depth. The class survives `dplyr` verbs such as
#' [dplyr::filter()] and [dplyr::mutate()] as long as the identifier and
#' depth columns are kept; otherwise the result falls back to a plain
#' tibble. Depth-aware functions such as [soil_spline()],
#' [soil_carbon_stock()] and [soil_depth_mean()] then find the columns
#' automatically.
#'
#' @param .data A data frame with one row per horizon or layer.
#' @param id <[`tidy-select`][tidyselect::language]> Profile identifier
#'   column.
#' @param top,bottom <[`tidy-select`][tidyselect::language]> Upper and lower
#'   depth columns, in centimetres by convention (positive downwards).
#' @param horizon <[`tidy-select`][tidyselect::language]> Optional horizon
#'   designation column.
#'
#' @return A `<soil_profile>` tibble.
#' @seealso [soil_profile_check()] to find gaps and overlaps.
#' @family soil classes
#' @export
#' @examples
#' prof <- soil_profile(png_profiles, id = profile, top = top, bottom = bottom,
#'   horizon = horizon)
#' prof
#'
#' # The class survives dplyr verbs that keep the depth columns
#' dplyr::filter(prof, top < 30)
soil_profile <- function(.data, id, top, bottom, horizon = NULL) {
  check_data_frame(.data)
  data <- as_plain_tibble(.data)
  id <- select_one(rlang::enquo(id), data, "id")
  top <- select_one(rlang::enquo(top), data, "top")
  bottom <- select_one(rlang::enquo(bottom), data, "bottom")
  horizon <- select_one(rlang::enquo(horizon), data, "horizon", required = FALSE)
  for (col in c(top, bottom)) {
    if (!is.numeric(data[[col]])) {
      cli::cli_abort("Depth column {.field {col}} must be numeric.")
    }
  }
  if (anyNA(data[[id]])) {
    cli::cli_abort("Profile identifier {.field {id}} must not contain missing values.")
  }
  bad <- which(!is.na(data[[top]]) & !is.na(data[[bottom]]) & data[[bottom]] <= data[[top]])
  if (length(bad)) {
    cli::cli_warn(c(
      "{length(bad)} layer{?s} {?has/have} a lower depth that is not below the upper depth.",
      "i" = "Use {.fn soil_profile_check} to list them."
    ))
  }
  data <- data[order(data[[id]], data[[top]]), , drop = FALSE]
  meta <- list(id = id, top = top, bottom = bottom, horizon = horizon,
    required = c(id, top, bottom))
  new_soil_tbl(data, meta, "soil_profile")
}

#' Soil spectra table
#'
#' `soil_spectra()` marks the spectral columns of a data frame (one row per
#' sample, one column per wavelength) as a `<soil_spectra>` tibble that
#' remembers the wavelengths, the measurement type and the preprocessing
#' applied so far. The preprocessing functions [soil_snv()],
#' [soil_savgol()], [soil_msc()], [soil_absorbance()] and
#' [soil_resample()] use this metadata, and non-spectral columns (sample
#' identifiers, laboratory values) are carried along untouched.
#'
#' @param .data A data frame in wide format.
#' @param cols <[`tidy-select`][tidyselect::language]> The spectral columns,
#'   in order of increasing wavelength.
#' @param wavelengths Numeric vector of wavelengths, one per spectral column.
#'   By default the numbers in the column names are used, so columns named
#'   `nm_400`, `nm_410` and so on become 400, 410 and so on.
#' @param type Measurement type: `"reflectance"` (0 to 1), `"absorbance"`
#'   (log10(1/R)) or `"preprocessed"`.
#'
#' @return A `<soil_spectra>` tibble.
#' @family soil classes
#' @export
#' @examples
#' spec <- soil_spectra(soil_nir, dplyr::starts_with("nm_"))
#' spec
soil_spectra <- function(.data, cols, wavelengths = NULL,
                         type = c("reflectance", "absorbance", "preprocessed")) {
  check_data_frame(.data)
  type <- rlang::arg_match(type)
  data <- as_plain_tibble(.data)
  cols <- select_many(rlang::enquo(cols), data, "cols")
  if (length(cols) < 3) {
    cli::cli_abort("{.arg cols} must select at least three spectral columns.")
  }
  num <- vapply(data[cols], is.numeric, logical(1))
  if (!all(num)) {
    cli::cli_abort("Spectral columns must be numeric; {.field {cols[!num]}} {?is/are} not.")
  }
  if (is.null(wavelengths)) {
    wavelengths <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", cols)))
    if (anyNA(wavelengths)) {
      cli::cli_abort(c(
        "Could not read wavelengths from the column names.",
        "i" = "Supply them with {.arg wavelengths}."
      ))
    }
  }
  if (!is.numeric(wavelengths) || length(wavelengths) != length(cols)) {
    cli::cli_abort("{.arg wavelengths} must be numeric with one value per spectral column ({length(cols)}).")
  }
  if (any(diff(wavelengths) <= 0)) {
    cli::cli_abort("Wavelengths must be strictly increasing.")
  }
  meta <- list(spectra = cols, wavelengths = wavelengths, type = type, history = character(),
    required = cols)
  new_soil_tbl(data, meta, "soil_spectra")
}

new_soil_retention <- function(fits, meta) {
  meta$required <- c("theta_r", "theta_s", "alpha", "n")
  new_soil_tbl(fits, meta, "soil_retention")
}

#' Test for tidysoilr classes
#'
#' @param x An object.
#' @return `TRUE` or `FALSE`.
#' @family soil classes
#' @export
#' @examples
#' prof <- soil_profile(png_profiles, profile, top, bottom)
#' is_soil_profile(prof)
#' is_soil_spectra(prof)
is_soil_profile <- function(x) inherits(x, "soil_profile")

#' @rdname is_soil_profile
#' @export
is_soil_spectra <- function(x) inherits(x, "soil_spectra")

#' @rdname is_soil_profile
#' @export
is_soil_retention <- function(x) inherits(x, "soil_retention")

# Resolve profile columns from a soil_profile or explicit arguments.
profile_cols <- function(.data, id_quo, top_quo, bottom_quo, call = rlang::caller_env()) {
  meta <- if (is_soil_profile(.data)) soil_meta(.data) else list()
  id <- select_one(id_quo, .data, "id", required = is.null(meta$id), call = call) %||% meta$id
  top <- select_one(top_quo, .data, "top", required = is.null(meta$top), call = call) %||% meta$top
  bottom <- select_one(bottom_quo, .data, "bottom", required = is.null(meta$bottom), call = call) %||%
    meta$bottom
  for (col in c(top, bottom)) {
    if (!is.numeric(.data[[col]])) {
      cli::cli_abort("Depth column {.field {col}} must be numeric.", call = call)
    }
  }
  list(id = id, top = top, bottom = bottom)
}

#' Convert a soil profile table to an aqp SoilProfileCollection
#'
#' Hands a [soil_profile()] over to the aqp package (Beaudette et al. 2013)
#' for its profile sketches, horizon logic and slab summaries. aqp must be
#' installed.
#'
#' @param .data A [soil_profile()].
#' @return An `aqp::SoilProfileCollection`.
#' @references Beaudette, D. E., Roudier, P. and O'Geen, A. T. (2013).
#'   Algorithms for quantitative pedology: a toolkit for soil scientists.
#'   *Computers and Geosciences*, 52, 258-268.
#'   \doi{10.1016/j.cageo.2012.10.020}
#' @family soil classes
#' @export
#' @examples
#' if (requireNamespace("aqp", quietly = TRUE)) {
#'   prof <- soil_profile(png_profiles, profile, top, bottom, horizon = horizon)
#'   spc <- soil_as_aqp(prof)
#'   spc
#' }
soil_as_aqp <- function(.data) {
  if (!is_soil_profile(.data)) {
    cli::cli_abort("{.arg .data} must be a {.cls soil_profile}; create one with {.fn soil_profile}.")
  }
  rlang::check_installed("aqp", reason = "to build a SoilProfileCollection.")
  meta <- soil_meta(.data)
  df <- as.data.frame(as_plain_tibble(.data))
  df[[meta$id]] <- as.character(df[[meta$id]])
  fml <- stats::as.formula(paste(meta$id, "~", meta$top, "+", meta$bottom))
  aqp::depths(df) <- fml
  if (!is.null(meta$horizon)) {
    aqp::hzdesgnname(df) <- meta$horizon
  }
  df
}
