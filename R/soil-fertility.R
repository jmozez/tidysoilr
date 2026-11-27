# Soil chemistry and fertility ---------------------------------------------

#' Exchangeable bases, base saturation and sodicity
#'
#' From exchangeable cations (cmol(+)/kg) and the cation exchange capacity
#' computes the sum of bases, base saturation, exchangeable sodium
#' percentage (ESP) and the Ca:Mg ratio. When exchangeable acidity (Al + H)
#' is given, the effective CEC (ECEC = bases + acidity) and the aluminium
#' (acidity) saturation of the ECEC are added; acidity saturation above
#' about 60 percent indicates likely aluminium toxicity for many crops.
#'
#' @param .data A data frame.
#' @param ca,mg,k,na <[`tidy-select`][tidyselect::language]> Exchangeable
#'   calcium, magnesium, potassium and sodium (cmol(+)/kg).
#' @param cec <[`tidy-select`][tidyselect::language]> Cation exchange
#'   capacity (cmol(+)/kg), usually at pH 7.
#' @param acidity <[`tidy-select`][tidyselect::language]> Optional
#'   exchangeable acidity (Al + H, cmol(+)/kg).
#'
#' @return `.data` with `sum_bases`, `base_saturation` (percent), `esp`
#'   (percent), `ca_mg_ratio` and, with `acidity`, `ecec` and
#'   `acid_saturation` (percent).
#' @references Hazelton, P. and Murphy, B. (2016). *Interpreting Soil Test
#'   Results: What Do All the Numbers Mean?*, 3rd edition. CSIRO Publishing,
#'   Melbourne.
#' @family soil fertility
#' @export
#' @examples
#' df <- data.frame(ca = c(8, 2.1), mg = c(3, 0.8), k = c(0.6, 0.2),
#'   na = c(0.3, 0.1), cec = c(18, 12), al_h = c(0.1, 2.5))
#' soil_base_saturation(df, ca, mg, k, na, cec, acidity = al_h)
soil_base_saturation <- function(.data, ca, mg, k, na, cec, acidity = NULL) {
  check_data_frame(.data)
  ca <- numeric_col(rlang::enquo(ca), .data, "ca", non_negative = TRUE)
  mg <- numeric_col(rlang::enquo(mg), .data, "mg", non_negative = TRUE)
  k <- numeric_col(rlang::enquo(k), .data, "k", non_negative = TRUE)
  na <- numeric_col(rlang::enquo(na), .data, "na", non_negative = TRUE)
  cec <- numeric_col(rlang::enquo(cec), .data, "cec", non_negative = TRUE)
  acid <- numeric_col(rlang::enquo(acidity), .data, "acidity", required = FALSE,
    non_negative = TRUE)
  sb <- ca + mg + k + na
  bs <- 100 * sb / cec
  if (any(bs > 100, na.rm = TRUE)) {
    cli::cli_inform(c("i" = "Some base saturations exceed 100%, common in calcareous or saline soils."))
  }
  out <- list(sum_bases = sb, base_saturation = bs, esp = 100 * na / cec,
    ca_mg_ratio = ca / mg)
  if (!is.null(acid)) {
    out$ecec <- sb + acid
    out$acid_saturation <- 100 * acid / (sb + acid)
  }
  add_cols(.data, out)
}

#' Soil reaction class from pH
#'
#' Classifies pH (in water) into the soil reaction classes of the USDA Soil
#' Survey Manual, from "ultra acid" (< 3.5) to "very strongly alkaline"
#' (> 9.0).
#'
#' @param .data A data frame.
#' @param ph <[`tidy-select`][tidyselect::language]> pH measured in water.
#' @param .name Name of the new column.
#'
#' @return `.data` with an ordered factor of reaction classes.
#' @references Soil Science Division Staff (2017). *Soil Survey Manual*.
#'   USDA Handbook 18. Government Printing Office, Washington, DC.
#' @family soil fertility
#' @export
#' @examples
#' png_profiles |>
#'   soil_ph_class(ph_water) |>
#'   dplyr::count(ph_class)
soil_ph_class <- function(.data, ph, .name = "ph_class") {
  check_data_frame(.data)
  check_name(.name)
  x <- numeric_col(rlang::enquo(ph), .data, "ph")
  if (any(x < 0 | x > 14, na.rm = TRUE)) {
    cli::cli_abort("{.arg ph} values must lie between 0 and 14.")
  }
  labs <- c("ultra acid", "extremely acid", "very strongly acid", "strongly acid",
    "moderately acid", "slightly acid", "neutral", "slightly alkaline",
    "moderately alkaline", "strongly alkaline", "very strongly alkaline")
  brks <- c(-Inf, 3.45, 4.45, 5.05, 5.55, 6.05, 6.55, 7.35, 7.85, 8.45, 9.05, Inf)
  .data[[.name]] <- cut(x, breaks = brks, labels = labs, right = FALSE, ordered_result = TRUE)
  .data
}

#' Default fertility rating thresholds
#'
#' Rating limits for common soil test values, after the widely used limits
#' of Metson (1961) as tabulated by Landon (1991) and Hazelton and Murphy
#' (2016). Each row gives the lower limits of the "low", "medium", "high"
#' and "very high" classes; values below `low` are "very low". Edit the
#' table or supply your own, with the same columns, to use local
#' calibrations.
#'
#' @return A tibble with columns `property`, `units`, `low`, `medium`,
#'   `high` and `very_high`.
#' @references Metson, A. J. (1961). *Methods of Chemical Analysis for Soil
#'   Survey Samples*. New Zealand Soil Bureau Bulletin 12. Government
#'   Printer, Wellington.
#'
#'   Landon, J. R. (ed.) (1991). *Booker Tropical Soil Manual*. Longman,
#'   Harlow.
#'
#'   Hazelton, P. and Murphy, B. (2016). *Interpreting Soil Test Results:
#'   What Do All the Numbers Mean?*, 3rd edition. CSIRO Publishing,
#'   Melbourne.
#' @family soil fertility
#' @export
#' @examples
#' soil_fertility_thresholds()
soil_fertility_thresholds <- function() {
  tibble::tribble(
    ~property, ~units, ~low, ~medium, ~high, ~very_high,
    "cec", "cmol(+)/kg", 6, 12, 25, 40,
    "ca", "cmol(+)/kg", 2, 5, 10, 20,
    "mg", "cmol(+)/kg", 0.3, 1, 3, 8,
    "k", "cmol(+)/kg", 0.2, 0.3, 0.7, 2,
    "na", "cmol(+)/kg", 0.1, 0.3, 0.7, 2,
    "base_saturation", "%", 20, 40, 60, 80,
    "total_n", "%", 0.1, 0.2, 0.5, 1
  )
}

#' Rate soil test values against fertility thresholds
#'
#' Rates each selected soil test value as "very low", "low", "medium",
#' "high" or "very high" against a threshold table such as
#' [soil_fertility_thresholds()]. Name each argument after the `property`
#' it should be rated as, for example `k = exch_k`.
#'
#' @param .data A data frame.
#' @param ... Name-value pairs: the name is a `property` in `thresholds`
#'   and the value the column to rate, for example `cec = cec_ph7`. An
#'   unnamed column is rated as the property of the same name.
#' @param thresholds A threshold table with columns `property`, `low`,
#'   `medium`, `high` and `very_high`.
#'
#' @return `.data` with one ordered factor column per rated property, named
#'   `<property>_rating`.
#' @inherit soil_fertility_thresholds references
#' @family soil fertility
#' @export
#' @examples
#' png_profiles |>
#'   dplyr::mutate(total_n_pct = total_n / 10) |>
#'   soil_fertility_rating(cec, total_n = total_n_pct) |>
#'   dplyr::select(profile, top, cec, cec_rating, total_n_rating)
soil_fertility_rating <- function(.data, ..., thresholds = soil_fertility_thresholds()) {
  check_data_frame(.data)
  check_data_frame(thresholds)
  need <- c("property", "low", "medium", "high", "very_high")
  miss <- setdiff(need, names(thresholds))
  if (length(miss)) {
    cli::cli_abort("{.arg thresholds} is missing column{?s} {.field {miss}}.")
  }
  sel <- tidyselect::eval_select(rlang::expr(c(...)), .data, allow_rename = TRUE)
  if (!length(sel)) {
    cli::cli_abort("Select at least one column to rate in {.arg ...}.")
  }
  labs <- c("very low", "low", "medium", "high", "very high")
  for (i in seq_along(sel)) {
    prop <- names(sel)[i]
    col <- names(.data)[sel[i]]
    row <- thresholds[thresholds$property == prop, , drop = FALSE]
    if (nrow(row) != 1) {
      cli::cli_abort(c(
        "No thresholds for property {.val {prop}}.",
        "i" = "Name the argument after one of {.val {thresholds$property}}."
      ))
    }
    x <- .data[[col]]
    if (!is.numeric(x)) {
      cli::cli_abort("Column {.field {col}} must be numeric.")
    }
    brks <- c(-Inf, row$low, row$medium, row$high, row$very_high, Inf)
    .data[[paste0(prop, "_rating")]] <- cut(x, breaks = brks, labels = labs, right = FALSE,
      ordered_result = TRUE)
  }
  .data
}

#' Parker's nutrient index
#'
#' Summarises the fertility status of an area from the ratings of many
#' samples with the nutrient index of Parker et al. (1951):
#' \deqn{NI = \frac{1 \times N_L + 2 \times N_M + 3 \times N_H}{N_L + N_M + N_H},}
#' where \eqn{N_L}, \eqn{N_M} and \eqn{N_H} are the numbers of samples
#' rated low, medium and high. Five-class ratings are collapsed ("very low"
#' counts as low and "very high" as high). An index below 1.67 is rated
#' low, 1.67 to 2.33 medium and above 2.33 high.
#'
#' @param .data A data frame.
#' @param rating <[`tidy-select`][tidyselect::language]> A rating column, for
#'   example from [soil_fertility_rating()].
#' @param .by <[`tidy-select`][tidyselect::language]> Optional grouping
#'   columns (for example district or land use).
#'
#' @return A tibble with one row per group: `n_low`, `n_medium`, `n_high`,
#'   `nutrient_index` and `status`.
#' @references Parker, F. W., Nelson, W. L., Winters, E. and Miles, I. E.
#'   (1951). The broad interpretation and application of soil test
#'   information. *Agronomy Journal*, 43, 105-112.
#'   \doi{10.2134/agronj1951.00021962004300030001x}
#' @family soil fertility
#' @export
#' @examples
#' png_profiles |>
#'   soil_fertility_rating(cec) |>
#'   soil_nutrient_index(cec_rating, .by = wrb_group)
soil_nutrient_index <- function(.data, rating, .by = NULL) {
  check_data_frame(.data)
  col <- select_one(rlang::enquo(rating), .data, "rating")
  by_cols <- resolve_by(rlang::enquo(.by), .data)
  x <- tolower(as.character(.data[[col]]))
  ok <- x %in% c("very low", "low", "medium", "moderate", "high", "very high") | is.na(x)
  if (!all(ok)) {
    cli::cli_abort("Column {.field {col}} must hold ratings such as {.val low}, {.val medium} or {.val high}.")
  }
  data <- as_plain_tibble(dplyr::ungroup(.data))
  data$.lvl <- ifelse(x %in% c("very low", "low"), 1,
    ifelse(x %in% c("medium", "moderate"), 2, ifelse(is.na(x), NA, 3)))
  by_apply(data, by_cols, function(d) {
    l <- d$.lvl[!is.na(d$.lvl)]
    ni <- if (length(l)) mean(l) else NA_real_
    tibble::tibble(
      n_low = sum(l == 1), n_medium = sum(l == 2), n_high = sum(l == 3),
      nutrient_index = ni,
      status = factor(ifelse(is.na(ni), NA, ifelse(ni < 1.67, "low", ifelse(ni <= 2.33, "medium", "high"))),
        levels = c("low", "medium", "high"))
    )
  })
}
