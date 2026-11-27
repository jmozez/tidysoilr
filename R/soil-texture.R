# Texture classification and compositional transforms ----------------------

usda_levels <- c(
  "sand", "loamy sand", "sandy loam", "loam", "silt loam", "silt",
  "sandy clay loam", "clay loam", "silty clay loam", "sandy clay", "silty clay", "clay"
)

fao_levels <- c("coarse", "medium", "medium fine", "fine", "very fine")

# Resolve sand, silt and clay columns into a closed (sums to 100) matrix.
texture_parts <- function(.data, sand_quo, silt_quo, clay_quo, normalise = TRUE,
                          tolerance = 5, call = rlang::caller_env()) {
  sand <- numeric_col(sand_quo, .data, "sand", non_negative = TRUE, call = call)
  silt <- numeric_col(silt_quo, .data, "silt", non_negative = TRUE, call = call)
  clay <- numeric_col(clay_quo, .data, "clay", non_negative = TRUE, call = call)
  total <- sand + silt + clay
  if (all(total <= 1.5, na.rm = TRUE) && any(total > 0, na.rm = TRUE)) {
    cli::cli_abort(c(
      "Sand, silt and clay look like fractions (they sum to about 1).",
      "i" = "Give them as percentages (0 to 100)."
    ), call = call)
  }
  off <- which(abs(total - 100) > tolerance)
  if (length(off)) {
    cli::cli_warn(c(
      "{length(off)} sample{?s} {?has/have} sand + silt + clay more than {tolerance} away from 100.",
      "i" = if (normalise) "They were rescaled to sum to 100." else "They were left as they are."
    ), call = call)
  }
  if (normalise) {
    k <- ifelse(total > 0, 100 / total, NA_real_)
    sand <- sand * k
    silt <- silt * k
    clay <- clay * k
  }
  list(sand = sand, silt = silt, clay = clay)
}

classify_usda <- function(sand, silt, clay) {
  # Rules are applied in order; the first match wins.
  rules <- list(
    "silty clay" = clay >= 40 & silt >= 40,
    "clay" = clay >= 40 & sand <= 45,
    "sandy clay" = clay >= 35 & sand > 45,
    "silty clay loam" = clay >= 27 & clay < 40 & sand <= 20,
    "clay loam" = clay >= 27 & clay < 40 & sand > 20 & sand <= 45,
    "sandy clay loam" = clay >= 20 & clay < 35 & silt < 28 & sand > 45,
    "loam" = clay >= 7 & clay < 27 & silt >= 28 & silt < 50 & sand <= 52,
    "silt" = silt >= 80 & clay < 12,
    "silt loam" = (silt >= 50 & clay >= 12 & clay < 27) | (silt >= 50 & silt < 80 & clay < 12),
    "sand" = silt + 1.5 * clay < 15,
    "loamy sand" = silt + 2 * clay < 30,
    "sandy loam" = clay < 20 & silt + 2 * clay >= 30
  )
  cls <- rep(NA_character_, length(sand))
  for (nm in names(rules)) {
    hit <- is.na(cls) & !is.na(rules[[nm]]) & rules[[nm]]
    cls[hit] <- nm
  }
  factor(cls, levels = usda_levels)
}

classify_fao <- function(sand, silt, clay) {
  cls <- ifelse(clay >= 60, "very fine",
    ifelse(clay >= 35, "fine",
      ifelse(sand < 15, "medium fine",
        ifelse(clay < 18 & sand > 65, "coarse", "medium")
      )
    )
  )
  factor(cls, levels = fao_levels)
}

#' Soil texture class from sand, silt and clay
#'
#' Assigns each sample to a texture class of the USDA texture triangle
#' (12 classes; Soil Science Division Staff 2017) or the FAO texture
#' triangle used for the FAO-UNESCO Soil Map of the World and the
#' Harmonized World Soil Database (5 classes: coarse, medium, medium fine,
#' fine and very fine). Fractions are rescaled to sum to 100 first, with a
#' warning when the sum is far from 100.
#'
#' The USDA classes use the sand (0.05 to 2 mm), silt (0.002 to 0.05 mm) and
#' clay (< 0.002 mm) limits of the USDA system. Data measured with the
#' 0.063 mm sand limit of the FAO and ISO systems give slightly different
#' classes near the boundaries.
#'
#' @param .data A data frame.
#' @param sand,silt,clay <[`tidy-select`][tidyselect::language]> Columns of
#'   percentages (0 to 100).
#' @param system `"usda"` or `"fao"`.
#' @param .name Name of the new column.
#' @param normalise Rescale the three fractions to sum to 100?
#'
#' @return `.data` with a new factor column (named by `.name`) of texture
#'   classes, in the order of the classes on the triangle.
#' @references Soil Science Division Staff (2017). *Soil Survey Manual*.
#'   USDA Handbook 18. Government Printing Office, Washington, DC.
#'
#'   FAO (2006). *Guidelines for Soil Description*, 4th edition. FAO, Rome.
#'
#'   FAO/IIASA/ISRIC/ISSCAS/JRC (2012). *Harmonized World Soil Database*
#'   (version 1.2). FAO, Rome and IIASA, Laxenburg.
#' @family soil texture
#' @export
#' @examples
#' png_profiles |>
#'   soil_texture_class(sand, silt, clay) |>
#'   soil_texture_class(sand, silt, clay, system = "fao", .name = "fao_texture") |>
#'   dplyr::select(profile, top, sand, silt, clay, texture, fao_texture)
soil_texture_class <- function(.data, sand, silt, clay, system = c("usda", "fao"),
                               .name = "texture", normalise = TRUE) {
  check_data_frame(.data)
  system <- rlang::arg_match(system)
  check_name(.name)
  p <- texture_parts(.data, rlang::enquo(sand), rlang::enquo(silt), rlang::enquo(clay),
    normalise = normalise)
  fn <- if (system == "usda") classify_usda else classify_fao
  .data[[.name]] <- fn(p$sand, p$silt, p$clay)
  .data
}

#' Texture triangle class boundaries
#'
#' Returns the vertices of the texture class polygons of the USDA or FAO
#' texture triangle, for plotting (see [soil_texture_triangle()]) or for
#' your own classification code. The polygons follow the same rules as
#' [soil_texture_class()].
#'
#' @param system `"usda"` or `"fao"`.
#' @return A tibble with columns `class`, `vertex`, `sand`, `silt` and
#'   `clay` (percentages).
#' @family soil texture
#' @export
#' @examples
#' polys <- soil_texture_polygons("usda")
#' head(polys)
#' unique(polys$class)
soil_texture_polygons <- function(system = c("usda", "fao")) {
  system <- rlang::arg_match(system)
  # Vertices as (sand, clay); silt is the remainder.
  v <- if (system == "usda") {
    list(
      "sand" = c(85, 0, 90, 10, 100, 0),
      "loamy sand" = c(70, 0, 85, 15, 90, 10, 85, 0),
      "sandy loam" = c(52, 20, 80, 20, 85, 15, 70, 0, 50, 0, 43, 7, 52, 7),
      "loam" = c(45, 27, 52, 20, 52, 7, 43, 7, 23, 27),
      "silt loam" = c(0, 27, 23, 27, 50, 0, 20, 0, 8, 12, 0, 12),
      "silt" = c(0, 12, 8, 12, 20, 0, 0, 0),
      "sandy clay loam" = c(45, 35, 65, 35, 80, 20, 52, 20, 45, 27),
      "clay loam" = c(20, 40, 20, 27, 45, 27, 45, 40),
      "silty clay loam" = c(0, 40, 0, 27, 20, 27, 20, 40),
      "sandy clay" = c(45, 55, 45, 35, 65, 35),
      "silty clay" = c(0, 60, 0, 40, 20, 40),
      "clay" = c(0, 100, 0, 60, 20, 40, 45, 40, 45, 55)
    )
  } else {
    list(
      "coarse" = c(65, 0, 65, 18, 82, 18, 100, 0),
      "medium" = c(15, 0, 65, 0, 65, 18, 82, 18, 65, 35, 15, 35),
      "medium fine" = c(0, 0, 15, 0, 15, 35, 0, 35),
      "fine" = c(0, 35, 65, 35, 40, 60, 0, 60),
      "very fine" = c(0, 60, 40, 60, 0, 100)
    )
  }
  lv <- names(v)
  out <- lapply(lv, function(cl) {
    m <- matrix(v[[cl]], ncol = 2, byrow = TRUE)
    tibble::tibble(class = cl, vertex = seq_len(nrow(m)), sand = m[, 1],
      silt = 100 - m[, 1] - m[, 2], clay = m[, 2])
  })
  out <- dplyr::bind_rows(out)
  out$class <- factor(out$class, levels = lv)
  out
}

#' Isometric log-ratio coordinates for soil texture
#'
#' Sand, silt and clay are compositional data: they are constrained to sum
#' to 100, so ordinary statistics and models applied to them separately can
#' give impossible predictions. `soil_texture_ilr()` maps each composition
#' to two unconstrained isometric log-ratio (ilr) coordinates (Egozcue et
#' al. 2003) using the pivot balances
#' \deqn{z_1 = \sqrt{2/3}\,\ln\left(\frac{sand}{\sqrt{silt \cdot clay}}\right),
#'   \quad z_2 = \sqrt{1/2}\,\ln\left(\frac{silt}{clay}\right).}
#' Model or average the coordinates, then map back with
#' [soil_texture_ilr_inv()]. Zeros are replaced by `zero` before the
#' transform and the composition is closed again.
#'
#' @inheritParams soil_texture_class
#' @param names Names of the two new columns.
#' @param zero Replacement value (percentage) for zeros and values below it.
#'
#' @return `.data` with two new numeric columns.
#' @references Egozcue, J. J., Pawlowsky-Glahn, V., Mateu-Figueras, G. and
#'   Barcelo-Vidal, C. (2003). Isometric logratio transformations for
#'   compositional data analysis. *Mathematical Geology*, 35, 279-300.
#'   \doi{10.1023/A:1023818214614}
#' @family soil texture
#' @export
#' @examples
#' tex <- soil_texture_ilr(png_profiles, sand, silt, clay)
#' tex[c("sand", "silt", "clay", "ilr_1", "ilr_2")]
#'
#' # Round trip
#' soil_texture_ilr_inv(tex, ilr_1, ilr_2)[c("sand_ilr", "silt_ilr", "clay_ilr")]
soil_texture_ilr <- function(.data, sand, silt, clay, names = c("ilr_1", "ilr_2"), zero = 0.5) {
  check_data_frame(.data)
  check_number_positive(zero)
  if (!is.character(names) || length(names) != 2) {
    cli::cli_abort("{.arg names} must be a character vector of length 2.")
  }
  p <- texture_parts(.data, rlang::enquo(sand), rlang::enquo(silt), rlang::enquo(clay))
  z <- texture_ilr(p$sand, p$silt, p$clay, zero = zero)
  .data[[names[1]]] <- z[, 1]
  .data[[names[2]]] <- z[, 2]
  .data
}

#' Back-transform isometric log-ratio coordinates to sand, silt and clay
#'
#' The inverse of [soil_texture_ilr()].
#'
#' @param .data A data frame.
#' @param ilr_1,ilr_2 <[`tidy-select`][tidyselect::language]> The two ilr
#'   coordinate columns.
#' @param names Names of the three new columns (sand, silt and clay).
#'
#' @return `.data` with three new columns of percentages summing to 100.
#' @family soil texture
#' @export
#' @examples
#' df <- data.frame(ilr_1 = c(-1, 0, 1), ilr_2 = c(0.5, 0, -0.5))
#' soil_texture_ilr_inv(df, ilr_1, ilr_2)
soil_texture_ilr_inv <- function(.data, ilr_1, ilr_2,
                                 names = c("sand_ilr", "silt_ilr", "clay_ilr")) {
  check_data_frame(.data)
  if (!is.character(names) || length(names) != 3) {
    cli::cli_abort("{.arg names} must be a character vector of length 3.")
  }
  z1 <- numeric_col(rlang::enquo(ilr_1), .data, "ilr_1")
  z2 <- numeric_col(rlang::enquo(ilr_2), .data, "ilr_2")
  m <- texture_ilr_inv(z1, z2)
  for (i in 1:3) .data[[names[i]]] <- m[, i]
  .data
}

texture_ilr <- function(sand, silt, clay, zero = 0.5) {
  x <- cbind(sand, silt, clay)
  x[!is.na(x) & x < zero] <- zero
  x <- x / rowSums(x)
  z1 <- sqrt(2 / 3) * log(x[, 1] / sqrt(x[, 2] * x[, 3]))
  z2 <- sqrt(1 / 2) * log(x[, 2] / x[, 3])
  cbind(z1, z2)
}

texture_ilr_inv <- function(z1, z2) {
  # Centred log-ratio parts (summing to zero) from the pivot coordinates.
  a <- sqrt(2 / 3) * z1
  b <- (-a + sqrt(2) * z2) / 2
  cc <- (-a - sqrt(2) * z2) / 2
  e <- cbind(exp(a), exp(b), exp(cc))
  100 * e / rowSums(e)
}
