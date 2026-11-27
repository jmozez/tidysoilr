# Particle-size analysis --------------------------------------------------

#' Particle-size distribution from hydrometer readings
#'
#' Converts hydrometer (Bouyoucos) readings into percentages of sand, silt
#' and clay (Gee and Or 2002). A standard ASTM 152H hydrometer reads grams
#' of soil per litre of suspension. Each reading is corrected with a blank
#' (dispersant solution at the same temperature) when one is given;
#' otherwise it is corrected for temperature by `temp_coef` g/L per degree
#' above 20 degrees C. The reading taken after `40 s` estimates the silt +
#' clay still in suspension, and the reading after `2 h` estimates the
#' clay.
#'
#' @param .data A data frame with one row per sample.
#' @param mass <[`tidy-select`][tidyselect::language]> Oven-dry mass of soil
#'   (g) dispersed in the 1 L cylinder.
#' @param reading_silt_clay,reading_clay <[`tidy-select`][tidyselect::language]>
#'   Hydrometer readings (g/L) for silt + clay (typically at 40 s) and for
#'   clay (typically at 2 h).
#' @param blank_silt_clay,blank_clay <[`tidy-select`][tidyselect::language]>
#'   Optional blank readings (g/L) taken at the same times.
#' @param temp_silt_clay,temp_clay <[`tidy-select`][tidyselect::language]>
#'   Optional suspension temperatures (degrees C), used when there is no
#'   blank.
#' @param temp_coef Temperature correction in g/L per degree C.
#' @param volume Suspension volume in litres.
#'
#' @return `.data` with new columns `sand`, `silt` and `clay` (percent of
#'   the oven-dry mass).
#' @references Gee, G. W. and Or, D. (2002). Particle-size analysis. In
#'   Dane, J. H. and Topp, G. C. (eds), *Methods of Soil Analysis, Part 4:
#'   Physical Methods*, 255-293. Soil Science Society of America, Madison.
#'   \doi{10.2136/sssabookser5.4.c12}
#' @family particle-size analysis
#' @export
#' @examples
#' soil_hydrometer(soil_hydrometer_lab, mass = mass, reading_silt_clay = r_40s,
#'   reading_clay = r_2h, blank_silt_clay = blank_40s, blank_clay = blank_2h)
soil_hydrometer <- function(.data, mass, reading_silt_clay, reading_clay,
                            blank_silt_clay = NULL, blank_clay = NULL,
                            temp_silt_clay = NULL, temp_clay = NULL,
                            temp_coef = 0.36, volume = 1) {
  check_data_frame(.data)
  check_number_positive(temp_coef, allow_zero = TRUE)
  check_number_positive(volume)
  m <- numeric_col(rlang::enquo(mass), .data, "mass", non_negative = TRUE)
  r1 <- numeric_col(rlang::enquo(reading_silt_clay), .data, "reading_silt_clay")
  r2 <- numeric_col(rlang::enquo(reading_clay), .data, "reading_clay")
  b1 <- numeric_col(rlang::enquo(blank_silt_clay), .data, "blank_silt_clay", required = FALSE)
  b2 <- numeric_col(rlang::enquo(blank_clay), .data, "blank_clay", required = FALSE)
  t1 <- numeric_col(rlang::enquo(temp_silt_clay), .data, "temp_silt_clay", required = FALSE)
  t2 <- numeric_col(rlang::enquo(temp_clay), .data, "temp_clay", required = FALSE)
  correct <- function(r, b, t) {
    if (!is.null(b)) return(r - b)
    if (!is.null(t)) return(r + temp_coef * (t - 20))
    r
  }
  if (is.null(b1) && is.null(t1)) {
    cli::cli_inform(c("i" = "No blank or temperature given; readings are used uncorrected."))
  }
  c1 <- correct(r1, b1, t1) * volume
  c2 <- correct(r2, b2, t2) * volume
  if (any(c2 > c1, na.rm = TRUE)) {
    cli::cli_warn("Some clay readings exceed the silt + clay readings; check the data.")
  }
  silt_clay <- pmin(pmax(100 * c1 / m, 0), 100)
  clay <- pmin(pmax(100 * c2 / m, 0), silt_clay)
  add_cols(.data, list(sand = 100 - silt_clay, silt = silt_clay - clay, clay = clay))
}

#' Particle-size distribution from pipette analysis
#'
#' Computes sand, silt and clay from the pipette method (Gee and Or 2002).
#' Sand is weighed after wet sieving; silt + clay and clay are estimated
#' from the oven-dry residue of aliquots pipetted from the sedimentation
#' cylinder at the times given by Stokes' law. The mass of dispersant in an
#' aliquot (from a blank) is subtracted from each residue.
#'
#' @param .data A data frame with one row per sample.
#' @param mass <[`tidy-select`][tidyselect::language]> Oven-dry mass of the
#'   sample (g).
#' @param sand_mass <[`tidy-select`][tidyselect::language]> Oven-dry mass of
#'   the sand retained on the sieve (g).
#' @param residue_silt_clay,residue_clay <[`tidy-select`][tidyselect::language]>
#'   Oven-dry residue (g) of the silt + clay and clay aliquots.
#' @param aliquot Aliquot volume (mL).
#' @param volume Cylinder volume (mL).
#' @param dispersant Mass of dispersant (g) in one aliquot, from a blank.
#' @param normalise Rescale the three fractions to sum to 100? The
#'   unscaled total is kept in `recovery`.
#'
#' @return `.data` with new columns `sand`, `silt`, `clay` (percent) and
#'   `recovery` (percent of the sample mass accounted for).
#' @inherit soil_hydrometer references
#' @family particle-size analysis
#' @export
#' @examples
#' df <- data.frame(id = 1:2, mass = c(10, 10), sand = c(4.1, 1.2),
#'   sc = c(0.152, 0.226), cl = c(0.081, 0.142))
#' soil_pipette(df, mass = mass, sand_mass = sand, residue_silt_clay = sc,
#'   residue_clay = cl, aliquot = 25, volume = 1000, dispersant = 0.0125)
soil_pipette <- function(.data, mass, sand_mass, residue_silt_clay, residue_clay,
                         aliquot = 25, volume = 1000, dispersant = 0, normalise = TRUE) {
  check_data_frame(.data)
  check_number_positive(aliquot)
  check_number_positive(volume)
  check_number_positive(dispersant, allow_zero = TRUE)
  m <- numeric_col(rlang::enquo(mass), .data, "mass", non_negative = TRUE)
  s <- numeric_col(rlang::enquo(sand_mass), .data, "sand_mass", non_negative = TRUE)
  r1 <- numeric_col(rlang::enquo(residue_silt_clay), .data, "residue_silt_clay")
  r2 <- numeric_col(rlang::enquo(residue_clay), .data, "residue_clay")
  k <- volume / aliquot
  silt_clay <- pmax(r1 - dispersant, 0) * k
  clay <- pmin(pmax(r2 - dispersant, 0) * k, silt_clay)
  sand <- 100 * s / m
  silt <- 100 * (silt_clay - clay) / m
  clay <- 100 * clay / m
  total <- sand + silt + clay
  if (normalise) {
    sand <- 100 * sand / total
    silt <- 100 * silt / total
    clay <- 100 * clay / total
  }
  add_cols(.data, list(sand = sand, silt = silt, clay = clay, recovery = total))
}
