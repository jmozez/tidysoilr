# Bulk density, porosity and pedotransfer functions -----------------------

#' Bulk density from core samples
#'
#' Computes dry bulk density from the oven-dry mass and volume of a core
#' and, when the mass of coarse fragments (> 2 mm) is known, the fine-earth
#' bulk density and the volumetric coarse fragment content. Fine-earth bulk
#' density is the value needed for carbon stocks (Poeplau et al. 2017).
#'
#' @param .data A data frame with one row per core.
#' @param dry_mass <[`tidy-select`][tidyselect::language]> Oven-dry mass of
#'   the whole sample (g).
#' @param volume <[`tidy-select`][tidyselect::language]> Core volume
#'   (cm3).
#' @param coarse_mass <[`tidy-select`][tidyselect::language]> Optional mass of
#'   coarse fragments (g).
#' @param coarse_density Density of the coarse fragments (g/cm3).
#'
#' @return `.data` with `bulk_density` (g/cm3) and, when `coarse_mass` is
#'   given, `fine_earth_bd` (g/cm3) and `coarse_vol` (volume fraction).
#' @references Poeplau, C., Vos, C. and Don, A. (2017). Soil organic carbon
#'   stocks are systematically overestimated by misuse of the parameters
#'   bulk density and rock fragment content. *SOIL*, 3, 61-66.
#'   \doi{10.5194/soil-3-61-2017}
#' @family soil physics
#' @export
#' @examples
#' cores <- data.frame(core = 1:3, mass = c(130, 142, 118), vol = 100,
#'   stones = c(0, 12, 5))
#' soil_bulk_density(cores, dry_mass = mass, volume = vol, coarse_mass = stones)
soil_bulk_density <- function(.data, dry_mass, volume, coarse_mass = NULL, coarse_density = 2.65) {
  check_data_frame(.data)
  check_number_positive(coarse_density)
  m <- numeric_col(rlang::enquo(dry_mass), .data, "dry_mass", non_negative = TRUE)
  v <- numeric_col(rlang::enquo(volume), .data, "volume", non_negative = TRUE)
  cm <- numeric_col(rlang::enquo(coarse_mass), .data, "coarse_mass", required = FALSE,
    non_negative = TRUE)
  out <- list(bulk_density = m / v)
  if (!is.null(cm)) {
    if (any(cm > m, na.rm = TRUE)) {
      cli::cli_abort("{.arg coarse_mass} must not exceed {.arg dry_mass}.")
    }
    cv <- cm / coarse_density
    out$fine_earth_bd <- (m - cm) / (v - cv)
    out$coarse_vol <- cv / v
  }
  add_cols(.data, out)
}

#' Total porosity, air-filled porosity and water-filled pore space
#'
#' Total porosity is \eqn{1 - \rho_b / \rho_s}. When volumetric water
#' content is given, air-filled porosity (porosity minus water content) and
#' water-filled pore space (water content divided by porosity) are added.
#'
#' @param .data A data frame.
#' @param bulk_density <[`tidy-select`][tidyselect::language]> Bulk density
#'   (g/cm3).
#' @param particle_density Particle density (g/cm3); 2.65 is the usual
#'   value for mineral soils.
#' @param water <[`tidy-select`][tidyselect::language]> Optional volumetric
#'   water content (m3/m3).
#'
#' @return `.data` with `porosity` and, when `water` is given, `air_filled`
#'   and `wfps` (all as fractions).
#' @family soil physics
#' @export
#' @examples
#' df <- data.frame(bd = c(1.1, 1.3, 1.5), theta = c(0.30, 0.35, 0.25))
#' soil_porosity(df, bulk_density = bd, water = theta)
soil_porosity <- function(.data, bulk_density, particle_density = 2.65, water = NULL) {
  check_data_frame(.data)
  check_number_positive(particle_density)
  bd <- numeric_col(rlang::enquo(bulk_density), .data, "bulk_density", non_negative = TRUE)
  w <- numeric_col(rlang::enquo(water), .data, "water", required = FALSE, non_negative = TRUE)
  if (any(bd > particle_density, na.rm = TRUE)) {
    cli::cli_warn("Some bulk densities exceed {.arg particle_density}; porosity will be negative.")
  }
  phi <- 1 - bd / particle_density
  out <- list(porosity = phi)
  if (!is.null(w)) {
    out$air_filled <- phi - w
    out$wfps <- w / phi
  }
  add_cols(.data, out)
}

#' Saxton and Rawls (2006) pedotransfer functions
#'
#' Estimates soil water characteristics and saturated hydraulic
#' conductivity from sand, clay and organic matter with the equations of
#' Saxton and Rawls (2006), fitted to about 2000 USDA soil samples. Results
#' are for soils at normal density without gravel.
#'
#' @param .data A data frame.
#' @param sand,clay <[`tidy-select`][tidyselect::language]> Sand and clay
#'   (percent).
#' @param om <[`tidy-select`][tidyselect::language]> Organic matter
#'   (percent by weight). Use 1.724 times organic carbon if only carbon is
#'   known. The equations were fitted for organic matter up to 8 percent.
#'
#' @return `.data` with new columns `theta_1500` (wilting point, 1500 kPa),
#'   `theta_33` (field capacity, 33 kPa), `theta_sat` (saturation), `paw`
#'   (plant-available water, `theta_33 - theta_1500`), all in m3/m3,
#'   `air_entry` (kPa), `lambda` (slope of the log tension curve),
#'   `bulk_density` (g/cm3) and `ksat` (mm/h).
#' @references Saxton, K. E. and Rawls, W. J. (2006). Soil water
#'   characteristic estimates by texture and organic matter for hydrologic
#'   solutions. *Soil Science Society of America Journal*, 70, 1569-1578.
#'   \doi{10.2136/sssaj2005.0117}
#' @family soil physics
#' @family pedotransfer functions
#' @export
#' @examples
#' df <- data.frame(sand = c(80, 40, 10), clay = c(5, 20, 45), om = c(1, 2.5, 3))
#' soil_ptf_saxton(df, sand, clay, om)
soil_ptf_saxton <- function(.data, sand, clay, om) {
  check_data_frame(.data)
  s <- numeric_col(rlang::enquo(sand), .data, "sand", non_negative = TRUE) / 100
  c <- numeric_col(rlang::enquo(clay), .data, "clay", non_negative = TRUE) / 100
  o <- numeric_col(rlang::enquo(om), .data, "om", non_negative = TRUE)
  if (any(o > 8, na.rm = TRUE)) {
    cli::cli_warn("Organic matter above 8% is outside the range used to fit the equations.")
  }
  t1500 <- -0.024 * s + 0.487 * c + 0.006 * o + 0.005 * s * o - 0.013 * c * o +
    0.068 * s * c + 0.031
  th1500 <- t1500 + (0.14 * t1500 - 0.02)
  t33 <- -0.251 * s + 0.195 * c + 0.011 * o + 0.006 * s * o - 0.027 * c * o +
    0.452 * s * c + 0.299
  th33 <- t33 + (1.283 * t33^2 - 0.374 * t33 - 0.015)
  ts33 <- 0.278 * s + 0.034 * c + 0.022 * o - 0.018 * s * o - 0.027 * c * o -
    0.584 * s * c + 0.078
  ths33 <- ts33 + (0.636 * ts33 - 0.107)
  pe <- -21.67 * s - 27.93 * c - 81.97 * ths33 + 71.12 * s * ths33 + 8.29 * c * ths33 +
    14.05 * s * c + 27.16
  pe <- pe + (0.02 * pe^2 - 0.113 * pe - 0.70)
  ths <- th33 + ths33 - 0.097 * s + 0.043
  bd <- (1 - ths) * 2.65
  b <- (log(1500) - log(33)) / (log(th33) - log(th1500))
  lambda <- 1 / b
  ks <- 1930 * pmax(ths - th33, 0)^(3 - lambda)
  add_cols(.data, list(
    theta_1500 = th1500, theta_33 = th33, theta_sat = ths, paw = th33 - th1500,
    air_entry = pe, lambda = lambda, bulk_density = bd, ksat = ks
  ))
}

#' Bulk density pedotransfer functions from organic carbon
#'
#' Bulk density is often missing from legacy profile data but is needed for
#' carbon stocks. These widely used pedotransfer functions predict it from
#' organic carbon (OC, percent):
#'
#' * `"manrique_jones"`: \eqn{1.660 - 0.318\sqrt{OC}} (Manrique and Jones
#'   1991).
#' * `"alexander"`: \eqn{1.72 - 0.294\sqrt{OC}} (Alexander 1980).
#' * `"honeysett_ratkowsky"`: \eqn{1 / (0.548 + 0.0588\,OC)} (Honeysett and
#'   Ratkowsky 1989).
#' * `"adams"`: \eqn{100 / (OM / 0.244 + (100 - OM) / \rho_m)}, with
#'   \eqn{OM = 1.724\,OC} and mineral bulk density \eqn{\rho_m} (Adams 1973).
#'
#' Their performance varies widely between regions (De Vos et al. 2005);
#' check them against local measurements where possible.
#'
#' @param .data A data frame.
#' @param oc <[`tidy-select`][tidyselect::language]> Organic carbon
#'   (percent, g/100 g). Divide g/kg by 10.
#' @param method One of the methods above.
#' @param mineral_bd Mineral bulk density (g/cm3) for `method = "adams"`.
#' @param .name Name of the new column.
#'
#' @return `.data` with a new column of predicted bulk density (g/cm3).
#' @references De Vos, B., Van Meirvenne, M., Quataert, P., Deckers, J. and
#'   Muys, B. (2005). Predictive quality of pedotransfer functions for
#'   estimating bulk density of forest soils. *Soil Science Society of
#'   America Journal*, 69, 500-510. \doi{10.2136/sssaj2005.0500}
#' @family soil physics
#' @family pedotransfer functions
#' @export
#' @examples
#' png_profiles |>
#'   dplyr::mutate(oc_pct = oc / 10) |>
#'   soil_ptf_bd(oc_pct) |>
#'   dplyr::select(profile, top, oc_pct, bulk_density, bd_ptf)
soil_ptf_bd <- function(.data, oc, method = c("manrique_jones", "alexander",
                                              "honeysett_ratkowsky", "adams"),
                        mineral_bd = 1.64, .name = "bd_ptf") {
  check_data_frame(.data)
  method <- rlang::arg_match(method)
  check_number_positive(mineral_bd)
  check_name(.name)
  x <- numeric_col(rlang::enquo(oc), .data, "oc", non_negative = TRUE)
  if (any(x > 60, na.rm = TRUE)) {
    cli::cli_warn("Organic carbon above 60% looks like g/kg; the functions expect percent.")
  }
  .data[[.name]] <- switch(method,
    manrique_jones = 1.660 - 0.318 * sqrt(x),
    alexander = 1.72 - 0.294 * sqrt(x),
    honeysett_ratkowsky = 1 / (0.548 + 0.0588 * x),
    adams = {
      om <- 1.724 * x
      100 / (om / 0.244 + (100 - om) / mineral_bd)
    }
  )
  .data
}
