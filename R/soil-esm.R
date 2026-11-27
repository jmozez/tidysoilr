# Equivalent soil mass carbon stocks ---------------------------------------

esm_profile <- function(top, bottom, oc, bd, cf, ref_mass, method) {
  o <- order(top)
  top <- top[o]
  bottom <- bottom[o]
  mass <- bd[o] * (bottom - top) * (1 - cf[o]) * 100
  soc <- oc[o] / 1000 * mass
  bad <- anyNA(mass) || anyNA(soc) || any(abs(top[-1] - bottom[-length(bottom)]) > 1e-8)
  n_ref <- length(ref_mass)
  if (bad || !length(top)) {
    return(list(soc = rep(NA_real_, n_ref), depth = rep(NA_real_, n_ref),
                extrapolated = rep(NA, n_ref), ok = FALSE))
  }
  cm <- c(0, cumsum(mass))
  cc <- c(0, cumsum(soc))
  cd <- c(top[1], bottom)
  max_m <- cm[length(cm)]
  inside <- ref_mass <= max_m + 1e-8
  est <- rep(NA_real_, n_ref)
  if (method == "linear" || length(cm) < 3L) {
    est[inside] <- stats::approx(cm, cc, ref_mass[inside])$y
  } else {
    f <- stats::splinefun(cm, cc, method = "monoH.FC")
    est[inside] <- f(ref_mass[inside])
  }
  # Beyond the sampled mass, extend the deepest layer's concentration.
  last_conc <- soc[length(soc)] / mass[length(mass)]
  est[!inside] <- cc[length(cc)] + (ref_mass[!inside] - max_m) * last_conc
  last_rate <- mass[length(mass)] / (bottom[length(bottom)] - top[length(top)])
  depth <- rep(NA_real_, n_ref)
  depth[inside] <- stats::approx(cm, cd, ref_mass[inside])$y
  depth[!inside] <- cd[length(cd)] + (ref_mass[!inside] - max_m) / last_rate
  list(soc = est, depth = depth, extrapolated = !inside, ok = TRUE,
       mass_at = function(d) stats::approx(cd, cm, d, rule = 1)$y,
       soc_at = function(d) stats::approx(cd, cc, d, rule = 1)$y)
}

#' Organic carbon stocks on an equivalent soil mass basis
#'
#' @description
#' Fixed-depth carbon stocks are biased when bulk density differs between
#' the soils being compared (for example after tillage, compaction or a
#' change in land use), because a fixed depth then holds different amounts
#' of soil. The equivalent soil mass (ESM) approach compares the carbon in
#' the same mass of fine earth instead (Ellert and Bettany 1995; Wendt and
#' Hauser 2013), and is recommended by the FAO (2020) protocol for
#' measuring soil organic carbon change.
#'
#' For each profile, cumulative fine-earth mass (Mg/ha) and cumulative
#' organic carbon (Mg C/ha) are computed down the layers, and the carbon
#' held in each reference mass is interpolated, either linearly (constant
#' concentration within each layer) or with a monotone cubic spline of
#' cumulative carbon against cumulative mass (von Haden et al. 2020). The
#' default reference masses are the smallest masses found at `depths`
#' among the `reference` profiles, so no profile needs extrapolation; when
#' a reference mass is deeper than a profile was sampled, the deepest
#' layer's concentration is extended and the result is flagged.
#'
#' @param .data A [soil_profile()] or a data frame of layers, contiguous
#'   within each profile.
#' @param oc <[`tidy-select`][tidyselect::language]> Organic carbon.
#' @param bulk_density <[`tidy-select`][tidyselect::language]> Fine-earth
#'   bulk density (g/cm3).
#' @param coarse <[`tidy-select`][tidyselect::language]> Optional volumetric
#'   coarse fragment content.
#' @param depths Depths (cm) whose soil masses define the reference masses.
#' @param ref_mass Optional reference masses (cumulative fine-earth mass,
#'   Mg/ha), used instead of `depths`.
#' @param reference Optional profile identifiers that define the reference
#'   masses (for example the baseline sampling); all profiles by default.
#' @param method `"cubic"` (monotone cubic spline) or `"linear"`.
#' @param oc_units Units of `oc`: `"g_kg"` or `"percent"`.
#' @param coarse_units Units of `coarse`: `"percent"` or `"fraction"`.
#' @param id,top,bottom <[`tidy-select`][tidyselect::language]> Profile
#'   identifier and depth columns (cm). Taken from a `<soil_profile>` when
#'   not given.
#'
#' @return A tibble with one row per profile and reference mass: the
#'   identifier, `depth` (the reference depth, when `depths` is used),
#'   `ref_mass` (Mg/ha), `esm_depth` (the depth at which the profile
#'   reaches that mass, cm), `soc_esm` (Mg C/ha in the reference mass),
#'   `soc_fixed` (the fixed-depth stock to `depth`, for comparison) and
#'   `extrapolated`.
#' @references Ellert, B. H. and Bettany, J. R. (1995). Calculation of
#'   organic matter and nutrients stored in soils under contrasting
#'   management regimes. Canadian Journal of Soil Science, 75, 529-538.
#'   \doi{10.4141/cjss95-075}
#'
#'   Wendt, J. W. and Hauser, S. (2013). An equivalent soil mass procedure
#'   for monitoring soil organic carbon in multiple soil layers. European
#'   Journal of Soil Science, 64, 58-65. \doi{10.1111/ejss.12002}
#'
#'   von Haden, A. C., Yang, W. H. and DeLucia, E. H. (2020). Soils' dirty
#'   little secret: depth-based comparisons can be inadequate for
#'   quantifying changes in soil organic carbon and other mineral soil
#'   properties. Global Change Biology, 26, 3759-3770.
#'   \doi{10.1111/gcb.15124}
#'
#'   FAO (2020). A protocol for measurement, monitoring, reporting and
#'   verification of soil organic carbon in agricultural landscapes. GSOC-MRV
#'   Protocol. Food and Agriculture Organization of the United Nations,
#'   Rome. \doi{10.4060/cb0509en}
#' @family soil carbon
#' @export
#' @examples
#' # Three NCSS profiles from Papua New Guinea with measured bulk density
#' prof <- png_profiles[png_profiles$profile %in% c("82P0793", "82P0794", "82P0800"), ]
#' prof$coarse_frag[is.na(prof$coarse_frag)] <- 0
#' soil_carbon_esm(prof, oc, bulk_density, coarse = coarse_frag,
#'   depths = c(30, 100), id = profile, top = top, bottom = bottom
#' )
soil_carbon_esm <- function(.data, oc, bulk_density, coarse = NULL, depths = 30, ref_mass = NULL,
                            reference = NULL, method = c("cubic", "linear"),
                            oc_units = c("g_kg", "percent"),
                            coarse_units = c("percent", "fraction"),
                            id = NULL, top = NULL, bottom = NULL) {
  check_data_frame(.data)
  method <- rlang::arg_match(method)
  oc_units <- rlang::arg_match(oc_units)
  coarse_units <- rlang::arg_match(coarse_units)
  cols <- profile_cols(.data, rlang::enquo(id), rlang::enquo(top), rlang::enquo(bottom))
  x <- numeric_col(rlang::enquo(oc), .data, "oc", non_negative = TRUE)
  bd <- numeric_col(rlang::enquo(bulk_density), .data, "bulk_density", non_negative = TRUE)
  cf <- numeric_col(rlang::enquo(coarse), .data, "coarse", required = FALSE, non_negative = TRUE)
  if (oc_units == "percent") x <- x * 10
  cf <- if (is.null(cf)) rep(0, nrow(.data)) else if (coarse_units == "percent") cf / 100 else cf
  cf[is.na(cf)] <- 0
  if (any(cf > 1)) {
    cli::cli_abort("Coarse fragment fractions must not exceed 1; check {.arg coarse_units}.")
  }
  if (!is.null(ref_mass)) {
    if (!is.numeric(ref_mass) || anyNA(ref_mass) || any(ref_mass <= 0)) {
      cli::cli_abort("{.arg ref_mass} must be positive numbers (Mg/ha).")
    }
    depths <- NULL
  } else if (!is.numeric(depths) || !length(depths) || anyNA(depths) || any(depths <= 0)) {
    cli::cli_abort("{.arg depths} must be positive numbers (cm).")
  }
  ids <- unique(.data[[cols$id]])
  idx <- split(seq_len(nrow(.data)), factor(.data[[cols$id]], levels = ids))
  tp <- .data[[cols$top]]
  bt <- .data[[cols$bottom]]
  profiles <- lapply(idx, function(r) {
    esm_profile(tp[r], bt[r], x[r], bd[r], cf[r], numeric(), method)
  })
  bad <- !vapply(profiles, `[[`, logical(1), "ok")
  if (any(bad)) {
    cli::cli_warn(c(
      "{sum(bad)} profile{?s} had missing values or gaps between layers and gave {.val NA}.",
      "i" = "Affected: {.val {as.character(ids[bad])}}."
    ))
  }
  if (is.null(ref_mass)) {
    ref_ids <- reference %||% ids
    if (!all(ref_ids %in% ids)) {
      cli::cli_abort("Every {.arg reference} profile must be in {.field {cols$id}}.")
    }
    use <- which(ids %in% ref_ids & !bad)
    if (!length(use)) cli::cli_abort("No complete reference profile.")
    masses <- vapply(profiles[use], function(p) p$mass_at(depths), numeric(length(depths)))
    masses <- matrix(masses, nrow = length(depths))
    ref_mass <- apply(masses, 1, function(m) if (all(is.na(m))) NA_real_ else min(m, na.rm = TRUE))
    if (anyNA(ref_mass)) {
      cli::cli_abort("No reference profile reaches {max(depths)} cm.")
    }
  }
  out <- lapply(seq_along(ids), function(i) {
    r <- idx[[i]]
    p <- esm_profile(tp[r], bt[r], x[r], bd[r], cf[r], ref_mass, method)
    fixed <- if (is.null(depths) || !p$ok) NA_real_ else p$soc_at(depths)
    tibble::tibble(
      !!cols$id := ids[[i]],
      depth = if (is.null(depths)) NA_real_ else depths,
      ref_mass = ref_mass, esm_depth = p$depth, soc_esm = p$soc, soc_fixed = fixed,
      extrapolated = p$extrapolated
    )
  })
  dplyr::bind_rows(out)
}
