#' tidysoilr: tidy tools for soil analysis
#'
#' tidysoilr brings everyday soil laboratory and profile calculations into
#' the tidyverse and tidymodels. Functions are grouped by topic:
#'
#' * **Classes**: [soil_profile()], [soil_spectra()] and the
#'   `<soil_retention>` tibble returned by [soil_vg_fit()].
#' * **Profiles**: [soil_profile_check()], [soil_slice()],
#'   [soil_depth_mean()].
#' * **Texture**: [soil_texture_class()], [soil_texture_polygons()],
#'   [soil_texture_ilr()], [soil_texture_ilr_inv()].
#' * **Particle-size analysis**: [soil_hydrometer()], [soil_pipette()].
#' * **Physics and pedotransfer functions**: [soil_bulk_density()],
#'   [soil_porosity()], [soil_ptf_saxton()], [soil_ptf_bd()].
#' * **Water retention**: [soil_vg_fit()], [soil_vg_theta()], [soil_awc()].
#' * **Carbon and depth harmonisation**: [soil_spline()],
#'   [soil_carbon_stock()], [soil_carbon_total()].
#' * **Fertility**: [soil_base_saturation()], [soil_ph_class()],
#'   [soil_fertility_rating()], [soil_fertility_thresholds()],
#'   [soil_nutrient_index()].
#' * **Spectroscopy**: [soil_absorbance()], [soil_snv()], [soil_msc()],
#'   [soil_savgol()], [soil_resample()], [soil_continuum_removal()].
#' * **Soil quality**: [soil_quality_index()], [sqi_more()], [sqi_less()],
#'   [sqi_optimum()], [soil_mds()].
#' * **tidymodels**: [step_soil_snv()], [step_soil_savgol()],
#'   [step_soil_ilr()], [step_soil_texture()], [soil_profile_cv()],
#'   [soil_spatial_block_cv()], [soil_kennard_stone()], [sep()],
#'   [texture_adjacent_accuracy()], [soil_pls()].
#' * **Plots**: [soil_texture_triangle()] and `autoplot()` methods.
#'
#' All calculations are written in base R and the recommended packages, so
#' no compiled spatial or chemometric libraries are needed.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom rlang .data := %||%
#' @importFrom generics tidy glance augment required_pkgs
#' @importFrom ggplot2 autoplot
#' @importFrom recipes prep bake
#' @importFrom pillar tbl_sum
#' @importFrom vctrs vec_restore vec_ptype2 vec_cast
#' @importFrom dplyr dplyr_reconstruct
#' @importFrom stats predict update
## usethis namespace: end
NULL

#' @export
generics::tidy

#' @export
generics::glance

#' @export
generics::augment

#' @export
generics::required_pkgs

#' @export
ggplot2::autoplot

#' @export
recipes::prep

#' @export
recipes::bake

.onLoad <- function(libname, pkgname) {
  current <- parsnip::get_model_env()
  if (!"soil_pls" %in% current$models) {
    make_soil_pls()
  }
  invisible()
}
