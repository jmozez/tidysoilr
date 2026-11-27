#' Soil profiles from Papua New Guinea (WoSIS)
#'
#' Real laboratory data for 150 layers of 32 soil profiles described in
#' Papua New Guinea between 1967 and 1983, extracted from the ISRIC World
#' Soil Information Service (WoSIS). The profiles come from two source
#' datasets: the ISRIC-WISE harmonised profile database (`dataset ==
#' "WD-WISE"`, 24 profiles, licence CC BY 3.0) and the USDA-NRCS National
#' Cooperative Soil Survey characterisation database (`dataset ==
#' "US-NCSS"`, 8 profiles, US public domain). Values are the standardised
#' WoSIS values; profiles without any analytical data are omitted. Positions
#' of the WISE profiles are uncertain by more than 10 km. See
#' `data-raw/make_data.R` for the download and processing code.
#'
#' @format A tibble with 150 rows and 25 columns:
#' \describe{
#'   \item{profile}{Profile code.}
#'   \item{dataset}{WoSIS source dataset code.}
#'   \item{licence}{Licence of the source dataset.}
#'   \item{wosis_id}{WoSIS profile identifier.}
#'   \item{year}{Year of description.}
#'   \item{lon, lat}{Longitude and latitude (decimal degrees, WGS 84).}
#'   \item{wrb_group}{WRB reference soil group, where classified.}
#'   \item{usda_order}{USDA Soil Taxonomy order, where classified.}
#'   \item{horizon}{Horizon designation.}
#'   \item{top, bottom}{Upper and lower depth of the layer (cm).}
#'   \item{sand, silt, clay}{Particle-size fractions (percent).}
#'   \item{coarse_frag}{Coarse fragments (volume percent).}
#'   \item{oc}{Organic carbon (g/kg).}
#'   \item{total_n}{Total nitrogen (g/kg).}
#'   \item{bulk_density}{Fine-earth bulk density, oven dry (g/cm3).}
#'   \item{ph_water, ph_kcl}{pH in water and in KCl.}
#'   \item{cec}{Cation exchange capacity at pH 7 (cmol(+)/kg).}
#'   \item{ecec}{Effective cation exchange capacity (cmol(+)/kg).}
#'   \item{theta_33, theta_1500}{Volumetric water content at 33 and 1500 kPa
#'     (m3/m3).}
#' }
#' @source ISRIC World Soil Information, WoSIS latest web feature service (see
#'   <https://isric.org/explore/wosis>), accessed 27 November 2026. Cite as: Batjes,
#'   N. H., Calisto, L. and de Sousa, L. M. (2024). Providing
#'   quality-assessed and standardised soil data to support global mapping
#'   and modelling (WoSIS snapshot 2023). *Earth System Science Data*, 16,
#'   4735-4765. \doi{10.5194/essd-16-4735-2024}. WISE data: Batjes, N. H.
#'   (2009). Harmonized soil profile data for applications at global and
#'   continental scales: updates to the WISE database. *Soil Use and
#'   Management*, 25, 124-127. \doi{10.1111/j.1475-2743.2009.00202.x}. NCSS
#'   data: Soil Survey Staff, National Cooperative Soil Survey Soil
#'   Characterization Database, USDA Natural Resources Conservation Service.
#' @family datasets
#' @examples
#' png_profiles
#' soil_profile(png_profiles, profile, top, bottom)
"png_profiles"

#' Laboratory water retention curves (simulated)
#'
#' Simulated water retention measurements for 12 samples, two for each of
#' six texture classes, generated from the class-average van Genuchten
#' parameters of Carsel and Parrish (1988) with random measurement error
#' (standard deviation 0.006 m3/m3). These are not real measurements.
#'
#' @format A tibble with 108 rows and 5 columns:
#' \describe{
#'   \item{sample}{Sample identifier.}
#'   \item{texture}{USDA texture class used to simulate the curve.}
#'   \item{h_cm}{Matric suction (cm of water).}
#'   \item{h_kpa}{Matric suction (kPa).}
#'   \item{theta}{Volumetric water content (m3/m3).}
#' }
#' @source Simulated in `data-raw/make_data.R`. Parameters from Carsel, R.
#'   F. and Parrish, R. S. (1988). Developing joint probability distributions
#'   of soil water retention characteristics. *Water Resources Research*, 24,
#'   755-769. \doi{10.1029/WR024i005p00755}
#' @family datasets
#' @examples
#' soil_retention_lab
"soil_retention_lab"

#' Hydrometer particle-size readings (simulated)
#'
#' Simulated ASTM 152H hydrometer readings for 24 samples of 50 g, with
#' blank readings, for use with [soil_hydrometer()]. These are not real
#' measurements.
#'
#' @format A tibble with 24 rows and 7 columns:
#' \describe{
#'   \item{sample}{Sample identifier.}
#'   \item{mass}{Oven-dry soil mass (g).}
#'   \item{r_40s, r_2h}{Hydrometer readings at 40 seconds and 2 hours (g/L).}
#'   \item{blank_40s, blank_2h}{Blank readings at the same times (g/L).}
#'   \item{temp_c}{Suspension temperature (degrees C).}
#' }
#' @source Simulated in `data-raw/make_data.R`.
#' @family datasets
#' @examples
#' soil_hydrometer_lab
"soil_hydrometer_lab"

#' Visible and near-infrared soil spectra (simulated)
#'
#' Simulated diffuse reflectance spectra (400 to 2480 nm every 20 nm) of 120
#' samples from 30 profiles, with absorption features near 1400 and 1900 nm
#' (water and hydroxyl) and 2200 nm (clay minerals) scaled by clay content,
#' a visible-range darkening scaled by organic carbon, multiplicative
#' scatter and noise. Coordinates are random points in Papua New Guinea.
#' These are not real measurements; they are meant for demonstrating
#' preprocessing and calibration workflows.
#'
#' @format A tibble with 120 rows and 112 columns:
#' \describe{
#'   \item{sample}{Sample identifier.}
#'   \item{profile}{Profile identifier (four samples per profile).}
#'   \item{depth}{Depth interval (cm).}
#'   \item{lon, lat}{Simulated coordinates (decimal degrees).}
#'   \item{clay}{Clay (percent).}
#'   \item{oc}{Organic carbon (g/kg).}
#'   \item{nm_400, ..., nm_2480}{Reflectance (0 to 1) at each wavelength (nm).}
#' }
#' @source Simulated in `data-raw/make_data.R`.
#' @family datasets
#' @examples
#' soil_spectra(soil_nir, dplyr::starts_with("nm_"))
"soil_nir"

#' Soil quality indicators under three land uses (simulated)
#'
#' Simulated topsoil indicators for 36 plots, 12 each under forest, food
#' gardens and grassland, for use with [soil_quality_index()] and
#' [soil_mds()]. These are not real measurements.
#'
#' @format A tibble with 36 rows and 9 columns:
#' \describe{
#'   \item{plot}{Plot identifier.}
#'   \item{land_use}{Land use (factor).}
#'   \item{oc}{Organic carbon (g/kg).}
#'   \item{bulk_density}{Bulk density (g/cm3).}
#'   \item{ph}{pH in water.}
#'   \item{avail_p}{Available phosphorus (mg/kg).}
#'   \item{total_n}{Total nitrogen (percent).}
#'   \item{cec}{Cation exchange capacity (cmol(+)/kg).}
#'   \item{agg_stability}{Water-stable aggregates (percent).}
#' }
#' @source Simulated in `data-raw/make_data.R`.
#' @family datasets
#' @examples
#' soil_quality_trial
"soil_quality_trial"
