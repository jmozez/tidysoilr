# tidysoilr design note

## Why another soil package?

R already has excellent soil software. tidysoilr does not replace it; it
fills the gaps that appear when soil data are handled in a tidyverse and
tidymodels workflow. The survey below covers the main tools (CRAN status as
of November 2026).

| Tool | What it does well | What is awkward in a tidy workflow |
|---|---|---|
| aqp (CRAN) | Profile collections, horizon logic, depth slicing, profile sketches | S4 `SoilProfileCollection` objects do not flow through dplyr verbs or tidymodels; GPL-3 |
| soilDB (CRAN) | Access to USDA KSSL, SSURGO, NASIS and SoilGrids web services | US-centred data access; results arrive as aqp objects or wide data frames |
| soiltexture (CRAN) | Very complete set of texture triangles and plots | Matrix input with fixed column names, base graphics only, AGPL-3 |
| SoilGrids access (geodata, soilDB, terra) | Download of global gridded soil layers | Raster handling, belongs with a spatial companion package |
| ithir (GitHub only) | Equal-area splines and digital soil mapping teaching code | Not on CRAN, older code, list outputs |
| mpspline2 (CRAN) | Mass-preserving splines | Returns nested lists; one property at a time unless reshaped |
| prospectr, simplerspec (GitHub) | Spectral preprocessing, Kennard-Stone, spectral workflows | Matrix-based functions; no recipes steps, so preprocessing cannot be tuned or applied consistently to new data in a workflow |
| plsmod (parsnip::pls) | PLS in tidymodels | Depends on mixOmics from Bioconductor, a heavy install |
| soilphysics, soilwater, SoilHyP (CRAN) | Water retention models, compaction, hydraulic functions | Vector interfaces, no grouped fitting, no broom methods |
| SoilTaxonomy, aqp classification helpers | Taxonomic lookups | Full WRB or Soil Taxonomy keys need morphological data rarely held in lab tables |

Common pain points:

1. **Depth-aware data lose their structure.** Plain data frames do not know
   which columns are the profile identifier and depths, and S4 profile
   collections do not work with dplyr.
2. **Compositional texture data are modelled as independent variables.**
   Sand, silt and clay sum to 100, yet models usually treat them as free
   predictors or responses.
3. **Spectral preprocessing sits outside the modelling pipeline.** SNV,
   Savitzky-Golay and friends are applied once to a matrix, so the same
   transformation is not guaranteed on new samples and cannot be tuned.
4. **Resampling ignores soil structure.** Layers of one profile, and nearby
   profiles, are correlated; random folds overstate accuracy.
5. **Everyday lab calculations are scattered.** Hydrometer sheets, bulk
   density, carbon stocks, base saturation and fertility ratings are often
   done in spreadsheets.

## Scope of version 0.1.0

* Tibble subclasses built on vctrs: `<soil_profile>` (identifier and depth
  columns), `<soil_spectra>` (wavelengths, measurement type and
  preprocessing history) and `<soil_retention>` (van Genuchten fits with
  their data). They survive dplyr verbs while the key columns are kept.
* Profiles: depth checks, slicing, depth-weighted means and conversion to
  aqp for those who need its tools.
* Texture: USDA (12 classes) and FAO (5 classes) triangles, polygons for
  plotting, isometric log-ratio coordinates and their inverse.
* Laboratory calculations: hydrometer and pipette particle-size analysis,
  bulk density, porosity, Saxton and Rawls (2006) and organic carbon based
  bulk density pedotransfer functions.
* Water retention: grouped van Genuchten fitting by bounded least squares,
  prediction and available water capacity.
* Carbon: equal-area quadratic splines (penalised least squares in a
  quadratic B-spline basis, giving the Bishop et al. 1999 solution),
  per-layer organic carbon stocks with coarse fragment correction, and
  fixed-depth totals with proration.
* Fertility: exchangeable bases, base saturation, ESP, ECEC and acid
  saturation, USDA reaction classes, threshold-based ratings with an
  editable table, and Parker's nutrient index.
* Spectroscopy: absorbance, SNV, MSC, Savitzky-Golay smoothing and
  derivatives, resampling and continuum removal.
* Soil quality: linear scoring functions, additive index and the principal
  component minimum data set of Andrews et al. (2002).
* tidymodels: recipes steps (SNV, Savitzky-Golay, texture log-ratios,
  texture class), profile-grouped and spatial block cross-validation and a
  Kennard-Stone split for rsample, yardstick metrics (bias-corrected SEP and
  texture accuracy allowing adjacent classes; RPIQ, RPD and Lin's CCC
  already exist in yardstick and are used alongside), and a base R PLS
  regression model for parsnip.
* broom `tidy()`, `glance()` and `augment()` methods and ggplot2
  `autoplot()` methods, plus a texture triangle.

## Design rules

* Data first: every user-facing function takes `.data` first and returns a
  tibble or the input with new columns, so it works with `|>`.
* Columns are chosen with tidy selection; errors and warnings use cli.
* Calculations are written in base R, stats and splines. sf and aqp are
  suggested, never required, so the package installs without GDAL or
  compilers.
* Units are stated in every argument description.

## Left for later

* Equivalent soil mass carbon stocks.
* Full taxonomic keys (WRB, Soil Taxonomy).
* Hydraulic conductivity functions (Mualem-van Genuchten) and more
  pedotransfer functions (Rosetta-style, tropical soil PTFs).
* Memory-based learning for large spectral libraries.
* SoilGrids access, which belongs with the spatial companion tidymapr.

## Data

* `png_profiles`: real WoSIS data for Papua New Guinea (ISRIC-WISE profiles
  under CC BY 3.0 and USDA NCSS profiles in the US public domain).
* `soil_retention_lab`, `soil_hydrometer_lab`, `soil_nir` and
  `soil_quality_trial`: clearly labelled simulated data, generated in
  `data-raw/make_data.R`.
