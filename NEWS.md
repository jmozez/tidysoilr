# tidysoilr (development version)

* New `soil_carbon_esm()` computes organic carbon stocks on an equivalent
  soil mass basis (Ellert and Bettany 1995; Wendt and Hauser 2013), with
  linear or monotone cubic spline interpolation of cumulative carbon
  against cumulative soil mass (von Haden et al. 2020), as recommended by
  the FAO GSOC-MRV protocol for monitoring carbon change.
* New `soil_vg_conductivity()` predicts unsaturated hydraulic conductivity
  from van Genuchten parameters with the Mualem (1976) model.
* New `soil_infiltration_fit()` fits the Philip, Horton and Kostiakov
  infiltration equations to infiltrometer readings and compares them by
  AIC.

# tidysoilr 0.1.0

* First version. Tidy tools for soil texture, particle-size analysis,
  physics and pedotransfer functions, water retention, equal-area depth
  splines and carbon stocks, fertility, spectroscopy preprocessing and soil
  quality indices, with recipes steps, rsample resampling, yardstick
  metrics, a parsnip PLS model, broom methods and ggplot2 plots.
* Datasets: `png_profiles` (real WoSIS profiles from Papua New Guinea) and
  four clearly labelled simulated datasets.
