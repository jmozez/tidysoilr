# Builds the datasets shipped in data/.
#
# Run from the package root with `Rscript data-raw/make_data.R`. Downloads
# are cached in data-raw/cache/ (not committed), so later runs work offline.
#
# png_profiles: real soil profile data for Papua New Guinea from the ISRIC
#   World Soil Information Service (WoSIS), "WoSIS latest" web feature
#   service (https://maps.isric.org/), filtered to country "Papua New
#   Guinea". The PNG profiles come from the WISE (WD-WISE) dataset and are
#   distributed under CC BY 3.0 (licence field of every WoSIS layer).
# soil_retention_lab, soil_hydrometer, soil_nir, soil_quality_trial:
#   clearly labelled simulated data, generated below with a fixed seed.

library(xml2)
options(timeout = 600)

cache_dir <- file.path("data-raw", "cache")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

wfs <- "https://maps.isric.org/mapserv?map=/map/wosis_latest.map&SERVICE=WFS&VERSION=1.1.0"
png_bbox <- "-12,140.8,-0.5,160,urn:ogc:def:crs:EPSG::4326"

fetch_layer <- function(layer) {
  path <- file.path(cache_dir, paste0("wosis_", layer, ".xml"))
  if (!file.exists(path)) {
    url <- paste0(wfs, "&REQUEST=GetFeature&TYPENAME=ms:wosis_latest_", layer,
      "&BBOX=", png_bbox)
    for (attempt in 1:5) {
      ok <- tryCatch({
        utils::download.file(url, path, mode = "wb", quiet = TRUE)
        TRUE
      }, error = function(e) FALSE)
      if (ok) break
      Sys.sleep(5)
    }
    if (!ok) stop("Could not download WoSIS layer ", layer)
  }
  doc <- read_xml(path)
  feats <- xml_find_all(doc, paste0("//ms:wosis_latest_", layer), xml_ns(doc))
  if (!length(feats)) return(NULL)
  rows <- lapply(feats, function(f) {
    kids <- xml_children(f)
    nm <- xml_name(kids)
    keep <- !nm %in% c("boundedBy", "msGeometry")
    stats::setNames(as.list(xml_text(kids[keep])), nm[keep])
  })
  cols <- unique(unlist(lapply(rows, names)))
  out <- lapply(cols, function(cl) vapply(rows, function(r) r[[cl]] %||% NA_character_, ""))
  out <- as.data.frame(stats::setNames(out, cols), stringsAsFactors = FALSE)
  out[out$country_name == "Papua New Guinea", , drop = FALSE]
}

`%||%` <- function(x, y) if (is.null(x)) y else x

profiles <- fetch_layer("profiles")
cat("Profiles:", nrow(profiles), "\n")
print(table(profiles$dataset_code))

props <- c(
  clay = "clay", silt = "silt", sand = "sand", orgc = "orgc", bdfiod = "bdfiod",
  phaq = "phaq", phkc = "phkc", cecph7 = "cecph7", ecec = "ecec",
  nitkjd = "nitkjd", cfvo = "cfvo", wv0033 = "wv0033", wv1500 = "wv1500"
)
layers <- lapply(props, function(p) {
  x <- fetch_layer(p)
  if (is.null(x) || !nrow(x)) return(NULL)
  cat(p, nrow(x), "\n")
  x
})

# One row per layer: join every property on profile and layer ids.
key_cols <- c("profile_id", "layer_id", "layer_name", "upper_depth", "lower_depth")
hz <- NULL
for (p in names(layers)) {
  x <- layers[[p]]
  if (is.null(x)) next
  x <- x[, c(key_cols, "value_avg", "licence")]
  x$value_avg <- as.numeric(x$value_avg)
  x <- x[!duplicated(x[key_cols]), ]
  lic <- unique(x$licence)
  cat(p, "licence:", lic, "\n")
  names(x)[names(x) == "licence"] <- paste0("licence_", p)
  names(x)[names(x) == "value_avg"] <- p
  hz <- if (is.null(hz)) x else merge(hz, x, by = key_cols, all = TRUE)
}
hz$upper_depth <- as.numeric(hz$upper_depth)
hz$lower_depth <- as.numeric(hz$lower_depth)

sites <- profiles[, c("profile_id", "profile_code", "dataset_code", "year", "longitude",
  "latitude", "wrb_reference_soil_group", "fao_major_group", "usda_order_name")]
png <- merge(sites, hz, by = "profile_id")
png <- png[order(png$profile_code, png$upper_depth), ]
png <- png[!is.na(png$upper_depth) & !is.na(png$lower_depth) &
  png$lower_depth > png$upper_depth, ]

blank_na <- function(x) {
  x[!nzchar(x)] <- NA_character_
  x
}
lic_cols <- grep("^licence_", names(png), value = TRUE)
lic <- apply(png[lic_cols], 1, function(v) v[!is.na(v)][1])
lic <- ifelse(grepl("CC BY 3.0", lic), "CC BY 3.0",
  ifelse(grepl("Public Domain", lic), "US Public Domain", NA_character_))
png_profiles <- tibble::tibble(
  profile = png$profile_code,
  dataset = png$dataset_code,
  licence = lic,
  wosis_id = as.integer(png$profile_id),
  year = as.integer(png$year),
  lon = as.numeric(png$longitude),
  lat = as.numeric(png$latitude),
  wrb_group = blank_na(png$wrb_reference_soil_group),
  usda_order = blank_na(png$usda_order_name),
  horizon = blank_na(png$layer_name),
  top = png$upper_depth,
  bottom = png$lower_depth,
  sand = png$sand,
  silt = png$silt,
  clay = png$clay,
  coarse_frag = png$cfvo,
  oc = png$orgc,
  total_n = png$nitkjd,
  bulk_density = png$bdfiod,
  ph_water = png$phaq,
  ph_kcl = png$phkc,
  cec = png$cecph7,
  ecec = png$ecec,
  theta_33 = png$wv0033 / 100,
  theta_1500 = png$wv1500 / 100
)
# Drop rows with no analytical data at all.
vals <- png_profiles[, c("sand", "silt", "clay", "oc", "ph_water", "cec")]
png_profiles <- png_profiles[rowSums(!is.na(vals)) > 0, ]
cat("png_profiles:", nrow(png_profiles), "layers in",
  length(unique(png_profiles$profile)), "profiles\n")
print(table(png_profiles$dataset, png_profiles$licence, useNA = "ifany"))

usethis_save <- function(...) {
  objs <- list(...)
  for (nm in names(objs)) {
    assign(nm, objs[[nm]])
    save(list = nm, file = file.path("data", paste0(nm, ".rda")), compress = "xz")
  }
}

# Simulated data ------------------------------------------------------------
set.seed(2026)

# soil_retention_lab: retention curves simulated from the class-average van
# Genuchten parameters of Carsel and Parrish (1988), with measurement noise.
cp <- data.frame(
  texture = c("sand", "loamy sand", "sandy loam", "loam", "silt loam", "clay"),
  theta_r = c(0.045, 0.057, 0.065, 0.078, 0.067, 0.068),
  theta_s = c(0.43, 0.41, 0.41, 0.43, 0.45, 0.38),
  alpha = c(0.145, 0.124, 0.075, 0.036, 0.020, 0.008),
  n = c(2.68, 2.28, 1.89, 1.56, 1.41, 1.09)
)
h <- c(0, 10, 30, 60, 100, 330, 1000, 3000, 15000)
ret <- lapply(seq_len(nrow(cp)), function(i) {
  lapply(1:2, function(r) {
    m <- 1 - 1 / cp$n[i]
    th <- cp$theta_r[i] + (cp$theta_s[i] - cp$theta_r[i]) / (1 + (cp$alpha[i] * h)^cp$n[i])^m
    data.frame(
      sample = sprintf("S%02d", (i - 1) * 2 + r),
      texture = cp$texture[i],
      h_cm = h,
      h_kpa = round(h * 0.0980665, 2),
      theta = round(pmax(th + stats::rnorm(length(h), 0, 0.006), 0.01), 3)
    )
  })
})
soil_retention_lab <- tibble::as_tibble(do.call(rbind, unlist(ret, recursive = FALSE)))

# soil_hydrometer_lab: hydrometer readings simulated from known fractions.
n_h <- 24
true_clay <- round(stats::runif(n_h, 5, 65), 1)
true_sand <- round(stats::runif(n_h, 5, 95 - true_clay), 1)
mass <- 50
blank40 <- round(stats::rnorm(n_h, 5, 0.3), 1)
blank2h <- round(blank40 - stats::runif(n_h, 0, 1), 1)
soil_hydrometer_lab <- tibble::tibble(
  sample = sprintf("H%02d", seq_len(n_h)),
  mass = mass,
  r_40s = round((100 - true_sand) / 100 * mass + blank40 + stats::rnorm(n_h, 0, 0.3), 1),
  r_2h = round(true_clay / 100 * mass + blank2h + stats::rnorm(n_h, 0, 0.3), 1),
  blank_40s = blank40,
  blank_2h = blank2h,
  temp_c = round(stats::runif(n_h, 22, 28), 1)
)

# soil_nir: vis-NIR reflectance spectra simulated with absorption features
# at 1400 nm and 1900 nm (water and hydroxyl), 2200 nm (clay minerals) and
# a broad visible slope linked to organic carbon.
n_s <- 120
wl <- seq(400, 2480, by = 20)
profile_id <- rep(sprintf("P%02d", 1:30), each = 4)
site_lon <- stats::runif(30, 143, 151)
site_lat <- stats::runif(30, -10, -3)
clay <- round(pmin(pmax(stats::rnorm(n_s, 35, 15), 3), 80), 1)
oc <- round(pmax(stats::rlnorm(n_s, log(15), 0.6), 1), 1)
gauss <- function(x, mu, s) exp(-0.5 * ((x - mu) / s)^2)
spec <- t(vapply(seq_len(n_s), function(i) {
  base <- 0.15 + 0.35 * (1 - exp(-(wl - 350) / 600)) - 0.004 * oc[i] * gauss(wl, 600, 300)
  feat <- (0.0012 * clay[i] + 0.03) * gauss(wl, 1410, 30) +
    (0.0025 * clay[i] + 0.06) * gauss(wl, 1910, 45) +
    0.0020 * clay[i] * gauss(wl, 2205, 25)
  scatter <- stats::runif(1, 0.85, 1.15)
  pmin(pmax(scatter * (base - feat) + stats::rnorm(length(wl), 0, 0.002), 0.01), 0.99)
}, numeric(length(wl))))
colnames(spec) <- paste0("nm_", wl)
soil_nir <- tibble::as_tibble(data.frame(
  sample = sprintf("N%03d", seq_len(n_s)),
  profile = profile_id,
  depth = rep(c("0-15", "15-30", "30-60", "60-100"), 30),
  lon = round(rep(site_lon, each = 4), 4),
  lat = round(rep(site_lat, each = 4), 4),
  clay = clay,
  oc = oc,
  round(spec, 4),
  check.names = FALSE
))

# soil_quality_trial: indicators for plots under three land uses.
lu <- rep(c("forest", "food garden", "grassland"), each = 12)
mu <- list(
  oc = c(45, 22, 30), bulk_density = c(0.95, 1.25, 1.15), ph = c(5.6, 5.9, 5.7),
  avail_p = c(10, 14, 9), total_n = c(0.42, 0.21, 0.28), cec = c(28, 18, 22),
  agg_stability = c(78, 45, 62)
)
sdv <- c(oc = 8, bulk_density = 0.08, ph = 0.5, avail_p = 5, total_n = 0.07, cec = 4,
  agg_stability = 9)
k <- match(lu, c("forest", "food garden", "grassland"))
soil_quality_trial <- tibble::tibble(plot = sprintf("Q%02d", seq_along(lu)), land_use = lu)
for (v in names(mu)) {
  soil_quality_trial[[v]] <- round(pmax(mu[[v]][k] + stats::rnorm(length(lu), 0, sdv[[v]]), 0), 2)
}
soil_quality_trial$land_use <- factor(soil_quality_trial$land_use,
  levels = c("forest", "food garden", "grassland"))

dir.create("data", showWarnings = FALSE)
usethis_save(png_profiles = png_profiles, soil_retention_lab = soil_retention_lab,
  soil_hydrometer_lab = soil_hydrometer_lab, soil_nir = soil_nir,
  soil_quality_trial = soil_quality_trial)
