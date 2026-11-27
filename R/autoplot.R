# ggplot2 autoplot methods and the texture triangle ---------------------------

ternary_xy <- function(sand, silt, clay) {
  tot <- sand + silt + clay
  list(x = 100 * (silt + clay / 2) / tot, y = 100 * clay * sqrt(3) / 2 / tot)
}

#' Texture triangle plot
#'
#' Draws the USDA or FAO texture triangle with its class boundaries and
#' plots samples on it. Sand increases to the bottom left, silt to the
#' bottom right and clay to the top, as in the USDA Soil Survey Manual.
#'
#' @param .data A data frame of samples, or `NULL` for an empty triangle.
#' @param sand,silt,clay <[`tidy-select`][tidyselect::language]> Columns of
#'   percentages.
#' @param system `"usda"` or `"fao"`.
#' @param colour <[`tidy-select`][tidyselect::language]> Optional column
#'   mapped to point colour.
#' @param labels Label the classes?
#'
#' @return A ggplot object.
#' @family plots
#' @export
#' @examples
#' soil_texture_triangle(png_profiles, sand, silt, clay, colour = wrb_group)
#' soil_texture_triangle(NULL, system = "fao")
soil_texture_triangle <- function(.data, sand = NULL, silt = NULL, clay = NULL,
                                  system = c("usda", "fao"), colour = NULL, labels = TRUE) {
  system <- rlang::arg_match(system)
  polys <- soil_texture_polygons(system)
  pxy <- ternary_xy(polys$sand, polys$silt, polys$clay)
  polys$x <- pxy$x
  polys$y <- pxy$y
  h <- sqrt(3) / 2
  ticks <- seq(10, 90, by = 10)
  grid <- rbind(
    data.frame(x = ticks / 2, y = ticks * h, xend = 100 - ticks / 2, yend = ticks * h),
    data.frame(x = 100 - ticks, y = 0, xend = (100 - ticks) / 2, yend = (100 - ticks) * h),
    data.frame(x = ticks, y = 0, xend = 50 + ticks / 2, yend = (100 - ticks) * h)
  )
  frame <- data.frame(x = c(0, 100, 50, 0), y = c(0, 0, 100 * h, 0))
  tick_lab <- rbind(
    data.frame(x = ticks / 2 - 3, y = ticks * h, label = ticks, angle = 0),
    data.frame(x = 50 + (100 - ticks) / 2 + 3, y = ticks * h, label = 100 - ticks, angle = 0),
    data.frame(x = 100 - ticks, y = -3.5, label = ticks, angle = 0)
  )
  axis_lab <- data.frame(
    x = c(14, 86, 50), y = c(50, 50, -9), label = c("Clay (%)", "Silt (%)", "Sand (%)"),
    angle = c(60, -60, 0)
  )
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend,
      yend = .data$yend), data = grid, colour = "grey88", linewidth = 0.3) +
    ggplot2::geom_polygon(ggplot2::aes(x = .data$x, y = .data$y, group = .data$class),
      data = polys, fill = NA, colour = "grey30", linewidth = 0.4) +
    ggplot2::geom_path(ggplot2::aes(x = .data$x, y = .data$y), data = frame, linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      data = tick_lab, size = 2.6, colour = "grey30") +
    ggplot2::geom_text(ggplot2::aes(x = .data$x, y = .data$y, label = .data$label,
      angle = .data$angle), data = axis_lab, size = 3.4) +
    ggplot2::coord_equal(clip = "off") +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(10, 10, 10, 10))
  if (labels) {
    cent <- dplyr::summarise(dplyr::group_by(polys, .data$class),
      x = mean(.data$x), y = mean(.data$y))
    p <- p + ggplot2::geom_text(ggplot2::aes(x = .data$x, y = .data$y, label = .data$class),
      data = cent, size = 2.5, colour = "grey45", fontface = "italic")
  }
  if (!is.null(.data)) {
    check_data_frame(.data)
    parts <- texture_parts(.data, rlang::enquo(sand), rlang::enquo(silt), rlang::enquo(clay))
    xy <- ternary_xy(parts$sand, parts$silt, parts$clay)
    pts <- data.frame(x = xy$x, y = xy$y)
    col <- select_one(rlang::enquo(colour), .data, "colour", required = FALSE)
    if (!is.null(col)) {
      pts$.colour <- .data[[col]]
      p <- p + ggplot2::geom_point(ggplot2::aes(x = .data$x, y = .data$y, colour = .data$.colour),
        data = pts[!is.na(pts$x), ], size = 1.8, alpha = 0.8) +
        ggplot2::labs(colour = col)
    } else {
      p <- p + ggplot2::geom_point(ggplot2::aes(x = .data$x, y = .data$y),
        data = pts[!is.na(pts$x), ], size = 1.8, alpha = 0.8, colour = "#8C510A")
    }
  }
  p
}

#' Plot fitted water retention curves
#'
#' Draws the measured points and fitted van Genuchten curves of a
#' `<soil_retention>`, with suction on a log scale (points at zero suction
#' are drawn at the smallest positive suction).
#'
#' @param object A `<soil_retention>` from [soil_vg_fit()].
#' @param ... Not used.
#' @return A ggplot object.
#' @family plots
#' @export
#' @examples
#' fits <- soil_vg_fit(soil_retention_lab, h_cm, theta, .by = sample)
#' autoplot(fits)
autoplot.soil_retention <- function(object, ...) {
  meta <- soil_meta(object)
  obs <- meta$data
  hmin <- min(c(1, obs[[meta$h]][obs[[meta$h]] > 0]), na.rm = TRUE)
  hmax <- max(obs[[meta$h]], na.rm = TRUE)
  curves <- stats::predict(object, h = 10^seq(log10(hmin), log10(hmax), length.out = 120))
  lab <- function(d) {
    if (!length(meta$groups)) return(rep("all", nrow(d)))
    do.call(paste, c(unname(as.list(d[meta$groups])), sep = " / "))
  }
  obs$.curve <- lab(obs)
  curves$.curve <- lab(curves)
  obs$.h <- pmax(obs[[meta$h]], hmin)
  curves$.h <- curves[[meta$h]]
  ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(x = .data$.h, y = .data[[meta$theta]], colour = .data$.curve),
      data = obs) +
    ggplot2::geom_line(ggplot2::aes(x = .data$.h, y = .data$.pred, colour = .data$.curve),
      data = curves) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(x = paste0("Suction (", meta$h, ")"), y = "Volumetric water content",
      colour = NULL) +
    ggplot2::theme_minimal()
}

#' Plot soil profiles
#'
#' Draws one panel per profile with a property plotted against depth, each
#' layer as a bar from its upper to its lower depth.
#'
#' @param object A [soil_profile()].
#' @param var <[`tidy-select`][tidyselect::language]> Numeric property to
#'   plot.
#' @param max_profiles Maximum number of profiles to draw (the first ones).
#' @param ... Not used.
#' @return A ggplot object.
#' @family plots
#' @export
#' @examples
#' prof <- soil_profile(png_profiles, profile, top, bottom)
#' autoplot(prof, clay, max_profiles = 6)
autoplot.soil_profile <- function(object, var, max_profiles = 12, ...) {
  meta <- soil_meta(object)
  check_integerish(max_profiles)
  v <- select_one(rlang::enquo(var), object, "var")
  if (!is.numeric(object[[v]])) {
    cli::cli_abort("{.arg var} must be a numeric column.")
  }
  d <- as_plain_tibble(object)
  ids <- utils::head(unique(d[[meta$id]]), max_profiles)
  d <- d[d[[meta$id]] %in% ids & !is.na(d[[v]]), , drop = FALSE]
  ggplot2::ggplot(d) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = .data[[v]], ymin = .data[[meta$top]],
      ymax = .data[[meta$bottom]]), fill = "#BF812D", colour = "white", linewidth = 0.2) +
    ggplot2::scale_y_reverse() +
    ggplot2::facet_wrap(ggplot2::vars(.data[[meta$id]])) +
    ggplot2::labs(x = v, y = "Depth") +
    ggplot2::theme_minimal()
}

#' Plot soil spectra
#'
#' Draws every spectrum of a `<soil_spectra>` against wavelength.
#'
#' @param object A [soil_spectra()] table.
#' @param colour <[`tidy-select`][tidyselect::language]> Optional column
#'   mapped to line colour (for example clay content).
#' @param ... Not used.
#' @return A ggplot object.
#' @family plots
#' @export
#' @examples
#' spec <- soil_spectra(soil_nir, dplyr::starts_with("nm_"))
#' autoplot(spec, colour = clay)
#' autoplot(soil_savgol(spec, deriv = 1))
autoplot.soil_spectra <- function(object, colour = NULL, ...) {
  meta <- soil_meta(object)
  long <- tidy.soil_spectra(object)
  col <- select_one(rlang::enquo(colour), object, "colour", required = FALSE)
  ylab <- switch(meta$type, reflectance = "Reflectance", absorbance = "Absorbance",
    "Preprocessed value")
  if (is.null(col)) {
    p <- ggplot2::ggplot(long, ggplot2::aes(x = .data$wavelength, y = .data$value,
      group = .data$.row)) +
      ggplot2::geom_line(alpha = 0.4, colour = "#8C510A", linewidth = 0.3)
  } else {
    long$.colour <- rep(object[[col]], times = length(meta$wavelengths))
    p <- ggplot2::ggplot(long, ggplot2::aes(x = .data$wavelength, y = .data$value,
      group = .data$.row, colour = .data$.colour)) +
      ggplot2::geom_line(alpha = 0.6, linewidth = 0.3) +
      ggplot2::labs(colour = col)
  }
  p + ggplot2::labs(x = "Wavelength", y = ylab) + ggplot2::theme_minimal()
}
