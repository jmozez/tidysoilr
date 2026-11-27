# rsample resampling schemes ------------------------------------------------

# Assign groups of given sizes to v folds, largest first, emptiest fold next.
balance_folds <- function(sizes, v) {
  ord <- names(sizes)[sample.int(length(sizes))]
  ord <- ord[order(-sizes[ord])]
  load <- rep(0, v)
  fold_of <- stats::setNames(integer(length(ord)), ord)
  for (g in ord) {
    f <- which(load == min(load))
    f <- if (length(f) > 1) sample(f, 1) else f
    fold_of[g] <- f
    load[f] <- load[f] + sizes[[g]]
  }
  fold_of
}

make_fold_rset <- function(fold, v, data, attrib, subclass) {
  splits <- lapply(seq_len(v), function(f) {
    rsample::make_splits(list(analysis = which(fold != f), assessment = which(fold == f)),
      data = data)
  })
  ids <- tibble::tibble(id = sprintf(paste0("Fold%0", nchar(v), "d"), seq_len(v)))
  rsample::new_rset(splits = splits, ids = ids, attrib = attrib,
    subclass = c(subclass, "rset"))
}

#' Profile-grouped cross-validation
#'
#' Splits layers into `v` folds so that all layers of a profile fall in the
#' same fold. Layers of one profile are strongly correlated, so ordinary
#' random folds leak information between analysis and assessment sets and
#' overstate the accuracy of soil property models. Profiles are assigned to
#' folds so that folds hold similar numbers of layers.
#'
#' @param .data A [soil_profile()] or data frame.
#' @param id <[`tidy-select`][tidyselect::language]> Profile identifier.
#'   Taken from a `<soil_profile>` when not given.
#' @param v Number of folds.
#'
#' @return An `rset` of class `soil_profile_cv` with columns `splits` and
#'   `id`, usable with `tune::fit_resamples()`.
#' @family resampling
#' @export
#' @examples
#' set.seed(1)
#' prof <- soil_profile(png_profiles, profile, top, bottom)
#' folds <- soil_profile_cv(prof, v = 5)
#' folds
#' # No profile appears in both sets
#' a <- rsample::analysis(folds$splits[[1]])
#' b <- rsample::assessment(folds$splits[[1]])
#' intersect(a$profile, b$profile)
soil_profile_cv <- function(.data, id = NULL, v = 5) {
  check_data_frame(.data)
  check_integerish(v, min = 2)
  meta <- if (is_soil_profile(.data)) soil_meta(.data) else list()
  id_col <- select_one(rlang::enquo(id), .data, "id", required = is.null(meta$id)) %||% meta$id
  g <- as.character(.data[[id_col]])
  if (anyNA(g)) {
    cli::cli_abort("Profile identifier {.field {id_col}} must not contain missing values.")
  }
  sizes <- table(g)
  if (length(sizes) < v) {
    cli::cli_abort("Only {length(sizes)} profile{?s} for {v} folds.")
  }
  fold_of <- balance_folds(sizes, v)
  make_fold_rset(unname(fold_of[g]), v, as_plain_tibble(.data),
    attrib = list(v = v, group = id_col), subclass = "soil_profile_cv")
}

#' @export
print.soil_profile_cv <- function(x, ...) {
  cat("#  Profile-grouped cross-validation: ", attr(x, "v"), " folds by ", attr(x, "group"),
    "\n", sep = "")
  class(x) <- c("tbl_df", "tbl", "data.frame")
  print(x, ...)
  invisible(x)
}

#' Spatial block cross-validation
#'
#' Splits samples into `v` folds by square spatial blocks, so that each
#' assessment set is spatially separated from its analysis set (Roberts et
#' al. 2017). This gives honest estimates of how well a digital soil
#' mapping model predicts at unsampled locations. Blocks are assigned to
#' folds so that folds hold similar numbers of samples.
#'
#' @param .data A data frame or an `sf` POINT layer.
#' @param lon,lat <[`tidy-select`][tidyselect::language]> Coordinate columns
#'   (ignored for `sf` input).
#' @param v Number of folds.
#' @param block_size Block edge length in coordinate units (decimal degrees
#'   or metres). Defaults to a size giving about `4 * v` blocks across the
#'   larger extent.
#'
#' @return An `rset` of class `soil_spatial_block_cv` with columns `splits`
#'   and `id`.
#' @references Roberts, D. R., Bahn, V., Ciuti, S. et al. (2017).
#'   Cross-validation strategies for data with temporal, spatial,
#'   hierarchical, or phylogenetic structure. *Ecography*, 40, 913-929.
#'   \doi{10.1111/ecog.02881}
#' @family resampling
#' @export
#' @examples
#' set.seed(1)
#' folds <- soil_spatial_block_cv(png_profiles, lon = lon, lat = lat, v = 3)
#' folds
soil_spatial_block_cv <- function(.data, lon = NULL, lat = NULL, v = 5, block_size = NULL) {
  check_data_frame(.data)
  check_integerish(v, min = 2)
  check_number_positive(block_size, allow_null = TRUE)
  if (inherits(.data, "sf")) {
    rlang::check_installed("sf", reason = "to read coordinates from an sf object.")
    xy <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(.data)))
    x <- xy[, 1]
    y <- xy[, 2]
    data <- .data
  } else {
    x <- numeric_col(rlang::enquo(lon), .data, "lon")
    y <- numeric_col(rlang::enquo(lat), .data, "lat")
    data <- as_plain_tibble(.data)
  }
  if (anyNA(x) || anyNA(y)) {
    cli::cli_abort("Coordinates must not contain missing values.")
  }
  if (is.null(block_size)) {
    ext <- max(diff(range(x)), diff(range(y)))
    block_size <- if (ext > 0) ext / ceiling(sqrt(4 * v)) else 1
  }
  block <- paste(floor((x - min(x)) / block_size), floor((y - min(y)) / block_size), sep = "_")
  sizes <- table(block)
  if (length(sizes) < v) {
    cli::cli_abort(c(
      "Only {length(sizes)} block{?s} for {v} folds.",
      "i" = "Use a smaller {.arg block_size} or fewer folds."
    ))
  }
  fold_of <- balance_folds(sizes, v)
  make_fold_rset(unname(fold_of[block]), v, data,
    attrib = list(v = v, block_size = block_size), subclass = "soil_spatial_block_cv")
}

#' @export
print.soil_spatial_block_cv <- function(x, ...) {
  cat("#  Spatial block cross-validation: ", attr(x, "v"), " folds, block size ",
    signif(attr(x, "block_size"), 4), "\n", sep = "")
  class(x) <- c("tbl_df", "tbl", "data.frame")
  print(x, ...)
  invisible(x)
}

#' Kennard-Stone calibration and validation split
#'
#' Selects a calibration (analysis) set that spans the predictor space with
#' the Kennard and Stone (1969) algorithm, widely used to build soil
#' spectral libraries: the two most distant samples are taken first, then
#' repeatedly the sample farthest from those already selected. Distances
#' are Euclidean on the principal component scores of the standardised
#' predictors, which equals the Mahalanobis distance when all components are
#' kept.
#'
#' @param .data A data frame or [soil_spectra()] table.
#' @param ... <[`tidy-select`][tidyselect::language]> Predictor columns.
#'   Defaults to the spectral columns of a `<soil_spectra>`.
#' @param prop Proportion of samples in the calibration set.
#' @param n_comp Number of principal components used for distances
#'   (`NULL` for enough to explain 99 percent of the variance).
#'
#' @return An `rsplit` (class `soil_kennard_stone_split`); use
#'   [rsample::training()] and [rsample::testing()] on it.
#' @references Kennard, R. W. and Stone, L. A. (1969). Computer aided design
#'   of experiments. *Technometrics*, 11, 137-148.
#'   \doi{10.1080/00401706.1969.10490666}
#' @family resampling
#' @export
#' @examples
#' spec <- soil_spectra(soil_nir, dplyr::starts_with("nm_"))
#' split <- soil_kennard_stone(spec, prop = 0.7)
#' split
#' nrow(rsample::training(split))
soil_kennard_stone <- function(.data, ..., prop = 0.75, n_comp = NULL) {
  check_data_frame(.data)
  check_proportion(prop)
  check_integerish(n_comp %||% 1, min = 1)
  vars <- select_dots(.data, ...)
  if (!length(vars) && is_soil_spectra(.data)) vars <- soil_meta(.data)$spectra
  if (length(vars) < 1) {
    cli::cli_abort("Select predictor columns in {.arg ...}.")
  }
  x <- as.matrix(as.data.frame(.data)[vars])
  if (!is.numeric(x) || anyNA(x)) {
    cli::cli_abort("Predictors must be numeric without missing values.")
  }
  n <- nrow(x)
  k <- round(prop * n)
  if (k < 2 || k >= n) {
    cli::cli_abort("{.arg prop} gives {k} calibration samples out of {n}; adjust it.")
  }
  sds <- apply(x, 2, stats::sd)
  x <- x[, sds > 0, drop = FALSE]
  pc <- stats::prcomp(x, center = TRUE, scale. = TRUE)
  if (is.null(n_comp)) {
    n_comp <- which(cumsum(pc$sdev^2) / sum(pc$sdev^2) >= 0.99)[1]
  }
  n_comp <- min(n_comp, ncol(pc$x))
  s <- pc$x[, seq_len(n_comp), drop = FALSE]
  s <- sweep(s, 2, pc$sdev[seq_len(n_comp)], "/")
  d <- as.matrix(stats::dist(s))
  first <- which(d == max(d), arr.ind = TRUE)[1, ]
  sel <- unname(first)
  mind <- pmin(d[, sel[1]], d[, sel[2]])
  while (length(sel) < k) {
    mind[sel] <- -Inf
    nxt <- which.max(mind)
    sel <- c(sel, nxt)
    mind <- pmin(mind, d[, nxt])
  }
  data <- as_plain_tibble(.data)
  out <- rsample::make_splits(list(analysis = sort(sel), assessment = setdiff(seq_len(n), sel)),
    data = data)
  class(out) <- c("soil_kennard_stone_split", class(out))
  out
}
