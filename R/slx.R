# Original functions from slx_score.R

S_score <- function(ob, phi, k = 0.1, A = 4) {
  if (ob > k) {
    if (phi < ob - k) return(phi / (ob - k))
    if (phi <= ob) return(1.0)
    return(max(1 - (phi - ob) / (A * ob), 0.0))
  }
  if (phi <= k) return(1.0)
  max(1 - (phi - k) / (A * k), 0.0)
}

local_extreme_indices <- function(field, mode = "max", tolerance = 0.0) {
  nrow_f <- nrow(field)
  ncol_f <- ncol(field)
  extrema <- matrix(numeric(0), 0, 3, dimnames = list(NULL, c("row", "col", "value")))

  for (i in 2:(nrow_f - 1)) {
    for (j in 2:(ncol_f - 1)) {
      val <- field[i, j]
      if (is.na(val)) next
      neighborhood <- field[(i - 1):(i + 1), (j - 1):(j + 1)]
      neighbors <- as.vector(neighborhood)[-5]

      if (mode == "max") {
        if (all(val >= neighbors - tolerance) && any(val > neighbors + tolerance)) {
          extrema <- rbind(extrema, c(i, j, val))
        }
      } else if (mode == "min") {
        if (all(val <= neighbors + tolerance) && any(val < neighbors - tolerance)) {
          extrema <- rbind(extrema, c(i, j, val))
        }
      }
    }
  }
  colnames(extrema) <- c("row", "col", "value")
  return(extrema)
}

neighbourhood_view <- function(field, i, j, L, mode = 'max') {
  i0 <- max(1, i - L)
  i1 <- min(nrow(field), i + L)
  j0 <- max(1, j - L)
  j1 <- min(ncol(field), j + L)
  neigh <- field[i0:i1, j0:j1]
  if (mode == 'max') return(max(neigh, na.rm = TRUE))
  min(neigh, na.rm = TRUE)
}

SLX_components <- function(analysis, forecast, L, delta = 0.0) {
  ob_max_pts <- local_extreme_indices(analysis, 'max', delta)
  ob_min_pts <- local_extreme_indices(analysis, 'min', delta)
  fc_max_pts <- local_extreme_indices(forecast, 'max', delta)
  fc_min_pts <- local_extreme_indices(forecast, 'min', delta)
  avg_score <- function(pts, ob_field, fc_field, mode) {
    if (nrow(pts) == 0) return(NA)
    scores <- sapply(1:nrow(pts), function(idx) {
      i <- pts[idx, 1]
      j <- pts[idx, 2]
      if (mode == 'ob_max') {
        ob <- ob_field[i, j]
        phi <- neighbourhood_view(fc_field, i, j, L, 'max')
      } else if (mode == 'ob_min') {
        ob <- ob_field[i, j]
        phi <- neighbourhood_view(fc_field, i, j, L, 'min')
      } else if (mode == 'fc_max') {
        ob <- neighbourhood_view(ob_field, i, j, L, 'max')
        phi <- fc_field[i, j]
      } else { # fc_min
        ob <- neighbourhood_view(ob_field, i, j, L, 'min')
        phi <- fc_field[i, j]
      }
      S_score(ob, phi)
    })
    mean(scores, na.rm = TRUE)
  }
  s_ob_max <- avg_score(ob_max_pts, analysis, forecast, 'ob_max')
  s_ob_min <- avg_score(ob_min_pts, analysis, forecast, 'ob_min')
  s_fc_max <- avg_score(fc_max_pts, analysis, forecast, 'fc_max')
  s_fc_min <- avg_score(fc_min_pts, analysis, forecast, 'fc_min')

  list(
    s_ob_max = s_ob_max,
    s_ob_min = s_ob_min,
    s_fc_max = s_fc_max,
    s_fc_min = s_fc_min,
    combined = mean(c(s_ob_max, s_ob_min, s_fc_max, s_fc_min), na.rm = TRUE)
  )
}


#' SLX Score Calculation wrapper for verify_spatial
#'
#' This function calculates the SLX score for a given forecast and observation field.
#' It acts as a wrapper around the original SLX_components function to be compatible
#' with the verify_spatial workflow.
#'
#' @param obfield A geofield object representing the observation field.
#' @param fcfield A geofield object representing the forecast field.
#' @param scales A vector of scales (in grid points) to compute the score for.
#' @param ... Other arguments, not used.
#' @return A tibble with columns for scale, S, L, and X.
#' @export
slx <- function(obfield, fcfield, scales, ...) {

  # Transformation from geofield to matrix, as in load_field_data_slx.R
  obs_dims <- dim(obfield)
  obs_values <- as.numeric(obfield)
  obs_matrix <- matrix(obs_values, nrow = obs_dims[1], ncol = obs_dims[2])

  fc_dims <- dim(fcfield)
  fc_values <- as.numeric(fcfield)
  fc_matrix <- matrix(fc_values, nrow = fc_dims[1], ncol = fc_dims[2])
  #Handle NA values
  obs_matrix[is.na(obs_matrix)] <- 0
  fc_matrix[is.na(fc_matrix)] <- 0

  # Calculate scores for each scale
  results <- lapply(scales, function(l) {
    res <- SLX_components(obs_matrix, fc_matrix, L = l)
    tibble::tibble(
      scale = l,
      S_OB_MAX = res$s_ob_max,
      S_OB_MIN = res$s_ob_min,
      S_FC_MAX = res$s_fc_max,
      S_FC_MIN = res$s_fc_min,
      SLX_COMB = res$combined
    )
  })

  # Combine results into a single tibble
  dplyr::bind_rows(results)
}
