# Functions to calculate agreement scales based on Dey et al. (2016)
# This file provides the function fo_agreement, which is intended to be
# called from other verification scripts like verify_spatial.R.

#' Efficiently calculate window mean using a summed-area table.
#'
#' This function uses the summed-area table approach (integral image) for
#' O(1) window mean calculation after an O(N) preprocessing step. It relies
#' on the `cumsum_2d` function from the `harpSpatial` package.
#'
#' @param mat A numeric matrix.
#' @param k The half-width of the square window (in grid points). The full
#'   window size will be (2k+1)x(2k+1).
#' @return A matrix of the same dimensions as `mat` with the window mean
#'   calculated at each point.
window_mean_cumsum2d <- function(mat, k) {
  if (k == 0L) return(mat)

  ny <- nrow(mat)
  nx <- ncol(mat)

  mat_clean <- mat
  na_mask <- is.na(mat)
  mat_clean[na_mask] <- 0

  # Calculate cumulative sum using harpSpatial::cumsum_2d
  #cumsum_mat <- harpSpatial::cumsum_2d(mat_clean)
  cumsum_mat <- cumsum_2d(mat_clean)

  result <- matrix(NA, ny, nx)

  for (i in 1:ny) {
    for (j in 1:nx) {
      i_min <- max(1, i - k)
      i_max <- min(ny, i + k)
      j_min <- max(1, j - k)
      j_max <- min(nx, j + k)

      sum_val <- cumsum_mat[i_max, j_max]

      if (i_min > 1) {
        sum_val <- sum_val - cumsum_mat[i_min - 1, j_max]
      }
      if (j_min > 1) {
        sum_val <- sum_val - cumsum_mat[i_max, j_min - 1]
      }
      if (i_min > 1 && j_min > 1) {
        sum_val <- sum_val + cumsum_mat[i_min - 1, j_min - 1]
      }

      n_points <- (i_max - i_min + 1) * (j_max - j_min + 1)
      result[i, j] <- sum_val / n_points
    }
  }

  result[na_mask] <- NA
  return(result)
}

#' Calculate the similarity score D.
#'
#' This function implements Equation 1 from Dey et al. (2016).
#'
#' @param a A numeric matrix.
#' @param b A numeric matrix of the same dimensions as `a`.
#' @return A matrix of similarity scores.
similarity_D <- function(a, b) {
  zero_both <- (a == 0 & b == 0)

  a_mat <- as.matrix(a)
  b_mat <- as.matrix(b)

  numerator <- (a_mat - b_mat)^2
  denominator <- (a_mat^2 + b_mat^2)

  D <- numerator / denominator
  D[denominator == 0 & numerator == 0] <- 0 # Should be covered by zero_both
  D[zero_both] <- 1

  return(D)
}

#' Calculate the agreement scale map.
#'
#' This function implements the core algorithm from Dey et al. (2016) to find
#' the smallest scale at which two fields agree at each grid point.
#'
#' @param f1 A numeric matrix for the first field.
#' @param f2 A numeric matrix for the second field.
#' @param alpha The agreement criterion parameter.
#' @param S_lim The maximum neighbourhood half-width.
#' @param verbose If TRUE, print progress messages.
#' @return A matrix containing the agreement scale at each grid point.
agreement_scale_map <- function(f1, f2, alpha = 0.5, S_lim = 80L, verbose = TRUE) {

  ny <- nrow(f1)
  nx <- ncol(f1)
  SA <- matrix(S_lim, ny, nx)
  valid_points_mask <- !is.na(f1) & !is.na(f2)

  if (verbose) {
    cat(sprintf("Calculating agreement scales for %d x %d grid...\n", ny, nx))
  }

  for (S in 0:S_lim) {
    if (verbose && S %% 10 == 0) {
      cat(sprintf("  Processing scale S = %d/%d\n", S, S_lim))
    }

    f1_bar <- window_mean_cumsum2d(f1, S)
    f2_bar <- window_mean_cumsum2d(f2, S)

    D <- similarity_D(f1_bar, f2_bar)
    D_crit <- alpha + (1 - alpha) * S / S_lim

    agreement_achieved <- (D <= D_crit) & !is.na(D)
    first_agreement <- agreement_achieved & (SA == S_lim)

    SA[first_agreement] <- S

    if (sum(SA[valid_points_mask] == S_lim) == 0) {
      if (verbose) {
        cat(sprintf("  All points converged at scale S = %d\n", S))
      }
      break
    }
  }

  if (verbose) {
    cat("Agreement scale calculation completed!\n")
  }
  return(SA)
}

#' Calculate agreement scales for a forecast and observation field.
#'
#' This function serves as the main entry point for calculating the
#' forecast-observation agreement scales (SA_fo) based on the method by
#' Dey et al. (2016).
#'
#' @param fc_field A numeric matrix representing the forecast field.
#' @param obs_field A numeric matrix representing the observation field,
#'   with dimensions matching `fc_field`.
#' @param alpha A parameter for the agreement criterion. Default is 0.5.
#' @param S_lim An integer for the maximum neighbourhood half-width. Default is 80.
#' @param verbose A boolean to control progress messages. Default is TRUE.
#'
#' @return A tibble with the agreement scale map (SA_fo) and its summary statistics.
#'   The tibble has one row, with the SA_fo matrix in a list-column.
#' @export
fo_agreement_scales <- function(fcfield, obfield, alpha = 0.5, S_lim = 80L, verbose=TRUE, return_full_matrix = TRUE, ...) {

  obs_dims <- dim(obfield)
  obs_values <- as.numeric(obfield)
  obs_matrix <- matrix(obs_values, nrow = obs_dims[1], ncol = obs_dims[2])

  fc_dims <- dim(fcfield)
  fc_values <- as.numeric(fcfield)
  fc_matrix <- matrix(fc_values, nrow = fc_dims[1], ncol = fc_dims[2])

  if (!is.matrix(fc_matrix) || !is.matrix(obs_matrix)) {
    stop("fc_field and obs_field must be matrices.")
  }
  if (!all(dim(fcfield) == dim(obfield))) {
    stop("fc_field and obs_field must have the same dimensions.")
  }

  SA_fo <- agreement_scale_map(
    fc_matrix,
    obs_matrix,
    alpha = alpha,
    S_lim = S_lim,
    verbose = verbose
  )

  summary_stats <- list(
    mean_agreescale = mean(SA_fo, na.rm = TRUE),
    min_agreescale  = min(SA_fo, na.rm = TRUE),
    max_agreescale  = max(SA_fo, na.rm = TRUE),
    sd_agreescale   = sd(SA_fo, na.rm = TRUE)
  )

  if (return_full_matrix) {
    tibble::tibble(
      SA_fo               = list(SA_fo),
      mean_agreescale     = summary_stats$mean_agreescale,
      min_agreescale      = summary_stats$min_agreescale,
      max_agreescale      = summary_stats$max_agreescale,
      sd_agreescale       = summary_stats$sd_agreescale
    )
  } else {
    tibble::tibble(
      SA_fo               = summary_stats$mean_agreescale,
      mean_agreescale     = summary_stats$mean_agreescale,
      min_agreescale      = summary_stats$min_agreescale,
      max_agreescale      = summary_stats$max_agreescale,
      sd_agreescale       = summary_stats$sd_agreescale
    )
  }
}
