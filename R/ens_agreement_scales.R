# Functions to calculate ensemble agreement scales based on Dey et al. (2016).
# This file provides the function ens_agreement, which is intended to be
# called from other verification scripts. The helper functions are adapted
# from the `enhanced_dey_functions.R` script and are used internally.

#' Calculate the similarity score D.
#'
#' This function implements Equation 1 from Dey et al. (2016).
#'
#' @param f1_bar A numeric matrix or vector for the first field.
#' @param f2_bar A numeric matrix or vector for the second field.
#' @return A matrix or vector of similarity scores.
.similarity_D <- function(f1_bar, f2_bar) {
  if (is.matrix(f1_bar) || is.array(f1_bar)) {
    f1 <- as.numeric(f1_bar)
    f2 <- as.numeric(f2_bar)
  } else {
    f1 <- f1_bar
    f2 <- f2_bar
  }

  D <- numeric(length(f1))
  both_zero <- (f1 == 0) & (f2 == 0)
  D[both_zero] <- 1

  not_both_zero <- !both_zero
  if (any(not_both_zero)) {
    numerator <- (f1[not_both_zero] - f2[not_both_zero])^2
    denominator <- f1[not_both_zero]^2 + f2[not_both_zero]^2
    # Handle cases where denominator is zero (e.g. one is zero, one is not)
    D[not_both_zero] <- ifelse(denominator == 0, 1, numerator / denominator)
  }

  if (is.matrix(f1_bar) || is.array(f1_bar)) {
    D <- array(D, dim = dim(f1_bar))
  }

  return(D)
}

#' Optimized window mean using cumulative sum (integral image).
#'
#' @param mat A numeric matrix.
#' @param k The half-width of the window.
#' @return A matrix of window means.
.window_mean <- function(mat, k) {
  if (k == 0) return(mat)

  ny <- nrow(mat)
  nx <- ncol(mat)

  mat_clean <- mat
  na_mask <- is.na(mat)
  mat_clean[na_mask] <- 0

  padded <- matrix(0, ny + 2*k, nx + 2*k)
  padded[(k+1):(ny+k), (k+1):(nx+k)] <- mat_clean

  cumsum_mat <- matrix(0, ny + 2*k + 1, nx + 2*k + 1)
  for (i in 2:(ny + 2*k + 1)) {
    for (j in 2:(nx + 2*k + 1)) {
      cumsum_mat[i, j] <- padded[i-1, j-1] +
                          cumsum_mat[i-1, j] +
                          cumsum_mat[i, j-1] -
                          cumsum_mat[i-1, j-1]
    }
  }

  result <- matrix(NA, ny, nx)
  window_area <- (2*k + 1)^2

  for (i in 1:ny) {
    for (j in 1:nx) {
      i1 <- i
      j1 <- j
      i2 <- i + 2*k + 1
      j2 <- j + 2*k + 1
      window_sum <- cumsum_mat[i2, j2] -
                    cumsum_mat[i1, j2] -
                    cumsum_mat[i2, j1] +
                    cumsum_mat[i1, j1]
      result[i, j] <- window_sum / window_area
    }
  }

  result[na_mask] <- NA
  return(result)
}

#' Calculate the agreement scale map.
#'
#' @param f1 A numeric matrix for the first field.
#' @param f2 A numeric matrix for the second field.
#' @param alpha The agreement criterion parameter.
#' @param S_lim The maximum neighbourhood half-width.
#' @param verbose If TRUE, print progress messages.
#' @return A matrix containing the agreement scale at each grid point.
.agreement_scale_map <- function(f1, f2, alpha = 0.5, S_lim = 80L, verbose = TRUE) {
  if (inherits(f1, "geofield")) f1 <- as.array(f1)
  if (inherits(f2, "geofield")) f2 <- as.array(f2)

  ny <- nrow(f1)
  nx <- ncol(f1)
  SA <- matrix(S_lim, ny, nx)

  if (verbose) cat(sprintf("Calculating agreement scales for %d x %d grid...\n", ny, nx))

  for (S in 0:S_lim) {
    if (verbose && S %% 10 == 0) {
      cat(sprintf("  Processing scale S = %d/%d\n", S, S_lim))
    }

    f1_bar <- .window_mean(f1, S)
    f2_bar <- .window_mean(f2, S)
    D <- .similarity_D(f1_bar, f2_bar)
    D_crit <- alpha + (1 - alpha) * S / S_lim

    agreement_achieved <- (D <= D_crit) & !is.na(D)
    first_agreement <- agreement_achieved & (SA == S_lim)
    SA[first_agreement] <- S
  }

  if (verbose) cat("Agreement scale calculation completed!\n")
  return(SA)
}

#' Calculate ensemble member-member agreement scales (SA_mm).
#'
#' @param ensemble_fields A list of numeric matrices for the ensemble members.
#' @param alpha The agreement criterion parameter.
#' @param S_lim The maximum neighbourhood half-width.
#' @param verbose If TRUE, print progress messages.
#' @return A matrix containing the SA_mm values.
.calculate_SA_mm <- function(ensemble_fields, alpha = 0.5, S_lim = 80L, verbose = TRUE) {
  n_members <- length(ensemble_fields)
  if (n_members < 2) stop("Need at least 2 ensemble members for SA(mm)")

  n_pairs <- n_members * (n_members - 1) / 2
  if (verbose) cat(sprintf("Calculating SA(mm) for %d members (%d pairs)...\n", n_members, n_pairs))

  agreement_scales <- list()
  pair_count <- 0
  for (i in 1:(n_members - 1)) {
    for (j in (i + 1):n_members) {
      pair_count <- pair_count + 1
      if (verbose) cat(sprintf("Processing member pair %d-%d (%d/%d)\n", i, j, pair_count, n_pairs))
      scale_map <- .agreement_scale_map(
        ensemble_fields[[i]],
        ensemble_fields[[j]],
        alpha = alpha,
        S_lim = S_lim,
        verbose = FALSE # Verbosity handled in this function
      )
      agreement_scales[[pair_count]] <- scale_map
    }
  }

  SA_mm <- Reduce("+", agreement_scales) / length(agreement_scales)
  if (verbose) cat(sprintf("SA(mm) calculation completed. Mean value: %.2f grid points\n", mean(SA_mm, na.rm = TRUE)))
  return(SA_mm)
}

#' Calculate ensemble member-observation agreement scales (SA_mo).
#'
#' @param ensemble_fields A list of numeric matrices for the ensemble members.
#' @param observations A numeric matrix for the observation field.
#' @param alpha The agreement criterion parameter.
#' @param S_lim The maximum neighbourhood half-width.
#' @param verbose If TRUE, print progress messages.
#' @return A matrix containing the SA_mo values.
.calculate_SA_mo <- function(ensemble_fields, observations, alpha = 0.5, S_lim = 80L, verbose = TRUE) {
  n_members <- length(ensemble_fields)
  if (verbose) cat(sprintf("Calculating SA(mo) for %d members vs observations...\n", n_members))

  agreement_scales <- list()
  for (i in 1:n_members) {
    if (verbose) cat(sprintf("Processing member %d/%d vs observations\n", i, n_members))
    scale_map <- .agreement_scale_map(
      ensemble_fields[[i]],
      observations,
      alpha = alpha,
      S_lim = S_lim,
      verbose = FALSE # Verbosity handled in this function
    )
    agreement_scales[[i]] <- scale_map
  }

  SA_mo <- Reduce("+", agreement_scales) / length(agreement_scales)
  if (verbose) cat(sprintf("SA(mo) calculation completed. Mean value: %.2f grid points\n", mean(SA_mo, na.rm = TRUE)))
  return(SA_mo)
}

#' Calculate ensemble agreement scales (SA_mm and SA_mo).
#'
#' This function serves as the main entry point for calculating the
#' ensemble agreement scales (SA_mm and SA_mo) based on the method by
#' Dey et al. (2016).
#'
#' @param ensemble_fields A list of numeric matrices representing the ensemble forecast fields.
#' @param observation_field A numeric matrix representing the observation field,
#'   with dimensions matching the ensemble fields.
#' @param alpha A parameter for the agreement criterion. Default is 0.5.
#' @param S_lim An integer for the maximum neighbourhood half-width. Default is 80.
#' @param verbose A boolean to control progress messages. Default is TRUE.
#'
#' @return A list with four elements:
#'   \item{SA_mm}{A matrix of member-member agreement scales.}
#'   \item{SA_mo}{A matrix of member-observation agreement scales.}
#'   \item{summary_stats_mm}{A list of summary statistics for SA_mm.}
#'   \item{summary_stats_mo}{A list of summary statistics for SA_mo.}
#' @export
ens_agreement_scales <- function(ensemble_fields, observation_field, alpha = 0.5, S_lim = 80L, verbose = TRUE) {

  # --- Input validation ---
  if (!is.list(ensemble_fields) || length(ensemble_fields) < 2) {
    stop("ensemble_fields must be a list of at least 2 matrices.")
  }
  if (!all(sapply(ensemble_fields, is.matrix))) {
    stop("All elements of ensemble_fields must be matrices.")
  }
  dims <- dim(ensemble_fields[[1]])
  if (!all(sapply(ensemble_fields, function(m) all(dim(m) == dims)))) {
    stop("All matrices in ensemble_fields must have the same dimensions.")
  }
  if (!is.matrix(observation_field)) {
    stop("observation_field must be a matrix.")
  }
  if (!all(dim(observation_field) == dims)) {
    stop("observation_field must have the same dimensions as the ensemble fields.")
  }

  # --- Calculations ---
  SA_mm <- .calculate_SA_mm(ensemble_fields, alpha = alpha, S_lim = S_lim, verbose = verbose)
  SA_mo <- .calculate_SA_mo(ensemble_fields, observation_field, alpha = alpha, S_lim = S_lim, verbose = verbose)

  # --- Summary statistics ---
  summary_stats_mm <- list(
    mean_agreescale = mean(SA_mm, na.rm = TRUE),
    min_agreescale  = min(SA_mm, na.rm = TRUE),
    max_agreescale  = max(SA_mm, na.rm = TRUE),
    sd_agreescale   = sd(SA_mm, na.rm = TRUE)
  )

  summary_stats_mo <- list(
    mean_agreescale = mean(SA_mo, na.rm = TRUE),
    min_agreescale  = min(SA_mo, na.rm = TRUE),
    max_agreescale  = max(SA_mo, na.rm = TRUE),
    sd_agreescale   = sd(SA_mo, na.rm = TRUE)
  )

  # --- Return results ---
  list(
    SA_mm              = SA_mm,
    SA_mo              = SA_mo,
    summary_stats_mm   = summary_stats_mm,
    summary_stats_mo   = summary_stats_mo
  )
}
