# Agreement Scales Implementation
# Based on Dey et al. (2016) - Equations 1, 2, 3

# Similarity measure D (Equation 1 from Dey et al. 2016)
similarity_D <- function(a, b) {
  # Handle the case where both values are zero
  zero_both <- (a == 0 & b == 0)
  #print("in similarity")
  # Calculate D for non-zero cases
  # Extract numeric matrices from geofield objects
  a_mat <- as.matrix(a)
  b_mat <- as.matrix(b)

  # Now do the arithmetic
  # Compute numerator and denominator
  numerator <- (a_mat - b_mat)^2
  denominator <- (a_mat^2 + b_mat^2)

  # Compute the expression element-wise
  D <- numerator / denominator
  #D <- (a - b)^2 / pmax(a^2 + b^2, 1e-12)  # Small epsilon to avoid division by zero


  # Set D = 1 where both values are zero (as per Dey et al. 2016)
  D[zero_both] <- 1
  return(D)
}


#--- Efficient window mean using cumsum_2d (integral image) -------------------
# This uses the summed-area table approach for O(1) window mean calculation
# after O(N) preprocessing, making the total complexity O(N*S_lim) instead of O(N*S_lim^3)

window_mean_cumsum2d <- function(mat, k) {
  if (k == 0L) return(mat)

  # Get dimensions
  ny <- nrow(mat)
  nx <- ncol(mat)

  # Handle NA values by setting them to 0 for cumsum calculation
  mat_clean <- mat
  na_mask <- is.na(mat)
  mat_clean[na_mask] <- 0

  # Calculate cumulative sum using harpSpatial::cumsum_2d
  cumsum_mat <- cumsum_2d(mat_clean)

  # Create output matrix
  result <- matrix(NA, ny, nx)

  # Calculate window means using integral image
  for (i in 1:ny) {
    for (j in 1:nx) {
      # Define window bounds
      i_min <- max(1, i - k)
      i_max <- min(ny, i + k)
      j_min <- max(1, j - k)
      j_max <- min(nx, j + k)

      # Calculate sum using integral image
      # Sum = cumsum[i_max, j_max] - cumsum[i_min-1, j_max] - cumsum[i_max, j_min-1] + cumsum[i_min-1, j_min-1]
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

      # Calculate number of points in window
      n_points <- (i_max - i_min + 1) * (j_max - j_min + 1)

      # Calculate mean
      result[i, j] <- sum_val / n_points
    }
  }

  # Restore NA values where original data had NAs in the center point
  result[na_mask] <- NA

  return(result)
}


# Agreement scale calculation
agreement_scale_map <- function(f1, f2, alpha = 0.5, S_lim = 80L) {
  # Ensure we work with regular matrices, not geofield objects
  if (inherits(f1, "geofield")) {
    f1 <- as.array(f1)
  }
  if (inherits(f2, "geofield")) {
    f2 <- as.array(f2)
  }

  ny <- nrow(f1)
  nx <- ncol(f1)

  # Initialize agreement scale matrix with maximum scale
  SA <- matrix(S_lim, ny, nx)

  cat(sprintf("Calculating agreement scales for %d x %d grid...\n", ny, nx))

  # Loop over scales from 0 to S_lim
  for (S in 0:S_lim) {
    if (S %% 10 == 0) {
      cat(sprintf("  Processing scale S = %d/%d\n", S, S_lim))
    }

    # Calculate neighborhood means
    f1_bar <- window_mean_cumsum2d(f1, S)
    f2_bar <- window_mean_cumsum2d(f2, S)

    # Calculate similarity measure D (Equation 1 from Dey et al. 2016)
    D <- similarity_D(f1_bar, f2_bar)

    # Calculate agreement criterion threshold (Equation 3 from Dey et al. 2016)
    D_crit <- alpha + (1 - alpha) * S / S_lim

    # Find points where agreement is achieved for the first time
    agreement_achieved <- (D <= D_crit) & !is.na(D)
    first_agreement <- agreement_achieved & (SA == S_lim)

    # Update agreement scale for points achieving agreement for first time
    SA[first_agreement] <- S

    # Early termination if all valid points have found their agreement scale
    remaining_points <- sum(SA == S_lim & !is.na(f1) & !is.na(f2))
    if (remaining_points == 0) {
      cat(sprintf("  All points converged at scale S = %d\n", S))
      break
    }
  }

  cat("Agreement scale calculation completed!\n")
  return(SA)
}


# Ensemble agreement scales (SA_mm)
calculate_SA_mm <- function(ensemble_fields, alpha = 0.5, S_lim = 80L) {
  n_members <- length(ensemble_fields)

  if (n_members < 2) {
    stop("Need at least 2 ensemble members to calculate SA(mm)")
  }

  # Calculate number of member pairs
  n_pairs <- n_members * (n_members - 1) / 2
  cat(sprintf("Calculating SA(mm) for %d members (%d pairs)...\n", n_members, n_pairs))

  # Store all pairwise agreement scales
  agreement_scales <- list()
  pair_count <- 0

  # Calculate agreement scales for all member pairs
  for (i in 1:(n_members - 1)) {
    for (j in (i + 1):n_members) {
      pair_count <- pair_count + 1
      cat(sprintf("Processing member pair %d-%d (%d/%d)\n", i, j, pair_count, n_pairs))

      scale_map <- agreement_scale_map(
        ensemble_fields[[i]],
        ensemble_fields[[j]],
        alpha = alpha,
        S_lim = S_lim
      )

      agreement_scales[[pair_count]] <- scale_map
    }
  }

  # Average over all pairs (Equation 6 from Dey et al. 2016)
  SA_mm <- Reduce("+", agreement_scales) / length(agreement_scales)

  cat(sprintf("SA(mm) calculation completed. Mean value: %.2f grid points\n",
              mean(SA_mm, na.rm = TRUE)))

  return(SA_mm)
}

# Member-observation agreement scales (SA_mo)
calculate_SA_mo <- function(ensemble_fields, observations, alpha = 0.5, S_lim = 80L) {
  n_members <- length(ensemble_fields)

  cat(sprintf("Calculating SA(mo) for %d members vs observations...\n", n_members))

  # Store agreement scales for all member-observation pairs
  agreement_scales <- list()

  # Calculate agreement scales for each member vs observations
  for (i in 1:n_members) {
    cat(sprintf("Processing member %d/%d vs observations\n", i, n_members))

    scale_map <- agreement_scale_map(
      ensemble_fields[[i]],
      observations,
      alpha = alpha,
      S_lim = S_lim
    )

    agreement_scales[[i]] <- scale_map
  }

  # Average over all members (Equation 7 from Dey et al. 2016)
  SA_mo <- Reduce("+", agreement_scales) / length(agreement_scales)

  cat(sprintf("SA(mo) calculation completed. Mean value: %.2f grid points\n",
              mean(SA_mo, na.rm = TRUE)))

  return(SA_mo)
}

ens_agreement_scales <- function(obfield, fcfield, scales, thresholds, ...) {

  if (!is.list(fcfield)) {
    stop("fcfield must be a list of ensemble members for ens_agreement_scales")
  }

  alpha <- thresholds[1]
  S_lim <- scales[1]

  SA_mm <- calculate_SA_mm(fcfield, alpha = alpha, S_lim = S_lim)
  SA_mo <- calculate_SA_mo(fcfield, obfield, alpha = alpha, S_lim = S_lim)

  result <- tibble::tibble(
    SA_mm = mean(SA_mm, na.rm = TRUE),
    SA_mo = mean(SA_mo, na.rm = TRUE)
  )
  result
}

fo_agreement_scales <- function(obfield, fcfield, scales, thresholds, ...) {

  if (is.list(fcfield)) {
    stop("fcfield must be a single forecast field for fo_agreement_scales")
  }

  alpha <- thresholds[1]
  S_lim <- scales[1]

  SA_fo <- agreement_scale_map(fcfield, obfield, alpha = alpha, S_lim = S_lim)

  result <- tibble::tibble(
    SA_fo = mean(SA_fo, na.rm = TRUE)
  )
  result
}