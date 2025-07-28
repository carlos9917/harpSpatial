#' SLX Score Calculation
#'
#' This function calculates the SLX score for a given forecast and observation field.
#'
#' @param obfield A 2D matrix representing the observation field.
#' @param fcfield A 2D matrix representing the forecast field.
#' @param scales A vector of scales (in grid points) to compute the score for.
#' @param ... Other arguments, not used.
#' @return A tibble with columns for scale, S, L, and X.
#' @export
slx <- function(obfield, fcfield, scales, ...) {

  if (!is.matrix(fcfield) || !is.matrix(obfield)) {
    stop("fcfield and obfield must be matrices.")
  }
  if (!all(dim(fcfield) == dim(obfield))) {
    stop("fcfield and obfield must have the same dimensions.")
  }

  fcfield[is.na(fcfield)] <- 0
  obfield[is.na(obfield)] <- 0

  S_scores <- numeric(length(scales))
  L_scores <- numeric(length(scales))
  X_scores <- numeric(length(scales))

  for (i in seq_along(scales)) {
    l <- scales[i]
    if (l > 1) {
      # Assuming fastMean is available from the package
      fc_avg <- harpSpatial::fastMean(fcfield, l)
      obs_avg <- harpSpatial::fastMean(obfield, l)
    } else {
      fc_avg <- fcfield
      obs_avg <- obfield
    }

    # S component (Amplitude)
    s_l <- sum(fc_avg)
    o_l <- sum(obs_avg)
    if ((s_l + o_l) == 0) {
      S_scores[i] <- 0
    } else {
      S_scores[i] <- (s_l - o_l) / (s_l + o_l)
    }

    # L component (Location)
    gX <- which(fc_avg > 0, arr.ind = TRUE)
    gY <- which(obs_avg > 0, arr.ind = TRUE)

    if (nrow(gX) == 0 || nrow(gY) == 0) {
        L_scores[i] <- NA
    } else {
        cX <- colMeans(gX)
        cY <- colMeans(gY)
        d <- sqrt(sum((cX - cY)^2))
        D <- sqrt(sum(dim(fcfield)^2))
        L_scores[i] <- d / D
    }

    # X component (Structure)
    denom <- sum(fc_avg^2) + sum(obs_avg^2)
    if (denom == 0) {
        X_scores[i] <- 1 # Perfect score if both fields are zero
    } else {
        X_scores[i] <- 2 * sum(fc_avg * obs_avg) / denom
    }
  }

  tibble::tibble(
    scale = scales,
    S     = S_scores,
    L     = L_scores,
    X     = X_scores
  )
}