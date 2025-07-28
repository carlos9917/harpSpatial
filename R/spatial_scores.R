# various spatial verification scores
#' Calculate spatial scores
#' @param score The score to calculate. If NULL, the function will return a list of available scores.
#' @param obfield A matrix containing the observation field. If NULL, the function will return
#'   the table structure for the specified score.
#' @param fcfield A matrix containing the forecast field. Must have the same dimension as obfield.
#' @param ... Other options that may depend on the score (like scale, threshold, ...)
#' @export
spatial_scores <- function(score = NULL, obfield = NULL, fcfield = NULL, ...) {
  # TODO: add score options, plot_func and plot_opt
  # FIXME: you MUST indicate the primary fields (e.g. threshold & scale) !
  score_list <- list(
                     "bias"   = list(fields = c("bias"), "func" = "scores_sp_basic", "plot_func" = "plot_basic"),
                     "mse"    = list(fields = c("mse"),  "func" = "scores_sp_basic", "plot_func" = "plot_basic"),
                     "mae"    = list(fields = c("mae"),  "func" = "scores_sp_basic", "plot_func" = "plot_basic"),
#                     "gridded" = list(fields = c("bias", "mse"), "func" = "score_sp_gridded"),
                     "SAL"     = list(fields = c("S", "A", "L"), "func" = "SAL", "plot_func" = "plot_sal"),
                     "FSS"     = list(fields = c("fss"), primary = c("threshold", "scale"),
                                      "func" = "scores_sp_neighborhood", "plot_func" = "plot_fss"),
                     "NACT"    = list(fields = c("hit", "fa", "miss", "cr"), primary = c("threshold", "scale"),
                                      "func" = "scores_sp_neighborhood", "plot_func" = "plot_nact"),
                     "SLX"     = list(fields = c("S", "L", "X"), primary = c("scale"), "func" = "slx")
#                     , "FSS_p"     = list(fields = c("percentile", "scale", "fss"), "func" = "score_fss", "plot_func" = "plot_fss")
                     )

  # if called without "score", return a list of all scores
  if (is.null(score)) return(score_list)
  else if (!is.element(score, names(score_list))) stop("Unknown score ", score, ".\n")

  # Derive table structure
  # table_structure <- spatial_score_table(score_list[[score]]$fields)
  # if called without "obfield" and "fcfield", just return the table structure for the given score
  if (is.null(obfield) && is.null(fcfield)) {
    return(score_list[[score]])
  }

  # FIXME: we may be calling with options that are not recognised/used by the score
#  arglist <- names(as.list(args(score_list[[score]]$func)))
#  message("score function: ", score_list[[score]]$func)
#  message("argument list: ", paste(arglist, collapse=" "))

  myargs <- c(list(obfield = obfield, fcfield = fcfield), list(...))
  do.call(score_list[[score]]$func, myargs)
}

# simple wrappers to deal with unwanted arguments
# Yes, there are probably nicer ways...

##' @export
scores_sp_basic <- function(obfield, fcfield, ...) {
  harpSpatial_basic_scores(obfield=obfield, fcfield=fcfield)
}

##' @export
scores_sp_neighborhood <- function(obfield, fcfield, thresholds, scales, ...) {
  message("obfield dimensions: ", paste(dim(obfield), collapse="x"))
  message("fcfield dimensions: ", paste(dim(fcfield), collapse="x"))
  message("thresholds", paste(thresholds, collapse=","))
  message("scales", paste(scales, scales=","))
  harpSpatial_neighborhood_scores(obfield=obfield, fcfield=fcfield,
                                  thresholds=thresholds, scales=scales
  )
}
