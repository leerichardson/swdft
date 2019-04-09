#' Simple high pass filter
#'
#' @param x the vector or time-series
#' @param order
#'
moving_average <- function(x, order) {
  stats::filter(x, rep(1 / order, order), sides=2)
}
