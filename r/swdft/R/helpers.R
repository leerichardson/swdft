#' Get the maximum DFT coefficient
#'
#' @param x numeric vector
#'
#' @export
#'
#' @return numeric of largest frequency. Will be between 0 and .5
#'
get_max_freq <- function(x) {
  periodogram <- Mod(stats::fft(x))^2
  freqs <- (0:(length(x) - 1)) / length(x)
  max_freq <- freqs[which.max(periodogram)]
  return(max_freq)
}

#' Simple high pass filter
#'
#' @param x the vector or time-series
#' @param order the order of the filter
#'
moving_average <- function(x, order) {
  stats::filter(x, rep(1 / order, order), sides=2)
}

#' Phase unwrapping
#'
#' @param p vector of phases fit by demodulation
#'
unwrap_phase <- function(p) {
  pdiff <- diff(p, na.rm=TRUE)
  pd_no_nas <- pdiff
  pd_no_nas[is.na(pdiff)] <- 0

  p_no_nas <- p
  p_no_nas[is.na(p)] <- 0

  p[] <- cumsum(c(p_no_nas[1], pd_no_nas - round(pd_no_nas)))
  return(p)
}

#' Convert the SWDFT to proportions of frequency
#'
#' @param a swdft
#'
swdft_to_props <- function(a) {
  n <- nrow(a)
  m <- floor(n / 2)
  if (class(a[1, 1]) == "complex") { amod <- Mod(a)^2 }

  return( apply(X=amod[2:m, ], MARGIN=2, FUN=function(x) x / sum(x)) )
}
