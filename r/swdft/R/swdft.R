#' Sliding Window Discrete Fourier Transform (SWDFT)
#'
#' @param x real or complex vector
#' @param n integer window size.
#' @param type algorithm to implement. defaults to "fftw"
#' @param pad optionally 0 pad the array to that the output
#' array has the same dimension as the original time-series
#' @param taper type of taper for each window position. defaults
#' to 'none', can
#' @param p Proportion to be tapered at each end of the series. Argument
#' copied from the spec.taper function in the default stats package
#'
#'
#' @return n x P array, where P = length(x) - n + 1
#'
#' @export
#'
#' @examples
#' x <- rnorm(n = 20)
#' a <- swdft(x, n = 2^3)
#'
swdft <- function(x, n, type="fftw", pad=TRUE, taper='none', p=.1) {
  ## Optionally pad the array x with 0's
  if (pad == TRUE) { x <- c(rep(0, n-1), x) }

  ## Optionally create a taper to apply at each window position
  taper <- get_taper(n, taper, p)

  ## Call the appropriate wrapper FFT function
  if (type == "fftw") {
    if (requireNamespace("fftwtools", quietly = TRUE)) {
      a <- swdft_fftw(x, n, taper) # Run with the 'fftwtools' library
    } else {
      a <- swdft_fft(x, n, taper) # Run with base R's 'fft' function
    }
  } else if (type == "multitaper") {
    a <- swdft_multitaper(x=x, n=n)
  } else {
    stop("Only works for type = 'fftw'")
  }

  return(a)
}

#' Sliding Window Discrete Fourier Transform with base R
#'
#' @param x real or complex vector
#' @param n window size
#' @param taper
#'
#' @return n x P array, where P = length(x) - n + 1
#'
swdft_fft <- function(x, n, taper) {
  N <- length(x)
  P <- N - n + 1
  a <- array(data = NA, dim = c(n, P))

  for (p in n:N) {
    a[, p - n + 1] <- stats::fft(z = x[(p - n + 1):p] * taper)
  }

  return(a)
}


swdft_multitaper <- function(x, n) {
  N <- length(x)
  P <- N - n + 1
  a <- array(data = NA, dim = c(n, P))

  for (p in n:N) {
    a_p <- multitaper::spec.mtm(timeSeries=ts(x[(p - n + 1):p]), plot=FALSE)
    a[, p - n + 1] <- a_p$spec[1:n]
  }

  return(a)
}

#' Sliding Window Discrete Fourier Transform using fftw
#'
#' @param x real or complex vector
#' @param n window size
#' @param taper
#'
#' @return n x P array, where P = length(x) - n + 1
#'
swdft_fftw <- function(x, n, taper) {
  N <- length(x)
  P <- N - n + 1
  a <- array(data = NA, dim = c(n, P))

  for (p in n:N) {
    a[, p - n + 1] <- fftwtools::fftw(data = x[(p - n + 1):p] * taper)
  }

  return(a)
}
