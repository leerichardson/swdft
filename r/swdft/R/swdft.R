#' Sliding Window Discrete Fourier Transform (SWDFT)
#'
#' @param x real or complex vector
#' @param n integer window size.
#' @param type algorithm to implement. defaults to "fftw"
#' @param normalize how to scale SWDFT coefficients. Defaults
#' to unnormalized (1), but many definitions include either
#' (1 / n) or (1 / sqrt(n))
#'
#' @return n x P array, where P = length(x) - n + 1
#'
#' @export
#'
#' @examples
#' x <- rnorm(n = 20)
#' a <- swdft(x, n = 2^3)
#'
#'
swdft <- function(x, n, type="fftw", normalize=1) {
  if (type == "fftw") {
    if (requireNamespace("fftwtools", quietly = TRUE)) {
      a <- swdft_fftw(x, n) # Run with the 'fftwtools' library
    } else {
      a <- swdft_fft(x, n) # Run with base R's 'fft' function
    }
  } else {
    stop("Only works for type = 'fftw'")
  }

  return(a * normalize)
}

#' Sliding Window Discrete Fourier Transform with base R
#'
#' @param x real or complex vector
#' @param n window size
#'
#' @return n x P array, where P = length(x) - n + 1
#'
swdft_fft <- function(x, n) {
  N <- length(x)
  P <- N - n + 1
  a <- array(data = NA, dim = c(n, P))

  for (p in n:N) {
    a[, p - n + 1] <- stats::fft(x[(p - n + 1):p])
  }

  return(a)
}

#' Sliding Window Discrete Fourier Transform using fftw
#'
#' @param x real or complex vector
#' @param n window size
#'
#' @return n x P array, where P = length(x) - n + 1
#'
swdft_fftw <- function(x, n) {
  N <- length(x)
  P <- N - n + 1
  a <- array(data = NA, dim = c(n, P))

  for (p in n:N) {
    a[, p - n + 1] <- fftwtools::fftw(x[(p - n + 1):p])
  }

  return(a)
}
