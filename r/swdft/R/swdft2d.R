#' 2D Sliding Window Discrete Fourier Transform
#'
#' @param x 2D input signal
#' @param n0 window size in row direction
#' @param n1 window size in column direction
#' @param type algorithm to implement. defaults to "fftw", other option 'fft' for R's base FFT function.
#' R's base fft function is used if 'fftwtools' library is not installed.
#'
#' @export
#'
#' @return An S3 'swdft2d' object. See ?new_swdft for details.
#'
swdft2d <- function(x, n0, n1, type="fftw") {
  if (type == "fftw") {
    if (requireNamespace("fftwtools", quietly = TRUE)) {
      a <- swdft2d_fftw(x, n0, n1) # Run with the 'fftwtools' library
    } else {
      print("Using base fft, since fftwtools package is not available")
      a <- swdft2d_fft(x, n0, n1) # Run with base R's 'fft' function
    }
  } else if (type == "fft") {
    a <- swdft2d_fft(x, n0, n1)
  } else {
    stop("type must be 'base' or 'fftw'")
  }

  ## Return a 'swdft2d' S3 object
  swdft2d_obj <- new_swdft2d(a=a, x=x, n0=n0, n1=n1, type=type)
  return( swdft2d_obj )

  return(a)
}

#' Constructor function for class 'swdft2d'
#'
#' @param a 4D complex-valued array of 2D SWDFT coefficients
#' @param x 2D real or complex valued signal
#' @param n0 window size in row direction
#' @param n1 window size in column direction
#' @param type algorithm to implement. defaults to "fftw", other option 'fft' for R's base FFT function.
#' R's base fft function is used if
#'
#' @return S3 object w/ the same elements as arguments to this constructor function
#'
new_swdft2d <- function(a, x, n0, n1, type) {
  structure(list(a=a, x=x, n0=n0, n1=n1, type=type), class=c("swdft2d", "swdft"))
}

#' 2D Sliding Window Discrete Fourier Transform using fftw
#'
#' @param x 2D input signal
#' @param n0 window size in row direction
#' @param n1 window size in column direction
#'
swdft2d_fftw <- function(x, n0, n1) {
  N0 <- dim(x)[1]
  N1 <- dim(x)[2]
  P0 <- N0 - n0 + 1
  P1 <- N1 - n1 + 1
  a <- array(data=NA, dim=c(n0, n1, P0, P1))

  for (p0 in n0:N0) {
    for (p1 in n1:N1) {
      a[, , p0 - n0 + 1, p1 - n1 + 1] <- fftwtools::fftw2d(x[(p0 - n0 + 1):p0, (p1 - n1 + 1):p1])
    }
  }

  return(a)
}

#' 2D Sliding Window Discrete Fourier Transform using base R
#'
#' @param x 2D input signal
#' @param n0 window size in row direction
#' @param n1 window size in column direction
#'
swdft2d_fft <- function(x, n0, n1) {
  N0 <- dim(x)[1]
  N1 <- dim(x)[2]
  P0 <- N0 - n0 + 1
  P1 <- N1 - n1 + 1
  a <- array(data=NA, dim=c(n0, n1, P0, P1))

  for (p0 in n0:N0) {
    for (p1 in n1:N1) {
      a[, , p0 - n0 + 1, p1 - n1 + 1] <- stats::fft(x[(p0 - n0 + 1):p0, (p1 - n1 + 1):p1])
    }
  }

  return(a)
}
