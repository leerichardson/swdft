#' Sliding Window Discrete Fourier Transform
#'
#' @param x real or complex vector
#' @param n integer window size.
#' @param type algorithm to implement. defaults to "fftw", other option 'fft' for R's base FFT function.
#' R's base fft function is used if
#' @param pad optionally zero-pad the array to that the output
#' array has the same dimension as the original time-series
#' @param taper_type type of taper for each window position. defaults to 'none', can also be 'cosine'.
#' @param p Proportion to be tapered at each end of the series. Argument
#' copied from the spec.taper function in the default stats package. Defaults to .1.
#' @param smooth Type of smoother. Defaults to 'none', can also be 'daniell' or 'modified daniell'.
#' If smooth is 'none', then the SWDFT returns the smoothed squared modulus coefficients, not the  complex numbers
#' @param m width of kernel. Defaults to 2
#' @param num_convs Number of times to convolve the kernel. Defailts to 1
#'
#' @return An S3 'swdft' object. See ?new_swdft for details.
#'
#' @export
#'
#' @examples
#' x <- rnorm(n = 20)
#' a <- swdft(x, n = 2^3)
#'
swdft <- function(x, n, type="fftw", pad=TRUE, taper_type='none', p=.1, smooth='none', m=2, num_convs=1) {
  ## Optionally pad the array x with 0's
  if (pad == TRUE) { x <- c(rep(0, n-1), x) }

  ## Generate taper to be used by the input data
  taper <- get_taper(n, taper_type, p)

  ## Call the appropriate SWDFT wrapper
  if (type == "fftw") {
    if (requireNamespace("fftwtools", quietly = TRUE)) {
      a <- swdft_fftw(x, n, taper) # Run with the 'fftwtools' library
    } else {
      type <- "fft" ## Update the type of FFT used
      a <- swdft_fft(x, n, taper) # Run with base R's 'fft' function if 'fftwtools' not available
    }
  } else if (type == "fft") {
    a <- swdft_fft(x, n, taper) # Run with base R's 'fft' function
  } else {
    stop("Only works for type = 'fftw' or 'fft'")
  }

  ## Optionally smooth the final coefficients
  if (smooth != 'none') { a <- smooth_swdft(a=a, ktype=smooth, m=m, num_convs=num_convs) }

  ## Return a 'swdft' S3 object
  swdft_obj <- new_swdft(a=a, x=x, n=n, type=type, pad=pad, taper_type=taper_type,
                                taper=taper, p=p, smooth=smooth, m=m, num_convs=num_convs)
  return( swdft_obj )
}

#' Constructor function for class 'swdft'
#'
#' @param a 2D complex array of SWDFT coefficients. If there is smoothing, then
#' this represents the smoothed squared modulus coefficients.
#' @param x numeric input signal
#' @param n window size
#' @param type 'fftw' or 'fft'
#' @param pad whether or not it was padded
#' @param taper_type type of taper
#' @param taper numeric values of the taper
#' @param p of cosine taper (if used)
#' @param smooth type of smoother
#' @param m width of kernel for smoothing (optional)
#' @param num_convs number of kernel convolutions (optional)
#'
#' @return list w/ the same elements as the arguments, an S3 object of class 'swdft'
#'
new_swdft <- function(a, x, n, type, pad, taper_type, taper, p, smooth, m, num_convs) {
  structure(list(a=a, x=x, n=n, type=type, pad=pad, taper_type=taper_type, taper=taper, p=p,
                 smooth=smooth, m=m, num_convs=num_convs), class="swdft")
}

#' Sliding Window Discrete Fourier Transform with base R
#'
#' @inheritParams swdft
#' @param taper length n vector to multiply against the input data for each window position
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

#' Sliding Window Discrete Fourier Transform using fftw
#'
#' @inheritParams swdft
#' @inheritParams swdft_fft
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
