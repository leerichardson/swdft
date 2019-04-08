#' Smooth the spectrogram
#'
#' @param a real or complex-valued swdft. If real-valued. this is the squared modlus. If
#' complex-valued, then all of the imaginary components must be 0
#' @param ktype
#' @param m
#' @param num_convs
#'
#' @return Smooth squared modulues SWDFT coefficients
#'
smooth_swdft <- function(a, ktype='daniell', m=2, num_convs=1) {
  n <- nrow(a)
  P <- ncol(a)

  ## If a is complex-valued, convert to using the squared modulus
  if (is.complex(a)) { a <- Mod(a)^2 / n }

  ## Create the Kernel
  if ( ( ktype %in% c('daniell', 'modified.daniell') ) == FALSE) {
    stop("ktype must be either 'daniell' or 'modified.daniell'")
  }

  kern <- stats::kernel(coef=ktype, m=c(m, num_convs))
  newm <- kern$m

  ## Pre-compute the FFT of the kernel used in the convolution
  weights <- c( kern[0:newm], rep_len(0, n - (2 * newm) - 1), kern[-newm:-1])
  fft_weights <- fftwtools::fftw(data=weights)

  ## Apply the kernel smoothing across each window position
  asmooth <- apply(X=a, MARGIN=2, FUN=smooth_pgram, fft_weight=fft_weights)

  return(asmooth)
}

#' Smooth SWDFT coefficients with a convolution
#'
#' @param a real-valued length n periodogram
#' @param fft_weight optionally specify the pre-computed FFT of the weights
#'
#' @return smoothed coefficients
#'
smooth_pgram <- function(a, fft_weight=NULL) {
  Re( fftwtools::fftw(data = fftwtools::fftw(data=a) * fft_weight, inverse=1) / length(a) )
}
