#' Complex Demodulation
#'
#' @param x vector (time-series)
#' @param f0 Fixed frequency of the cosine function
#' @param smooth Type of smoothing method. Only accepts 'spline'
#' @param spar sparsity parameter in smooth.spline
#'
#' @return length(x) demodulated time-series
#'
#' @references Bloomfield: Fourier Analysis of Time-Series, Chapter 7
#'
complex_demod <- function(x, f0, smooth="spline", spar=.8) {
  ## Demodulate the original time-series
  N <- length(x)
  t <- 0:(N-1)
  demod <- complex(real=cos(2 * pi * f0 * t), imaginary=-sin(2 * pi * f0 * t))
  y <- x * demod

  ## Smooth the demodulated series
  if (smooth == "spline") {
    re_smooth <- smooth.spline(x=Re(y), spar=spar)$y
    im_smooth <- smooth.spline(x=Im(y), spar=spar)$y
    y_smooth <- complex(real=re_smooth, imaginary=im_smooth)
  } else {
    stop("smooth only accepts 'spline' for now")
  }

  ## Extract the amplitude and phase, and fit the time-varying function
  A_t <- 2 * Mod(y_smooth)
  Phi_t <- Arg(y_smooth)

  return( A_t * cos((2 * pi * f0 * t) + Phi_t) )
}
