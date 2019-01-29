#' Complex Demodulation
#'
#' @param x vector (time-series)
#' @param f0 Fixed frequency of the cosine function
#' @param smooth Type of smoothing method. Only accepts 'spline'
#' @param spar sparsity parameter in smooth.spline
#' @param order moving average parameter if 'smooth' argument equals 'ma'
#'
#' @return length(x) demodulated time-series
#'
#' @references Bloomfield: Fourier Analysis of Time-Series, Chapter 7
#'
complex_demod <- function(x, f0, smooth="spline", spar=.8, order=5, match_swdft=FALSE) {
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

  } else if (smooth == "ma") {
    re_smooth <- swdft::moving_average(x=Re(y), order=order)
    im_smooth <- swdft::moving_average(x=Im(y), order=order)
    y_smooth <- complex(real=re_smooth, imaginary=im_smooth)

    ## Optiona: Apply shift factor that sets the demodulated series equal
    ## to the SWDFT outputs
    if (match_swdft == TRUE) {
      l <- floor( order / 2)
      k <- f0 * n
      s <- (t - l) * k
      shift <- swdft::prou(n=order)^(s)
      y_smooth <- y_smooth * shift
    }

  } else {
    stop("smooth only accepts 'spline' for now")
  }

  ## Extract the amplitude and phase, and fit the time-varying function
  A_t <- 2 * Mod(y_smooth)
  Phi_t <- Arg(y_smooth)
  fitted <- A_t * cos((2 * pi * f0 * t) + Phi_t)

  return( list(fitted=fitted, y=y, y_smooth=y_smooth, amp=A_t, phase=Phi_t) )
}

#' Demodulate a Fourier Frequency with the SWDFT
#'
#' @param a swdft
#' @param k frequency to demodulate
#'
demod_swdft <- function(a, k) {
  N <- ncol(a)
  n <- nrow(a)
  l <- floor( n / 2 )
  t <- 0:(N-1)

  ## Shift the frequency band to math the demodulated series
  s <- ( (1:N) %% n) * k
  shift <- swdft::prou(n=n)^(-s)
  demod_k <- a[k + 1, ] * shift
  demod_signal <- c(rep(NA, l), demod_k[n:N], rep(NA, l))

  ## Extract amplitude, phase, and fitted values
  A_t <- 2 * Mod(demod_signal)
  Phi_t <- Arg(demod_signal)
  fitted <- A_t * cos((2 * pi * (k / n) * t) + Phi_t)

  return( list(fitted=fitted, demod=demod_signal, amp=A_t, phase=Phi_t) )
}
