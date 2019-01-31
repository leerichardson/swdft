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

    ## Optional shift factor that sets the demodulated series equal
    ## to the SWDFT outputs. Used primarily for testing
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
#' @param resmooth logical. Resmooth the components?
#' @param spar defaults to .8
#'
demod_swdft <- function(a, k, resmooth=FALSE, spar=.8) {
  N <- ncol(a)
  n <- nrow(a)
  l <- floor( n / 2 )
  t <- 0:(N-1)

  ## Shift the frequency band to math the demodulated series
  s <- ( (1:N) %% n) * k
  shift <- swdft::prou(n=n)^(-s)
  demod_k <- a[k + 1, ] * shift
  demod_signal <- c(rep(NA, l), demod_k[n:N], rep(NA, l))

  ## Optionally resmooth the phase and amplitude with a smoothing spline
  if (resmooth == TRUE) {
    non_na_inds <- !is.na(demod_signal)
    demod_smooth <- demod_signal[non_na_inds]
    re_smooth <- stats::smooth.spline(x=Re(demod_smooth), spar=spar)$y
    im_smooth <- stats::smooth.spline(x=Im(demod_smooth), spar=spar)$y
    demod_signal[non_na_inds] <- complex(real=re_smooth, imaginary=im_smooth)
  }

  ## Extract amplitude, phase, and fitted values
  A_t <- 2 * Mod(demod_signal)
  Phi_t <- Arg(demod_signal)
  fitted <- A_t * cos((2 * pi * (k / n) * t) + Phi_t)

  return( list(fitted=fitted, demod=demod_signal, amp=A_t, phase=Phi_t) )
}

#' Sequentially Demodulate a signal
#'
#' @param x signal to demodulate
#' @param n window size
#' @param max_cycles maximum number of iterations
#' @param debug
#' @param thresh
#'
demod_signal <- function(x, n, max_cycles=5, order_prop=.02, debug=FALSE, thresh=.05) {
  ## Create an array to store the demodulations
  N <- length(x)
  demods <- array(data=NA_real_, dim=c(max_cycles, N))
  amps <- array(data=NA_real_, dim=c(max_cycles, N))
  phases <- array(data=NA_real_, dim=c(max_cycles, N))
  khats <- c()
  total_fit <- rep(0, N)

  ## Iteratively demodulate the signal until the stopping criteria is reached
  cycle <- 1
  while (cycle <= max_cycles) {
    ## Select the frequency that explains the maximum variance proportion
    a <- swdft::swdft(x=x-total_fit, n=n, taper='cosine') * (1 / n)
    aprop <- swdft::swdft_to_props(a=a)
    maxval <- max(Mod(a)^2)
    cat("Max SqMod: ", maxval, " Max Prop: ", max(aprop, na.rm=TRUE), " \n");

    if (maxval < thresh) {
      cat("No more large components! \n")
      break;
    }

    ## Demodulate the maximum frequency and add to the total fit
    khat <- which(Mod(a) == max(Mod(a)), arr.ind=TRUE)[1,1]
    khat_demod <- swdft::complex_demod(x=x-total_fit, f0=(khat-1)/n, smooth='ma', order=(N * order_prop))
    total_fit <- total_fit + khat_demod$fitted
    total_fit[is.na(total_fit)] <- 0

    ## Store the parameters of this cycles demodulated signal
    demods[cycle, ] <- khat_demod$fitted
    amps[cycle, ] <- khat_demod$amp
    phases[cycle, ] <- khat_demod$phase
    khats <- c(khats, khat)

    ## Optionally plot what the current fit looks like
    if (debug == TRUE) {
      par(mfrow=c(2, 1))
      plot(x=x, cex=.4, pch=19)
      lines(total_fit, lwd=2, col="red")
      plot_swdft(a=a, col="other", hertz=TRUE, fs=1000, hertz_range=c(0, 100))
      par(mfrow=c(1, 1))
    }

    cycle <- cycle + 1
  }

  return_rows <- !(apply(X=demods, MARGIN=1, FUN=function(x) all(is.na(x))))
  return(list(demods=demods[return_rows, ],
              amps=amps[return_rows, ],
              phases=phases[return_rows, ],
              khats=khats,
              total_fit=total_fit))
}

#' Convert the SWDFT to proportions of frequency
#'
#' @param a swdft
swdft_to_props <- function(a) {
  n <- nrow(a)
  m <- floor(n / 2)
  amod <- Mod(a)^2

  return( apply(X=amod[2:m, ], MARGIN=2, FUN=function(x) x / sum(x)) )
}
