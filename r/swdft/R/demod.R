#' Complex Demodulation
#'
#' @param x numeric vector
#' @param f0 numeric scalar. Frequency to demodulate
#' @param smooth. character. Type of smoothing to use, accepts either 'ma', 'double_ma',
#' or 'butterworth' (the default)
#' @param order moving average parameter if 'smooth' argument equals 'ma' or 'double_ma'. Defaults to 5
#' @param passfreq numeric scalar. Pass frequency used in butterworth low-pass filter. Defaults to 2 * f0
#' @param match_swdft logical. Only used to demonstrate equivalence w/ SWDFT when
#' a moving average filter is used. Otherwise, never used.
#' @param window_size defaults to NULL, only used when match_swdft=TRUE, so can ignore.
#'
#' @return Object of class 'swdft_demod'
#'
#' @references Chapter 7 of 'Fourier Analysis of Time-Series' by Peter Bloomfield
#' The following blog-post for the idea of a butterworth filter: #' https://dankelley.github.io/r/2014/02/17/demodulation.html
#'
complex_demod <- function(x, f0, smooth="butterworth", order=5, passfreq=2*f0, match_swdft=FALSE, window_size=NULL) {
  N <- length(x)
  t <- 0:(N-1)

  ## Demodulate the original time-series
  demod <- complex(real=cos(2 * pi * f0 * t), imaginary=-sin(2 * pi * f0 * t))
  y <- x * demod

  ## Smooth the demodulated series
  if (smooth == "ma") {
    re_smooth <- swdft::moving_average(x=Re(y), order=order)
    im_smooth <- swdft::moving_average(x=Im(y), order=order)
    y_smooth <- complex(real=re_smooth, imaginary=im_smooth)

    ## Optional shift factor that sets the demodulated series equal
    ## to the SWDFT outputs. Primarily used for testing.
    if (match_swdft == TRUE) {
      l <- floor( order / 2)
      k <- f0 * window_size
      s <- (t - l) * k
      shift <- swdft::prou(n=order)^(s)
      y_smooth <- y_smooth * shift
    }

  } else if (smooth == "double_ma") {
    re_smooth <- swdft::moving_average(x=swdft::moving_average(x=Re(y), order=order), order=order)
    im_smooth <- swdft::moving_average(x=swdft::moving_average(x=Im(y), order=order), order=order)
    y_smooth <- complex(real=re_smooth, imaginary=im_smooth)

  } else if (smooth == "butterworth") {
    filter <- signal::butter(n=order, W=freqcut, type="low")
    re_smooth <- signal::filtfilt(filt=filter, Re(y))
    im_smooth <- signal::filtfilt(filt=filter, Im(y))
    y_smooth <- complex(real=re_smooth, imaginary=im_smooth)

  } else {
    stop("smooth only accepts 'ma', 'double_ma', or 'butterworth'")
  }

  ## Extract the amplitude and phase, and fit the time-varying function
  A_t <- 2 * Mod(y_smooth)
  Phi_t <- Arg(y_smooth)
  fitted <- A_t * cos((2 * pi * f0 * t) + Phi_t)

  ## Return a 'swdft_demod' S3 object
  swdft_demod_obj <- swdft::new_swdft_demod(x=x, f0=f0, A_t=A_t, Phi_t=Phi_t, fitted=fitted, y=y,
                                            y_smooth=y_smooth, smooth=smooth, order=order,
                                            passfreq=passfreq)
  return( swdft_demod_obj )
}

#' Constructor function for 'swdft_demod'
#'
#' @param x input signal
#' @param f0 center frequency to demodulate
#' @param A_t extracted amplitude from y_smooth
#' @param Phi_t extracted phase from y_smooth
#' @param fitted fitted values
#' @param y non-smoothed demodulated signal
#' @param y_smooth smoothed demodulated signal
#' @param smooth. character. Type of smoothing to use, accepts either 'ma', 'double_ma',
#' or 'butterworth' (the default)
#' @param order moving average parameter if 'smooth' argument equals 'ma' or 'double_ma'. Defaults to 5
#' @param passfreq numeric scalar. Pass frequency used in butterworth low-pass filter. Defaults to 2 * f0
#'
#' @return list with the following elements
#' \itemize{
#'   \item coefficients. A matrix of parameters, the three columns are: 1. amplitude 2. phase, and 3. frequency.
#'   There is only more that one row used when multiple frequencies are fit sequentially.
#'   \item fitted. fitted values of cosine regression model
#'   \item residuals. residuals of cosine regression model
#'   \item data. original signal used to fit cosine regression
#'   \item list with the filter used ('smooth') and parameters ('order' for 'ma' or 'double_ma', 'passfreq' for butterworth)
#'   \item list w/ the demodulated signal, and smoothed demodulated signal
#' }
#'
new_swdft_demod <- function(x, f0, A_t, Phi_t, fitted, y, y_smooth, smooth, order, passfreq) {
  structure(list(coefficients=list(f0=f0, inst_amp=A_t, inst_phase=Phi_t),
                 fitted=fitted,
                 residuals=x-fitted,
                 data=x,
                 filter=list(smooth=smooth, order=order, passfreq=passfreq),
                 demod=list(y=y, y_smooth=y_smooth)),
            class=c("swdft_demod", "swdft_mod"))
}

#' Matching Demodulation
#'
#' @param x
#' @param n
#' @param filter
#' @param thresh
#' @param max_cycles
#' @param order_prop
#' @param freqcut_scale
#' @param thresh
#' @param debug
#'
matching_demod <- function(x, n, filter="ma", thresh=.05, max_cycles=5,
                           order_prop=.02, freqcut_scale=2,
                           debug=FALSE) {
  N <- length(x)
  demods <- array(data=NA_real_, dim=c(max_cycles, N))
  amps <- array(data=NA_real_, dim=c(max_cycles, N))
  phases <- array(data=NA_real_, dim=c(max_cycles, N))
  resids <- array(data=NA_real_, dim=c(max_cycles, N))
  total_fit <- rep(0, N)
  khats <- c()

  ## Iteratively demodulate the signal until the stopping criteria is reached
  cycle <- 1
  while (cycle <= max_cycles) {
    ## Calculate the residuals of this cycle
    cycle_resids <- x - total_fit

    ## Select the frequency that explains the maximum variance proportion
    a <- swdft::swdft(x=cycle_resids, n=n, taper='cosine') * (1 / n)
    maxval <- max(Mod(a)^2)
    cat("Max SqMod: ", maxval, " \n");

    ## Break if no more large components in the SWDFT exist
    if (maxval < thresh) {
      cat("No large spectral components remain! \n")
      break;
    }

    ## Calculate the maximum frequency plus the residuals
    khat <- which(Mod(a) == max(Mod(a)), arr.ind=TRUE)[1,1]
    f0 <- (khat-1) / n

    ## Smooth the demodulated series with optional different filters
    if (filter == "ma") {
      khat_demod <- swdft::complex_demod(x=x, f0=f0, smooth='ma', order=N*order_prop)
    } else if (filter == "butterworth") {
      khat_demod <- swdft::complex_demod(x=x, f0=f0, smooth='butterworth', freqcut=freqcut_scale*f0)
    } else {
      stop("filter must be 'ma' or 'butterworth'")
    }

    ## Update the total fitter values with the newly demodulated series
    total_fit <- total_fit + khat_demod$fitted
    if (any(is.na(total_fit))) {
      warning("adding zero's to total fit, likely because MA filter used")
      total_fit[is.na(total_fit)] <- 0
    }

    ## Store the parameters of this cycles demodulated signal
    demods[cycle, ] <- khat_demod$fitted
    amps[cycle, ] <- khat_demod$amp
    phases[cycle, ] <- khat_demod$phase
    resids[cycle, ] <- cycle_resids
    khats <- c(khats, khat)

    ## Optionally plot what the current fit looks like
    if (debug == TRUE) {
      par(mfrow=c(2, 1))
      plot(x=x, cex=.4, pch=19)
      lines(total_fit, lwd=2, col="red")
      plot_swdft(a=a, freq_type="hertz", fs=1000, hertz_range=c(0, 100), col="tim.colors")
      par(mfrow=c(1, 1))
    }

    cycle <- cycle + 1
  }

  return_rows <- !(apply(X=demods, MARGIN=1, FUN=function(x) all(is.na(x))))
  return(list(demods=demods[return_rows, ],
              amps=amps[return_rows, ],
              phases=phases[return_rows, ],
              resids=resids[return_rows, ],
              khats=khats,
              total_fit=total_fit))
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
