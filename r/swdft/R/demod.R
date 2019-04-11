#' Complex Demodulation
#'
#' @param x numeric vector
#' @param f0 numeric scalar. Frequency to demodulate
#' @param smooth character. Type of smoothing to use, accepts either 'ma', 'double_ma',
#' or 'butterworth' (the default)
#' @param order moving average parameter if 'smooth' argument equals 'ma' or 'double_ma'. Defaults to 5
#' @param passfreq numeric scalar. Pass frequency used in butterworth low-pass filter. Defaults to .1
#' which corresponds to a pass frequency of 2 * f0.
#' @param match_swdft logical. Only used to demonstrate equivalence w/ SWDFT when
#' a moving average filter is used. Otherwise, never used.
#' @param window_size defaults to NULL, only used when match_swdft=TRUE, so can ignore.
#'
#' @export
#'
#' @return An S3 'swdft_demod' object. See ?new_swdft_matching_demod for details.
#'
#' @references Chapter 7 of 'Fourier Analysis of Time-Series' by Peter Bloomfield
#' The following blog-post for the idea of a butterworth filter: #' https://dankelley.github.io/r/2014/02/17/demodulation.html
#'
complex_demod <- function(x, f0, smooth="butterworth", order=5, passfreq=.1, match_swdft=FALSE, window_size=NULL) {
  N <- length(x)
  t <- 0:(N-1)

  ## Demodulate the original time-series
  demod <- complex(real=cos(2 * pi * f0 * t), imaginary=-sin(2 * pi * f0 * t))
  y <- x * demod

  ## Smooth the demodulated series
  if (smooth == "ma") {
    re_smooth <- moving_average(x=Re(y), order=order)
    im_smooth <- moving_average(x=Im(y), order=order)
    y_smooth <- complex(real=re_smooth, imaginary=im_smooth)

    ## Optional shift factor that sets the demodulated series equal
    ## to the SWDFT outputs. Primarily used for testing.
    if (match_swdft == TRUE) {
      l <- floor( order / 2)
      k <- f0 * window_size
      s <- (t - l) * k
      shift <- prou(n=order)^(s)
      y_smooth <- y_smooth * shift
    }

  } else if (smooth == "double_ma") {
    re_smooth <- moving_average(x=moving_average(x=Re(y), order=order), order=order)
    im_smooth <- moving_average(x=moving_average(x=Im(y), order=order), order=order)
    y_smooth <- complex(real=re_smooth, imaginary=im_smooth)

  } else if (smooth == "butterworth") {
    filter <- signal::butter(n=order, W=passfreq, type="low")
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
  swdft_demod_obj <- new_swdft_demod(x=x, f0=f0, A_t=A_t, Phi_t=Phi_t, fitted=fitted, y=y,
                                            y_smooth=y_smooth, smooth=smooth, order=order,
                                            passfreq=passfreq)
  return( swdft_demod_obj )
}

#' Constructor function for class 'swdft_demod'
#'
#' @inheritParams complex_demod
#' @param A_t extracted amplitude from y_smooth
#' @param Phi_t extracted phase from y_smooth
#' @param fitted fitted values
#' @param y non-smoothed demodulated signal
#' @param y_smooth smoothed demodulated signal
#' @param passfreq numeric frequency used as the passfreq in the low-pass filter
#'
#' @export
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
                 smooth=list(smooth=smooth, order=order, passfreq=passfreq),
                 demod=list(y=y, y_smooth=y_smooth)),
            class=c("swdft_demod", "swdft_mod"))
}

#' Matching Demodulation
#'
#' @param x numeric. Signal to demodulate
#' @param n integer. Window size for SWDFT
#' @param thresh numeric. Threshold to determine whether to continue demodulating
#' @param max_cycles maximum number of demodulation cycles
#' @param smooth character. Type of smoothing to use, accepts either 'ma', 'double_ma',
#' or 'butterworth' (the default)
#' @param order moving average parameter if 'smooth' argument equals 'ma' or 'double_ma'. Defaults to 5
#' @param passfreq numeric scalar. Pass frequency used in butterworth low-pass filter. defaults to .1
#' @param debug Logical. Whether to print out intermediate output.
#'
#' @export
#'
#' @return An S3 'swdft_matching_demod' object. See ?new_swdft_matching_demod for details.
#'
matching_demod <- function(x, n, thresh=.05, max_cycles=5, smooth="butterworth", order=5, passfreq=.1, debug=FALSE) {
  # --- Generate arrays and vectors to store the final results --
  N <- length(x)
  demods <- list(y=array(data=NA_real_, dim=c(max_cycles, N)), y_smooth=array(data=NA_real_, dim=c(max_cycles, N)))
  amps <- array(data=NA_real_, dim=c(max_cycles, N))
  phases <- array(data=NA_real_, dim=c(max_cycles, N))
  resids <- array(data=NA_real_, dim=c(max_cycles, N))
  fits <- array(data=NA_real_, dim=c(max_cycles, N))
  fitted <- rep(0, N)
  maxvals <- c()
  passfreqs <- c()
  freqs <- c()

  # --- Iteratively demodulate the signal until stopping criteria is reached ---
  cycle <- 0
  while (cycle <= max_cycles) {
    ## Select the frequency that corresponds to the largest SWDFT frequency
    cycle_resids <- x - fitted

    ## Take the SWDFT based on whether we are debugging or not
    if (debug == TRUE) {
      a_debug <- swdft(x=cycle_resids, n=n, taper_type='cosine')
      a <- a_debug$a * (1 / n)
    } else {
      a <- swdft(x=cycle_resids, n=n, taper_type='cosine')$a * (1 / n)
    }

    ## Get the largest SWDFT coefficient
    maxval <- max(Mod(a)^2)
    maxvals <- c(maxvals, maxval)
    if (debug == TRUE) { cat("Max SqMod: ", maxval, "in iteration ", cycle, " \n") }

    ## Break if no SWDFT coefficients exceed a threshold
    if (maxval < thresh) {
      if (debug == TRUE) { cat("No large spectral components remain! \n") }
      break;
    } else {
      cycle <- cycle + 1
    }

    ## Calculate the maximum frequency plus the residuals
    khat <- which(Mod(a) == max(Mod(a)), arr.ind=TRUE)[1,1]
    f0 <- (khat-1) / n

    ## Apply complex demodulation at the frequency w/ the largest SWDFT coefficient
    if (smooth == "ma") {
      khat_demod <- complex_demod(x=x, f0=f0, smooth='ma', order=order)
    } else if (smooth == "butterworth") {
      khat_demod <- complex_demod(x=cycle_resids, f0=f0, smooth='butterworth', passfreq=passfreq)
    } else {
      stop("filter must be 'ma' or 'butterworth'")
    }

    ## Update the total fitter values with the newly demodulated series
    fitted <- fitted + khat_demod$fitted
    if (any(is.na(fitted))) {
      if (debug == TRUE) {  warning("adding zero's to total fit, likely because MA filter used") }
      fitted[is.na(fitted)] <- 0
    }

    ## Store the parameters from this cycle of complex demodulation
    demods$y[cycle, ] <- khat_demod$demod$y
    demods$y_smooth[cycle, ] <- khat_demod$demod$y_smooth
    amps[cycle, ] <- khat_demod$coefficients$inst_amp
    phases[cycle, ] <- khat_demod$coefficients$inst_phase
    resids[cycle, ] <- cycle_resids
    fits[cycle, ] <- khat_demod$fitted
    freqs <- c(freqs, f0)
    passfreqs <- c(passfreqs, passfreq)

    ## Optionally plot the current fit and SWDFT of the resids
    if (debug == TRUE) {
      graphics::par(mfrow=c(2, 1))
      graphics::plot(x=x, cex=.4, pch=19)
      graphics::lines(fitted, lwd=2, col="red")
      plot(x=a_debug, freq_type="fraction", col="tim.colors")
      graphics::par(mfrow=c(1, 1))
    }
  }

  ## Return an S3 object of class 'swdft_matching_demod'
  return_rows <- !(apply(X=demods$y, MARGIN=1, FUN=function(x) all(is.na(x))))
  swdft_matching_demod_obj <- new_swdft_matching_demod(x, n, fitted, thresh, max_cycles, smooth,
                                                       order, passfreqs, maxvals, freqs, amps, phases,
                                                       demods, cycle, resids, fits, return_rows)

  return( swdft_matching_demod_obj )
}

#' Constructor function for class 'swdft_matching_demod'
#'
#' @inheritParams matching_demod
#' @param fitted fitted values
#' @param passfreqs pass frequency used in each iteration
#' @param maxvals Maximum SWDFT coefficient for each iteration
#' @param freqs Frequencies used in each iteration
#' @param amps Instantaneous amplitudes for each iteration
#' @param phases Instantaneous phases for each iteration
#' @param demods List of demodulated signal and smoothed demodulated signal for each iteration
#' @param cycle Number of cycles used
#' @param resids Residuals for each iteration
#' @param fits Fitted values for each iteration
#' @param return_rows Logival vector indicating which iterations occured. Used for subsetting.
#'
#' @export
#'
#' @return list with the following elements
#' \itemize{
#'   \item coefficients. coefficients from the R local signals with time-varying ampltidue and phase model.
#'   \item fitted. fitted values of cosine regression model
#'   \item residuals. residuals of cosine regression model
#'   \item data. original signal used to fit cosine regression
#'   \item smooth. list with the filter used ('smooth') and parameters ('order' for 'ma' or 'double_ma', 'passfreq' for butterworth)
#'   \item demod. list w/ the demodulated signal, and smoothed demodulated signal
#'   \item thresh. Threshold used.
#'   \item iterations. List of fits, residuals, and maximum values for each iteration
#' }
#'
new_swdft_matching_demod <- function(x, n, fitted, thresh, max_cycles, smooth, order, passfreqs,
                                     maxvals, freqs, amps, phases, demods, cycle, resids, fits,
                                     return_rows) {
  structure(list(coefficients=list(R=cycle, f0=freqs, inst_amp=amps[return_rows, ], inst_phase=phases[return_rows, ]),
                 fitted=fitted,
                 residuals=x-fitted,
                 data=x,
                 smooth=list(smooth=smooth, order=order, passfreq=passfreqs),
                 demod=list(y=demods$y[return_rows, ], y_smooth=demods$y_smooth[return_rows, ]),
                 thresh=thresh,
                 iterations=list(num_iters=cycle, iter_fits=fits[return_rows, ],
                                 iter_resids=resids[return_rows, ], maxvals=maxvals)),
            class=c("swdft_matching_demod", "swdft_demod", "swdft_mod"))
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
  shift <- prou(n=n)^(-s)
  demod_k <- a[k + 1, ] * shift
  demod_signal <- c(rep(NA, l), demod_k[n:N], rep(NA, l))

  ## Extract amplitude, phase, and fitted values
  A_t <- 2 * Mod(demod_signal)
  Phi_t <- Arg(demod_signal)
  fitted <- A_t * cos((2 * pi * (k / n) * t) + Phi_t)

  return( list(fitted=fitted, demod=demod_signal, amp=A_t, phase=Phi_t) )
}
