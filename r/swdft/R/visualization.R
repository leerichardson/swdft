#' Plot method for 'swdft' object
#'
#' @param x Object of class 'swdft'. If x$a is complex-valued, it is converted to the squared
#' modulus. If x$a is real-valued, then we assume that it represents the squared
#' @param y not used, but required by plot generic function
#' @param freq_type Specify how to display the frequency axis. Either 'cycles' (default), 'angular', or 'hertz'
#' @param fs sample rate. Used if freq_type='hertz'
#' @param hertz_range integer vector, given by (low, high). Specifies the range of hertz to display and
#' is only used when freq_type='hertz'
#' @param take_log logical. Whether to take the log before plotting
#' @param log_thresh numeric. Threshold for smallest possible value. Defaults to .000001, and is
#' used to keep plots from displaying of ~ -40.
#' @param use_fields logical. Determines whether we use image.plot from the fields package, or 'image'
#' from the graphics package. The advantage of image.plot is that we get a color scale, so the default is TRUE
#' @param scale_shrink Proportion between 0 and 1 to shrink the scale
#' @param zlim Custom z range
#' @param xlab Custom x-label
#' @param ylab Custom y-label
#' @param title Custom title
#' @param cex_main how large to make the title
#' @param cex_lab how large to make the labels
#' @param custom_xaxis Defaults to NULL. Otherwise, used to change the x-axis
#' @param custom_yaxis Defaults to NULL. Otherwise, used to change the y-axis
#' @param col defauts to grayscale, can also be 'tim.colors' from fields package
#' @param display logical. Defaults to TRUE, only used for testing purposes, so it should always be TRUE.
#'
#' @return NULL
#'
plot.swdft <- function(x, y=NULL, freq_type="cycles", fs=NULL, hertz_range=NULL,
                       take_log=FALSE, log_thresh=.00001, use_fields=TRUE, scale_shrink=.9,
                       zlim=NULL, xlab="Window Position", ylab="Frequency (Cycles/Window)", title="SWDFT",
                       cex_main=1, cex_lab=1, cex_axis=1, custom_xaxis=NULL, custom_yaxis=NULL,
                       col="grayscale", display=TRUE) {
  a <- x$a
  n <- nrow(a)
  P <- ncol(a)

  ## If passed the complex-valued SWDFT, convert to the squared modulus
  if (class(a[1, 1]) == "complex") { a <- Mod(a)^2 }

  ## Optionally take the logarithm of the coefficients
  if (take_log) { a[which(a < log_thresh)] <- log_thresh; a <- log(a) }

  ## Optionally set a custom x or y axis. If none, set as the defaults
  if ( is.null(custom_xaxis) ) {
    windows <- 0:(P-1)
  } else {
    windows <- custom_xaxis
  }

  if ( is.null(custom_yaxis) ) {
    freqs <- 0:(n-1)
  } else {
    freqs <- custom_yaxis
  }

  ## Optionally specify the frequency output
  if (freq_type == "hertz") {
    if (is.null(fs)) { stop("If plotting in Hertz, must specify the sampling rate parameter 'fs'") }
    if (ylab == "Frequency (Cycles/Window)") { ylab <- "Frequency (Hz)" }

    ## Convert the Fourier Frequencies into Hertz based on the window size
    hertz <- (0:(n - 1)) * (fs / n)

    ## Optonally specify a frequency band in Hertz
    if (is.null(hertz_range)) {
      freqs <- hertz
    } else {
      freq_inds <- hertz >= hertz_range[1] & hertz <= hertz_range[2]
      freqs <- hertz[freq_inds]
      a <- a[freq_inds, ]
    }

  } else if (freq_type == "angular") {
    ## Only keep frequencies between 0 and .5
    freqs <- freqs / n
    freq_inds <- which(freqs <= .5)
    a <- a[freq_inds, ]
    freqs <- freqs[freq_inds]

    ## Change the ylab to remove cycles/window
    if (ylab == "Frequency (Cycles/Window)") { ylab <- "Frequency" }

  } else if (freq_type != "cycles") {
    stop("freq_type must be 'cycles', 'hertz', or 'angular'")
  }

  ## Optionally determine the color of swdft output
  if (col == "grayscale") {
    color <- grDevices::grey(seq(1, 0, length = 256))
  } else if (col == "tim.colors") {
    color <- fields::tim.colors(n=256)
  } else {
    stop("col must be 'grayscale' or 'tim.colors'")
  }

  # Time-frequency plot using either the fields image.plot function or base R's image function.
  if (display == TRUE) {
    if (requireNamespace("fields", quietly = TRUE) & use_fields == TRUE) {
      fields::image.plot(x=windows, y=freqs, z=t(a),
                         col=color, xlab=xlab, ylab=ylab,
                         main=title, cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis,
                         zlim=zlim, legend.shrink=scale_shrink)
    } else {
      graphics::image(x=windows, y=freqs, z=t(a),
                      col=color, xlab=xlab, ylab=ylab, main=title,
                      cex.main=cex_main, cex.lab=cex_lab, cex.axis=cex_axis,
                      zlim=zlim)
    }
  } else {
    return("display set to false, 'plot.swdft' runs without errors")
  }
}

#' Plot method for swdft_mod object
#'
#' @param x A swdft_cosreg object
#' @param y not used, but required by plot generic function
#'
plot.swdft_mod <- function(x, y, ...) {
  N <- length(x$data)
  t <- 0:(N-1)
  plot(t, x$data, main="Fitted values for 'swdft_mod' object", xlab="", ylab="", pch=19)
  lines(t, x$data)
  lines(t, x$fitted, col="red", lty=2)
}
