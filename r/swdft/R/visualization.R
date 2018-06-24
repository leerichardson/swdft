#' Time-Frequency Plot of the SWDFT
#'
#' @param a 2D complex array of output of the 'swdft' function
#' @param type type of complex-number ourput. Either 'Re',
#' 'Im', Arg', or 'Mod'. Defaults to 'Mod'
#' @param take_log logical. whether to take the log before plotting
#' @param use_fields logical that determines whether we use image or
#' image.plot from the fields package. The key advantage of fields is
#' that it automatically provides a legend
#' @param zlim Custom z range
#' @param xlab Custom x-label
#' @param ylab Custom y-label
#' @param title Custom title
#'
#' @export
#'
#' @examples
#' x <- rnorm(n = 100)
#' a <- swdft(x, n = 2^5)
#' plot_swdft(a)
#'
plot_swdft <- function(a, type="Mod", take_log=FALSE, log_thresh=.01,
                       only_unique=FALSE, use_fields=TRUE, pad_array=FALSE,
                       zlim=NULL, xlab="Window Position",
                       ylab="Frequency (Cycles/Window)", title="SWDFT") {
  # Construct the correct complex-number output to plot
  if (class(a[1, 1]) != "complex") {
    warning("Expecting complex array from 'swdft'")
  } else if (type == 'Mod') {
    a <- Mod(a)^2
  } else if (type == "Re") {
    a <- Re(a)
  } else if (type == "Im") {
    a <- Im(a)
  } else if (type == "Arg") {
    a <- Arg(a)
  }

  # Optionally convert to log scale
  if (take_log) {
    a[which(a < log_thresh)] <- log_thresh
    a <- log(a)
  }

  if (only_unique) {
    m <- floor(nrow(a) / 2) + 1
    a <- a[2:m, ]
  }

  # Extract parameters used in plotting
  P <- ncol(a)
  n <- nrow(a)
  N <- P + n - 1
  if (only_unique) { freqs <- 1:n } else { freqs <- 0:(n - 1) }
  if (pad_array) {
    apad <- matrix(data = NA, nrow = n, ncol = n - 1)
    a <- cbind(apad, a)
    windows <- 0:(N - 1)
  } else {
    windows <- (n - 1):(N - 1)
  }

  grayscale <- grey(seq(0, 1, length = 256))

  # Time-frequency plot using either the fields image.plot function or
  # base R's image function. T
  if (requireNamespace("fields", quietly = TRUE) & use_fields == TRUE) {
    # If fields package is available, use 'image.plot' to add a legend
    fields::image.plot(x=windows, y=freqs, z=t(a),
                       col=grayscale, xlab=xlab, ylab=ylab,
                       main=title, cex.main = 1, zlim=zlim)
  } else {
    image(x=windows, y=freqs, z=t(a),
           col=grayscale, xlab=xlab, ylab=ylab,
           main=title, cex.main = 1)
  }

}



#' Multivariate Time-Series plot of the SWDFT
#'
#' @param a 2D complex array of output of the 'swdft' function
#' @param legend logical whether to include a legend
#' @param title Custom title
#' @param cex_leg size of the legend
#' @param ylim Custom y-limits
#' @param colors custom vector of colots
#'
#' @examples
#' x <- rnorm(n = 100)
#' a <- swdft(x, n = 2^5)
#' plot_mvts(a)
#'
plot_mvts <- function(a, legend=FALSE, title="Title Left Blank", cex_leg=1, ylim=NULL, colors=NULL) {
  n <- nrow(a)
  m <- floor(n / 2) + 1
  P <- ncol(a)
  N <- P + n - 1

  if (class(a[1,1]) == "complex") { cat("Converting to Squared Modulus \n"); a <- Mod(a)^2 }
  if (is.null(ylim)) { ylim = c(min(a) - .1, max(a) + .1) }
  if (is.null(colors)) { colors <-randomcoloR::distinctColorPalette(k = n, altCol = TRUE) }

  plot((n - 1):(N - 1), a[1, ], ylim = ylim, col = colors[1], type = "l",
       xlim = c(-20, N - 1), xlab = "Window Position",
       xaxt = 'n', ylab = "", main = title, cex.main = 1, lwd = 1.2, lty = 1)
  axis(1) #, at = seq(from = n - 1, to = N - 1, by = 10))

  for (i in 2:m) {
    lines((n - 1):(N - 1), a[i, ], col = colors[i], lwd = 1.2, lty = i)
  }

  if (legend) {
    legend(x=-20, y=ylim[2], paste0("Freq: ", 0:(m - 1)), col = colors,
           cex = cex_leg, lwd = 2, lty = 1:m)
  }
}
