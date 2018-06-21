#' Time-Frequency Plot of the SWDFT
#'
#' @param a 2D complex array of output of the 'swdft' function
#' @param type type of complex-number ourput. Either 'Re',
#' 'Im', Arg', or 'Mod'. Defaults to 'Mod'
#' @param use_fields logical that determines whether we use image or
#' image.plot from the fields package. The key advantage of fields is
#' that it automatically provides a legend
#'
#' @export
#'
#' @examples
#' x <- rnorm(n = 100)
#' a <- swdft(x, n = 2^5)
#' plot_swdft(a)
#'
plot_swdft <- function(a, type="Mod", use_fields=TRUE, zlim=NULL,
                    xlab="Window Position", ylab="Frequency (Cycles/Window)",
                    title="SWDFT") {
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

  # Extract parameters used in plotting
  P <- ncol(a)
  n <- nrow(a)
  N <- P + n - 1
  freqs <- 0:(n - 1)
  windows <- (n - 1):(N - 1)
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
          main=title, cex.main = 1, zlim=zlim)
  }

}
