#' 3D Sliding Window Discrete Fourier Transform
#'
#' @param x 3D real or complex-valued array
#' @param n0 window size in dimension 0
#' @param n1 window size in dimension 1
#' @param n2 window size in dimension 2
#' @param type detaults to 'base', which is the only option
#'
#' @export
#'
#' @return An S3 'swdft3d' object. See ?new_swdft for details.
#'
swdft3d <- function(x, n0, n1, n2, type="base") {
  if (type == "base") {
    a <- swdft_base_3d(x, n0, n1, n2)
  } else {
    stop("type for 3D algorithm must be 'base' (fftwtools doesn't have a 3D program)")
  }

  swdft3d_obj <- new_swdft3d(a=a, x=x, n0=n0, n1=n1, n2=n2, type=type)
  return( swdft3d_obj )
}

#' Constructor function for class 'swdft3d'
#'
#' @param a 4D complex-valued array of 2D SWDFT coefficients
#' @param x 3D real or complex-valued array
#' @param n0 window size in dimension 0
#' @param n1 window size in dimension 1
#' @param n2 window size in dimension 2
#' @param type detaults to 'base', which is the only option
#'
#' @return S3 object w/ the same elements as arguments to this constructor function
#'
new_swdft3d <- function(a, x, n0, n1, n2, type) {
  structure(list(a=a, x=x, n0=n0, n1=n1, n2=n2, type=type), class=c("swdft3d", "swdft"))
}

#' 3D SWDFT using base R
#'
#' @param x 3D real or complex-valued array
#' @param n0 window size in dimension 0
#' @param n1 window size in dimension 1
#' @param n2 window size in dimension 2
swdft_base_3d <- function(x, n0, n1, n2) {
  N0 <- dim(x)[1]
  N1 <- dim(x)[2]
  N2 <- dim(x)[3]
  P0 <- N0 - n0 + 1
  P1 <- N1 - n1 + 1
  P2 <- N2 - n2 + 1
  a <- array(data=NA, dim=c(n0, n1, n2, P0, P1, P2))

  for (p0 in n0:N0) {
    for (p1 in n1:N1) {
      for (p2 in n2:N2) {
        a[, , , p0 - n0 + 1, p1 - n1 + 1, p2 - n2 + 1] <- stats::fft(x[(p0 - n0 + 1):p0, (p1 - n1 + 1):p1, (p2 - n2 + 1):p2])
      }
    }
  }

  return(a)
}
