#' Dirichlet Kernel
#'
#' @param x
#' @param n size of Dirichlet kernel
#' @param weight logical whether to add the Dirichlet Weight factor
#'
#' @return evaluation of the Dirichlet Kernel
#'
dirichlet_kernel <- function(x, n, weight=FALSE) {
  if (x %% ((2 * pi)) == 0) {
    if ( (x / (2 * pi) ) %% 2  == 1 && zero_start == FALSE) { return (-n) }
    return(n)
  }

  val <- sin( (n * x) / 2) / sin( x / 2 )

  # Optionally convert to the Dirichlet Weight (DW)
  if (weight == TRUE) {
    i <- complex(length.out = 1, real = 0, imaginary = 1)
    shift <- exp(i * ((n - 1) / 2) * x)
    val <- shift * val
  }

  return(val)
}
