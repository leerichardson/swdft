#' Dirichlet Kernel
#'
#' @param x variable evaluated by dirichlet kernel
#' @param n size of Dirichlet kernel
#' @param dw logical whether to add the Dirichlet Weight (DW) factor
#'
#' @return evaluation of the Dirichlet Kernel (D_n(x))
#'
dirichlet_kernel <- function(x, n, dw=FALSE) {
  if (x %% ((2 * pi)) == 0) {
    if ( (x / (2 * pi) ) %% 2  == 1 ) { return (-n) }
    return(n)
  }

  val <- sin( (n * x) / 2) / sin( x / 2 )

  # Optionally convert to the Dirichlet Weight (DW)
  if (dw == TRUE) {
    i <- complex(length.out = 1, real = 0, imaginary = 1)
    shift <- exp(i * ((n - 1) / 2) * x)
    val <- shift * val
  }

  return(val)
}
