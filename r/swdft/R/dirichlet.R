#' Dirichlet Kernel (Weight) for arbitrary summation indices
#'
#' @param x numeric to evaluate
#' @param phase defaults to 0
#' @param a start of summation index
#' @param b end of summation index
#'
#' @return sum of a complex exponential sum
#'
dirichlet <- function(x, phase=0, a=0, b=length(x)-1) {
  i <- complex(length.out=1, real=0, imaginary=1)
  exp(x = i * (( ((a + b) * x) / 2 ) + phase) ) * swdft::dirichlet_kernel(x=x, n=b-a+1)
}

#' Dirichlet Kernel
#'
#' @param x variable evaluated by dirichlet kernel
#' @param n size of Dirichlet kernel
#' @param dw logical whether to add the Dirichlet Weight (DW) factor
#'
#' @return evaluation of the Dirichlet Kernel (D_n(x))
#'
dirichlet_kernel <- function(x, n, dw=FALSE) {
  # Special cases to make the Dirichlet Kernel continuous.
  if ( (x %% pi) == 0) {
    x_nopi <- round(x / pi, digits=0)
    if (x_nopi %% 2 != 0) {
      return(0)
    } else if ( (x_nopi %% 4) == 0) {
      return(n)
    } else {
      return(-n)
    }
  }

  val <- sin( (n * x) / 2) / sin( x / 2 )

  # Optionally add the Dirchlet Weight factor
  if (dw == TRUE) {
    i <- complex(length.out = 1, real = 0, imaginary = 1)
    shift <- exp(i * ((n - 1) / 2) * x)
    val <- shift * val
  }

  return(val)
}
