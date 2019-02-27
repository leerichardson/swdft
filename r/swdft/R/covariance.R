#' Covariance between two complex-numbered outputs
#'
#' @param k
#' @param l
#' @param delta
#' @param n
#' @param sigma
#'
#' @return complex-valued number of the covariance
#'
cov_swdft_cnum <- function(k, l, delta, n, sigma) {
  shift <- swdft::dirichlet(x=(2*pi*(l-k))/n, phase=(-2*pi*delta*l)/n, a=delta, b=n-1)
  return( (sigma^2 / n) * shift)
}

cov_swdft_rr <- function(k, l, delta, n, sigma) {
  shift1 <- Re( swdft::dirichlet(x=(-2*pi*(k-l))/n, phase=(-2*pi*delta*l)/n, a=delta, b=n-1) )
  shift2 <- Re( swdft::dirichlet(x=(-2*pi*(k+l))/n, phase=(2*pi*delta*l)/n, a=delta, b=n-1) )
  return( (sigma^2 / (2*n)) * (shift1 + shift2) )
}

cov_swdft_ii <- function(k, l, delta, n, sigma) {
  shift1 <- Re( swdft::dirichlet(x=(-2*pi*(k-l))/n, phase=(-2*pi*delta*l)/n, a=delta, b=n-1) )
  shift2 <- Re( swdft::dirichlet(x=(-2*pi*(k+l))/n, phase=(2*pi*delta*l)/n, a=delta, b=n-1) )
  return( (sigma^2 / (2*n)) * (shift1 - shift2) )
}

cov_swdft_ri <- function(k, l, delta, n, sigma) {
  shift1 <- Im( swdft::dirichlet(x=(-2*pi*(k+l))/n, phase=(2*pi*delta*l)/n, a=delta, b=n-1) )
  shift2 <- Im( swdft::dirichlet(x=(-2*pi*(k-l))/n, phase=(-2*pi*delta*l)/n, a=delta, b=n-1) )
  return( (sigma^2 / (2*n)) * (shift1 - shift2) )
}

cov_swdft_ir <- function(k, l, delta, n, sigma) {
  shift1 <- Im( swdft::dirichlet(x=(-2*pi*(k+l))/n, phase=(2*pi*delta*l)/n, a=delta, b=n-1) )
  shift2 <- Im( swdft::dirichlet(x=(-2*pi*(k-l))/n, phase=(-2*pi*delta*l)/n, a=delta, b=n-1) )
  return( (sigma^2 / (2*n)) * (shift1 + shift2) )
}

#' #' Covariance between two complex-numbered outputs
#' #'
#' #' @param k
#' #' @param l
#' #' @param delta
#' #' @param n
#' #' @param sigma
#' #'
#' #' @return real-valued covariance between squared modulus SWDT coefficients
#' #'
#' cov_swdft_sqmod <- function(k, l, delta, n, sigma) {
#'   case1 <- 3 * (n - delta)
#'   case2 <- (n^2 - (n - delta))
#'
#'   case3 <- complex(length.out=1, real=0, imaginary=0)
#'   case4 <- complex(length.out=1, real=0, imaginary=0)
#'   srange <- delta:(n-1)
#'   vrange <- delta:(n-1)
#'   for (s in srange) {
#'     for (v in vrange) {
#'       indicator <- ifelse(s==v, 0, 1)
#'       case3_pow <- (k + l) * (v - s)
#'       case4_pow <- (k * (v - s)) + (l * (s - v))
#'       case3 <- case3 + (indicator * prou(n=n)^( -case3_pow ))
#'       case4 <- case4 + (indicator * prou(n=n)^( -case4_pow ))
#'     }
#'   }
#'
#'   part1 <- (sigma^4 / n^2) * (case1 + case2 + Re(case3) + Re(case4))
#'   return(part1 - sigma^4)
#' }
