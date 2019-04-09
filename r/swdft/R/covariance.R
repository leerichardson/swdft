#' Covariance between two complex-numbered outputs
#'
#' @param k frequency of first coefficient
#' @param l frequency of second coefficient
#' @param delta window position shift of second coefficient
#' @param n window size
#' @param sigma white noise standard error
#'
#' @return complex-valued number of the covariance
#'
cov_swdft_cnum <- function(k, l, delta, n, sigma) {
  shift <- dirichlet(x=(2*pi*(l-k))/n, phase=(-2*pi*delta*l)/n, a=delta, b=n-1)
  return( (sigma^2 / n) * shift)
}

cov_swdft_rr <- function(k, l, delta, n, sigma) {
  shift1 <- Re( dirichlet(x=(-2*pi*(k-l))/n, phase=(-2*pi*delta*l)/n, a=delta, b=n-1) )
  shift2 <- Re( dirichlet(x=(-2*pi*(k+l))/n, phase=(2*pi*delta*l)/n, a=delta, b=n-1) )
  return( (sigma^2 / (2*n)) * (shift1 + shift2) )
}

cov_swdft_ii <- function(k, l, delta, n, sigma) {
  shift1 <- Re( dirichlet(x=(-2*pi*(k-l))/n, phase=(-2*pi*delta*l)/n, a=delta, b=n-1) )
  shift2 <- Re( dirichlet(x=(-2*pi*(k+l))/n, phase=(2*pi*delta*l)/n, a=delta, b=n-1) )
  return( (sigma^2 / (2*n)) * (shift1 - shift2) )
}

cov_swdft_ri <- function(k, l, delta, n, sigma) {
  shift1 <- Im( dirichlet(x=(-2*pi*(k+l))/n, phase=(2*pi*delta*l)/n, a=delta, b=n-1) )
  shift2 <- Im( dirichlet(x=(-2*pi*(k-l))/n, phase=(-2*pi*delta*l)/n, a=delta, b=n-1) )
  return( (sigma^2 / (2*n)) * (shift1 - shift2) )
}

cov_swdft_ir <- function(k, l, delta, n, sigma) {
  shift1 <- Im( dirichlet(x=(-2*pi*(k+l))/n, phase=(2*pi*delta*l)/n, a=delta, b=n-1) )
  shift2 <- Im( dirichlet(x=(-2*pi*(k-l))/n, phase=(-2*pi*delta*l)/n, a=delta, b=n-1) )
  return( (sigma^2 / (2*n)) * (shift1 + shift2) )
}
