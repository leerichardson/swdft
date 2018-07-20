#' Estimate parameters of a local cosine signal
#'
#' @param b 2D complex array from the 'swdft' function
#' @param lmin the smallest local-signal to search for
#' @param full_estimation optionally return all parameter
#' estimates for every grid position, defaults to false.
#'
#' @return A list with three components
#' \item{ls_params}{A named numeric vector with the 5 least squars parameter
#' estimates, mean squared error, mean squared error of the linear fit, mean squared
#' error of fitting the mean, and the change in mse between the two}
#' \item{fitted}{length N time-series with the fitted values corresponding to
#' our least squares fit }
#' \item{full_estimation}{if full_estimation is TRUE, then return all parameter estimates
#' corresponding to each grid position}
#'
#' @export
#'
#' @examples
#' # First generate a local periodic signal
#' N <- 40
#' window_size <- 2^4
#' A <- 1
#' Fr <- 4
#' f <- (Fr * window_size) / N
#' phase <- 1
#' S <- 10
#' L <- 10
#' x <- local_signal(N=N, A=A, Fr=Fr, phase=phase, S=S, L=L)
#'
#' # Take the SWDFT and estimate parameters with fit_local_cosine
#' b <- swdft(x=x, n=window_size, normalize=(1/sqrt(window_size)))
#' bhat <- fit_local_cosine(b)
#'
fit_local_cosine <- function(b, lmin=8, full_estimation=FALSE) {
  P <- ncol(b)
  n <- nrow(b)
  N <- P + n - 1
  p_range <- (n - 1):(N - 1)
  i <- complex(length.out = 1, real = 0, imaginary = 1)

  # Select which frequency time-series to search
  k <- coef_to_search(b=b, n=n)
  bk <- b[k + 1, ]

  # Estimate the parameters based on a grig search
  grid_mat <- get_grid(N, n, lmin)
  params <- t(apply(X = grid_mat, MARGIN = 1, FUN = grid_search,
                    bk=bk, N=N, n=n, i=i, k=k, p_range=p_range))
  colnames(params) <- c("S", "L", "f", "A", "phase", "mse_b", "mse_a", "mse_c")
  ls_params <- params[which.min(params[, "mse_b"]), ]

  # Calculate the fitted values of the time-series using the SWDFT
  fitted <- swdft::local_signal(N=N,
                                A=ls_params["A"],
                                Fr=(ls_params["f"] * N) / n,
                                phase=ls_params["phase"],
                                S=ls_params["S"],
                                L=ls_params["L"])

  # Return least squares parameters, fitted values, and the full estimation results
  if (!full_estimation) { params <- NULL }
  return(list(ls_params=c(ls_params, k=as.numeric(k)), fitted=fitted, full_estimation=params))
}

#' Compute possible start/length values of local signal
#'
#' @param N Length of original signal
#' @param n Window Size of SWDFT
#' @param lmin Minimum signal to search for
#'
get_grid <- function(N, n, lmin) {
  num_basis <- sum((N - lmin + 1):1)
  grid_mat <- matrix(data=NA, nrow=num_basis, ncol=2)

  count <- 0
  for (S in 0:(N - lmin)) {
    for (L in (lmin:(N - S))) {
      count <- count + 1
      grid_mat[count, ] <- c(S, L)
    }
  }

  return(grid_mat)
}

#' Choose frequency time-series in SWDFT to search
#'
#' @param b 2D complex array from the 'swdft' function
#' @param n window size
#'
#' @return k: the frequency time-series to search
#'
coef_to_search <- function(b, n) {
  m <- ifelse(n %% 2 == 0, floor(n / 2) - 1, floor(n / 2))
  bmod_elig <- Mod(b[2:(m + 1), ])^2

  # Take the first maximum value (Note: Could do better here, but works for now)
  largest_k <- which(bmod_elig == max(bmod_elig), arr.ind = TRUE)[1, 1]
  return(largest_k)
}

#' Estimate Parameters for all Start/Length Combinations of Local Signal
#'
#' @param SL Vector of S and L paramers
#' @param bk frequency time-series to search
#' @param N lenth of original time-series
#' @param n window size
#' @param p_range range of window positions
#' @param k fourier frequency we are searching
#' @param i imaginary number
#'
grid_search <- function(SL, bk, N, n, p_range, k, i) {
  S <- SL[1]
  L <- SL[2]

  optimal_freq <- stats::optimize(f=compute_mse, interval=c(k - .5, k + .5),
                                   S=S, L=L, bk=bk, n=n, N=N,
                                   p_range=p_range, i=i, k=k)
  optimal_params <- compute_mse(f=optimal_freq$minimum,
                                S=S, L=L, bk=bk, n=n, N=N,
                                p_range=p_range, i=i, k=k,
                                return_amp_phase=TRUE)

  return(c(S, L, optimal_freq$minimum, optimal_params[4], optimal_params[5],
           optimal_params[1], optimal_params[2], optimal_params[3]))
}

#' Compute MSE of local fit given frequency F and grid-point
#'
#' @param f frequency: number of cycles in length n window
#' @param S start of local signal
#' @param L length of local signal
#' @param bk frequency time-series to search
#' @param n window size
#' @param N length of original time-series
#' @param p_range range of window positions
#' @param k fourier frequency we are searching
#' @param i imaginary number
#' @param return_amp_phase should we return  *just* the mse, or the mse plus
#' the parameters estimates. Used because we need a single number output
#' for optimization functions.
#'
compute_mse <- function(f, S, L, bk, n, N, p_range, k, i, return_amp_phase=FALSE) {
  Fr <- (f * N) / n
  C1 <- sapply(X = p_range, FUN = c1kp, n=n, N=N, f=Fr, k=k, S=S, L=L, i=i)
  C2 <- sapply(X = p_range, FUN = c2kp, n=n, N=N, f=Fr, k=k, S=S, L=L, i=i)

  # Fit the linearized model to the real-part of the signal
  X_re <- cbind(Re(C1), Re(C2))
  Y_re <- Re(bk)
  fit_re <- stats::.lm.fit(x = X_re, y = Y_re)
  beta1_re <- stats::coefficients(fit_re)[1]
  beta2_re <- stats::coefficients(fit_re)[2]
  fitted_re <- (beta1_re * Re(C1)) + (beta2_re * Re(C2))
  mse_fitted_re <- sum( ( Re(bk) - fitted_re )^2 )

  # If we are optimizing the function, only return delta MSE
  if (return_amp_phase == FALSE) { return(mse_fitted_re) }

  # Optionally return parameter estimates for phase and amplitude
  amp_re <- sqrt(beta1_re^2 + beta2_re^2)
  phase_re <- atan2(y = -beta2_re, x = beta1_re)
  if (phase_re < 0){ phase_re <- pi + (pi - abs(phase_re)) }

  mse_xbar_re <- sum(  ( Re(bk) - mean(Re(bk)) )^2 )
  delta_mse_re <- mse_fitted_re - mse_xbar_re

  return(c(mse_fitted_re, mse_xbar_re, delta_mse_re, amp_re, phase_re))
}

#' Basis C_{k, 1, p} in fitting procedure
#'
#' @param p window position
#' @param n window size
#' @param N lenth of original time-series
#' @param f frequency: number of cycles in length n window
#' @param k fourier frequency we are searching
#' @param S start of local signal
#' @param L length of local signal
#' @param i imaginary number
#'
c1kp <- function(p, n, N, f, k, S, L, i) {
  j <- 0:(n - 1)
  phat <- p - n + 1 + j
  inds <- phat >= S & phat <= S + L - 1

  cos_part <- cos((2 * pi * phat[inds] * f) / N)
  twiddle <- prou(n = n)^(-j[inds] * k)
  twiddle <- exp((-2 * pi * i * j[inds] * k) / n)

  (1 / sqrt(n)) * sum(twiddle * cos_part)
}

#' Basis C_{k, 1, p} in fitting procedure
#'
#' @param p window position
#' @param n window size
#' @param N lenth of original time-series
#' @param f frequency: number of cycles in length n window
#' @param k fourier frequency we are searching
#' @param S start of local signal
#' @param L length of local signal
#' @param i imaginary number
#'
c2kp <- function(p, n, N, f, k, S, L, i) {
  j <- 0:(n - 1)
  phat <- p - n + 1 + j
  inds <- phat >= S & phat <= S + L - 1

  sin_part <- sin((2 * pi * phat[inds] * f) / N)
  twiddle <- exp((-2 * pi * i * j[inds] * k) / n)

  (1 / sqrt(n)) * sum(twiddle * sin_part)
}
