#' MLE Estimate of Local Signal
#'
#' @param x 1D Vector
#' @param n window size
#'
fit_burst <- function(x, n) {
  # Extract the initial parameters
  N <- length(x)

  # Pad the time-series  w/ 0's for east time-domain translation
  xpad <- c(rep(0, n-1), x)

  # Compute the Spectrogram ---
  a <- swdft::swdft(x=xpad, n=n)
  num_freqs <- floor( (n - 1) / 2)
  freq_range <- (1:num_freqs) + 1
  amod <- Mod(a[freq_range, ])^2

  # Estimate starting values and numerical range for S, L, and F ---
  max_amod_inds <- which(amod==max(amod), arr.ind=TRUE)

  ## Compute the range of F to numerically search
  khat <- max_amod_inds[1, 1]
  Fmin <- (khat - 1) * (N / n); cat("Fmin: ", Fmin, " \n")
  Fmax <- (khat + 1) * (N / n); cat("Fmax: ", Fmax, " \n")

  ## Estimate the starting locations (S and L) of the periodic signal
  phat <- max_amod_inds[1, 2]
  mod_thresh <- as.numeric(quantile(x=c(amod), probs = c(.9)))
  seq_above_thresh <- get_seq_above(amod[khat, ], phat, mod_thresh)

  S_range <- max(seq_above_thresh[1] - n, 0):seq_above_thresh[1]
  num_above_thresh <- (seq_above_thresh[2] - seq_above_thresh[1])
  L_range <- (num_above_thresh - n):num_above_thresh

  # Maximize the Log-Likelood w/ a numerical search ---
  SL_grid <- expand.grid(S=S_range, L=L_range)
  SL_grid$Fhat <- NA
  SL_grid$Ahat <- NA
  SL_grid$Phihat <- NA
  SL_grid$Sigma_hat <- NA
  SL_grid$Loglik <- NA

  for (i in 1:nrow(SL_grid)) {
    max_sl <- max_loglik_sl(x=x, S=SL_grid[i, 1], L=SL_grid[i, 2], Fmin=Fmin, Fmax=Fmax, n=n, N=N)
    SL_grid[i, "Fhat"] <- max_sl$Fhat
    SL_grid[i, "Ahat"] <- max_sl$Ahat
    SL_grid[i, "Phihat"] <- max_sl$Phihat
    SL_grid[i, "Sigma_hat"] <- max_sl$Sigma_hat
    SL_grid[i, "Loglik"] <- max_sl$Loglik
  }

  return(SL_grid)
}

#' Compute the log likelihood
#'
#' @param x
#' @param S
#' @param L
#' @param Fmin
#' @param Fmax
#' @param n
#' @param N
#'
max_loglik_sl <- function(x, S, L, Fmin, Fmax, n, N) {
  ## Compute the maximum log likelihood over this interval of frequencies
  max_loglik_F <- stats::optimize(f=compute_burst_loglik, interval=c(Fmin, Fmax), maximum=TRUE,
                                  x=x, S=S, L=L, n=n, N=N)

  ## Recover the exact parameters at the maximum
  analytic_maxF <- get_analytic_aps(x=x, Fr=max_loglik_F$maximum, S=S, L=L, n=n, N=N)

  return(list(Fhat=max_loglik_F$maximum,
              Ahat=analytic_maxF$A_hat,
              Phihat=analytic_maxF$Phi_hat,
              Sigma_hat=analytic_maxF$sigma_hat,
              Loglik=max_loglik_F$objective))
}

#' Compute the log likeligood given S, L, and Ft
#'
#' @param Fr
#' @param x
#' @param S
#' @param L
#' @param n
#' @param N
#'
compute_burst_loglik <- function(Fr, x, S, L, n, N) {
  ## Compute the analytic solutions for A, phi, and sigma
  analytic_params <- get_analytic_aps(x=x, Fr=Fr, S=S, L=L, n=n, N=N)

  ## Use analytic solutions to compute the log likelihood value
  loglik <- burst_loglik(x=x, Fr=Fr, S=S, L=L, N=N, A=analytic_params$A_hat,
                         Phi=analytic_params$Phi_hat, sigma=analytic_params$sigma_hat,
                         fitted=analytic_params$fitted)

  ## Output: numeric log likelihood
  return(loglik)
}

#' Analytic SOlutipn for A, phi, and sigma given S, L, and F
#'
#' @param x
#' @param Fr
#' @param S
#' @param L
#' @param n
#' @param N
#'
get_analytic_aps <- function(x, Fr, S, L, n, N) {
  # Implement the MLE estimation of this signal w/ known (S, L, and F)
  t_range <- S:(S+L-1)

  ## Get Ahat, Phihat with least squares
  indicator <- rep(0, N)
  indicator[(S + 1):(S + L)] <- 1

  Cos <- cos_vec(Freq=Fr, N=N, S=S, L=L) * indicator
  Sin <- sin_vec(Freq=Fr, N=N, S=S, L=L) * indicator
  Hat <- cbind(Cos, Sin)
  betas <- solve(t(Hat) %*% Hat) %*% t(Hat) %*% x
  A_hat <- sqrt(betas[1]^2 + betas[2]^2)
  Phi_hat <- atan(-betas[2] / betas[1])

  ## Compute sigma hat
  fitted <- swdft::local_signal(N=N, A=A_hat, Fr=Fr, phase=Phi_hat, S=S, L=L)
  sigma_hat <- sqrt((1 / N) * sum( (x - fitted)^2 ))

  ## Output A, Phi, and Sigma
  return(list(A_hat=A_hat, Phi_hat=Phi_hat, sigma_hat=sigma_hat, fitted=fitted))
}

#' Calculate log likelihoodfor a set of parameters
#'
#' @param x
#' @param Fr
#' @param S
#' @param L
#' @param A
#' @param Phi
#' @param sigma
#' @param fitted
#'
burst_loglik <- function(x, Fr, S, L, A, Phi, sigma, N, fitted) {
  sum_val <- (-1 / (2 * sigma^2)) * sum( (x - fitted)^2 )
  -N * log(sigma) + sum_val
}

#' Extract the sequence around the maximum above a threshold
#'
#' @param x sequence from the maximum frequency
#' @param phat maximum index
#' @param thresh threshold
#'
get_seq_above <- function(a, phat, thresh) {
  under_thresh_before <- phat - which(a[(phat-1):1] < thresh)[1] + 1
  under_thresh_after <- phat + which(a[(phat + 1):length(a)] < thresh)[1] - 1

  stopifnot( all( a[under_thresh_before:under_thresh_after] ) )

  return(c(under_thresh_before, under_thresh_after))
}

#' Cosine basis vector for MLE
#'
#' @param Freq
#' @param N
#' @param S,
#' @param L
#'
cos_vec <- function(Freq, N, S, L) {
  cos(x = ( 2 * pi * Freq * (0:(N-1)) ) / N)
}

#' Sin basis vector for MLE
#'
#' @param Freq
#' @param N
#' @param S,
#' @param L
#'
sin_vec <- function(Freq, N, S, L) {
  sin(x = ( 2 * pi * Freq * (0:(N-1)) ) / N)
}
