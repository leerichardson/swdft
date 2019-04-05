#' Local cosine regression
#'
#' @param x
#' @param lmin
#' @param pwidth
#' @param kwidth
#' @param verbose
#'
#' @return S3 object of class 'swdft_local_cosreg'
#'
local_cosreg <- function(x, lmin=6, pwidth=5, kwidth=1, verbose=FALSE) {
  N <- length(x)

  ## Set up matrix to store parameters for each window size
  grid_names <- c("n", "f", "p", "S", "L", "A", "Phi", "sigma", "loglik")
  param_grid <- matrix(data=NA_real_, ncol=length(grid_names), nrow=length(lmin:N) * ((2*pwidth) + 1) )
  colnames(param_grid) <- grid_names

  ## Extract the parameters for each possible window size
  iter <- 0
  for (n in lmin:N) {
    if (verbose) { cat("Window Size: ", n, " \n") }

    ## Compute the range of frequencies and window positions to search
    a <- swdft::swdft(x=x, n=n)
    freq_range <- swdft::get_freq_range(a=a, kwidth=kwidth)
    phat <- which(Mod(a)^2 == max(Mod(a)^2), arr.ind=TRUE)[1,2]
    prange <- swdft::get_p_range(phat=phat, n=n, N=N, pwidth=pwidth)

    ## Optimize the frequuency selection for each window position
    for (p in prange) {
      iter <- iter + 1
      SL <- get_sl(n=n, p=p-1)

      ## Search for the global maxima of the frequency
      maxfreq <- nloptr::nloptr(x0=mean(freq_range), eval_f=lcr_loglik, lb=freq_range[1], ub=freq_range[2],
                                 x=x, S=SL[1], L=SL[2], ftype="negoptim",
                                 opts=list("algorithm"="NLOPT_GN_DIRECT_L", "maxeval"=100))
      loglik <- lcr_loglik(f=maxfreq$solution, x=x, S=SL[1], L=SL[2], ftype="full")

      ## Store parameter estimates for this window size
      param_grid[iter, ] <- c(n, maxfreq$solution, p, SL[1], SL[2], loglik[4:7])
    }
  }

  ## Return an S3 object w/ results
  max_ind <- which.max(param_grid[, "loglik"])
  fitted <- swdft::local_signal(N=N, A=param_grid[max_ind,"A"], Fr=param_grid[max_ind,"f"],
                                phase=param_grid[max_ind,"Phi"], S=param_grid[max_ind,"S"],
                                L=param_grid[max_ind, "L"])

  local_cosreg_obj <- structure(list(coefficients=param_grid[max_ind, c("f", "S", "L", "A", "Phi", "sigma")],
                                     fitted=fitted,
                                     residuals=x-fitted,
                                     data=x,
                                     window_params=param_grid[complete.cases(param_grid), ]),
                                class=c("swdft_local_cosreg", "swdft_cosreg"))
  return(local_cosreg_obj)
}

#' Get range of frequencies to search
#'
#' @param a swdft to search
get_freq_range <- function(a, kwidth) {
  khat <- which(Mod(a)^2 == max(Mod(a)^2), arr.ind=TRUE)[1,1] - 1
  fmin <- max(0, (khat - kwidth) / nrow(a))
  fmax <- min((khat + kwidth) / nrow(a), .5)

  return( c(fmin, fmax) )
}

#' Get range of P's to search
#'
#' @param phat
#' @param n
#' @param N
#' @param prange
#' @param type
#'
get_p_range <- function(phat, n, N, pwidth, type="around_max") {
  fullp <- 1:N
  minp <- min( fullp[(fullp - n + 1) > 0] )
  maxp <- max( fullp )

  if (type == "around_max") {
    phat_min <- max(phat-pwidth, minp)
    phat_max <- ifelse(test=phat+pwidth>minp, yes=min(maxp, phat+pwidth), no=maxp)
    prange <- phat_min:phat_max

  } else if (type == "fullp") {
    prange <- minp:maxp
  }

  return(prange)
}

#' Extract signal parameters
#'
#' @param n
#' @param p
#'
get_sl <- function(n, p) {
  L <- n
  S <- p - n + 1
  return( c(S, L) )
}

#' Log Likelihood
#'
#' @param f
#' @param x
#' @param S
#' @param L
#' @param ftype
#'
lcr_loglik <- function(f, x, S, L, ftype="full") {
  A_Phi <- swdft::get_aphi(x=x, S=S, L=L, f=f)
  fitted <- swdft::local_signal(N=length(x), A=A_Phi[1], Fr=f, phase=A_Phi[2], S=S, L=L)
  sigma <- swdft::get_sigma(x=x, fitted=fitted, N=length(x))
  loglik <- swdft::get_loglik(x=x, fitted=fitted, sigma=sigma, N=length(x))

  ## Optionally return either all the parameters or just the log likelihood
  ## in case we are optimizing the function
  if (ftype == "full") {
    return(c(S, L, f, A_Phi[1], A_Phi[2], sigma, loglik))
  } else if (ftype == "optim") {
    return(loglik)
  } else if (ftype == "negoptim") {
    return(-loglik)
  }
}

#' Extract amplitude and phase
#'
#' @param x
#' @param f
#' @param S,
#' @param L
#'
get_aphi <- function(x, S, L, f) {
  N <- length(x)
  indicator <- rep(0, N)
  indicator[(S+1):(S+L)] <- 1
  U <- matrix(data=NA, nrow=N, ncol=2)
  U[, 1] <- swdft::cosine(N=N, Fr=f) * indicator
  U[, 2] <- swdft::sine(N=N, Fr=f) * indicator
  fitted <- lm(x ~ U - 1)
  beta <- coefficients(fitted)
  A <- sqrt( sum(beta^2) )
  Phi <- atan2(y=-beta[2], x=beta[1])

  return( c(A, Phi) )
}

#' Extract estimator of sigma
#'
#' @param x
#' @param fitted
#' @param N
#'
get_sigma <- function(x, fitted, N) {
  sqrt((1 / N) * sum( (x - fitted)^2 ))
}

#' Compute the log likelihood
#'
#' @param x
#' @param fitted
#' @param sigma
#' @param N
#'
get_loglik <- function(x, fitted, sigma, N) {
  sum_val <- (-1 / (2 * sigma^2)) * sum( (x - fitted)^2 )
  return( -N * log(sigma) + sum_val )
}

#' Evaluate the localized periodogram function
#'
#' @param f
#' @param x
#' @param n
#' @param t
#' @param normalize
#' @param ftype
#'
eval_swdft <- function(f, x, n, t, normalize=sqrt(2/n), ftype="optim") {
  twiddle <- swdft::prou(n=n)^(-f * t )
  if (length(x) != length(twiddle)) { browser() }

  val <- normalize * Mod( sum(x * twiddle) )^2

  if (ftype=="optim") {
    return(val)
  } else if (ftype=="negoptim") {
    return(-val)
  }
}
