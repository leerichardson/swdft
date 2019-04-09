#' Local cosine regression
#'
#' @param x numeric signal to apply local cosine regression on
#' @param lmin integer. minimum signal length (L parameter) to search
#' @param pwidth integer. the range of window positions to search for each window size
#' @param kwidth integer. the width of frequencies to search
#' @param verbose logical. whether or not to print intermediate results
#'
#' @export
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
    a <- swdft(x=x, n=n)$a
    freq_range <- get_freq_range(a=a, kwidth=kwidth)
    phat <- which(Mod(a)^2 == max(Mod(a)^2), arr.ind=TRUE)[1,2]
    prange <- get_p_range(phat=phat, n=n, N=N, pwidth=pwidth)

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
  fitted <- local_signal(N=N, A=param_grid[max_ind,"A"], Fr=param_grid[max_ind,"f"],
                                phase=param_grid[max_ind,"Phi"], S=param_grid[max_ind,"S"],
                                L=param_grid[max_ind, "L"])

  local_cosreg_obj <- new_swdft_local_cosreg(coefficients=param_grid[max_ind, c("f", "S", "L", "A", "Phi", "sigma")],
                                                    fitted=fitted, residuals=x-fitted, data=x,
                                                    window_params=param_grid[stats::complete.cases(param_grid),])

  return(local_cosreg_obj)
}

#' Constructor function for class 'swdft_local_cosreg'
#'
#' @param coefficients matrix of coefficients for cosine regression model
#' @param fitted fitted values of cosine regression model
#' @param residuals residuals of cosine regression model
#' @param data original signal used to fit cosine regression
#' @param window_params data frame of fitted coefficients for each window size
#'
#' @return list with the following elements
#' \itemize{
#'   \item coefficients. A matrix of parameters, the three columns are: 1. amplitude 2. phase, and 3. frequency.
#'   There is only more that one row used when multiple frequencies are fit sequentially.
#'   \item fitted. fitted values of cosine regression model
#'   \item residuals. residuals of cosine regression model
#'   \item data. original signal used to fit cosine regression
#'   \item window_params. data frame of fitted coefficients for each window size
#' }
#'
new_swdft_local_cosreg <- function(coefficients, fitted, residuals, data, window_params) {
  structure(list(coefficients=coefficients,
                 fitted=fitted,
                 residuals=residuals,
                 data=data,
                 window_params=window_params),
            class=c("swdft_local_cosreg", "swdft_mod"))
}

#' Get range of frequencies to search
#'
#' @param a 2D complex-valued array. The SWDFT to search
#' @param kwidth integer. the width of frequencies to search
#'
get_freq_range <- function(a, kwidth) {
  maxk <- floor(nrow(a) / 2) + 1

  khat <- which(Mod(a[1:maxk,])^2 == max(Mod(a[1:maxk,])^2), arr.ind=TRUE)[1,1] - 1
  fmin <- max(0, (khat - kwidth) / nrow(a))
  fmax <- min((khat + kwidth) / nrow(a), .5)

  if (fmax < fmin) { stop("fmax greater than fmin") }

  return( c(fmin, fmax) )
}

#' Get range of P's to search
#'
#' @param phat integer. Window position with largest SWDFT coefficient
#' @param n integer. window size
#' @param N integer. Signal length
#' @param pwidth integer. the range of window positions to search for each window size
#' @param type character. either 'around max' or 'fullp'.
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
#' @param n window size
#' @param p window position
#'
get_sl <- function(n, p) {
  L <- n
  S <- p - n + 1
  return( c(S, L) )
}

#' Log Likelihood
#'
#' @param f frequency
#' @param x signal
#' @param S start parameter
#' @param L length pe
#' @param ftype what to return
#'
lcr_loglik <- function(f, x, S, L, ftype="full") {
  A_Phi <- get_aphi(x=x, S=S, L=L, f=f)
  fitted <- local_signal(N=length(x), A=A_Phi[1], Fr=f, phase=A_Phi[2], S=S, L=L)
  sigma <- get_sigma(x=x, fitted=fitted, N=length(x))
  loglik <- get_loglik(x=x, fitted=fitted, sigma=sigma, N=length(x))

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
#' @inheritParams lcr_loglik
#'
get_aphi <- function(x, S, L, f) {
  N <- length(x)
  indicator <- rep(0, N)
  indicator[(S+1):(S+L)] <- 1
  U <- matrix(data=NA, nrow=N, ncol=2)
  U[, 1] <- cosine(N=N, Fr=f) * indicator
  U[, 2] <- sine(N=N, Fr=f) * indicator
  fitted <- stats::lm(x ~ U - 1)
  beta <- stats::coefficients(fitted)
  A <- sqrt( sum(beta^2) )
  Phi <- atan2(y=-beta[2], x=beta[1])

  return( c(A, Phi) )
}

#' Extract estimator of sigma
#'
#' @param x signal
#' @param fitted fitted values
#' @param N length of x
#'
get_sigma <- function(x, fitted, N) {
  sqrt((1 / N) * sum( (x - fitted)^2 ))
}

#' Compute the log likelihood
#'
#' @param sigma estimated standard deviation
#' @inheritParams get_sigma
#'
get_loglik <- function(x, fitted, sigma, N) {
  sum_val <- (-1 / (2 * sigma^2)) * sum( (x - fitted)^2 )
  return( -N * log(sigma) + sum_val )
}

#' #' Evaluate the localized periodogram function
#' #'
#' #' @param f
#' #' @param x
#' #' @param n
#' #' @param t
#' #' @param normalize
#' #' @param ftype
#' #'
#' eval_swdft <- function(f, x, n, t, normalize=sqrt(2/n), ftype="optim") {
#'   twiddle <- prou(n=n)^(-f * t )
#'   if (length(x) != length(twiddle)) { browser() }
#'
#'   val <- normalize * Mod( sum(x * twiddle) )^2
#'
#'   if (ftype=="optim") {
#'     return(val)
#'   } else if (ftype=="negoptim") {
#'     return(-val)
#'   }
#' }
