#' Local Cosine Regression
#'
#' @param x
#' @param slf_type either "grid" or "window"
#' @param ptype
#' @param ktype
#'
local_cosreg <- function(x, slf_type="grid", lmin=6, ptype="around_max", ftype="optim",
                         freq_grid_len=100, verbose=FALSE) {
  N <- length(x)

  ## Fit local cosine regression using a grid search of S and L
  if (slf_type == "grid") {
    ## Create the grid of all potential window positions
    SL_grid <- swdft::get_grid(N=length(x), lmin=lmin)

    ## Get the overall of frequenciees to consider
    a <- swdft::swdft(x=x, n=floor(length(x)/2), taper='none') *  (2 / floor(length(x)/2))
    freq_range <- swdft::get_freq_range(a=a)
    freq_grid <- seq(from=freq_range[1], to=freq_range[2], length=freq_grid_len)

    ## Compute the log likelihood at each window position
    for (i in 1:nrow(SL_grid)) {
      S <- SL_grid[i, 1]
      L <- SL_grid[i, 2]
      if (verbose == TRUE) { cat("S: ", SL_grid[i, 1], " L: ", SL_grid[i, 2], " \n") }

      if (ftype == "grid") {
        logliks <- sapply(X=freq_grid, FUN=lcr_loglik, x=x, S=S, L=L)
        SL_grid[i, ] <- logliks[, which.max(logliks[7,])]

        ## Find the global optimum within the specified frequency range
      } else if (ftype == "optim") {
        maxfreq <- nloptr::nloptr(x0=mean(freq_range), eval_f=lcr_loglik, lb=freq_range[1], ub=freq_range[2],
                                  x=x, S=S, L=L, ftype="negoptim",
                                  opts=list("algorithm"="NLOPT_GN_DIRECT_L", "maxeval"=50))
        SL_grid[i, ] <- lcr_loglik(f=maxfreq$solution, x=x, S=SL_grid[i, 1], L=SL_grid[i,2], ftype="full")
      }
    }

  ## Fit the model by optimizing over window size
  } else if (slf_type == "window") {
    grid_names <- c("n", "f", "p", "S", "L", "A", "Phi", "sigma", "loglik")
    SL_grid <- data.frame(matrix(data=NA, nrow=(length(lmin:N))*10, ncol=length(grid_names)))
    names(SL_grid) <- grid_names

    iter <- 0
    for (n in lmin:N) {
      if (verbose == TRUE) { cat("Window Size: ", n, " \n") }
      t <- 0:(n-1)
      a <- swdft::swdft(x=x, n=n)
      freq_range <- swdft::get_freq_range(a=a)

      phat <- which(Mod(a)^2 == max(Mod(a)^2), arr.ind=TRUE)[1,2]
      prange <- swdft::get_p_range(phat=phat, n=n, N=N)

      for (p in prange) {
        iter <- iter + 1
        xwin <- x[(p-n+1):p]

         if (ftype == "optim") {
           SL <- get_sl(n=n, p=p-1)

           ## Search for the global maxima of the frequency
           maxfreq <- nloptr::nloptr(x0=mean(freq_range), eval_f=lcr_loglik, lb=freq_range[1], ub=freq_range[2],
                                     x=x, S=SL[1], L=SL[2], ftype="negoptim",
                                     opts=list("algorithm"="NLOPT_GN_DIRECT_L", "maxeval"=100))
          # ## Secondary local optimization around the global optimia (reccomended by Steven Johnson in NLOpt library)
          # maxfreq_local <- nloptr::nloptr(x0=maxfreq$solution, eval_f=lcr_loglik, lb=maxfreq$solution-.01, ub=maxfreq$solution+.01,
          #                                 x=x, S=SL[1], L=SL[2], ftype="negoptim",
          #                                 opts=list("algorithm"="NLOPT_LN_COBYLA", maxval=50))
          ## Get final analytic estimates
          loglik <- lcr_loglik(f=maxfreq$solution, x=x, S=SL[1], L=SL[2], ftype="full")
          SL_grid[iter, "n"] <- n
          SL_grid[iter, "f"] <- maxfreq$solution
          SL_grid[iter, "p"] <- p
          SL_grid[iter, "S"] <- SL[1]
          SL_grid[iter, "L"] <- SL[2]
          SL_grid[iter, 6:9] <- loglik[4:7]
         }
      }
    }
  } else {
    stop("slf_type must be 'grid' or 'window'!")
  }

  return(SL_grid)
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

#' Get range of frequencies to search
#'
#' @param a swdft to search
get_freq_range <- function(a, freq_width=1) {
  khat <- which(Mod(a)^2 == max(Mod(a)^2), arr.ind=TRUE)[1,1] - 1
  fmin <- max(0, (khat - freq_width) / nrow(a))
  fmax <- min((khat + freq_width) / nrow(a), .5)

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
get_p_range <- function(phat, n, N, pwidth=5, type="around_max") {
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

#' Compute possible start/length values of local signal
#'
#' @param N Length of original signal
#' @param n Window Size of SWDFT
#' @param lmin Minimum signal to search for
get_grid <- function(N, lmin) {
  num_basis <- sum((N - lmin + 1):1)
  grid_df <- data.frame(matrix(data=NA, nrow=num_basis, ncol=7))
  names(grid_df) <- c("S", "L", "f", "A", "Phi", "sigma", "loglik")

  count <- 0
  for (S in 0:(N - lmin)) {
    for (L in (lmin:(N - S))) {
      count <- count + 1
      grid_df[count, 1:2] <- c(S, L)
    }
  }

  return(grid_df)
}
