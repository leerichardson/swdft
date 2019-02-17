#' Local Cosine Regression
#'
#' @param x 1D Signal
#' @param ptype
#' @param ktype
#'
local_cosreg <- function(x, ptype="around_max", ktype="grid") {
  cat("Local Cosine Regression \n")
}

#' Refined optimization of Frequency and Window Position for fixed window size
#'
#' @param x
#' @param n
#' @param phat
#' @param ptype
#' @param prange
#' @param kgrid_len
#'
optimize_kp <- function(x, a, phat, khat, pwidth=5, ptype="around_max",
                        krange=NULL, ktype="grid", kgrid_len=200) {
  N <- ncol(a)
  n <- nrow(a)
  t <- 0:(n-1)
  browser()
  ## Extract the range of window positions to search
  prange <- get_p_range(phat=phat, n=n, N=N, pwidth=pwidth, type=ptype)

  ## Optionally specify the grid to search of that's the way we're selecting k
  if (ktype == "grid") {
    kpmat <- matrix(data=NA, nrow=length(prange)* kgrid_len, ncol=3)
  }

  iter <- 0
  for (p in prange) {
    iter <- iter + 1
    iter_inds <- (((iter - 1) * kgrid_len) + 1):(iter * kgrid_len)
    if (p-n+1 < 1) { stop("FOUND AN INELIGIBLE WINDOW SIZE COMBINATION!") }
    xwin <- x[(p-n+1):p]

    ## If we're using a grid search,
    if (ktype == "grid") {
      if (is.null(krange)) { krange <- get_k_range() }
      kwin_grid <- seq(from=krange[1], to=krange[2], length=kgrid_len)
      maxval <- sapply(X=kwin_grid * n, FUN=eval_kp, x=xwin, n=n, t=t)
    }

    kpmat[iter_inds, 1] <- kwin_grid
    kpmat[iter_inds, 2] <- p
    kpmat[iter_inds, 3] <- maxval
  }

  return(kpmat)
}

get_k_range <- function(khat, n) {
  browser()
}

#' Evaluate the SWDFT periodogram function
#'
#' @param
#'
eval_kp <- function(k, x, n, t) {
  twiddle <- swdft::prou(n=n)^(-k * t)
  stopifnot( length(twiddle) == length(x))
  return( sqrt(2 / n) * Mod(sum(x * twiddle))^2 )
}

#' Get range of P's to search
#'
#' @param phat
#' @param n
#' @param N
#' @param prange
#' @param type
#'
get_p_range <- function(phat, n, N, pwidth, type="around_max") { browser();
  fullp <- 1:N
  minp <- min( fullp[(fullp - n + 1) > 0] )
  maxp <- max( fullp )

  if (type == "around_max") {
    prange <- max(minp, phat-pwidth):min(maxp, (phat+pwidth))
  } else if (type == "fullp") {
    prange <- minp:maxp
  }

  return(prange)
}
