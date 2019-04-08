#' Create taper for the SWDFT
#'
#' @param n window size
#' @param taper taper type. Can be either 'none' (default) or 'cosine'
#' @param p proportion to taper on each end, if cosine taper is used
#'
#' @param A length n taper
#'
get_taper <- function(n, taper, p) {
  if (taper == 'none') {
    w <- rep(1, n)
  } else if (taper == 'cosine') {
    w <- cosine_taper(n=n)
  } else {
    stop("taper must be 'none' or 'cosine'")
  }

  return(w)
}

#' Cosine bell data taper
#'
#' @param n length of time-series to taper
#' @param p proportion of ends to taper
#'
#' @return length n cosine bell taper w/ proportion p
#'
cosine_taper <- function(n, p=.1) {
  ## Verify the proportion of series to taper is chosen correctly
  if (p < 0 | p > .5) { stop("'p' must  be between 0 and .5") }

  ## Construct the cosine taper
  m <- floor(n * p)
  w <- 0.5 * (1 - cos(pi * seq.int(1, 2 * m - 1, by = 2)/ (2 * m)) )
  taper <- c(w, rep_len(1, n - 2 * m), rev(w))

  return(taper)
}
