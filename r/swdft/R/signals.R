#' Local Periodic Signal
#'
#' @param N signal length
#' @param A Amplitude
#' @param Fr Frequency: Number of cycles in a length N period
#' @param phase phase
#' @param S start of local signal
#' @param L length of local signal
#'
#' @export
#'
#' @return length N local periodic signal
#'
local_signal <- function(N, A=1, Fr=1, phase=0, S=0, L=N) {
  periodic_signal <- cosine(N, A=A, Fr=Fr, phase=phase)

  ## Indicator, adjusting for the fact the math notation uses a 0-index
  indicator <- rep(0, N)
  indicator[(S + 1):(S + L)] <- 1

  periodic_signal * indicator
}

#' The principal nth root of unity
#'
#' @param n integer root
#'
#' @export
#'
#' @return complex number
#'
prou <- function(n) {
  i <- complex(length.out = 1, real = 0, imaginary = 1)
  exp((2 * pi * i) / n)
}

#' Cosine signal with adjustable parameters
#'
#' @param N signal length
#' @param A Amplitude
#' @param Fr Frequency: Number of cycles in a length N period
#' @param phase phase
#'
#' @export
#'
#' @return numeric vector with cosine function of x
#'
cosine <- function(N, A=1, Fr=1, phase=0) {
  A * cos( (2 * pi * Fr * (0:(N - 1)) ) + phase )
}

#' Sine signal with adjustable parameters
#'
#' @param N length signal
#' @param A Amplitude
#' @param Fr Frequency: Number of cycles in a length N period
#' @param phase phase
#'
#' @export
#'
#' @return numeric vector with sine
#'
sine <- function(N, A=1, Fr=1, phase=0) {
  A * sin( (2 * pi * Fr * (0:(N - 1)) ) + phase )
}
