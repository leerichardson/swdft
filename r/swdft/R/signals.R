#' The principal nth root of unity
#'
#' @param n integer root
#'
#' @return the resulting complex number
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
#' @return numeric vector with cosine function of x
#'
cosine <- function(N, A=1, Fr=1, phase=0) {
  x <- 0:(N - 1)
  A * cos( ((2 * pi * Fr * x) / N) + phase )
}

#' Sine signal with adjustable parameters
#'
#' @param N length signal
#' @param A Amplitude
#' @param Fr Frequency: Number of cycles in a length N period
#' @param phase phase
#'
#' @return numeric vector with sine
#'
sine <- function(N, A=1, Fr=1, phase=0) {
  x <- 0:(N - 1)
  A * sine( ((2 * pi * Fr * x) / N) + phase )
}

#' Local Periodic Signal
#'
#' @param N signal length
#' @param A Amplitude
#' @param Fr Frequency: Number of cycles in a length N period
#' @param phase phase
#' @param S start of local signal
#' @param L length of local signal
#'
#' @return length N local periodic signal
#'
local_signal <- function(N, A=1, Fr=1, phase=0, S=0, L=N) {
  periodic_signal <- cosine(N, A=A, Fr=Fr, phase=phase)

  indicator <- rep(0, N)
  indicator[(S + 1):(S + L)] <- 1

  periodic_signal * indicator
}
