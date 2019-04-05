context("Cosine and local cosine regression")

test_that("Cosine regression works", {
  ## Generate a noiseless cosine signal
  N <- 15
  window_size <- 2^3
  A <- 1
  Fr <- 2 / window_size
  phase <- 1
  S <- 2
  L <- 10
  signal <- swdft::cosine(N=N, A=A, Fr=Fr, phase=phase)

  ## Fit the cosine regression at the true frequency
  cosreg_fit <- swdft::cosreg(x=signal, f=Fr)

  ## Check that the parameters match
  expect_true( class(cosreg_fit) == "swdft_cosreg" )
  expect_true( round(coefficients(cosreg_fit)[1, "A"] - A, digits=5) == 0 )
  expect_true( round(coefficients(cosreg_fit)[1, "phi"] - phase, digits=5) == 0 )
  expect_true( round(coefficients(cosreg_fit)[1, "f"] - Fr, digits=5) == 0 )
})

test_that("Local cosine regression works for noiseless signal", {
  ## Generate a local-in-time periodic signal
  N <- 15
  window_size <- 2^3
  A <- 1
  Fr <- 2 / window_size
  phase <- 1
  S <- 2
  L <- 10
  x <- swdft::local_signal(N=N, A=A, Fr=Fr, phase=phase, S=S, L=L)

  ## Fit local cosine regression
  slgrid_fit_window <- swdft::local_cosreg(x=x, slf_type="window", ftype="optim", verbose=TRUE)
  maxind <- which.max(slgrid_fit_window$loglik)
  maxparams <- slgrid_fit_window[maxind, ]
  maxparams
  # ## Check that the parameters match for the noiseless case
  # expect_true( round( bhat$ls_params["A"] - A, digits = 1)  == 0)
  # expect_true( round( bhat$ls_params["S"] - S, digits = 1)  == 0)
  # expect_true( round( bhat$ls_params["L"] - L, digits = 1)  == 0)
  # expect_true( round( bhat$ls_params["f"] - f, digits = 1)  == 0)
  # expect_true( round( bhat$ls_params["phase"] - phase, digits = 1)  == 0)
})
