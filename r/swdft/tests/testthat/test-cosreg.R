context("Cosine and local cosine regression")

test_that("Cosine regression works", {
  ## Generate a noiseless cosine signal
  N <- 15
  window_size <- 2^3
  A <- 1
  Fr <- 2 / window_size
  phase <- 1
  signal <- swdft::cosine(N=N, A=A, Fr=Fr, phase=phase)

  ## Fit the cosine regression at the true frequency
  cosreg_fit <- swdft::cosreg(x=signal, f=Fr)

  ## Check that the parameters match
  expect_true( all(class(cosreg_fit) == c("swdft_cosreg", "swdft_mod")) )
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
  local_cosreg_fit <- swdft::local_cosreg(x=x)

  ## Verify we get the right fit in the noiseless case
  expect_true(all(class(local_cosreg_fit) == c("swdft_local_cosreg", "swdft_mod")))
  expect_true(all(round(residuals(local_cosreg_fit), digits=3) == 0))
})
