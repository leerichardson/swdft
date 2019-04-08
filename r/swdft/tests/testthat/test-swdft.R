context("1D SWDFT")

test_that("1D SWDFT Algorithms give the same answer", {
  N <- 100
  window_size <- 2^4
  x <- complex(N, rnorm(N), rnorm(N))

  a_ref <- swdft::swdft(x=x, n=window_size, type="fft")
  a_fftw <- swdft::swdft(x=x, n=window_size, type="fftw")

  expect_true(all(round( Re(a_ref$a - a_fftw$a) , digits = 5) == 0))
  expect_true(all(round( Im(a_ref$a - a_fftw$a) , digits = 5) == 0))
})
