context("2D and 3D SWDFTs")

test_that("2D, 3D SWDFT Algorithms give the same answer", {
  # --- 2D SWDFT ----
  N0 <- 40
  N1 <- 40
  n0 <- 16
  n1 <- 16

  data <- complex(N0 * N1, rnorm(N0 * N1), rnorm(N0 * N1))
  x <- array(data = data, dim = c(N0, N1))

  a_base <- swdft2d(x=x, n0=n0, n1=n1, type="fft")
  a_fftw <- swdft2d(x=x, n0=n0, n1=n1, type="fftw")

  expect_true( all( round( Re(a_base$a - a_fftw$a), digits = 5) == 0) )
  expect_true( all( round( Im(a_base$a - a_fftw$a), digits = 5) == 0) )

  # --- 3D SWDFT ---
  N0 <- 20
  N1 <- 20
  N2 <- 20
  n0 <- 8
  n1 <- 8
  n2 <- 8
  data <- complex(length.out= N0 * N1 * N2, real=rnorm(N0 * N1 * N2), imaginary=rnorm(N0 * N1 * N2))
  x <- array(data=data, dim=c(N0, N1, N2))
  a_base <- swdft3d(x=x, n0=n0, n1=n1, n2=n2)

  expect_true( class(a_base)[1] == "swdft3d" )
  expect_true( all( dim(a_base$a) == c(8, 8, 8, 13, 13, 13) ) )
})
