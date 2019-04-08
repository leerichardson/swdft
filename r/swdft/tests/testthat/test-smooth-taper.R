context("Smoothing and Tapering")

test_that("Verify that smoothing and tapering the SWDFT works", {
  N <- 100
  window_size <- 2^5
  x <- rnorm(n=N)

  ## Take smoothed and tapered SWDFTs
  asmooth <- swdft::swdft(x=x, n=window_size, type="fftw", smooth="daniell")
  ataper <- swdft::swdft(x=x, n=window_size, type="fftw", taper_type='cosine')

  ## Check that the class and dimensions work out properly
  expect_true(typeof(asmooth$a) == "double")
  expect_true(typeof(ataper$a) == "complex")
  expect_true(all(dim(asmooth$a) == c(window_size, N)))
  expect_true(all(dim(ataper$a) == c(window_size, N)))
})
