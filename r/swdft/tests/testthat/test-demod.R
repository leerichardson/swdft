context("Complex and matching demodulation")

test_that("Complex demodulation matches SWDFT w/ moving average filter", {
  # --- Generate white-noise signal and corresponding SWDFT ---
  N <- 15
  window_size <- 5
  x <- rnorm(n=N, mean=0, sd=1)
  a <- swdft::swdft(x=x, n=window_size, type="fftw") * (1 / window_size)
  f0 <- sample(x=1:3, size=1) / window_size

  # --- Tests ---

  ## Verify that the shifted moving average filter matches the SWDFT after shifting
  x_demod_match <- swdft::complex_demod(x=x, f0=f0, smooth='ma', order=window_size, match_swdft=TRUE, window_size=window_size)
  expect_true(
    all( round(x_demod_match$demod$y_smooth[3:13], digits=5) == round(a[as.integer(f0 * window_size) + 1, 5:N], digits=5) )
  )

  ## Verify that we can shift the SWDFT to match the demodulation
  x_demod <- swdft::complex_demod(x=x, f0=f0, smooth='ma', order=window_size)
  k_swdft_demod <- swdft::demod_swdft(a=a, k=round(f0 * window_size))
  expect_true(
    all.equal(target=round(k_swdft_demod$demod, digits=3), current=round(x_demod$demod$y_smooth,digits=3))
  )

  ## Other tests on class and length of outputs
  expect_true(all(class(x_demod_match) == c("swdft_demod", "swdft_mod")))
  expect_true(length(x_demod_match$coefficients$inst_amp) == length(x))
})

test_that("Matching demodulation works", {
  print("Add matching demodulation tests!")
})
