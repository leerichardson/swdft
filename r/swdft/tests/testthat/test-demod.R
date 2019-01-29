context("Complex Demodulation")

test_that("Demodulation Functions Work", {
  # --- Generate a random signal and SWDFT ---
  N <- 15
  n <- 5
  x <- rnorm(n=N, mean=0, sd=1)
  a <- swdft::swdft(x=x, n=n, type="fftw") * (1 / n)
  f0 <- 3 / n

  ## Verify that the shifted moving average filter matches the SWDFT after shifting
  x_demod_match <- swdft::complex_demod(x=x, f0=f0, smooth='ma', order=n, match_swdft=TRUE)
  expect_true(
    all( round(x_demod_match$y_smooth[3:13], digits=5) == round(a[as.integer(f0 * n) + 1, 5:N], digits=5) )
  )

  ## Verify that we can shift the SWDFT to match the demodulation
  x_demod <- swdft::complex_demod(x=x, f0=f0, smooth='ma', order=n)
  k_swdft_demod <- swdft::demod_swdft(a=a, k=round(f0 * n))
  expect_true(
    all(round(x_demod$y_smooth[3:13], digits=4) == round(k_swdft_demod$k_demod[5:N], digits=4))
  )

})

