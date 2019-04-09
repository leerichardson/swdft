context("Visualization")

test_that("Visualizations work", {
  N <- 60
  n <- 32
  signal <- swdft::local_signal(N=N, A=1, Fr=4/32, phase=1, S=10, L=40)
  noise <- rnorm(n=N, mean=0, sd=.2)
  x <- signal + noise
  a <- swdft::swdft(x=x, n=32, type="fftw")
  expect_true( plot(a, display=FALSE) == "display set to false, 'plot.swdft' runs without errors" )
})
