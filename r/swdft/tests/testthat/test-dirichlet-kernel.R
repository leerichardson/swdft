context("Dirichlet Kernel")

test_that("Dirichlet Kernel works", {
  # Tests for all of the derivations in the Trigonometrix Identities Section of the Thesis
  n <- 10
  i <- complex(length.out=1, real=0, imaginary=1)
  expect_true(dirichlet_kernel(x=(1*pi), n=n) == 0)
  expect_true(dirichlet_kernel(x=(-4*pi), n=n) == n)
  expect_true(dirichlet_kernel(x=(6*pi), n=n) == -n)
  expect_true(dirichlet_kernel(x=(1), n=n) == sin((n * 1) / 2)  / sin(1 / 2))

  # Cosine and Sine identities
  x <- 2
  j <- 0:(n-1)
  phase <- 1.2
  expect_true( round(sum( cos(x = (j * x)) ), digits=4) == round(cos(x = ((n - 1) / 2) * x) * dirichlet_kernel(x=x, n=n), digits=4) )
  expect_true( round(sum( sin(x = (j * x)) ), digits=4) == round(sin(x = ((n - 1) / 2) * x) * dirichlet_kernel(x=x, n=n), digits=4) )

  sum( exp( x = i * ((j * x) + phase) ) )
  exp( x = i * ( (((n - 1) / 2) * x) + phase) ) * dirichlet_kernel(x=x, n=n)

  sum( cos(x = (( (j * x) + phase))) )
  cos(x = (((n - 1) / 2) * x) + phase) * dirichlet_kernel(x=x, n=n)

  sum( sin(x = (( (j * x) + phase))) )
  sin(x = (((n - 1) / 2) * x) + phase) * dirichlet_kernel(x=x, n=n)

  # The Dirichlet Weight
  expect_true(round(sum( exp(i * j * x) ), digits=4) == round(dirichlet_kernel(x=x, n=n, dw=TRUE), digits=4))

  # Complex Summation Identity ---
  a <- 2
  b <- 10
  k <- a:b
  sum_vals <- exp(x = i * k * x)
  complex_sum <- sum(sum_vals)
  complex_sum_identity <- exp(x = i * a * x) * exp(x = i * ((b - a) / 2) * x) * dirichlet_kernel(x=x, n=b-a+1)
  expect_true( round(complex_sum, digits=4) == round(complex_sum_identity, digits=4) )

  # Cosine and Sine Partial Sum identities
  expect_true(round(sum(cos(x = k * x)), digits=4) ==  round(cos(x = ((a+b) / 2) * x) * dirichlet_kernel(x=x, n=b-a+1), digits=4))
  expect_true(round(sum(sin(x = k * x)), digits=4) ==  round(sin(x = ((a+b) / 2) * x) * dirichlet_kernel(x=x, n=b-a+1), digits=4))
})
