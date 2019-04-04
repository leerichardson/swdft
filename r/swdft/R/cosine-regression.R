#' Local Cosine Regression
#'
#' @param x numeric. Signal
#' @param f numeric. Frequency or vector of frequencies
#'
cosreg <- function(x, f) {
  ## Create design matrix from frequencies
  N <- length(x)
  design_matrix <- matrix(data=NA, nrow=N, ncol = 2 * length(f))

  iter <- 0
  for (freq in f) {
    design_matrix[, (2 * iter) + 1] <- swdft::cosine(N=N, Fr=freq)
    design_matrix[, (2 * iter) + 2] <- swdft::sine(N=N, Fr=freq)
    iter <- iter + 1
  }

  ## Fit the model and return the parameters
  cosreg_fit <- lm(x ~ design_matrix - 1)

  ## Extract the amplitudes and phases
  coefs <- coefficients(object=cosreg_fit)
  amps <- vector(mode="numeric", length=length(f))
  phases <- vector(mode="numeric", length=length(f))
  fitted <- vector(mode="numeric", length=N)

  iter <- 0
  for (freq in f) {
    cval <- complex(length.out=1, real=coefs[(2 * iter) + 1], imaginary=coefs[(2 * iter) + 2])

    amps[iter + 1] <- Mod(cval)
    phases[iter + 1] <- Arg(cval)
    fitted <- fitted + (Re(cval) * design_matrix[, (2 * iter) + 1] + Im(cval) * design_matrix[, (2 * iter) + 2])

    iter <- iter + 1
  }

  return(fitted)
}
