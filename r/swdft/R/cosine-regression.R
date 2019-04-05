#' Cosine regression
#'
#' @param x numeric. Signal
#' @param f numeric. Single of vector of frequenies to fit. Defaults to NULL and
#' in this case the maximum DFT coefficient in will be used.
#'
#' @return S3 object of class 'swdft_cosreg'
#'
cosreg <- function(x, f=NULL) {
  N <- length(x)

  ## Construct the design matrix for frequency f
  design_matrix <- matrix(data=NA, nrow=N, ncol = 2 * length(f))
  iter <- 0
  for (freq in f) {
    design_matrix[, (2 * iter) + 1] <- swdft::cosine(N=N, Fr=freq)
    design_matrix[, (2 * iter) + 2] <- swdft::sine(N=N, Fr=freq)
    iter <- iter + 1
  }

  ## Fit the design matrix with least squares
  cosreg_fit <- lm(x ~ design_matrix - 1)
  lm_coefs <- coefficients(object=cosreg_fit)

  ## Extract the amplitudes and phases
  coef_mat <- matrix(data=NA_real_, ncol=3, nrow=length(f))
  fitted <- vector(mode="numeric", length=N)

  iter <- 0
  for (freq in f) {
    ## Store the amplitude, phase, and frequency in the covariance matrix
    cval <- complex(length.out=1, real=lm_coefs[(2 * iter) + 1], imaginary=lm_coefs[(2 * iter) + 2])
    coef_mat[iter + 1, 1] <- sqrt(Mod(cval))
    coef_mat[iter + 1, 2] <- atan2(-Im(cval), Re(cval))
    coef_mat[iter + 1, 3] <- freq

    ## Update the fitted values w/ the next cosine term
    fitted <- fitted + (Re(cval) * design_matrix[, (2 * iter) + 1] + Im(cval) * design_matrix[, (2 * iter) + 2])

    iter <- iter + 1
  }

  ## Compute the residuals
  residuals <- fitted - x

  ## Compute the residuals and create a 'swdft_cosreg' object
  cosreg_obj <- structure(list(coefficients=coef_mat,
                               fitted=fitted,
                               residuals=residuals,
                               data=x),
                          class="swdft_cosreg")
  return(cosreg_obj)
}
