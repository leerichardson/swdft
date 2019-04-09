#' Cosine regression
#'
#' @param x numeric. Signal.
#' @param f numeric. scalar or vector of frequencies to fit.
#'
#' @export
#'
#' @importFrom stats lm coefficients
#'
#' @return S3 object of class 'swdft_cosreg'. See ?new_swdft_cosreg for details.
#'
cosreg <- function(x, f) {
  N <- length(x)

  ## Construct the design matrix for frequency f
  design_matrix <- matrix(data=NA, nrow=N, ncol = 2 * length(f))
  iter <- 0
  for (freq in f) {
    design_matrix[, (2 * iter) + 1] <- cosine(N=N, Fr=freq)
    design_matrix[, (2 * iter) + 2] <- sine(N=N, Fr=freq)
    iter <- iter + 1
  }

  ## Fit the design matrix with least squares
  cosreg_fit <- stats::lm(x ~ design_matrix - 1)
  lm_coefs <- stats::coefficients(object=cosreg_fit)

  ## Extract the amplitudes and phases
  coef_mat <- matrix(data=NA_real_, ncol=3, nrow=length(f))
  colnames(coef_mat) <- c("A", "phi", "f")
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

  ## Return an S3 'swdft_cosreg' object w/ results
  cosreg_obj <- new_swdft_cosreg(coefficients=coef_mat, fitted=fitted, residuals=x-fitted, data=x)
  return(cosreg_obj)
}

#' Constructor function for class swdft_mod
#'
#' @param coefficients matrix of coefficients for cosine regression model
#' @param fitted fitted values of cosine regression model
#' @param residuals residuals of cosine regression model
#' @param data original signal used to fit cosine regression
#'
#' @export
#'
#' @return list with the following elements
#' \itemize{
#'   \item coefficients. A matrix of parameters, the three columns are: 1. amplitude 2. phase, and 3. frequency.
#'   There is only more that one row used when multiple frequencies are fit sequentially.
#'   \item fitted. fitted values of cosine regression model
#'   \item residuals. residuals of cosine regression model
#'   \item data. original signal used to fit cosine regression
#' }
#'
new_swdft_cosreg <- function(coefficients, fitted, residuals, data) {
  structure(list(coefficients=coefficients,
                 fitted=fitted,
                 residuals=residuals,
                 data=data),
            class=c("swdft_cosreg", "swdft_mod"))
}

#' Coefficients method for swdft_cosreg objects
#'
#' @param object A swdft_cosreg object
#' @param ... optional arguments to match generic function
#'
#' @export
#'
coefficients.swdft_mod <- function(object, ...) {
  object$coefficients
}

#' Fitted values method for swdft_cosreg objects
#'
#' @inheritParams coefficients.swdft_mod
#'
#' @export
#'
fitted.swdft_mod <- function(object, ...) {
  object$fitted
}

#' Residuals method for swdft_cosreg objects
#'
#' @inheritParams coefficients.swdft_mod
#'
#' @export
#'
residuals.swdft_mod <- function(object, ...) {
  object$residuals
}
