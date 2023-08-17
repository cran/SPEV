#' unsmoothed_penalized_EV
#'
#' This function takes a matrix (m), a lambda value (lambda), and the number
#' of desired eigenvectors (k) as input. It then computes
#' eigenvectors 1 to k, penalized by the supplied lambda.
#'
#' @param m A matrix generated from a large dataset.
#' @param lambda A numeric vector of lambda values to use for the penalty.
#' @param k The number of eigenvectors we consider in the analysis.
#' @return Returns eigenvectors 1 to k for the specified lambda value.
#' @importFrom stats optim runif
#' @export
#' @examples
#' # Generate a small matrix for testing
#' m <- matrix(rnorm(100), nrow = 10)
#' # Call function (using matrix, lambda, and k)
#' unsmoothed_penalized_EV(
#'   m = m,
#'   lambda = 1,
#'   k = 2
#' )

unsmoothed_penalized_EV <- function(m, lambda, k) {

  pcs <- matrix(nrow = nrow(m), ncol = k)  # initializing pcs

  gev3 <- function(m, lambda) {
    v <- optim(par = runif(nrow(m)), fn = function(w) {
      w <- w/sqrt(sum(w*w))
      result <- -as.numeric(w %*% m %*% w)
      penalty <- lambda * sum(abs(w))
      return(result + penalty)
    }, method = "BFGS")$par
    v <- v/sqrt(sum(v*v))
    return(v)
  }

  # compute the first eigenvector
  pcs[,1] <- gev3(m, lambda)

  # loop to compute the remaining eigenvectors
  if (k > 1) {
    for (i in 2:k) {
      m_adjusted <- m
      for (j in 1:(i-1)) {
        m_adjusted <- m_adjusted - as.numeric(pcs[,j] %*% m_adjusted %*% pcs[,j]) * outer(pcs[,j], pcs[,j])
      }
      pcs[,i] <- gev3(m_adjusted, lambda)
    }
  }

  return(pcs)
}


####################################################################################################################

#' smoothed_penalized_EV
#'
#' This function takes a matrix (m), a lambda value (lambda), the number
#' of desired eigenvectors (k), and a mu value (mu) as input. It then computes
#' eigenvectors 1 to k, penalized by the supplied lambda and smoothed by the
#' Nesterov smoothing function.
#'
#' @param m A matrix generated from a large dataset.
#' @param lambda A numeric vector of lambda values to use for the penalty.
#' @param k The number of eigenvectors we consider in the analysis.
#' @param mu A number assigned to mu; we are typically using 0.1.
#' @return Returns smoothed eigenvectors 1 to k for the specified lambda value.
#' @importFrom stats optim runif
#' @export
#' @examples
#' # Generate a small matrix for testing
#' m <- matrix(rnorm(100), nrow = 10)
#' # Call function (using matrix, lambda, mu, and k)
#' smoothed_penalized_EV(
#'   m = m,
#'   lambda = 1,
#'   k = 2,
#'   mu = 0.1
#' )

smoothed_penalized_EV <- function(m, lambda, k, mu) {

  pcs <- matrix(nrow = nrow(m), ncol = k)  # initializing pcs

  entropy_prox <- function(z, mu) {
    mu * log((0.5 * exp(-z/mu)) + (0.5 * exp(z/mu)))
  }

  gev4 <- function(m, lambda) {
    v <- optim(par = runif(nrow(m)), fn = function(w) {
      w <- w/sqrt(sum(w*w))
      result <- -as.numeric(w %*% m %*% w)
      penalty <- lambda * sum(entropy_prox(w, mu))
      return(result + penalty)
    }, method = "BFGS")$par
    v <- v/sqrt(sum(v*v))
    return(v)
  }

  # compute the first eigenvector
  pcs[,1] <- gev4(m, lambda)

  # loop to compute the remaining eigenvectors
  if (k > 1) {
    for (i in 2:k) {
      m_adjusted <- m
      for (j in 1:(i-1)) {
        m_adjusted <- m_adjusted - as.numeric(pcs[,j] %*% m_adjusted %*% pcs[,j]) * outer(pcs[,j], pcs[,j])
      }
      pcs[,i] <- gev4(m_adjusted, lambda)
    }
  }

  return(pcs)
}


