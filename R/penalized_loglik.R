#' Function to compute the negative log-likelihood
#'
#' @param params A vector of parameters (gamma + b + sigma)
#' @param Y Survival time
#' @param delta Censoring indicator (1 = event)
#' @param X Design matrix
#' @param family "lognormal" or "loglogistic"
#' @param lambda A chosen nonnegative smoothing parameter
#' @param Pen A chosen penalty matrix
#'
#' @returns Negative penalized log-likelihood

penalized_loglik <- function(params, Y, delta, X, family = c("lognormal", "loglogistic"), lambda, Pen) {
  # Extract parameters
  b_coef <- params[-length(params)]
  sigma <- params[length(params)]
  family <- match.arg(family)

  # Ensure sigma > 0
  if (sigma <= 0 || !is.finite(sigma)) return(Inf)

  mu <- X %*% b_coef
  z <- (log(Y) - mu) / sigma

  if (family == "lognormal"){
    log_f <- delta * (-log(Y) - log(sigma) - 0.5 * z^2 - 0.5 * log(2 * pi))
    log_S <- (1 - delta) * pnorm(z, lower.tail = FALSE, log.p = TRUE)
  } else if (family == "loglogistic"){
    p_z <- 1 / (1 + exp(-z)) # logistic cdf
    log_f <- delta * log((1 / (sigma * Y)) * p_z * (1 - p_z))
    log_S <- (1 - delta) * log(1 - p_z)
  }

  # Penalized log-likelihood
  penalty <- lambda * crossprod(b_coef, Pen) %*% b_coef
  loglik <- sum(log_f + log_S) - penalty

  # Return negative log-likelihood for minimization
  if (!is.finite(loglik)) return(Inf)
  return(-loglik)
}
