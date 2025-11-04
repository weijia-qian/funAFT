#' Function to compute the gradient (score equations)
#'
#' @param params A vector of parameters (gamma + b + sigma)
#' @param Y Survival time
#' @param delta Censoring indicator (1 = event)
#' @param X Design matrix
#' @param family "lognormal" or "loglogistic"
#' @param lambda A chosen nonnegative smoothing parameter
#' @param Pen A chosen penalty matrix
#'
#' @returns A vector of negative penalized scores

penalized_score <- function(params, Y, delta, X, family = c("lognormal", "loglogistic"), lambda, Pen) {
  # Extract parameters
  b_coef <- params[-length(params)]
  sigma <- params[length(params)]
  family <- match.arg(family)

  mu <- X %*% b_coef
  z <- (log(Y) - mu) / sigma

  if (family == "lognormal"){
    f_z <- dnorm(z) # PDF of z
    S_z <- pnorm(-z)  # survival function of z

    # Score for b
    score_b <- t(X) %*% (delta * z / sigma + (1 - delta) * (f_z / (S_z * sigma))) - 2 * lambda * Pen %*% b_coef

    # Score for sigma
    score_sigma <- sum(delta * (-1 / sigma + z^2 / sigma) + (1 - delta) * f_z * z / (sigma * S_z))

  } else if (family == "loglogistic"){
    # Logistic CDF
    p_z <- 1 / (1 + exp(-z))

    # Score for b
    score_b <-  crossprod(X, delta * (2 * p_z - 1) + (1 - delta) * p_z) / sigma - 2 * lambda * Pen %*% b_coef

    # Score for sigma
    score_sigma <- sum(- delta + delta * (2 * p_z - 1) * z + (1 - delta) * p_z * z) / sigma
  }

  return(-c(score_b, score_sigma)) # negative score for minimization
}
