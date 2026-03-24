#' Function to compute generalized cross-validation score
#'
#' @param Y Survival time
#' @param delta Censoring indicator (1 = event)
#' @param C Design matrix
#' @param mu_hat Estimated linear predictor
#' @param sigma_hat Estimated scale parameter
#' @param lambda A chosen nonnegative smoothing parameter
#' @param Pen A chosen penalty matrix
#' @param family "lognormal" or "loglogistic"
#'
#' @returns Generalized cross-validation score
#'
compute_gcv <- function(Y, delta, C, mu_hat, sigma_hat, lambda, Pen, family = c("lognormal", "loglogistic")) {
  family <- match.arg(family)
  n <- length(Y)

  # Log-likelihood
  loglik <- switch(family,
                   "lognormal" = {
                     z <- (log(Y) - mu_hat) / sigma_hat
                     log_f <- dnorm(z, log = TRUE) - log(Y * sigma_hat)
                     log_S <- pnorm(z, lower.tail = FALSE, log.p = TRUE)
                     sum(delta * log_f + (1 - delta) * log_S)
                   },
                   "loglogistic" = {
                     z <- (log(Y) - mu_hat) / sigma_hat
                     log_f <- log(1 / sigma_hat) - log(Y) - 2 * log(1 + exp(-z))
                     log_S <- -log(1 + exp(z))
                     sum(delta * log_f + (1 - delta) * log_S)
                   }
  )

  # Weights for DF calculation (approximate second derivatives)
  w <- switch(family,
              "lognormal" = {
                z <- (log(Y) - mu_hat) / sigma_hat
                phi_z <- dnorm(z)
                S_z <- pnorm(z, lower.tail = FALSE)
                w <- numeric(n)
                w[delta == 1] <- 1 / sigma_hat^2
                w[delta == 0] <- (phi_z[delta == 0]^2) / (S_z[delta == 0]^2 * sigma_hat^2)
                w
              },
              "loglogistic" = {
                z <- (log(Y) - mu_hat) / sigma_hat
                # Approximate weight as derivative of log-hazard:
                p <- exp(z) / (1 + exp(z))  # derivative of log(1 + exp(z)) is logistic
                w <- numeric(n)
                w[delta == 1] <- 1 / sigma_hat^2  # constant approx
                w[delta == 0] <- (p[delta == 0]^2) / sigma_hat^2  # rough approx
                w
              }
  )

  # Clip weights for stability
  w <- pmin(pmax(w, 1e-6), 1e6)

  # Degrees of freedom
  # IMPROVEMENT: crossprod(C, w * C) computes X^T W X without building an N x N diagonal matrix
  XtWX <- crossprod(C, w * C)

  df <- tryCatch({
    # Primary attempt: Fast standard solve
    H <- solve(XtWX + lambda * Pen, XtWX)
    sum(diag(H))
  }, error = function(e) {
    # Fallback attempt: Generalized inverse handles singular matrices (empty factor levels)
    tryCatch({
      suppressPackageStartupMessages(require(MASS))
      H <- ginv(XtWX + lambda * Pen) %*% XtWX
      sum(diag(H))
    }, error = function(e2) {
      warning("Both standard and generalized matrix inversion failed during GCV calculation.")
      return(NA)
    })
  })

  # Guard against invalid df
  if (is.na(df) || df <= 0 || df / n > 0.95) {
    return(Inf)
  }

  # GCV
  gcv <- -loglik / (1 - df / n)^2
  return(gcv)
}
