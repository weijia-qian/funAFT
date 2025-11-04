#' Linear Functional Accelerated Failure Time (lfAFT) Models
#'
#' Fits a linear functional accelerated failure time (lfAFT) model
#' for right-censored survival data with functional and scalar covariates.
#'
#' @param data Optional \code{data.frame} containing all variables used in the model.
#'             If \code{NULL}, variables are searched in the calling environment.
#' @param y Character. Name of the variable for survival time.
#' @param delta Character. Name of the censoring indicator
#'              (1 = observed, 0 = censored).
#' @param x Character vector or single regex/prefix pattern.
#'          - If \code{x_as_regex = FALSE} (default), supply the exact column names
#'            for the functional covariate.
#'          - If \code{x_as_regex = TRUE} and \code{length(x) == 1}, all columns in
#'            \code{data} (or objects in the environment) matching the pattern are used.
#' @param x_as_regex Logical. If \code{TRUE} and \code{length(x) == 1},
#'                   treat \code{x} as a regular expression or prefix to match
#'                   functional covariate columns.
#' @param z Optional character vector of scalar covariate names.
#'          If \code{NULL}, the model includes no scalar covariates.
#' @param family Character; either \code{"lognormal"} or \code{"loglogistic"}.
#'               Determines the error distribution of the AFT model.
#' @param k Integer; number of spline basis functions used to represent
#'          \eqn{\beta(s)} (effective degrees of freedom for the functional coefficient).
#' @param lambda_grid Numeric vector of candidate smoothing parameters.
#'                    Defaults to \code{exp(seq(log(1000), log(10000), length.out = 100))},
#'                    which creates 100 exponentially spaced values from 1e3 to 1e4.
#'                    The function evaluates the generalized cross-validation (GCV)
#'                    score for each value and selects the one minimizing GCV.
#' @param s_grid Optional numeric vector of the same length as the number of
#'               functional measurement points (\code{ncol(X)}). If \code{NULL},
#'               assumes an equally spaced grid over \[0, 1\].
#' @param basis Character; one of \code{c("bs", "ns")}. Determines the spline
#'              basis type used for the functional coefficient. Defaults to
#'              \code{"bs"} (B-spline).
#' @param basis_args List of additional arguments passed to the spline basis
#'                   constructor (e.g., \code{degree = 3} for \code{bs},
#'                   or \code{intercept = FALSE}).
#' @param se Logical. If \code{TRUE}, computes Wald standard errors and 95\%
#'           confidence intervals based on the inverse Hessian.
#' @param bootstrap Logical. If \code{TRUE}, computes percentile-based
#'                  bootstrap confidence intervals via subject-level
#'                  resampling.
#' @param B Integer; number of bootstrap replicates (default 500).
#' @param boot_seed Optional integer seed for reproducibility in bootstrapping.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{beta0_hat}}{Estimated intercept.}
#'   \item{\code{betaZ_hat}}{Estimated scalar covariate coefficients (if any).}
#'   \item{\code{betaX_hat}}{Estimated functional coefficient curve
#'         \eqn{\hat{\beta}(s)} evaluated on the supplied grid.}
#'   \item{\code{sigma_hat}}{Estimated scale parameter.}
#'   \item{\code{lp}}{Linear predictor values (\eqn{\hat{\mu}_i}).}
#'   \item{\code{GCV}}{Generalized cross-validation score for the optimal \eqn{\lambda}.}
#'   \item{\code{lambda_grid}}{Grid of candidate lambda values.}
#'   \item{\code{lambda_gcv}}{Vector of GCV scores across all grid values.}
#'   \item{\code{se_*}}{Wald standard errors (if \code{se=TRUE}).}
#'   \item{\code{*_ci_lower}, \code{*_ci_upper}}{Wald or bootstrap 95%
#'         confidence intervals.}
#'   \item{\code{family, lambda, basis, kbasis, s_grid}}{Model specifications.}
#'   \item{\code{convergence, value, message}}{Optimizer diagnostics.}
#' }
#'
#' @export
#'
#' @examples
#' # Example with synthetic data
#' set.seed(1)
#' n <- 100; s <- seq(0, 1, length.out = 50)
#' beta_true <- sin(2 * pi * s)
#' X <- matrix(rnorm(n * 50), n, 50)
#' mu <- 1 + X %*% beta_true * 0.5
#' Y <- exp(mu + rnorm(n))
#' Delta <- rbinom(n, 1, 0.8)
#' dat <- data.frame(Y = Y, Delta = Delta, X)
#'
#' fit <- fit_lfAFT(data = dat,
#'                  y = "Y", delta = "Delta",
#'                  x = "^X", x_as_regex = TRUE,
#'                  family = "lognormal",
#'                  k = 10,
#'                  lambda_grid = 10^seq(-3, 3, length.out = 10),
#'                  se = TRUE)
#'
#' plot(fit$s_grid, fit$betaX_hat, type = "l",
#'      xlab = "s", ylab = "Estimated beta(s)")


fit_lfAFT <- function(data = NULL,
                      y,
                      delta,
                      x,
                      x_as_regex = FALSE,
                      z = NULL,
                      family = c("lognormal", "loglogistic"),
                      k = 20,
                      lambda_grid = exp(seq(log(1000), log(10000), length.out = 100)),
                      s_grid = NULL,
                      basis = c("bs", "ns"),
                      basis_args = list(),
                      se = FALSE,
                      bootstrap = FALSE,
                      B = 500,
                      boot_seed = NULL) {

  # 1) GCV-based tuning of lambda
  opt <- optimize_lambda(
    lambda_grid = lambda_grid,
    data   = data,
    y      = y,
    delta  = delta,
    x      = x,
    x_as_regex = x_as_regex,
    z      = z,
    family = family,
    k      = k,
    s_grid = s_grid,
    basis  = basis,
    basis_args = basis_args
  )

  best_lambda <- opt$lambda

  # 2) Refit at the optimal lambda, now allowing se / bootstrap
  fit <- fit_lfAFT_single_lambda(
    data   = data,
    y      = y,
    delta  = delta,
    x      = x,
    x_as_regex = x_as_regex,
    z      = z,
    family = family,
    k      = k,
    lambda = best_lambda,
    s_grid = s_grid,
    basis  = basis,
    basis_args = basis_args,
    se        = se,
    bootstrap = bootstrap,
    B         = B,
    boot_seed = boot_seed
  )

  # 3) Attach tuning info
  fit$lambda_grid       <- lambda_grid
  fit$lambda_gcv        <- opt$gcv

  fit
}
