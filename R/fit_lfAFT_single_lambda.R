#' Fit Linear Functional Accelerated Failure Time (lfAFT) Models with A Given Lambda
#'
#' Fits a linear functional accelerated failure time (lfAFT) model
#' for right-censored survival outcomes with functional and scalar covariates
#' with a given smoothing parameter (\eqn{\lambda}).
#'
#' @param data Optional \code{data.frame} containing all variables in the model.
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
#' @param s_grid Optional numeric vector of the same length as the number of
#'               functional measurement points (\code{ncol(X)}). If \code{NULL},
#'               assumes an equally spaced grid over \[0, 1\].
#' @param k Integer; number of spline basis functions used to represent
#'          \eqn{\beta(s)} (effective degrees of freedom for the functional coefficient).
#' @param lambda Numeric; smoothing parameter controlling the penalty on
#'               roughness of \eqn{\beta(s)}.
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
#'   \item{\code{GCV}}{Generalized cross-validation score.}
#'   \item{\code{se_*}}{Wald standard errors (if \code{se=TRUE}).}
#'   \item{\code{*_ci_lower}, \code{*_ci_upper}}{Wald or bootstrap
#'         95 % confidence intervals.}
#'   \item{\code{family, lambda, basis, kbasis, s_grid}}{Model specifications.}
#'   \item{\code{convergence, value, message}}{Optimizer diagnostics.}
#' }
#'
#' @importFrom stats as.formula dnorm optim optimHess pnorm qnorm
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
#'                  lambda = 10, k = 10, se = TRUE)
#'
#' plot(fit$s_grid, fit$betaX_hat, type = "l",
#'      xlab = "s", ylab = "Estimated beta(s)")


fit_lfAFT_single_lambda <- function(data = NULL, y, delta,x, x_as_regex = FALSE, z = NULL, family = c("lognormal", "loglogistic"), s_grid = NULL,
                      lambda, k = 20, basis = c("bs", "ns"), basis_args = list(), se = FALSE, bootstrap = FALSE, B = 500, boot_seed = NULL) {

  family <- match.arg(family)
  basis  <- match.arg(basis)

  .env <- parent.frame()
  if (!is.null(data) && !.is_df(data)) stop("`data` must be a data.frame or NULL.")

  # ---- pull Y, Delta ----
  Y <- .get_vec(y, data, .env, required_n = NULL, label = "y")
  n <- length(Y)
  Delta <- .get_vec(delta, data, .env, required_n = n, label = "delta")

  # ---- select/build X (functional, wide) ----
  if (length(x) == 1L && isTRUE(x_as_regex)) {
    discovered <- .discover_x_by_regex(x, data, .env)
    if (!length(discovered)) stop("No X candidates matched regex/pattern: ", x)
    # If matches are data columns, bind those; otherwise bind env objects
    if (!is.null(data) && all(discovered %in% names(data))) {
      X <- as.matrix(data[, discovered, drop = FALSE])
      x_cols <- colnames(X)
    } else {
      X <- .bind_cols(discovered, data = NULL, env = .env, required_n = n, label = "X")
      x_cols <- colnames(X)
    }
  } else {
    # Names explicitly given: each can be a vector (1 col) or matrix/df (multi-col)
    X <- .bind_cols(x, data = data, env = .env, required_n = n, label = "X")
    x_cols <- colnames(X)
  }
  X <- as.matrix(X)
  nS <- ncol(X)

  # ---- Z (scalar) ----
  if (!is.null(z)) {
    Z_mat <- .bind_cols(z, data = data, env = .env, required_n = n, label = "Z")
    Z_names <- colnames(Z_mat)
  } else {
    Z_mat <- NULL
    Z_names <- NULL
  }
  nZ <- if (is.null(Z_mat)) 0L else ncol(Z_mat)

  # ---- s-grid ----
  if (is.null(s_grid)) {
    s <- seq(0, 1, length.out = nS)
  } else {
    if (length(s_grid) != nS) stop("`s_grid` must have length = ncol(X) = ", nS, ".")
    s <- as.numeric(s_grid)
  }

  # ---- build bs/ns basis matrix on s (df = k) ----
  if (basis == "bs") {
    Bspl <- do.call(splines::bs, c(list(x = s, df = k), basis_args))
  } else {
    Bspl <- do.call(splines::ns, c(list(x = s, df = k), basis_args))
  }
  kb <- ncol(Bspl)  # typically == k

  # ---- projected functional design (trapezoid weights) ----
  # q <- (s[nS] - s[1]) / (nS - 1)     # quadrature weight
  # X_mat <- X %*% Bspl * q         # n x kb
  d <- diff(s)
  w <- c(d[1]/2, (d[-length(d)] + d[-1])/2, d[length(d)]/2)
  Bs_w <- sweep(Bspl, 1, w, "*")  # apply weights to rows
  X_mat <- X %*% Bs_w

  # ---- design matrix: intercept + Z + X_mat ----
  C <- cbind(1, Z_mat, X_mat)
  p <- ncol(C)

  # ---- penalty matrix on functional block ----
  Pen <- matrix(0, nrow = p, ncol = p)
  Pen_block <- penalty_matrix(kp = kb, nS = nS, a = 0.001)
  Pen[(2 + nZ):p, (2 + nZ):p] <- Pen_block

  # ---- optimization ----
  init_params <- c(rep(0, p), 1)   # last is scale parameter
  fit <- optim(par     = init_params,
               fn      = penalized_loglik,
               gr      = penalized_score,
               method  = "BFGS",
               Y       = Y,
               delta   = Delta,
               X       = C,
               family  = family,
               lambda  = lambda,
               Pen     = Pen,
               control = list(maxit = 2000))

  # ---- estimates ----
  coef_all  <- fit$par[1:p]
  beta0_hat <- coef_all[1]
  betaZ_hat <- if (nZ > 0) coef_all[2:(1 + nZ)] else NULL
  betaX_hat <- as.numeric(coef_all[(2 + nZ):p] %*% t(Bspl))
  sigma_hat <- fit$par[p + 1]
  mu_hat    <- as.vector(C %*% coef_all)

  # ---- GCV ----
  GCV <- compute_gcv(Y, Delta, C, mu_hat, sigma_hat, lambda, Pen, family)

  out <- list(
    beta0_hat = beta0_hat,
    betaZ_hat = betaZ_hat,
    betaX_hat = betaX_hat,
    sigma_hat = sigma_hat,
    lp        = mu_hat,
    GCV       = GCV,
    family    = family,
    lambda    = lambda,
    Z_names   = Z_names,
    X_names   = x_cols,
    s_grid    = s,
    basis     = basis,
    kbasis    = kb,
    convergence = fit$convergence,
    value        = fit$value,
    message      = fit$message %||% NULL
  )

  # ---- Wald SEs / CIs (optional) ----
  if (isTRUE(se)) {
    hessian  <- optimHess(fit$par, fn = penalized_loglik, gr = penalized_score,
                          Y = Y, delta = Delta, X = C, family = family,
                          lambda = lambda, Pen = Pen)
    cov_beta <- tryCatch(solve(hessian), error = function(e) NULL)

    if (is.null(cov_beta)) {
      warning("Hessian inversion failed; Wald SEs unavailable.")
    } else {
      se_beta0 <- sqrt(diag(cov_beta)[1])
      se_betaZ <- if (nZ > 0) sqrt(diag(cov_beta)[2:(1 + nZ)]) else NULL
      se_betaX <- sqrt(rowSums((Bspl %*% cov_beta[(2 + nZ):p, (2 + nZ):p]) * Bspl))
      se_sigma <- sqrt(diag(cov_beta)[p + 1])

      zq <- qnorm(0.975)
      out$beta0_ci_lower <- beta0_hat - zq * se_beta0
      out$beta0_ci_upper <- beta0_hat + zq * se_beta0

      if (nZ > 0) {
        out$betaZ_ci_lower <- as.numeric(betaZ_hat - zq * se_betaZ)
        out$betaZ_ci_upper <- as.numeric(betaZ_hat + zq * se_betaZ)
        names(out$betaZ_ci_lower) <- Z_names
        names(out$betaZ_ci_upper) <- Z_names
      }

      out$betaX_ci_lower <- betaX_hat - zq * se_betaX
      out$betaX_ci_upper <- betaX_hat + zq * se_betaX
      out$sigma_ci_lower <- sigma_hat - zq * se_sigma
      out$sigma_ci_upper <- sigma_hat + zq * se_sigma

      out$se_beta0 <- se_beta0
      out$se_betaZ <- se_betaZ
      out$se_betaX <- se_betaX
      out$se_sigma <- se_sigma
    }
  }

  # ---- Bootstrap CIs (optional) ----
  if (isTRUE(bootstrap)) {
    if (!is.null(boot_seed)) set.seed(boot_seed)

    beta0_boot <- numeric(B)
    sigma_boot <- numeric(B)
    betaX_boot <- matrix(NA_real_, nrow = B, ncol = nS)
    betaZ_boot <- if (nZ > 0) matrix(NA_real_, nrow = B, ncol = nZ) else NULL

    for (bb in seq_len(B)) {
      idx <- sample.int(n, size = n, replace = TRUE)
      Yb     <- Y[idx]
      Deltab <- Delta[idx]
      Xb     <- X[idx, , drop = FALSE]
      Zb     <- if (nZ > 0) Z_mat[idx, , drop = FALSE] else NULL

      # reuse same Bspl and penalty
      X_mat_b <- Xb %*% Bs_w
      Cb      <- cbind(1, Zb, X_mat_b)

      fit_b <- tryCatch(
        optim(par     = init_params,
              fn      = penalized_loglik,
              gr      = penalized_score,
              method  = "BFGS",
              Y       = Yb,
              delta   = Deltab,
              X       = Cb,
              family  = family,
              lambda  = lambda,
              Pen     = Pen,
              control = list(maxit = 2000)),
        error = function(e) NULL
      )

      if (is.null(fit_b) || fit_b$convergence != 0 || any(!is.finite(fit_b$par))) next

      coef_b <- fit_b$par[1:p]
      beta0_boot[bb] <- coef_b[1]
      if (nZ > 0) betaZ_boot[bb, ] <- coef_b[2:(1 + nZ)]
      betaX_boot[bb, ] <- as.numeric(coef_b[(2 + nZ):p] %*% t(Bspl))
      sigma_boot[bb]   <- fit_b$par[p + 1]
    }

    qfun <- function(x) stats::quantile(x, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE)

    out$bootstrap_info <- list(
      B    = B,
      used = sum(is.finite(beta0_boot) & is.finite(sigma_boot))
    )

    if (isTRUE(bootstrap) && out$bootstrap_info$used < B/2)
      warning("Fewer than half of bootstrap fits converged; CIs may be unreliable.")

    beta0_ci <- qfun(beta0_boot)
    out$beta0_boot_ci_lower <- beta0_ci[1]; out$beta0_boot_ci_upper <- beta0_ci[2]

    if (nZ > 0) {
      Z_ci <- apply(betaZ_boot, 2, qfun)
      out$betaZ_boot_ci_lower <- as.numeric(Z_ci[1, ])
      out$betaZ_boot_ci_upper <- as.numeric(Z_ci[2, ])
      names(out$betaZ_boot_ci_lower) <- Z_names
      names(out$betaZ_boot_ci_upper) <- Z_names
    }

    betaX_ci <- apply(betaX_boot, 2, qfun)
    out$betaX_boot_ci_lower <- as.numeric(betaX_ci[1, ])
    out$betaX_boot_ci_upper <- as.numeric(betaX_ci[2, ])

    sigma_ci <- qfun(sigma_boot)
    out$sigma_boot_ci_lower <- sigma_ci[1]; out$sigma_boot_ci_upper <- sigma_ci[2]
  }

  out
}

