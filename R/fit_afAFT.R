#' Additive Functional Accelerated Failure Time (afAFT) Models
#'
#' Fit an additive functional regression model for survival data using
#' \pkg{mgcv}, where the functional predictor \eqn{X(s)} enters the model
#' through an additive term of the form
#' \deqn{
#'   \int F\{X(s), s\}\, ds,
#' }
#' with an unknown bivariate smooth function \eqn{F(\cdot,\cdot)}.  This is
#' implemented using a tensor-product smooth over the functional covariate and
#' its domain.
#'
#' Currently only the log-normal AFT model is supported (see \url{https://rdrr.io/cran/mgcv/man/cnorm.html}).
#'
#' @param data Optional \code{data.frame} containing variables referenced by
#'   \code{y}, \code{delta}, \code{x}, and \code{z}. If \code{NULL}, the
#'   variables are looked up in the calling environment.
#' @param y Survival time variable. Either a character string (column name in
#'   \code{data}) or an object name in the calling environment.
#' @param delta Censoring indicator (1 = event, 0 = right-censored).
#' @param x Functional covariate specification. Either:
#'   \itemize{
#'     \item A character vector of names to be column-bound into the wide
#'           functional matrix \code{X}, or
#'     \item A single regex/prefix (with \code{x_as_regex = TRUE}) used to
#'           discover matching columns or objects.
#'   }
#' @param x_as_regex Logical; if \code{TRUE} and \code{length(x) == 1}, treat
#'   \code{x} as a regex/prefix and discover matching functional columns.
#' @param z Optional scalar covariates (may be \code{NULL}).  Typically a
#'   character vector of column/object names.
#' @param family Character. Determines the error distribution of the AFT model;
#'   currently only support \code{"lognormal"}.
#' @param s_grid Optional numeric vector giving the functional grid
#'   \eqn{s_1,\ldots,s_{n_S}} (must match \code{ncol(X)}). If \code{NULL}, an
#'   equally spaced grid on \[0, 1\] is used.
#' @param k Integer vector of length 1 or 2. Basis dimensions for the
#'   tensor-product smooth \code{ti(S, X, ...)}; recycled to length 2 if a
#'   single integer is provided.
#' @param basis Character vector of length 1 or 2 specifying the marginal
#'   spline bases for \code{\link[mgcv]{ti}}.  Defaults to
#'   \code{c("cr","cr")}.
#'
#' @return A fitted \code{\link[mgcv]{gam}} object with additional elements:
#' \itemize{
#'   \item \code{X_names}, \code{Z_names}: variable names used.
#'   \item \code{s_grid}: the functional grid employed.
#'   \item \code{family_name}: chosen family label.
#' }
#'
#' @examples
#' \dontrun{
#' library(mgcv)
#' set.seed(1)
#' n <- 100; nS <- 50
#' dat <- data.frame(
#'   time = rexp(n, 0.1),
#'   status = rbinom(n, 1, 0.7),
#'   age = rnorm(n, 60, 8)
#' )
#' Xfun <- matrix(rnorm(n*nS), n, nS)
#' colnames(Xfun) <- paste0("X", 1:nS)
#' dat <- cbind(dat, Xfun)
#'
#' fit <- fit_afAFT(
#'   data = dat, y = "time", delta = "status",
#'   x = "^X", x_as_regex = TRUE,
#'   z = "age")
#' summary(fit)
#' }
#'
#' @export

fit_afAFT <- function(data = NULL,
                      y,
                      delta,
                      x,
                      x_as_regex = FALSE,
                      z = NULL,
                      family = c("lognormal"),
                      k = c(20, 20),
                      s_grid = NULL,
                      basis = c("cr", "cr")) {

  family <- match.arg(family)

  .env <- parent.frame()
  if (!is.null(data) && !.is_df(data)) {
    stop("`data` must be a data.frame or NULL.")
  }

  # ---- pull Y, Delta ----
  Y <- .get_vec(y, data, .env, required_n = NULL, label = "y")
  n <- length(Y)
  Delta <- .get_vec(delta, data, .env, required_n = n, label = "delta")

  # ---- select/build X (functional, wide) ----
  if (length(x) == 1L && isTRUE(x_as_regex)) {
    discovered <- .discover_x_by_regex(x, data, .env)
    if (!length(discovered)) stop("No X candidates matched regex/pattern: ", x)

    if (!is.null(data) && all(discovered %in% names(data))) {
      X <- as.matrix(data[, discovered, drop = FALSE])
      x_cols <- colnames(X)
    } else {
      X <- .bind_cols(discovered, data = NULL, env = .env,
                      required_n = n, label = "X")
      x_cols <- colnames(X)
    }
  } else {
    X <- .bind_cols(x, data = data, env = .env,
                    required_n = n, label = "X")
    x_cols <- colnames(X)
  }
  X <- as.matrix(X)
  nS <- ncol(X)

  # ---- Z (scalar) ----
  if (!is.null(z)) {
    Z_mat <- .bind_cols(z, data = data, env = .env,
                        required_n = n, label = "Z")
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
    if (length(s_grid) != nS) {
      stop("`s_grid` must have length = ncol(X) = ", nS, ".")
    }
    s <- as.numeric(s_grid)
  }

  # ---- quadrature weights (trapezoidal rule) ----
  d <- diff(s)
  w <- c(d[1] / 2,
         (d[-length(d)] + d[-1]) / 2,
         d[length(d)] / 2)  # length nS

  # L: n x nS matrix of weights; S: n x nS matrix of s-grid
  L <- matrix(rep(w, times = n), nrow = n, byrow = TRUE)
  S <- matrix(rep(s, times = n), nrow = n, byrow = TRUE)

  # ---- tidy k and basis ----
  if (length(k) == 1L) k <- rep(k, 2L)
  if (length(basis) == 1L) basis <- rep(basis, 2L)
  if (length(k) != 2L) stop("`k` must have length 1 or 2.")
  if (length(basis) != 2L) stop("`basis` must have length 1 or 2.")

  # ---- build RHS string: Z terms (by name) + tensor term ----
  tp_term <- "ti(S, X, by = L, bs = basis, k = k, mc = c(FALSE, TRUE))"

  if (nZ > 0) {
    z_term <- paste(Z_names, collapse = " + ")
    rhs <- paste(z_term, tp_term, sep = " + ")
  } else {
    rhs <- tp_term
  }

  if (family == "lognormal") {
    # censored log-time response for cnorm()
    logY <- cbind(log(Y), log(Y))
    logY[Delta == 0, 2] <- Inf

    form <- as.formula(paste("logY ~", rhs))

    if (!is.null(data)) {
      fit <- mgcv::gam(
        formula = form,
        family  = mgcv::cnorm(),
        method  = "REML",
        data    = data
      )
    } else {
      fit <- mgcv::gam(
        formula = form,
        family  = mgcv::cnorm(),
        method  = "REML"
      )
    }
  } else {
    stop("Unsupported `family`: ", family)
  }

  # ---- attach meta info ----
  fit$X_names     <- x_cols
  fit$Z_names     <- Z_names
  fit$s_grid      <- s
  fit$family_name <- family

  fit
}
