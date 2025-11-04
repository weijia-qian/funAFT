#' Optimize the Smoothing Parameter (Lambda) via Generalized Cross-Validation (GCV)
#'
#' Evaluates a user-specified grid of smoothing parameters (\eqn{\lambda})
#' for the linear functional accelerated failure time (lfAFT) model and
#' selects the value that minimizes the generalized cross-validation (GCV) score.
#'
#' @param lambda_grid Numeric vector of candidate smoothing parameters.
#'                    Must have length ≥ 1. Each value is evaluated separately
#'                    by fitting alfAFT model and computing its GCV.
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
#' @param k Integer; number of spline basis functions used to represent
#'          \eqn{\beta(s)} (effective degrees of freedom for the functional coefficient).
#' @param s_grid Optional numeric vector of the same length as the number of
#'               functional measurement points (\code{ncol(X)}). If \code{NULL},
#'               assumes an equally spaced grid over \[0, 1\].
#' @param basis Character; one of \code{c("bs", "ns")}. Determines the spline
#'              basis type used for the functional coefficient. Defaults to
#'              \code{"bs"} (B-spline).
#' @param basis_args List of additional arguments passed to the spline basis
#'                   constructor (e.g., \code{degree = 3} for \code{bs},
#'                   or \code{intercept = FALSE}).
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{lambda}}{Optimal smoothing parameter minimizing GCV.}
#'   \item{\code{gcv}}{Numeric vector of GCV scores corresponding to
#'         each candidate \code{lambda_grid} value (may contain \code{NA}s).}
#'   \item{\code{grid}}{The input \code{lambda_grid}.}
#'   \item{\code{best_index}}{Integer index of the selected optimal \code{lambda}.}
#' }
#'
#' @seealso \code{\link{fit_lfAFT}} for model fitting with automatic
#'          lambda selection, and \code{\link{fit_lfAFT_single_lambda}}
#'          for fitting at a fixed \eqn{\lambda}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Define a candidate grid of smoothing parameters
#' lam_grid <- 10^seq(1, 4, length.out = 10)
#'
#' # Evaluate GCV for each lambda (example with synthetic data)
#' opt <- optimize_lambda(lambda_grid = lam_grid,
#'                        data = mydata,
#'                        y = "Y", delta = "Delta",
#'                        x = "^X", x_as_regex = TRUE,
#'                        family = "lognormal",
#'                        k = 10)
#'
#' plot(lam_grid, opt$gcv, type = "b", log = "x",
#'      xlab = "lambda", ylab = "GCV")
#' abline(v = opt$lambda, col = "red", lty = 2)
#' }

optimize_lambda <- function(lambda_grid,
                            data = NULL,
                            y,
                            delta,
                            x,
                            x_as_regex = FALSE,
                            z = NULL,
                            family = c("lognormal", "loglogistic"),
                            k = 20,
                            s_grid = NULL,
                            basis = c("bs", "ns"),
                            basis_args = list()) {

  family <- match.arg(family)
  basis  <- match.arg(basis)

  if (!is.numeric(lambda_grid) || length(lambda_grid) < 1L)
    stop("`lambda_grid` must be a non-empty numeric vector.")

  gcv_values <- sapply(lambda_grid, function(lam) {
    fit_i <- try(
      fit_lfAFT_single_lambda(
        data   = data,
        y      = y,
        delta  = delta,
        x      = x,
        x_as_regex = x_as_regex,
        z      = z,
        family = family,
        k      = k,
        lambda = lam,
        s_grid = s_grid,
        basis  = basis,
        basis_args = basis_args,
        se        = FALSE,
        bootstrap = FALSE
      ),
      silent = TRUE
    )

    if (inherits(fit_i, "try-error") ||
        is.null(fit_i$GCV) ||
        !is.finite(fit_i$GCV)) {
      return(NA_real_)
    }
    fit_i$GCV
  })

  finite_idx <- which(is.finite(gcv_values))
  if (length(finite_idx) == 0L) {
    stop("All GCV evaluations failed or returned non-finite values.")
  }

  best_idx <- finite_idx[which.min(gcv_values[finite_idx])]
  optimal_lambda <- lambda_grid[best_idx]

  names(gcv_values) <- format(lambda_grid, digits = 6)

  list(
    lambda     = optimal_lambda,
    gcv        = gcv_values,
    grid       = lambda_grid,
    best_index = best_idx
  )
}




