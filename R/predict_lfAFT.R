#' Predict from a Linear Functional AFT (lfAFT) Model
#'
#' Compute predictions for new observations from a fitted linear functional
#' AFT model, using the estimated functional
#' coefficient \eqn{\hat{\beta}(s)} and scalar covariates.
#'
#' @param fit A fitted lfAFT model object returned by \code{\link{fit_lfAFT}}.
#' @param newdata A \code{data.frame} (or similar) containing new covariate
#'   information. It must contain:
#'   \itemize{
#'     \item Columns whose names match \code{fit$X_names}, representing the
#'           functional covariate grid used in fitting.
#'     \item If the original model included scalar covariates, columns whose
#'           names match \code{fit$Z_names}.
#'   }
#' @param type Character; one of \code{c("link", "response", "survival")}.
#'   \describe{
#'     \item{\code{"link"}}{Return the linear predictor \eqn{\mu_i}.}
#'     \item{\code{"response"}}{Return \code{exp(mu_i)}, which corresponds to the
#'       predicted median survival time under the lognormal and loglogistic AFT
#'       parameterizations used in \code{fit_lfAFT}.}
#'     \item{\code{"survival"}}{Return the predicted survival probabilities
#'       \eqn{S(t \mid X_i, Z_i)} at user-specified time points \code{times}.}
#'   }
#' @param times Optional numeric vector of time points (on the original time
#'   scale) at which to evaluate the survival function. Required when
#'   \code{type = "survival"}. All values must be strictly positive.
#'
#' @return
#' A numeric object whose form depends on \code{type}:
#' \itemize{
#'   \item \code{type = "link"}: numeric vector of linear predictors
#'         \eqn{\hat{\mu}_i}.
#'   \item \code{type = "response"}: numeric vector of predicted median survival
#'         times \eqn{\exp(\hat{\mu}_i)}.
#'   \item \code{type = "survival"}: numeric matrix of survival probabilities
#'         with one row per observation and one column per time point in
#'         \code{times}. Row names are taken from \code{newdata} (if present);
#'         column names are set to \code{times}.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Suppose 'fit' is an lfAFT model fitted via fit_lfAFT()
#' # newdata must contain columns matching fit$X_names and fit$Z_names (if any):
#'
#' # 1) Linear predictor
#' mu_new <- predict_lfAFT(fit, newdata, type = "link")
#'
#' # 2) Predicted median survival times
#' t_med_new <- predict_lfAFT(fit, newdata, type = "response")
#'
#' # 3) Survival probabilities at t = 1, 2, 5
#' S_new <- predict_lfAFT(fit, newdata,
#'                        type  = "survival",
#'                        times = c(1, 2, 5))
#' head(S_new)
#' }

predict_lfAFT <- function(fit,
                          newdata,
                          type = c("link", "response", "survival"),
                          times = NULL) {
  type <- match.arg(type)

  # ---- basic checks on fit ----
  if (is.null(fit$beta0_hat) || is.null(fit$betaX_hat)) {
    stop("`fit` does not appear to be a valid lfAFT fit (missing beta0_hat or betaX_hat).")
  }
  if (is.null(fit$X_names)) {
    stop("`fit$X_names` is NULL; cannot locate functional covariates in `newdata`.")
  }

  # ---- extract new functional covariates by name ----
  X_names <- fit$X_names
  if (!all(X_names %in% names(newdata))) {
    missing <- X_names[!X_names %in% names(newdata)]
    stop("The following functional covariate columns are missing in `newdata`: ",
         paste(missing, collapse = ", "))
  }
  X_new <- as.matrix(newdata[, X_names, drop = FALSE])

  beta0 <- fit$beta0_hat
  betaX <- fit$betaX_hat     # estimated beta(s) on s_grid
  betaZ <- fit$betaZ_hat     # may be NULL if no scalar covariates

  # grid & dimensional consistency
  s <- fit$s_grid
  if (is.null(s)) stop("`fit$s_grid` is missing; cannot compute quadrature weights for prediction.")
  nS_fit <- length(s)
  nS_new <- ncol(X_new)

  if (nS_fit != nS_new) {
    warning(sprintf(
      "Number of grid points in fit (length(s_grid) = %d) does not match ncol(X_new) = %d.",
      nS_fit, nS_new
    ))
  }

  if (length(betaX) != nS_fit) {
    stop(sprintf(
      "Length of betaX_hat (%d) does not match length of s_grid (%d).",
      length(betaX), nS_fit
    ))
  }

  # ---- trapezoid weights from s_grid (same as used in fitting) ----
  d <- diff(s)
  w <- c(d[1] / 2,
         (d[-length(d)] + d[-1]) / 2,
         d[length(d)] / 2)

  if (length(w) != nS_fit) {
    stop("Length of trapezoid weights does not match length of s_grid.")
  }

  # precompute weights for functional term: sum_j w_j * X_ij * betaX_j
  beta_w <- betaX * w  # elementwise

  # ---- scalar covariates, if present ----
  Z_names <- fit$Z_names
  if (!is.null(betaZ)) {
    if (is.null(Z_names)) {
      stop("Scalar covariates were estimated, but `fit$Z_names` is NULL.")
    }
    if (!all(Z_names %in% names(newdata))) {
      missing <- Z_names[!Z_names %in% names(newdata)]
      stop("The following scalar covariates are missing in `newdata`: ",
           paste(missing, collapse = ", "))
    }
    Z_mat <- as.matrix(newdata[, Z_names, drop = FALSE])
    functional_part <- as.vector(X_new %*% beta_w)
    mu_new <- as.vector(beta0 + Z_mat %*% betaZ + functional_part)
  } else {
    # no scalar covariates
    functional_part <- as.vector(X_new %*% beta_w)
    mu_new <- as.vector(beta0 + functional_part)
  }

  # ---- output by type ----
  if (type == "link") {
    pred <- mu_new

  } else if (type == "response") {
    # For both lognormal and loglogistic AFT parameterizations, exp(mu)
    # corresponds to the median survival time on the original time scale.
    pred <- exp(mu_new)

  } else {  # type == "survival"
    if (is.null(times)) {
      stop("`times` must be supplied when type = \"survival\".")
    }
    if (!is.numeric(times) || any(times <= 0)) {
      stop("`times` must be a numeric vector of strictly positive values.")
    }

    if (is.null(fit$family) || is.null(fit$sigma_hat)) {
      stop("`fit$family` and `fit$sigma_hat` must be present to compute survival probabilities.")
    }

    family <- fit$family
    sigma  <- fit$sigma_hat

    n_obs   <- length(mu_new)
    n_times <- length(times)
    log_t   <- log(times)

    if (family == "lognormal") {
      z <- outer(mu_new, log_t, function(m, lt) (lt - m) / sigma)
      S <- 1 - pnorm(z)
    } else if (family == "loglogistic") {
      z <- outer(mu_new, log_t, function(m, lt) (lt - m) / sigma)
      S <- 1 / (1 + exp(z))
    } else {
      stop("Unknown family in `fit$family`: ", family)
    }

    rownames(S) <- rownames(newdata)
    colnames(S) <- as.character(times)
    pred <- S
  }

  # preserve row names for vector outputs
  if (type != "survival" && !is.null(rownames(newdata))) {
    names(pred) <- rownames(newdata)
  }

  pred
}
