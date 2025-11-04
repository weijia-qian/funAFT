#' Predict from an Additive Functional AFT (afAFT) Model
#'
#' Compute predictions for new observations from a fitted additive
#' functional AFT model obtained via \code{\link{fit_afAFT}}.
#'
#' @param fit A fitted model object returned by \code{\link{fit_afAFT}}
#'   (i.e., a \code{\link[mgcv]{gam}} object with additional components).
#' @param newdata A \code{data.frame} (or compatible object) containing
#'   new covariate information. It must contain:
#'   \itemize{
#'     \item Columns whose names match \code{fit$X_names}, representing the
#'           functional covariate grid used in fitting.
#'     \item If scalar covariates were included in the fit,
#'           columns whose names match \code{fit$Z_names}.
#'   }
#' @param type Character; one of \code{c("link", "response", "survival")}.
#'   \describe{
#'     \item{\code{"link"}}{
#'       the linear predictor \eqn{\mu_i = E[\log T_i \mid X_i, Z_i]}. \cr
#'     }
#'     \item{\code{"response"}}{
#'        \eqn{\exp(\mu_i)}, the predicted median survival time on
#'       the original time scale. \cr
#'     }
#'     \item{\code{"survival"}}{
#'       returns predicted survival probabilities \eqn{S(t \mid X_i, Z_i)} at user-specified
#'       time points \code{times}.
#'     }
#'   }
#' @param times Optional numeric vector of time points (on the original time
#'   scale) at which to evaluate the survival function. Required when
#'   \code{type = "survival"}. All values must be strictly positive.
#'
#'
#' @return
#' A numeric object whose form depends on \code{type}:
#' \itemize{
#'   \item \code{type = "link"}:
#'     numeric vector of linear predictors \eqn{\hat{\mu}_i}.
#'   \item \code{type = "response"}:
#'     numeric vector of median survival times \eqn{\exp(\hat{\mu}_i)}.
#'   \item \code{type = "survival"}:
#'     numeric matrix of survival probabilities with one row per new
#'     observation and one column per time point in \code{times}. Row names
#'     are taken from \code{newdata} (if present) and column names are set
#'     to \code{as.character(times)}.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Suppose 'fit' is an afAFT model fitted via fit_afAFT()
#' # newdata must contain columns matching fit$X_names and fit$Z_names (if any)
#'
#' # 1) Linear predictor
#' mu_new <- predict_afAFT(fit, newdata, type = "link")
#'
#' # 2) Predicted median survival times
#' t_med_new <- predict_afAFT(fit, newdata, type = "response")
#'
#' # 3) Survival probabilities at t = 1, 2, 5
#' S_new <- predict_afAFT(fit, newdata,
#'                        type  = "survival",
#'                        times = c(1, 2, 5))
#' head(S_new)
#' }

predict_afAFT <- function(fit,
                          newdata,
                          type  = c("link", "response", "survival"),
                          times = NULL) {
  type <- match.arg(type)

  if (is.null(fit$X_names) || is.null(fit$s_grid)) {
    stop("`fit` does not appear to be a valid afAFT fit (missing X_names or s_grid).")
  }
  family_name <- fit$family_name

  # Coerce newdata to data.frame for safe subsetting/augmentation
  newdata <- as.data.frame(newdata)

  # ---- functional covariate: X_new, S, L ----
  X_names <- fit$X_names
  if (!all(X_names %in% names(newdata))) {
    missing <- X_names[!X_names %in% names(newdata)]
    stop("The following functional covariate columns are missing in `newdata`: ",
         paste(missing, collapse = ", "))
  }
  X_new <- as.matrix(newdata[, X_names, drop = FALSE])

  s   <- fit$s_grid
  n   <- nrow(X_new)
  nS  <- length(s)
  if (nS != ncol(X_new)) {
    warning(sprintf("Number of grid points in fit (length(s_grid) = %d) ",
                    nS),
            sprintf("does not match ncol(newdata[, X_names]) = %d.", ncol(X_new)))
  }

  # Trapezoid weights and S,L matrices
  d <- diff(s)
  w <- c(d[1] / 2,
         (d[-length(d)] + d[-1]) / 2,
         d[length(d)] / 2)

  if (length(w) != nS) {
    stop("Length of quadrature weights does not match length of s_grid.")
  }

  S_mat <- matrix(rep(s, times = n), nrow = n, byrow = TRUE)
  L_mat <- matrix(rep(w, times = n), nrow = n, byrow = TRUE)

  # add matrix covariates to newdata for predict.gam
  newdata$X <- X_new
  newdata$S <- S_mat
  newdata$L <- L_mat

  # ---- check scalar covariates ----
  Z_names <- fit$Z_names
  if (!is.null(Z_names)) {
    if (!all(Z_names %in% names(newdata))) {
      missing <- Z_names[!Z_names %in% names(newdata)]
      stop("The following scalar covariates are missing in `newdata`: ",
           paste(missing, collapse = ", "))
    }
  }

  # ---- get linear predictor from mgcv ----
  eta <- as.vector(mgcv::predict.gam(fit, newdata = newdata, type = "link"))

  # ---- handle type & family ----
  if (type == "link") {
    pred <- eta

  } else if (type == "response") {
      # median survival time
      pred <- exp(eta)

  } else {  # type == "survival"
    if (is.null(times)) {
      stop("`times` must be supplied when type = \"survival\".")
    }
    if (!is.numeric(times) || any(times <= 0)) {
      stop("`times` must be a numeric vector of strictly positive values.")
    }

    # extract sigma (scale) from cnorm family if available
    if (!is.null(fit$family$family)) {
      sigma <- as.numeric(gsub(".*\\(([^)]+)\\).*", "\\1", fit$family$family))
    } else {
      warning("Could not find `family$theta` in fit; using sigma = 1.")
      sigma <- 1
    }

    log_t <- log(times)
    z     <- outer(eta, log_t, function(m, lt) (lt - m) / sigma)
    S     <- 1 - pnorm(z)

    rownames(S) <- rownames(newdata)
    colnames(S) <- as.character(times)
    pred <- S
  }

  # preserve row names for vector outputs
  if (!is.matrix(pred) && !is.null(rownames(newdata))) {
    names(pred) <- rownames(newdata)
  }

  pred
}
