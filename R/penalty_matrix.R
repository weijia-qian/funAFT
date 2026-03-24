#' Function to construct the penalty matrix
#'
#' @param kp Number of basis coefficients
#' @param nS Density of the grid
#' @param a Parameter to ensure the penalty matrix is full rank (a < 0.01)
#' @param basis Character string specifying the basis type (e.g., "bs", "ns", "cc", "cr"). Defaults to "bs".
#'
#' @returns A kp x kp penalty matrix
#'
#'
penalty_matrix <- function(kp, nS, a, basis = "bs") {
  D <- nS
  s <- seq(0, 1, length.out = nS)

  if (basis %in% c("bs", "ns")) {

    if (basis == "bs") {
      spline_basis <- splines::bs(s, df = kp, intercept = TRUE)
    } else {
      spline_basis <- splines::ns(s, df = kp, intercept = TRUE)
    }

    # Construct the second-order difference matrix
    diff2 <- diff(diag(D), differences = 2)

    P2 <- t(spline_basis) %*% t(diff2) %*% diff2 %*% spline_basis # not full rank
    Pen <- a * diag(kp) + (1 - a) * P2

    return(Pen)

  }
}
