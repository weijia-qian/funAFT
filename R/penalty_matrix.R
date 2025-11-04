#' Function to construct the penalty matrix
#'
#' @param kp Number of basis coefficients
#' @param nS Density of the grid
#' @param a Parameter to ensure the penalty matrix is full rank (a < 0.01)
#'
#' @returns A kp x kp penalty matrix
#'
#' @examples
#' Pen_block <- penalty_matrix(kp = 30, nS = 100, a = 0.001)
#'
penalty_matrix <- function(kp, nS, a){
  D <- nS
  s <- seq(0, 1, length.out = nS)
  spline_basis <- splines::bs(s, df = kp, intercept = TRUE)
  diff2 <- matrix(rep(c(1, -2, 1, rep(0, D-2)), D - 2)[1:((D-2) * D)], D-2, D, byrow = TRUE)
  P2 <- t(spline_basis) %*% t(diff2) %*% diff2 %*% spline_basis # not full rank
  Pen <- a * diag(kp) + (1-a) * P2
  return(Pen)
}
