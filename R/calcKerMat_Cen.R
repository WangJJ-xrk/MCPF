#' Similarity Matrix Calculation.
#'
#' \code{calcKerMat_Cen} calculates centered similarity matrix for a data
#' matrix.
#'
#' @param X the data matrix, each row represents an observation.
#'
#' @return the centered similarity matrix.
#'
#' @importFrom stats dist
#'
#' @export
#'
#' @examples
#' library(MASS)
#' n = 50
#' p = 100
#' k = 15
#' sigmax = diag(rep(0.5,k)) + matrix(0.5,k,k)
#' x = mvrnorm(n, rep(0,k), sigmax)
#' Kx = calcKerMat_Cen(x)
calcKerMat_Cen <- function(X){

  dis.X = as.matrix(dist(X))

  K.X = -0.5*dis.X^2

  n = dim(K.X)[1]

  H = diag(rep(1,n)) - matrix(1/n,n,n)

  cK.X = H %*% K.X %*% H

  return(cK.X)

}


