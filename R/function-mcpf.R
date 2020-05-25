#' Calculation of Pseudo F Statistic Matrix.
#'
#' \code{calct1t2} calculates the pseudo F statistics of the divided similarity matrixes
#' at given dividing point.
#'
#' @param Y the projection of pseudo outcome data.
#' @param kcut a number which indicates that the dividing points are
#'  1:kcut.
#' @param H1 the difference between two projection matrixes.
#' @param H2 the difference between the identity matirx and a projection
#'  matirx.
#'
#' @return a pseudo F statistic matrix.

calct1t2 <- function(Y, kcut, H1, H2){

  k = dim(Y)[2]

  tMat = matrix(NA, kcut,2)

  for (i in 1:kcut){

    G1 = as.matrix(Y[,1:i])%*%as.matrix(t(Y[,1:i]))

    G2 = as.matrix(Y[,(i+1):k])%*%as.matrix(t(Y[,(i+1):k]))

    tMat[i,1] = sum(H1*G1)/sum(H2*G1)

    tMat[i,2] = sum(H1*G2)/sum(H2*G2)

  }

  return(tMat)

}


#' Significance Calculation of MCPF in Distance Regression Model
#'
#' \code{mcpf} calculates the p value of MCPF in distance regression model.
#'
#' @param simi.mat a similarity matrix.
#' @param dis.mat a distance matrix. Only one of these two matrixes needs
#'  to be given. If dis.mat is given, it is then transformed to the
#'  similarity matrix.
#' @param null.space an index vector containing the index of ${X}_1$,
#'  the covariate matrix to be adjusted for.
#' @param x.mat a predictor matrix which contains the covariates and
#'  the predictors of interest.
#' @param n.monterCarlo the numbers of permutation.
#'
#' @return the p value of MCPF
#'
#' @importFrom MASS ginv
#'
#' @export
#'
#' @examples
#' library(MASS)
#' n = 50
#' p = 100
#' k = 15
#' sigmax = diag(rep(0.5,k)) + matrix(0.5,k,k)
#' sigmay = diag(rep(1,p))
#' for(i in 1:p){
#'   for(j in 1:p){
#'     sigmay[i,j] = 0.5^abs(i-j)
#'   }
#' }
#' r1 = 0.05
#' beta0 = r1*matrix(rbinom(k*p,1,0.9), k, p)
#' x = mvrnorm(n, rep(0,k), sigmax)
#' y = x%*%beta0 + mvrnorm(n, rep(0,p), sigmay)
#' Ky = calcKerMat_Cen(y)
#' mcpf(Ky,null.space = 1:5,x.mat = x)

mcpf <- function(simi.mat, dis.mat = NULL, null.space,
                 x.mat, n.monterCarlo = 1000){

  if(is.null(dis.mat) == FALSE){

    simi.mat = -0.5*dis.mat^2

  }

  x1 = x.mat[,null.space]

  if(length(null.space)==1){

    x1 = matrix(x1, ncol=1)

  }


  Hx = x.mat %*% ginv(t(x.mat)%*%x.mat) %*% t(x.mat)

  Hx1 = x1 %*% ginv(t(x1)%*%x1) %*% t(x1)

  n = nrow(simi.mat)

  In = diag(rep(1,n))

  H1 = Hx - Hx1

  H2 = In - Hx

  H3 = In - Hx1

  eigG = eigen(simi.mat)

  Gval = eigG$values

  Gvec = eigG$vectors

  k = sum(Gval>0.001)

  GY = Gvec[,1:k] %*% diag(sqrt(Gval[1:k]))

  QY = H3 %*% GY


  kcut = k-1

  difGval = diff(Gval)

  for(vv in 1:(k-2)){

    if(difGval[vv]>difGval[vv+1]){

      kcut = vv

      break

    }

  }

  B2 = n.monterCarlo

  T0 = calct1t2(QY, kcut, H1, H2)

  mat1 = matrix(NA, kcut, B2)

  mat2 = matrix(NA, kcut, B2)

  U = 1:n

  for(j in 1:B2){

    newQY = QY[sample(U, replace = FALSE),]

    newT = calct1t2(newQY,kcut,H1,H2)

    mat1[,j] = newT[,1]

    mat2[,j] = newT[,2]

  }

  rank1 = matrix(NA, kcut, B2)

  rank2 = matrix(NA, kcut, B2)

  combvec = rep(NA, kcut)

  combmat = matrix(NA, kcut, B2)

  for(l in 1:kcut){

    p1 = (sum(mat1[l,] > T0[l,1])+1)/(B2+1)

    p2 = (sum(mat2[l,] > T0[l,2])+1)/(B2+1)

    combvec[l] = -2*log(p1) - 2*log(p2)

    rank1[l,] = (B2 -rank(mat1[l,]) + 1)/B2

    rank2[l,] = (B2 -rank(mat2[l,]) + 1)/B2

    for(w in 1:B2){

      combmat[l,w] = -2*log(rank1[l,w]) - 2*log(rank2[l,w])

    }

  }

  mcpf = max(combvec)

  maxvec = apply(combmat,2,max)

  pvalue = sum(maxvec > mcpf)/B2

  return(pvalue)

}
