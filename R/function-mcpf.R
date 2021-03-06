#' Calculation of Pseudo F Statistic Matrix.
#'
#' \code{calct1t2} calculates the pseudo F statistics of the divided similarity matrixes
#' at given dividing point.
#'
#' @param Y the projection of pseudo outcome data.
#' @param point.seq a vector containing the dividing points.
#' @param H1 the difference between two projection matrixes.
#' @param H2 the difference between the identity matirx and a projection
#'  matirx.
#'
#' @return a length(point.seq)*2 matrix of pseudo F statistics.

calct1t2 <- function(Y, point.seq, H1, H2){

  k = dim(Y)[2]

  kcut = length(point.seq)

  tMat = matrix(NA, kcut,2)

  for (i in 1:kcut){

    ipoint = point.seq[i]

    G1 = as.matrix(Y[,1:ipoint])%*%as.matrix(t(Y[,1:ipoint]))

    G2 = as.matrix(Y[,(ipoint+1):k])%*%as.matrix(t(Y[,(ipoint+1):k]))

    tMat[i,1] = sum(H1*G1)/sum(H2*G1)

    tMat[i,2] = sum(H1*G2)/sum(H2*G2)

  }

  return(tMat)

}



#' Generate Fibonacci Sequence to be Used as Dividing Points
#'
#' \code{Fibonacci} generates a Fibonacci sequence to be used as dividing
#'  points in MCPF. So it deletes the first element 1 in the original
#'  Fibonacci sequence.
#'
#' @param n a number indicating the maximum of the generated sequence is
#'  less than n.
#'
#' @return a vector.
Fibonacci <- function(n){
  fab = NA
  fab[1] = 1
  fab[2] = 1
  if(n == 1|n == 2){
    return(fab[1])
  }else{
    i = 2
    while(fab[i] < n){
      i = i + 1
      fab[i] = fab[i-1] + fab[i-2]
    }
    return(fab[2:(i-1)])
  }
}




#' Significance Calculation of MCPF in Distance-based Regression Model
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
#' @param n.monterCarlo the number of permutation replicates.
#' @param rMethod the method used to generate the vector of dividing point.
#'  It should be one of "interpolation" and "fibonacci", where
#'  "interpolation" indicates that the vector is generated via interpolation
#'  strategy, while "fibonacci" indicates that it is generated through
#'  Fibonacci sequence.
#'
#' @return the p value of MCPF.
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
#' mcpf(Ky,null.space = 1:5,x.mat = x,rMethod = "interpolation")

mcpf <- function(simi.mat, dis.mat = NULL, null.space,
                 x.mat, n.monterCarlo = 1000,
                 rMethod = c("interpolation","fibonacci")){

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

  if(rMethod == "interpolation"){

    kcut = k-1

    difGval = diff(Gval)

    for(vv in 1:(k-2)){

      if(difGval[vv]>difGval[vv+1]){

        kcut = vv

        break

      }

    }

    seq0 = seq(1:kcut)

  }else if(rMethod == "fibonacci"){

    seq0 = Fibonacci(k)

    kcut = length(seq0)

  }else{

    return(message("'arg' should be one of \"interpolation\" and \"fibonacci\"."))


  }

  B2 = n.monterCarlo

  T0 = calct1t2(QY, seq0, H1, H2)

  mat1 = matrix(NA, kcut, B2)

  mat2 = matrix(NA, kcut, B2)

  U = 1:n

  for(j in 1:B2){

    newQY = QY[sample(U, replace = FALSE),]

    newT = calct1t2(newQY,seq0,H1,H2)

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

