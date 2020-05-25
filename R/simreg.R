#' Significance Calculation of the Pseudo F Statistic in Distance
#' Regression Model.
#'
#' \code{simreg} calculates the p value of the pseudo F statistic in distance
#' regression model.
#'
#' @param simi.mat a similarity matrix.
#' @param dis.mat a distance matrix. Only one of these two matrixes needs
#'  to be given. If dis.mat is given, it is then transformed to the
#'  similarity matrix.
#' @param null.space an index vector containing the index of ${X}_1$,
#'  the covariate matrix to be adjusted for.
#' @param x.mat a predictor matrix which contains the covariates and
#'  the predictors of interest.
#' @param permute logical.Should sampling with replacement.
#' @param n.monterCarlo the numbers of permutation.
#' @param seed the seed for permutation.
#'
#' @return the p value of pseudo F test statistic.

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
#' simreg(Ky,null.space = 1:5,x.mat = x)

simreg <- function(simi.mat, dis.mat = NULL, null.space, x.mat,
                   permute=TRUE, n.monterCarlo=1000, seed=NULL){

  if(is.null(dis.mat) == FALSE){
    simi.mat = -0.5*dis.mat^2
  }

  if(!is.null(seed)){
    set.seed(seed)
  }

  x1 = x.mat[,null.space]

  if(length(null.space)==1){
    x1 = matrix(x1, ncol=1)
  }


  x.hat = x.mat %*% ginv(t(x.mat)%*%x.mat) %*% t(x.mat)

  x1.hat = x1 %*% ginv(t(x1)%*%x1) %*% t(x1)

  n = nrow(simi.mat)

  I.n = diag(n)

  cent = I.n - matrix(1,nrow=n,ncol=n)/n

  i.x  = I.n - x.hat

  i.x1 = I.n - x1.hat

  Q = i.x1 %*% cent %*% simi.mat %*% cent %*% i.x1

  alter.hat = x.hat - x1.hat

  F.obs = sum(alter.hat*Q) / sum(i.x*Q)

  U = 1:n

  F.star = rep(NA, n.monterCarlo)

  for (i in 1:n.monterCarlo){
    id.sam = sample(U, replace=!permute)

    Q.star = Q[id.sam, id.sam]

    F.star[i] = sum(alter.hat*Q.star) / sum(i.x*Q.star)

  }

  pvalue = sum(F.star > F.obs)/n.monterCarlo
  return(pvalue)
}
