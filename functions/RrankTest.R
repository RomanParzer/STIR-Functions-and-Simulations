#' Test to Estimate Rank of Random Matrix (Bura and Yang 2011, non weighted version)
#'

#' @param Bhat estimated random matrix d1 x d2
#' @param V consistent estimate of the asymptotic covariance-matrix of vec(Bhat)
#' @param n number of observations
#' @param alpha test-level for sequential test rank = j vs > j  for j = 1, ... , min(d1,d2)
#'
#' @returns a list with components
#'  k: estimated rank
#'  p_vals: vector of p-values for dimensions up to k
#'
TestRMRank <- function(Bhat, V, n, alpha=0.05) {
  d1 <- nrow(Bhat)
  d2 <- ncol(Bhat)
  maxd <- min(d1,d2)
  k <- maxd
  
  mysvd <- svd(Bhat,nu=d1,nv=d2)
  pvals <- NULL
  for (j in 0:(maxd-1)) {
    U0 <- mysvd$u[,(j+1):d1]
    R0 <- mysvd$v[,(j+1):d2]
    K0 <- matrix(c(0), d1-j,d2-j)
    if (maxd - j > 1) {
      K0[1:(maxd - j),1:(maxd - j)] <- diag(mysvd$d[(j+1):maxd])
    } else {
      K0[1,1] <- mysvd$d[(j+1):maxd]
    }
    
    s <- min(rankMatrix(V),(d1-j)*(d2-j))
    krRU <- kronecker(R0,U0)
    Qhat <- crossprod(krRU,V)%*%krRU
    Lambda1 <- n*crossprod(c(K0),ginv(Qhat))%*%c(K0)
    
    # Qsvd <- svd(Qhat)
    # rkQ <- rankMatrix(Qhat)
    # Qsvd$d[1:rkQ] <- 1/Qsvd$d[1:rkQ]
    # Qsvd$d[-(1:rkQ)] <- 0
    # Qinv <- Qsvd$u%*%tcrossprod(diag(Qsvd$d),Qsvd$v)
    # Lambda2 <- n*crossprod(c(K0),Qinv)%*%c(K0)
    
    pvals <- c(pvals,1 - pchisq(Lambda1,s))
    
    if (tail(pvals,1) > alpha) {
      k <- j
      break
    }
  }
  return(list(k=k,pvals=pvals))
}

pacman::p_load(MASS,Matrix)

# # small test
# n <- 100
# res <- STIR(X=array(rnorm(n*10*5),dim=c(5,10,n)), y=rbinom(n,1,1/2), t=5, p=10,basis_S = "Polynomial",H=1,d=2,k=2,time_pts = 1:5,
#                  calc_VC_vecB=TRUE, calc_covX=TRUE, comp_Xred = FALSE,comp_B_ML = TRUE)
# TestRMRank(res$B,res$VC_vecB*n, n, 0.05)

