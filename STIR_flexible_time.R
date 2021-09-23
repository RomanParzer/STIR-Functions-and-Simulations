require(MASS) # for ginv

#' Sufficient Dimension Reduction for Longitudinal Data
#'
#' @param X array of dim \eqn{t \times p \times n} with each \eqn{t \times p} matrix  representing a
#'  observation (assumed to be centered over observations).
#' @param y vector of \eqn{n} elements as responses
#' @param t,p,d,k dimensions, \eqn{d \leq t. k \leq \min(p,dH)}.
#' @param time_pts either vector of t ordered timepoints or array 
#' of dim \eqn{t \times n}, with each col giving t ordered time points for individual i
#' @param S_t either \eqn{t \times d} matrix (indicated with ind_tpts = FALSE) or a list of length n with 
#' each entry being a \eqn{t \times d} matrix (indicated with ind_tpts = TRUE) for one individual. 
#' Each \eqn{t \times d} matrix should have centered columns and full rank d. Interpreted as d basis functions applied to t timepoints.
#' 
#' @returns a list with components
#'  B: \eqn{dH \times p} matrix of coefficients, where \eqn{H+1} is the number of levels of y, OLS based estimate
#'  X_fitted: \eqn{t \times p \times n} array of fitted values
#'  delta_p, delta_t, delta_ls: estimates of the covariance matrix 
#'  X_red: reduced predictors, array of dim \eqn{t k \times n} with each row representing one vectorized \eqn{t \times k} observation
#'  
STIR <- function(X, y, t, p, d=1, H = 1, cont_y = FALSE, k = NULL, ind_tpts = FALSE, S_t=scale(sapply(1:d,function(di) (1:t/t)^di),scale=FALSE),
                 calc_VC_vecB=FALSE, calc_covX=FALSE, comp_Xred = FALSE,comp_B_ML=FALSE) {
  # Check and transform parameters
  if (!is.array(X)) X <- as.array(X)
  
  n <- dim(X)[3] ### add that here, it was further down before
  
  if (cont_y == FALSE) {
    if (!is.factor(y)) y <- factor(y)
    Gy <- model.matrix(~y)[,-1]
  } else {
    Gy <- do.call(cbind, lapply(1:(H%/%2), function(s, z) {
      cbind(cos(s * z), sin(s * z))
    }, z = 2 * pi * y))
    if (H%%2 == 1) {
      Gy <- cbind(Gy,cos((H%/%2+1)*2*pi*y))
    }
  }
  Gy <- scale(Gy,center=TRUE,scale=FALSE) # center columns
  
  stopifnot(nrow(Gy) == n, dim(X)[1:2] == c(t,p), H == ncol(Gy), d <= t)
  if (!is.null(k)) stopifnot(k <= min(p,d*H))
  
  # X: t p n
  X_n <- aperm(X, c(3, 1, 2))
  Sig_p <- matrix(rowMeans(apply(X_n, 2, cov)), p, p)
  Sig_t <- matrix(rowMeans(apply(X_n, 3, cov)), t, t)
  dim(X_n) <- c(n*t, p)
  
  # find Bhat
  Z <- NULL
  if (ind_tpts) {
    myX <- aperm(X,c(1,3,2))
    dim(myX) <- c(t*n,p)
    Z <- do.call(rbind,lapply(as.list(1:n),function(i) kronecker(S_t[[i]],matrix(Gy[i,],nrow=1))))
    
    Bhat <- solve(crossprod(Z),crossprod(Z,myX))
    Xhat <- sapply(1:n,function(i) c(kronecker(S_t[[i]],matrix(Gy[i,],nrow=1))%*%Bhat) )
  } else {
    M_S <- tcrossprod(ginv(crossprod(S_t)),S_t)
    M_G <- solve(crossprod(Gy),t(Gy))
    
    Bhat <- as.matrix(kronecker(M_S,M_G)%*%X_n)
    Xhat <- sapply(1:n,function(i) c(kronecker(S_t,matrix(Gy[i,],nrow=1))%*%Bhat) )
  }
  
  dim(Xhat) <- c(t,p,n)
  U <- X-Xhat # U already has mean 0
  dim(U) <- c(t*p, n)
  delta_ls <- tcrossprod(U) / (n-p*min(d,t)*H)
  delta_ml <- delta_ls * (n-p*min(d,t)*H) / n
  
  dim(U) <- c(t,p,n)
  U_n <- aperm(U, c(3, 1, 2))
  delta_p <- matrix(rowMeans(apply(U_n, 2, cov)), p, p)
  delta_t <- matrix(rowMeans(apply(U_n, 3, cov)), t, t)
  
  sd_B <- VC_vecB <- X_red <- lsvs <- X_red_kron <- covX <- NULL
  if (calc_VC_vecB==TRUE) {
    # C <- kronecker(diag(1,p),kronecker(M_S,M_G))
    # VC_vecB <- C %*% tcrossprod(kronecker(delta_ml,diag(1,n)),C)
    #### faster implementation
    VC_vecB <- matrix(c(0),d*H*p,d*H*p)
    if (ind_tpts) {
      ZZinv <- solve(crossprod(Z))
      for (r in 1:p) {
        for (s in 1:p) {
          tmp_help <- Reduce(`+`, lapply(as.list(1:n),function(i) {
            kronecker(crossprod(S_t[[i]],delta_ls[1:t+(r-1)*t,1:t+(s-1)*t])%*%S_t[[i]],crossprod(matrix(Gy[i,],nrow=1)))
          }))
          VC_vecB[1:(d*H)+(r-1)*d*H,1:(d*H)+(s-1)*d*H] <- ZZinv %*% tmp_help %*% ZZinv
        }
      }
    } else {
      GG_inv <- solve(crossprod(Gy))
      for (r in 1:p) {
        for (s in 1:p) {
          VC_vecB[1:(d*H)+(r-1)*d*H,1:(d*H)+(s-1)*d*H] <- kronecker(M_S%*%tcrossprod(delta_ls[1:t+(r-1)*t,1:t+(s-1)*t],M_S),GG_inv)
        }
      }
    }
    
    sd_B <- matrix(sqrt(diag( VC_vecB )),nrow=d*H,ncol=p)
  }
  if (calc_covX) {
    dim(X) <- c(t*p, n)
    covX <- cov(t(X))
  } 
  
  if (!is.null(k)) {
    dim(X) <- c(t*p, n)
    tmp <- svd(t(Bhat),nu=k,nv=0)
    lsvs <- tmp$u
    svd_d <- tmp$d
    
    if (comp_Xred) {
      X_red <- apply(X,2,function(xi) {
        matrix(solve(delta_ls,xi),t,p) %*% lsvs
      })
      dim(X_red) <- c(t*k,n)
    }
  }
  ## ML estimate with cov
  B_ML <- NULL
  if (comp_B_ML) {
    dim(X) <- c(t*p, n)
    del_inv <- solve(delta_ml)
    if (ind_tpts) {
      MDelM <- Reduce(`+`, lapply(as.list(1:n),function(i) {
        M_i <- kronecker(diag(1,p),kronecker(S_t[[i]],matrix(Gy[i,],nrow=1)))
        crossprod(M_i,del_inv%*%M_i)
      }))
      MDelX <- Reduce(`+`, lapply(as.list(1:n),function(i) {
        M_i <- kronecker(diag(1,p),kronecker(S_t[[i]],matrix(Gy[i,],nrow=1)))
        crossprod(M_i,del_inv%*%X[,i])
      }))
      B_ML <- matrix(solve(MDelM, MDelX),d*H,p)
    } else {
      MDelM <- Reduce(`+`, lapply(as.list(1:n),function(i) {
        M_i <- kronecker(diag(1,p),kronecker(S_t,matrix(Gy[i,],nrow=1)))
        crossprod(M_i,del_inv%*%M_i)
      }))
      MDelX <- Reduce(`+`, lapply(as.list(1:n),function(i) {
        M_i <- kronecker(diag(1,p),kronecker(S_t,matrix(Gy[i,],nrow=1)))
        crossprod(M_i,del_inv%*%X[,i])
      }))
      B_ML <- matrix(solve(MDelM, MDelX),d*H,p)
    }
  }
  
  result <- list(B = Bhat, sd_B = sd_B, VC_vecB=VC_vecB, X_fitted = Xhat, B_ML = B_ML,
                 delta_p = delta_p, delta_t = delta_t, delta_ls = delta_ls,delta_ml=delta_ml,
                 X_red = X_red, k=k, covX=covX, lsvs=lsvs, svd_d=svd_d,
                 cont_y=cont_y, S_t=S_t, Sig_p=Sig_p, Sig_t=Sig_t, ind_tpts=ind_tpts, Z=Z, Gy=Gy)
  return(result)
}


