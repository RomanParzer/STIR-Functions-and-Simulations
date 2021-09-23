require(MASS) # for ginv

#' Sufficient Dimension Reduction for Longitudinal Data
#'
#' @param X array of dim \eqn{t \times p \times n} with each \eqn{t \times p} matrix  representing a
#'  observation (assumed to be centered over observations).
#' @param y vector of \eqn{n} elements as responses
#' @param t,p,d,k dimensions, \eqn{d \leq t. k \leq \min(p,dH)}.
#' @param time_pts either vector of t ordered timepoints or array 
#' of dim \eqn{t \times n}, with each col giving t ordered time points for individual i
#'
#' @returns a list with components
#'  B: \eqn{dH \times p} matrix of coefficients, where \eqn{H+1} is the number of levels of y, OLS based estimate
#'  X_fitted: \eqn{t \times p \times n} array of fitted values
#'  delta_p, delta_t, delta_ls: estimates of the covariance matrix 
#'  X_red: reduced predictors, array of dim \eqn{t k \times n} with each row representing one vectorized \eqn{t \times k} observation
#'  
STIR <- function(X, y, t, p, d=1, H = 1, cont_y = FALSE, k = NULL, time_pts = 1:t/t, basis_S = c("Polynomial","Fourier"),
                 center_S = TRUE,calc_VC_vecB=FALSE, calc_covX=FALSE, comp_Xred = FALSE,comp_B_ML=FALSE) {
  # Check and transform parameters
  if (!is.array(X)) X <- as.array(X)
  
  n <- dim(X)[3] ### add that here, it was further down before
  
  time_pts <- as.matrix(time_pts)
  dim_tpts <- dim(time_pts)[2]
  stopifnot(dim_tpts == 1 | dim_tpts==n)
  if (dim_tpts==1) {
    ind_tpts <- FALSE
  } else {
    ind_tpts <- TRUE
  }
  
  basis_type <- match.arg(basis_S)
  if (basis_type == "Fourier") {
    if (d>1) {
      S_t <- do.call(cbind, lapply(1:(d%/%2), function(s, z) {
        cbind(cos(s * z), sin(s * z))
      }, z = 2 * pi * c(time_pts)))
    } else {
      S_t <- NULL
    }
    if (d%%2 == 1) {
      S_t <- cbind(S_t,cos((d%/%2+1)*2*pi*c(time_pts)))
    }
  } else {
    S_t <- sapply(1:(d),function(di) c(time_pts)^di) # polynomial basis on times in [0,1]
  }
  S_t <- scale(S_t,center = center_S,scale=FALSE) # center columns (over timepoints)
  
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
  Gy <- scale(Gy,scale=FALSE)
  
  stopifnot(nrow(Gy) == n, dim(X)[1:2] == c(t,p), dim(time_pts)[1]==t, H == ncol(Gy), d <= t)
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
    Z <- do.call(rbind,lapply(as.list(1:n),function(i) kronecker(S_t[((i-1)*t+1):(i*t),],matrix(Gy[i,],nrow=1))))
    
    Bhat <- solve(crossprod(Z),crossprod(Z,myX))
    Xhat <- sapply(1:n,function(i) c(kronecker(S_t[((i-1)*t+1):(i*t),],matrix(Gy[i,],nrow=1))%*%Bhat) )
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
  
  sd_B <- VC_vecB <- X_red <- lsvs <- svd_d <- X_red_kron <- covX <- NULL
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
            kronecker(crossprod(S_t[((i-1)*t+1):(i*t),],delta_ml[1:t+(r-1)*t,1:t+(s-1)*t])%*%S_t[((i-1)*t+1):(i*t),],crossprod(matrix(Gy[i,],nrow=1)))
            }))
          VC_vecB[1:(d*H)+(r-1)*d*H,1:(d*H)+(s-1)*d*H] <- ZZinv %*% tmp_help %*% ZZinv
        }
      }
    } else {
      GG_inv <- solve(crossprod(Gy))
      for (r in 1:p) {
        for (s in 1:p) {
          VC_vecB[1:(d*H)+(r-1)*d*H,1:(d*H)+(s-1)*d*H] <- kronecker(M_S%*%tcrossprod(delta_ml[1:t+(r-1)*t,1:t+(s-1)*t],M_S),GG_inv)
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
        matrix(solve(delta_ml,xi),t,p) %*% lsvs
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
        M_i <- kronecker(diag(1,p),kronecker(S_t[((i-1)*t+1):(i*t),],matrix(Gy[i,],nrow=1)))
        crossprod(M_i,del_inv%*%M_i)
      }))
      MDelX <- Reduce(`+`, lapply(as.list(1:n),function(i) {
        M_i <- kronecker(diag(1,p),kronecker(S_t[((i-1)*t+1):(i*t),],matrix(Gy[i,],nrow=1)))
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

#' compute Grassman-metric distance between B_hat and true B (of possibly different dimensions)
#'
#' @param B_hat coefficient matrix of dim \eqn{tp \times kr} (LSIR) or \eqn{dH \times p} (STIR)
#' @param B true coefficient matrix of \eqn{dH \times p} (STIR)
#' 
#' see Ye and Lim (2016): Schubert varieties and distances between subspaces of different dimensions
#' 
grass_metric <- function(B_hat,B) {
  # get ONBs of the spans
  svd1 <- svd(B_hat)
  U1 <- matrix(svd1$u[,svd1$d>1e-13])
  
  svd2 <- svd(B)
  U2 <- matrix(svd2$u[,svd2$d>1e-13])
  
  stopifnot(nrow(U1) >= ncol(U1),
            nrow(U2) >= ncol(U2))
  
  if (nrow(U1) > nrow(U2)) { # make U2 the "bigger" one
    tmp <- U1
    U1 <- U2
    U2 <- tmp
  }
  
  stopifnot(ncol(U1)<= ncol(U2),
            ncol(U2)-ncol(U1) <=  nrow(U2)-nrow(U1))
  
  U1_emb <- rbind( U1, matrix( rep(0,ncol(U1)*(nrow(U2)-nrow(U1))) ,ncol = ncol(U1) ) )
  
  d <- svd(crossprod(U2,U1_emb), nu=0, nv=0)$d

  return(list(dist=sqrt( pi^2/4 *(ncol(U2)-ncol(U1)) + sum(acos(pmin(pmax(d,-1.0),1.0))^2) ),
              d=d))
}


