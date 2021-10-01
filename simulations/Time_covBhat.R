######################################################################
######################################################################
### computing times Cov Bhat same
source('../functions/STIR.R')
source('../functions/multi_assign.R')
source('../functions/lsir_adapted.R')


#require(dplyr)
#require(stats)
require(MASS) # for ginv and mvrnorm
require(Matrix)
require(mvtnorm)
require(ROCR)
require(mgcv)
require(microbenchmark)

generateDataSTIR <- function(n,B,t,p,rho,H=1L,d=1L,what_time = "1:t/t", basis_S = c("Fourier","Polynomial"),
                             center_S = TRUE, ind_tpts = FALSE, delta_y = FALSE) {
  y <- rnorm(n,0,0.1)
  if (ind_tpts) {
    time_pts <- sapply(1:n,function(i) {
      sort(eval(str2expression(what_time)))
    })
  } else {
    time_pts <- sort(eval(str2expression(what_time)))
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
  
  Gy <- do.call(cbind, lapply(1:(H%/%2), function(s, z) {
    cbind(cos(s * z), sin(s * z))
  }, z = 2 * pi * y))
  if (H%%2 == 1) {
    Gy <- cbind(Gy,cos((H%/%2+1)*2*pi*y))
  }
  Gy <- scale(Gy,scale=FALSE)
  
  delta <- rho^abs(outer(1:(t*p), 1:(t*p), FUN = `-`))
  
  X <- sapply(1:n,function(i) {
    if (ind_tpts) {
      mu <- kronecker(S_t[((i-1)*t+1):(i*t),],matrix(Gy[i,],nrow = 1))%*%B
    } else {
      mu <- kronecker(S_t,matrix(Gy[i,],nrow = 1))%*%B
    }
    if (delta_y) {
      myscale <- exp(-min(sqrt(10)*abs(y[i]),2))
    } else {
      myscale <- 1
    }
    mvrnorm(1, c(mu), myscale*delta)
  })
  dim(X) <- c(t,p,n)
  X <- sweep(X,1:2,apply(X,1:2,mean))
  return(list(X=X, y=y, B=B, S=t(S_t), Gy=Gy, delta=delta, time_pts = time_pts))
}

myB12 <- matrix(c(rep(c(1,0.5),2),rep(c(0.7,0.3),2),rep(0.1,4)),4,3)

#                       n,         B, t, p, H, d, rho, k,   list(r_lsir,m_lsir),  time, ind_tpts, delta_y
#              --------------------------------------------------------
param <- list( 500,myB12,10, 3, 2, 2,  0.8, 3,list(c(1,1),c(3,1),c(3,2),c(7,2)),"exp(1:t/6)/exp(t/6)", FALSE, TRUE)
c(n,B,t,p,H,d,rho,k,list_rm_lsir,what_time,ind_tpts,delta_y) %<-% param

set.seed(123456)
datas <- generateDataSTIR(n,B,t,p,rho,basis_S = "Polynomial",H=H,d=d,what_time=what_time,ind_tpts=ind_tpts, delta_y=delta_y)
X <- datas$X  
y <- datas$y
# Check and transform parameters
if (!is.array(X)) X <- as.array(X)

n <- dim(X)[3] ### add that here, it was further down before

time_pts <- datas$time_pts
time_pts <- as.matrix(time_pts)
dim_tpts <- dim(time_pts)[2]
stopifnot(dim_tpts == 1 | dim_tpts==n)
if (dim_tpts==1) {
  ind_tpts <- FALSE
} else {
  ind_tpts <- TRUE
}

basis_type <- "Poly"
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
S_t <- scale(S_t,center = TRUE,scale=FALSE) # center columns (over timepoints)
cont_y <- TRUE
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

VC1 <- function() {
  C <- kronecker(diag(1,p),kronecker(M_S,M_G))
  VC_vecB <- C %*% tcrossprod(kronecker(delta_ls,diag(1,n)),C)
}

#### faster implementation
VC2 <- function() {
  GG_inv <- solve(crossprod(Gy))
  VC_vecB <- matrix(c(0),d*H*p,d*H*p)
  for (r in 1:p) {
    for (s in 1:p) {
      VC_vecB[1:(d*H)+(r-1)*d*H,1:(d*H)+(s-1)*d*H] <- kronecker(M_S%*%tcrossprod(delta_ls[1:t+(r-1)*t,1:t+(s-1)*t],M_S),GG_inv)
    }
  }
}

timesn500 <- microbenchmark(VC1(),VC2(),times = 10,unit = "s")

######### n 1000

myB12 <- matrix(c(rep(c(1,0.5),2),rep(c(0.7,0.3),2),rep(0.1,4)),4,3)

#                       n,         B, t, p, H, d, rho, k,   list(r_lsir,m_lsir),  time, ind_tpts, delta_y
#              --------------------------------------------------------
param <- list( 1000,myB12,10, 3, 2, 2,  0.8, 3,list(c(1,1),c(3,1),c(3,2),c(7,2)),"exp(1:t/6)/exp(t/6)", FALSE, TRUE)
c(n,B,t,p,H,d,rho,k,list_rm_lsir,what_time,ind_tpts,delta_y) %<-% param

set.seed(123456)
datas <- generateDataSTIR(n,B,t,p,rho,basis_S = "Polynomial",H=H,d=d,what_time=what_time,ind_tpts=ind_tpts, delta_y=delta_y)
X <- datas$X  
y <- datas$y
# Check and transform parameters
if (!is.array(X)) X <- as.array(X)

n <- dim(X)[3] ### add that here, it was further down before

time_pts <- datas$time_pts
time_pts <- as.matrix(time_pts)
dim_tpts <- dim(time_pts)[2]
stopifnot(dim_tpts == 1 | dim_tpts==n)
if (dim_tpts==1) {
  ind_tpts <- FALSE
} else {
  ind_tpts <- TRUE
}

basis_type <- "Poly"
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
S_t <- scale(S_t,center = TRUE,scale=FALSE) # center columns (over timepoints)
cont_y <- TRUE
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

VC1 <- function() {
  C <- kronecker(diag(1,p),kronecker(M_S,M_G))
  VC_vecB <- C %*% tcrossprod(kronecker(delta_ls,diag(1,n)),C)
}

#### faster implementation
VC2 <- function() {
  GG_inv <- solve(crossprod(Gy))
  VC_vecB <- matrix(c(0),d*H*p,d*H*p)
  for (r in 1:p) {
    for (s in 1:p) {
      VC_vecB[1:(d*H)+(r-1)*d*H,1:(d*H)+(s-1)*d*H] <- kronecker(M_S%*%tcrossprod(delta_ls[1:t+(r-1)*t,1:t+(s-1)*t],M_S),GG_inv)
    }
  }
}

timesn1000 <- microbenchmark(VC1(),VC2(),times = 10,unit = "s")

######### n 100

myB12 <- matrix(c(rep(c(1,0.5),2),rep(c(0.7,0.3),2),rep(0.1,4)),4,3)

#                       n,         B, t, p, H, d, rho, k,   list(r_lsir,m_lsir),  time, ind_tpts, delta_y
#              --------------------------------------------------------
param <- list( 100,myB12,10, 3, 2, 2,  0.8, 3,list(c(1,1),c(3,1),c(3,2),c(7,2)),"exp(1:t/6)/exp(t/6)", FALSE, TRUE)
c(n,B,t,p,H,d,rho,k,list_rm_lsir,what_time,ind_tpts,delta_y) %<-% param

set.seed(123456)
datas <- generateDataSTIR(n,B,t,p,rho,basis_S = "Polynomial",H=H,d=d,what_time=what_time,ind_tpts=ind_tpts, delta_y=delta_y)
X <- datas$X  
y <- datas$y
# Check and transform parameters
if (!is.array(X)) X <- as.array(X)

n <- dim(X)[3] ### add that here, it was further down before

time_pts <- datas$time_pts
time_pts <- as.matrix(time_pts)
dim_tpts <- dim(time_pts)[2]
stopifnot(dim_tpts == 1 | dim_tpts==n)
if (dim_tpts==1) {
  ind_tpts <- FALSE
} else {
  ind_tpts <- TRUE
}

basis_type <- "Poly"
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
S_t <- scale(S_t,center = TRUE,scale=FALSE) # center columns (over timepoints)
cont_y <- TRUE
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

VC1 <- function() {
  C <- kronecker(diag(1,p),kronecker(M_S,M_G))
  VC_vecB <- C %*% tcrossprod(kronecker(delta_ls,diag(1,n)),C)
}

#### faster implementation
VC2 <- function() {
  GG_inv <- solve(crossprod(Gy))
  VC_vecB <- matrix(c(0),d*H*p,d*H*p)
  for (r in 1:p) {
    for (s in 1:p) {
      VC_vecB[1:(d*H)+(r-1)*d*H,1:(d*H)+(s-1)*d*H] <- kronecker(M_S%*%tcrossprod(delta_ls[1:t+(r-1)*t,1:t+(s-1)*t],M_S),GG_inv)
    }
  }
}

timesn100 <- microbenchmark(VC1(),VC2(),times = 10,unit = "s")


timesn100
timesn500
timesn1000

