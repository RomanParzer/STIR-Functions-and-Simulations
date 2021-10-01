#source('./random.R')
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

## returns X: T x p x n
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

simulateSDR <- function(nsim, n, B, t, p, rho,H=1, d=1, k=1, npred = 100, what_time = "1:t/t",list_rm_lsir=list(c(1,1)),nlev_ydis = 10,
                        delta_y=FALSE, ind_tpts = FALSE) {
  res <- list(STIR=NULL,LSIR=NULL)
  B_est <- list(STIR=array(c(0),dim=c(nsim,d*H,p)),LSIR=array(c(0),dim=c(nsim,p*t,p*t)))
  B_ML_est <- list(STIR=array(c(0),dim=c(nsim,d*H,p)))
  delta_ls_est <- list(STIR=array(c(0),dim=c(nsim,p*t,p*t)))
  
  for (i in 1:nsim) {
    if ((10*i)%%nsim == 0) {
      cat(sprintf('Starting rep %d/%d\n', i, nsim))
    }
    
    set.seed(i+1000)
    datas <- generateDataSTIR(n,B,t,p,rho,basis_S = "Polynomial",H=H,d=d,what_time=what_time,ind_tpts=ind_tpts, delta_y=delta_y)
    if (i==1) {
      delta <- datas$delta
    }
    if (ind_tpts) {
      res$STIR <- STIR(X=datas$X, y=datas$y, t=t, p=p,basis_S = "Polynomial",H=H,d=d,k=k,time_pts = datas$time_pts,
                       calc_VC_vecB=FALSE, calc_covX=FALSE, comp_Xred = FALSE,cont_y = TRUE,comp_B_ML = TRUE)
    } else {
      res$STIR <- STIR(X=datas$X, y=datas$y, t=t, p=p,basis_S = "Polynomial",H=H,d=d,k=k,time_pts = datas$time_pts,
                       calc_VC_vecB=TRUE, calc_covX=FALSE, comp_Xred = FALSE,cont_y = TRUE,comp_B_ML = TRUE)
    }
    
    for (j in 1:1) {
      B_est[[j]][i,,] <- res[[j]]$B
    }
    B_ML_est[[1]][i,,] <- res[[1]]$B_ML
    delta_ls_est[[1]][i,,] <- res[[1]]$delta_ls
  }
  
  
  return(list(B_est=B_est,B_ML_est=B_ML_est,B=B,delta=delta,delta_ls_est=delta_ls_est,n=n,t=t,p=p,H=H,d=d)
  )
}

## dH >= T,p !!! 
myB12 <- matrix(c(rep(c(1,0.5),2),rep(c(0.7,0.3),2),rep(0.1,4)),4,3)
myB3 <- matrix(c(rep(c(1,0.5),5),rep(c(0.7,0.3),5),rep(0.1,10)),10,3)
myB4 <- matrix(c(rep(c(1,0.5,1),1),rep(c(0.7,0.3,0.7),1),rep(0.1,3)),3,3)
myB5 <- 2*matrix(c(1,  0,1/2,  0,1/2,  0,1,
                 0,1/4,1/2,1/4,1/2,1/4,0),2,7,byrow = TRUE)
#                       nrep, n,         B, t, p, H, d, rho, k,   list(r_lsir,m_lsir),  time, ind_tpts, delta_y
#              --------------------------------------------------------
params <- list( list( 2000, 500, myB12,10, 3, 2, 2,  0.8,   3,list(c(1,1),c(3,1),c(3,2),c(7,2)),"exp(1:t/6)/exp(t/6)", FALSE, FALSE),
                list( 2000, 500,    myB5,2, 7, 2, 1,  0.8,   2,list(c(1,1),c(1,3),c(2,4)),       "1:t/t", FALSE, FALSE)
)

# i <- 1
# what_time <- "exp(1:t/6)/exp(t/6)"
# npred <- 100
# nlev_ydis <- 10
# j <- 10
# traceback()
# warnings()

# k means now for all k=1,...,k
counter <- 1
nsettings <- length(params)
for (param in params) {
  cat(sprintf('Parameter setting %d / %d:\n',counter,nsettings))
  counter <- counter + 1
  c(nsim,n,B,t,p,H,d,rho,k,list_rm_lsir,time_sel,ind_tpts,delta_y) %<-% param
  sim <- simulateSDR(nsim, n, B, t, p, rho,H, d,k=k,list_rm_lsir=list_rm_lsir,npred=100,nlev_ydis = 10,what_time = time_sel,
                     ind_tpts = ind_tpts,delta_y = delta_y)

  saveRDS(sim, file = sprintf("../results/Delta_est_meas_contSTIR_%d_%d_%d_%d_%d_%d_%d_%.2f_%d_time_%.1s_ind_tpts_%s_delta_y_%s.rds",
                            nsim, n, round(sum(abs(B))), t, p, H, d, rho,k,time_sel,ind_tpts,delta_y))
}
