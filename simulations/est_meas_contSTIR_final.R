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
  
  psi_est <- list(LSIR=array(c(0),dim=c(nsim,t,t)))
  phi_est <- list(LSIR=array(c(0),dim=c(nsim,p,p)))
  
  delta_p_est <- list(LSIR=array(c(0),dim=c(nsim,p,p)))
  delta_t_est <- list(LSIR=array(c(0),dim=c(nsim,t,t)))
  delta_ls_est <- list(STIR=array(c(0),dim=c(nsim,p*t,p*t)))
  
  rkB <- rankMatrix(B)[1]
  #scaling_est_acc <- sqrt(rkB) + sqrt(min(d*H,p))
  
  BminBhats <- array(c(0),dim=c(nsim))
  est_accs <- array(c(0),dim=c(nsim))
  princ_angs <- array(c(0),dim=c(nsim))
  
  BminBhats_ML <- array(c(0),dim=c(nsim))
  est_accs_ML <- array(c(0),dim=c(nsim))
  princ_angs_ML <- array(c(0),dim=c(nsim))
  
  chisqu_vals <- array(c(0),dim=c(nsim))
  rk_VC_Bhat <- array(c(0),dim=c(nsim))
  
  delta_Fs <- array(c(0),dim=c(nsim))
  delta_Fs_ML <- array(c(0),dim=c(nsim))
  
  for (i in 1:nsim) {
    if ((10*i)%%nsim == 0) {
      cat(sprintf('Starting rep %d/%d\n', i, nsim))
    }
    
    set.seed(i+1000)
    datas <- generateDataSTIR(n,B,t,p,rho,basis_S = "Polynomial",H=H,d=d,what_time=what_time,ind_tpts=ind_tpts, delta_y=delta_y)
    newdatas <- generateDataSTIR(npred,B,t,p,rho,basis_S = "Polynomial",d=d,H=H,what_time = what_time,ind_tpts=ind_tpts, delta_y=delta_y)
    if (i==1) {
      if (delta_y) {
        delta <- datas$delta*mean(exp(-pmin(sqrt(10)*abs(datas$y),2)))
      } else {
        delta <- datas$delta
      }
      
    }
    
    res$STIR <- STIR(X=datas$X, y=datas$y, t=t, p=p,basis_S = "Polynomial",H=H,d=d,k=k,time_pts = datas$time_pts,
                       calc_VC_vecB=TRUE, calc_covX=FALSE, comp_Xred = FALSE,cont_y = TRUE,comp_B_ML = TRUE)
    
    y_dis <- cut(pnorm(datas$y,0,0.1),breaks=(0:(nlev_ydis))/(nlev_ydis))
    new_y_dis <- cut(pnorm(newdatas$y,0,0.1),breaks=(0:(nlev_ydis))/(nlev_ydis))
    
    res$LSIR <- LSIR(datas$X, y=y_dis, p=p, t=t,m=p,r=t)
    
    for (j in 1:2) {
      B_est[[j]][i,,] <- res[[j]]$B
    }
    B_ML_est[[1]][i,,] <- res[[1]]$B_ML
    phi_est[[1]][i,,] <- res[[2]]$beta
    psi_est[[1]][i,,] <- res[[2]]$alpha
    delta_p_est[[1]][i,,] <- res[[2]]$delta_p
    delta_t_est[[1]][i,,] <- res[[2]]$delta_t
    delta_ls_est[[1]][i,,] <- res[[1]]$delta_ls
    
    
    BminBhats[i] <- norm(datas$B-res$STIR$B,'F')
    
    est_accs[i] <- norm(datas$B%*%tcrossprod(ginv(crossprod(datas$B)),datas$B) -
                          res[[1]]$B%*%tcrossprod(ginv(crossprod(res[[1]]$B)),res[[1]]$B), '2')
    
    princ_angs[i] <-  1-norm(crossprod(datas$B,ginv(tcrossprod(datas$B)))%*%datas$B%*%crossprod(res[[1]]$B,ginv(tcrossprod(res[[1]]$B)))%*%res[[1]]$B, '2')
    
    BminBhats_ML[i] <- norm(datas$B-res$STIR$B_ML,'F')
    
    est_accs_ML[i] <- norm(datas$B%*%tcrossprod(ginv(crossprod(datas$B)),datas$B) -
                          res[[1]]$B_ML%*%tcrossprod(ginv(crossprod(res[[1]]$B_ML)),res[[1]]$B_ML), '2')
    
    princ_angs_ML[i] <-  1-norm(crossprod(datas$B,ginv(tcrossprod(datas$B)))%*%datas$B%*%crossprod(res[[1]]$B_ML,ginv(tcrossprod(res[[1]]$B_ML)))%*%res[[1]]$B_ML, '2')
    
    delta_Fs[i] <- norm(delta-res[[1]]$delta_ls,'F')
    delta_Fs_ML[i] <- norm(delta-res[[1]]$delta_ls*(n-min(t,d)*d*H)/n,'F')
    
    
    if (TRUE) {
      rk_VC_Bhat[i] <- rankMatrix(res$STIR$VC_vecB)
      if (rk_VC_Bhat[i] < d*H*p) {
        chisqu_vals[i] <- crossprod(c(res$STIR$B-datas$B),ginv(res$STIR$VC_vecB)%*%c(res$STIR$B-datas$B))
      } else {
        chisqu_vals[i] <- crossprod(c(res$STIR$B-datas$B),solve(res$STIR$VC_vecB,c(res$STIR$B-datas$B)))
      }
    }
  }
  
  B_sd <- B_mean <- list(STIR=NULL,LSIR=NULL)
  phi_sd <- phi_mean <- psi_sd <- psi_mean <- delta_p_mean <- delta_t_mean <- list(LSIR=NULL) 
  delta_ls_mean <- B_ML_mean <- B_ML_sd <- list(STIR=NULL)
  
  if (dim(B_est[[1]])[2] == 1) {
    B_sd[[1]] <- matrix(apply(B_est[[1]],3,function(Bi) apply(Bi,2,sd)),nrow=1)
    B_mean[[1]] <- matrix(apply(B_est[[1]],3,function(Bi) apply(Bi,2,mean)),nrow=1)
    B_ML_sd[[1]] <- matrix(apply(B_ML_est[[1]],3,function(Bi) apply(Bi,2,sd)),nrow=1)
    B_ML_mean[[1]] <- matrix(apply(B_ML_est[[1]],3,function(Bi) apply(Bi,2,mean)),nrow=1)
  } else {
    B_sd[[1]] <- apply(B_est[[1]],3,function(Bi) apply(Bi,2,sd))
    B_mean[[1]] <- apply(B_est[[1]],3,function(Bi) apply(Bi,2,mean))
    B_ML_sd[[1]] <- apply(B_ML_est[[1]],3,function(Bi) apply(Bi,2,sd))
    B_ML_mean[[1]] <- apply(B_ML_est[[1]],3,function(Bi) apply(Bi,2,mean))
  }
  delta_p_mean[[1]] <- apply(delta_p_est[[1]],3,function(Si) apply(Si,2,mean))
  delta_t_mean[[1]] <- apply(delta_t_est[[1]],3,function(Si) apply(Si,2,mean))
  delta_ls_mean[[1]] <- apply(delta_ls_est[[1]],3,function(Si) apply(Si,2,mean))
  
  
  B_sd[[2]] <- apply(B_est[[2]],3,function(Bi) apply(Bi,2,sd))
  B_mean[[2]] <- apply(B_est[[2]],3,function(Bi) apply(Bi,2,mean))
  phi_sd[[1]] <- apply(phi_est[[1]],3,function(Bi) apply(Bi,2,sd))
  phi_mean[[1]] <- apply(phi_est[[1]],3,function(Bi) apply(Bi,2,mean))
  psi_sd[[1]] <- apply(psi_est[[1]],3,function(Bi) apply(Bi,2,sd))
  psi_mean[[1]] <- apply(psi_est[[1]],3,function(Bi) apply(Bi,2,mean))
  
  ks_test_res <- ks.test(chisqu_vals,pchisq,df=d*H*p)
  
  return(list(B_mean=B_mean, B_sd = B_sd,delta_p_mean=delta_p_mean, delta_t_mean=delta_t_mean, delta_ls_mean = delta_ls_mean,
              phi_mean=phi_mean, phi_sd=phi_sd, psi_mean=psi_mean, psi_sd=psi_sd, B_ML_mean=B_ML_mean,B_ML_sd=B_ML_sd,
              what_time=what_time,B=B,delta=delta,
              BminBhat = mean(BminBhats), BminBhat_sd = sd(BminBhats), 
              BminBhat_ML = mean(BminBhats_ML), BminBhat_ML_sd = sd(BminBhats_ML), 
              est_acc = mean(est_accs), est_acc_sd = sd(est_accs),
              princ_ang = mean(princ_angs), princ_ang_sd = sd(princ_angs),
              delta_F = mean(delta_Fs), delta_F_sd = sd(delta_Fs),
              delta_F_ML = mean(delta_Fs_ML), delta_F_ML_sd = sd(delta_Fs_ML),
              chisqu_vals=chisqu_vals, rk_VC_Bhat=rk_VC_Bhat, ks_test_res= ks_test_res,
              unb_delta = norm(delta_ls_mean[[1]]-delta,'F'),unb_delta_ML = norm(delta_ls_mean[[1]]*(n-min(d,t)*H*p)/n-delta,'F')
  ))
}

## dH >= T,p !!! 
myB12 <- matrix(c(rep(c(1,0.5),2),rep(c(0.7,0.3),2),rep(0.1,4)),4,3)
myB3 <- matrix(c(rep(c(1,0.5),5),rep(c(0.7,0.3),5),rep(0.1,10)),10,3)
myB4 <- matrix(c(rep(c(1,0.5,1),1),rep(c(0.7,0.3,0.7),1),rep(0.1,3)),3,3)
myB5 <- 2*matrix(c(1,  0,1/2,  0,1/2,  0,1,
                 0,1/4,1/2,1/4,1/2,1/4,0),2,7,byrow = TRUE)
#                       n,         B, t, p, H, d, rho, k,   list(r_lsir,m_lsir),  time, ind_tpts, delta_y
#              --------------------------------------------------------
params <- list( #list( 200,myB12,10, 3, 2, 2,  0.8, 3,list(c(1,1),c(3,1),c(3,2),c(7,2)),"exp(1:t/6)/exp(t/6)", FALSE, TRUE),
                # list( 200, myB12,10, 3, 2, 2,  0.8,   3,list(c(1,1),c(3,1),c(3,2),c(7,2)),"exp(1:t/6)/exp(t/6)", FALSE, FALSE),
                list( 500,myB12,10, 3, 2, 2,  0.8, 3,list(c(1,1),c(3,1),c(3,2),c(7,2)),"exp(1:t/6)/exp(t/6)", FALSE, TRUE),
                list( 500, myB12,10, 3, 2, 2,  0.8,   3,list(c(1,1),c(3,1),c(3,2),c(7,2)),"exp(1:t/6)/exp(t/6)", FALSE, FALSE),
                list( 2000,myB12,10, 3, 2, 2,  0.8, 3,list(c(1,1),c(3,1),c(3,2),c(7,2)),"exp(1:t/6)/exp(t/6)", FALSE, TRUE),
                list( 2000, myB12,10, 3, 2, 2,  0.8,   3,list(c(1,1),c(3,1),c(3,2),c(7,2)),"exp(1:t/6)/exp(t/6)", FALSE, FALSE),
                list( 10000,myB12,10, 3, 2, 2,  0.8, 3,list(c(1,1),c(3,1),c(3,2),c(7,2)),"exp(1:t/6)/exp(t/6)", FALSE, TRUE),
                list( 10000, myB12,10, 3, 2, 2,  0.8,   3,list(c(1,1),c(3,1),c(3,2),c(7,2)),"exp(1:t/6)/exp(t/6)", FALSE, FALSE),
                # list( 200, myB12,10, 3, 2, 2,  0.8,  3,list(c(1,1),c(3,1),c(3,2),c(7,2)),       "runif(t)",  TRUE, TRUE),
                # list( 200, myB12,10, 3, 2, 2,  0.8,   3,list(c(1,1),c(3,1),c(3,2),c(7,2)),       "runif(t)",  TRUE, FALSE),
                list( 500, myB12,10, 3, 2, 2,  0.8,  3,list(c(1,1),c(3,1),c(3,2),c(7,2)),       "runif(t)",  TRUE, TRUE),
                list( 500, myB12,10, 3, 2, 2,  0.8,   3,list(c(1,1),c(3,1),c(3,2),c(7,2)),       "runif(t)",  TRUE, FALSE),
                list( 2000, myB12,10, 3, 2, 2,  0.8,  3,list(c(1,1),c(3,1),c(3,2),c(7,2)),       "runif(t)",  TRUE, TRUE),
                list( 2000, myB12,10, 3, 2, 2,  0.8,   3,list(c(1,1),c(3,1),c(3,2),c(7,2)),       "runif(t)",  TRUE, FALSE),
                list( 10000, myB12,10, 3, 2, 2,  0.8,  3,list(c(1,1),c(3,1),c(3,2),c(7,2)),       "runif(t)",  TRUE, TRUE),
                list( 10000, myB12,10, 3, 2, 2,  0.8,   3,list(c(1,1),c(3,1),c(3,2),c(7,2)),       "runif(t)",  TRUE, FALSE),
                # list( 200, cbind(myB4,myB4,myB4,myB4,myB4),3, 15, 3, 1,  0.8,   3,list(c(1,1),c(1,3),c(2,3),c(2,10)),"exp(1:t/6-t/6)", FALSE, TRUE),
                # list( 200, cbind(myB4,myB4,myB4,myB4,myB4),3, 15, 3, 1,  0.8,   3,list(c(1,1),c(1,3),c(2,3),c(2,10)),"exp(1:t/6-t/6)", FALSE, FALSE),
                list( 500, cbind(myB4,myB4,myB4,myB4,myB4),3, 15, 3, 1,  0.8,   3,list(c(1,1),c(1,3),c(2,3),c(2,10)),"exp(1:t/6-t/6)", FALSE, TRUE),
                list( 500, cbind(myB4,myB4,myB4,myB4,myB4),3, 15, 3, 1,  0.8,   3,list(c(1,1),c(1,3),c(2,3),c(2,10)),"exp(1:t/6-t/6)", FALSE, FALSE),
                list( 2000, cbind(myB4,myB4,myB4,myB4,myB4),3, 15, 3, 1,  0.8,   3,list(c(1,1),c(1,3),c(2,3),c(2,10)),"exp(1:t/6-t/6)", FALSE, TRUE),
                list( 2000, cbind(myB4,myB4,myB4,myB4,myB4),3, 15, 3, 1,  0.8,   3,list(c(1,1),c(1,3),c(2,3),c(2,10)),"exp(1:t/6-t/6)", FALSE, FALSE),
                list( 10000, cbind(myB4,myB4,myB4,myB4,myB4),3, 15, 3, 1,  0.8,   3,list(c(1,1),c(1,3),c(2,3),c(2,10)),"exp(1:t/6-t/6)", FALSE, TRUE),
                list( 10000, cbind(myB4,myB4,myB4,myB4,myB4),3, 15, 3, 1,  0.8,   3,list(c(1,1),c(1,3),c(2,3),c(2,10)),"exp(1:t/6-t/6)", FALSE, FALSE),
                # list( 200,       myB5,2, 7, 2, 1,  0.8,   2,list(c(1,1),c(1,3),c(2,4)),       "1:t/t", FALSE, TRUE),
                # list( 200,       myB5,2, 7, 2, 1,  0.8,   2,list(c(1,1),c(1,3),c(2,4)),       "1:t/t", FALSE, FALSE),
                list( 500,       myB5,2, 7, 2, 1,  0.8,   2,list(c(1,1),c(1,3),c(2,4)),       "1:t/t", FALSE, TRUE),
                list( 500,       myB5,2, 7, 2, 1,  0.8,   2,list(c(1,1),c(1,3),c(2,4)),       "1:t/t", FALSE, FALSE),
                list( 2000,       myB5,2, 7, 2, 1,  0.8,   2,list(c(1,1),c(1,3),c(2,4)),       "1:t/t", FALSE, TRUE),
                list( 2000,       myB5,2, 7, 2, 1,  0.8,   2,list(c(1,1),c(1,3),c(2,4)),       "1:t/t", FALSE, FALSE),
                list( 10000,       myB5,2, 7, 2, 1,  0.8,   2,list(c(1,1),c(1,3),c(2,4)),       "1:t/t", FALSE, TRUE),
                list( 10000,       myB5,2, 7, 2, 1,  0.8,   2,list(c(1,1),c(1,3),c(2,4)),       "1:t/t", FALSE, FALSE),
                # list( 200, cbind(myB5,myB5,myB5),2, 21, 2, 1, 0.8,  2,list(c(1,1),c(1,4),c(2,10)), "1:t/t", FALSE, TRUE),
                # list( 200, cbind(myB5,myB5,myB5),2, 21, 2, 1,  0.8,  2,list(c(1,1),c(1,4),c(2,10)),"1:t/t", FALSE, FALSE),
                list( 500, cbind(myB5,myB5,myB5),2, 21, 2, 1, 0.8,  2,list(c(1,1),c(1,4),c(2,10)), "1:t/t", FALSE, TRUE),
                list( 500, cbind(myB5,myB5,myB5),2, 21, 2, 1,  0.8,  2,list(c(1,1),c(1,4),c(2,10)),"1:t/t", FALSE, FALSE),
                list( 2000, cbind(myB5,myB5,myB5),2, 21, 2, 1, 0.8,  2,list(c(1,1),c(1,4),c(2,10)), "1:t/t", FALSE, TRUE),
                list( 2000, cbind(myB5,myB5,myB5),2, 21, 2, 1,  0.8,  2,list(c(1,1),c(1,4),c(2,10)),"1:t/t", FALSE, FALSE),
                list( 10000, cbind(myB5,myB5,myB5),2, 21, 2, 1, 0.8,  2,list(c(1,1),c(1,4),c(2,10)), "1:t/t", FALSE, TRUE),
                list( 10000, cbind(myB5,myB5,myB5),2, 21, 2, 1,  0.8,  2,list(c(1,1),c(1,4),c(2,10)),"1:t/t", FALSE, FALSE),
                # list( 200,       3*myB5,2, 7, 2, 1,  0.8,   2,list(c(1,1),c(1,3),c(2,4)),       "1:t/t", FALSE, TRUE),
                # list( 200,       3*myB5,2, 7, 2, 1,  0.8,   2,list(c(1,1),c(1,3),c(2,4)),       "1:t/t", FALSE, FALSE)
                list( 500,       3*myB5,2, 7, 2, 1,  0.8,   2,list(c(1,1),c(1,3),c(2,4)),       "1:t/t", FALSE, TRUE),
                list( 500,       3*myB5,2, 7, 2, 1,  0.8,   2,list(c(1,1),c(1,3),c(2,4)),       "1:t/t", FALSE, FALSE),
                list( 2000,       3*myB5,2, 7, 2, 1,  0.8,   2,list(c(1,1),c(1,3),c(2,4)),       "1:t/t", FALSE, TRUE),
                list( 2000,       3*myB5,2, 7, 2, 1,  0.8,   2,list(c(1,1),c(1,3),c(2,4)),       "1:t/t", FALSE, FALSE),
                list( 10000,       3*myB5,2, 7, 2, 1,  0.8,   2,list(c(1,1),c(1,3),c(2,4)),       "1:t/t", FALSE, TRUE),
                list( 10000,       3*myB5,2, 7, 2, 1,  0.8,   2,list(c(1,1),c(1,3),c(2,4)),       "1:t/t", FALSE, FALSE)
)

# i <- 1
# what_time <- "exp(1:t/6)/exp(t/6)"
# npred <- 100
# nlev_ydis <- 10
# j <- 10
# traceback()
# warnings()

# k means now for all k=1,...,k
nsim <- 500
counter <- 1
nsettings <- length(params)
for (param in params) {
  cat(sprintf('Parameter setting %d / %d:\n',counter,nsettings))
  counter <- counter + 1
  c(n,B,t,p,H,d,rho,k,list_rm_lsir,time_sel,ind_tpts,delta_y) %<-% param
  sim <- simulateSDR(nsim, n, B, t, p, rho,H, d,k=k,list_rm_lsir=list_rm_lsir,npred=100,nlev_ydis = 10,what_time = time_sel,
                     ind_tpts = ind_tpts,delta_y = delta_y)

  saveRDS(sim, file = sprintf("../results/Final_est_meas_contSTIR_%d_%d_%d_%d_%d_%d_%d_%.2f_%d_time_%.1s_ind_tpts_%s_delta_y_%s.rds",
                            nsim, n, round(sum(abs(B))), t, p, H, d, rho,k,time_sel,ind_tpts,delta_y))
}
