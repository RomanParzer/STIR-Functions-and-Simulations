#source('./random.R')
source('../functions/STIR.R')
source('../functions/multi_assign.R')
source('../functions/lsir_adapted.R')
source('../functions/RrankTest.R')

#require(dplyr)
#require(stats)
require(MASS) # for ginv and mvrnorm
require(Matrix)
require(mvtnorm)
require(ROCR)
require(mgcv)

## returns X: T x p x n
generateDataLSIR <- function(n,phi,psi,t,p,rho,what_time = "1:t/t",ind_tpts=FALSE, delta_y=FALSE) {
  y <- rnorm(n,0,0.1)
  
  if (ind_tpts) {
    time_pts <- sapply(1:n,function(i) {
      sort(eval(str2expression(what_time)))
    })
  } else {
    time_pts <- sort(eval(str2expression(what_time)))
  }
  
  H <- ncol(phi)*ncol(psi)
  Gy <- do.call(cbind, lapply(1:(H%/%2), function(s, z) {
    cbind(cos(s * z), sin(s * z))
  }, z = 2 * pi * y))
  if (H%%2 == 1) {
    Gy <- cbind(Gy,cos((H%/%2+1)*2*pi*y))
  }
  Gy <- scale(Gy,scale=FALSE)
  
  B <- as.matrix(kronecker(phi,psi))
  
  delta <- rho^abs(outer(1:(t*p), 1:(t*p), FUN = `-`))
  
  X <- sapply(1:n,function(i) {
    mu <- tcrossprod(B,matrix(Gy[i,],nrow = 1)) 
    if (delta_y) {
      myscale <- exp(-min(sqrt(10)*abs(y[i]),2))
    } else {
      myscale <- 1
    }
    mvrnorm(1, c(mu), myscale*delta)
  })
  dim(X) <- c(t,p,n)
  X <- sweep(X,1:2,apply(X,1:2,mean))
  return(list(X=X, y=y, B=B, delta=delta, time_pts = time_pts))
}

simulateSDR <- function(nsim, n, phi, psi, t, p, rho,H=1, d=1, k=1, npred = 100, what_time = "1:t/t",nlev_ydis = 10,
                        delta_y=FALSE, ind_tpts = FALSE) {
  res <- list(STIR=NULL,LSIR=NULL)
  
  X_red <- as.list(rep(0,5*k+2+k))
  names(X_red) <- c(paste(c("STIR_Del","STIR_Del_kron","STIR_Sig","STIR_Sig_kron","STIR_noCOV"),rep(1:k,each=5),sep="_k"),
                    "LSIR","orig",paste("STIR_ML",1:k,sep="_k"))
  new_X_red <- X_red 
  
  predmodBIC <- y_mse <- y_cor <- array(c(0),dim=c(5*k + 2+k,3,nsim),dimnames=list(names(X_red),c("lm","lm w.2.ord.terms","gam")))
  rkTest <- matrix(c(0),2,nsim)
  rownames(rkTest) <- c("Est. Rank","p-value")
  
  for (i in 1:nsim) {
    if ((10*i)%%nsim == 0) {
      cat(sprintf('Starting rep %d/%d\n', i, nsim))
    }
    
    set.seed(i+1000)
    datas <- generateDataLSIR(n,phi,psi,t,p,rho,what_time = what_time,ind_tpts=ind_tpts, delta_y=delta_y)
    newdatas <- generateDataLSIR(npred,phi,psi,t,p,rho,what_time = what_time,ind_tpts=ind_tpts, delta_y=delta_y)
    res$STIR <- STIR(X=datas$X, y=datas$y, t=t, p=p,basis_S = "Polynomial",H=H,d=d,k=k,time_pts = datas$time_pts,
                     calc_VC_vecB=TRUE, calc_covX=TRUE, comp_Xred = FALSE,cont_y = TRUE,comp_B_ML = TRUE)
    y_dis <- cut(pnorm(datas$y,0,0.1),breaks=(0:(nlev_ydis))/(nlev_ydis))
    new_y_dis <- cut(pnorm(newdatas$y,0,0.1),breaks=(0:(nlev_ydis))/(nlev_ydis))
    
    res$LSIR <- LSIR(datas$X, y=y_dis, p=p, t=t,m=p,r=t)
    
    Rktestres <- TestRMRank(res$STIR$B,res$STIR$VC_vecB*n, n, 0.05)
    rkTest[1,i] <- Rktestres$k
    rkTest[2,i] <- tail(Rktestres$pvals,1)
    
    ### y estimation
    X <- datas$X
    dim(X) <- c(t*p,n)
    new_X <- newdatas$X
    dim(new_X) <- c(t*p, npred)
    svd_B_ML <- svd(t(res$STIR$B_ML),nu=k,nv=0)
    
    for (myk in 1:k) {
      X_red[[1+5*(myk-1)]] <- apply(X,2,function(xi) {
        matrix(solve(res$STIR$delta_ls,xi),t,p) %*% res$STIR$lsvs[,1:myk]
      })
      new_X_red[[1+5*(myk-1)]] <- apply(new_X,2,function(xi) {
        matrix(solve(res$STIR$delta_ls,xi),t,p) %*% res$STIR$lsvs[,1:myk]
      })
      X_red[[2+5*(myk-1)]] <- apply(X,2,function(xi) {
        matrix(kronecker(solve(res$STIR$delta_p),solve(res$STIR$delta_t))%*%xi,t,p) %*% res$STIR$lsvs[,1:myk]
      })
      new_X_red[[2+5*(myk-1)]] <- apply(new_X,2,function(xi) {
        matrix(kronecker(solve(res$STIR$delta_p),solve(res$STIR$delta_t))%*%xi,t,p) %*% res$STIR$lsvs[,1:myk]
      })
      X_red[[3+5*(myk-1)]] <- apply(X,2,function(xi) {
        matrix(solve(res$STIR$covX,xi),t,p) %*% res$STIR$lsvs[,1:myk]
      })
      new_X_red[[3+5*(myk-1)]] <- apply(new_X,2,function(xi) {
        matrix(solve(res$STIR$covX,xi),t,p) %*% res$STIR$lsvs[,1:myk]
      })
      X_red[[4+5*(myk-1)]] <- apply(X,2,function(xi) {
        matrix(kronecker(solve(res$STIR$Sig_p),solve(res$STIR$Sig_t))%*%xi,t,p) %*% res$STIR$lsvs[,1:myk]
      })
      new_X_red[[4+5*(myk-1)]] <- apply(new_X,2,function(xi) {
        matrix(kronecker(solve(res$STIR$Sig_p),solve(res$STIR$Sig_t))%*%xi,t,p) %*% res$STIR$lsvs[,1:myk]
      })
      X_red[[5+5*(myk-1)]] <- apply(X,2,function(xi) {
        matrix(xi,t,p) %*% res$STIR$lsvs[,1:myk]
      })
      new_X_red[[5+5*(myk-1)]] <- apply(new_X,2,function(xi) {
        matrix(xi,t,p) %*% res$STIR$lsvs[,1:myk]
      })
      X_red[[2+5*k + myk]] <- apply(X,2,function(xi) {
        matrix(solve(res$STIR$delta_ml,xi),t,p) %*% svd_B_ML$u[,1:myk]
      })
      new_X_red[[2+5*k +myk]] <- apply(new_X,2,function(xi) {
        matrix(solve(res$STIR$delta_ml,xi),t,p) %*% svd_B_ML$u[,1:myk]
      })
    }
    X_red[[5*k+1]] <- crossprod(res$LSIR$B,X)
    new_X_red[[5*k+1]] <- crossprod(res$LSIR$B,new_X)
    X_red[[5*k+2]] <- X
    new_X_red[[5*k+2]] <- new_X
    
    for (j in 1:(5*k + 2 + k)) {
      dim <- nrow(X_red[[j]])
      ### lm 
      if (dim==1) {
        lm_fit <- lm(y~.,data = data.frame(X1=t(X_red[[j]]),y=datas$y))
        preds <- predict.lm(lm_fit,newdata = data.frame(X1=t(new_X_red[[j]])))
      } else {
        lm_fit <- lm(y~.,data = data.frame(t(X_red[[j]]),y=datas$y))
        preds <- predict.lm(lm_fit,newdata = data.frame(t(new_X_red[[j]])))
      }
      #summary(lm_fit)
      
      y_mse[j,1,i] <- sum((preds-newdatas$y)^2) / sum((newdatas$y-mean(newdatas$y))^2)
      y_cor[j,1,i] <- cor(preds,newdatas$y)
      predmodBIC[j,1,i] <- BIC(lm_fit)
      
      ### lm + 2nd order
      myform2 <- formula(paste("y ~ (.)^2 + ",paste(paste("I(X",1:(dim),sep=""),collapse = "^2) + ",sep=""),"^2)"))
      if (dim==1) {
        fit_lm2 <- try(lm(myform2,data = data.frame(X1=t(X_red[[j]]),y=datas$y)),
                       silent = FALSE)
      } else {
        fit_lm2 <- try(lm(myform2,data = data.frame(t(X_red[[j]]),y=datas$y)),
                       silent=FALSE)
      }
      
      if (!inherits(fit_lm2,"try-error")) {
        if (dim==1) {
          preds2 <- predict.lm(fit_lm2,data.frame(X1=t(new_X_red[[j]])))
        } else {
          preds2 <- predict.lm(fit_lm2,data.frame(t(new_X_red[[j]])))
        }
        y_mse[j,2,i] <- sum((preds2-newdatas$y)^2) / sum((newdatas$y-mean(newdatas$y))^2)
        y_cor[j,2,i] <- cor(preds2,newdatas$y)
        predmodBIC[j,2,i] <- BIC(fit_lm2)
      }
      
      ## gam
      myform <- formula(paste("y ~ ",paste(paste("s(X",1:(dim),sep=""),collapse = ", bs = 'cr', k = 3) + "),")"))
      #myform <- formula(paste("y ~ ",paste(paste("s(X",1:(dim),sep=""),collapse = ") + "),")"))
      
      if (dim==1) {
        gam_fit <- try(mgcv::gam(myform,data = data.frame(X1=t(X_red[[j]]),y=datas$y,method="REML")),silent = FALSE)
      } else {
        gam_fit <- try(mgcv::gam(myform,data = data.frame(t(X_red[[j]]),y=datas$y,method="REML")),silent = FALSE)
      }
      # if (dim==1) {
      #   gam_fit <- try(mgcv::gam(myform,data = data.frame(X1=t(X_red[[j]]),y=datas$y)),silent = FALSE)
      # } else {
      #   gam_fit <- try(mgcv::gam(myform,data = data.frame(t(X_red[[j]]),y=datas$y)),silent = FALSE)
      # }
      
      #summary(gam_fit)
      if (!inherits(gam_fit, "try-error")) {
        if (dim==1) {
          preds <- predict(gam_fit,data.frame(X1=t(new_X_red[[j]])))
        } else {
          preds <- predict(gam_fit,data.frame(t(new_X_red[[j]])))
        }
        
        y_mse[j,3,i] <- sum((preds-newdatas$y)^2) / sum((newdatas$y-mean(newdatas$y))^2)
        y_cor[j,3,i] <- cor(preds,newdatas$y)
      }
    }
  }
  
  return(list(what_time=what_time, phi=phi,psi=psi,
              y_mse_mean = apply(y_mse,1:2,mean), y_mse_sd = apply(y_mse,1:2,sd),y_mse=y_mse,
              y_cor_mean = apply(y_cor,1:2,mean), y_cor_sd = apply(y_cor,1:2,sd),y_cor=y_cor,
              BICmean = apply(predmodBIC,1:2,mean), BICsd=apply(predmodBIC,1:2,sd),predmodBIC=predmodBIC,rkTest=rkTest
              ))
}

## dH >= rk(B) !!! 
myB12 <- matrix(c(rep(c(1,0.5),2),rep(c(0.7,0.3),2),rep(0.1,4)),4,3)
myB3 <- matrix(c(rep(c(1,0.5),5),rep(c(0.7,0.3),5),rep(0.1,10)),10,3)
myB4 <- matrix(c(rep(c(1,0.5,1),1),rep(c(0.7,0.3,0.7),1),rep(0.1,3)),3,3)
myB5 <- matrix(c(1,  0,1/2,  0,1/2,  0,1,
                 0,1/4,1/2,1/4,1/2,1/4,0),2,7,byrow = TRUE)
#                       n,          phi p x m,                            psi txr, t, p, H, d, rho, k,  time, ind_tpts, delta_y
#              --------------------------------------------------------
params <- list( list( 500, matrix(c(rep(c(0.1,-0.05),4),rep(c(1,-1),4)),nrow=8),matrix(c(1,0,1),nrow=3), 3, 8, 3, 2,   0.8, 3,"exp(1:t/6)/exp(t/6)", FALSE, TRUE),
                list( 500, matrix(c(rep(c(0.1,-0.05),4),rep(c(1,-1),4)),nrow=8),matrix(c(1,0,1),nrow=3), 3, 8, 3, 2,   0.8, 3,"exp(1:t/6)/exp(t/6)", FALSE, FALSE),
                list( 500, matrix(c(rep(c(0.1,-0.05),4),rep(c(1,-1),4)),nrow=8),matrix(c(1,0,1),nrow=3), 3, 8, 3, 2,   0.8, 3,"runif(t)", TRUE, TRUE),
                list( 500, matrix(c(rep(c(0.1,-0.05),4),rep(c(1,-1),4)),nrow=8),matrix(c(1,0,1),nrow=3), 3, 8, 3, 2,   0.8, 3,"runif(t)", TRUE, FALSE)
  # list( 2000, matrix(c(rep(c(0.1,-0.05),4),rep(c(1,-1),4)),nrow=8),matrix(c(1,0,1),nrow=3), 3, 8, 3, 2,   0.8, 1,"exp(1:t/6)/exp(t/6)", FALSE, TRUE),
  # list( 2000, matrix(c(rep(c(0.1,-0.05),4),rep(c(1,-1),4)),nrow=8),matrix(c(1,0,1),nrow=3), 3, 8, 3, 2,   0.8, 1,"exp(1:t/6)/exp(t/6)", FALSE, FALSE),
  # list( 2000, matrix(c(rep(c(0.1,-0.05),4),rep(c(1,-1),4)),nrow=8),matrix(c(1,0,1),nrow=3), 3, 8, 3, 2,   0.8, 1,"runif(t)", TRUE, TRUE),
  # list( 2000, matrix(c(rep(c(0.1,-0.05),4),rep(c(1,-1),4)),nrow=8),matrix(c(1,0,1),nrow=3), 3, 8, 3, 2,   0.8, 1,"runif(t)", TRUE, FALSE)
)

# i <- 1
# what_time <- "exp(1:t/6)/exp(t/6)"
# npred <- 100
# nlev_ydis <- 10
# j <- 10
# traceback()
# warnings()

# k means now for all k=1,...,k
nsim <- 100
counter <- 1
nsettings <- length(params)
for (param in params) {
  cat(sprintf('Parameter setting %d / %d:\n',counter,nsettings))
  counter <- counter + 1
  c(n,phi,psi,t,p,H,d,rho,k,time_sel,ind_tpts,delta_y) %<-% param
  sim <- simulateSDR(nsim, n, phi, psi, t, p, rho,H, d,k=k,npred=100,nlev_ydis = 10,what_time = time_sel,
                     ind_tpts = ind_tpts,delta_y = delta_y)

  saveRDS(sim, file = sprintf("../results/rank_contLSIR_%d_%d_%d_%d_%d_%d_%d_%.2f_%d_time_%.1s_ind_tpts_%s_delta_y_%s.rds",
                            nsim, n, round(sum(abs(phi))+sum(abs(psi))), t, p, H, d, rho,k,time_sel,ind_tpts,delta_y))
}

