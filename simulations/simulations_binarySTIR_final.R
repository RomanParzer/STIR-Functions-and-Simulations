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
require(glmnet)
require(mgcv)

## returns X: T x p x n
generateDataSTIR <- function(n,B,t,p,rho,H=1L,d=1L,what_time = "1:t/t", basis_S = c("Fourier","Polynomial"),
                             center_S = TRUE, ind_tpts = FALSE, delta_y = FALSE) {
  y <- sample(0:1,n,replace=TRUE)
  
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
  
  if (!is.factor(y)) y <- factor(y)
  Gy <- model.matrix(~y)[,-1]
  Gy <- scale(Gy,scale=FALSE)
  
  if (delta_y) {
    delta1 <- 0.1*(rho^abs(outer(1:(t*p), 1:(t*p), FUN = `-`)))
    delta2 <- rho^abs(outer(1:(t*p), 1:(t*p), FUN = `-`))
  } else {
    delta1 <- delta2 <- rho^abs(outer(1:(t*p), 1:(t*p), FUN = `-`))
  }
  X <- sapply(1:n,function(i) {
    if (ind_tpts) {
      mu <- kronecker(S_t[((i-1)*t+1):(i*t),],matrix(Gy[i,],nrow = 1))%*%B
    } else {
      mu <- kronecker(S_t,matrix(Gy[i,],nrow = 1))%*%B
    }
    if (y[i]==0) {
      return(mvrnorm(1, c(mu), delta1))
    } else {
      return(mvrnorm(1, c(mu), delta2))
    }
  })
  dim(X) <- c(t,p,n)
  X <- sweep(X,1:2,apply(X,1:2,mean))
  return(list(X=X, y=y, B=B, S_t=S_t, Gy=Gy, delta1=delta1,delta2=delta2, time_pts = time_pts))
}

simulateSDR <- function(nsim, n, B, t, p, rho,H=1, d=1, k=1, npred = n, what_time = "1:t/t",
                        delta_y=FALSE, ind_tpts = FALSE) {
  res <- list(STIR=NULL,LSIR=NULL)
  
  X_red <- as.list(rep(0,5*k+2+k))
  names(X_red) <- c(paste(c("STIR_Del","STIR_Del_kron","STIR_Sig","STIR_Sig_kron","STIR_noCOV"),rep(1:k,each=5),sep="_k"),
                    "LSIR","orig",paste("STIR_ML",1:k,sep="_k"))
  new_X_red <- X_red 
  AUC_est <- array(c(0),dim=c(5*k + 2+k,3,nsim),dimnames=list(names(X_red),c("logit","logit w.2.ord.terms","gam")))
  
  
  for (i in 1:nsim) {
    if ((10*i)%%nsim == 0) {
      cat(sprintf('Starting rep %d/%d\n', i, nsim))
    }
    
    set.seed(i+1000)
    datas <- generateDataSTIR(n,B,t,p,rho,basis_S = "Polynomial",H=H,d=d,what_time=what_time,ind_tpts=ind_tpts, delta_y=delta_y)
    newdatas <- generateDataSTIR(npred,B,t,p,rho,basis_S = "Polynomial",d=d,H=H,what_time = what_time,ind_tpts=ind_tpts, delta_y=delta_y)
    res$STIR <- STIR(X=datas$X, y=datas$y, t=t, p=p,basis_S = "Polynomial",H=H,d=d,k=k,time_pts = datas$time_pts,
                     calc_VC_vecB=FALSE, calc_covX=TRUE, comp_Xred = FALSE,comp_B_ML = TRUE)
    res$LSIR <- LSIR(datas$X, y=datas$y, p=p, t=t,m=1,r=1)
    
    ### y estimation
    X <- datas$X
    dim(X) <- c(t*p,n)
    new_X <- newdatas$X
    dim(new_X) <- c(t*p, npred)
    svd_B_ML <- svd(t(res$STIR$B_ML),nu=k,nv=0)
    
    for (myk in 1:k) {
      X_red[[1+5*(myk-1)]] <- apply(X,2,function(xi) {
        matrix(solve(res$STIR$delta_ml,xi),t,p) %*% res$STIR$lsvs[,1:myk]
      })
      new_X_red[[1+5*(myk-1)]] <- apply(new_X,2,function(xi) {
        matrix(solve(res$STIR$delta_ml,xi),t,p) %*% res$STIR$lsvs[,1:myk]
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
    X_red[[5*k+2]] <- X
    new_X_red[[5*k+1]] <- crossprod(res$LSIR$B,new_X)
    new_X_red[[5*k+2]] <- new_X
    
    for (j in 1:(5*k + 2+k)) {
      dim <- nrow(X_red[[j]])
      ## logit
      if (dim==1) {
        fit_logit <- glm(y~.,family = binomial(link = "logit"),data = data.frame(X1=t(X_red[[j]]),y=datas$y))
        probs <- predict(fit_logit,data.frame(X1=t(new_X_red[[j]])))
      } else {
        fit_logit <- glm(y~.,family = binomial(link = "logit"),data = data.frame(t(X_red[[j]]),y=datas$y))
        probs <- predict(fit_logit,data.frame(t(new_X_red[[j]])))
      }
      # summary(fit_logit)
      AUC_est[j,1,i] <- performance(prediction(probs,newdatas$y),measure="auc")@y.values[[1]]
      
      myform <- formula(paste("y ~ (.)^2 + ",paste(paste("I(X",1:(dim),sep=""),collapse = "^2) + ",sep=""),"^2)"))
      if (dim==1) {
        fit_logit2 <- try(glm(myform,family = binomial(link = "logit"),data = data.frame(X1=t(X_red[[j]]),y=datas$y)),
                         silent = FALSE)
      } else {
        fit_logit2 <- try(glm(myform,family = binomial(link = "logit"),data = data.frame(t(X_red[[j]]),y=datas$y)),
                          silent = FALSE)
      }
      # summary(fit_logit2)
      if (!inherits(fit_logit2,"try-error")) {
        if (dim==1) {
          probs2 <- predict(fit_logit2,data.frame(X1=t(new_X_red[[j]])))
        } else {
          probs2 <- predict(fit_logit2,data.frame(t(new_X_red[[j]])))
        }
        AUC_est[j,2,i] <- performance(prediction(probs2,newdatas$y),measure="auc")@y.values[[1]]
      }
      
      
      # if (dim>1) {
      #   fit_glmnet <- cv.glmnet(x=t(X_red[[j]]),y=datas$y,family = "binomial",type.measure = "auc")
      #   preds <- predict(fit_glmnet,newx =t(new_X_red[[j]]),s="lambda.min")
      #   AUC_est[j,2,i] <- performance(prediction(c(preds),newdatas$y),measure="auc")@y.values[[1]]
      # }
      
      ### only for sc 5 6
      ## gam
      # myform <- formula(paste("y ~ ",paste(paste("s(X",1:(dim),sep=""),collapse = ") + "),")")) ## takes too much time
      # myform <- formula(paste("y ~ ",paste(paste("s(X",1:(dim),sep=""),collapse = ", bs = 'cr', k = 3) + "),")"))
      # if (dim==1) {
      #   gam_fit <- try(mgcv::gam(myform,data = data.frame(X1=t(X_red[[j]]),y=datas$y),family = binomial,method = "REML"),
      #                     silent = FALSE)
      # } else {
      #   gam_fit <- try(mgcv::gam(myform,data = data.frame(t(X_red[[j]]),y=datas$y),family = binomial,method = "REML"),
      #                     silent = FALSE)
      # }
      # # summary(gam_fit)
      # if (!inherits(gam_fit,"try-error")) {
      #   if (dim==1) {
      #     probs_gam <- predict(gam_fit,data.frame(X1=t(new_X_red[[j]])))
      #   } else {
      #     probs_gam <- predict(gam_fit,data.frame(t(new_X_red[[j]])))
      #   }
      #   AUC_est[j,3,i] <- performance(prediction(c(probs_gam),newdatas$y),measure="auc")@y.values[[1]]
      # }
    }
  }

  
  return(list(yAUC_mean = apply(AUC_est,1:2,mean), yAUC_sd=apply(AUC_est,1:2,sd),
              what_time=what_time,B=B
              ))
}

## dH >= rk(B) !!! 
myB12 <- matrix(c(rep(c(1,0.5),2),rep(c(0.7,0.3),2),rep(0.1,4)),4,3)
myB3 <- matrix(c(rep(c(1,0.5),5),rep(c(0.7,0.3),5),rep(0.1,10)),10,3)
myB4 <- matrix(c(rep(c(1,0.5,1),1),rep(c(0.7,0.3,0.7),1),rep(0.1,3)),3,3)
myB5 <- matrix(c(1,  0,1/2,  0,1/2,  0,1,
                 0,1/4,1/2,1/4,1/2,1/4,0),2,7,byrow = TRUE)
#                       n,         B, t, p, H, d, rho, k,                time, ind_tpts, delta_y
#              --------------------------------------------------------
params <- list( # list( 500,0.1*myB12,10, 3, 1, 4,  0.8, 3,"exp(1:t/6)/exp(t/6)", FALSE, TRUE),
                # list( 500, 0.1*myB12,10, 3, 1, 4, 0.8, 3,"exp(1:t/6)/exp(t/6)", FALSE, FALSE),
                # list( 500, 0.1*myB12,10, 3, 1, 4,  0.8, 3,           "runif(t)",  TRUE, TRUE),
                # list( 500, 0.1*myB12,10, 3, 1, 4,  0.8, 3,           "runif(t)",  TRUE, FALSE),
                # list( 500, 0.1*cbind(myB4,myB4,myB4,myB4,myB4),4, 15, 1, 3,  0.8, 3,"exp(1:t/6-t/6)", FALSE, TRUE),
  #               list( 500, 0.1*cbind(myB4,myB4,myB4,myB4,myB4),4, 15, 1, 3,   0.8, 3,"exp(1:t/6-t/6)", FALSE, FALSE),
  #               list( 500,       myB5,3, 7, 1, 2,   0.8, 2,             "1:t/t", FALSE, TRUE),
  #               list( 500,       myB5,3, 7, 1, 2,   0.8, 2,             "1:t/t", FALSE, FALSE),
  #               list( 500, cbind(myB5,myB5,myB5),3, 21, 1, 2, 0.8, 2,"1:t/t", FALSE, TRUE),
  #               list( 500, cbind(myB5,myB5,myB5),3, 21, 1, 2,    0.8, 2,"1:t/t", FALSE, FALSE),
  #               list( 500,       3*myB5,3, 7, 1, 2,   0.8, 2,             "1:t/t", FALSE, TRUE),
  #               list( 500,       3*myB5,3, 7, 1, 2,   0.8, 2,             "1:t/t", FALSE, FALSE),
  list( 2000,0.1*myB12,10, 3, 1, 4,  0.8, 1,"exp(1:t/6)/exp(t/6)", FALSE, TRUE),
  list( 2000, 0.1*myB12,10, 3, 1, 4, 0.8, 1,"exp(1:t/6)/exp(t/6)", FALSE, FALSE),
  list( 2000, 0.1*myB12,10, 3, 1, 4,  0.8, 1,           "runif(t)",  TRUE, TRUE),
  list( 2000, 0.1*myB12,10, 3, 1, 4,  0.8, 1,           "runif(t)",  TRUE, FALSE),
  list( 2000, 0.1*cbind(myB4,myB4,myB4,myB4,myB4),4, 15, 1, 3, 0.8, 1,"exp(1:t/6-t/6)", FALSE, TRUE),
  list( 2000, 0.1*cbind(myB4,myB4,myB4,myB4,myB4),4, 15, 1, 3,  0.8, 1,"exp(1:t/6-t/6)", FALSE, FALSE),
  list( 2000,       myB5,3, 7, 1, 2,0.8, 1,             "1:t/t", FALSE, TRUE),
  list( 2000,       myB5,3, 7, 1, 2 , 0.8, 1,             "1:t/t", FALSE, FALSE),
  list( 2000, cbind(myB5,myB5,myB5),3, 21, 1, 2, 0.8, 1,"1:t/t", FALSE, TRUE),
  list( 2000, cbind(myB5,myB5,myB5),3, 21, 1, 2,  0.8, 1,"1:t/t", FALSE, FALSE),
  list( 2000,       3*myB5,3, 7, 1, 2,0.8, 1,             "1:t/t", FALSE, TRUE),
  list( 2000,       3*myB5,3, 7, 1, 2 , 0.8, 1,             "1:t/t", FALSE, FALSE)
)

# i <- 1
# what_time <- "exp(1:t/6)/exp(t/6)"
# npred <- 100
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
  c(n,B,t,p,H,d,rho,k,time_sel,ind_tpts,delta_y) %<-% param
  sim <- simulateSDR(nsim, n, B, t, p, rho,H, d,k=k,npred=100,what_time = time_sel,
                     ind_tpts = ind_tpts,delta_y = delta_y)

  saveRDS(sim, file = sprintf("../results/simulation_binarySTIR_%d_%d_%d_%d_%d_%d_%d_%.2f_%d_time_%.1s_ind_tpts_%s_delta_y_%s.rds",
                            nsim, n, round(sum(abs(B))), t, p, H, d, rho,k,time_sel,ind_tpts,delta_y))
}

