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
require(glmnet)

## returns X: T x p x n

generateDataLSIR <- function(n,phi,psi,t,p,rho=0,what_time="1:t/t",ind_tpts = FALSE,delta_y = FALSE) {
  y <- sample(0:1,n,replace=TRUE)
  
  y <- as.factor(y)
  Gy <- model.matrix(~y)[,-1]
  Gy <- scale(Gy,scale=FALSE)
  
  B <- as.matrix(kronecker(phi,psi))
  
  if (delta_y) {
    delta1 <- 0.1*(rho^abs(outer(1:(t*p), 1:(t*p), FUN = `-`)))
    delta2 <- rho^abs(outer(1:(t*p), 1:(t*p), FUN = `-`))
  } else {
    delta1 <- delta2 <- rho^abs(outer(1:(t*p), 1:(t*p), FUN = `-`))
  }
  
  X <- sapply(1:n,function(i) {
    mu <- tcrossprod(B,Gy[i,]) 
    if (y[i]==0) {
      return(mvrnorm(1, c(mu), delta1))
    } else {
      return(mvrnorm(1, c(mu), delta2))
    }
  })
  dim(X) <- c(t,p,n)
  X <- sweep(X,1:2,apply(X,1:2,mean))
  if (ind_tpts) {
    time_pts <- sapply(1:n,function(i) {
      sort(eval(str2expression(what_time)))
    })
  } else {
    time_pts <- sort(eval(str2expression(what_time)))
  }
  return(list(X=X, y=y, B=B, delta1=delta1,delta2=delta2, time_pts=time_pts))
}


simulateSDR <- function(nsim, n, phi,psi, t, p, rho,H=1, d=1, k=1, npred = n, what_time = "1:t/t",
                        delta_y=FALSE, ind_tpts = FALSE) {
  res <- list(STIR=NULL,LSIR=NULL)
  
  X_red <- as.list(rep(0,5*k+2+k))
  names(X_red) <- c(paste(c("STIR_Del","STIR_Del_kron","STIR_Sig","STIR_Sig_kron","STIR_noCOV"),rep(1:k,each=5),sep="_k"),
                    "LSIR","orig",paste("STIR_ML",1:k,sep="_k"))
  new_X_red <- X_red 
  predmodBIC <- AUC_est <- array(c(0),dim=c(5*k + 2+k,2,nsim),dimnames=list(names(X_red),c("logit","logit w.2.ord.terms")))
  rkTest <- matrix(c(0),2,nsim)
  rownames(rkTest) <- c("Est. Rank","p-value")
  
  for (i in 1:nsim) {
    if ((10*i)%%nsim == 0) {
      cat(sprintf('Starting rep %d/%d\n', i, nsim))
    }
    
    set.seed(i+1000)
    datas <- generateDataLSIR(n,phi,psi,t,p,rho,what_time=what_time,ind_tpts=ind_tpts, delta_y=delta_y)
    newdatas <- generateDataLSIR(npred,phi,psi,t,p,rho,what_time=what_time,ind_tpts=ind_tpts, delta_y=delta_y)
    res$STIR <- STIR(X=datas$X, y=datas$y, t=t, p=p,basis_S = "Polynomial",H=H,d=d,k=k,time_pts = datas$time_pts,
                     calc_VC_vecB=TRUE, calc_covX=TRUE, comp_Xred = FALSE,comp_B_ML = TRUE)
    res$LSIR <- LSIR(datas$X, y=datas$y, p=p, t=t,m=1,r=1)
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
    X_red[[5*k+2]] <- X
    new_X_red[[5*k+1]] <- crossprod(res$LSIR$B,new_X)
    new_X_red[[5*k+2]] <- new_X
    
    for (j in 1:(5*k + 2 + k)) {
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
      predmodBIC[j,1,i] <- BIC(fit_logit)
      
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
        predmodBIC[j,2,i] <- BIC(fit_logit2)
      }
      
    }
  }
  
  return(list(yAUC_mean = apply(AUC_est,1:2,mean), yAUC_sd=apply(AUC_est,1:2,sd),AUC_est=AUC_est,
              BICmean = apply(predmodBIC,1:2,mean), BICsd=apply(predmodBIC,1:2,sd),predmodBIC=predmodBIC,
              what_time=what_time,phi=phi,psi=psi,rkTest=rkTest
  ))
}

## dH >= T,p !!! 
myB12 <- matrix(c(rep(c(1,0.5),2),rep(c(0.7,0.3),2),rep(0.1,4)),4,3)
myB3 <- matrix(c(rep(c(1,0.5),5),rep(c(0.7,0.3),5),rep(0.1,10)),10,3)
myB4 <- matrix(c(rep(c(1,0.5,1),1),rep(c(0.7,0.3,0.7),1),rep(0.1,3)),3,3)
myB5 <- matrix(c(1,  0,1/2,  0,1/2,  0,1,
                 0,1/4,1/2,1/4,1/2,1/4,0),2,7,byrow = TRUE)
#                       n,          phi,     psi,   t, p, H, d, rho, k,                time, ind_tpts, delta_y
#              --------------------------------------------------------
params <- list( list( 500, c(1,0.5,0.1),(1:8)/8,      8, 3, 1, 4,   0.8, 3, "exp(1:t/6-t/6)", FALSE, TRUE),
                list( 500, c(1,0.5,0.1),(1:8)/8,      8, 3, 1, 4,   0.8, 3, "exp(1:t/6-t/6)", FALSE, FALSE),
  list( 500, c(1,0.5,0.1),(1:8)/8,      8, 3, 1, 4,   0.8, 3, "runif(t)", TRUE, TRUE),
  list( 500, c(1,0.5,0.1),(1:8)/8,      8, 3, 1, 4,   0.8, 3, "runif(t)", TRUE, FALSE)
                
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
  c(n,phi,psi,t,p,H,d,rho,k,time_sel,ind_tpts,delta_y) %<-% param
  sim <- simulateSDR(nsim, n, phi,psi, t, p, rho,H, d,k=k,npred=100,what_time = time_sel,
                     ind_tpts = ind_tpts,delta_y = delta_y)

  saveRDS(sim, file = sprintf("../results/rank_binaryLSIR_%d_%d_%d_%d_%d_%d_%d_%.2f_%d_time_%.1s_ind_tpts_%s_delta_y_%s.rds",
                            nsim, n, round(sum(abs(phi))+sum(abs(psi))), t, p, H, d, rho,k,time_sel,ind_tpts,delta_y))
}

# param <- params[[1]]
