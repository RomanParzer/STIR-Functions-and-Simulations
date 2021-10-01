#source('./random.R')
source('../functions/STIR.R')
source('../functions/multi_assign.R')
source('../functions/lsir_adapted.R')


#require(dplyr)
#require(stats)
require(MASS) # for ginv and mvrnorm
require(Matrix)
require(mvtnorm)
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

simulateSDR_wald <- function(nsim, n, B, t, p, rho,H=1, d=1,what_time = "1:t/t",cvec = 0:5/5,delta_y=FALSE) {
  len_c <- length(cvec) 
  res <- NULL
  rej_rates <- list(STIR=numeric(len_c),lm=numeric(len_c),lm2=numeric(len_c))
  rkwald <- wald <- myF <- myF2 <- matrix(c(0),len_c,nsim)
  rkS <- NULL
  rkVCB <- matrix(c(0),len_c,nsim)
  
  for (c_count in 1:len_c) {
    rej_count <- list(STIR=numeric(nsim),lm=numeric(nsim),lm2=numeric(nsim))
    
    scaledB <- B
    scaledB[,3] <- cvec[c_count]*B[,3]
    
    for (i in 1:nsim) {
      if ((10*i)%%nsim == 0) {
        cat(sprintf('Starting rep %d/%d for c %d/%d \n', i, nsim,c_count,len_c))
      }
      
      set.seed(i+1000)
      datas <- generateDataSTIR(n,scaledB,t,p,rho,basis_S = "Polynomial",H=H,d=d,what_time =what_time,ind_tpts=FALSE, delta_y=delta_y)
      res <- STIR(X=datas$X, y=datas$y, t=t, p=p, d=d,basis_S = "Polynomial",H=H,cont_y = TRUE,k=NULL,time_pts = datas$time_pts,calc_VC_vecB=TRUE)
      #A3 <- cbind(matrix(0,d*H,d*H*(2)),diag(1,d*H),matrix(0,d*H,d*H*(p-3)))
      rkwald[c_count,i] <- rankMatrix(res$VC_vecB[2*d*H+1:(d*H),2*d*H+1:(d*H)])
      
      if (i==1 && c_count==1) {rkS <- rankMatrix(datas$S)}
      rkVCB[c_count,i] <- rankMatrix(res$VC_vecB)
      
      if (rkwald[c_count,i] < d*H) {
        wald[c_count,i] <- crossprod(res$B[,3],ginv(res$VC_vecB[2*d*H+1:(d*H),2*d*H+1:(d*H)])%*%res$B[,3])
      } else {
        wald[c_count,i] <- crossprod(res$B[,3],solve(res$VC_vecB[2*d*H+1:(d*H),2*d*H+1:(d*H)],res$B[,3]))
      }
      if (pchisq(wald[c_count,i],rkwald[c_count,i],lower.tail = FALSE) <= 0.05) {
        rej_count$STIR[i] <- 1
      }
      X <- datas$X
      dim(X) <- c(t*p,n)
      lm_fit <- lm(y~.,data = data.frame(t(X),y=datas$y))
      
      lmX <- model.matrix(lm_fit)
      XXinv <- solve(crossprod(lmX))
      sub_delta <- res$delta_ls[2*t+1:t,]
      RXXR <- tcrossprod(sub_delta%*%XXinv[-1,-1],sub_delta)
      myF[c_count,i] <- crossprod(sub_delta%*%coef(lm_fit)[-1],solve(RXXR,sub_delta%*%coef(lm_fit)[-1])) / (t*summary(lm_fit)$sigma^2)
      if (pf(myF[c_count,i],t,n-t*p-1,lower.tail=FALSE) <= 0.05) {
        rej_count$lm[i] <- 1
      }
      if (delta_y) {
        sub_delta2 <- (datas$delta* mean(exp(-min(sqrt(10)*abs(datas$y),2))) )[2*t+1:t,]
      } else {
        sub_delta2 <- datas$delta[2*t+1:t,]
      }
      
      RXXR <- tcrossprod(sub_delta2%*%XXinv[-1,-1],sub_delta2)
      myF2[c_count,i] <- crossprod(sub_delta2%*%coef(lm_fit)[-1],solve(RXXR,sub_delta2%*%coef(lm_fit)[-1])) / (t*summary(lm_fit)$sigma^2)
      if (pf(myF2[c_count,i],t,n-t*p-1,lower.tail=FALSE) <= 0.05) {
        rej_count$lm2[i] <- 1
      }
      
    }
    rej_rates$STIR[c_count] <- mean(rej_count$STIR)
    rej_rates$lm[c_count] <- mean(rej_count$lm)
    rej_rates$lm2[c_count] <- mean(rej_count$lm2)
  }
  return(list(what_time=what_time, cvec = cvec , rej_rates=rej_rates,
              rkwald = rkwald, wald=wald,rkS=rkS, 
              rkVCB=rkVCB,myF=myF))
}

## dH >= T,p !!! 
## dH >= rk(B) !!! 
myB12 <- matrix(c(rep(c(1,0.5),2),rep(c(0.7,0.3),2),rep(0.1,4)),4,3)
myB3 <- matrix(c(rep(c(1,0.5),5),rep(c(0.7,0.3),5),rep(0.1,10)),10,3)
myB4 <- matrix(c(rep(c(1,0.5,1),1),rep(c(0.7,0.3,0.7),1),rep(0.1,3)),3,3)
myB5 <- 2*matrix(c(1,  0,1/2,  0,1/2,  0,1,
                   0,1/4,1/2,1/4,1/2,1/4,0),2,7,byrow = TRUE)
myB5flipped <- 2*matrix(c(0,1/4,1/2,1/4,1/2,1/4,0,
                          1,  0,1/2,  0,1/2,  0,1),2,7,byrow = TRUE)
myB5ad <- 2*matrix(c(1,  0,1/2,  0,1/2,  0,1,
                     1,1/4,1/2,1/4,1/2,1/4,0),2,7,byrow = TRUE)
myB12flipped <- matrix(c(rep(c(0.5,1),2),rep(c(0.7,0.3),2),rep(0.1,4)),4,3)
#                       n,         B, t, p, H, d, rho,  time, ind_tpts, delta_y
#              --------------------------------------------------------
params <- list( list( 500,myB12,10, 3, 2, 2,  0.8,"exp(1:t/6)/exp(t/6)", 0:5/2,0:5, TRUE),
  list( 500, myB12,10, 3, 2, 2,  0.8, "exp(1:t/6)/exp(t/6)", 0:5/2,0:5, FALSE),
  list( 2000,myB12,10, 3, 2, 2,  0.8,"exp(1:t/6)/exp(t/6)", 0:5/2,0:5, TRUE),
  list( 2000, myB12,10, 3, 2, 2,  0.8, "exp(1:t/6)/exp(t/6)", 0:5/2,0:5, FALSE),
  list( 500,       myB5,2, 7, 2, 1,  0.8,     "1:t/t", 0:5/5,0:5/5, TRUE),
  list( 500,       myB5,2, 7, 2, 1,  0.8,     "1:t/t", 0:5/5,0:5/5, FALSE),
  list( 2000,       myB5,2, 7, 2, 1,  0.8,     "1:t/t", 0:5/5,0:5/5, TRUE),
  list( 2000,       myB5,2, 7, 2, 1,  0.8,     "1:t/t", 0:5/5,0:5/5, FALSE)
)


# what_time <- "exp(1:t/6)/exp(t/6)"
# # what_time <- "1:t/t"
# npred <- 100
# nlev_ydis <- 10
# c_count <- 1
# j <- 10
# i <- 1
# traceback()
# warnings()

# no k or rm_list needed
## do 1 and 4
nsim <- 500
counter <- 1
nsettings <- length(params)
for (param in params) {
  cat(sprintf('Parameter setting %d / %d:\n',counter,nsettings))
  counter <- counter + 1
  c(n,B,t,p,H,d,rho,time_sel,cvec,cvec2,delta_y) %<-% param
  sim <- simulateSDR_wald(nsim, n, B, t, p, rho,H, d,what_time = time_sel,cvec = cvec,delta_y = delta_y)

  saveRDS(sim, file = sprintf("../results/wald_contSDR_%d_%d_%d_%d_%d_%d_%d_%.2f_time_%.1s_delta_y_%s_maxc_%.2f_3rdMarker_delta_beta.rds",
                            nsim, n, round(sum(abs(B))), t, p, H, d, rho,time_sel,delta_y,max(cvec)))
}

