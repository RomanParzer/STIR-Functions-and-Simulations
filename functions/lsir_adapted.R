source('../functions/matpow.R')

#' Longitudinal Sliced Inverse Regression
#'
#' @param X array of dim \eqn{t \times p \times n} with each \eqn{t \times p} matrix  representing a
#'  observation (assumed to be centered over observations).
#' @param y vector of \eqn{n} elements as factors. (can be coersed to factors)
#' @param p,t,k,r dimensions.
#'
#' @returns a list with components
#'  alpha: matrix of \eqn{t \times r}
#'  beta: matrix of \eqn{p \times k}
#'  X_fitted: \eqn{t \times p \times n} array of fitted values
#'
LSIR <- function(X, y, p, t, m = 1L, r = 1L) {
    # Check and transform parameters.
    if (!is.array(X)) X <- as.array(X)
    n <- dim(X)[3]
    stopifnot(
        dim(X)[1] == t,
        dim(X)[2] == p,
        n == length(y)
    )
    if (!is.factor(y)) y <- factor(y)
    
    X_perm <- aperm(X,c(3,2,1))
    # the code assumes:
    # alpha: T x r, beta: p x k, X_i: p x T, for ith observation

    # Restructure X into a 3D tensor with axis (observations, predictors, time).
    dim(X_perm) <- c(n, p, t)

    # Estimate predictor/time covariance matrices \hat{Sigma}_1, \hat{Sigma}_2.
    sigma_p <- matrix(rowMeans(apply(X_perm, 3, cov)), p, p)
    sigma_t <- matrix(rowMeans(apply(X_perm, 2, cov)), t, t)

    # Normalize X as vec(Z) = Sigma^-1/2 (vec(X) - E(vec(X)))
    dim(X_perm) <- c(n, p * t)
    sigma_p_isqrt <- matpow(sigma_p, -0.5)
    sigma_t_isqrt <- matpow(sigma_t, -0.5)
    Z <- scale(X_perm, scale = FALSE) %*% kronecker(sigma_t_isqrt, sigma_p_isqrt)
    # Both as 3D tensors.
    dim(X_perm) <- dim(Z) <- c(n, p, t)

    # Estimate the conditional predictor/time covariance matrix Omega = cov(E(Z|Y)).
    omega_p <- matrix(Reduce(`+`, lapply(levels(y), function(l) {
        rowMeans(apply(Z[y == l, , ], 3, function(z) {
            (nrow(z) / n) * tcrossprod(colMeans(z))
        }))
    })), p, p)
    omega_t <- matrix(Reduce(`+`, lapply(levels(y), function(l) {
        rowMeans(apply(Z[y == l, , ], 2, function(z) {
            (nrow(z) / n) * tcrossprod(colMeans(z))
        }))
    })), t, t)

    # Compute seperate SVD of estimated omega's and use that for an estimate of
    # a central subspace basis.
    svd_p <- La.svd(omega_p,nu=m,nv=0)
    svd_t <- La.svd(omega_t,nu=r,nv=0)
    beta  <- as.matrix(sigma_p_isqrt %*% svd_p$u)
    alpha <- as.matrix(sigma_t_isqrt %*% svd_t$u)
    
    Bhat <- as.matrix(kronecker(beta, alpha))
    # alpha <- alpha*beta[1]
    # beta <- beta/beta[1]
    delta_p <- sigma_p
    delta_t <- sigma_t
    
    return(list(sigma_p = sigma_p, sigma_t = sigma_t,
                sigma = kronecker(sigma_p, sigma_t),
                alpha = alpha, beta = beta,delta_p=delta_p,delta_t=delta_t,
                B = Bhat,m=m,r=r))
}
