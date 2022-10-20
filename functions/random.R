#' Multivariate Normal Distribution.
#'
#' Random generation for the multivariate normal distribution.
#' \deqn{X \sim N_p(\mu, \Sigma)}{X ~ N_p(\mu, \Sigma)}
#'
#' @param n number of samples.
#' @param mu mean
#' @param sigma covariance matrix.
#'
#' @return a \eqn{n\times p}{n x p} matrix with samples in its rows.
#'
#' @examples
#' \dontrun{
#' rmvnorm(20, sigma = matrix(c(2, 1, 1, 2), 2))
#' rmvnorm(20, mu = c(3, -1, 2))
#' }
#' @keywords internal
rmvnorm <- function(n = 1, mu = rep(0, p), sigma = diag(p)) {
    if (!missing(sigma)) {
        p <- nrow(sigma)
    } else if (!missing(mu)) {
        mu <- matrix(mu, ncol = 1)
        p <- nrow(mu)
    } else {
        stop("At least one of 'mu' or 'sigma' must be supplied.")
    }

    # See: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
    return(rep(mu, each = n) + matrix(rnorm(n * p), n) %*% chol(sigma))
}
