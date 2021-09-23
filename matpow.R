#' Generalized matrix power function for symmetric matrices.
#'
#' Using the SVD of the matrix \eqn{A = U D V'} where \eqn{D} is the
#' diagonal matrix with the singular values of \eqn{A}, the powers are then
#' computed as \deqn{A^p = U D^p V'} using the symmetry of \eqn{A}.
#'
#' @details
#' It is assumed that the argument \code{A} is symmetric and it is not checked
#' for symmerie, the result will just be wrong. The actual formula for negative
#' powers is \deqn{A^-p = V D^-p U'}. For symmetric matrices \eqn{U = V} which
#' gives the formula used by this function.
#' The reason is for speed, using the symmetrie propertie as described avoids
#' two transpositions in the algorithm (one in \code{svd} using \code{La.svd}).
#'
#' @param A Matrix.
#' @param pow numeric power.
#' @param tol relative tolerance to detect zero singular values as well as the
#'  \code{qr} factorization tolerance.
#'
#' @return a matrix.
#'
#' @seealso \code{\link{solve}}, \code{\link{qr}}, \code{\link{svd}}.
#'
#' @examples
#' # General full rank square matrices.
#' A <- matrix(rnorm(121), 11, 11)
#' all.equal(matpow(A,  1), A)
#' all.equal(matpow(A,  0), diag(nrow(A)))
#' all.equal(matpow(A, -1), solve(A))
#'
#' # Roots of full rank symmetric matrices.
#' A <- crossprod(A)
#' B <- matpow(A, 0.5)
#' all.equal(B %*% B, A)
#' all.equal(matpow(A, -0.5), solve(B))
#' C <- matpow(A, -0.5)
#' all.equal(C %*% C %*% A, diag(nrow(A)))
#' all.equal(A %*% C %*% C, diag(nrow(A)))
#'
#' # General singular matrices.
#' A <- matrix(rnorm(72), 12, 12) # rank(A) = 6
#' B <- matpow(A, -1) # B = A^+
#' # Check generalized inverse properties.
#' all.equal(A %*% B %*% A, A)
#' all.equal(B %*% A %*% B, B)
#' all.equal(B %*% A, t(B %*% A))
#' all.equal(A %*% B, t(A %*% B))
#'
#' # Roots of singular symmetric matrices.
#' A <- crossprod(matrix(rnorm(72), 12, 12)) # rank(A) = 6
#' B <- matpow(A, -0.5) # B = (A^+)^1/2
#' # Check generalized inverse properties.
#' all.equal(A %*% B %*% B %*% A, A)
#' all.equal(B %*% B %*% A %*% B %*% B, B %*% B)
#' all.equal(B %*% A, t(B %*% A))
#' all.equal(A %*% B, t(A %*% B))
#'
matpow <- function(A, pow, tol = 1e-7) {
    if (nrow(A) != ncol(A)) {
        stop("Expected a square matix, but 'A' is ", nrow(A), " by ", ncol(A))
    }
    # Case study for negative, zero or positive power.
    if (pow > 0) {
        if (pow == 1) { return(A) }
        # Perform SVD and return power as A^pow = U diag(d^pow) V'.
        svdA <- La.svd(A)
        return(svdA$u %*% ((svdA$d^pow) * svdA$vt))
    } else if (pow == 0) {
        return(diag(nrow(A)))
    } else {
        # make QR decomposition.
        qrA <- qr(A, tol = tol)
        # Check rank of A.
        if (qrA$rank == nrow(A)) {
            # Full rank, calc inverse the classic way using A's QR decomposition
            return(matpow(solve(qrA), abs(pow), tol = tol))
        } else {
            # For singular matrices use the SVD decomposition for the power
            svdA <- svd(A)
            # Get (numerically) positive singular values.
            positives <- svdA$d > tol * svdA$d[1]
            # Apply the negative power to positive singular values and augment
            # the rest with zero.
            d <- c(svdA$d[positives]^pow, rep(0, sum(!positives)))
            # The pseudo invers as A^pow = V diag(d^pow) U' for pow < 0.
            return(svdA$v %*% (d * t(svdA$u)))
        }
    }
}
