#' Multi-Value assigne operator.
#'
#' @param lhs vector of variables (or variable names) to assign values.
#' @param rhs object that can be coersed to list of the same length as
#'  \code{lhs} storing the values of the variables defined in \code{lhs}.
#'
#' @details The parameter names \code{lhs} and \code{rhs} stand for "Left Hand
#'  Side" and "Right Hand Side", respectively.
#'
#' @examples
#' c(a, b) %<-% list(1, 2)
#' # is equivalent to
#' ## a <- 1
#' ## b <- 2
#'
#' # Switching the values of a, b could be done by.
#' c(a, b) %<-% list(b, a)
#' note the usage of 'list' on the right side, otherwise an exraction of the
#' first two values of the concatenated object is performed. See next:
#'
#' # Extract values.
#' c(d1, d2, d3) %<-% 1:10
#' extracting the first three valus from the vector of length 10.
#'
"%<-%" <- function(lhs, rhs) {
    var.names <- make.names(as.list(substitute(lhs))[-1])
    values <- as.list(rhs)
    env <- parent.frame()
    for (i in seq_along(var.names)) {
        assign(var.names[i], values[[i]], envir = env)
    }
}
