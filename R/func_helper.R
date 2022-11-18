#' if.null
#'
#' This function is wrapper for ifelse(is.null(.),.,.)
if.null <- function(vec, val) {
  ifelse(is.null(vec), val, vec)
}

alt_chol <- function(S) {
  eig_decomp <- eigen(S)
  eig_decomp$values <- Re(ifelse(Re(eig_decomp$values) < 0, 0, eig_decomp$values))
  return(eig_decomp$vectors %*% diag(sqrt(eig_decomp$values)) %*% t(eig_decomp$vectors))
}
