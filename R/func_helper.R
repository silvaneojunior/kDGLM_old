#' if.null
#'
#' This function is wrapper for ifelse(is.null(.),.,.)
if.null <- function(vec, val) {
  ifelse(is.null(vec), val, vec)
}

#' alt_chol
#'
alt_chol <- function(S) {
  eig_decomp <- eigen(S)
  eig_decomp$values <- Re(ifelse(Re(eig_decomp$values) < 0, 0, eig_decomp$values))
  return(eig_decomp$vectors %*% diag(sqrt(eig_decomp$values)) %*% t(eig_decomp$vectors))
}

#' colProd
#'
colProd <- function(x) {
  apply(x, 2, prod)
}

#' calcula_max
#'
#' Auxiliary function to calculate the axis limits and gradation for plots.
#'
#' @param pre_max Numeric: A vector/matrix from which to calculate the axis limits and gradation.
#'
#' @return A list contaning the gradation for the axis, the number of ticks in the axis and the maximum value.
calcula_max <- function(pre_max) {
  if (length(pre_max) == 0 | sum(pre_max**2) < 10**-20) {
    pre_max <- 1
  } else {
    pre_max <- max(pre_max)
  }
  scaled_max <- log10(pre_max)
  category <- scaled_max %% 1
  value <- 10**(floor(log10(max(pre_max))))
  if (category < 0.1) {
    value <- value / 10
  } else {
    if (category < 0.25) {
      value <- value / 5
    } else {
      if (category < 0.5) {
        value <- value / 2
      }
    }
  }
  interval_size <- (pre_max %/% value) + 2
  max_value <- value * interval_size

  return(list(value, interval_size, max_value))
}
