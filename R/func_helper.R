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

#' colQuantile
#'
#' A function that calculates the column-wise quantiles of a matrix.
#'
#' @param X Matrix.
#' @param q Numeric: A number between 0 and 1.
#'
#' @importFrom Rfast colnth
#'
#' @return Vector: The chosen quatile for each column of X.
colQuantile <- function(X, q) {
  n <- dim(X)[1]
  k <- dim(X)[2]
  min_index <- floor(n * q)
  max_index <- ceiling(n * q)
  (colnth(X, rep(min_index, k)) + colnth(X, rep(max_index, k))) / 2
}

#' rowQuantile
#'
#' A function that calculates the row-wise quantiles of a matrix.
#'
#' @param X Matrix.
#' @param q Numeric: A number between 0 and 1.
#'
#' @importFrom Rfast rownth
#'
#' @return Vector: The chosen quatile for each row of X.
rowQuantile <- function(X, q) {
  n <- dim(X)[1]
  k <- dim(X)[2]
  min_index <- floor(k * q)
  max_index <- ceiling(k * q)
  (rownth(X, rep(min_index, n)) + rownth(X, rep(max_index, n))) / 2
}


#' f_root
#'
#' Calculates the root of a function given an initial value and a function to calculate it's derivatives.
#'
#' @param f function: A function that receives a vector and return a vector of the same size.
#' @param df function: A function that receives a vector and return the derivatives of fx in relation to it's argmuments (must return a matrix, if fx returns a vector).
#' @param start vector: The initial value to start the algorithm.
#' @param tol numeric: The tolerance for the solution error.
#' @param n_max numeric: The maximum number of iterations allowed.
#'
#' @return A list contaning:
#' \itemize{
#'    \item root vector: The solution for the system.
#'    \item f.root vector: The function fx evaluated at the root.
#'    \item iter numeric: The number of steps taken.
#' }
f_root <- function(f, df, start, tol = 1e-8, n_max = 1000) {
  x_root <- start
  fx <- f(x_root)
  dfx <- df(x_root)
  error <- max(abs(fx))
  count <- 0
  while (error >= tol & count < n_max) {
    count <- count + 1
    change <- solve(dfx, -fx)
    x_root <- x_root + change
    fx <- f(x_root)
    dfx <- df(x_root)
    error <- max(abs(fx))
  }
  if (count >= n_max) {
    warning("Steady state not reached.\n")
  }
  return(list("root" = x_root, "f.root" = fx, "inter." = count))
}
