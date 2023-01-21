#' Creates the structure for a polynomial block with desired order.
#'
#' @param ... Named values for the planing matrix.
#' @param order Positive integer: The order of the polimial structure.
#' @param name String: An optional argument providing the name for this block. Can be useful to identify the models with meaningful labels, also, the name used will be used in some auxiliary functions.
#' @param D Array, Matrix, vector or  scalar: The values for the discount factors at each time. If D is a array, it's dimensions should be nxnxt, where n is the order of the polynomial block and t is the length of the outcomes. If D is a matrix, it's dimesions should be nxn and it's values will be used for each time. If D is a vector or scalar, a discount factor matrix will be created as a diagonal matrix with the values of D in the diagonal.
#' @param W Array, Matrix, vector or  scalar: The values for the covariance matrix for the noise factor at each time. If W is a array, it's dimensions should be nxnxt, where n is the order of the polynomial block and t is the length of the outcomes. If W is a matrix, it's dimesions should be nxn and it's values will be used for each time. If W is a vector or scalar, a discount factor matrix will be created as a diagonal matrix with the values of W in the diagonal.
#' @param m0 Vector or scalar: The prior mean for the latent variables associated with this block. If m0 is a vector, it's dimesion should be equal to the order of the polynomial block. If m0 is a scalar, it's value will be used for all latent variables.
#' @param C0 Matrix, vector or scalar: The prior covariance matrix for the latent variables associated with this block. If C0 is a matrix, it's dimesions should be nxn. If W is a vector or scalar, a covariance matrix will be created as a diagonal matrix with the values of C0 in the diagonal.
#'
#' @return An object of the class dlm_block containing the following values:
#' \itemize{
#'    \item FF Array: A 3D-array containing the regression matrix for each time. It's dimension should be n x m x T, where n is the number of latent variables, m is the number of outcomes in the model and T is the time series length.
#'    \item G Matrix: The state evolution matrix.
#'    \item D Array: A 3D-array containing the discount factor matrix for each time. It's dimension should be n x n x T, where n is the number of latent variables and T is the time series length.
#'    \item W Array: A 3D-array containing the covariance matrix of the noise for each time. It's dimension should be the same as D.
#'    \item m0 Vector: The prior mean for the latent vector.
#'    \item C0 Matrix: The prior covariance matrix for the latent vector.
#'    \item names list: A list containg the variables indexes by their name.
#'    \item order Positive integer: Same as argument.
#'    \item n Positive integer: The number of latent variables associated with this block (same value as order).
#'    \item t Positive integer: The number of time steps associated with this block. If 1, the block is compatible with blocks of any time length, but if t is greater than 1, this block can only be used with blocks of the same time length.
#'    \item k Positive integer: The number of outcomes associated with this block. This block can only be used with blocks with the same outcome length.
#'    \item var_names Vector: The name of the linear predictors associated with this block,
#'    \item type Character: The type of block (polynomial).
#' }
#'
#' @export
#' @examples
#' # Creating a first order structure for a model with 2 outcomes.
#' # One block is created for each outcome, with each block being associated with only one of the outcomes.
#' level_1 <- polynomial_block(alpha1 = 1, order = 1)
#' level_2 <- polynomial_block(alpha2 = 1, order = 1)
#'
#' # Creating a block with shared effect between the oucomes
#' level_3 <- polynomial_block(alpha1 = 1, alpha2 = 1, order = 2)
polynomial_block <- function(..., order = 1, name = "Var_Poly", D = 1, W = 0, m0 = 0, C0 = 1) {
  values <- list(...)
  vars <- names(values)
  if (any(vars == "")) {
    stop("Error: One or more outcomes are unnamed. Please, name all arguments.")
  }
  var_len <- sapply(values, length)
  k <- length(values)
  t <- max(var_len)
  if (any(var_len != t & var_len > 1)) {
    stop(paste0("Error: Outcomes have mismatching lengths. Expected 1 or ", t, ", got: ", paste(var_len, collapse = ", "), "."))
  }

  if (length(D) == 1) {
    D <- array(1, c(order, order, t)) * D
  } else if (is.vector(D)) {
    D <- array(D, c(length(D), order, order)) %>% aperm(c(3, 2, 1))
  } else if (is.matrix(D)) {
    D <- array(D, c(dim(D)[1], dim(D)[2], t))
  }
  t <- if (t == 1) {
    dim(D)[3]
  } else {
    t
  }

  if (length(dim(D)) > 3 | any(dim(D)[1:2] != order) | (dim(D)[3] != t & t > 1)) {
    stop(paste0("Error: Invalid shape for D. Expected ", order, "x", order, "x", t, ". Got ", paste(dim(D), collapse = "x"), "."))
  }

  if (length(W) == 1) {
    W <- array(diag(order), c(order, order, t)) * W
  } else if (is.vector(W)) {
    W_vals <- W
    W <- array(diag(order), c(order, order, length(W_vals)))
    for (i in 1:length(W_vals)) {
      W[, , i] <- W[, , i] * W_vals[i]
    }
  } else if (is.matrix(W)) {
    W <- array(W, c(dim(W)[1], dim(W)[2], t))
  }
  t <- if (t == 1) {
    dim(W)[3]
  } else {
    t
  }
  D <- array(D, c(order, order, t))

  if (length(dim(W)) > 3 | any(dim(W)[1:2] != order) | (dim(W)[3] != t & t > 1)) {
    stop(paste0("Error: Invalid shape for W. Expected ", order, "x", order, "x", t, ". Got ", paste(dim(W), collapse = "x"), "."))
  }

  for (name_var in vars[var_len == 1]) {
    var_len[[name_var]] <- t
    values[[name_var]] <- rep(values[[name_var]], t)
  }

  FF <- array(0, c(order, k, t))
  FF[1, , ] <- matrix(sapply(values, c), k, t, byrow = TRUE)
  D[, , apply(is.na(FF), 3, any)] <- 1
  W[, , apply(is.na(FF), 3, any)] <- 0
  FF <- ifelse(is.na(FF), 0, FF)

  G <- diag(order)
  if (order == 2) {
    G[1, 2] <- 1
  } else if (order > 2) {
    diag(G[1:(order - 1), 2:order]) <- 1
  }

  m0 <- if (length(m0) == 1) {
    rep(m0, order)
  } else {
    m0
  }
  C0 <- if (length(C0) == 1) {
    diag(order) * C0
  } else if (is.vector(C0)) {
    diag(C0)
  } else {
    C0
  }

  if (length(dim(C0)) > 2) {
    stop(paste0("Error: C0 must be a matrix, but it has ", length(dim(C0)), " dimensions."))
  }
  if (any(dim(C0) != order)) {
    stop(paste0("Error: C0 must have dimensions ", order, "x", order, ". Got ", dim(C0)[1], "x", dim(C0)[2], "."))
  }

  names <- list()
  names[[name]] <- c(1:order)
  block <- list(
    "FF" = FF,
    "G" = G,
    "D" = D,
    "W" = W,
    "m0" = m0,
    "C0" = C0,
    "names" = names,
    "order" = order,
    "n" = order,
    "t" = t,
    "k" = k,
    "var_names" = vars,
    "type" = "Polynomial"
  )
  class(block) <- "dlm_block"
  return(block)
}


#' harmonic_block
#'
#' Creates the structure for a harmonic block with desired periodicity.
#'
#' @param ... Named values for the planing matrix.
#' @param period Positive integer: The size of the harmonic cycle.
#' @param name String: An optional argument providing the name for this block. Can be useful to identify the models with meaningful labels, also, the name used will be used in some auxiliary functions.
#' @param D Array, Matrix, vector or  scalar: The values for the discount factors at each time. If D is a array, it's dimensions should be nxnxt, where n is the order of the polynomial block and t is the length of the outcomes. If D is a matrix, it's dimesions should be nxn and it's values will be used for each time. If D is a vector or scalar, a discount factor matrix will be created as a diagonal matrix with the values of D in the diagonal.
#' @param W Array, Matrix, vector or  scalar: The values for the covariance matrix for the noise factor at each time. If W is a array, it's dimensions should be nxnxt, where n is the order of the polynomial block and t is the length of the outcomes. If W is a matrix, it's dimesions should be nxn and it's values will be used for each time. If W is a vector or scalar, a discount factor matrix will be created as a diagonal matrix with the values of W in the diagonal.
#' @param m0 Vector or scalar: The prior mean for the latent variables associated with this block. If m0 is a vector, it's dimesion should be equal to the order of the polynomial block. If m0 is a scalar, it's value will be used for all latent variables.
#' @param C0 Matrix, vector or scalar: The prior covariance matrix for the latent variables associated with this block. If C0 is a matrix, it's dimesions should be nxn. If W is a vector or scalar, a covariance matrix will be created as a diagonal matrix with the values of C0 in the diagonal.
#'
#' @return An object of the class dlm_block containing the following values:
#' \itemize{
#'    \item FF Array: A 3D-array containing the regression matrix for each time. It's dimension should be n x m x T, where n is the number of latent variables, m is the number of outcomes in the model and T is the time series length.
#'    \item G Matrix: The state evolution matrix.
#'    \item D Array: A 3D-array containing the discount factor matrix for each time. It's dimension should be n x n x T, where n is the number of latent variables and T is the time series length.
#'    \item W Array: A 3D-array containing the covariance matrix of the noise for each time. It's dimension should be the same as D.
#'    \item m0 Vector: The prior mean for the latent vector.
#'    \item C0 Matrix: The prior covariance matrix for the latent vector.
#'    \item names list: A list containg the variables indexes by their name.
#'    \item period Positive integer: Same as argument.
#'    \item n Positive integer: The number of latent variables associated with this block (2).
#'    \item t Positive integer: The number of time steps associated with this block. If 1, the block is compatible with blocks of any time length, but if t is greater than 1, this block can only be used with blocks of the same time length.
#'    \item k Positive integer: The number of outcomes associated with this block. This block can only be used with blocks with the same outcome length.
#'    \item var_names Vector: The name of the linear predictors associated with this block,
#'    \item type Character: The type of block (Harmonic).
#' }
#'
#' @export
#' @examples
#' # Creating seasonal structure for a model with 2 outcomes.
#' # One block is created for each outcome, with each block being associated with only one of the outcomes.
#' season_1 <- harmonic_block(alpha1 = 1, period = 3)
#' season_2 <- harmonic_block(alpha2 = 1, period = 6)
#'
#' # Creating a block with shared effect between the oucomes
#' season_3 <- harmonic_block(alpha = 1, alpha2 = 1, period = 12)
harmonic_block <- function(..., period, name = "Var_Sazo", D = 1, W = 0, m0 = 0, C0 = 1) {
  w <- 2 * pi / period
  order <- 2
  block <- polynomial_block(..., order = order, name = name, D = D, W = W, m0 = m0, C0 = C0)

  G <- matrix(c(cos(w), -sin(w), sin(w), cos(w)), order, order)
  block$G <- G
  block$order <- NULL
  block$period <- period
  block$type <- "Harmonic"
  return(block)
}

#' AR_block
#'
#' DESCRIPTION
#'
#' @param ... Named values for the planing matrix.
#' @param order Positive integer: DESCRIPTION.
#' @param name String: An optional argument providing the name for this block. Can be useful to identify the models with meaningful labels, also, the name used will be used in some auxiliary functions.
#' @param D Array, Matrix, vector or  scalar: The values for the discount factors at each time. If D is a array, it's dimensions should be nxnxt, where n is the order of the polynomial block and t is the length of the outcomes. If D is a matrix, it's dimesions should be nxn and it's values will be used for each time. If D is a vector or scalar, a discount factor matrix will be created as a diagonal matrix with the values of D in the diagonal.
#' @param W Array, Matrix, vector or  scalar: The values for the covariance matrix for the noise factor at each time. If W is a array, it's dimensions should be nxnxt, where n is the order of the polynomial block and t is the length of the outcomes. If W is a matrix, it's dimesions should be nxn and it's values will be used for each time. If W is a vector or scalar, a discount factor matrix will be created as a diagonal matrix with the values of W in the diagonal.
#' @param m0 Vector or scalar: The prior mean for the latent variables associated with this block. If m0 is a vector, it's dimesion should be equal to the order of the polynomial block. If m0 is a scalar, it's value will be used for all latent variables.
#' @param C0 Matrix, vector or scalar: The prior covariance matrix for the latent variables associated with this block. If C0 is a matrix, it's dimesions should be nxn. If W is a vector or scalar, a covariance matrix will be created as a diagonal matrix with the values of C0 in the diagonal.
#'
#' @return An object of the class dlm_block containing the following values:
#' \itemize{
#'    \item FF Array: A 3D-array containing the regression matrix for each time. It's dimension should be n x m x T, where n is the number of latent variables, m is the number of outcomes in the model and T is the time series length.
#'    \item G Matrix: The state evolution matrix.
#'    \item D Array: A 3D-array containing the discount factor matrix for each time. It's dimension should be n x n x T, where n is the number of latent variables and T is the time series length.
#'    \item W Array: A 3D-array containing the covariance matrix of the noise for each time. It's dimension should be the same as D.
#'    \item m0 Vector: The prior mean for the latent vector.
#'    \item C0 Matrix: The prior covariance matrix for the latent vector.
#'    \item names list: A list containg the variables indexes by their name.
#'    \item order Positive integer: Same as argument.
#'    \item n Positive integer: The number of latent variables associated with this block (2).
#'    \item t Positive integer: The number of time steps associated with this block. If 1, the block is compatible with blocks of any time length, but if t is greater than 1, this block can only be used with blocks of the same time length.
#'    \item k Positive integer: The number of outcomes associated with this block. This block can only be used with blocks with the same outcome length.
#'    \item var_names Vector: The name of the linear predictors associated with this block,
#'    \item type Character: The type of block (AR).
#' }
#'
#' @export
#' @examples
#' # EXAMPLE
AR_block <- function(..., order, name = "Var_AR", D = 1, W = 0, m0 = 0, C0 = 1) {
  block <- polynomial_block(..., order = 2 * order, name = name, D = D, W = W, m0 = m0, C0 = C0)
  k <- block$k

  if (order == 1) {
    G <- diag(2)
    G[1, 1] <- NA
  } else {
    G <- matrix(0, 2 * order, 2 * order)
    G[1, 1:order] <- NA
    G[2:order, -(order:(2 * order))] <- diag(order - 1)
    G[(order + 1):(2 * order), (order + 1):(2 * order)] <- diag(order)
    index <- sort(c(c(1:order), c(1:order)))
    G <- G[index + c(0, order), index + c(0, order)]
  }
  block$G <- G
  block$order <- order
  block$type <- "AR"
  return(block)
}

#' correlation_block
#'
#' DESCRIPTION
#'
#' @param var_names Vector: Name of the linear predictors associated with this block.
#' @param name String: An optional argument providing the name for this block. Can be useful to identify the models with meaningful labels, also, the name used will be used in some auxiliary functions.
#' @param D Array, Matrix, vector or  scalar: The values for the discount factors at each time. If D is a array, it's dimensions should be nxnxt, where n is the order of the polynomial block and t is the length of the outcomes. If D is a matrix, it's dimesions should be nxn and it's values will be used for each time. If D is a vector or scalar, a discount factor matrix will be created as a diagonal matrix with the values of D in the diagonal.
#' @param W Array, Matrix, vector or  scalar: The values for the covariance matrix for the noise factor at each time. If W is a array, it's dimensions should be nxnxt, where n is the order of the polynomial block and t is the length of the outcomes. If W is a matrix, it's dimesions should be nxn and it's values will be used for each time. If W is a vector or scalar, a discount factor matrix will be created as a diagonal matrix with the values of W in the diagonal.
#' @param m0 Vector or scalar: The prior mean for the latent variables associated with this block. If m0 is a vector, it's dimesion should be equal to the order of the polynomial block. If m0 is a scalar, it's value will be used for all latent variables.
#' @param C0 Matrix, vector or scalar: The prior covariance matrix for the latent variables associated with this block. If C0 is a matrix, it's dimesions should be nxn. If W is a vector or scalar, a covariance matrix will be created as a diagonal matrix with the values of C0 in the diagonal.
#'
#' @return An object of the class dlm_block containing the following values:
#' \itemize{
#'    \item FF Array: A 3D-array containing the regression matrix for each time. It's dimension should be n x m x T, where n is the number of latent variables, m is the number of outcomes in the model and T is the time series length.
#'    \item G Matrix: The state evolution matrix.
#'    \item D Array: A 3D-array containing the discount factor matrix for each time. It's dimension should be n x n x T, where n is the number of latent variables and T is the time series length.
#'    \item W Array: A 3D-array containing the covariance matrix of the noise for each time. It's dimension should be the same as D.
#'    \item m0 Vector: The prior mean for the latent vector.
#'    \item C0 Matrix: The prior covariance matrix for the latent vector.
#'    \item names list: A list containg the variables indexes by their name.
#'    \item n Positive integer: The number of latent variables associated with this block (2).
#'    \item t Positive integer: The number of time steps associated with this block. If 1, the block is compatible with blocks of any time length, but if t is greater than 1, this block can only be used with blocks of the same time length.
#'    \item k Positive integer: The number of outcomes associated with this block. This block can only be used with blocks with the same outcome length.
#'    \item var_names Vector: The name of the linear predictors associated with this block,
#'    \item type Character: The type of block (Correlation).
#' }
#'
#' @importFrom Matrix bdiag
#'
#' @export
#' @examples
#' # EXAMPLE
correlation_block <- function(var_names, name = "Var_cor", D = 1, W = 0, m0 = 0, C0 = 1) {
  k <- length(var_names)
  arg_list <- c(rep(0, k), list(order = k, name = name, D = D, W = W, m0 = m0, C0 = C0))
  names(arg_list) <- c(var_names, "order", "name", "D", "W", "m0", "C0")
  block <- do.call(polynomial_block, arg_list)

  block$FF <- simplify2array(apply(block$FF, 3, function(x) {
    diag(x) <- NA
    rbind(x, NA)
  },
  simplify = FALSE
  ))

  coef <- 0
  W <- 1
  block$G <- diag(c(rep(1, k), coef))
  block$D <- simplify2array(apply(block$D, 3, function(x) {
    as.matrix(bdiag(x, 1))
  },
  simplify = FALSE
  ))
  block$W <- simplify2array(apply(block$W, 3, function(x) {
    as.matrix(bdiag(x, W))
  },
  simplify = FALSE
  ))
  block$m0 <- c(block$m0, 1 / (1 - coef))
  block$C0 <- as.matrix(bdiag(block$C0, W))
  block$k <- k

  n <- k + 1
  block$names[[name]] <- c(block$names[[name]], n)
  block$order <- NULL
  block$n <- n


  block$type <- "Correlation"
  return(block)
}

#' block_merge
#'
#' An auxiliary function to merge blocks.
#'
#' @param ... dlm_block: A sequence of block to be merged.
#'
#' @return The merged block as a dlm_block.
#' @export
#'
#' @examples
#' level_1 <- polynomial_block(alpha1 = 1, order = 1)
#' level_2 <- polynomial_block(alpha2 = 1, order = 2)
#' season_2 <- harmonic_block(alpha2 = 1, period = 20)
#'
#' final_block <- block_merge(level_1, level_2, season_2)
block_merge <- function(...) {
  blocks <- list(...)
  for (block in blocks) {
    if (class(block) != "dlm_block") {
      stop(paste0("Error: Expected all arguments to be dlm_block's, but got a ", class(block), "."))
    }
  }

  n <- 0
  t <- 1
  k <- 1

  names <- list()
  var_names <- c()
  for (block in blocks) {
    ref_names <- block$names
    for (name in names(ref_names)) {
      ref_names[[name]] <- ref_names[[name]] + n
    }
    names <- c(names, ref_names)
    var_names <- c(var_names, block$var_names)
    if (block$t > 1) {
      if (block$t != t & t > 1) {
        stop(paste("Error: Blocks should have same length or length equal 1. Got", block$t, "and", t))
      }
      t <- block$t
    }
    n <- n + block$n
  }
  var_names <- sort(unique(var_names))
  k <- length(var_names)
  for (name in names(names)) {
    ref_idx <- which(names(names) == name)
    n_names <- length(ref_idx)
    if (n_names > 1) {
      names(names)[ref_idx] <- paste0(names(names)[ref_idx], "_", c(1:n_names))
    }
  }

  FF <- array(0, c(n, k, t), dimnames = list(NULL, var_names, NULL))
  G <- matrix(0, n, n)
  D <- array(0, c(n, n, t))
  W <- array(0, c(n, n, t))
  m0 <- c()
  C0 <- matrix(0, n, n)
  position <- 1
  for (block in blocks) {
    k_i <- length(block$var_names)
    current_range <- position:(position + block$n - 1)
    FF[current_range, var_names %in% block$var_names, ] <- block$FF[, (1:k_i)[order(block$var_names)], ]
    G[current_range, current_range] <- block$G
    D[current_range, current_range, ] <- block$D
    W[current_range, current_range, ] <- block$W
    m0 <- c(m0, block$m0)
    C0[current_range, current_range] <- block$C0
    position <- position + block$n
  }
  block <- list("FF" = FF, "G" = G, "D" = D, "W" = W, "m0" = m0, "C0" = C0, "n" = n, "t" = t, "k" = k, "names" = names, "var_names" = var_names)
  class(block) <- "dlm_block"
  return(block)
}
