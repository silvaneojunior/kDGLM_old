#' polynomial_block
#'
#' Creates the structure for a polynomial block with desired order.
#'
#' @param order Positive integer: The order of the polimial structure.
#' @param values Matrix, vector or scalar: The values to be used on the first row of the regression matrix. If values is a matrix, it's dimensions should be k x t, where k is the number of outcomes of the model and t is the length of the outcome. If values is a vector and it's dimesions is equal to k (or k is Null), then it's values will be repeated for each times.  If values is a vector and it's dimesions is not equal to k (and k is not Null), then it's values will be repeated for each outcome (it's length will be used as time length). If values is a scalar, it's value will be used for all series and all times.
#' @param name String: An optional argument providing the name for this block. Can be useful to identify the models with meaningful labels, also, the name used will be used in some auxiliary functions.
#' @param D Array, Matrix, vector or  scalar: The values for the discount factors at each time. If D is a array, it's dimensions should be nxnxt, where n is the order of the polynomial block and t is the length of the outcomes. If D is a matrix, it's dimesions should be nxn and it's values will be used for each time. If D is a vector or scalar, a discount factor matrix will be created as a diagonal matrix with the values of D in the diagonal.
#' @param W Array, Matrix, vector or  scalar: The values for the covariance matrix for the noise factor at each time. If W is a array, it's dimensions should be nxnxt, where n is the order of the polynomial block and t is the length of the outcomes. If W is a matrix, it's dimesions should be nxn and it's values will be used for each time. If W is a vector or scalar, a discount factor matrix will be created as a diagonal matrix with the values of W in the diagonal.
#' @param m0 Vector or scalar: The prior mean for the latent variables associated with this block. If m0 is a vector, it's dimesion should be equal to the order of the polynomial block. If m0 is a scalar, it's value will be used for all latent variables.
#' @param C0 Matrix, vector or scalar: The prior covariance matrix for the latent variables associated with this block. If C0 is a matrix, it's dimesions should be nxn. If W is a vector or scalar, a covariance matrix will be created as a diagonal matrix with the values of C0 in the diagonal.
#' @param by_time Bool: A flag indicating if values should be interpreted as time or as outcomes when passing a vector.
#' @param k Positive integer: The number of outcomes in the model. Must be consistent with the dimension of values.
#'
#' @return A list containing the following values:
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
#' }
#'
#' @export
#' @examples
#' # Creating a first order structure for a model with 2 outcomes.
#' # One block is created for each outcome, with each block being associated with only one of the outcomes.
#' level_1 <- polynomial_block(order = 1, values = c(1, 0), by_time = FALSE)
#' level_2 <- polynomial_block(order = 1, values = c(1, 0), by_time = FALSE)
#'
#' # Creating a block with shared effect between the oucomes
#' level_3 <- polynomial_block(order = 2, values = c(1, 1), by_time = FALSE)
polynomial_block <- function(order, values = 1, name = "Var_Poly", D = 1, W = 0, m0 = 0, C0 = 1 / sqrt(mean(values**2, na.rm = TRUE)), by_time = TRUE, k = NULL) {
  multiple_block <- (k %>% if.null(1)) > 1 & by_time & (dim(values) %>% is.null())
  if (!is.null(k) & (!is.null(dim(values)) | !by_time)) {
    warning("The number of outcomes is being induced by values argument (values is a matrix or by_time is FALSE), but k has been passed. Ignoring k.")
  }

  if (!is.null(dim(values))) {
    k <- dim(values)[1]
    t <- dim(values)[2]
  } else {
    if (by_time) {
      t <- length(values)
    } else {
      k <- length(values)
      t <- 1
    }
  }
  k <- k %>% if.null(1)

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

  FF <- array(0, c(order, k, t))
  FF[1, , ] <- matrix(values, k, t, byrow = by_time)
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
    "type" = "Polynomial"
  )

  if (multiple_block) {
    block <- multiple_block(block, k, values = FF[1, 1, ])
  }


  return(block)
}


#' harmonic_block
#'
#' Creates the structure for a harmonic block with desired periodicity.
#'
#' @param period Positive integer: The size of the harmonic cycle.
#' @param values Matrix, vector or scalar: The values to be used on the first row of the regression matrix. If values is a matrix, it's dimensions should be k x t, where k is the number of outcomes of the model and t is the length of the outcome. If values is a vector and it's dimesions is equal to k (or k is Null), then it's values will be repeated for each times.  If values is a vector and it's dimesions is not equal to k (and k is not Null), then it's values will be repeated for each outcome (it's length will be used as time length). If values is a scalar, it's value will be used for all series and all times.
#' @param name String: An optional argument providing the name for this block. Can be useful to identify the models with meaningful labels, also, the name used will be used in some auxiliary functions.
#' @param D Array, Matrix, vector or  scalar: The values for the discount factors at each time. If D is a array, it's dimensions should be nxnxt, where n is the order of the polynomial block and t is the length of the outcomes. If D is a matrix, it's dimesions should be nxn and it's values will be used for each time. If D is a vector or scalar, a discount factor matrix will be created as a diagonal matrix with the values of D in the diagonal.
#' @param W Array, Matrix, vector or  scalar: The values for the covariance matrix for the noise factor at each time. If W is a array, it's dimensions should be nxnxt, where n is the order of the polynomial block and t is the length of the outcomes. If W is a matrix, it's dimesions should be nxn and it's values will be used for each time. If W is a vector or scalar, a discount factor matrix will be created as a diagonal matrix with the values of W in the diagonal.
#' @param m0 Vector or scalar: The prior mean for the latent variables associated with this block. If m0 is a vector, it's dimesion should be equal to the order of the polynomial block. If m0 is a scalar, it's value will be used for all latent variables.
#' @param C0 Matrix, vector or scalar: The prior covariance matrix for the latent variables associated with this block. If C0 is a matrix, it's dimesions should be nxn. If W is a vector or scalar, a covariance matrix will be created as a diagonal matrix with the values of C0 in the diagonal.
#' @param by_time Bool: A flag indicating if values should be interpreted as time or as outcomes when passing a vector.
#' @param k Positive integer: The number of outcomes in the model. Must be consistent with the dimension of values.
#'
#' @return A list containing the following values:
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
#' }
#'
#' @export
#' @examples
#' # Creating seasonal structure for a model with 2 outcomes.
#' # One block is created for each outcome, with each block being associated with only one of the outcomes.
#' season_1 <- harmonic_block(period = 3, values = c(1, 0), by_time = FALSE)
#' season_2 <- harmonic_block(order = 6, values = c(1, 0), by_time = FALSE)
#'
#' # Creating a block with shared effect between the oucomes
#' season_3 <- harmonic_block(order = 12, values = c(1, 1), by_time = FALSE)
harmonic_block <- function(period, values = 1, name = "Var_Sazo", D = 1, W = 0, m0 = 0, C0 = 1, by_time = TRUE, k = NULL) {
  w <- 2 * pi / period
  order <- 2
  block <- polynomial_block(order, values, name, D, W, m0, C0, by_time, k = NULL)

  G <- matrix(c(cos(w), -sin(w), sin(w), cos(w)), order, order)
  block$G <- G
  block$order <- NULL
  block$period <- period
  if (!is.null(k) & length(dim(values)) == 0) {
    block <- multiple_block(block, k, values = values)
  }
  block$type <- "Harmonic"
  return(block)
}

#' AR_block
#'
#' DESCRIPTION
#'
#' @param order Positive integer: DESCRIPTION.
#' @param values Matrix, vector or scalar: The values to be used on the first row of the regression matrix. If values is a matrix, it's dimensions should be k x t, where k is the number of outcomes of the model and t is the length of the outcome. If values is a vector and it's dimesions is equal to k (or k is Null), then it's values will be repeated for each times.  If values is a vector and it's dimesions is not equal to k (and k is not Null), then it's values will be repeated for each outcome (it's length will be used as time length). If values is a scalar, it's value will be used for all series and all times.
#' @param name String: An optional argument providing the name for this block. Can be useful to identify the models with meaningful labels, also, the name used will be used in some auxiliary functions.
#' @param D Array, Matrix, vector or  scalar: The values for the discount factors at each time. If D is a array, it's dimensions should be nxnxt, where n is the order of the polynomial block and t is the length of the outcomes. If D is a matrix, it's dimesions should be nxn and it's values will be used for each time. If D is a vector or scalar, a discount factor matrix will be created as a diagonal matrix with the values of D in the diagonal.
#' @param W Array, Matrix, vector or  scalar: The values for the covariance matrix for the noise factor at each time. If W is a array, it's dimensions should be nxnxt, where n is the order of the polynomial block and t is the length of the outcomes. If W is a matrix, it's dimesions should be nxn and it's values will be used for each time. If W is a vector or scalar, a discount factor matrix will be created as a diagonal matrix with the values of W in the diagonal.
#' @param m0 Vector or scalar: The prior mean for the latent variables associated with this block. If m0 is a vector, it's dimesion should be equal to the order of the polynomial block. If m0 is a scalar, it's value will be used for all latent variables.
#' @param C0 Matrix, vector or scalar: The prior covariance matrix for the latent variables associated with this block. If C0 is a matrix, it's dimesions should be nxn. If W is a vector or scalar, a covariance matrix will be created as a diagonal matrix with the values of C0 in the diagonal.
#' @param by_time Bool: A flag indicating if values should be interpreted as time or as outcomes when passing a vector.
#' @param k Positive integer: The number of outcomes in the model. Must be consistent with the dimension of values.
#'
#' @return A list containing the following values:
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
#' }
#'
#' @export
#' @examples
#' # EXAMPLE
AR_block <- function(order, values = 1, name = "Var_AR", D = 1, W = 0, m0 = 0, C0 = 1, by_time = TRUE, k = NULL) {
  block <- polynomial_block(2*order, values, name, D, W, m0, C0, by_time)

  if(order==1){
    G <- diag(2)
    G[1,1]=NA
  }else{
    G <- matrix(0, 2*order, 2*order)
    G[1,1:order]=NA
    G[2:order,-(order:(2*order))]=diag(order-1)
    G[(order+1):(2*order),(order+1):(2*order)]=diag(order)
    index=sort(c(c(1:order),c(1:order)))
    G=G[index+c(0,order),index+c(0,order)]
  }
  block$G <- G
  block$order <- order
  if (!is.null(k) & length(dim(values)) == 0) {
    block <- multiple_block(block, k, values = values)
  }
  block$type <- "AR"
  return(block)
}

#' block_merge
#'
#' An auxiliary function to merge blocks.
#'
#' @param ... <undefined class object> or list: A sequence of block to be merged.
#'
#' @return The merged block as a <undefined class object>.
#' @export
#'
#' @examples
#' level_1 <- polynomial_block(order = 1, values = matrix(c(rep(1, T), rep(0, T)), 2, T, byrow = TRUE))
#' level_2 <- polynomial_block(order = 2, values = matrix(c(rep(0, T), rep(1, T)), 2, T, byrow = TRUE))
#' season_2 <- harmonic_block(period = 20, values = matrix(c(rep(0, T), rep(1, T)), 2, T, byrow = TRUE))
#'
#' final_block <- block_merge(level_1, level_2, season_2)
block_merge <- function(...) {
  blocks <- list(...)

  n <- 0
  t <- 1
  k <- 1

  n_AR <- 0
  k_AR <- 0
  names <- list()
  for (block in blocks) {
    ref_names <- block$names
    # ref_AR_names <- block$AR_names
    for (name in names(ref_names)) {
      ref_names[[name]] <- ref_names[[name]] + n
    }
    # for (AR_name in names(ref_AR_names)) {
    #   ref_AR_names[[AR_name]] <- ref_AR_names[[AR_name]] + n_AR
    # }
    names <- c(names, ref_names)
    # AR_names <- c(AR_names, ref_AR_names)
    if (block$t > 1) {
      if (block$t != t & t > 1) {
        stop(paste("Error: Blocks should have same length or length equal 1. Got", block$t, "and", t))
      }
      t <- block$t
    }
    n <- n + block$n
    k <- max(block$k, k)
    n_AR <- n + if.null(block$AR_n, 0)
    k_AR <- k + if.null(block$AR_k, 0)
  }
  for (name in names(names)) {
    ref_idx <- which(names(names) == name)
    n_names <- length(ref_idx)
    if (n_names > 1) {
      names(names)[ref_idx] <- paste0(names(names)[ref_idx], "_", c(1:n_names))
    }
  }

  FF <- array(0, c(n, k, t))
  G <- matrix(0, n, n)
  D <- array(0, c(n, n, t))
  W <- array(0, c(n, n, t))
  m0 <- c()
  C0 <- matrix(0, n, n)
  position <- 1
  for (block in blocks) {
    current_range <- position:(position + block$n - 1)
    FF[current_range, , ] <- block$FF
    G[current_range, current_range] <- block$G
    D[current_range, current_range, ] <- block$D
    W[current_range, current_range, ] <- block$W
    m0 <- c(m0, block$m0)
    C0[current_range, current_range] <- block$C0
    position <- position + block$n
  }
  return(list("FF" = FF, "G" = G, "D" = D, "W" = W, "m0" = m0, "C0" = C0, "n" = n, "t" = t, "k" = k, "names" = names))
}

#' multiple_block
#'
#' An auxiliary function to repeat the same block for multiple outcomes. Useful if a regressor will be used in all outcomes, but with different effects in each one.
#'
#' @param ref_block <undefined class object> or list: The block to be copied.
#' @param k positive integer: The number of outcomes.
#' @param values vector: The values of the regressor at each time. If not provided, the first row of the FF array from the ref_block will be used.
#'
#' @return The merged block as a <undefined class object>.
#' @export
#'
#' @examples
#' level_i <- polynomial_block(order = 1, values = 0)
#' level <- multiple_block(level_i)
multiple_block <- function(ref_block, k, values = ref_block$FF[1, 1, ]) {
  if (ref_block$t > 1 & ref_block$t != length(values)) {
    stop("ERROR: ref_block have time length greater than 1, but not equal to the length of values")
  }
  n <- ref_block$n
  t <- length(values)
  if (ref_block$t == 1 | ref_block$k == 1) {
    ref_block$FF <- array(ref_block$FF, c(n, k, t))
  }
  aux_func <- function(i) {
    ref_block_i <- ref_block
    ref_block_i$FF[1, -i, ] <- 0
    ref_block_i$FF[1, i, ] <- values[]
    ref_block_i$k <- k
    ref_block_i$t <- t
    ref_block_i$name <- paste(ref_block_i$name, i, sep = "_")
    return(ref_block_i)
  }
  block_list <- lapply(1:k, aux_func)

  return(do.call(block_merge, block_list))
}
