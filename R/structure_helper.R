#' Structural blocks for polynomial trends and regressions
#'
#' Creates the structure for a polynomial block with desired order.
#'
#' @param ... Named values for the planing matrix.
#' @param order Positive integer: The order of the polynomial structure.
#' @param name String: An optional argument providing the name for this block. Can be useful to identify the models with meaningful labels, also, the name used will be used in some auxiliary functions.
#' @param D Array, Matrix, vector or  scalar: The values for the discount factors at each time. If D is a array, it's dimensions should be n x n x T, where n is the order of the polynomial block and T is the length of the serie If D is a matrix, it's dimensions should be n x n and it's values will be used for each time. If D is a vector or scalar, a discount factor matrix will be created as a diagonal matrix with the values of D in the diagonal.
#' @param W Array, Matrix, vector or  scalar: The values for the covariance matrix for the noise factor at each time. If W is a array, it's dimensions should be n x n x T, where n is the order of the polynomial block and T is the length of the serie. If W is a matrix, it's dimensions should be n x n and it's values will be used for each time. If W is a vector or scalar, a discount factor matrix will be created as a diagonal matrix with the values of W in the diagonal.
#' @param m0 Vector or scalar: The prior mean for the latent variables associated with this block. If m0 is a vector, it's dimension should be equal to the order of the polynomial block. If m0 is a scalar, it's value will be used for all latent variables.
#' @param C0 Matrix, vector or scalar: The prior covariance matrix for the latent variables associated with this block. If C0 is a matrix, it's dimensions should be n x n. If W is a vector or scalar, a covariance matrix will be created as a diagonal matrix with the values of C0 in the diagonal. If the first element of C0 is equal to NA, then a proper variance is calculated for the first latent state based on the scale of the effect of that state on the linear predictor.
#' @param drift Matrix, vector or scalar: A drift to be add after the temporal evoltion (can be interpreted as the mean of the random noise at each time). If a matrix, it's dimension should be n x T, where n is the number of latent variables (i.e., the order) and T is the length of the serie. If a vector, it should have size T, and each value will be applied to the first latent variable (the one which affects the linear predictors) in their respective time. If a scalar, the passed value will be used for the first latent variable at each time.
#'
#' @return An object of the class dlm_block containing the following values:
#' \itemize{
#'    \item FF Array: A 3D-array containing the regression matrix for each time. It's dimension should be n x m x T, where n is the number of latent variables, m is the number of outcomes in the model and T is the time series length.
#'    \item FF_labs Array: DESCRIPTION
#'    \item G Matrix: The state evolution matrix.
#'    \item D Array: A 3D-array containing the discount factor matrix for each time. It's dimension should be n x n x T, where n is the number of latent variables and T is the time series length.
#'    \item W Array: A 3D-array containing the covariance matrix of the noise for each time. It's dimension should be the same as D.
#'    \item m0 Vector: The prior mean for the latent vector.
#'    \item C0 Matrix: The prior covariance matrix for the latent vector.
#'    \item drift Matrix: The mean for the random noise of the temporal evolution. It's dimension should be n x T.
#'    \item var_names list: A list containing the variables indexes by their name.
#'    \item order Positive integer: Same as argument.
#'    \item n Positive integer: The number of latent variables associated with this block (same value as order).
#'    \item t Positive integer: The number of time steps associated with this block. If 1, the block is compatible with blocks of any time length, but if t is greater than 1, this block can only be used with blocks of the same time length.
#'    \item k Positive integer: The number of outcomes associated with this block. This block can only be used with blocks with the same outcome length.
#'    \item pred_names Vector: The name of the linear predictors associated with this block,
#'    \item type Character: The type of block (polynomial).
#' }
#'
#' @export
#' @examples
#' # Creating a first order structure for a model with 2 outcomes.
#' # One block is created for each outcome
#' # with each block being associated with only one of the outcomes.
#' level_1 <- polynomial_block(alpha1 = 1, order = 1)
#' level_2 <- polynomial_block(alpha2 = 1, order = 1)
#'
#' # Creating a block with shared effect between the outcomes
#' level_3 <- polynomial_block(alpha1 = 1, alpha2 = 1, order = 2)
#'
#' @details
#'
#' For the ..., D, W, m0 and C0 arguments, the user may set one or more of it's values as a string.
#' By doing so, the user will leave the block partially undefined and it can no longer be used in the \code{\link{fit_model}} function.
#' Instead, the user must use the \code{\link{search_model}} function to search the best hyper parameters among a defined range of possible values.
#' See the \code{\link{search_model}} function for details on it's usage.
#'
#' For the details about the implementation see \insertCite{ArtigoPacote;textual}{kDGLM}.
#'
#' For the details about polynomial trend in the context of DLM's, see \insertCite{WestHarr-DLM;textual}{kDGLM}, chapter 7.
#'
#' For the details about dynamic regression models in the context of DLM's, see \insertCite{WestHarr-DLM;textual}{kDGLM}, chapters 6 and 9.
#'
#' @seealso \code{\link{fit_model}}
#' @seealso \code{\link{search_model}}
#' @family {auxiliary functions for structural blocks}
#'
#'
#' @references
#'    \insertAllCited{}
polynomial_block <- function(..., order = 1, name = "Var_Poly", D = 1, W = 0, m0 = 0, C0 = diag(order), drift = 0) {
  if (any(D > 1 | D < 0) & is.numeric(D)) {
    stop("Error: The discount factor D must be a value between 0 and 1 (included).")
  }
  if (any(D == 0)) {
    warning("Some value of D are equal to 0. Those values will be treated as 1.")
  }
  values <- list(...)
  vars <- names(values)
  if (any(vars == "")) {
    stop("Error: One or more linear predictors are unnamed. Please, name all arguments.")
  }
  if (any(vars == "const")) {
    stop("Error: Cannot create a linear predictor named 'const' (reserved name). Choose another label.")
  }
  var_len <- sapply(values, length)
  k <- length(values)
  t <- max(var_len)
  if (any(var_len != t & var_len > 1)) {
    stop(paste0("Error: Outcomes have mismatching lengths. Expected 1 or ", t, ", got: ", paste(var_len, collapse = ", "), "."))
  }

  if (length(D) == 1) {
    D <- array(D, c(order, order, t))
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
    pre_W <- diag(order)
    diag(pre_W) <- W
    W <- array(pre_W, c(order, order, t))
  } else if (is.vector(W)) {
    W_vals <- W
    pre_W <- diag(order)
    W <- array(0, c(order, order, length(W_vals)))
    for (i in 1:length(W_vals)) {
      diag(pre_W) <- W_vals[i]
      W[, , i] <- pre_W
    }
  } else if (is.matrix(W)) {
    W <- array(W, c(dim(W)[1], dim(W)[2], t))
  }
  t <- if (t == 1) {
    dim(W)[3]
  } else {
    t
  }

  if (length(dim(drift)) < 2) {
    if (t == 1 & length(drift) > 1) {
      t <- length(drift)
    }
    placeholder <- drift
    drift <- matrix(0, order, t)
    drift[1, ] <- placeholder
  }

  if (any(dim(drift) != c(order, t))) {
    stop(paste0("Error: Invalid shape for drift. Expected ", order, "x", t, ". Got ", paste(dim(drift), collapse = "x"), "."))
  }

  D <- array(D, c(order, order, t))
  W <- array(W, c(order, order, t))

  if (length(dim(W)) > 3 | any(dim(W)[1:2] != order) | (dim(W)[3] != t & t > 1)) {
    stop(paste0("Error: Invalid shape for W. Expected ", order, "x", order, "x", t, ". Got ", paste(dim(W), collapse = "x"), "."))
  }

  FF <- array(0, c(order, k, t), dimnames = list(NULL, vars, NULL))
  FF_labs <- matrix("const", order, k, dimnames = list(NULL, vars))
  for (i in 1:k) {
    name_var <- vars[i]
    if (typeof(values[[name_var]]) == "character") {
      FF[1, i, ] <- NA
      FF_labs[1, i] <- values[[name_var]]
      if (values[[name_var]] == "const") {
        stop("Error: Predictor value is equal to 'const', but 'const' is a reserved name. Choose another label.")
      }
    } else {
      FF[1, i, ] <- values[[name_var]]
    }
  }

  not_observed_flag <- is.na(FF) & (array(FF_labs, c(order, k, t)) == "const")
  D[, , apply(not_observed_flag, 3, any)] <- 1
  W[, , apply(not_observed_flag, 3, any)] <- 0
  FF <- ifelse(not_observed_flag, 0, FF)

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
  if (length(C0) == 1 | is.vector(C0)) {
    pre_C0 <- diag(order)
    diag(pre_C0) <- C0
    C0 <- pre_C0
  } else {
    C0
  }


  if (length(dim(C0)) > 2) {
    stop(paste0("Error: C0 must be a matrix, but it has ", length(dim(C0)), " dimensions."))
  }
  if (any(dim(C0) != order)) {
    stop(paste0("Error: C0 must have dimensions ", order, "x", order, ". Got ", dim(C0)[1], "x", dim(C0)[2], "."))
  }

  var_names <- list()
  var_names[[name]] <- c(1:order)
  block <- list(
    "FF" = FF,
    "FF_labs" = FF_labs,
    "drift" = drift,
    "G" = array(G, c(order, order, t)),
    "G_labs" = matrix("const", order, order),
    "D" = D,
    "W" = W,
    "m0" = m0,
    "C0" = C0,
    "var_names" = var_names,
    "order" = order,
    "n" = order,
    "t" = t,
    "k" = k,
    "pred_names" = vars,
    "type" = "Polynomial"
  )
  class(block) <- "dlm_block"
  block$status <- check.block.status(block)
  return(block)
}


#' Structural blocks for seasonal trends and regressions
#'
#' Creates the structure for a harmonic block with desired periodicity.
#'
#' @param ... Named values for the planing matrix.
#' @param period Positive integer: The size of the harmonic cycle.
#' @param name String: An optional argument providing the name for this block. Can be useful to identify the models with meaningful labels, also, the name used will be used in some auxiliary functions.
#' @param D Array, Matrix, vector or  scalar: The values for the discount factors at each time. If D is a array, it's dimensions should be n x n x t, where n is the order of the polynomial block and t is the length of the outcomes. If D is a matrix, it's dimensions should be n x n and it's values will be used for each time. If D is a vector or scalar, a discount factor matrix will be created as a diagonal matrix with the values of D in the diagonal.
#' @param W Array, Matrix, vector or  scalar: The values for the covariance matrix for the noise factor at each time. If W is a array, it's dimensions should be n x n x t, where n is the order of the polynomial block and t is the length of the outcomes. If W is a matrix, it's dimensions should be n x n and it's values will be used for each time. If W is a vector or scalar, a discount factor matrix will be created as a diagonal matrix with the values of W in the diagonal.
#' @param m0 Vector or scalar: The prior mean for the latent variables associated with this block. If m0 is a vector, it's dimension should be equal to the order of the polynomial block. If m0 is a scalar, it's value will be used for all latent variables.
#' @param C0 Matrix, vector or scalar: The prior covariance matrix for the latent variables associated with this block. If C0 is a matrix, it's dimensions should be n x n. If W is a vector or scalar, a covariance matrix will be created as a diagonal matrix with the values of C0 in the diagonal.
#' @param drift Matrix, vector or scalar: A drift to be add after the temporal evoltion (can be interpreted as the mean of the random noise at each time). If a matrix, it's dimension should be 2 x T, where T is the length of the serie. If a vector, it should have size T, and each value will be applied to the first latent variable (the one which affects the linear predictors) in their respective time. If a scalar, the passed value will be used for the first latent variable at each time.
#'
#' @return An object of the class dlm_block containing the following values:
#' \itemize{
#'    \item FF Array: A 3D-array containing the regression matrix for each time. It's dimension should be n x m x T, where n is the number of latent variables, m is the number of outcomes in the model and T is the time series length.
#'    \item FF_labs Array: DESCRIPTION
#'    \item drift Matrix: The mean for the random noise of the temporal evolution. It's dimension should be n x T.
#'    \item G Matrix: The state evolution matrix.
#'    \item D Array: A 3D-array containing the discount factor matrix for each time. It's dimension should be n x n x T, where n is the number of latent variables and T is the time series length.
#'    \item W Array: A 3D-array containing the covariance matrix of the noise for each time. It's dimension should be the same as D.
#'    \item m0 Vector: The prior mean for the latent vector.
#'    \item C0 Matrix: The prior covariance matrix for the latent vector.
#'    \item var_names list: A list containing the variables indexes by their name.
#'    \item period Positive integer: Same as argument.
#'    \item n Positive integer: The number of latent variables associated with this block (2).
#'    \item t Positive integer: The number of time steps associated with this block. If 1, the block is compatible with blocks of any time length, but if t is greater than 1, this block can only be used with blocks of the same time length.
#'    \item k Positive integer: The number of outcomes associated with this block. This block can only be used with blocks with the same outcome length.
#'    \item pred_names Vector: The name of the linear predictors associated with this block,
#'    \item type Character: The type of block (Harmonic).
#' }
#'
#' @export
#' @examples
#' # Creating seasonal structure for a model with 2 outcomes.
#' # One block is created for each outcome
#' # with each block being associated with only one of the outcomes.
#' season_1 <- harmonic_block(alpha1 = 1, period = 3)
#' season_2 <- harmonic_block(alpha2 = 1, period = 6)
#'
#' # Creating a block with shared effect between the outcomes
#' season_3 <- harmonic_block(alpha = 1, alpha2 = 1, period = 12)
#'
#' @details
#'
#' For the ..., D, W, m0 and C0 arguments, the user may set one or more of it's values as a string.
#' By doing so, the user will leave the block partially undefined and it can no longer be used in the \code{\link{fit_model}} function.
#' Instead, the user must use the \code{\link{search_model}} function to search the best hyper parameters among a defined range of possible values.
#' See the \code{\link{search_model}} function for details on it's usage.
#'
#' For the details about the implementation see \insertCite{ArtigoPacote;textual}{kDGLM}.
#'
#' For the details about the modelling of seasonal trends using harmonics in the context of DLM's, see \insertCite{WestHarr-DLM;textual}{kDGLM}, chapter 8.
#'
#' For the details about dynamic regression models in the context of DLM's, see \insertCite{WestHarr-DLM;textual}{kDGLM}, chapters 6 and 9.
#'
#' @seealso \code{\link{fit_model}}
#' @seealso \code{\link{search_model}}
#' @family {auxiliary functions for structural blocks}
#'
#'
#' @references
#'    \insertAllCited{}
harmonic_block <- function(..., period, name = "Var_Sazo", D = 1, W = 0, m0 = 0, C0 = diag(2), drift = 0) {
  w <- 2 * pi / period
  order <- 2
  block <- polynomial_block(..., order = order, name = name, D = D, W = W, m0 = m0, C0 = C0, drift = drift)

  G <- matrix(c(cos(w), -sin(w), sin(w), cos(w)), order, order)
  block$G <- array(G, c(order, order, block$t))
  block$order <- NULL
  block$period <- period
  block$type <- "Harmonic"
  return(block)
}

#' Structural blocks for auto regressive trends and regressions
#'
#' Creates the structure for a Auto Regressive (AR) block (see West & Harrison (1997), chapter 9) with desired order.
#' As the package suppose that the structure of the model is linear, a linearization is applied to the evolution equation, as described in West & Harrison (1997), chapter 13.
#' This block also supports Transfer Functions, being necessary to specify the associated pulse when calling the AR_block function (see arg.).
#'
#' @param ... Named values for the planing matrix.
#' @param order Positive integer: The order of the AR block.
#' @param noise_var Non negative scalar: The variance of the white noise added to the latent state.
#' @param noise_var Non negative scalar: The variance of the white noise added to the latent state.
#' @param pulse Vector or scalar: An optional argument providing the values for the pulse for a Transfer Function. Default is 0.
#' @param name String: An optional argument providing the name for this block. Can be useful to identify the models with meaningful labels, also, the name used will be used in some auxiliary functions.
#' @param AR_support String: Either "constrained" or "free". If AR_support is "constrained", then the AR coefficients will be forced to be on the interval (-1,1), otherwise, the coefficients will be unrestricted. Beware that, under no restriction on the coefficients, there is no guarantee that the estimated coefficients will imply in a stationary process, furthermore, if the order of the AR block is greater than 1, then the restriction imposed when AR_support is equal to "constrained" does NOT guarantee that the process will be stationary (although it may help).
#' @param m0 Vector or scalar: The prior mean for the coefficients associated with this block. If m0 is a vector, it's dimension should be equal to the order of the AR block. If m0 is a scalar, it's value will be used for all coefficients. If the coefficients are restricted to the interval (-1,1), the m0 is interpreted as the mean for logit((rho+1)/2), where rho is the AR coefficient.
#' @param C0 Matrix, vector or scalar: The prior covariance matrix for the coefficients associated with this block. If C0 is a matrix, it's dimensions should be n x n, where n is the order of the AR block. If C0 is a vector or scalar, a covariance matrix will be created as a diagonal matrix with the values of C0 in the diagonal. If the coefficients are restricted to the interval (-1,1), the C0 is interpreted as the covariance matrix for logit((rho+1)/2), where rho is the AR coefficient.
#' @param D Array, Matrix, vector or  scalar: The values for the discount factors associated with the AR coefficients at each time. If D is a array, it's dimensions should be n x n x t, where n is the order of the AR block and t is the length of the outcomes. If D is a matrix, it's dimensions should be n x n and it's values will be used for each time. If D is a vector or scalar, a discount factor matrix will be created as a diagonal matrix with the values of D in the diagonal.
#' @param W Array, Matrix, vector or  scalar: The values for the covariance matrix for the noise factor associated with the AR coefficients at each time. If W is a array, it's dimensions should be n x n x t, where n is the order of the AR block and t is the length of the outcomes. If W is a matrix, it's dimensions should be n x n and it's values will be used for each time. If W is a vector or scalar, a discount factor matrix will be created as a diagonal matrix with the values of W in the diagonal.
#' @param m0 Vector or scalar: The prior mean for the coefficients associated with this block. If m0 is a vector, it's dimension should be equal to the order of the AR block. If m0 is a scalar, it's value will be used for all coefficients. If the coefficients are restricted to the interval (-1,1), the m0 is interpreted as the mean for logit((rho+1)/2), where rho is the AR coefficient.
#' @param C0 Matrix, vector or scalar: The prior covariance matrix for the coefficients associated with this block. If C0 is a matrix, it's dimensions should be n x n, where n is the order of the AR block. If C0 is a vector or scalar, a covariance matrix will be created as a diagonal matrix with the values of C0 in the diagonal. If the coefficients are restricted to the interval (-1,1), the C0 is interpreted as the covariance matrix for logit((rho+1)/2), where rho is the AR coefficient.
#' @param drift Matrix, vector or scalar: A drift to be add in the AR coefficients after the temporal evoltion (can be interpreted as the mean of the random noise at each time). If a matrix, it's dimension should be n x T, where n is the order of the AR block and T is the length of the serie. If a scalar, the passed value will be used for all coefficients at each time.
#' @param m0_states Vector or scalar: The prior mean for the states associated with this block. If m0_states is a vector, it's dimension should be equal to the order of the AR block. If m0_states is a scalar, it's value will be used for all coefficients.
#' @param C0_states Matrix, vector or scalar: The prior covariance matrix for the states associated with this block. If C0_states is a matrix, it's dimensions should be n x n, where n is the order of the AR block. If C0_state is a vector or scalar, a covariance matrix will be created as a diagonal matrix with the values of C0_state in the diagonal.
#' @param D_states Array, Matrix, vector or  scalar: The values for the discount factors for the states associated with this block. If D_states is a array, it's dimensions should be n x n x t, where n is the order of the AR block and t is the length of the outcomes. If D_states is a matrix, it's dimensions should be n x n and it's values will be used for each time. If D_states is a vector or scalar, a discount factor matrix will be created as a diagonal matrix with the values of D_states in the diagonal.
#' @param drift_states Matrix, vector or scalar: A drift to be add in the states after the temporal evoltion (can be interpreted as the mean of the random noise at each time). If a matrix, it's dimension should be n x T, where n is the order of the AR block and T is the length of the serie. If a vector, it should have size T, and each value will be applied to the first latent variable (the one which affects the linear predictors) in their respective time. If a scalar, the passed value will be used for the first latent variable at each time.
#' @param m0_pulse Vector or scalar: The prior mean for the coefficients associated with the pulses. If m0_pulse is a vector, it's dimension should be equal to the number of pulses. If m0_pulse is a scalar, it's value will be used for all coefficients.
#' @param C0_pulse  Matrix, vector or scalar: The prior covariance matrix for the coefficients associated with the pulses. If C0_pulse is a matrix, it's dimensions should be n x n, where n is the number of pulses. If C0_pulse is a vector or scalar, a covariance matrix will be created as a diagonal matrix with the values of C0_pulse in the diagonal.
#' @param D_pulse Array, Matrix, vector or  scalar: The values for the discount factors associated with pulse coefficients at each time. If D_pulse is a array, it's dimensions should be n x n x t, where n is the number of pulses and t is the length of the outcomes. If D_pulse is a matrix, it's dimensions should be n x n and it's values will be used for each time. If D_pulse is a vector or scalar, a discount factor matrix will be created as a diagonal matrix with the values of D_pulse in the diagonal.
#' @param W_pulse Array, Matrix, vector or  scalar: The values for the covariance matrix for the noise factor associated with pulse coefficients at each time. If W_pulse is a array, it's dimensions should be n x n x t, where n is the number of pulses and t is the length of the outcomes. If W_pulse is a matrix, it's dimensions should be n x n and it's values will be used for each time. If W_pulse is a vector or scalar, a covariance matrix will be created as a diagonal matrix with the values of W_pulse in the diagonal.
#' @param drift_pulse Matrix, vector or scalar: A drift to be add in the pulse effect after the temporal evoltion (can be interpreted as the mean of the random noise at each time). If a matrix, it's dimension should be n x T, where n is the number of pulses and T is the length of the serie. If a scalar, the passed value will be used for all latent variable at each time.
#'
#' @return An object of the class dlm_block containing the following values:
#' \itemize{
#'    \item FF Array: A 3D-array containing the regression matrix for each time. It's dimension should be n x m x T, where n is the number of latent variables, m is the number of outcomes in the model and T is the time series length.
#'    \item FF_labs Array: DESCRIPTION
#'    \item G Matrix: The state evolution matrix.
#'    \item D Array: A 3D-array containing the discount factor matrix for each time. It's dimension should be n x n x T, where n is the number of latent variables and T is the time series length.
#'    \item W Array: A 3D-array containing the covariance matrix of the noise for each time. It's dimension should be the same as D.
#'    \item m0 Vector: The prior mean for the latent vector.
#'    \item C0 Matrix: The prior covariance matrix for the latent vector.
#'    \item var_names list: A list containing the variables indexes by their name.
#'    \item order Positive integer: Same as argument.
#'    \item n Positive integer: The number of latent variables associated with this block (2).
#'    \item t Positive integer: The number of time steps associated with this block. If 1, the block is compatible with blocks of any time length, but if t is greater than 1, this block can only be used with blocks of the same time length.
#'    \item k Positive integer: The number of outcomes associated with this block. This block can only be used with blocks with the same outcome length.
#'    \item pred_names Vector: The name of the linear predictors associated with this block,
#'    \item type Character: The type of block (AR).
#'    \item AR_support Character: Same as argument.
#' }
#'
#' @export
#' @examples
#'
#' AR_block(mu = 1, pulse = rnorm(200), order = 3, noise_var = 0.1)
#'
#' @details
#'
#' For the ..., noise_var, D, W, m0, C0, m0_states, C0_states, D_states, m0_pulse, C0_pulse, D_pulse, W_pulse arguments, the user may set one or more of it's values as a string.
#' By doing so, the user will leave the block partially undefined and it can no longer be used in the \code{\link{fit_model}} function.
#' Instead, the user must use the \code{\link{search_model}} function to search the best hyper parameters among a defined range of possible values.
#' See the \code{\link{search_model}} function for details on it's usage.
#'
#' For the details about the implementation see \insertCite{ArtigoPacote;textual}{kDGLM}.
#'
#' For the details about Auto regressive models in the context of DLM's, see \insertCite{WestHarr-DLM;textual}{kDGLM}, chapter 9.
#'
#' For the details about the linearization of non-linear evolution equations in the context of DLM's, see \insertCite{WestHarr-DLM;textual}{kDGLM}, chapter 13.
#'
#' For the details about dynamic regression models in the context of DLM's, see \insertCite{WestHarr-DLM;textual}{kDGLM}, chapters 6 and 9.
#'
#' @seealso \code{\link{fit_model}}
#' @family {auxiliary functions for structural blocks}
#'
#' @references
#'    \insertAllCited{}
AR_block <- function(..., order, noise_var, pulse = 0, name = "Var_AR", AR_support = "free",
                     D = 1, W = 0, m0 = c(1, rep(0, order - 1)), C0 = 1, drift = 0,
                     m0_states = 0, C0_states = diag(order), D_states = 1, drift_states = 0,
                     m0_pulse = 0, C0_pulse = 1, D_pulse = 1, W_pulse = 0, drift_pulse = 0) {
  W_states <- array(0, c(order, order, length(noise_var)))
  W_states[1, 1, ] <- noise_var
  block_state <-
    polynomial_block(..., order = order, name = paste0(name, "_State"), m0 = m0_states, C0 = C0_states, D = D_states, W = W_states, drift = drift_states)


  dummy_var <- list()
  dummy_var[[names(list(...))[1]]] <- rep(0, block_state$t)
  if (length(drift_states) > 1 & length(dim(drift_states)) < 2) {
    drift_states <- matrix(drift_states, order, length(drift_states))
  }
  block_coeff <-
    do.call(
      function(...) {
        polynomial_block(..., order = order, name = paste0(name, "_Coeff"), m0 = m0, C0 = C0, D = D, W = W, drift = drift_states)
      },
      dummy_var
    )

  block <- block_state + block_coeff
  k <- block$k

  if (order == 1) {
    G <- diag(2)
    G_labs <- matrix("const", 2, 2)
    G[1, 1] <- NA
    G_labs[1, 1] <- tolower(AR_support)
  } else {
    G <- matrix(0, 2 * order, 2 * order)
    G_labs <- matrix("const", 2 * order, 2 * order)
    G[1, 1:order] <- NA
    G_labs[1, 2 * (1:order) - 1] <- tolower(AR_support)
    G[2:order, -(order:(2 * order))] <- diag(order - 1)
    G[(order + 1):(2 * order), (order + 1):(2 * order)] <- diag(order)
    index <- sort(c(c(1:order), c(1:order)))
    G <- G[index + c(0, order), index + c(0, order)]
  }
  block$G <- array(G, c(2 * order, 2 * order, block$t))
  block$G_labs <- G_labs
  block$order <- order
  block$type <- "AR"
  block$AR_support <- AR_support
  true_order <- c(rbind(1:order, 1:order + order))
  block$m0 <- block$m0[true_order]
  block$C0 <- block$C0[true_order, true_order]
  block$D <- block$D[true_order, true_order, ]
  block$W <- block$W[true_order, true_order, ]
  block$var_names[[1]] <- seq(1, 2 * order, 2)[block$var_names[[1]]]
  block$var_names[[2]] <- seq(2, 2 * order, 2)[block$var_names[[2]] - order]
  if (any(pulse != 0)) {
    k <- if.null(dim(pulse)[2], 1)
    t <- if.null(dim(pulse)[1], length(pulse))
    dummy_var <- list()
    dummy_var[[names(list(...))[1]]] <- rep(0, t)
    if (length(drift_pulse) > 1 & length(dim(drift_pulse)) < 2) {
      drift_pulse <- matrix(drift_pulse, order, length(drift_pulse))
    }
    block_pulse <-
      do.call(
        function(...) {
          polynomial_block(..., order = k, name = paste0(name, "_Pulse"), m0 = m0_pulse, C0 = C0_pulse, D = D_pulse, W = W_pulse, drift = drift_pulse)
        },
        dummy_var
      )
    block_pulse$G <- diag(k)
    block <- block + block_pulse
    block$G[1, (2 * order + 1):(2 * order + k), ] <- t(matrix(pulse, block$t, k))
  }
  return(block)
}

#' Auxiliary function to merge blocks
#'
#' An auxiliary function to merge blocks.
#'
#' @param ... dlm_block: A sequence of block to be merged.
#'
#' @return The merged block as a dlm_block.
#' @export
#'
#' @examples
#'
#' # Long way
#' level_1 <- polynomial_block(alpha1 = 1, order = 1)
#' level_2 <- polynomial_block(alpha2 = 1, order = 2)
#' season_2 <- harmonic_block(alpha2 = 1, period = 20)
#'
#' final_block <- block_merge(level_1, level_2, season_2)
#'
#' # Short way
#' final_block <- polynomial_block(alpha1 = 1, order = 1) +
#'   polynomial_block(alpha2 = 1, order = 2) +
#'   harmonic_block(alpha2 = 1, period = 20)
#'
#' @family {auxiliary functions for structural blocks}
block_merge <- function(...) {
  blocks <- list(...)
  if (length(blocks) == 1) {
    return(blocks[[1]])
  }
  for (block in blocks) {
    if (!inherits(block, "dlm_block")) {
      stop(paste0("Error: Expected all arguments to be dlm_block's, but got a ", class(block), "."))
    }
  }

  n <- 0
  t <- 1
  k <- 1

  var_names <- list()
  pred_names <- c()
  for (block in blocks) {
    ref_names <- block$var_names
    for (name in names(ref_names)) {
      ref_names[[name]] <- ref_names[[name]] + n
      var_names[[name]] <- c(var_names[[name]], ref_names[[name]])
    }
    pred_names <- c(pred_names, block$pred_names)
    if (block$t > 1) {
      if (block$t != t & t > 1) {
        stop(paste("Error: Blocks should have same length or length equal 1. Got", block$t, "and", t))
      }
      t <- block$t
    }
    n <- n + block$n
  }
  pred_names <- sort(unique(pred_names))
  k <- length(pred_names)

  FF <- array(0, c(n, k, t), dimnames = list(NULL, pred_names, NULL))
  FF_labs <- matrix("const", n, k, dimnames = list(NULL, pred_names))
  drift <- matrix(0, n, t)
  G <- array(0, c(n, n, t))
  G_labs <- matrix("const", n, n)
  D <- array(0, c(n, n, t))
  W <- array(0, c(n, n, t))
  m0 <- c()
  C0 <- matrix(0, n, n)
  position <- 1
  status <- "defined"
  for (block in blocks) {
    k_i <- length(block$pred_names)
    current_range <- position:(position + block$n - 1)
    FF[current_range, pred_names %in% block$pred_names, ] <- block$FF[, (1:k_i)[order(block$pred_names)], ]
    FF_labs[current_range, pred_names %in% block$pred_names] <- block$FF_labs[, (1:k_i)[order(block$pred_names)]]

    drift[current_range, ] <- block$drift
    G[current_range, current_range, ] <- block$G
    G_labs[current_range, current_range] <- block$G_labs
    D[current_range, current_range, ] <- block$D
    W[current_range, current_range, ] <- block$W
    m0 <- c(m0, block$m0)
    C0[current_range, current_range] <- block$C0
    position <- position + block$n
  }
  block <- list(
    "FF" = FF,
    "FF_labs" = FF_labs,
    "drift" = drift,
    "G" = G,
    "G_labs" = G_labs,
    "D" = D,
    "W" = W,
    "m0" = m0,
    "C0" = C0,
    "n" = n,
    "t" = t,
    "k" = k,
    "status" = status,
    "var_names" = var_names,
    "pred_names" = pred_names,
    "type" = "Mixed"
  )
  class(block) <- "dlm_block"
  block$status <- check.block.status(block)
  return(block)
}

#' Auxiliary function to replicate blocks
#'
#' An auxiliary function to merge blocks.
#'
#' @param block dlm_block: A block to be multiplied.
#' @param k Integer: The number of blocks to generate.
#'
#' @return The merged block as a dlm_block.
#' @export
#'
#' @examples
#' # Long way
#' level <- polynomial_block(alpha = 1, order = 1)
#'
#' final_block <- block_mult(level, 5)
#'
#' # Short way
#' final_block <- 5 * polynomial_block(alpha = 1, order = 1)
#'
#' @seealso \code{\link{block_rename}}
#' @family {auxiliary functions for structural blocks}
block_mult <- function(block, k) {
  block_list <- list()
  size_total <- floor(log10(k)) + 1
  if (k > 1) {
    block_ref <- block
    block$pred_names <- paste0(block$pred_names, "_", paste0(rep("0", size_total - 1), collapse = ""), "1")
    block_list[[1]] <- block
    for (i in 2:k) {
      size_i <- floor(log10(i)) + 1
      block_clone <- block_ref
      block_clone$pred_names <- paste0(block_ref$pred_names, "_", paste0(rep("0", size_total - size_i), collapse = ""), i)
      block_list[[i]] <- block_clone
    }
    block <- do.call(block_merge, block_list)
  }
  return(block)
}

#' block_rename
#'
#' @param block A dlm_block object.
#' @param pred_names A vector of string with names for each linear predictor in block.
#'
#' @return A dlm_block with the linear predictors renamed to the values passed in names.
#' @export
#'
#' @examples
#'
#' base_block <- polynomial_block(
#'   eta = 1,
#'   order = 1,
#'   name = "Poly",
#'   D = 0.95
#' )
#'
#' final_block <- block_rename(2 * base_block, c("mu", "sigma"))
#'
#' @family {auxiliary functions for structural blocks}
block_rename <- function(block, pred_names) {
  if (!inherits(block, "dlm_block")) {
    stop("Error: The block argument is not a dlm_block object.")
  }
  if (length(pred_names) != length(block$pred_names)) {
    stop(paste0("Error: The number of names provided does not match the number of linear predictor in the block. Expected ", length(block$pred_names), ", got ", length(names), "."))
  }
  if (length(pred_names) != length(unique(pred_names))) {
    stop(paste0("Error: Repeated names are not allowed."))
  }

  block$pred_names <- pred_names
  colnames(block$FF) <- pred_names
  block$status <- check.block.status(block)
  return(block)
}
