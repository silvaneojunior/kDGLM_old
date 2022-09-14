#' generic_smoother
#'
#' Generic smoother for all models.
#'
#' @param mt Matrix: A matrix containing the filtered mean of the latent variables at each time. Each row should represent one variable.
#' @param Ct Array: A 3D-array representing the filtered covariance matrix of the latent variables at each time. The third dimension should represent the time index.
#' @param at Matrix: A matrix containing the one-step-ahead mean of the latent variables at each time based upon the filtered mean. Each row should represent one variable.
#' @param Rt Array: A 3D-array representing the one-step-ahead covariance matrix of the latent variables at each time based upon the filtered covariance matrix. The third dimension should represent the time index.
#' @param G  Matrix: A matrix representing the transition matrix of the model.
#'
#' @return List: The smoothed mean (mts) and covariance (Cts) of the latent variables at each time. Their dimension follows, respectivelly, the dimensions of mt and Ct.
#' @export
#' @importFrom MASS ginv
#'
#' @examples
#' T <- 20
#'
#' mt <- matrix(c(cumsum(rnorm(T) + 1), rep(1, T)), 2, T, byrow = TRUE)
#' Ct <- array(diag(c(1, 1)), c(2, 2, T))
#' G <- matrix(c(1, 0, 1, 1), 2, 2)
#' at <- G %*% mt
#' Rt <- array(G %*% t(G) + diag(c(0.1, 0.1)), c(2, 2, T))
#'
#' smoothed_values <- generic_smoother(mt, Ct, at, Rt, G)
generic_smoother <- function(mt, Ct, at, Rt, G) {
  T <- dim(mt)[2]
  n <- dim(mt)[1]
  mts <- mt
  Cts <- Ct

  var_index <- matrix(apply(Ct, 3, diag), n, T) != 0

  for (t in (T - 1):1) {
    var_ref <- var_index[, t]
    restricted_Rt <- Rt[var_ref, var_ref, t + 1]
    restricted_Ct <- Ct[var_ref, var_ref, t]

    simple_Rt_inv <- restricted_Ct %*% t(G[var_ref, var_ref]) %*% ginv(restricted_Rt)

    mts[var_ref, t] <- mt[var_ref, t] + simple_Rt_inv %*% (mts[var_ref, t + 1] - at[var_ref, t + 1])
    Cts[var_ref, var_ref, t] <- restricted_Ct - simple_Rt_inv %*% (restricted_Rt - Cts[var_ref, var_ref, t + 1]) %*% t(simple_Rt_inv)
  }
  return(list("mts" = mts, "Cts" = Cts))
}
