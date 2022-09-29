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

#' normal_gamma_fit
#'
#' Fit the  normal model given the observed value and the model parameters.
#'
#' @param outcome Matrix: The observed data. It's dimension shoulb be T x m, where T is the length of the time series and m is the number of outcomes at each time.
#' @param m0 Vector: The prior mean for the latent vector.
#' @param C0 Matrix: The prior covariance matrix for the latent vector.
#' @param FF Array: A 3D-array containing the regression matrix for each time. It's dimension should be n x m x T, where n is the number of latent variables, m is the number of outcomes in the model and T is the time series length.
#' @param G Matrix: The state evolution matrix.
#' @param D Array: A 3D-array containing the discount factor matrix for each time. It's dimension should be n x n x T, where n is the number of latent variables and T is the time series length.
#' @param W Array: A 3D-array containing the covariance matrix of the noise for each time. It's dimension should be the same as D.
#' @param offset Matrix: The offset of the model. It's dimension should be the same as outcome.
#' @param parms list: a list contaning extra arguments. In this model, extra parameters are not used.
#' @param family function: the method for filtering.
#'
#' @return A list containing the following values:
#' \itemize{
#'    \item mt Matrix: The filtered mean of the latent variables for each time. Dimensions are n x T.
#'    \item Ct Array: A 3D-array containing the filtered covariance matrix of the latent variable for each time. Dimensions are n x n x T.
#'    \item ft Matrix: The one-step-ahead linear predictor for each time. Dimensions are m x T.
#'    \item qt Array: A 3D-array containing the one-step-ahead covariance matrix for the linear predictor for each time. Dimensions are m x T.
#'    \item FF Array: The same as the argument (same values).
#'    \item G Matrix: The same as the argument (same values).
#'    \item D Array: The same as the argument (same values).
#'    \item W Array: The same as the argument (same values).
#'    \item offset Vector: The same as the argument (same values).
#'    \item outcome Matrix: The same as the argument outcome (same values).
#'    \item parms: The same as the argument.
#' }
#' @export
#'
#' @examples
analytic_filter <- function(outcome, m0 = 0, C0 = 1, FF, G, D, W, offset, family, parms = list()) {

  # Definindo quantidades
  T <- nrow(outcome)
  r <- ncol(outcome)
  n <- dim(FF)[1]
  k <- dim(FF)[2]

  D <- ifelse(D == 0, 1, D)
  m0 <- matrix(m0, n, 1)
  C0 <- C0
  mt <- matrix(NA, nrow = n, ncol = T)
  Ct <- array(NA, dim = c(n, n, T))
  Rt <- array(NA, dim = c(n, n, T))
  ft <- matrix(NA, nrow = k, ncol = T)
  at <- matrix(NA, nrow = n, ncol = T)
  Qt <- array(NA, dim = c(k, k, T))

  param_names=family$param_names(outcome)
  conj_prior_param=matrix(NA,T,length(param_names)) %>% as.data.frame
  names(conj_prior_param)=param_names
  conj_post_param=conj_prior_param

  last_m <- m0
  last_C <- C0

  for (t in 1:T) {
    FF_step=matrix(FF[, , t], n, k)
    offset_step=offset[t,]
    na.flag=any(is.null(offset_step) | any(offset_step==0) | any(is.na(outcome[t,])))
    next_step=one_step_evolve(last_m,last_C,FF_step,G,D[, , t],W[, , t])

    at_step <- next_step$at
    Rt_step <- next_step$Rt

    next_step=family$offset(next_step$ft,next_step$Qt,if(na.flag){1}else{offset_step})

    ft_step <- next_step$ft
    Qt_step <- next_step$Qt

    conj_prior=family$conj_prior(ft_step,Qt_step,parms)
    conj_prior_param[t,]=conj_prior
    if(na.flag){
      conj_post_param[t,]=conj_prior
      norm_post=next_step
    }else{
      conj_post=family$update(conj_prior,outcome[t,],parms=parms)

      conj_post_param[t,]=conj_post
      norm_post=family$conj_post(conj_post,parms)
      }
    ft_star=norm_post$ft
    Qt_star=norm_post$Qt

    At <- Rt_step %*% FF_step %*% ginv(Qt_step)
    mt_step <- at_step + At %*% (ft_star - ft_step)
    Ct_step <- Rt_step + At %*% (Qt_star - Qt_step) %*% t(At)

    at[, t] <- at_step
    Rt[, , t] <- Rt_step

    ft[, t] <- ft_step
    Qt[, , t] <- Qt_step

    mt[, t] <- last_m <- mt_step
    Ct[, , t] <- last_C <- Ct_step
  }

  result <- list(
    "mt" = mt, "Ct" = Ct,
    "at" = at, "Rt" = Rt,
    "ft" = ft, "Qt" = Qt,
    "conj_prior_param"=conj_prior_param,"conj_post_param"=conj_post_param,
    "FF" = FF, "G" = G, "D" = D, "W" = W,
    "outcome" = outcome,'offset'=offset, "parms" = parms
  )
  return(result)
}

one_step_evolve=function(m0,C0,FF,G,D,W){
  at <- (G %*% m0)
  Pt <- G %*% C0 %*% (t(G))
  Rt <- as.matrix(D * Pt) + W

  ft <- (t(FF) %*% at)
  Qt <- as.matrix(t(FF) %*% Rt %*% FF)
  list('at'=at,'Rt'=Rt,'ft'=ft,'Qt'=Qt)
}

ident_offset=function(ft,Qt,offset){
  list('ft'=ft*offset,'Qt'=Qt*offset**2)
}

log_offset=function(ft,Qt,offset){
  list('ft'=ft+log(offset),'Qt'=Qt)
}

logit_offset=function(ft,Qt,offset){
  r=length(ft)
  list('ft'=ft+log(offset[1:r]/offset[r+1]),'Qt'=Qt)
}
