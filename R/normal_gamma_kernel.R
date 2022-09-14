#' normal_gamma_filter
#'
#' Filtering function for the Normal model with unknown variance.
#'
#' @param outcome Vector: The observed values at time t. Must have length greater than 1.
#' @param m0 Vector: The prior mean for the latent vector at time t.
#' @param C0 Matrix: The prior covariance matrix for the latent vector at time t.
#' @param FF Matrix: The regression matrix at time t. Should be compatible with the dimensions of m0 and outcome, i.e., if m0 has dimension n and outcome has dimension m, FF should be n x m.
#' @param G Matrix: The state evolution matrix.
#' @param D Matrix: The discount factor matrix at time t.
#' @param W Matrix: The noise matrix at time t.
#' @param offset Vector: Same dimension as outcome. A vector contaning the offset at time t.
#' @param parms list: a list contaning extra arguments. In this model, extra parameters are not used.
#'
#' @return A list containing:
#' \itemize{
#'  \item at: One-step-ahead mean for the latent vectors.
#'  \item Rt: One-step-ahead covariance matrix for the latent vectors.
#'  \item ft: One-step-ahead linear predictor.
#'  \item Qt: One-step-ahead covariance matrix for the linear predictior.
#'  \item tau: BLANK (see Ref. Raíra).
#'  \item tau_star: BLANK (see Ref. Raíra).
#'  \item mt: The filtered mean of the latent vector at time t.
#'  \item Ct: The filtered covariance matrix of the latent vector at time t.
#'  \item outcome: The observed value at time t.
#'  \item params: The same as the argument.
#'  }
#' @export
#' @importFrom MASS ginv
#' @importFrom rootSolve multiroot
#'
#' @examples
#' outcome <- c(rnorm(1, 0, 1), 0)
#' m0 <- c(2, -log(1.1))
#' C0 <- diag(2)
#' G <- matrix(c(1, 0, 0, 1), 2, 2)
#' FF <- matrix(c(1, 0, 0, 1), 2, 2)
#' D <- diag(c(0.1, 0.1), 2, 2) + 1
#' W <- diag(c(0, 0), 2, 2)
#' offset <- c(0)
#' filtered_data <- GDLM::normal_gamma_filter(outcome, m0, C0, FF, G, D, W, offset)
normal_gamma_filter <- function(outcome, m0, C0, FF, G, D, W, offset = NULL, parms = list()) {
  offset=if(is.null(offset)){1}else{offset}
  r <- 2

  y=outcome[1]

  at <- (G %*% m0)
  Pt <- G %*% C0 %*% (t(G))
  Rt <- as.matrix(D * Pt) + W

  ft <- (t(FF) %*% at)*offset
  Qt <- as.matrix(t(FF) %*% Rt %*% FF)*(offset**2)

  # Normal prior
  f1 <- ft[1,]
  f2 <- ft[2,]
  q1 <- Qt[1, 1]
  q2 <- Qt[2, 2]
  q12 <- Qt[1, 2]

  ###########################################
  # IGNORE ESSE SISTEMA
  # Esse trecho serve para testar se o problema no ajuste
  # vem da aproximação da função digamma (o que não é o caso).
  #
  # system <- function(x, parms) {
  #   outcome=exp(x)
  #   out <- digamma(outcome)-x + parms$q2 / 2
  #   return(out)
  # }
  #
  # s <- multiroot(
  #   f = system,
  #   start = c(0),
  #   parms = list("q2" = q2)
  # )
  ###########################################

  # Conjugated prior
  mu0 <- f1+q12
  c0 <- 1/(exp(f2 + q2 / 2)*q1)
  helper=-3+3*sqrt(1+2*q2/3)
  alpha=1/helper
  beta <- alpha*exp(-f2 - q2 / 2)

  # Updating parameter
  mu0_star <- (c0 * mu0 + y) / (c0 + 1)
  c0_star <- c0 + 1
  alpha_star <- alpha + 0.5
  beta_star <- beta + 0.5 * c0 * ((mu0 - y)**2) / (c0 + 1)

  # Normal posterior
  f1star <- mu0_star
  f2star <- digamma(alpha_star)-log(beta_star)
  q1star <- beta_star / (c0_star * (alpha_star-1))
  q2star <- trigamma(alpha_star)
  q12star <- 0

  ####################################################
  # Se estas equações forem usadas no lugar das que estão acima o ajuste fica bom, mesmo em casos relativamente extremos (priori vaga).
  #
  # q1star <- beta_star / (c0_star * (alpha_star))
  # q2star <- 2*log(alpha_star)-2*digamma(alpha_star)
  ####################################################

  fstar <- c(f1star, f2star)
  Qstar <- matrix(c(q1star, q12star, q12star, q2star), byrow = F, ncol = 2)

  At <- Rt %*% FF %*% ginv(Qt)
  mt <- at + At %*% (fstar - ft)
  Ct <- Rt + At %*% (Qstar - Qt) %*% t(At)

  return(list(
    "at" = at, "Rt" = Rt,
    "ft" = ft, "Qt" = Qt,
    "mu0" = mu0, "c0" = c0, "alpha" = alpha, "beta" = beta,
    "mu0_star" = mu0_star, "c0_star" = c0_star, "alpha_star" = alpha_star, "beta_star" = beta_star,
    "tau0" = c0/2, "tau1" = c0*mu0, "tau2" = -(c0*(mu0**2)/2+beta), "tau3" = alpha-1,
    "tau0_star" = c0_star/2, "tau1_star" = c0_star*mu0_star, "tau2_star" = -(c0_star*(mu0_star**2)/2+beta_star), "tau3_star" = alpha_star-1,
    "mt" = mt, "Ct" = Ct,
    "outcome" = outcome, "parms" = parms
  ))
}


#' normal_gamma_pred
#'
#' Calculate the values for the predictive distribuition given the values of the parameter of the conjugated distribuition of the linear predictor.
#'
#' @param model List: contaning the parameters of the conjugated distribuition for the linear predictor (may be based on the prior, filtered or smoothed distribuition).
#' @param pred_cred Numeric: the desired credibility for the credibility interval.
#'
#' @return A list containing the following values:
#' \itemize{
#'    \item pred vector/matrix: the mean of the predictive distribuition of a next observation. Same type and shape as the parameter in model.
#'    \item var.pred vector/matrix: the variance of the predictive distribuition of a next observation. Same type and shape as the parameter in model.
#'    \item icl.pred vector/matrix: the percentile of 100*((1-pred_cred)/2)% of the predictive distribuition of a next observation. Same type and shape as the parameter in model.
#'    \item icu.pred vector/matrix: the percentile of 100*(1-(1-pred_cred)/2)% of the predictive distribuition of a next observation. Same type and shape as the parameter in model.
#' }
#' @export
#' @importFrom extraDistr qbetapr
#'
#' @examples
#' # A fitted model shoulb be used as argument, but you can also pass only the parameter themselves.
#'
#' model <- list(
#'   "tau0" = 0.5,
#'   "tau1" = -1 / 2,
#'   "tau2" = 0,
#'   "tau3" = -2.0,
#'   "Qt" = matrix(1, 1, 1),
#'   "outcome" = 0
#' )
#'
#' normal_gamma_pred(model)
normal_gamma_pred <- function(model, pred_cred = 0.95) {
  c0=model$c0
  mu0=model$mu0
  alpha=model$alpha
  beta=model$beta

  mu=mu0
  nu=2*alpha
  sigma2=(beta/alpha)*(1+1/c0)

  pred=mu %>% as.matrix %>% t
  var.pred=beta/((alpha-1))*(1+1/c0) %>% as.matrix %>% t

  icl.pred=qt((1-pred_cred)/2,nu)*sqrt(sigma2)+mu
  icu.pred=qt(1-(1-pred_cred)/2,nu)*sqrt(sigma2)+mu

  list(
    "pred"     = pred,
    "var.pred" = var.pred,
    "icl.pred" = icl.pred,
    "icu.pred" = icu.pred,
    "log.like" = 0
  )
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
#' @param pred_cred Numeric: the desired credibility for the credibility interval.
#' @param parms list: a list contaning extra arguments. In this model, extra parameters are not used.
#' @param kernel_filter function: the method for filtering.
#'
#' @return A list containing the following values:
#' \itemize{
#'    \item mt Matrix: The filtered mean of the latent variables for each time. Dimensions are n x T.
#'    \item Ct Array: A 3D-array containing the filtered covariance matrix of the latent variable for each time. Dimensions are n x n x T.
#'    \item ft Matrix: The one-step-ahead linear predictor for each time. Dimensions are m x T.
#'    \item qt Array: A 3D-array containing the one-step-ahead covariance matrix for the linear predictor for each time. Dimensions are m x T.
#'    \item a Matrix: The alpha parameter for the gamma prior for each time. Dimensions are m x T.
#'    \item b Matrix: The beta parameter for the gamma prior for each time. Dimensions are m x T.
#'    \item a.post Matrix: The alpha parameter for the gamma posteior for each time. Dimensions are m x T.
#'    \item b.post Matrix: The beta parameter for the gamma posteior for each time. Dimensions are m x T.
#'    \item FF Array: The same as the argument (same values).
#'    \item G Matrix: The same as the argument (same values).
#'    \item D Array: The same as the argument (same values).
#'    \item W Array: The same as the argument (same values).
#'    \item pred Matrix: The one-step-ahead predictions for each time. Dimensions are m x T.
#'    \item var.pred Matrix: The variance for the one-step-ahead predictions for each time. Dimensions are m x T. Note that, in the multivariate Poisson case, the series are supossed independent, so, in particular, they are uncorrelated.
#'    \item icl.pred Matrix: The lower credibility interval for the prediction at each time. Dimensions are m x T.
#'    \item icu.pred Matrix: The upper credibility interval for the prediction at each time. Dimensions are m x T.
#'    \item pred_cred Numeric: Deprecated
#'    \item offset Vector: The same as the argument (same values).
#'    \item log.like Matrix: log likelyhood for the filtered parameters based on the one-step ahead prediction.
#'    \item outcome Matrix: The same as the argument outcome (same values).
#'    \item parms: The same as the argument.
#' }
#' @export
#'
#' @examples
#' # Not ideal way: should use fit_model function.
#' T <- 200
#' mean <- 2 * 1:T / 200
#' var <- 3 * 1:T / 200
#' outcome <- matrix(rt(T, 10) * 2 * sqrt(var) + mean)
#' m0 <- c(0, 0, 0, 0)
#' C0 <- diag(c(1, 1, 1, 1))
#' G <- as.matrix(Matrix::bdiag(
#'   matrix(c(1, 0, 1, 1), 2, 2),
#'   matrix(c(1, 0, 1, 1), 2, 2)
#' ))
#' FF <- array(matrix(c(1, 0, 0, 0, 0, 0, 1, 0), 4, 2), c(4, 2, T))
#' D <- array(diag(c(0.1, 0, 0.1, 0)), c(4, 4, T)) + 1
#' W <- array(diag(c(0, 0, 0)), c(4, 4, T))
#' offset <- matrix(1, T, 1)
#'
#' fitted_data <- GDLM::normal_gamma_fit(outcome = outcome, m0 = m0, C0 = C0, FF = FF, G = G, D = D, W = W, offset = offset, pred_cred = 0.95)
#'
#' plot(outcome)
#' lines(fitted_data$pred[1, ])
normal_gamma_fit <- function(outcome, m0 = 0, C0 = 1, FF, G, D, W, offset, pred_cred = 0.95, parms = list(), kernel_filter = normal_gamma_filter) {

  # Definindo quantidades
  T <- nrow(outcome)
  n <- dim(FF)[1]

  D.aux <- D
  D <- ifelse(D.aux == 0, 1, D.aux)
  r <- 2
  m0 <- matrix(m0, n, 1)
  C0 <- C0
  mt <- matrix(0, nrow = n, ncol = T)
  Ct <- array(rep(diag(n), T), dim = c(n, n, T))
  Rt <- array(rep(0, T), dim = c(n, n, T))
  ft <- matrix(0, nrow = r, ncol = T)
  at <- matrix(0, nrow = n, ncol = T)
  Qt <- array(0, dim = c(r, r, T))
  # At <- array(0,dim=c(n,r,T))

  pred <- matrix(0, r - 1, T)
  var.pred <- array(0, c(r - 1, r - 1, T))
  icl.pred <- matrix(0, r - 1, T)
  icu.pred <- matrix(0, r - 1, T)
  log.like <- matrix(0, r - 1, T)

  # f_star <- matrix(0,nrow=T,ncol=r)
  # Q_star <- array(0,c(r,r,T))

  mt <- matrix(0, nrow = n, ncol = T)
  Ct <- array(rep(diag(n), T), dim = c(n, n, T))
  mu0 <- c0 <- alpha <- beta <- rep(NA, T)
  mu0_star <- c0_star <- alpha_star <- beta_star <- rep(NA, T)
  tau0 <- tau1 <- tau2 <- tau3 <- rep(NA, T)
  tau0_star <- tau1_star <- tau2_star <- tau3_star <- rep(NA, T)

  D <- ifelse(D == 0, 1, D)

  last_m <- m0
  last_C <- C0

  for (t in 1:T) {
    filter <- kernel_filter(outcome[t,], last_m, last_C, matrix(FF[, , t], n, r), G, D[, , t], W[, , t])

    at[, t] <- filter$at
    Rt[, , t] <- filter$Rt
    ft[, t] <- filter$ft
    Qt[, , t] <- filter$Qt

    mu0[t] <- filter$mu0
    c0[t] <- filter$c0
    alpha[t] <- filter$alpha
    beta[t] <- filter$beta
    mu0_star[t] <- filter$mu0_star
    c0_star[t] <- filter$c0_star
    alpha_star[t] <- filter$alpha_star
    beta_star[t] <- filter$beta_star

    tau0[t] <- filter$tau0
    tau1[t] <- filter$tau1
    tau2[t] <- filter$tau2
    tau3[t] <- filter$tau3
    tau0_star[t] <- filter$tau0_star
    tau1_star[t] <- filter$tau1_star
    tau2_star[t] <- filter$tau2_star
    tau3_star[t] <- filter$tau3_star
    # f_star[t,]     <- filter$f_star
    # Q_star[,,t]    <- filter$Q_star
    # At[,,t]        <- filter$At
    mt[, t] <- filter$mt
    Ct[, , t] <- filter$Ct

    last_m <- mt[, t]
    last_C <- Ct[, , t]
  }

  result <- list(
    "mt" = mt, "Ct" = Ct,
    "at" = at, "Rt" = Rt,
    "ft" = ft, "Qt" = Qt,
    "mu0" = mu0, "c0" = c0, "alpha" = alpha, "beta" = beta,
    "mu0_star" = mu0_star, "c0_star" = c0_star, "alpha_star" = alpha_star, "beta_star" = beta_star,
    "tau0" = tau0, "tau1" = tau1, "tau2" = tau2, "tau3" = tau3,
    "tau0_star" = tau0_star, "tau1_star" = tau1_star, "tau2_star" = tau2_star, "tau3_star" = tau3_star,
    "FF" = FF, "G" = G, "D" = D, "W" = W,
    "outcome" = outcome, "parms" = parms
  )
  return(result)
}

#' @export
normal_gamma_kernel <- list(
  "fit" = normal_gamma_fit,
  "filter" = normal_gamma_filter,
  "smoother" = generic_smoother,
  "pred" = normal_gamma_pred,
  "multi_var" = FALSE
)
