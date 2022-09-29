#' convert_NG_Normal
#'
#' Calculate the parameters of the Normal-Gamma that best approximates the given Normal distribuition.
#' The approximation is the best in the sense that it minimizes the KL divergence from the Normal to the Normal-Gamma.
#' In this approach, we suppose that the first entry of the multivariate normal represents the mean of the observed data and the second represent the log variance.
#'
#' @param ft vector: A vector representing the means from the normal distribution.
#' @param Qt matrix: A matrix representing the covariance matrix of the normal distribution.
#' @param parms list: A list of extra known parameters of the distribuition. Not used in this kernel.
#'
#' @return
#' @export
convert_NG_Normal=function(ft,Qt,parms=list()){
  mu0 <- ft[1,]+Qt[1,2]
  c0 <- exp(-ft[2,] - Qt[2,2] / 2)/Qt[1,1]
  helper=-3+3*sqrt(1+2*Qt[2,2]/3)
  alpha=1/helper
  beta <- alpha*exp(-ft[2,] - Qt[2,2] / 2)
  return(list('mu0'=mu0,'c0'=c0,'alpha'=alpha,'beta'=beta))
}

#' convert_Normal_NG
#'
#' Calculate the parameters of the Normal that best approximates the given Normal-Gamma distribuition.
#' The approximation is the best in the sense that it minimizes the KL divergence from the Normal-Gamma to the Normal.
#' In this approach, we suppose that the first entry of the multivariate normal represents the mean of the observed data and the second represent the log variance.
#'
#' @param conj_distr list: A list containg the parameters of the Normal-Gamma distribuition (mu0,c0,alpha,beta).
#' @param parms list: A list of extra known parameters of the distribuition. Not used in this kernel.
#'
#' @return
#' @export
convert_Normal_NG=function(conj_distr,parms=list()){
  f1 <- conj_distr$mu0
  f2 <- digamma(conj_distr$alpha)-log(conj_distr$beta)
  q1 <- conj_distr$beta / (conj_distr$c0 * (conj_distr$alpha-1))
  q2 <- trigamma(conj_distr$alpha)
  q12 <- 0
  ####################################################
  # Se estas equações forem usadas no lugar das que estão acima o ajuste fica bom, mesmo em casos relativamente extremos (priori vaga).
  #
  # q1 <- conj_distr$beta / (conj_distr$c0 * (conj_distr$alpha))
  # q2 <- 2*log(conj_distr$alpha)-2*digamma(conj_distr$alpha)
  ####################################################

  ft <- c(f1, f2)
  Qt <- matrix(c(q1, q12, q12, q2), byrow = F, ncol = 2)
  return(list('ft'=ft,'Qt'=Qt))
}

#' update_NG
#'
#' Calculate posterior parameter for the Normal-Gamma, assuming that the observed values came from a Normal model from which the prior distribuition for the mean and the precision have joint distribuition Normal-Gamma
#'
#' @param conj_prior list: A vector containing the parameters of the Normal-Gamma (mu0,c0,alpha,beta).
#' @param y vector: A vector containing the observations.
#' @param parms list: A list of extra known parameters of the distribuition. Not used in this kernel.
#'
#' @return
#' @export
update_NG=function(conj_prior,y,parms=list()){
  mu0   <- (conj_prior$c0 * conj_prior$mu0 + y) / (conj_prior$c0 + 1)
  c0    <- conj_prior$c0 + 1
  alpha <- conj_prior$alpha + 0.5
  beta  <- conj_prior$beta + 0.5 * conj_prior$c0 * ((conj_prior$mu0 - y)**2) / (conj_prior$c0 + 1)
  return(list('mu0'=mu0,'c0'=c0,'alpha'=alpha,'beta'=beta))
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
normal_gamma_pred <- function(conj_param,outcome=NULL,parms=list(), pred_cred = 0.95) {
  c0=conj_param$c0
  mu0=conj_param$mu0
  alpha=conj_param$alpha
  beta=conj_param$beta

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
    "log.like" = dt((outcome-mu)/sqrt(sigma2),nu,log=TRUE)
  )
}

#' @export
normal_gamma_kernel <- list(
  "conj_prior" = convert_NG_Normal,
  "conj_post" = convert_Normal_NG,
  "update" = update_NG,
  "smoother" = generic_smoother,
  "pred" = normal_gamma_pred,
  "offset" = ident_offset,
  "param_names" = function(y){c('mu0','c0','alpha','beta')},
  "multi_var" = FALSE
)
