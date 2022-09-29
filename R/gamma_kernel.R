#' convert_IGamma_Normal
#'
#' Calculate the parameters of the Inverse-Gamma that best approximates the given log-Normal distribuition.
#' The approximation is the best in the sense that it minimizes the KL divergence from the log-Normal to the Inverse-Gamma
#'
#' @param ft vector: A vector representing the means from the normal distribution.
#' @param Qt matrix: A matrix representing the covariance matrix of the normal distribution.
#' @param parms list: A list of extra known parameters of the distribuition. Not used in this function.
#'
#' @return
#' @export
convert_IGamma_Normal=function(ft,Qt,parms){
  alpha <- 1 / (-3+3*sqrt(1+2*Qt/3))
  beta <- alpha * exp(ft - Qt / 2)
  return(list('alpha'=alpha,'beta'=beta))
}

#' convert_Normal_IGamma
#'
#' Calculate the parameters of the log-Normal that best approximates the given Inverse-Gamma distribuition.
#' The approximation is the best in the sense that it minimizes the KL divergence from the Inverse-Gamma to the log-Normal
#'
#' @param conj_prior list: A vector containing the parameters of the Inverse-Gamma (alpha,beta).
#' @param parms list: A list of extra known parameters of the distribuition. Not used in this function.
#'
#' @return
#' @export
convert_Normal_IGamma=function(conj_prior,parms){
  ft <- -digamma(conj_prior$alpha) + log(conj_prior$beta)
  Qt <- trigamma(conj_prior$alpha)
  return(list('ft'=ft,'Qt'=Qt))
}

#' update_IGamma
#'
#' Calculate posterior parameter for the Inverser-Gamma, assuming that the observed values came from a Gamma model from which the shape parameter (phi) is known and the mean (mu) have prior distribuition Inverse-Gamma.
#'
#' @param conj_prior list: A vector containing the parameters of the Inverse-Gamma (alpha,beta).
#' @param y vector: A vector containing the observations.
#' @param parms list: A list of extra known parameters of the distribuition. For this kernel, parms should containg the shape parameter (phi) for the observational gamma model.
#'
#' @return
#' @export
update_IGamma=function(conj_prior,y,parms){
  alpha <- conj_prior$alpha + parms$phi
  beta <- conj_prior$beta + y * parms$phi
  return(list('alpha'=alpha,'beta'=beta))
}

#' convert_IGamma_Normal_LB
#'
#' Calculate the parameters of the Inverse-Gamma that approximates the given log-Normal distribuition following the approach proposed in ref: Migon.
#'
#' @param ft vector: A vector representing the means from the normal distribution.
#' @param Qt matrix: A matrix representing the covariance matrix of the normal distribution.
#' @param parms list: A list of extra known parameters of the distribuition. Not used in this function.
#'
#' @return
#' @export
convert_IGamma_Normal_LB=function(ft,Qt,parms){
  mu <- exp(ft + Qt / 2)
  var <- exp(Qt - 1) * exp(2 * ft + Qt)

  alpha <- (mu**2) / var + 2
  beta <- mu * (alpha - 1)
  return(list('alpha'=alpha,'beta'=beta))
}

#' convert_Normal_IGamma_LB
#'
#' Calculate the parameters of the log-Normal that best approximates the given Inverse-Gamma distribuition following the approach proposed in ref: Migon.
#'
#' @param conj_prior list: A vector containing the parameters of the Inverse-Gamma (alpha,beta).
#' @param parms list: A list of extra known parameters of the distribuition. Not used in this function.
#'
#' @return
#' @export
convert_Normal_IGamma_LB=function(conj_prior,parms){
  mu <- conj_prior$beta / (conj_prior$alpha - 1)
  var <- (conj_prior$beta**2) / (((conj_prior$alpha - 1)**2) * (conj_prior$alpha - 2))

  Qt <- log(var / (mu**2)) + 1
  ft <- log(mu) - Qt / 2
  return(list('ft'=ft,'Qt'=Qt))
}

#' gamma_pred
#'
#' Calculate the values for the predictive distribuition given the values of the parameter of the conjugated distribuition of the linear predictor.
#' It's worth noting that the credibility interval is calculated as +-2 standard deviations.
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
#'    \item log.like vector/matrix: the log likelyhood of the observed data based on the predictive distribuition.
#' }
#' @export
#' @importFrom extraDistr qbetapr
#' @importFrom extraDistr dbetapr
#'
#' @examples
#' # A fitted model shoulb be used as argument, but you can also pass only the parameter themselves.
#'
#' model <- list(
#'   "outcome" = 2,
#'   "tau0" = 3,
#'   "tau1" = 4,
#'   "parms" = list("alpha" = 1)
#' )
#'
#' gamma_pred(model)
gamma_pred <- function(conj_param,outcome=NULL,parms=list(), pred_cred = 0.95) {
  phi <- parms$phi

  alpha <- (phi * conj_param$alpha) %>% t
  beta <- (phi * conj_param$beta) %>% t
  outcome_list <- list(
    "pred"     = beta / (alpha - 1),
    "var.pred" = ((alpha / (beta - 1))**2) * (alpha + phi - 1) / ((alpha - 2) * phi),
    "icl.pred" = qbetapr((1 - pred_cred) / 2, phi, alpha, beta / phi),
    "icu.pred" = qbetapr(1 - (1 - pred_cred) / 2, phi, alpha, beta / phi)
  )
  if (all(!is.null(outcome))) {
    outcome_list$log.like <- dbetapr(outcome, phi, alpha, beta / phi, log = TRUE)
  }
  return(outcome_list)
}

#' @export
gamma_kernel <- list(
  "conj_prior" = convert_IGamma_Normal,
  "conj_post" = convert_Normal_IGamma,
  "update" = update_IGamma,
  "smoother" = generic_smoother,
  "pred" = gamma_pred,
  "offset" = log_offset,
  "param_names" = function(y){c('alpha','beta')},
  "multi_var" = FALSE
)

#' @export
gamma_lb_kernel <- list(
  "conj_prior" = convert_IGamma_Normal_LB,
  "conj_post" = convert_Normal_IGamma_LB,
  "update" = update_IGamma,
  "smoother" = generic_smoother,
  "pred" = gamma_pred,
  "offset" = log_offset,
  "param_names" = function(y){c('alpha','beta')},
  "multi_var" = FALSE
)
