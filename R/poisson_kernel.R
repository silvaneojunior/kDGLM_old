#' convert_Gamma_Normal
#'
#' Calculate the parameters of the Gamma that best approximates the given log-Normal distribuition.
#' The approximation is the best in the sense that it minimizes the KL divergence from the log-Normal to the Gamma
#'
#' @param ft vector: A vector representing the means from the normal distribution.
#' @param Qt matrix: A matrix representing the covariance matrix of the normal distribution.
#' @param parms list: A list of extra known parameters of the distribuition. Not used in this kernel.
#'
#' @return
#' @export
convert_Gamma_Normal=function(ft,Qt,parms){
  # diag(Qt)=ifelse(diag(Qt)<0,0,diag(Qt))
  h=-3+3*sqrt(1+2*Qt/3)
  alpha <- (1 / h)
  beta <- alpha*exp(-ft - 0.5 * Qt)
  return(list('alpha'=alpha,'beta'=beta))
}

#' convert_Normal_Gamma
#'
#' Calculate the parameters of the log-Normal that best approximates the given Gamma distribuition.
#' The approximation is the best in the sense that it minimizes the KL divergence from the Gamma to the log-Normal
#'
#' @param conj_prior list: A vector containing the parameters of the Gamma (alpha,beta).
#' @param parms list: A list of extra known parameters of the distribuition. Not used in this kernel.
#'
#' @return
#' @export
convert_Normal_Gamma=function(conj_prior,parms){
  ft <- digamma(conj_prior$alpha) - log(conj_prior$beta)
  Qt <- trigamma(conj_prior$alpha)
  return(list('ft'=ft,'Qt'=Qt))
}

#' update_Gamma
#'
#' Calculate posterior parameter for the Gamma, assuming that the observed values came from a Poisson model from which the rate parameter (lambda) have prior distribuition Gamma.
#'
#' @param conj_prior list: A vector containing the parameters of the Gamma (alpha,beta).
#' @param y vector: A vector containing the observations.
#' @param parms list: A list of extra known parameters of the distribuition. Not used in this kernel.
#'
#' @return
#' @export
update_Gamma=function(conj_prior,y,parms){
  alpha <- conj_prior$alpha + y
  beta <- conj_prior$beta + 1
  return(list('alpha'=alpha,'beta'=beta))
}

#' convert_Gamma_Normal_LB
#'
#' Calculate the parameters of the Gamma that approximates the given log-Normal distribuition following the approach proposed in ref: Migon.
#'
#' @param ft vector: A vector representing the means from the normal distribution.
#' @param Qt matrix: A matrix representing the covariance matrix of the normal distribution.
#' @param parms list: A list of extra known parameters of the distribuition. Not used in this function.
#'
#' @return
#' @export
convert_Gamma_Normal_LB=function(ft,Qt,parms){
  h=-3+3*sqrt(1+2*Qt/3)
  alpha <- (1 / h)
  beta <- alpha*exp(-ft + 0.5 * Qt)
  return(list('alpha'=alpha,'beta'=beta))
}

#' poisson_pred
#'
#' Calculate the values for the predictive distribuition given the values of the parameter of the conjugated distribuition of the linear predictor.
#' It's worth noting that, since the conjugated distribuition of the linear predictior is Gamma (see Ref. Raíra), then the predictive distribuition is Negative Binomial.
#'
#' @param model List: contaning the parameters of the conjugated distribuition for the linear predictor (may be based on the prior, filtered or smoothed distribuition).
#' @param pred_cred Numeric: the desired credibility for the credibility interval
#'
#' @return A list containing the following values:
#' \itemize{
#'    \item pred vector/matrix: the mean of the predictive distribuition of a next observation. Same type and shape as the parameter in model.
#'    \item var.pred vector/matrix: the variance of the predictive distribuition of a next observation. Same type and shape as the parameter in model.
#'    \item icl.pred vector/matrix: the percentile of 100*((1-pred_cred)/2)% of the predictive distribuition of a next observation. Same type and shape as the parameter in model.
#'    \item icu.pred vector/matrix: the percentile of 100*(1-(1-pred_cred)/2)% of the predictive distribuition of a next observation. Same type and shape as the parameter in model.
#'    \item log.like vector/matrix: the log likelyhood of the observed data based on the predictive distribuition.
#' }
#'
#' @export
#'
#' @examples
#' # A fitted model shoulb be used as argument, but you can also pass only the parameter themselves.
#'
#' model <- list(
#'   "a" = c(1:3),
#'   "b" = c(3:1)
#' )
#'
#' poisson_pred(model)
poisson_pred <- function(conj_param,outcome=NULL,parms=list(), pred_cred = 0.95) {
  a <- conj_param$a %>% t
  b <- conj_param$b %>% t

  if(any(b < 10**-40)){
    b=ifelse(b >= 10**-40,b,NA)
    warning('The beta parameter for the predictive distribution is very low (<1e-40) at some times. Predicition for those times are unviable.')
  }

  pred <- a / b
  var.pred <- a * (b + 1) / (b)^2

  icl.pred <- qnbinom((1 - pred_cred) / 2, a, (b / (b + 1)))
  icu.pred <- qnbinom(1 - (1 - pred_cred) / 2, a, (b / (b + 1)))

  if (all(!is.null(outcome))) {
    log.like <- dnbinom(outcome, a, (b / (b + 1)), log = TRUE)
  } else {
    log.like <- NA
  }
  return(list(
    "pred" = pred,
    "var.pred" = var.pred,
    "icl.pred" = icl.pred,
    "icu.pred" = icu.pred,
    "log.like" = log.like
  ))
}

#' @export
poisson_kernel <- list(
  "conj_prior" = convert_Gamma_Normal,
  "conj_post" = convert_Normal_Gamma,
  "update" = update_Gamma,
  "smoother" = generic_smoother,
  "pred" = poisson_pred,
  "offset" = log_offset,
  "param_names" = function(y){c('alpha','beta')},
  "multi_var" = FALSE
)

#' @export
poisson_lb_kernel <- list(
  "conj_prior" = convert_Gamma_Normal_LB,
  "conj_post" = convert_Normal_Gamma,
  "update" = update_Gamma,
  "smoother" = generic_smoother,
  "pred" = poisson_pred,
  "offset" = log_offset,
  "param_names" = function(y){c('alpha','beta')},
  "multi_var" = FALSE
)
