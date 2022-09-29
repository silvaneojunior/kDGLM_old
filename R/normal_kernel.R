#' convert_dummy_prior
#'
#' Calculate the parameters of the Normal that best approximates the given Normal distribuition.
#' The approximation is the best in the sense that it minimizes the KL divergence from the Normal to the Normal.
#' This funciton is a dummy used for compatibility between kernels.
#'
#' @param ft vector: A vector representing the means from the normal distribution.
#' @param Qt matrix: A matrix representing the covariance matrix of the normal distribution.
#' @param parms list: A list of extra known parameters of the distribuition. Not used in this function.
#'
#' @importFrom rootSolve multiroot
#' @return
#' @export
convert_dummy_prior=function(ft,Qt,parms){
  return(do.call(c,list(ft,Qt)))
}

#' convert_dummy_post
#'
#' Calculate the parameters of the Normal that best approximates the given Normal distribuition.
#' The approximation is the best in the sense that it minimizes the KL divergence from the Normal to the Normal
#' This funciton is a dummy used for compatibility between kernels.
#'
#' @param conj_prior list: A vector containing the parameters of the Normal
#' @param parms list: A list of extra known parameters of the distribuition. Not used in this function.
#'
#' @return
#' @export
convert_dummy_post=function(conj_prior,parms){
  r_star=length(conj_prior)
  r=(-1+sqrt(1+4*r_star))/2
  ft <- conj_prior[1:r]
  Qt <- conj_prior[(r+1):(r*r+r)] %>% matrix(r,r)
  return(list('ft'=ft,'Qt'=Qt))
}

#' update_Normal_dummy
#'
#' Calculate posterior parameter for the Normal, assuming that the observed values came from a Normal model from which the covariance is known and the prior distribuition for the mean vector have Normal distribuition
#'
#' @param conj_prior list: A vector containing the concentration parameters of the Normal.
#' @param y vector: A vector containing the observations.
#' @param parms list: A list of extra known parameters of the distribuition. For this kernel, parms should containg the covariance matrix parameter (Sigma) for the observational Normal model.
#'
#' @return
#' @export
update_Normal_dummy=function(conj_prior,y,parms){
  r_star=length(conj_prior)
  r=(-1+sqrt(1+4*r_star))/2
  ft <- conj_prior[1:r]
  Qt <- conj_prior[(r+1):(r*r+r)] %>% matrix(r,r)

  Tau0=ginv(Qt)
  Tau1=ginv(parms$Sigma)

  Qt <- ginv(Tau0+Tau1)
  ft <- Qt%*%(Tau0%*%ft + Tau1 %*% y)
  return(do.call(c,list(ft,Qt)))
}

#' normal_pred
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
#' normal_pred(model)
normal_pred <- function(conj_param,outcome,parms=list(), pred_cred = 0.95) {
  if(is.null(dim(conj_param))){
    conj_param=conj_param %>% matrix(1,length(conj_param))
  }
  if(is.null(dim(conj_param))){
    outcome=outcome %>% matrix(1,length(outcome))
  }
  r=dim(outcome)[2]
  T=dim(outcome)[1]

  pred <- conj_param[,1:r] %>% as.matrix
  var.pred <- conj_param[(r+1):(r*r+r)] %>% t %>% array(c(r,r,T))
  icl.pred <- matrix(NA,T,r)
  icu.pred <- matrix(NA,T,r)
  log.like <- rep(NA,T)
  for(t in 1:T){
    mu=pred[t,]
    sigma2=var.pred[,,t]+parms$Sigma
    icl.pred[t,]=qnorm((1 - pred_cred) / 2)*sqrt(diag(sigma2))+mu
    print(icl.pred[t,])
    icu.pred[t,]=qnorm(1 - (1 - pred_cred) / 2)*sqrt(diag(sigma2))+mu
    # Placeholder: O valor está incorreto!
    log.like[t]=sum(dnorm(outcome[t,],mu,sqrt(diag(sigma2)),log=TRUE))
  }

  list(
    "pred"     = pred,
    "var.pred" = var.pred,
    "icl.pred" = icl.pred,
    "icu.pred" = icu.pred,
    "log.like" = log.like
  )
}

create_normal_names=function(y){
  c(paste('ft_',1:dim(y)[2],sep=''),
    paste('Qt_',
          c(matrix(1:dim(y)[2],dim(y)[2],dim(y)[2])),
          c(matrix(1:dim(y)[2],dim(y)[2],dim(y)[2],byrow=TRUE)),
          sep=''))}

#' @export
normal_kernel <- list(
  "conj_prior" = convert_dummy_prior,
  "conj_post" = convert_dummy_post,
  "update" = update_Normal_dummy,
  "smoother" = generic_smoother,
  "pred" = normal_pred,
  "offset" = ident_offset,
  "param_names" = create_normal_names,
  "multi_var" = TRUE
)
