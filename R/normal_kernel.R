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
#' @return The parameters of the conjugated distribuition (which is Normal) of the linear predictor.
#' @export
convert_dummy_prior <- function(ft, Qt, parms) {
  return(do.call(c, list(ft, Qt)))
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
#' @return The parameters of the Normal distribuition of the linear predictor.
#' @export
convert_dummy_post <- function(conj_prior, parms) {
  r_star <- length(conj_prior)
  r <- (-1 + sqrt(1 + 4 * r_star)) / 2
  ft <- conj_prior[1:r]
  Qt <- conj_prior[(r + 1):(r * r + r)] %>% matrix(r, r)
  return(list("ft" = ft, "Qt" = Qt))
}

#' update_Normal_dummy
#'
#' Calculate posterior parameter for the Normal, assuming that the observed values came from a Normal model from which the covariance is known and the prior distribuition for the mean vector have Normal distribuition
#'
#' @param conj_prior list: A vector containing the concentration parameters of the Normal.
#' @param y vector: A vector containing the observations.
#' @param parms list: A list of extra known parameters of the distribuition. For this kernel, parms should containg the covariance matrix parameter (Sigma) for the observational Normal model.
#'
#' @return The parameters of the posterior distribution.
#' @export
update_Normal_dummy <- function(conj_prior, y, parms) {
  r_star <- length(conj_prior)
  r <- (-1 + sqrt(1 + 4 * r_star)) / 2
  ft <- conj_prior[1:r]
  Qt <- conj_prior[(r + 1):(r * r + r)] %>% matrix(r, r)

  Sigma=parms$Sigma

  if(all(diag(Sigma)==0)){
    Qt=Sigma*0
    ft=y
  }else{
    Tau0 <- ginv(Qt)
    Tau1 <- ginv(parms$Sigma)

    Qt <- ginv(Tau0 + Tau1)
    ft <- Qt %*% (Tau0 %*% ft + Tau1 %*% y)
  }
  return(do.call(c, list(ft, Qt)))
}

#' normal_pred
#'
#' Calculate the values for the predictive distribuition given the values of the parameter of the conjugated distribuition of the linear predictor.
#' The data is assumed to have Normal distribuition with known variance and it's mean having distribuition Normal.
#' In this scenario, the marginal distribuition of the data is also Normal.
#'
#' @param conj_param List or data.frame: The parameters of the conjugated distribuitions of the linear predictor.
#' @param outcome Vector or matrix (optional): The observed values at the current time. Not used in this function.
#' @param parms List: A list of extra parameters for the model. For this function, it must contain the observational covariance matrix, Sigma.
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
#'
#' @examples
#'
# ft=c(0.5, -1 / 2, 0, -2.0)
# Qt=diag(4)
#' #
# conj_param <- c(ft,Qt)
#' #
# normal_pred(conj_param,parms=list('Sigma'=diag(4)))
normal_pred <- function(conj_param, outcome = NULL, parms = list(), pred_cred = 0.95) {
  if (is.null(dim(conj_param))) {
    conj_param <- conj_param %>% matrix(1, length(conj_param))
  }
  r <- dim(conj_param)[2]
  r <- (-1 + sqrt(1 + 4 * r)) / 2
  T <- dim(conj_param)[1]

  pred <- conj_param[, 1:r] %>% as.matrix()
  var.pred <- conj_param[(r + 1):(r * r + r)] %>%
    t() %>%
    array(c(r, r, T))
  icl.pred <- matrix(NA, T, r)
  icu.pred <- matrix(NA, T, r)
  for (t in 1:T) {
    mu <- pred[t, ]
    sigma2 <- var.pred[, , t] + parms$Sigma
    icl.pred[t, ] <- qnorm((1 - pred_cred) / 2) * sqrt(sigma2) + mu
    icu.pred[t, ] <- qnorm(1 - (1 - pred_cred) / 2) * sqrt(sigma2) + mu
  }

  list(
    "pred"     = pred,
    "var.pred" = var.pred,
    "icl.pred" = icl.pred,
    "icu.pred" = icu.pred
  )
}

#' normal_log_like
#'
#' Calculate the predictive log-likelyhood of the data given the values of the parameter of the conjugated distribuition of the linear predictor.
#' The data is assumed to have Normal distribuition with it's mean and precision having joint distribuition Normal-Gamma with parameters mu0, c0, alpha, beta.
#' In this scenario, the marginal distribuition of the data is t-student with 2*alpha degrees of freedoom, mu0 as the location parameter and (beta/alpha)*(1+1/c0) as the scale parameter.
#'
#' @param conj_param List or data.frame: The parameters of the conjugated distribuition (Normal-Gamma) of the linear predictor.
#' @param outcome Vector or matrix: The observed values at the current time.
#' @param parms List: A list of extra parameters for the model. For this function, it must contain the observational covariance matrix, Sigma.
#'
#' @return The predictive log-likelyhood of the observed data.
#' @export
#'
#' @examples
#'
#' ft <- c(0.5, -1 / 2, 0, -2.0)
#' Qt <- diag(4)
#'
#' conj_param <- c(ft, Qt)
#'
#' normal_log_like(conj_param, rnorm(1), parms = list("Sigma" = diag(4)))
normal_log_like <- function(conj_param, outcome, parms = list()) {
  if (is.null(dim(conj_param))) {
    conj_param <- conj_param %>% matrix(1, length(conj_param))
  }
  if (is.null(dim(outcome))) {
    outcome <- outcome %>% matrix(1, length(outcome))
  }
  r <- dim(outcome)[2]
  T <- dim(outcome)[1]

  log.like <- rep(NA, T)
  pred <- conj_param[, 1:r] %>% as.matrix()
  var.pred <- conj_param[(r + 1):(r * r + r)] %>%
    t() %>%
    array(c(r, r, T))
  for (t in 1:T) {
    mu <- pred[t, ]
    sigma2 <- var.pred[, , t] + parms$Sigma
    log.like[t] <- sum(dnorm(outcome[t, ], mu, sqrt(diag(sigma2)), log = TRUE))
  }
  log.like
}

create_normal_names <- function(y) {
  c(
    paste("ft_", 1:dim(y)[2], sep = ""),
    paste("Qt_",
      c(matrix(1:dim(y)[2], dim(y)[2], dim(y)[2])),
      c(matrix(1:dim(y)[2], dim(y)[2], dim(y)[2], byrow = TRUE)),
      sep = ""
    )
  )
}

#' @export
normal_kernel <- list(
  "conj_prior" = convert_dummy_prior,
  "conj_post" = convert_dummy_post,
  "update" = update_Normal_dummy,
  "smoother" = generic_smoother,
  "pred" = normal_pred,
  "log.like" = normal_log_like,
  "offset" = ident_offset,
  "link_function" = function(x) {
    x
  },
  "inv_link_function" = function(x) {
    x
  },
  "param_names" = create_normal_names,
  "multi_var" = TRUE
)
