#' Normal
#'
#' Creates an outcome with Normal distribuition with the chosen parameters (can only specify 2,).
#'
#' @param mu character: Name of the linear preditor associated with the mean parameter of the Normal distribuition. The parameter is treated as unknowed and equal to the associated linear preditor.
#' @param tau character: Name of the linear preditor associated with the precision parameter of the Normal distribuition. The parameter is treated as unknowed and equal to the exponential of the associated linear preditor. It cannot be specified with sigma or sigma2
#' @param sigma character: Name of the linear preditor associated with the scale parameter of the Normal distribuition. The parameter is treated as unknowed and equal to the exponential of the associated linear preditor. It cannot be specified with tau or sigma2.
#' @param sigma2 character or numeric: Name of the linear preditor associated with the variance parameter of the Normal distribuition. If numeric, this parameter is treated as knowed and equal to the value passed. If a character, the parameter is treated as unknowed and equal to the exponential of the associated linear preditor. It cannot be specified with sigma or tau.
#' @param outcome vector: Values of the observed data.
#' @param offset vector: The offset at each observation. Must have the same shape as outcome.
#'
#' @return A object of the class dlm_distr
#' @export
#'
#' @examples
#'
#' # Normal case
#' T <- 200
#' mu <- rnorm(T, 0, 0.1)
#' data <- rnorm(T, cumsum(mu))
#'
#' level <- polynomial_block(
#'   mu = 1,
#'   D = 1 / 0.95
#' )
#' variance <- polynomial_block(
#'   sigma2 = 1
#' )
#'
#' # Known variance
#' outcome <- Normal(mu = "mu", sigma2 = 1, outcome = data)
#'
#' fitted_data <- fit_model(level, outcomes = outcome)
#' summary(fitted_data)
#'
#' show_fit(fitted_data, smooth = TRUE)$plot
#'
#' # Unknown variance
#' outcome <- Normal(mu = "mu", sigma2 = "sigma2", outcome = data)
#'
#' fitted_data <- fit_model(level, variance, outcomes = outcome)
#' summary(fitted_data)
#'
#' show_fit(fitted_data, smooth = TRUE)$plot
Normal <- function(mu, tau = NA, sigma = NA, sigma2 = NA, outcome, offset = outcome**0) {
  t <- length(outcome)
  r <- 1
  flags <- is.na(c(tau, sigma, sigma2))
  if (all(flags)) {
    stop("Error: Scale not specified.")
  }
  if (sum(!flags) > 1) {
    stop("Error: Scale specified in more than one value.")
  }
  val_type <- c("tau", "sigma", "sigma2")[!flags]
  val <- c(tau, sigma, sigma2)[!flags]
  if (is.numeric(val)) {
    var_names <- c(mu)
    names(var_names) <- c("mu")
    family <- normal_kernel
    parms <- list(Sigma = if (val_type == "tau") {
      1 / val
    } else if (val_type == "sigma") {
      val**2
    } else {
      val
    })
    convert_mat_canom <- convert_mat_default <- diag(1)
  } else {
    var_names <- c(mu, val)
    names(var_names) <- c("mu", val_type)
    family <- normal_gamma_cor_kernel
    parms <- list()
    convert_mat_canom <- if (val_type == "tau") {
      diag(2)
    } else if (val_type == "sigma") {
      matrix(c(1, 0, 0, -2), 2, 2)
    } else {
      matrix(c(1, 0, 0, -1), 2, 2)
    }
    convert_mat_default <- if (val_type == "tau") {
      diag(2)
    } else if (val_type == "sigma") {
      matrix(c(1, 0, 0, -0.5), 2, 2)
    } else {
      matrix(c(1, 0, 0, -1), 2, 2)
    }
  }
  distr <- list(
    var_names = var_names,
    family = family,
    r = 1,
    k = length(var_names),
    t = t,
    offset = matrix(offset, t, r),
    outcome = matrix(outcome, t, r),
    convert_mat_canom = convert_mat_canom,
    convert_mat_default = convert_mat_default,
    parms = parms,
    name = "Normal"
  )
  class(distr) <- "dlm_distr"

  return(distr)
}

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
#' @return The parameters of the conjugated distribuition of the linear predictor.
#' @export
convert_NG_Normal <- function(ft, Qt, parms = list()) {
  ft <- matrix(ft, 2, 1)
  mu0 <- ft[1, ] + Qt[1, 2]
  c0 <- exp(-ft[2, ] - Qt[2, 2] / 2) / (Qt[1, 1])
  helper <- -3 + 3 * sqrt(1 + 2 * Qt[2, 2] / 3)
  # helper=Qt[2,2]
  alpha <- 1 / helper
  beta <- alpha * exp(-ft[2, ] - Qt[2, 2] / 2)
  return(list("mu0" = mu0, "c0" = c0, "alpha" = alpha, "beta" = beta))
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
#' @return The parameters of the Normal distribuition of the linear predictor.
#' @export
convert_Normal_NG <- function(conj_distr, parms = list()) {
  # q12 <- conj_distr$mu0*(digamma(conj_distr$alpha)-log(conj_distr$beta))/(1-log(conj_distr$alpha)+log(conj_distr$beta))
  f1 <- conj_distr$mu0 #-q12
  f2 <- digamma(conj_distr$alpha) - log(conj_distr$beta)
  # f2 <- log(conj_distr$alpha)-log(conj_distr$beta)
  q1 <- conj_distr$beta / (conj_distr$c0 * (conj_distr$alpha - 1))
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
  return(list("ft" = ft, "Qt" = Qt))
}

convert_Normal_NG_cor <- function(conj_distr, parms = list()) {
  # q12 <- conj_distr$mu0*(digamma(conj_distr$alpha)-log(conj_distr$beta))/(1-log(conj_distr$alpha)+log(conj_distr$beta))
  f1 <- conj_distr$mu0 #-q12
  f2 <- digamma(conj_distr$alpha) - log(conj_distr$beta)
  # f2 <- log(conj_distr$alpha)-log(conj_distr$beta)
  q1 <- conj_distr$beta / (conj_distr$c0 * conj_distr$alpha)
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
  return(list("ft" = ft, "Qt" = Qt))
}

#' update_NG
#'
#' Calculate posterior parameter for the Normal-Gamma, assuming that the observed values came from a Normal model from which the prior distribuition for the mean and the precision have joint distribuition Normal-Gamma
#'
#' @param conj_prior list: A vector containing the parameters of the Normal-Gamma (mu0,c0,alpha,beta).
#' @param y vector: A vector containing the observations.
#' @param parms list: A list of extra known parameters of the distribuition. Not used in this kernel.
#'
#' @return The parameters of the posterior distribution.
#' @export
update_NG <- function(conj_prior, y, parms = list()) {
  mu0 <- (conj_prior$c0 * conj_prior$mu0 + y) / (conj_prior$c0 + 1)
  c0 <- conj_prior$c0 + 1
  alpha <- conj_prior$alpha + 0.5
  beta <- conj_prior$beta + 0.5 * conj_prior$c0 * ((conj_prior$mu0 - y)**2) / (conj_prior$c0 + 1)
  return(list("mu0" = mu0, "c0" = c0, "alpha" = alpha, "beta" = beta))
}

#' normal_gamma_pred
#'
#' Calculate the values for the predictive distribuition given the values of the parameter of the conjugated distribuition of the linear predictor.
#' The data is assumed to have Normal distribuition with it's mean and precision having joint distribuition Normal-Gamma with parameters mu0, c0, alpha, beta.
#' In this scenario, the marginal distribuition of the data is t-student with 2*alpha degrees of freedoom, mu0 as the location parameter and (beta/alpha)*(1+1/c0) as the scale parameter.
#'
#' @param conj_param List or data.frame: The parameters of the conjugated distribuition (Normal-Gamma) of the linear predictor.
#' @param outcome Vector or matrix (optional): The observed values at the current time.
#' @param parms List (optional): A list of extra parameters for the model. Not used in this function.
#' @param pred_cred Numeric: the desired credibility for the credibility interval.
#'
#' @return A list containing the following values:
#' \itemize{
#'    \item pred vector/matrix: the mean of the predictive distribuition of a next observation. Same type and shape as the parameter in model.
#'    \item var.pred vector/matrix: the variance of the predictive distribuition of a next observation. Same type and shape as the parameter in model.
#'    \item icl.pred vector/matrix: the percentile of 100*((1-pred_cred)/2)% of the predictive distribuition of a next observation. Same type and shape as the parameter in model.
#'    \item icu.pred vector/matrix: the percentile of 100*(1-(1-pred_cred)/2)% of the predictive distribuition of a next observation. Same type and shape as the parameter in model.
#' }
#'
#' @export
#'
#' @examples
#'
#' conj_param <- data.frame(
#'   "mu0" = rep(0, 3),
#'   "c0" = rep(1, 3),
#'   "alpha" = rep(3, 3),
#'   "beta" = rep(3, 3)
#' )
#'
#' normal_gamma_pred(conj_param)
normal_gamma_pred <- function(conj_param, outcome = NULL, parms = list(), pred_cred = 0.95) {
  c0 <- conj_param$c0
  mu0 <- conj_param$mu0
  alpha <- conj_param$alpha
  beta <- conj_param$beta

  mu <- mu0
  nu <- 2 * alpha
  sigma2 <- (beta / alpha) * (1 + 1 / c0)

  pred <- mu %>%
    as.matrix() %>%
    t()
  var.pred <- (beta / (alpha - 1)) * (1 + 1 / c0) %>%
    as.matrix() %>%
    t()

  icl.pred <- qt((1 - pred_cred) / 2, nu) * sqrt(sigma2) + mu
  icu.pred <- qt(1 - (1 - pred_cred) / 2, nu) * sqrt(sigma2) + mu

  list(
    "pred"     = pred,
    "var.pred" = var.pred,
    "icl.pred" = icl.pred,
    "icu.pred" = icu.pred
  )
}

#' normal_gamma_log_like
#'
#' Calculate the predictive log-likelyhood of the data given the values of the parameter of the conjugated distribuition of the linear predictor.
#' The data is assumed to have Normal distribuition with it's mean and precision having joint distribuition Normal-Gamma with parameters mu0, c0, alpha, beta.
#' In this scenario, the marginal distribuition of the data is t-student with 2*alpha degrees of freedoom, mu0 as the location parameter and (beta/alpha)*(1+1/c0) as the scale parameter.
#'
#' @param conj_param List or data.frame: The parameters of the conjugated distribuition (Normal-Gamma) of the linear predictor.
#' @param outcome Vector or matrix: The observed values at the current time.
#' @param parms List (optional): A list of extra parameters for the model. Not used in this function.
#'
#' @return The predictive log-likelyhood of the observed data.
#' @export
#'
#' @examples
#'
#' conj_param <- data.frame(
#'   "mu0" = rep(0, 3),
#'   "c0" = rep(1, 3),
#'   "alpha" = rep(3, 3),
#'   "beta" = rep(3, 3)
#' )
#'
#' normal_gamma_log_like(conj_param, rnorm(3))
normal_gamma_log_like <- function(conj_param, outcome, parms = list()) {
  c0 <- conj_param$c0
  mu0 <- conj_param$mu0
  alpha <- conj_param$alpha
  beta <- conj_param$beta

  mu <- mu0
  nu <- 2 * alpha
  sigma2 <- (beta / alpha) * (1 + 1 / c0)

  return(dt((outcome - mu) / sqrt(sigma2), nu, log = TRUE))
}

#' @export
normal_gamma_kernel <- list(
  "conj_prior" = convert_NG_Normal,
  "conj_post" = convert_Normal_NG,
  "update" = update_NG,
  "smoother" = generic_smoother,
  "pred" = normal_gamma_pred,
  "log.like" = normal_gamma_log_like,
  "offset" = ident_log_offset,
  "link_function" = function(x) {
    rbind(x[1, ], log(x[2, ]))
  },
  "inv_link_function" = function(x) {
    rbind(x[1, ], exp(x[2, ]))
  },
  "param_names" = function(y) {
    c("mu0", "c0", "alpha", "beta")
  },
  "multi_var" = FALSE
)

#' @export
normal_gamma_cor_kernel <- list(
  "conj_prior" = convert_NG_Normal,
  "conj_post" = convert_Normal_NG_cor,
  "update" = update_NG,
  "smoother" = generic_smoother,
  "pred" = normal_gamma_pred,
  "log.like" = normal_gamma_log_like,
  "offset" = ident_log_offset,
  "link_function" = function(x) {
    rbind(x[1, ], log(x[2, ]))
  },
  "inv_link_function" = function(x) {
    rbind(x[1, ], exp(x[2, ]))
  },
  "param_names" = function(y) {
    c("mu0", "c0", "alpha", "beta")
  },
  "multi_var" = FALSE
)
