#' Poisson
#'
#' Creates an outcome with Poisson distribuition with the chosen parameter.
#'
#' @param lambda character: Name of the linear preditor associated with the rate (mean) parameter of the Poisson distribuition. The parameter is treated as unknowed and equal to the exponential of the associated linear preditor.
#' @param outcome vector: Values of the observed data.
#' @param offset vector: The offset at each observation. Must have the same shape as outcome.
#'
#' @return A object of the class dlm_distr
#'
#' @export
#'
#' @examples
#'
#' # Poisson case
#' T <- 200
#' w <- (200 / 40) * 2 * pi
#' data <- rpois(T, 20 * (sin(w * 1:T / T) + 2))
#'
#' level <- polynomial_block(rate = 1, D = 1 / 0.95)
#' season <- harmonic_block(rate = 1, period = 40, D = 1 / 0.98)
#'
#' outcome <- Poisson(lambda = "rate", outcome = data)
#'
#' fitted_data <- fit_model(level, season, outcomes = outcome)
#' summary(fitted_data)
#'
#' show_fit(fitted_data, smooth = TRUE)$plot
Poisson <- function(lambda, outcome, offset = outcome**0) {
  family <- poisson_kernel
  t <- length(outcome)
  r <- 1
  convert_mat_default <- convert_mat_canom <- diag(r)

  distr <- list(
    var_names = c(lambda),
    family = family,
    r = r,
    t = t,
    offset = matrix(offset, t, r),
    outcome = matrix(outcome, t, r),
    convert_mat_canom = convert_mat_canom,
    convert_mat_default = convert_mat_default,
    parms = list(),
    name = "Poisson"
  )
  class(distr) <- "dlm_distr"

  return(distr)
}

#' convert_Gamma_Normal
#'
#' Calculate the parameters of the Gamma that best approximates the given log-Normal distribuition.
#' The approximation is the best in the sense that it minimizes the KL divergence from the log-Normal to the Gamma
#'
#' @param ft vector: A vector representing the means from the normal distribution.
#' @param Qt matrix: A matrix representing the covariance matrix of the normal distribution.
#' @param parms list: A list of extra known parameters of the distribuition. Not used in this kernel.
#'
#' @return The parameters of the conjugated distribuition of the linear predictor.
#' @export
convert_Gamma_Normal <- function(ft, Qt, parms) {
  # diag(Qt)=ifelse(diag(Qt)<0,0,diag(Qt))
  if (length(ft) > 1) {
    Qt <- diag(Qt)
    ft <- c(ft)
  }
  h <- -3 + 3 * sqrt(1 + 2 * Qt / 3)
  alpha <- (1 / h)
  beta <- alpha * exp(-ft - 0.5 * Qt)
  return(c(alpha, beta))
}

convert_Gamma_Normal_ord_1 <- function(ft, Qt, parms) {
  # diag(Qt)=ifelse(diag(Qt)<0,0,diag(Qt))
  h <- Qt
  alpha <- (1 / h)
  beta <- alpha * exp(-ft - 0.5 * Qt)
  return(list("alpha" = alpha, "beta" = beta))
}

convert_Gamma_Normal_RN_default <- function(ft, Qt, parms) {
  # diag(Qt)=ifelse(diag(Qt)<0,0,diag(Qt))
  root <- multiroot(function(alpha) {
    digamma(exp(alpha)) - log(exp(alpha)) + Qt / 2
  }, start = 0)
  alpha <- exp(root$root)
  beta <- alpha * exp(-ft - 0.5 * Qt)
  return(list("alpha" = alpha, "beta" = beta))
}

convert_Gamma_Normal_true <- function(ft, Qt, parms) {
  # diag(Qt)=ifelse(diag(Qt)<0,0,diag(Qt))
  root <- multiroot(function(alpha) {
    digamma(exp(alpha)) - log(exp(alpha)) + Qt / 2
  }, start = 0, rtol = 10**-20, maxiter = 1000)
  alpha <- exp(root$root)
  beta <- alpha * exp(-ft - 0.5 * Qt)
  return(list("alpha" = alpha, "beta" = beta))
}

#' convert_Normal_Gamma
#'
#' Calculate the parameters of the log-Normal that best approximates the given Gamma distribuition.
#' The approximation is the best in the sense that it minimizes the KL divergence from the Gamma to the log-Normal
#'
#' @param conj_prior list: A vector containing the parameters of the Gamma (alpha,beta).
#' @param parms list: A list of extra known parameters of the distribuition. Not used in this kernel.
#'
#' @return The parameters of the Normal distribuition of the linear predictor.
#' @export
convert_Normal_Gamma <- function(conj_prior, parms) {
  if (is.null(dim(conj_prior))) {
    r <- length(conj_prior) / 2
  } else {
    r <- dim(conj_prior)[2] / 2
  }
  alpha <- conj_prior[1:r]
  beta <- conj_prior[(r + 1):(2 * r)]
  ft <- digamma(alpha) - log(beta)
  Qt <- trigamma(alpha)
  if (length(alpha) > 1) {
    Qt <- diag(Qt)
    ft <- matrix(ft, r, 1)
  }
  # print(Qt)
  return(list("ft" = ft, "Qt" = Qt))
}

#' update_Gamma
#'
#' Calculate posterior parameter for the Gamma, assuming that the observed values came from a Poisson model from which the rate parameter (lambda) have prior distribuition Gamma.
#'
#' @param conj_prior list: A vector containing the parameters of the Gamma (alpha,beta).
#' @param y vector: A vector containing the observations.
#' @param parms list: A list of extra known parameters of the distribuition. Not used in this kernel.
#'
#' @return The parameters of the posterior distribution.
#' @export
update_Gamma <- function(conj_prior, y, parms) {
  if (is.null(dim(conj_prior))) {
    r <- length(conj_prior) / 2
  } else {
    r <- dim(conj_prior)[2] / 2
  }
  alpha <- conj_prior[1:r]
  beta <- conj_prior[(r + 1):(2 * r)]
  alpha <- alpha + y
  beta <- beta + 1
  return(c(alpha, beta))
}

#' convert_Gamma_Normal_LB
#'
#' Calculate the parameters of the Gamma that approximates the given log-Normal distribuition following the approach proposed in ref: Migon.
#'
#' @param ft vector: A vector representing the means from the normal distribution.
#' @param Qt matrix: A matrix representing the covariance matrix of the normal distribution.
#' @param parms list: A list of extra known parameters of the distribuition. Not used in this function.
#'
#' @return The parameters of the posterior distribution.
#' @export
convert_Gamma_Normal_LB <- function(ft, Qt, parms) {
  h <- -3 + 3 * sqrt(1 + 2 * Qt / 3)
  alpha <- (1 / h)
  beta <- alpha * exp(-ft + 0.5 * Qt)
  return(c(alpha, beta))
}

#' poisson_pred
#'
#' Calculate the values for the predictive distribuition given the values of the parameter of the conjugated distribuition of the linear predictor.
#' The data is assumed to have Poisson distribuition with it's mean having distribuition Gamma with shape parameter a e rate parameter b.
#' In this scenario, the marginal distribuition of the data is Negative Binomial with a as the dispersion parameter and b/(b+1) as the probability.
#'
#' @param conj_param List or data.frame: The parameters of the conjugated distribuitions of the linear predictor.
#' @param outcome Vector or matrix (optional): The observed values at the current time. Not used in this function.
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
#'   "a" = c(1:3),
#'   "b" = c(3:1)
#' )
#'
#' poisson_pred(conj_param)
poisson_pred <- function(conj_param, outcome = NULL, parms = list(), pred_cred = 0.95) {
  r <- 1
  a <- conj_param[1:r] %>% t()
  b <- conj_param[(1 + r):(2 * r)] %>% t()
  t <- length(a)
  pred <- matrix(NA, r, t)
  var.pred <- array(NA, c(1, 1, t))
  icl.pred <- matrix(NA, r, t)
  icu.pred <- matrix(NA, r, t)
  log.like <- rep(NA, t)
  flags <- b > 1e-40

  N <- 5000
  pred[, flags] <- a[flags] / b[flags]
  var.pred[, , flags] <- a[flags] * (b[flags] + 1) / (b[flags]**2)

  icl.pred[, flags] <- qnbinom((1 - pred_cred) / 2, a[flags], (b[flags] / (b[flags] + 1)))
  icu.pred[, flags] <- qnbinom(1 - (1 - pred_cred) / 2, a[flags], (b[flags] / (b[flags] + 1)))

  for (i in (1:t)[!flags]) {
    sample_lambda <- rgamma(N, a[i], b[i])
    sample_y <- rpois(N, sample_lambda)

    pred[, i] <- mean(sample_y)
    var.pred[, , i] <- var(sample_y)
    icl.pred[, i] <- quantile(sample_y, (1 - pred_cred) / 2)
    icu.pred[, i] <- quantile(sample_y, 1 - (1 - pred_cred) / 2)
  }

  return(list(
    "pred" = pred,
    "var.pred" = var.pred,
    "icl.pred" = icl.pred,
    "icu.pred" = icu.pred
  ))
}

#' poisson_log_like
#'
#' Calculate the predictive log-likelyhood of the data given the values of the parameter of the conjugated distribuition of the linear predictor.
#' The data is assumed to have Poisson distribuition with it's mean having distribuition Gamma with shape parameter a e rate parameter b.
#' In this scenario, the marginal distribuition of the data is Negative Binomial with a as the dispersion parameter and b/(b+1) as the probability.
#'
#' @param conj_param List or data.frame: The parameters of the conjugated distribuition (Gamma) of the linear predictor.
#' @param outcome Vector or matrix: The observed values at the current time.
#' @param parms List (optional): A list of extra parameters for the model. Not used in this function.
#'
#' @return The predictive log-likelyhood of the observed data.
#' @export
#'
#' @examples
#'
#' conj_param <- data.frame(
#'   "a" = c(1:3),
#'   "b" = c(3:1)
#' )
#'
#' poisson_log_like(conj_param, rpois(3, 1))
poisson_log_like <- function(conj_param, outcome, parms = list()) {
  r <- 1
  a <- conj_param[1:r] %>% t()
  b <- conj_param[(1 + r):(2 * r)] %>% t()
  t <- length(a)
  pred <- matrix(NA, r, t)
  var.pred <- array(NA, c(1, 1, t))
  icl.pred <- matrix(NA, r, t)
  icu.pred <- matrix(NA, r, t)
  log.like <- rep(NA, t)
  flags <- b > 1e-40

  N <- 5000
  log.like[flags] <- dnbinom(outcome[flags], a[flags], (b[flags] / (b[flags] + 1)), log = TRUE)

  for (i in (1:t)[!flags]) {
    sample_lambda <- rgamma(N, a[i], b[i])
    sample_y <- rpois(N, sample_lambda)
    log.like.list <- dpois(outcome[i], sample_lambda, log = TRUE)
    max.log.like <- max(log.like.list)
    like.list <- exp(log.like.list - max.log.like)
    log.like[i] <- log(mean(like.list)) + max.log.like
  }
  return(log.like)
}

#' @export
poisson_kernel <- list(
  "conj_prior" = convert_Gamma_Normal,
  "conj_post" = convert_Normal_Gamma,
  "update" = update_Gamma,
  "smoother" = generic_smoother,
  "pred" = poisson_pred,
  "log.like" = poisson_log_like,
  "offset" = log_offset,
  "link_function" = log,
  "inv_link_function" = exp,
  "param_names" = function(y) {
    c(paste0("alpha_", 1:dim(y)[2]), paste("beta_", 1:dim(y)[2]))
  },
  "multi_var" = FALSE
)

#' @export
poisson_kernel_ord1 <- list(
  "conj_prior" = convert_Gamma_Normal_ord_1,
  "conj_post" = convert_Normal_Gamma,
  "update" = update_Gamma,
  "smoother" = generic_smoother,
  "pred" = poisson_pred,
  "log.like" = poisson_log_like,
  "link_function" = log,
  "inv_link_function" = exp,
  "offset" = log_offset,
  "param_names" = function(y) {
    c("alpha", "beta")
  },
  "multi_var" = FALSE
)

#' @export
poisson_kernel_RN_default <- list(
  "conj_prior" = convert_Gamma_Normal_RN_default,
  "conj_post" = convert_Normal_Gamma,
  "update" = update_Gamma,
  "smoother" = generic_smoother,
  "pred" = poisson_pred,
  "log.like" = poisson_log_like,
  "offset" = log_offset,
  "link_function" = log,
  "inv_link_function" = exp,
  "param_names" = function(y) {
    c("alpha", "beta")
  },
  "multi_var" = FALSE
)

#' @export
poisson_kernel_true <- list(
  "conj_prior" = convert_Gamma_Normal_true,
  "conj_post" = convert_Normal_Gamma,
  "update" = update_Gamma,
  "smoother" = generic_smoother,
  "pred" = poisson_pred,
  "log.like" = poisson_log_like,
  "offset" = log_offset,
  "link_function" = log,
  "inv_link_function" = exp,
  "param_names" = function(y) {
    c("alpha", "beta")
  },
  "multi_var" = FALSE
)

#' @export
poisson_lb_kernel <- list(
  "conj_prior" = convert_Gamma_Normal_LB,
  "conj_post" = convert_Normal_Gamma,
  "update" = update_Gamma,
  "smoother" = generic_smoother,
  "pred" = poisson_pred,
  "log.like" = poisson_log_like,
  "offset" = log_offset,
  "link_function" = log,
  "inv_link_function" = exp,
  "param_names" = function(y) {
    c("alpha", "beta")
  },
  "multi_var" = FALSE
)
