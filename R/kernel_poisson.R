#### Default Method ####

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
  distr <- list()
  t <- length(outcome)
  r <- k <- 1
  convert_mat_default <- convert_mat_canom <- diag(r)

  distr <- list(
    conj_prior = convert_Poisson_Normal,
    conj_post = convert_Normal_Poisson,
    update = update_Poisson,
    smoother = generic_smoother,
    calc_pred = poisson_pred,
    apply_offset = function(ft, Qt, offset) {
      t <- if.null(dim(ft)[2], 1)
      offset <- matrix(offset, t, r)

      list("ft" = ft + log(t(offset)), "Qt" = Qt)
    },
    link_function = log,
    inv_link_function = exp,
    param_names = function(y) {
      c(paste0("alpha_", 1:dim(y)[2]), paste("beta_", 1:dim(y)[2]))
    },
    var_names = c(lambda),
    family = family,
    r = r,
    k = k,
    l = k,
    t = t,
    offset = matrix(offset, t, r),
    outcome = matrix(outcome, t, r),
    convert_mat_canom = convert_mat_canom,
    convert_mat_default = convert_mat_default,
    convert_canom_flag = FALSE,
    parms = list(),
    name = "Poisson"
  )
  class(distr) <- "dlm_distr"
  distr$alt_method <- FALSE

  return(distr)
}

#' convert_Poisson_Normal
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
convert_Poisson_Normal <- function(ft, Qt, parms) {
  # diag(Qt)=ifelse(diag(Qt)<0,0,diag(Qt))
  if (length(ft) > 1) {
    Qt <- diag(Qt)
    ft <- c(ft)
  }
  h <- -3 + 3 * sqrt(1 + 2 * Qt / 3)
  alpha <- (1 / h)
  beta <- alpha * exp(-ft - 0.5 * Qt)
  return(list("alpha" = alpha, "beta" = beta))
}

#' convert_Normal_Poisson
#'
#' Calculate the parameters of the log-Normal that best approximates the given Gamma distribuition.
#' The approximation is the best in the sense that it minimizes the KL divergence from the Gamma to the log-Normal
#'
#' @param conj_prior list: A vector containing the parameters of the Gamma (alpha,beta).
#' @param parms list: A list of extra known parameters of the distribuition. Not used in this kernel.
#'
#' @return The parameters of the Normal distribuition of the linear predictor.
#' @export
convert_Normal_Poisson <- function(conj_prior, parms) {
  alpha <- conj_prior$alpha
  beta <- conj_prior$beta
  ft <- digamma(alpha) - log(beta)
  Qt <- trigamma(alpha)
  if (length(alpha) > 1) {
    Qt <- diag(Qt)
    ft <- matrix(ft, r, 1)
  }
  # print(Qt)
  return(list("ft" = ft, "Qt" = Qt))
}

#' update_Poisson
#'
#' Calculate posterior parameter for the Gamma, assuming that the observed values came from a Poisson model from which the rate parameter (lambda) have prior distribuition Gamma.
#'
#' @param conj_prior list: A vector containing the parameters of the Gamma (alpha,beta).
#' @param ft vector: A vector representing the means from the normal distribution. Not used in the default method.
#' @param Qt matrix: A matrix representing the covariance matrix of the normal distribution. Not used in the default method.
#' @param y vector: A vector containing the observations.
#' @param parms list: A list of extra known parameters of the distribuition. Not used in this kernel.
#'
#' @return The parameters of the posterior distribution.
#' @export
update_Poisson <- function(conj_prior, ft, Qt, y, parms) {
  alpha <- conj_prior$alpha
  beta <- conj_prior$beta
  alpha <- alpha + y
  beta <- beta + 1
  return(list("alpha" = alpha, "beta" = beta))
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
#'    \item log.like vector: the The log likelihood for the outcome given the conjugated parameters.
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
  pred.flag <- !is.na(pred_cred)
  like.flag <- !is.null(outcome)
  if (!like.flag & !pred.flag) {
    return(list())
  }

  r <- 1
  a <- conj_param$alpha
  b <- conj_param$beta
  t <- length(a)
  pred <- matrix(NA, r, t)
  var.pred <- array(NA, c(1, 1, t))
  icl.pred <- matrix(NA, r, t)
  icu.pred <- matrix(NA, r, t)
  log.like <- rep(NA, t)
  flags <- b > 1e-20
  if (like.flag) {
    outcome <- matrix(outcome, t, r)
  }

  if (pred.flag) {
    pred[, flags] <- a[flags] / b[flags]
    var.pred[, , flags] <- a[flags] * (b[flags] + 1) / (b[flags]**2)

    icl.pred[, flags] <- qnbinom((1 - pred_cred) / 2, a[flags], (b[flags] / (b[flags] + 1)))
    icu.pred[, flags] <- qnbinom(1 - (1 - pred_cred) / 2, a[flags], (b[flags] / (b[flags] + 1)))
  }
  if (like.flag) {
    log.like[flags] <- dnbinom(outcome[flags, 1], a[flags], (b[flags] / (b[flags] + 1)), log = TRUE)
  }

  N <- 5000
  for (i in (1:t)[!flags]) {
    sample_lambda <- rgamma(N, a[i], b[i])
    sample_y <- rpois(N, sample_lambda)

    if (pred.flag) {
      pred[, i] <- mean(sample_y)
      var.pred[, , i] <- var(sample_y)
      icl.pred[, i] <- quantile(sample_y, (1 - pred_cred) / 2)
      icu.pred[, i] <- quantile(sample_y, 1 - (1 - pred_cred) / 2)
    }
    if (like.flag) {
      log.like.list <- dpois(outcome[i, 1], sample_lambda, log = TRUE)
      max.like.list <- max(like.list)

      log.like[i] <- log(mean(exp(log.like.list - max.like.list))) + max.like.list
    }
  }
  if (!pred.flag) {
    pred <- NULL
    var.pred <- NULL
    icl.pred <- NULL
    icu.pred <- NULL
  }
  if (!like.flag) {
    log.like <- NULL
  }

  return(list(
    "pred" = pred,
    "var.pred" = var.pred,
    "icl.pred" = icl.pred,
    "icu.pred" = icu.pred,
    "log.like" = log.like
  ))
}

#### Alternative Method ####

#' Poisson_alt
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
#' outcome <- Poisson_alt(lambda = "rate", outcome = data)
#'
#' fitted_data <- fit_model(level, season, outcomes = outcome)
#' summary(fitted_data)
#'
#' show_fit(fitted_data, smooth = TRUE)$plot
Poisson_alt <- function(lambda, outcome, offset = outcome**0) {
  distr <- Poisson(lambda, outcome, offset)

  distr$update <- update_Poisson_alt
  distr$alt_method <- TRUE
  return(distr)
}

#' update_Poisson_alt
#'
#' Calculate posterior parameter for the Gamma, assuming that the observed values came from a Poisson model from which the rate parameter (lambda) have prior distribuition Gamma.
#'
#' @param conj_prior list: A vector containing the parameters of the Gamma (alpha,beta). Not used in the alternative method.
#' @param y vector: A vector containing the observations.
#' @param ft vector: A vector representing the means from the normal distribution.
#' @param Qt matrix: A matrix representing the covariance matrix of the normal distribution.
#' @param parms list: A list of extra known parameters of the distribuition. Not used in this kernel.
#'
#' @importFrom cubature cubintegrate
#'
#' @return The parameters of the posterior distribution.
#' @export
update_Poisson_alt <- function(conj_prior, ft, Qt, y, parms) {
  # f0 <- ft
  # Q0 <- Qt
  # S0 <- ginv(Qt)

  # val_const=lgamma(y+1)
  c_val <- -Inf

  f <- function(x) {
    log.prob <- y * log(x) - x + dlnorm(x, ft, sqrt(Qt), log = TRUE)
    max.prob <- max(log.prob)
    if (max.prob > c_val) {
      c_val <- max.prob
    }

    prob <- exp(log.prob - c_val)

    rbind(
      prob,
      log(x) * prob,
      (log(x)**2) * prob
    )
  }

  val <- cubintegrate(f, c(0), c(Inf), fDim = 3, nVec = 1000)$integral
  ft <- matrix(val[2] / val[1], 1, 1)
  Qt <- matrix(val[3], 1, 1) / val[1] - ft**2

  return(list("ft" = ft, "Qt" = Qt))
}
