#### Default Method ####

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
    k <- 1
    distr <- list(
      "conj_prior" = convert_dummy_prior,
      "conj_post" = convert_dummy_post,
      "update" = update_Normal_dummy,
      "smoother" = generic_smoother,
      "calc_pred" = normal_pred,
      "apply_offset" = function(ft, Qt, offset) {
        t <- dim(ft)[2]
        offset <- t(matrix(offset, t, r))
        Qt <- array(Qt, c(r, r, t))

        ft_off <- ft * offset
        Qt_off <- sapply(1:t, function(i) {
          diag(offset[, i]) %*% Qt[, , i] %*% diag(offset[, i])
        }, simplify = "array")
        if (t == 1) {
          Qt_off <- matrix(Qt_off, r, r)
        }

        return(
          list(
            "ft" = ft_off,
            "Qt" = Qt_off
          )
        )
      },
      "link_function" = function(x) {
        x
      },
      "inv_link_function" = function(x) {
        x
      },
      "param_names" = function(y) {
        c(
          paste("ft_", 1:dim(y)[2], sep = ""),
          paste("Qt_",
            c(matrix(1:dim(y)[2], dim(y)[2], dim(y)[2])),
            c(matrix(1:dim(y)[2], dim(y)[2], dim(y)[2], byrow = TRUE)),
            sep = ""
          )
        )
      }
    )
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
    k <- 2
    distr <- list(
      "conj_prior" = convert_NG_Normal,
      "conj_post" = convert_Normal_NG,
      "update" = update_NG,
      "smoother" = generic_smoother,
      "calc_pred" = normal_gamma_pred,
      "apply_offset" = function(ft, Qt, offset) {
        ft[1, ] <- ft[1, ] * offset
        ft[2, ] <- ft[2, ] - log(offset)

        Qt[1, ] <- Qt[1, ] * offset
        Qt[, 1] <- Qt[, 1] * offset
        return(
          list("ft" = ft, "Qt" = Qt)
        )
      },
      "link_function" = function(x) {
        rbind(x[1, ], log(x[2, ]))
      },
      "inv_link_function" = function(x) {
        rbind(x[1, ], exp(x[2, ]))
      },
      "param_names" = function(y) {
        c("mu0", "c0", "alpha", "beta")
      }
    )
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
  distr$var_names <- var_names
  distr$family <- family
  distr$r <- 1
  distr$k <- k
  distr$t <- t
  distr$offset <- matrix(offset, t, r)
  distr$outcome <- matrix(outcome, t, r)
  distr$convert_mat_canom <- convert_mat_canom
  distr$convert_mat_default <- convert_mat_default
  distr$parms <- parms
  distr$name <- "Normal"

  class(distr) <- "dlm_distr"

  distr$alt_method <- FALSE

  return(distr)
}

##### Normal with unknown variance #####

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
  c0 <- exp(-ft[2, ] - Qt[2, 2] / 2) / (Qt[1, 1] + 1e-40)
  helper <- -3 + 3 * sqrt(1 + 2 * Qt[2, 2] / 3)
  # helper=Qt[2,2]
  alpha <- 1 / helper
  beta <- alpha * exp(-ft[2, ] - Qt[2, 2] / 2)
  return(list("mu0" = mu0, "c0" = c0, "alpha" = alpha, "beta" = beta))
}

convert_Normal_NG <- function(conj_distr, parms = list()) {
  # q12 <- conj_distr$mu0*(digamma(conj_distr$alpha)-log(conj_distr$beta))/(1-log(conj_distr$alpha)+log(conj_distr$beta))
  f1 <- conj_distr$mu0 #-q12
  f2 <- digamma(conj_distr$alpha) - log(conj_distr$beta + 1e-40)
  # f2 <- log(conj_distr$alpha)-log(conj_distr$beta)
  q1 <- conj_distr$beta / (conj_distr$c0 * conj_distr$alpha)
  q2 <- trigamma(conj_distr$alpha)
  q12 <- 0

  ft <- c(f1, f2)
  Qt <- matrix(c(q1, q12, q12, q2), byrow = F, ncol = 2)
  return(list("ft" = ft, "Qt" = Qt))
}

#' update_NG
#'
#' Calculate posterior parameter for the Normal-Gamma, assuming that the observed values came from a Normal model from which the prior distribuition for the mean and the precision have joint distribuition Normal-Gamma
#'
#' @param conj_prior list: A vector containing the parameters of the Normal-Gamma (mu0,c0,alpha,beta).
#' @param ft vector: A vector representing the means from the normal distribution. Not used in the default method.
#' @param Qt matrix: A matrix representing the covariance matrix of the normal distribution. Not used in the default method.
#' @param y vector: A vector containing the observations.
#' @param parms list: A list of extra known parameters of the distribuition. Not used in this kernel.
#'
#' @return The parameters of the posterior distribution.
#' @export
update_NG <- function(conj_prior, ft, Qt, y, parms = list()) {
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
#'    \item log.like vector: the The log likelihood for the outcome given the conjugated parameters.
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
  pred.flag <- !is.na(pred_cred)
  like.flag <- !is.null(outcome)
  if (!like.flag & !pred.flag) {
    return(list())
  }

  c0 <- conj_param$c0
  mu0 <- conj_param$mu0
  alpha <- conj_param$alpha
  beta <- conj_param$beta

  mu <- mu0
  nu <- 2 * alpha
  sigma2 <- (beta / alpha) * (1 + 1 / c0)

  if (pred.flag) {
    pred <- mu %>%
      as.matrix() %>%
      t()
    var.pred <- (beta / (alpha - 1)) * (1 + 1 / c0) %>%
      as.matrix() %>%
      t()

    icl.pred <- qt((1 - pred_cred) / 2, nu) * sqrt(sigma2) + mu
    icu.pred <- qt(1 - (1 - pred_cred) / 2, nu) * sqrt(sigma2) + mu
  } else {
    pred <- NULL
    var.pred <- NULL
    icl.pred <- NULL
    icu.pred <- NULL
  }

  if (like.flag) {
    log.like <- dt((outcome - mu) / sqrt(sigma2), nu, log = TRUE)
  } else {
    log.like <- NULL
  }

  list(
    "pred"     = pred,
    "var.pred" = var.pred,
    "icl.pred" = icl.pred,
    "icu.pred" = icu.pred,
    "log.like" = log.like
  )
}

##### Normal with known variance #####

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
#' @param ft vector: A vector representing the means from the normal distribution. Not used in the default method.
#' @param Qt matrix: A matrix representing the covariance matrix of the normal distribution. Not used in the default method.
#' @param y vector: A vector containing the observations.
#' @param parms list: A list of extra known parameters of the distribuition. For this kernel, parms should containg the covariance matrix parameter (Sigma) for the observational Normal model.
#'
#' @return The parameters of the posterior distribution.
#' @export
update_Normal_dummy <- function(conj_prior, ft, Qt, y, parms) {
  Sigma <- parms$Sigma
  if (length(Sigma) == 1) {
    null.flag <- Sigma == 0
  } else {
    null.flag <- all(diag(Sigma) == 0)
  }

  if (null.flag) {
    Qt <- Sigma * 0
    ft <- y
  } else {
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
#'    \item log.like vector: the The log likelihood for the outcome given the conjugated parameters.
#' }
#'
#' @importFrom Rfast dmvnorm
#' @export
#'
#' @examples
#'
#' ft <- c(0.5, -1 / 2, 0, -2.0)
#' Qt <- diag(4)
#'
#' conj_param <- c(ft, Qt)
#'
#' normal_pred(conj_param, parms = list("Sigma" = diag(4)))
normal_pred <- function(conj_param, outcome = NULL, parms = list(), pred_cred = 0.95) {
  pred.flag <- !is.na(pred_cred)
  like.flag <- !is.null(outcome)
  if (!like.flag & !pred.flag) {
    return(list())
  }

  if (is.null(dim(conj_param))) {
    conj_param <- conj_param %>% matrix(1, length(conj_param))
  }
  r <- dim(conj_param)[2]
  r <- (-1 + sqrt(1 + 4 * r)) / 2
  t <- dim(conj_param)[1]

  pred <- t(conj_param[, 1:r] %>% as.matrix())
  var.pred <- conj_param[(r + 1):(r * r + r)] %>%
    t() %>%
    array(c(r, r, t))
  icl.pred <- matrix(NA, t, r)
  icu.pred <- matrix(NA, t, r)
  log.like <- rep(NA, t)
  if (like.flag) {
    outcome <- matrix(outcome, t, r)
  }
  for (i in 1:t) {
    mu <- pred[, i]
    sigma2 <- (var.pred[, , i] + parms$Sigma) %>% matrix(r, r)
    if (pred.flag) {
      var.pred[, , i] <- sigma2
      if (length(sigma2) > 1) {
        sigma2 <- diag(sigma2)
      }
      icl.pred[i, ] <- qnorm((1 - pred_cred) / 2) * sqrt(sigma2) + mu
      icu.pred[i, ] <- qnorm(1 - (1 - pred_cred) / 2) * sqrt(sigma2) + mu
    }
    if (like.flag) {
      log.like[i] <- dmvnorm(outcome[i, ], mu, sigma2, log = TRUE)
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

  list(
    "pred" = pred,
    "var.pred" = var.pred,
    "icl.pred" = icl.pred,
    "icu.pred" = icu.pred,
    "log.like" = log.like
  )
}

#### Alternative Method ####

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
#' outcome <- Normal_alt(mu = "mu", sigma2 = 1, outcome = data)
#'
#' fitted_data <- fit_model(level, outcomes = outcome)
#' summary(fitted_data)
#'
#' show_fit(fitted_data, smooth = TRUE)$plot
#'
#' # Unknown variance
#' outcome <- Normal_alt(mu = "mu", sigma2 = "sigma2", outcome = data)
#'
#' fitted_data <- fit_model(level, variance, outcomes = outcome)
#' summary(fitted_data)
#'
#' show_fit(fitted_data, smooth = TRUE)$plot
Normal_alt <- function(mu, tau = NA, sigma = NA, sigma2 = NA, outcome, offset = outcome**0) {
  distr <- Normal(mu, tau, sigma, sigma2, outcome, offset)

  if (!is.numeric(sigma2)) {
    distr$update <- update_NG_alt
    distr$alt_method <- TRUE
  }
  return(distr)
}

##### Normal with unknown variance #####

#' update_NG_alt
#'
#' Calculate posterior parameter for the Normal-Gamma, assuming that the observed values came from a Normal model from which the prior distribuition for the mean and the precision have joint distribuition Normal-Gamma
#'
#' @param conj_prior list: A vector containing the parameters of the Normal-Gamma (mu0,c0,alpha,beta). Not used in the alternative method.
#' @param ft vector: A vector representing the means from the normal distribution.
#' @param Qt matrix: A matrix representing the covariance matrix of the normal distribution.
#' @param y vector: A vector containing the observations.
#' @param parms list: A list of extra known parameters of the distribuition. Not used in this kernel.
#'
#' @return The parameters of the posterior distribution.
#' @export
update_NG_alt <- function(conj_prior, ft, Qt, y, parms = list()) {
  f0 <- ft
  Q0 <- Qt
  S0 <- ginv(Qt)

  mu1 <- ft[1]
  mu2 <- ft[2]
  sigma1 <- Qt[1, 1]
  sigma2 <- Qt[2, 2]
  rho <- Qt[1, 2] / sqrt(sigma1 * sigma2)

  f <- function(x2) {
    mu1_cond <- mu1 + rho * sqrt(sigma1 / sigma2) * (log(x2) - mu2)
    sigma1_cond <- (1 - rho**2) * sigma1

    prec_obs <- x2
    prec_mu <- 1 / sigma1_cond
    mu1_cond_update <- (prec_obs * y + prec_mu * mu1_cond) / (prec_obs + prec_mu)
    sigma1_cond_update <- 1 / (prec_obs + prec_mu)

    prob <- exp(dnorm(y, mu1_cond, sqrt(sigma1_cond + 1 / x2), log = TRUE) + dlnorm(x2, mu2, sqrt(sigma2), log = TRUE))


    rbind(
      prob,
      mu1_cond_update * prob,
      log(x2) * prob,
      (sigma1_cond_update + mu1_cond_update**2) * prob,
      (mu1_cond_update * log(x2)) * prob,
      (mu1_cond_update * log(x2)) * prob,
      (log(x2)**2) * prob
    )
  }

  val <- cubintegrate(f, c(0), c(Inf), fDim = 7, nVec = 1000)$integral
  ft <- matrix(val[2:3] / val[1], 2, 1)
  Qt <- matrix(val[4:7], 2, 2) / val[1] - ft %*% t(ft)

  return(list("ft" = ft, "Qt" = Qt))
}
