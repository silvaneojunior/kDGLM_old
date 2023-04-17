#### Default Method ####

#' Gamma
#'
#' Creates an outcome with gamma distribuition with the chosen parameters (can only specify 2).
#'
#' @param phi character or numeric: Name of the linear preditor associated with the shape parameter of the gamma distribuition. If numeric, this parameter is treated as knowed and equal to the value passed. If a character, the parameter is treated as unknowed and equal to the exponential of the associated linear preditor. It cannot be specified with alpha.
#' @param mu character: Name of the linear preditor associated with the mean parameter of the gamma distribuition. The parameter is treated as unknowed and equal to the exponential of the associated linear preditor.
#' @param alpha character: Name of the linear preditor associated with the shape parameter of the gamma distribuition. The parameter is treated as unknowed and equal to the exponential of the associated linear preditor. It cannot be specified with phi.
#' @param beta character: Name of the linear preditor associated with the rate parameter of the gamma distribuition. The parameter is treated as unknowed and equal to the exponential of the associated linear preditor. It cannot be specified with sigma.
#' @param sigma character: Name of the linear preditor associated with the scale parameter of the gamma distribuition. The parameter is treated as unknowed and equal to the exponential of the associated linear preditor. It cannot be specified with beta.
#' @param outcome vector: Values of the observed data.
#' @param offset vector: The offset at each observation. Must have the same shape as outcome.
#'
#' @return A object of the class dlm_distr
#' @export
#'
#' @examples
#'
#' # Gamma case
#' T <- 200
#' w <- (200 / 40) * 2 * pi
#' phi <- 2.5
#' data <- matrix(rgamma(T, phi, phi / (20 * (sin(w * 1:T / T) + 2))), T, 1)
#'
#' level <- polynomial_block(mu = 1, D = 1 / 0.95)
#' season <- harmonic_block(mu = 1, period = 40, D = 1 / 0.98)
#' scale <- polynomial_block(phi = 1, D = 1 / 1)
#'
#' # Known shape
#' outcome <- Gamma(phi = phi, mu = "mu", outcome = data)
#'
#' fitted_data <- fit_model(level, season, outcomes = outcome)
#' summary(fitted_data)
#'
#' show_fit(fitted_data, smooth = TRUE)$plot
#'
#' # Unknown shape
#' outcome <- Gamma(phi = "phi", mu = "mu", outcome = data)
#'
#' fitted_data <- fit_model(level, season, scale, outcomes = outcome)
#' summary(fitted_data)
#'
#' show_fit(fitted_data, smooth = TRUE)$plot
Gamma <- function(phi = NA, mu = NA, alpha = NA, beta = NA, sigma = NA, outcome, offset = outcome**0) {
  distr <- list(
    "conj_prior" = convert_FGamma_Normal,
    "conj_post" = convert_Normal_FGamma,
    "update" = update_FGamma,
    "smoother" = generic_smoother,
    "calc_pred" = Fgamma_pred,
    "apply_offset" = function(ft, Qt, offset) {
      list("ft" = ft + matrix(c(0, log(offset)), 2, dim(ft)[2]), "Qt" = Qt)
    },
    "link_function" = log,
    "inv_link_function" = exp,
    "param_names" = function(y) {
      c("n", "k", "tau", "theta")
    }
  )
  t <- length(outcome)
  r <- 1

  if (is.numeric(phi)) {
    k <- 1
    if (any(!is.na(c(alpha, beta, sigma)))) {
      stop("Error: When phi is known only mu can be estimated.")
    }
    var_names <- c(mu)
    convert_mat_default <- convert_mat_canom <- diag(1)
    distr <- list(
      "conj_prior" = convert_Gamma_Normal,
      "conj_post" = convert_Normal_Gamma,
      "update" = update_Gamma,
      "smoother" = generic_smoother,
      "calc_pred" = gamma_pred,
      "apply_offset" = function(ft, Qt, offset) {
        t <- if.null(dim(ft)[2], 1)
        offset <- matrix(offset, t, r)

        list("ft" = ft + log(t(offset)), "Qt" = Qt)
      },
      "link_function" = log,
      "inv_link_function" = exp,
      "param_names" = function(y) {
        c("alpha", "beta")
      }
    )

    parms <- list(phi = phi)
  } else {
    k <- 2
    flags <- !is.na(c(phi, mu, alpha, beta, sigma))
    if (sum(flags) < 2) {
      stop("Error: Parameters not fully specified. You must specified exactly 2 non reduntant parameters.")
    }
    if (sum(flags) > 2) {
      stop("Error: Parameters over specified. You must specified exactly 2 non reduntant parameters.")
    }
    if (flags[4] & flags[5]) {
      top("Error: Scale specified in more than one value.")
    }
    convert_mat_default <- matrix(c(1, 0, 0, 1, 1, 0, 1, -1, -2, 2), 2, 5)[, flags]
    convert_mat_canom <- solve(convert_mat_default)
    parms <- list()
    var_names <- c(phi, mu, alpha, beta, sigma)[flags]
  }

  distr$var_names <- var_names
  distr$r <- r
  distr$k <- k
  distr$t <- t
  distr$offset <- matrix(offset, t, r)
  distr$outcome <- matrix(outcome, t, r)
  distr$convert_mat_canom <- convert_mat_canom
  distr$convert_mat_default <- convert_mat_default
  distr$parms <- parms
  distr$name <- "Gamma"

  class(distr) <- "dlm_distr"
  distr$alt_method <- FALSE
  return(distr)
}

##### Gamma with unknown shape and mean #####

#' system_full_gamma
#'
#' Evaluate the compatibilizing equation for the full gamma model (see Ref. Raíra).
#'
#' @param x vector: current tau values.
#' @param parms list: auxiliary values for the system.
#'
#' @importFrom cubature cubintegrate
#'
#' @return a vector with the values of the system (see Ref. Raíra).
system_full_gamma <- function(x, parms) {
  n <- exp(x) # exp(x[1])
  k <- n # exp(x[2])
  tau <- (n * parms$Hq1 + 1) / parms$Hq2
  theta <- n - k + n * log(tau / n) - (k + 1) / (2 * parms$Hq1)

  a <- (k + 1) / 2
  b <- (n - k + n * log(tau / n) - theta)

  # print((parms$Hq3 + parms$Hq4))
  if (a <= 5) {
    # Densidade marginal aproximada de alpha (uso opcional).
    # f_densi=function(x){dgamma(a,b))}
    # c_val=1
    # Densidade marginal exata de phi.
    f_densi_raw <- function(x) {
      exp(k * (x + 1) * log(x) + lgamma(n * x + 1) + theta * x - k * lgamma(x + 1) - (n * x + 1) * log(x * tau))
    }
    lim_sup <- Inf
    c_val <- cubintegrate(f_densi_raw, 0, lim_sup, nVec = 200)$integral
    f_densi <- function(x) {
      f_densi_raw(x) / c_val
    }
    # print('a')
    f <- function(x) {
      (x * digamma(x * n + 1) - x * log(x) - x * log(tau)) * f_densi(x)
    }
    Hp3 <- cubintegrate(f, 0, lim_sup, nVec = 200)$integral

    # print('b')
    f <- function(x) {
      (x * log(x) - lgamma(x)) * f_densi(x)
    }
    Hp4 <- cubintegrate(f, 0, lim_sup, nVec = 200)$integral

    # print('c')
    # f <- function(x) {
    #   (x * digamma(x * n + 1)  - lgamma(x)- x*log(tau)) * f_densi(x)
    # }
    # Hp5 <- cubintegrate(f, 0, Inf, nVec = 200)$integral
    # print('sd')
    Hp5 <- Hp3 + Hp4


    # f <- function(x) {
    #   (lgamma(x)-x*digamma(n*x+1)) * f_densi(x)
    # }
    # Hp5 <- integrate(f, 0, Inf)$value+parms$Hq1*log(tau)
  } else {
    c_val <- 1
    # Hp3 <- log(tau / n) * a / b - 1 / n + b / (12 * (n**2) * (a - 1))
    # Hp4 <- a / b + 0.5 * (digamma(a) - log(b)) - b / (12 * (a - 1)) - 11 / 12
    # Hp5=Hp3+Hp4
    Hp5 <- parms$Hq1 * (log(tau / n) - 1) - 0.5 * digamma((n + 1) / 2) + 0.5 * log(n * log(tau / n) - theta) + (log(tau / n) - theta / n) / 6 + 11 / 12 - 1 / (2 * n)
  }

  f_all <- c(
    (Hp5 - (parms$Hq3 + parms$Hq4))
  )
  # print(f_all)
  return(f_all)
}

#' convert_FGamma_Normal
#'
#' DESCRIPTION
#'
#' @param ft vector: A vector representing the means from the normal distribution.
#' @param Qt matrix: A matrix representing the covariance matrix of the normal distribution.
#' @param parms list: A list of extra known parameters of the distribuition. Not used in this function.
#'
#' @return The parameters of the conjugated distribuition of the linear predictor.
#' @export
convert_FGamma_Normal <- function(ft, Qt, parms) {
  # print(ft)
  # print(Qt)


  # ft=matrix(c(-1.457271, 3.533958),2,1)
  # Qt=matrix(c(0.17878138,0.03349385,0.03349385,0.15984907),2,2)

  s <- exp(ft[2, ] + -1)
  f1 <- ft[1, ]
  f2 <- ft[2, ] - log(s)
  q1 <- Qt[1, 1]
  q2 <- Qt[2, 2]
  q12 <- Qt[1, 2]

  Hq1 <- exp(f1 + q1 / 2)
  Hq2 <- exp(f1 - f2 + (q1 + q2 - 2 * q12) / 2)
  Hq3 <- -(f2 + q12) * Hq1

  Hq4 <- cubintegrate(function(x) {
    (x * log(x) - lgamma(x)) * dlnorm(x, f1, sqrt(q1))
  }, 0, Inf, nVec = 200)$integral

  parms <- list(
    "Hq1" = Hq1,
    "Hq2" = Hq2,
    "Hq3" = Hq3,
    "Hq4" = Hq4
  )

  ss1 <- multiroot(f = system_full_gamma, start = c(0), parms = parms, maxiter = 2000, atol = 10**-20)
  x <- as.numeric(ss1$root)
  n <- exp(x) # exp(x[1])
  # print(q1)
  # if(q1<1.5){
  #   h=-1+sqrt(1+2*q1)
  #   n=2*(1/h)-1
  # }else{
  #   n=1
  # }

  k <- n # exp(x[2])

  # Calculando tau e theta dado n e k
  tau <- ((n * parms$Hq1 + 1) / parms$Hq2)
  # tau=exp(x[3])
  theta <- (n - k + n * log(tau / n) - (k + 1) / (2 * parms$Hq1))
  tau <- tau * s
  theta <- theta + n * log(s)
  return(list("n" = n, "k" = k, "tau" = tau, "theta" = theta))
}

#' convert_Normal_FGamma
#'
#' Calculate the parameters of the log-Normal that best approximates the given Inverse-Gamma distribuition.
#' The approximation is the best in the sense that it minimizes the KL divergence from the Inverse-Gamma to the log-Normal
#'
#' @param conj_prior list: A vector containing the parameters of the Inverse-Gamma (alpha,beta).
#' @param parms list: A list of extra known parameters of the distribuition. Not used in this function.
#'
#' @importFrom cubature cubintegrate
#'
#' @return The parameters of the Normal distribuition of the linear predictor.
#' @export
convert_Normal_FGamma <- function(conj_prior, parms) {
  n <- conj_prior$n
  tau <- conj_prior$tau
  theta <- conj_prior$theta
  k <- conj_prior$k

  s <- (10 * tau / n)

  tau <- tau / s
  theta <- theta - n * log(s)

  # Parâmetros da densidade aproximada de alpha
  a <- (k + 1) / 2
  b <- (n - k + n * log(tau / n) - theta)

  # Comentar essas linhas caso a densidade aproximada seja usada.
  f_densi <- function(x) {
    exp(k * (x + 1) * log(x) + lgamma(n * x + 1) + theta * x - k * lgamma(x + 1) - (n * x + 1) * log(x * tau))
  }
  c_val <- cubintegrate(f_densi, 0, Inf, nVec = 200)$integral

  # Média 1 calculada com a densidade exata.
  f <- function(x) {
    log(x) * f_densi(x)
  }
  f1 <- cubintegrate(f, 0, Inf, nVec = 200)$integral / c_val
  # Média 1 calculada com a densidade aproximada para phi.
  # f1=digamma(a)-log(b)

  # Média 2 calculada com a densidade exata.
  # Lembremos que mu|phi ~ IG(n*phi+1,phi*tau), logo E[log(mu)]=E[E[log(mu)|phi]]=E[digamma(n*phi+1)-log(tau*phi)]
  f <- function(x) {
    (-digamma(n * x + 1) + log(x * tau)) * f_densi(x)
  }
  f2 <- cubintegrate(f, 0, Inf, nVec = 200)$integral / c_val
  # Média 2 calculada com a densidade aproximada
  # f2=log(tau/n)-(1/(2*n))*(b/(a-1))+(1/(12*n**2))*(b**2)/((a-1)*(a-2))

  # Variância 1 calculada com a densidade exata.
  f <- function(x) {
    (log(x)**2) * f_densi(x)
  }
  Q1 <- cubintegrate(f, 0, Inf, nVec = 200)$integral / c_val - f1**2
  # Variância 1 calculada com a densidade aproximada.
  # Q1=trigamma(a)

  # Variância 2 calculada com a densidade exata.
  # O mesmo argumento para a média foi usado para o segundo momento.
  f <- function(x) {
    ((-digamma(n * x + 1) + log(x * tau))**2) * f_densi(x)
  }
  Q2 <- cubintegrate(f, 0, Inf, nVec = 200)$integral / c_val - f2**2
  # Variância 2 calculada com a densidade aproximada.

  # Covariância calculada com a densidade exata.
  # O mesmo argumento para a média e para o segundo momento foi usado para a covariância.
  f <- function(x) {
    (log(x) - f1) * (-digamma(n * x + 1) + log(x * tau) - f2) * f_densi(x)
  }
  Q12 <- cubintegrate(f, 0, Inf, nVec = 200)$integral / c_val
  # Covariância calculada com a densidade aproximada

  ft <- matrix(c(f1, f2 + log(s)), 2, 1)
  Qt <- matrix(c(Q1, Q12, Q12, Q2), 2, 2)
  # print("c")
  return(list("ft" = ft, "Qt" = Qt))
}

#' update_FGamma
#'
#' Calculate posterior parameter for the Inverser-Gamma, assuming that the observed values came from a Gamma model from which the shape parameter (phi) is known and the mean (mu) have prior distribuition Inverse-Gamma.
#'
#' @param conj_prior list: A vector containing the parameters of the Inverse-Gamma (alpha,beta).
#' @param ft vector: A vector representing the means from the normal distribution. Not used in the default method.
#' @param Qt matrix: A matrix representing the covariance matrix of the normal distribution. Not used in the default method.
#' @param y vector: A vector containing the observations.
#' @param parms list: A list of extra known parameters of the distribuition. For this kernel, parms should containg the shape parameter (phi) for the observational gamma model.
#'
#' @return The parameters of the posterior distribution.
#' @export
update_FGamma <- function(conj_prior, ft, Qt, y, parms) {
  n0 <- conj_prior$n
  k0 <- conj_prior$k
  tau0 <- conj_prior$tau
  theta0 <- conj_prior$theta

  a <- (k0 + 1) / 2
  b <- (n0 - k0 + n0 * log(tau0 / n0) - theta0)

  n1 <- n0 + 1
  k1 <- k0 + 1
  tau1 <- tau0 + y
  theta1 <- theta0 + log(y)

  a <- (k1 + 1) / 2
  b <- (n1 - k1 + n1 * log(tau1 / n1) - theta1)

  return(list("n" = n1, "k" = k1, "tau" = tau1, "theta" = theta1))
}

#' Fgamma_pred
#'
#' Calculate the values for the predictive distribuition given the values of the parameter of the distribuition of the linear predictor.
#' The data is assumed to have Gamma distribuition with unknown shape phi and unknown mean having log-Normal distribuition.
#' In this scenario, the marginal distribuition of the data is obtained via Monte Carlo.
#'
#' @param conj_param List or data.frame: The parameters of the distribuition of the linear predictor.
#' @param outcome Vector or matrix (optional): The observed values at the current time. Not used in this function.
#' @param parms List: A list of extra parameters for the model. For this function, it must contain the shape parameter phi of the observational model.
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
#' params <- data.frame(
#'   "f1" = c(1:3),
#'   "f2" = c(3:1),
#'   "Q11" = rep(1, 3),
#'   "Q12" = rep(0, 3),
#'   "Q21" = rep(0, 3),
#'   "Q22" = rep(1, 3),
#' )
#'
#' Fgamma_pred(params)
Fgamma_pred <- function(conj_param, outcome = NULL, parms = list(), pred_cred = 0.95) {
  pred.flag <- !is.na(pred_cred)
  like.flag <- !is.null(outcome)
  if (!like.flag & !pred.flag) {
    return(list())
  }
  r <- 1
  t <- length(conj_param$k)
  k <- 2

  if (pred.flag) {
    pred <- matrix(NA, r, t)
    var.pred <- array(NA, c(r, r, t))
    icl.pred <- matrix(NA, r, t)
    icu.pred <- matrix(NA, r, t)
  } else {
    pred <- NULL
    var.pred <- NULL
    icl.pred <- NULL
    icu.pred <- NULL
  }
  if (like.flag) {
    outcome <- matrix(outcome, r, t)
    log.like <- rep(NA, t)
  } else {
    log.like <- NULL
  }
  N <- 5000

  for (i in 1:t) {
    n <- conj_param$n[i]
    k <- conj_param$k[i]
    tau <- conj_param$tau[i]
    theta <- conj_param$theta[i]


    a <- (k + 1) / 2
    b <- (n - k + n * log(tau / n) - theta)

    alpha_i <- rgamma(N, a, b)
    mu_i <- 1 / rgamma(N, n * alpha_i + 1, alpha_i * tau)

    sample_y <- rgamma(N, alpha_i, alpha_i / mu_i)
    if (pred.flag) {
      pred[, i] <- mean(sample_y)
      var.pred[, , i] <- var(sample_y)
      icl.pred[, i] <- quantile(sample_y, (1 - pred_cred) / 2)
      icu.pred[, i] <- quantile(sample_y, 1 - (1 - pred_cred) / 2)
    }
    if (like.flag) {
      log.like.list <- dgamma(outcome[i, ], alpha_i, alpha_i / mu_i, log = TRUE)
      max.log.like <- max(log.like.list)
      like.list <- exp(log.like.list - max.log.like)
      log.like[i] <- log(mean(like.list)) + max.log.like
    }
  }

  outcome_list <- list(
    "pred"     = pred,
    "var.pred" = var.pred,
    "icl.pred" = icl.pred,
    "icu.pred" = icu.pred,
    "log.like" = log.like
  )
  return(outcome_list)
}

##### Gamma with known shape but unknown mean #####

#' convert_Gamma_Normal
#'
#' Calculate the parameters of the Inverse-Gamma that best approximates the given log-Normal distribuition.
#' The approximation is the best in the sense that it minimizes the KL divergence from the log-Normal to the Inverse-Gamma
#'
#' @param ft vector: A vector representing the means from the normal distribution.
#' @param Qt matrix: A matrix representing the covariance matrix of the normal distribution.
#' @param parms list: A list of extra known parameters of the distribuition. Not used in this function.
#'
#' @return The parameters of the conjugated distribuition of the linear predictor.
#' @export
convert_Gamma_Normal <- function(ft, Qt, parms) {
  alpha <- 1 / (-3 + 3 * sqrt(1 + 2 * Qt / 3))
  beta <- alpha * exp(ft - Qt / 2)
  return(list("alpha" = alpha, "beta" = beta))
}

#' convert_Normal_Gamma
#'
#' Calculate the parameters of the log-Normal that best approximates the given Inverse-Gamma distribuition.
#' The approximation is the best in the sense that it minimizes the KL divergence from the Inverse-Gamma to the log-Normal
#'
#' @param conj_prior list: A vector containing the parameters of the Inverse-Gamma (alpha,beta).
#' @param parms list: A list of extra known parameters of the distribuition. Not used in this function.
#'
#' @return The parameters of the Normal distribuition of the linear predictor.
#' @export
convert_Normal_Gamma <- function(conj_prior, parms) {
  ft <- -digamma(conj_prior$alpha) + log(conj_prior$beta)
  Qt <- trigamma(conj_prior$alpha)
  return(list("ft" = ft, "Qt" = Qt))
}

#' update_Gamma
#'
#' Calculate posterior parameter for the Inverse-Gamma, assuming that the observed values came from a Gamma model from which the shape parameter (phi) is known and the mean (mu) have prior distribuition Inverse-Gamma.
#'
#' @param conj_prior list: A vector containing the parameters of the Inverse-Gamma (alpha,beta).
#' @param ft vector: A vector representing the means from the normal distribution. Not used in the default method.
#' @param Qt matrix: A matrix representing the covariance matrix of the normal distribution. Not used in the default method.
#' @param y vector: A vector containing the observations.
#' @param parms list: A list of extra known parameters of the distribution. For this kernel, parms should containg the shape parameter (phi) for the observational gamma model.
#'
#' @return The parameters of the posterior distribution.
#' @export
update_Gamma <- function(conj_prior, ft, Qt, y, parms) {
  alpha <- conj_prior$alpha + parms$phi
  beta <- conj_prior$beta + y * parms$phi
  return(list("alpha" = alpha, "beta" = beta))
}

#' gamma_pred
#'
#' Calculate the values for the predictive distribuition given the values of the parameter of the conjugated distribuition of the linear predictor.
#' The data is assumed to have Gamma distribuition with known shape parameter phi and it's mean having distribuition Inverse-Gamma with shape parameter a e rate parameter b.
#' In this scenario, the marginal distribuition of the data is Beta prime with parameters phi, alpha, beta / phi.
#'
#' @param conj_param List or data.frame: The parameters of the conjugated distribuitions of the linear predictor.
#' @param outcome Vector or matrix (optional): The observed values at the current time. Not used in this function.
#' @param parms List: A list of extra parameters for the model. For this function, it must contain the shape parameter phi of the observational model.
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
#' @importFrom extraDistr qbetapr dbetapr
#' @export
#'
#' @examples
#'
#' conj_param <- data.frame(
#'   "alpha" = c(1:3),
#'   "beta" = c(3:1)
#' )
#'
#' gamma_pred(conj_param, parms = list("phi" = 1))
gamma_pred <- function(conj_param, outcome = NULL, parms = list(), pred_cred = 0.95) {
  pred.flag <- !is.na(pred_cred)
  like.flag <- !is.null(outcome)
  if (!like.flag & !pred.flag) {
    return(list())
  }

  phi <- parms$phi

  alpha <- (phi * conj_param$alpha) %>% t()
  beta <- (phi * conj_param$beta) %>% t()
  pred <- NULL
  var.pred <- NULL
  icl.pred <- NULL
  icu.pred <- NULL
  log.like <- NULL

  if (pred.flag) {
    pred <- beta / (alpha - 1)
    var.pred <- ((alpha / (beta - 1))**2) * (alpha + phi - 1) / ((alpha - 2) * phi)

    icl.pred <- qbetapr((1 - pred_cred) / 2, phi, alpha, beta / phi)
    icu.pred <- qbetapr(1 - (1 - pred_cred) / 2, phi, alpha, beta / phi)
  }
  if (like.flag) {
    log.like <- dbetapr(outcome, phi, alpha, beta / phi, log = TRUE)
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

#' Gamma_alt
#'
#' Creates an outcome with gamma distribuition with the chosen parameters (can only specify 2).
#'
#' @param phi character or numeric: Name of the linear preditor associated with the shape parameter of the gamma distribuition. If numeric, this parameter is treated as knowed and equal to the value passed. If a character, the parameter is treated as unknowed and equal to the exponential of the associated linear preditor. It cannot be specified with alpha.
#' @param mu character: Name of the linear preditor associated with the mean parameter of the gamma distribuition. The parameter is treated as unknowed and equal to the exponential of the associated linear preditor.
#' @param alpha character: Name of the linear preditor associated with the shape parameter of the gamma distribuition. The parameter is treated as unknowed and equal to the exponential of the associated linear preditor. It cannot be specified with phi.
#' @param beta character: Name of the linear preditor associated with the rate parameter of the gamma distribuition. The parameter is treated as unknowed and equal to the exponential of the associated linear preditor. It cannot be specified with sigma.
#' @param sigma character: Name of the linear preditor associated with the scale parameter of the gamma distribuition. The parameter is treated as unknowed and equal to the exponential of the associated linear preditor. It cannot be specified with beta.
#' @param outcome vector: Values of the observed data.
#' @param offset vector: The offset at each observation. Must have the same shape as outcome.
#'
#' @return A object of the class dlm_distr
#' @export
#'
#' @examples
#'
#' # Gamma case
#' T <- 200
#' w <- (200 / 40) * 2 * pi
#' phi <- 2.5
#' data <- matrix(rgamma(T, phi, phi / (20 * (sin(w * 1:T / T) + 2))), T, 1)
#'
#' level <- polynomial_block(mu = 1, D = 1 / 0.95)
#' season <- harmonic_block(mu = 1, period = 40, D = 1 / 0.98)
#' scale <- polynomial_block(phi = 1, D = 1 / 1)
#'
#' # Known shape
#' outcome <- Gamma_alt(phi = phi, mu = "mu", outcome = data)
#'
#' fitted_data <- fit_model(level, season, outcomes = outcome)
#' summary(fitted_data)
#'
#' show_fit(fitted_data, smooth = TRUE)$plot
#'
#' # Unknown shape
#' outcome <- Gamma_alt(phi = "phi", mu = "mu", outcome = data)
#'
#' fitted_data <- fit_model(level, season, scale, outcomes = outcome)
#' summary(fitted_data)
#'
#' show_fit(fitted_data, smooth = TRUE)$plot
Gamma_alt <- function(phi = NA, mu = NA, alpha = NA, beta = NA, sigma = NA, outcome, offset = outcome**0) {
  distr <- Gamma(phi, mu, alpha, beta, sigma, outcome, offset)

  if (is.numeric(phi)) {
    distr$update <- update_Gamma_alt
  } else {
    distr$conj_prior <- format_ft
    distr$conj_post <- format_param
    distr$update <- update_FGamma_alt
    distr$calc_pred <- Fgamma_pred_alt
    distr$param_names <- function(y) {
      c("f1", "f2", "Q11", "Q12", "21", "Q22")
    }
  }
  distr$alt_method <- TRUE
  return(distr)
}

##### Gamma with unknown shape and mean #####

#' update_FGamma_alt
#'
#' Calculate posterior parameter for the Inverser-Gamma, assuming that the observed values came from a Gamma model from which the shape parameter (phi) is known and the mean (mu) have prior distribuition Inverse-Gamma.
#'
#' @param conj_prior list: A vector containing the parameters of the Inverse-Gamma (alpha,beta). Not used in the alternative method.
#' @param ft vector: A vector representing the means from the normal distribution.
#' @param Qt matrix: A matrix representing the covariance matrix of the normal distribution.
#' @param y vector: A vector containing the observations.
#' @param parms list: A list of extra known parameters of the distribuition. For this kernel, parms should containg the shape parameter (phi) for the observational gamma model.
#'
#' @importFrom MASS ginv
#'
#' @return The parameters of the posterior distribution.
#' @export
update_FGamma_alt <- function(conj_prior, ft, Qt, y, parms) {
  # y=1
  # ft=matrix(0,2,1)
  # Qt=diag(2)*1

  f0 <- ft
  S0 <- ginv(Qt)

  # log.like=function(x){
  #   phi=exp(x[1])
  #   mu=exp(x[2])
  #
  #   phi*log(phi)-phi*log(mu)-lgamma(phi)+(phi-1)*log(y)-phi*y/mu-0.5*t(x-f0)%*%S0%*%(x-f0)
  # }

  d1.log.like <- function(x) {
    phi <- exp(x[1])
    mu <- exp(x[2])

    c(
      (log(phi) + 1 - log(mu) - digamma(phi) + log(y) - y / mu) * phi,
      (-phi / mu + phi * y / (mu**2)) * mu
    ) - S0 %*% (x - f0)
  }

  d2.inv <- function(x) {
    phi <- exp(x[1])
    mu <- exp(x[2])

    cross <- (-1 + y / mu) * phi

    mat <- matrix(
      c(
        (log(phi) + 1 - log(mu) - digamma(phi) + log(y) - y / mu + 1 - phi * trigamma(phi)) * phi,
        cross,
        cross,
        (-1 + y / mu + 1 - 2 * y / mu) * phi
      ),
      2, 2
    ) - S0

    # mat_inv=mat
    # mat_inv[1,1]=mat[2,2]
    # mat_inv[2,2]=mat[1,1]
    # mat_inv[1,2]=mat_inv[2,1]=-mat[1,2]
    # det_mat=mat[1,1]*mat[2,2]+mat[1,2]**2
    # mat_inv=mat_inv/det_mat

    return(mat)
  }

  mode <- f_root(d1.log.like, d2.inv, f0)$root
  H <- d2.inv(mode)
  S <- spdinv(-H)

  return(list("ft" = matrix(mode, 2, 1), "Qt" = S))
}

#' Fgamma_pred_alt
#'
#' Calculate the values for the predictive distribuition given the values of the parameter of the distribuition of the linear predictor.
#' The data is assumed to have Gamma distribuition with unknown shape phi and unknown mean having log-Normal distribuition.
#' In this scenario, the marginal distribuition of the data is obtained via Monte Carlo.
#'
#' @param conj_param List or data.frame: The parameters of the distribuition of the linear predictor.
#' @param outcome Vector or matrix (optional): The observed values at the current time. Not used in this function.
#' @param parms List: A list of extra parameters for the model. For this function, it must contain the shape parameter phi of the observational model.
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
#' @importFrom Rfast cholesky
#' @export
#'
#' @examples
#'
#' params <- data.frame(
#'   "f1" = c(1:3),
#'   "f2" = c(3:1),
#'   "Q11" = rep(1, 3),
#'   "Q12" = rep(0, 3),
#'   "Q21" = rep(0, 3),
#'   "Q22" = rep(1, 3),
#' )
#'
#' Fgamma_pred_alt(params)
Fgamma_pred_alt <- function(conj_param, outcome = NULL, parms = list(), pred_cred = 0.95) {
  pred.flag <- !is.na(pred_cred)
  like.flag <- !is.null(outcome)

  norm_param <- format_param(conj_param)
  ft <- norm_param$ft
  Qt <- norm_param$Qt

  k <- 2
  t <- dim(ft)[2]
  r <- 1

  if (pred.flag) {
    pred <- matrix(NA, r, t)
    var.pred <- array(NA, c(r, r, t))
    icl.pred <- matrix(NA, r, t)
    icu.pred <- matrix(NA, r, t)
  } else {
    pred <- NULL
    var.pred <- NULL
    icl.pred <- NULL
    icu.pred <- NULL
  }
  if (like.flag) {
    outcome <- matrix(outcome, r, t)
    log.like <- rep(NA, t)
  } else {
    log.like <- NULL
  }

  N <- 5000
  sample <- matrix(rnorm(k * N), N, k)
  for (i in 1:t) {
    ft_i <- sample %*% cholesky(Qt[, , i]) + matrix(ft[, i], N, k, byrow = TRUE)
    sample_y <- rgamma(N, exp(ft_i[, 1]), exp(ft_i[, 1] - ft_i[, 2]))
    if (pred.flag) {
      pred[, i] <- mean(sample_y)
      var.pred[, , i] <- var(sample_y)
      icl.pred[, i] <- quantile(sample_y, (1 - pred_cred) / 2)
      icu.pred[, i] <- quantile(sample_y, 1 - (1 - pred_cred) / 2)
    }
    if (like.flag) {
      log.like.list <- dgamma(outcome[, i], exp(ft_i[1, ]), exp(ft_i[1, ] - ft_i[2, ]), log = TRUE)
      max.log.like <- max(log.like.list)
      like.list <- exp(log.like.list - max.log.like)
      log.like[i] <- log(mean(like.list)) + max.log.like
    }
  }

  outcome_list <- list(
    "pred"     = pred,
    "var.pred" = var.pred,
    "icl.pred" = icl.pred,
    "icu.pred" = icu.pred,
    "log.like" = log.like
  )
  return(outcome_list)
}

##### Gamma with known shape but unknown mean #####

#' update_Gamma_alt
#'
#' Calculate posterior parameter for the Inverse-Gamma, assuming that the observed values came from a Gamma model from which the shape parameter (phi) is known and the mean (mu) have prior distribuition Inverse-Gamma.
#'
#' @param conj_prior list: A vector containing the parameters of the Inverse-Gamma (alpha,beta). Not used in the alternative method.
#' @param ft vector: A vector representing the means from the normal distribution.
#' @param Qt matrix: A matrix representing the covariance matrix of the normal distribution.
#' @param y vector: A vector containing the observations.
#' @param parms list: A list of extra known parameters of the distribution. For this kernel, parms should containg the shape parameter (phi) for the observational gamma model.
#'
#' @importFrom MASS ginv
#' @importFrom cubature cubintegrate
#'
#' @return The parameters of the posterior distribution.
#' @export
update_Gamma_alt <- function(conj_prior, ft, Qt, y, parms) {
  f0 <- ft
  Q0 <- Qt
  S0 <- ginv(Qt)

  f <- function(x) {
    prob <- exp(dgamma(y, parms$phi, parms$phi / x, log = TRUE) + dlnorm(x, ft, sqrt(Qt), log = TRUE))

    rbind(
      prob,
      log(x) * prob,
      (log(x)**2) * prob
    )
  }

  val <- cubintegrate(f, c(0), c(Inf), fDim = 3, nVec = 1000)$integral
  ft <- matrix(val[2] / val[1], 1, 1)
  Qt <- matrix(val[3], 1, 1) / val[1] - ft %*% t(ft)

  return(list("ft" = ft, "Qt" = Qt))
}