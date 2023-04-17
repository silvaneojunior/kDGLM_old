#### Default Method ####

#' Multinom
#'
#' Creates an outcome with Multinomial distribuition with the chosen parameters.
#'
#' @param p character: a vector with the name of the linear preditor associated with the probality of each category (except the base one, which is assumed to be the last).
#' @param outcome vector: Values of the observed data.
#' @param offset vector: The offset at each observation. Must have the same shape as outcome.
#'
#' @return A object of the class dlm_distr
#' @export
#'
#' @examples
#'
#' # Multinomial case
#' T <- 200
#' y1 <- rpois(T, exp(5 + (-T:T / T) * 5))
#' y2 <- rpois(T, exp(6 + (-T:T / T) * 5 + sin((-T:T) * (2 * pi / 12))))
#' y3 <- rpois(T, exp(5))
#'
#' y <- cbind(y1, y2, y3)
#'
#' level1 <- polynomial_block(p1 = 1)
#' level2 <- polynomial_block(p2 = 1)
#' season <- harmonic_block(p2 = 1, period = 12)
#' outcome <- Multinom(p = c("p1", "p2"), outcome = y)
#'
#' fitted_data <- fit_model(level1, level2, season, outcomes = outcome)
#' summary(fitted_data)
#'
#' show_fit(fitted_data, smooth = TRUE)$plot
Multinom <- function(p, outcome, offset = outcome**0) {
  t <- dim(outcome)[1]
  r <- dim(outcome)[2]
  k <- dim(outcome)[2] - 1
  if (length(p) != k) {
    stop(paste0("Error: Incorrect number of parameter, expected ", k, " got ", length(p), "."))
  }
  convert_mat_default <- convert_mat_canom <- diag(k)
  parms <- list()
  # var_names=p
  # names(var_names)=paste0('p',c(1:r))
  distr <- list(
    var_names = p,
    family = family,
    r = r,
    k = k,
    t = t,
    offset = matrix(offset, t, r),
    outcome = matrix(outcome, t, r),
    convert_mat_canom = convert_mat_canom,
    convert_mat_default = convert_mat_default,
    parms = parms,
    name = "Multinomial",
    conj_prior = convert_Multinom_Normal,
    conj_post = convert_Normal_Multinom,
    update = update_Multinom,
    smoother = generic_smoother,
    calc_pred = multnom_pred,
    apply_offset = function(ft, Qt, offset) {
      t <- dim(ft)[2]
      offset <- matrix(offset, r, t)
      offset_class <- offset[-r, ]
      offset_ref <- matrix(offset[r, ], r - 1, t, byrow = TRUE)
      return(list("ft" = ft + log(offset_class / offset_ref), "Qt" = Qt))
    },
    link_function = function(x) {
      T <- length(x)
      return(log(x[-T, ] / x[T, ]))
    },
    inv_link_function = function(x) {
      y <- exp(x)
      x_last <- 1 / (1 + colSums(y))
      return(rbind(y * x_last, x_last))
    },
    param_names = function(y) {
      paste0("alpha_", 1:dim(y)[2])
    }
  )
  class(distr) <- "dlm_distr"
  distr$alt_method <- FALSE
  return(distr)
}

#' system_multinom
#'
#' Evaluate the compatibilizing equation for the multinomial model (see Ref. Raíra).
#'
#' @param x vector: current tau values.
#' @param parms list: auxiliary values for the system.
#'
#' @return a vector with the values of the system (see Ref. Raíra).
system_multinom <- function(x, parms) {
  x <- exp(x)
  digamma_last <- digamma(x[length(x)] - sum(x[-length(x)]))
  digamma_vec <- digamma(x)

  f_all <- parms$ft - digamma_vec[-length(x)] + digamma_last
  last_guy <- parms$media.log - digamma_last + digamma_vec[length(x)]

  f_all <- c(f_all, last_guy)

  return(f_all)
}

#' convert_Multinom_Normal
#'
#' Calculate the parameters of the Dirichlet that best approximates the given log-Normal distribuition.
#' The approximation is the best in the sense that it minimizes the KL divergence from the log-Normal to the Dirichlet.
#'
#' @param ft vector: A vector representing the means from the normal distribution.
#' @param Qt matrix: A matrix representing the covariance matrix of the normal distribution.
#' @param parms list: A list of extra known parameters of the distribuition. Not used in this kernel.
#'
#' @importFrom rootSolve multiroot
#' @return The parameters of the conjugated distribuition of the linear predictor.
#' @export
convert_Multinom_Normal <- function(ft, Qt, parms = list()) {
  calc_helper <- 1 + sum(exp(ft))
  r <- length(ft)

  H <- exp(ft) %*% t(exp(ft)) / (calc_helper**2)
  diag(H) <- diag(H) - exp(ft) / calc_helper

  media.log <-
    -log(calc_helper) + sum(diag(0.5 * (H %*% Qt)))

  parms <- list("ft" = ft, "media.log" = media.log)

  ss1 <- multiroot(f = system_multinom, start = log(c(rep(0.01, r), 0.01 * (r + 1))), parms = parms, maxiter = 1000, rtol = 10**-20)

  tau <- exp(as.numeric(ss1$root))

  alpha <- tau
  alpha[r + 1] <- tau[r + 1] - sum(tau[-r - 1])
  return(alpha)
}

#' convert_Normal_Multinom
#'
#' Calculate the parameters of the log-Normal that best approximates the given Dirichlet distribuition.
#' The approximation is the best in the sense that it minimizes the KL divergence from the Dirichlet to the log-Normal
#'
#' @param conj_prior list: A vector containing the concentration parameters of the Dirichlet.
#' @param parms list: A list of extra known parameters of the distribuition. Not used in this kernel.
#'
#' @return The parameters of the Normal distribuition of the linear predictor.
#' @export
convert_Normal_Multinom <- function(conj_prior, parms = list()) {
  alpha <- conj_prior
  r <- length(alpha) - 1
  ft <- digamma(alpha[-r - 1]) - digamma(alpha[r + 1])
  Qt <- matrix(trigamma(alpha[r + 1]), r, r)
  diag(Qt) <- trigamma(alpha[-r - 1]) + trigamma(alpha[r + 1])
  return(list("ft" = ft, "Qt" = Qt))
}

#' update_Multinom
#'
#' Calculate posterior parameter for the Dirichlet, assuming that the observed values came from a Multinomial model from which the number of trials is known and the prior distribuition for the probabilities of each category have joint distribuition Dirichlet.
#'
#' @param conj_prior list: A vector containing the concentration parameters of the Dirichlet.
#' @param ft vector: A vector representing the means from the normal distribution. Not used in the default method.
#' @param Qt matrix: A matrix representing the covariance matrix of the normal distribution. Not used in the default method.
#' @param y vector: A vector containing the observations.
#' @param parms list: A list of extra known parameters of the distribuition. Not used in this kernel.
#'
#' @return The parameters of the posterior distribution.
#' @export
update_Multinom <- function(conj_prior, ft, Qt, y, parms = list()) {
  r <- length(y)
  alpha <- conj_prior + y
  return(alpha)
}

#' multnom_pred
#'
#' Calculate the values for the predictive distribuition given the values of the parameter of the conjugated distribuition of the linear predictor.
#' The data is assumed to have Multinomial distribuition with known number of trial N and the probability vector having distribuition Dirichlet with parameters alpha_i.
#' In this scenario, the marginal distribuition of the data is Dirichlet-Multinomial with parameters N and alpha_i.
#'
#' @param conj_param List or data.frame: The parameters of the conjugated distribuitions of the linear predictor.
#' @param outcome Vector or matrix: The observed values at the current time. The value passed is used to compute N.
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
#' conj_param <- matrix(c(1:15), 5, 3, byrow = TRUE)
#'
#' multnom_pred(conj_param, matrix(5, 5, 3))
multnom_pred <- function(conj_param, outcome, parms = list(), pred_cred = 0.95) {
  pred.flag <- !is.na(pred_cred)
  like.flag <- !is.null(outcome)
  if (!like.flag & !pred.flag) {
    return(list())
  }

  if (is.null(dim(conj_param))) {
    conj_param <- conj_param %>%
      data.frame.to_matrix() %>%
      t()
  }

  t <- nrow(conj_param)
  k <- ncol(conj_param) - 1
  r <- ncol(conj_param)

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
    outcome <- matrix(outcome, t, r)
    log.like <- rep(NA, t)
  } else {
    outcome <- matrix(1 / r, t, r)
    log.like <- NULL
  }

  for (t_i in 1:t) {
    outcome_t <- outcome[t_i, ]
    N <- sum(outcome_t)
    N <- max(N, 1)

    alpha <- conj_param[t_i, ] %>% as.numeric()
    alpha0 <- sum(alpha)

    const <- lgamma(alpha0) + lgamma(N + 1) - lgamma(N + alpha0)

    if (pred.flag) {
      p <- alpha / alpha0
      p_var <- p * (1 - p) / (alpha0 + 1)

      pred[, t_i] <- N * p
      var.pred[, , t_i] <- diag(N * p * (1 - p) * (N + alpha0) / (alpha0 + 1))
      for (i in 2:r) {
        for (j in 1:(i - 1)) {
          var.pred[i, j, t_i] <- var.pred[j, i, t_i] <- -N * p[i] * p[j] * (N + alpha0) / (alpha0 + 1)
        }
      }

      x_mat <- matrix(0:N, N + 1, r)
      alpha_mat <- matrix(alpha, N + 1, r, byrow = TRUE)
      x_alpha_mat <- x_mat + alpha_mat

      prob_mat <- lgamma(x_alpha_mat) - lgamma(x_mat + 1) - lgamma(alpha_mat) + lgamma(N + alpha0 - x_alpha_mat) - lgamma(N - x_mat + 1) - lgamma(alpha0 - alpha_mat)
      prob_mat <- exp(const + prob_mat)
      for (i in 1:r) {
        probs_acum <- cumsum(prob_mat[, i])

        icl.pred[i, t_i] <- sum(probs_acum <= ((1 - pred_cred) / 2)) - 1
        icu.pred[i, t_i] <- sum(probs_acum <= (1 - (1 - pred_cred) / 2))

        icl.pred[i, t_i] <- max(0, icl.pred[i])
        icu.pred[i, t_i] <- min(N, icu.pred[i])
      }
    }
    if (like.flag) {
      log.like[t_i] <- const + sum(lgamma(outcome_t + alpha) - lgamma(outcome_t + 1) - lgamma(alpha))
    }
  }

  list(
    "pred"     = pred,
    "var.pred" = var.pred,
    "icl.pred" = icl.pred,
    "icu.pred" = icu.pred,
    "log.like" = log.like
  )
}

#### Alternative Method ####

#' Multinom_alt
#'
#' Creates an outcome with Multinomial distribuition with the chosen parameters.
#'
#' @param p character: a vector with the name of the linear preditor associated with the probality of each category (except the base one, which is assumed to be the last).
#' @param outcome vector: Values of the observed data.
#' @param offset vector: The offset at each observation. Must have the same shape as outcome.
#'
#' @return A object of the class dlm_distr
#' @export
#'
#' @examples
#'
#' # Multinomial case
#' T <- 200
#' y1 <- rpois(T, exp(5 + (-T:T / T) * 5))
#' y2 <- rpois(T, exp(6 + (-T:T / T) * 5 + sin((-T:T) * (2 * pi / 12))))
#' y3 <- rpois(T, exp(5))
#'
#' y <- cbind(y1, y2, y3)
#'
#' level1 <- polynomial_block(p1 = 1)
#' level2 <- polynomial_block(p2 = 1)
#' season <- harmonic_block(p2 = 1, period = 12)
#' outcome <- Multinom_alt(p = c("p1", "p2"), outcome = y)
#'
#' fitted_data <- fit_model(level1, level2, season, outcomes = outcome)
#' summary(fitted_data)
#'
#' show_fit(fitted_data, smooth = TRUE)$plot
Multinom_alt <- function(p, outcome, offset = outcome**0) {
  distr <- Multinom(p, outcome, offset)

  distr$alt_method <- TRUE
  distr$update <- update_Multinom_alt
  return(distr)
}

#' update_Multinom_alt
#'
#' Calculate posterior parameter for the Dirichlet, assuming that the observed values came from a Multinomial model from which the number of trials is known and the prior distribuition for the probabilities of each category have joint distribuition Dirichlet.
#'
#' @param conj_prior list: A vector containing the concentration parameters of the Dirichlet. Not used in the alternative method.
#' @param ft vector: A vector representing the means from the normal distribution.
#' @param Qt matrix: A matrix representing the covariance matrix of the normal distribution.
#' @param y vector: A vector containing the observations.
#' @param parms list: A list of extra known parameters of the distribuition. Not used in this kernel.
#'
#' @return The parameters of the posterior distribution.
#' @export
update_Multinom_alt <- function(conj_prior, ft, Qt, y, parms = list()) {
  f0 <- ft
  S0 <- ginv(Qt)
  n <- sum(y)
  r <- length(y)

  log.like <- function(x) {
    p0 <- c(exp(x), 1)
    p <- p0 / sum(p0)

    sum(y * log(p)) - 0.5 * t(x - f0) %*% S0 %*% (x - f0)
  }

  d1.log.like <- function(x) {
    p0 <- c(x, 0)
    # p0 <- p0-max(p0)
    p0 <- exp(p0)
    p <- p0 / sum(p0)

    y[-r] - n * p[-r] +
      -S0 %*% (x - f0)
  }

  # d2.log.like <- function(x) {
  #   p0 <- c(x, 0)
  #   # p0 <- p0-max(p0)
  #   p0=exp(p0)
  #   p <- p0 / sum(p0)
  #
  #   mat <- -n * p[-r] %*% t(p[-r]) +
  #     -S0
  #   mat
  # }

  d2.log.like <- function(x) {
    mat <- matrix(NA, r - 1, r - 1)
    dev.ref <- d1.log.like(x)
    for (i in 1:(r - 1)) {
      x.dev <- x
      x.dev[i] <- x[i] + 1e-8
      mat[, i] <- (d1.log.like(x.dev) - dev.ref) / (1e-8)
    }
    mat
  }

  mode <- f_root(d1.log.like, d2.log.like, start = f0)$root
  H <- d2.log.like(mode)
  S <- spdinv(-H)
  return(list("ft" = matrix(mode, length(mode), 1), "Qt" = S))
}