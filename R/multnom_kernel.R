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

#' convert_Dir_Normal
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
convert_Dir_Normal <- function(ft, Qt, parms = list()) {
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

#' convert_Normal_Dir
#'
#' Calculate the parameters of the log-Normal that best approximates the given Dirichlet distribuition.
#' The approximation is the best in the sense that it minimizes the KL divergence from the Dirichlet to the log-Normal
#'
#' @param conj_prior list: A vector containing the concentration parameters of the Dirichlet.
#' @param parms list: A list of extra known parameters of the distribuition. Not used in this kernel.
#'
#' @return The parameters of the Normal distribuition of the linear predictor.
#' @export
convert_Normal_Dir <- function(conj_prior, parms = list()) {
  alpha <- conj_prior
  r <- length(alpha) - 1
  ft <- digamma(alpha[-r - 1]) - digamma(alpha[r + 1])
  Qt <- matrix(trigamma(alpha[r + 1]), r, r)
  diag(Qt) <- trigamma(alpha[-r - 1]) + trigamma(alpha[r + 1])
  return(list("ft" = ft, "Qt" = Qt))
}

#' update_Dir
#'
#' Calculate posterior parameter for the Dirichlet, assuming that the observed values came from a Multinomial model from which the number of trials is known and the prior distribuition for the probabilities of each category have joint distribuition Dirichlet.
#'
#' @param conj_prior list: A vector containing the concentration parameters of the Dirichlet.
#' @param y vector: A vector containing the observations.
#' @param parms list: A list of extra known parameters of the distribuition. Not used in this kernel.
#'
#' @return The parameters of the posterior distribution.
#' @export
update_Dir <- function(conj_prior, y, parms = list()) {
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
  conj_param <- if (is.null(dim(conj_param))) {
    conj_param %>% t()
  } else {
    conj_param
  }
  outcome <- if (is.null(dim(outcome))) {
    outcome %>% t()
  } else {
    outcome
  }
  T <- nrow(conj_param)
  r <- ncol(conj_param) - 1

  pred <- matrix(NA, r + 1, T)
  var.pred <- array(NA, c(r + 1, r + 1, T))
  icl.pred <- matrix(NA, r + 1, T)
  icu.pred <- matrix(NA, r + 1, T)
  for (t in 1:T) {
    outcome_t <- outcome[t, ]
    N <- sum(outcome_t)
    N <- max(N, 1)

    alpha <- conj_param[t, ] %>% as.numeric()
    alpha0 <- sum(alpha)

    p <- alpha / alpha0
    p_var <- p * (1 - p) / (alpha0 + 1)

    pred[, t] <- N * p
    var.pred[, , t] <- diag(N * p * (1 - p) * (N + alpha0) / (alpha0 + 1))
    for (i in 2:(r + 1)) {
      for (j in 1:(i - 1)) {
        var.pred[i, j, t] <- var.pred[j, i, t] <- -N * p[i] * p[j] * (N + alpha0) / (alpha0 + 1)
      }
    }

    const <- lgamma(alpha0) + lgamma(N + 1) - lgamma(N + alpha0)

    x_mat <- matrix(0:N, N + 1, r + 1)
    alpha_mat <- matrix(alpha, N + 1, r + 1, byrow = TRUE)
    x_alpha_mat <- x_mat + alpha_mat

    prob_mat <- lgamma(x_alpha_mat) - lgamma(x_mat + 1) - lgamma(alpha_mat) + lgamma(N + alpha0 - x_alpha_mat) - lgamma(N - x_mat + 1) - lgamma(alpha0 - alpha_mat)
    prob_mat <- exp(const + prob_mat)
    for (i in 1:(r + 1)) {
      probs_acum <- cumsum(prob_mat[, i])

      icl.pred[i, t] <- sum(probs_acum <= ((1 - pred_cred) / 2)) - 1
      icu.pred[i, t] <- sum(probs_acum <= (1 - (1 - pred_cred) / 2))

      icl.pred[i, t] <- max(0, icl.pred[i])
      icu.pred[i, t] <- min(N, icu.pred[i])
    }
  }

  list(
    "pred"     = pred,
    "var.pred" = var.pred,
    "icl.pred" = icl.pred,
    "icu.pred" = icu.pred
  )
}

#' multnom_log_like
#'
#' Calculate the predictive log-likelyhood of the data given the values of the parameter of the conjugated distribuition of the linear predictor.
#' The data is assumed to have Multinomial distribuition with known number of trial N and the probability vector having distribuition Dirichlet with parameters alpha_i.
#' In this scenario, the marginal distribuition of the data is Dirichlet-Multinomial with parameters N and alpha_i.
#'
#' @param conj_param List or data.frame: The parameters of the conjugated distribuition (Dirichlet) of the linear predictor.
#' @param outcome Vector or matrix: The observed values at the current time.
#' @param parms List (optional): A list of extra parameters for the model. Not used in this function.
#'
#' @return The predictive log-likelyhood of the observed data.
#' @export
#'
#' @examples
#'
#' conj_param <- matrix(c(1:15), 5, 3, byrow = TRUE)
#'
#' multnom_log_like(conj_param, matrix(5, 5, 3))
multnom_log_like <- function(conj_param, outcome, parms = list()) {
  conj_param <- if (is.null(dim(conj_param))) {
    conj_param %>% t()
  } else {
    conj_param
  }
  outcome <- if (is.null(dim(outcome))) {
    outcome %>% t()
  } else {
    outcome
  }
  T <- nrow(conj_param)
  r <- ncol(conj_param) - 1

  log.like <- rep(NA, T)
  for (t in 1:T) {
    outcome_t <- outcome[t, ]
    N <- sum(outcome_t)
    N <- max(N, 1)

    alpha <- conj_param[t, ] %>% as.numeric()
    alpha0 <- sum(alpha)

    const <- lgamma(alpha0) + lgamma(N + 1) - lgamma(N + alpha0)

    log.like[t] <- const + sum(lgamma(outcome_t + alpha) - lgamma(outcome_t + 1) - lgamma(alpha))
    # lgamma(N + alpha0 - outcome_t) - lgamma(N - outcome_t + 1) - lgamma(alpha0 - alpha)
  }
  log.like
}

#' @export
multnom_kernel <- list(
  "conj_prior" = convert_Dir_Normal,
  "conj_post" = convert_Normal_Dir,
  "update" = update_Dir,
  "smoother" = generic_smoother,
  "pred" = multnom_pred,
  "log.like" = multnom_log_like,
  "offset" = logit_offset,
  "link_function" = function(x) {
    T <- length(x)
    return(log(x[-T, ] / x[T, ]))
  },
  "inv_link_function" = function(x) {
    y <- exp(x)
    x_last <- 1 / (1 + colSums(y))
    return(rbind(y * x_last, x_last))
  },
  "param_names" = function(y) {
    paste0("alpha.", 1:dim(y)[2])
  },
  "multi_var" = TRUE
)
