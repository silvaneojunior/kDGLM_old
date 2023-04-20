#' Dirichlet
#'
#' Creates an outcome with Dirichlet distribuition with the chosen parameters.
#'
#' @param alpha Vector: A vector of names fir the linear preditor associated with the concentration parameter of the gamma distribuition. The parameter is treated as unknowed and equal to the exponential of the associated linear preditor.
#' @param outcome vector: Values of the observed data.
#' @param offset vector: The offset at each observation. Must have the same shape as outcome.
#'
#' @return A object of the class dlm_distr
#' @export
#'
#' @examples
#'
#' T <- 200
#' a <- 10
#' b <- 5
#' data <- rdirichlet(T, c(a, b))
#'
#' alpha_block <- polynomial_block(alpha = 1) * 2
#'
#' outcome <- Dirichlet(alpha = c("alpha_1", "alpha_2"), outcome = data)
#' fitted_data <- fit_model(alpha_block, outcomes = outcome)
#' summary(fitted_data)
#' plot(fitted_data)
#'
Dirichlet <- function(alpha, outcome, offset = outcome**0) {
  t <- dim(outcome)[1]
  r <- dim(outcome)[2]
  k <- dim(outcome)[2]
  if (length(alpha) != k) {
    stop(paste0("Error: Incorrect number of parameter, expected ", k, " got ", length(alpha), "."))
  }
  convert_mat_default <- convert_mat_canom <- diag(k)
  parms <- list()

  distr <- list(
    var_names = alpha,
    r = r,
    k = k,
    l = k,
    t = t,
    offset = matrix(offset, t, r),
    outcome = matrix(outcome, t, r),
    convert_mat_canom = convert_mat_canom,
    convert_mat_default = convert_mat_default,
    parms = parms,
    name = "Dirichlet",
    conj_prior = convert_Dirichlet_Normal,
    update = update_Dirichlet,
    calc_pred = Dirichlet_pred,
    smoother = generic_smoother,
    apply_offset = function(ft, Qt, offset) {
      t <- if.null(dim(ft)[2], 1)
      offset <- matrix(offset, t, r)

      list("ft" = ft + log(t(offset)), "Qt" = Qt)
    },
    link_function = log,
    inv_link_function = exp,
    param_names = function(y) {
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
  class(distr) <- "dlm_distr"
  distr$alt_method <- TRUE
  return(distr)
}

#' convert_Dirichlet_Normal
#'
#' DESCRIPTION
#'
#' @param ft vector: A vector representing the means from the normal distribution.
#' @param Qt matrix: A matrix representing the covariance matrix of the normal distribution.
#' @param parms list: A list of extra known parameters of the distribuition. Not used in this function.
#'
#' @return The parameters of the conjugated distribuition of the linear predictor.
#' @export
convert_Dirichlet_Normal <- function(ft, Qt, parms) {
  return(do.call(c, list(ft, Qt)))
}

#' update_Dirichlet
#'
#' Dirichlet
#'
#' @param conj_prior list: A vector containing the parameters of the Inverse-Gamma (alpha,beta). Not used in the alternative method.
#' @param ft vector: A vector representing the means from the prior distribution.
#' @param Qt matrix: A matrix representing the covariance matrix of the prior distribution.
#' @param y vector: A vector containing the observations.
#' @param parms list: A list of extra known parameters of the distribuition. For this kernel, parms should containg the shape parameter (phi) for the observational gamma model.
#'
#' @importFrom Rfast spdinv
#' @importFrom MASS ginv
#'
#' @return The parameters of the posterior distribution.
#' @export
update_Dirichlet <- function(conj_prior, ft, Qt, y, parms) {
  f0 <- ft
  S0 <- ginv(Qt)

  # log.like=function(x){
  #   alpha=exp(x)
  #
  #   lgamma(sum(alpha))+sum((alpha-1)*log(y)-lgamma(alpha))-0.5*t(x-f0)%*%S0%*%(x-f0)
  # }

  d1.log.like <- function(x) {
    alpha <- exp(c(x))

    (digamma(sum(alpha)) +
      log(y) - digamma(alpha)) * alpha +
      -S0 %*% (x - f0)
  }

  d2.log.like <- function(x) {
    alpha <- exp(c(x))

    mat <- trigamma(sum(alpha)) * alpha %*% t(alpha) +
      diag(digamma(sum(alpha)) * alpha) +
      diag(log(y) * alpha) +
      diag(-trigamma(alpha) * (alpha**2) - digamma(alpha) * alpha) +
      -S0
    mat
  }

  mode <- f_root(d1.log.like, d2.log.like, start = f0)$root
  H <- d2.log.like(mode)
  S <- spdinv(-H)

  return(list("ft" = matrix(mode, length(mode), 1), "Qt" = S))
  # return(list("ft" = ft, "Qt" = Qt))
}

#' Dirichlet_pred
#'
#' Calculate the values for the predictive distribuition given the values of the parameter of the distribuition of the linear predictor.
#' The data is assumed to have Dirichlet distribuition with unknown concentration parameters having log-Normal distribuition.
#' In this scenario, the marginal distribuition of the data is obtained via Monte Carlo.
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
#' }
#'
#' @importFrom extraDistr rdirichlet ddirichlet
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
#' Dirichlet_pred(params)
Dirichlet_pred <- function(conj_param, outcome = NULL, parms = list(), pred_cred = 0.95) {
  pred.flag <- !is.na(pred_cred)
  like.flag <- !is.null(outcome)
  if (!like.flag & !pred.flag) {
    return(list())
  }

  norm_param <- format_param(conj_param)
  ft <- norm_param$ft
  Qt <- norm_param$Qt

  r <- k <- dim(ft)[1]
  t <- dim(ft)[2]

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
    log.like <- rep(NA, t)
  } else {
    log.like <- NULL
  }
  N <- 5000

  outcome <- matrix(outcome, r, t)
  sample <- matrix(rnorm(r * N), N, r)
  for (i in 1:t) {
    ft_i <- sample %*% cholesky(Qt[, , i]) + matrix(ft[, i], N, r, byrow = TRUE)
    sample_y <- rdirichlet(N, alpha = exp(ft_i))
    if (pred.flag) {
      pred[, i] <- colmeans(sample_y)
      var.pred[, , i] <- var(sample_y)
      icl.pred[, i] <- colQuantile(sample_y, (1 - pred_cred) / 2)
      icu.pred[, i] <- colQuantile(sample_y, 1 - (1 - pred_cred) / 2)
    }
    if (like.flag) {
      log.like.list <- ddirichlet(outcome[, i], alpha = exp(t(ft_i)), log = TRUE)
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
