#' Laplace outcome for kDGLM models
#'
#' Creates an outcome with Laplace distribution with the chosen parameters.
#'
#' @param alpha Vector: A vector of names for the linear predictor associated with the mean parameter of the Laplace distribution. The parameter is treated as unknown and equal to the associated linear predictor.
#' @param outcome Vector: Values of the observed data.
#' @param offset Vector: The offset at each observation. Must have the same shape as outcome.
#'
#' @return A object of the class dlm_distr
#' @importFrom extraDistr dlaplace
#' @export
#'
#' @details
#'
#' For evaluating the posterior parameters, we use a modified version of the method proposed in \insertCite{ArtigokParametrico;textual}{kDGLM}.
#'
#' For computational efficiency, we also use a Laplace approximations to obtain the first and second moments of the posterior \insertCite{@see @TierneyKadane1 and @TierneyKadane2 }{kDGLM}.
#'
#' For the details about the implementation see  \insertCite{ArtigoPacote;textual}{kDGLM}.
#'
#' For the detail about the modification of the method proposed in \insertCite{ArtigokParametrico;textual}{kDGLM}, see \insertCite{ArtigoAltMethod;textual}{kDGLM}.
#'
#' @seealso \code{\link{fit_model}}
#' @family {auxiliary functions for a creating outcomes}
#'
#' @examples
#' @references
#'    \insertAllCited{}
Laplace <- function(mu, outcome, offset = outcome**0) {
  t <- length(outcome)
  r <- 1
  k <- 1
  convert_mat_default <- convert_mat_canom <- diag(k)
  parms <- list()

  distr <- list(
    var_names = c(mu),
    r = r,
    k = k,
    l = k,
    t = t,
    offset = matrix(offset, t, r),
    outcome = matrix(outcome, t, r),
    convert_mat_canom = convert_mat_canom,
    convert_mat_default = convert_mat_default,
    convert_canom_flag = FALSE,
    parms = parms,
    name = "Laplace",
    conj_prior = convert_Laplace_Normal,
    conj_prior = convert_Laplace_Normal,
    update = update_Laplace,
    conj_prior = format_ft,
    conj_post = format_param,
    calc_pred = Laplace_pred,
    smoother = generic_smoother,
    apply_offset = function(ft, Qt, offset) {
      list("ft" = ft, "Qt" = Qt)
    },
    link_function = function(x) {
      x
    },
    inv_link_function = function(x) {
      x
    },
    param_names = generic_param_names(k)
  )
  class(distr) <- "dlm_distr"
  distr$alt_method <- TRUE
  return(distr)
}

#' update_Laplace
#'
#' Calculate the (approximated) posterior parameter for the linear predictors, assuming that the observed values came from a Laplace model from which the mean parameter have prior distribution in the Normal family.
#'
#' @param conj_prior list: A vector containing the parameters of the Inverse-Gamma (alpha,beta). Not used in the alternative method.
#' @param ft vector: A vector representing the means from the prior distribution.
#' @param Qt matrix: A matrix representing the covariance matrix of the prior distribution.
#' @param y vector: A vector containing the observations.
#' @param parms list: A list of extra known parameters of the distribution. Not used in this kernel.
#'
#'
#' @return The parameters of the posterior distribution.
#' @keywords internal
#' @family {auxiliary functions for a Laplace outcome}
#' @details
#'
#' For evaluating the posterior parameters, we use a modified version of the method proposed in \insertCite{ArtigokParametrico;textual}{kDGLM}.
#'
#' For computational efficiency, we also use a Laplace approximations to obtain the first and second moments of the posterior \insertCite{@see @TierneyKadane1 and @TierneyKadane2 }{kDGLM}.
#'
#' For the details about the implementation see  \insertCite{ArtigoPacote;textual}{kDGLM}.
#'
#' For the detail about the modification of the method proposed in \insertCite{ArtigokParametrico;textual}{kDGLM}, see \insertCite{ArtigoAltMethod;textual}{kDGLM}.
#'
#' @references
#'    \insertAllCited{}
update_Laplace <- function(conj_prior, ft, Qt, y, parms) {
  c_val <- -Inf
  f <- function(x) {
    log.prob <- -abs(y - x) + dnorm(x, ft, sqrt(Qt), log = TRUE)
    max.prob <- max(log.prob)
    if (max.prob > c_val) {
      c_val <- max.prob
    }

    prob <- exp(log.prob - c_val)

    rbind(
      prob,
      x * prob,
      (x**2) * prob
    )
  }
  # f(seq(ft-3*sqrt(Qt),ft+3*sqrt(Qt),l=100))
  val <- cubintegrate(f, c(-Inf), c(Inf), fDim = 3, nVec = 1000)$integral
  ft <- matrix(val[2] / val[1], 1, 1)
  Qt <- matrix(val[3], 1, 1) / val[1] - ft**2

  return(list("ft" = ft, "Qt" = Qt))
  # return(list("ft" = ft, "Qt" = Qt))
}

#' Laplace_pred
#'
#' Calculate the values for the predictive distribution given the values of the parameter of the distribution of the linear predictor.
#' The data is assumed to have Laplace distribution with unknown concentration parameters having log-Normal distribution.
#' In this scenario, the marginal distribution of the data is obtained via Monte Carlo.
#'
#' @param conj_param List or data.frame: The parameters of the conjugated distributions of the linear predictor.
#' @param outcome Vector or matrix (optional): The observed values at the current time. Not used in this function.
#' @param parms List: A list of extra parameters for the model. For this function, it must contain the shape parameter phi of the observational model.
#' @param pred_cred Numeric: the desired credibility for the credibility interval.
#'
#' @return A list containing the following values:
#' \itemize{
#'    \item pred vector/matrix: the mean of the predictive distribution of a next observation. Same type and shape as the parameter in model.
#'    \item var.pred vector/matrix: the variance of the predictive distribution of a next observation. Same type and shape as the parameter in model.
#'    \item icl.pred vector/matrix: the percentile of 100*((1-pred_cred)/2)% of the predictive distribution of a next observation. Same type and shape as the parameter in model.
#'    \item icu.pred vector/matrix: the percentile of 100*(1-(1-pred_cred)/2)% of the predictive distribution of a next observation. Same type and shape as the parameter in model.
#' }
#'
#' @importFrom extraDistr rlaplace dlaplace
#' @importFrom stats rnorm var
#' @keywords internal
#' @family {auxiliary functions for a Laplace outcome}
Laplace_pred <- function(conj_param, outcome = NULL, parms = list(), pred_cred = 0.95) {
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

  outcome <- matrix(outcome, t, r)
  sample <- matrix(rnorm(r * N), N, r)
  for (i in 1:t) {
    ft_i <- sample %*% var_decomp(Qt[, , i]) + matrix(ft[, i], N, r, byrow = TRUE)
    sample_y <- rlaplace(N, mu = ft_i)
    if (pred.flag) {
      pred[, i] <- mean(sample_y)
      var.pred[, , i] <- var(sample_y)
      icl.pred[, i] <- quantile(sample_y, (1 - pred_cred) / 2)
      icu.pred[, i] <- quantile(sample_y, 1 - (1 - pred_cred) / 2)
    }
    if (like.flag) {
      log.like.list <- dlaplace(outcome[i, ], mu = t(ft_i), log = TRUE)
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
