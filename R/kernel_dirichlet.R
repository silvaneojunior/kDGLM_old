#' Dirichlet outcome for kDGLM models
#'
#' Creates an outcome with Dirichlet distribution with the chosen parameters.
#'
#' @param alpha Vector: A vector of names for the linear predictor associated with the concentration parameter of the Dirichlet distribution. The parameter is treated as unknown and equal to the exponential of the associated linear predictor.
#' @param outcome Vector: Values of the observed data.
#' @param offset Vector: The offset at each observation. Must have the same shape as outcome.
#'
#' @return A object of the class dlm_distr
#' @importFrom extraDistr ddirichlet
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
#'
#' T <- 200
#' a <- 10
#' b <- 5
#' data <- extraDistr::rdirichlet(T, c(a, b))
#'
#' alpha_block <- polynomial_block(alpha = 1) * 2
#'
#' outcome <- Dirichlet(alpha = c("alpha_1", "alpha_2"), outcome = data)
#' fitted_data <- fit_model(alpha_block, outcomes = outcome)
#' summary(fitted_data)
#' plot(fitted_data)
#'
#' @references
#'    \insertAllCited{}
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
    convert_canom_flag = FALSE,
    parms = parms,
    name = "Dirichlet",
    conj_prior = convert_Dirichlet_Normal,
    update = update_Dirichlet,
    log.like.cond = function(param, outcome) {
      ddirichlet(outcome, param, log = TRUE)
    },
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
#' This is a dummy function, since, for an Dirichelt outcome, the conjugated prior is not used.
#'
#' @param ft vector: A vector representing the means from the normal distribution.
#' @param Qt matrix: A matrix representing the covariance matrix of the normal distribution.
#' @param parms list: A list of extra known parameters of the distribution. Not used in this function.
#'
#' @return The parameters of the conjugated distribution of the linear predictor.
#' @keywords internal
#' @family {auxiliary functions for a Dirichlet outcome}
convert_Dirichlet_Normal <- function(ft, Qt, parms) {
  return(do.call(c, list(ft, Qt)))
}

#' update_Dirichlet
#'
#' Calculate the (approximated) posterior parameter for the linear predictors, assuming that the observed values came from a Dirichlet model from which the concentration parameters have prior distribution in the log-Normal family.
#'
#' @param conj_prior list: A vector containing the parameters of the Inverse-Gamma (alpha,beta). Not used in the alternative method.
#' @param ft vector: A vector representing the means from the prior distribution.
#' @param Qt matrix: A matrix representing the covariance matrix of the prior distribution.
#' @param y vector: A vector containing the observations.
#' @param parms list: A list of extra known parameters of the distribution. For this kernel, parms should containg the shape parameter (phi) for the observational gamma model.
#'
#'
#' @return The parameters of the posterior distribution.
#' @keywords internal
#' @family {auxiliary functions for a Dirichlet outcome}
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
update_Dirichlet <- function(conj_prior, ft, Qt, y, parms) {
  f0 <- ft
  S0 <- ginv(Qt)

  # log.like=function(x){
  #   alpha=exp(x)
  #
  #   lgamma(sum(alpha))+sum((alpha-1)*log(y)-lgamma(alpha))-0.5*crossprod(x-f0,S0)%*%(x-f0)
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
  S <- ginv(-H)

  return(list("ft" = matrix(mode, length(mode), 1), "Qt" = S))
  # return(list("ft" = ft, "Qt" = Qt))
}

#' Dirichlet_pred
#'
#' Calculate the values for the predictive distribution given the values of the parameter of the distribution of the linear predictor.
#' The data is assumed to have Dirichlet distribution with unknown concentration parameters having log-Normal distribution.
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
#' @importFrom extraDistr rdirichlet ddirichlet
#' @importFrom stats rnorm var
#' @keywords internal
#' @family {auxiliary functions for a Dirichlet outcome}
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
    ft_i <- sample %*% var_decomp(Qt[, , i]) + matrix(ft[, i], N, r, byrow = TRUE)
    sample_y <- rdirichlet(N, alpha = exp(ft_i))
    if (pred.flag) {
      pred[, i] <- colMeans(sample_y)
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
