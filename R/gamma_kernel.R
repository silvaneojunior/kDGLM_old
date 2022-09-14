#' gamma_filter
#'
#' Filtering function for the Gamma model with alpha known.
#'
#' @param outcome Vector: The observed value at time t.
#' @param m0 Vector: The prior mean for the latent vector at time t.
#' @param C0 Matrix: The prior covariance matrix for the latent vector at time t.
#' @param FF Matrix: The regression matrix at time t. Should be compatible with the dimensions of m0 and outcome, i.e., if m0 has dimension n and outcome has dimension m, FF should be n x m.
#' @param G Matrix: The state evolution matrix.
#' @param D Matrix: The discount factor matrix at time t.
#' @param W Matrix: The noise matrix at time t.
#' @param offset Vector: Same dimension as outcome. A vector contaning the offset at time t.
#' @param parms list: a list contaning extra arguments. In this model, parms must contain a variable named phi, which is equal to alpha.
#'
#' @return A list containing:
#' \itemize{
#'  \item at: One-step-ahead mean for the latent vectors.
#'  \item Rt: One-step-ahead covariance matrix for the latent vectors.
#'  \item ft: One-step-ahead linear predictor.
#'  \item qt: One-step-ahead covariance matrix for the linear predictior.
#'  \item tau0: The alpha parameter of the compatibilized gamma prior (see Ref. Raíra).
#'  \item tau1: The beta parameter of the compatibilized gamma prior (see Ref. Raíra).
#'  \item tau0_star: The alpha parameter of the gamma posteirior (see Ref. Raíra).
#'  \item tau1_star: The beta parameter of the gamma posteirior (see Ref. Raíra).
#'  \item mt: The filtered mean of the latent vector at time t.
#'  \item Ct: The filtered covariance matrix of the latent vector at time t.
#'  \item outcome: The observed value at time t.
#'  \item parms: The same as the argument.
#'  }
#' @export
#'
#' @examples
#' outcome <- 5
#' m0 <- c(4, 1)
#' C0 <- diag(c(1, 1))
#' G <- matrix(c(1, 0, 1, 1), 2, 2)
#' FF <- matrix(c(1, 0), 2, 1)
#' D <- diag(c(0.1, 0.1), 2, 2) + 1
#' W <- diag(c(0, 0), 2, 2)
#' offset <- 1
#'
#' filtered_data <- gamma_filter(outcome, m0, C0, FF, G, D, W, offset, parms = list("alpha" = 1))
gamma_filter <- function(outcome, m0, C0, FF, G, D, W, offset = 1, parms) {
  na.flag=is.na(outcome) | offset==0
  na.mult=if(na.flag){0}else{1}
  at <- G %*% m0
  Rt <- G %*% C0 %*% (t(G)) * D**na.mult + W*na.mult

  reduc_RFF <- Rt %*% FF

  # One-step-ahead prediction
  ft <- t(FF) %*% at + if(na.flag){0}else{log(offset)}
  qt <- t(FF) %*% reduc_RFF

  # Compatibilizing priors

  tau0 <- 1 / (-3+3*sqrt(1+2*qt/3))
  tau1 <- tau0 * exp(ft - qt / 2)

  # Calculating posterior
  if(!na.flag){
    tau0_star <- tau0 + phi
    tau1_star <- tau1 + outcome * phi

    # Compatibilizing posterior

    ft_star <- -digamma(tau0_star) + log(tau1_star)
    qt_star <- trigamma(tau0_star)

    mt <- at + reduc_RFF * as.vector((ft_star - ft) * (1 / (qt)))
    if (length(qt) > 1) {
      Ct <- Rt - (reduc_RFF %*% t(reduc_RFF)) %*% diag((1 - qt_star / qt) * (1 / qt))
    } else {
      Ct <- Rt - (reduc_RFF %*% t(reduc_RFF)) * as.vector((1 - qt_star / qt) * (1 / qt))
    }
  }else{
    tau0_star <- tau0
    tau1_star <- tau1

    ft_star <- ft
    qt_star <- qt

    mt <- at
    Ct <- Rt
  }

  return(list(
    "at" = at, "Rt" = Rt,
    "ft" = ft, "Qt" = qt,
    "tau0" = tau0, "tau1" = tau1,
    "tau0_star" = tau0_star, "tau1_star" = tau1_star,
    "ft_star" = ft_star, "qt_star" = qt_star,
    "mt" = mt, "Ct" = Ct,
    "outcome" = outcome, "parms" = parms
  ))
}

#' gamma_lb_filter
#'
#' Filtering function for the Gamma model with alpha known and using Linear Bayes.
#'
#' @param outcome Vector: The observed value at time t.
#' @param m0 Vector: The prior mean for the latent vector at time t.
#' @param C0 Matrix: The prior covariance matrix for the latent vector at time t.
#' @param FF Matrix: The regression matrix at time t. Should be compatible with the dimensions of m0 and outcome, i.e., if m0 has dimension n and outcome has dimension m, FF should be n x m.
#' @param G Matrix: The state evolution matrix.
#' @param D Matrix: The discount factor matrix at time t.
#' @param W Matrix: The noise matrix at time t.
#' @param offset Vector: Same dimension as outcome. A vector contaning the offset at time t.
#' @param parms list: a list contaning extra arguments. In this model, parms must contain a variable named phi, which is equal to alpha.
#'
#' @return A list containing:
#' \itemize{
#'  \item at: One-step-ahead mean for the latent vectors.
#'  \item Rt: One-step-ahead covariance matrix for the latent vectors.
#'  \item ft: One-step-ahead linear predictor.
#'  \item qt: One-step-ahead covariance matrix for the linear predictior.
#'  \item tau0: The alpha parameter of the compatibilized gamma prior (see Ref. Raíra).
#'  \item tau1: The beta parameter of the compatibilized gamma prior (see Ref. Raíra).
#'  \item tau0_star: The alpha parameter of the gamma posteirior (see Ref. Raíra).
#'  \item tau1_star: The beta parameter of the gamma posteirior (see Ref. Raíra).
#'  \item mt: The filtered mean of the latent vector at time t.
#'  \item Ct: The filtered covariance matrix of the latent vector at time t.
#'  \item outcome: The observed value at time t.
#'  \item parms: The same as the argument.
#'  }
#' @export
#'
#' @examples
#' outcome <- 5
#' m0 <- c(4, 1)
#' C0 <- diag(c(1, 1))
#' G <- matrix(c(1, 0, 1, 1), 2, 2)
#' FF <- matrix(c(1, 0), 2, 1)
#' D <- diag(c(0.1, 0.1), 2, 2) + 1
#' W <- diag(c(0, 0), 2, 2)
#' offset <- 1
#'
#' filtered_data <- gamma_filter(outcome, m0, C0, FF, G, D, W, offset, parms = list("alpha" = 1))
gamma_lb_filter <- function(outcome, m0, C0, FF, G, D, W, offset = 1, parms) {
  at <- G %*% m0
  Rt <- G %*% C0 %*% (t(G)) * D + W
  reduc_RFF <- Rt %*% FF
  # One-step-ahead prediction
  ft <- t(FF) %*% at + log(offset)
  qt <- t(FF) %*% reduc_RFF
  phi <- parms$phi

  mu <- exp(ft + qt / 2)
  var <- exp(qt - 1) * exp(2 * ft + qt)

  # Compatibilizing priors

  tau0 <- (mu**2) / var + 2
  tau1 <- mu * (tau0 - 1)

  # Calculating posterior

  tau0_star <- tau0 + phi
  tau1_star <- tau1 + phi * outcome

  # Compatibilizing posterior

  mu_star <- tau1_star / (tau0_star - 1)
  var_star <- (tau1_star**2) / (((tau0_star - 1)**2) * (tau0_star - 2))

  qt_star <- log(var_star / (mu_star**2)) + 1
  ft_star <- log(mu_star) - qt_star / 2

  mt <- at + reduc_RFF * as.vector((ft_star - ft) * (1 / (qt)))
  if (length(qt) > 1) {
    Ct <- Rt - (reduc_RFF %*% t(reduc_RFF)) %*% diag((1 - qt_star / qt) * (1 / qt))
  } else {
    Ct <- Rt - (reduc_RFF %*% t(reduc_RFF)) * as.vector((1 - qt_star / qt) * (1 / qt))
  }

  return(list(
    "at" = at, "Rt" = Rt,
    "ft" = ft, "Qt" = qt,
    "tau0" = tau0, "tau1" = tau1,
    "tau0_star" = tau0_star, "tau1_star" = tau1_star,
    "ft_star" = ft_star, "qt_star" = qt_star,
    "mt" = mt, "Ct" = Ct,
    "outcome" = outcome, "parms" = parms
  ))
}

#' gamma_pred
#'
#' Calculate the values for the predictive distribuition given the values of the parameter of the conjugated distribuition of the linear predictor.
#' It's worth noting that the credibility interval is calculated as +-2 standard deviations.
#'
#' @param model List: contaning the parameters of the conjugated distribuition for the linear predictor (may be based on the prior, filtered or smoothed distribuition).
#' @param pred_cred Numeric: the desired credibility for the credibility interval.
#'
#' @return A list containing the following values:
#' \itemize{
#'    \item pred vector/matrix: the mean of the predictive distribuition of a next observation. Same type and shape as the parameter in model.
#'    \item var.pred vector/matrix: the variance of the predictive distribuition of a next observation. Same type and shape as the parameter in model.
#'    \item icl.pred vector/matrix: the percentile of 100*((1-pred_cred)/2)% of the predictive distribuition of a next observation. Same type and shape as the parameter in model.
#'    \item icu.pred vector/matrix: the percentile of 100*(1-(1-pred_cred)/2)% of the predictive distribuition of a next observation. Same type and shape as the parameter in model.
#'    \item log.like vector/matrix: the log likelyhood of the observed data based on the predictive distribuition.
#' }
#' @export
#' @importFrom extraDistr qbetapr
#' @importFrom extraDistr dbetapr
#'
#' @examples
#' # A fitted model shoulb be used as argument, but you can also pass only the parameter themselves.
#'
#' model <- list(
#'   "outcome" = 2,
#'   "tau0" = 3,
#'   "tau1" = 4,
#'   "parms" = list("alpha" = 1)
#' )
#'
#' gamma_pred(model)
gamma_pred <- function(model, pred_cred = 0.95) {
  phi <- model$parms$alpha

  tau0 <- phi * model$tau0
  tau1 <- phi * model$tau1
  outcome <- list(
    "pred"     = tau1 / (tau0 - 1),
    "var.pred" = ((tau0 / (tau1 - 1))**2) * (tau0 + phi - 1) / ((tau0 - 2) * phi),
    "icl.pred" = qbetapr((1 - pred_cred) / 2, phi, tau0, tau1 / phi),
    "icu.pred" = qbetapr(1 - (1 - pred_cred) / 2, phi, tau0, tau1 / phi)
  )
  if (all(!is.null(model$outcome))) {
    outcome$log.like <- dbetapr(model$outcome, phi, tau0, tau1 / phi, log = TRUE)
  }
  return(outcome)
}

#' gamma_fit
#'
#' Fit the  gamma model giver the observed value and the model parameters.
#'
#' @param outcome Matrix: The observed data. It's dimension shoulb be T x m, where T is the length of the time series and m is the number of outcomes at each time.
#' @param m0 Vector: The prior mean for the latent vector.
#' @param C0 Matrix: The prior covariance matrix for the latent vector.
#' @param FF Array: A 3D-array containing the regression matrix for each time. It's dimension should be n x m x T, where n is the number of latent variables, m is the number of outcomes in the model and T is the time series length.
#' @param G Matrix: The state evolution matrix.
#' @param D Array: A 3D-array containing the discount factor matrix for each time. It's dimension should be n x n x T, where n is the number of latent variables and T is the time series length.
#' @param W Array: A 3D-array containing the covariance matrix of the noise for each time. It's dimension should be the same as D.
#' @param offset Matrix: The offset of the model. It's dimension should be the same as outcome.
#' @param parms list: a list contaning extra arguments. In this model, parms must contain a variable named phi, which is equal to alpha.
#' @param kernel_filter function: the method for filtering.
#'
#' @return A list containing the following values:
#' \itemize{
#'    \item mt Matrix: The filtered mean of the latent variables for each time. Dimensions are n x T.
#'    \item Ct Array: A 3D-array containing the filtered covariance matrix of the latent variable for each time. Dimensions are n x n x T.
#'    \item ft Matrix: The one-step-ahead linear predictor for each time. Dimensions are m x T.
#'    \item qt Array: A 3D-array containing the one-step-ahead covariance matrix for the linear predictor for each time. Dimensions are m x T.
#'    \item tau0 Matrix: (see Ref. Raíra).
#'    \item tau1 Matrix: (see Ref. Raíra).
#'    \item tau0_star Matrix: (see Ref. Raíra).
#'    \item tau1_star Matrix: (see Ref. Raíra).
#'    \item FF Array: The same as the argument (same values).
#'    \item G Matrix: The same as the argument (same values).
#'    \item D Array: The same as the argument (same values).
#'    \item W Array: The same as the argument (same values).
#'    \item offset Vector: The same as the argument (same values).
#'    \item outcome Matrix: The same as the argument outcome (same values).
#'    \item parms: The same as the argument.
#' }
#' @export
#'
#' @examples
#' # Not ideal way: should use fit_model function.
#' T <- 200
#' w <- (200 / 40) * 2 * pi
#' phi <- 2.5
#' outcome <- matrix(rgamma(T, phi, phi / 20 * (sin(w * 1:T / T) + 2)), T, 1)
#' m0 <- c(0, 0, 0)
#' C0 <- diag(c(1, 1, 1)) * 10
#' G <- as.matrix(Matrix::bdiag(1, matrix(c(cos(w), sin(w), -sin(w), cos(w)), 2, 2)))
#' FF <- array(matrix(c(1, 1, 0), 3, 1), c(3, 1, T))
#' D <- array(diag(c(0.1, 0, 0)), c(3, 3, T)) + 1
#' W <- array(diag(c(0, 0, 0)), c(3, 3, T))
#' offset <- matrix(1, T, 1)
#'
#' fitted_data <- gamma_fit(outcome = outcome, m0 = m0, C0 = C0, FF = FF, G = G, D = D, W = W, offset = offset, parms = list("phi" = phi))
#'
#' plot(outcome)
#' lines(fitted_data$pred[1, ])
gamma_fit <- function(outcome, m0 = 0, C0 = 1, FF, G, D, W, offset, parms = list(), kernel_filter = gamma_filter) {
  if (!("alpha" %in% names(parms)) & ("phi" %in% names(parms))) {
    parms$alpha <- parms$phi
  }
  if (("alpha" %in% names(parms)) & !("phi" %in% names(parms))) {
    parms$phi <- parms$alpha
  }
  if (!("alpha" %in% names(parms)) & !("phi" %in% names(parms))) {
    stop("Erro: alpha/phi must be passed in parms.")
  }

  T <- dim(outcome)[1]
  n <- dim(FF)[1]
  r <- dim(FF)[2]

  D.aux <- D
  D <- ifelse(D.aux == 0, 1, D.aux)

  m0 <- matrix(m0, n, 1)
  C0 <- C0

  # r=1
  at <- matrix(0, nrow = n, ncol = T)
  mt <- matrix(0, nrow = n, ncol = T)
  ft <- matrix(0, nrow = r, ncol = T)
  qt <- matrix(0, nrow = r, ncol = T)
  Ct <- array(rep(diag(n), T), dim = c(n, n, T))
  Rt <- array(rep(diag(n), T), dim = c(n, n, T))
  tau0 <- tau1 <- matrix(0, r, T)
  tau0_star <- tau1_star <- matrix(0, r, T)

  # Prior

  last_m <- m0
  last_C <- C0

  for (t in 1:T) {
    filter <- kernel_filter(outcome[t, ], last_m, last_C, FF[, , t], G, D[, , t], W[, , t], offset[t, ], parms)

    at[, t] <- filter$at
    Rt[, , t] <- filter$Rt
    ft[, t] <- filter$ft
    qt[, t] <- filter$Qt
    tau0[, t] <- filter$tau0
    tau1[, t] <- filter$tau1
    tau0_star[, t] <- filter$tau0_star
    tau1_star[, t] <- filter$tau1_star
    mt[, t] <- filter$mt
    Ct[, , t] <- filter$Ct

    last_m <- mt[, t]
    last_C <- Ct[, , t]
  }

  result <- list(
    "mt" = mt, "Ct" = Ct,
    "ft" = ft, "Qt" = qt,
    "at" = at, "Rt" = Rt,
    "tau0" = tau0, "tau1" = tau1,
    "tau0_star" = tau0_star, "tau1_star" = tau1_star,
    "FF" = FF, "G" = G, "D" = D, "W" = W,
    "offset" = offset, "outcome" = outcome, "parms" = parms
  )
  return(result)
}


#' @export
gamma_kernel <- list(
  "fit" = gamma_fit,
  "filter" = gamma_filter,
  "smoother" = generic_smoother,
  "pred" = gamma_pred,
  "multi_var" = FALSE
)

#' @export
gamma_lb_kernel <- list(
  "fit" = gamma_fit,
  "filter" = gamma_lb_filter,
  "smoother" = generic_smoother,
  "pred" = gamma_pred,
  "multi_var" = FALSE
)
