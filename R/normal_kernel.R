#' normal_filter
#'
#' Filtering function for the Normal model with known variance.
#'
#' @param outcome Vector: The observed values at time t. Must have length greater than 1.
#' @param m0 Vector: The prior mean for the latent vector at time t.
#' @param C0 Matrix: The prior covariance matrix for the latent vector at time t.
#' @param FF Matrix: The regression matrix at time t. Should be compatible with the dimensions of m0 and outcome, i.e., if m0 has dimension n and outcome has dimension m, FF should be n x m.
#' @param G Matrix: The state evolution matrix.
#' @param D Matrix: The discount factor matrix at time t.
#' @param W Matrix: The noise matrix at time t.
#' @param offset Vector: Same dimension as outcome. A vector contaning the offset at time t.
#' @param parms list: a list contaning extra arguments. In this model, sigma2 (covariance matrix for the noise associated withy) should be passed.
#'
#' @return A list containing:
#' \itemize{
#'  \item at: One-step-ahead mean for the latent vectors.
#'  \item Rt: One-step-ahead covariance matrix for the latent vectors.
#'  \item ft: One-step-ahead linear predictor.
#'  \item Qt: One-step-ahead covariance matrix for the linear predictior.
#'  \item mt: The filtered mean of the latent vector at time t.
#'  \item Ct: The filtered covariance matrix of the latent vector at time t.
#'  \item outcome: The observed value at time t.
#'  \item params: The same as the argument.
#'  }
#' @export
#' @importFrom MASS ginv
#' @importFrom rootSolve multiroot
#'
#' @examples
#' outcome <- c(rnorm(1, 0, 1), 0)
#' m0 <- c(2, -log(1.1))
#' C0 <- diag(2)
#' G <- matrix(c(1, 0, 0, 1), 2, 2)
#' FF <- matrix(c(1, 0, 0, 1), 2, 2)
#' D <- diag(c(0.1, 0.1), 2, 2) + 1
#' W <- diag(c(0, 0), 2, 2)
#' offset <- c(0)
#' filtered_data <- GDLM::normal_filter(outcome, m0, C0, FF, G, D, W, offset, parms = list("sigma2" = 1))
normal_filter <- function(outcome, m0, C0, FF, G, D, W, offset = 1, parms = list()) {

  na.flag=any(is.na(outcome)) | offset[r + 1]==0
  na.mult=if(na.flag){0}else{1}

  at <- (G %*% m0)
  Pt <- G %*% C0 %*% (t(G))
  Rt <- as.matrix(D * Pt) + W

  ft <- (t(FF) %*% at) * offset
  Qt <- (as.matrix(t(FF) %*% Rt %*% FF) + parms$sigma2) %*% diag(offset**2)

  At <- Rt %*% FF %*% ginv(Qt)
  mt <- at + At %*% (outcome - ft)
  Ct <- Rt - At %*% t(Rt %*% FF)

  return(list(
    "at" = at, "Rt" = Rt,
    "ft" = ft, "Qt" = Qt,
    "mt" = mt, "Ct" = Ct,
    "outcome" = outcome, "parms" = parms
  ))
}


#' normal_pred
#'
#' Calculate the values for the predictive distribuition given the values of the parameter of the conjugated distribuition of the linear predictor.
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
#' }
#' @export
#' @importFrom extraDistr qbetapr
#'
#' @examples
#' # A fitted model shoulb be used as argument, but you can also pass only the parameter themselves.
#'
#' model <- list(
#'   "tau0" = 0.5,
#'   "tau1" = -1 / 2,
#'   "tau2" = 0,
#'   "tau3" = -2.0,
#'   "Qt" = matrix(1, 1, 1),
#'   "outcome" = 0
#' )
#'
#' normal_pred(model)
normal_pred <- function(model, pred_cred = 0.95) {
  mu <- model$ft
  sigma2 <- model$Qt

  pred <- mu
  var.pred <- as.matrix(sigma2)

  icl.pred <- qnorm((1 - pred_cred) / 2) * sqrt(diag(sigma2)) + mu
  icu.pred <- qnorm(1 - (1 - pred_cred) / 2) * sqrt(diag(sigma2)) + mu

  list(
    "pred"     = pred,
    "var.pred" = var.pred,
    "icl.pred" = icl.pred,
    "icu.pred" = icu.pred,
    "log.vero" = dnorm((model$outcome[1] - mu) / sqrt(sigma2), log = TRUE)
  )
}

#' normal_fit
#'
#' Fit the  normal model given the observed value and the model parameters.
#'
#' @param outcome Matrix: The observed data. It's dimension shoulb be T x m, where T is the length of the time series and m is the number of outcomes at each time.
#' @param m0 Vector: The prior mean for the latent vector.
#' @param C0 Matrix: The prior covariance matrix for the latent vector.
#' @param FF Array: A 3D-array containing the regression matrix for each time. It's dimension should be n x m x T, where n is the number of latent variables, m is the number of outcomes in the model and T is the time series length.
#' @param G Matrix: The state evolution matrix.
#' @param D Array: A 3D-array containing the discount factor matrix for each time. It's dimension should be n x n x T, where n is the number of latent variables and T is the time series length.
#' @param W Array: A 3D-array containing the covariance matrix of the noise for each time. It's dimension should be the same as D.
#' @param offset Matrix: The offset of the model. It's dimension should be the same as outcome.
#' @param pred_cred Numeric: the desired credibility for the credibility interval.
#' @param parms list: a list contaning extra arguments. In this model, sigma2 (covariance matrix for the noise associated withy) should be passed.
#' @param kernel_filter function: the method for filtering.
#'
#' @return A list containing the following values:
#' \itemize{
#'    \item mt Matrix: The filtered mean of the latent variables for each time. Dimensions are n x T.
#'    \item Ct Array: A 3D-array containing the filtered covariance matrix of the latent variable for each time. Dimensions are n x n x T.
#'    \item ft Matrix: The one-step-ahead linear predictor for each time. Dimensions are m x T.
#'    \item qt Array: A 3D-array containing the one-step-ahead covariance matrix for the linear predictor for each time. Dimensions are m x T.
#'    \item FF Array: The same as the argument (same values).
#'    \item G Matrix: The same as the argument (same values).
#'    \item D Array: The same as the argument (same values).
#'    \item W Array: The same as the argument (same values).
#'    \item pred Matrix: The one-step-ahead predictions for each time. Dimensions are m x T.
#'    \item var.pred Matrix: The variance for the one-step-ahead predictions for each time. Dimensions are m x T. Note that, in the multivariate Poisson case, the series are supossed independent, so, in particular, they are uncorrelated.
#'    \item icl.pred Matrix: The lower credibility interval for the prediction at each time. Dimensions are m x T.
#'    \item icu.pred Matrix: The upper credibility interval for the prediction at each time. Dimensions are m x T.
#'    \item pred_cred Numeric: Deprecated
#'    \item offset Vector: The same as the argument (same values).
#'    \item log.vero Matrix: log likelyhood for the filtered parameters based on the one-step ahead prediction.
#'    \item outcome Matrix: The same as the argument outcome (same values).
#'    \item parms: The same as the argument.
#' }
#' @export
#'
#' @examples
#' # Not ideal way: should use fit_model function.
#' T <- 200
#' mean <- 0
#' var <- 1
#' outcome <- matrix(rnorm(T) * sqrt(var) + mean)
#' m0 <- c(0, 0)
#' C0 <- diag(c(1, 1))
#' G <- matrix(c(1, 0, 1, 1), 2, 2)
#' FF <- array(matrix(c(1, 0), 2, 1), c(2, 1, T))
#' D <- array(diag(c(0, 0)), c(2, 2, T)) + 1
#' W <- array(diag(c(0, 0)), c(2, 2, T))
#'
#' fitted_data <- normal_fit(outcome = outcome, m0 = m0, C0 = C0, FF = FF, G = G, D = D, W = W, pred_cred = 0.95, parms = list("sigma2" = 1))
#'
#' plot(outcome)
#' lines(fitted_data$pred[1, ])
normal_fit <- function(outcome, m0 = 0, C0 = 1, FF, G, D, W, offset = outcome * 0, pred_cred = 0.95, parms = list(), kernel_filter = normal_gamma_filter) {

  # Definindo quantidades
  T <- nrow(outcome)
  n <- dim(FF)[1]

  D.aux <- D
  D <- ifelse(D.aux == 0, 1, D.aux)
  r <- dim(FF)[2]
  m0 <- matrix(m0, n, 1)
  C0 <- C0
  mt <- matrix(0, nrow = n, ncol = T)
  Ct <- array(rep(diag(n), T), dim = c(n, n, T))
  Rt <- array(rep(0, T), dim = c(n, n, T))
  ft <- matrix(0, nrow = r, ncol = T)
  at <- matrix(0, nrow = n, ncol = T)
  Qt <- array(0, dim = c(r, r, T))
  # At <- array(0,dim=c(n,r,T))

  pred <- matrix(0, r, T)
  var.pred <- array(0, c(r, r, T))
  icl.pred <- matrix(0, r, T)
  icu.pred <- matrix(0, r, T)
  log.vero <- matrix(0, r, T)

  # f_star <- matrix(0,nrow=T,ncol=r)
  # Q_star <- array(0,c(r,r,T))

  mt <- matrix(0, nrow = n, ncol = T)
  Ct <- array(rep(diag(n), T), dim = c(n, n, T))
  tau <- matrix(NA, nrow = 4, ncol = T)
  tau_star <- matrix(NA, nrow = 4, ncol = T)

  D <- ifelse(D == 0, 1, D)

  last_m <- m0
  last_C <- C0

  for (t in 1:T) {
    filter <- kernel_filter(outcome[t, ], last_m, last_C, matrix(FF[, , t], n, r), G, D[, , t], W[, , t], offset = offset[t, ], parms = parms)

    at[, t] <- filter$at
    Rt[, , t] <- filter$Rt
    ft[, t] <- filter$ft
    Qt[, , t] <- filter$Qt
    mt[, t] <- filter$mt
    Ct[, , t] <- filter$Ct

    last_m <- mt[, t]
    last_C <- Ct[, , t]
  }

  result <- list(
    "mt" = mt, "Ct" = Ct,
    "at" = at, "Rt" = Rt,
    "ft" = ft, "Qt" = Qt,
    "tau" = tau, "tau_star" = tau_star,
    "FF" = FF, "G" = G, "D" = D, "W" = W,
    "pred" = pred, "var.pred" = var.pred, "icl.pred" = icl.pred, "icu.pred" = icu.pred,
    "log.vero" = log.vero,
    "outcome" = outcome, "parms" = parms
  )
  return(result)
}

#' normal_mcmc_fit
#'
#' Fit the  normal model given the observed value and the model parameters using MCMC.
#'
#' @param outcome Matrix: The observed data. It's dimension shoulb be T x m, where T is the length of the time series and m is the number of outcomes at each time.
#' @param m0 Vector: The prior mean for the latent vector.
#' @param C0 Matrix: The prior covariance matrix for the latent vector.
#' @param FF Array: A 3D-array containing the regression matrix for each time. It's dimension should be n x m x T, where n is the number of latent variables, m is the number of outcomes in the model and T is the time series length.
#' @param G Matrix: The state evolution matrix.
#' @param D Array: A 3D-array containing the discount factor matrix for each time. It's dimension should be n x n x T, where n is the number of latent variables and T is the time series length.
#' @param W Array: A 3D-array containing the covariance matrix of the noise for each time. It's dimension should be the same as D.
#' @param offset Matrix: The offset of the model. It's dimension should be the same as outcome.
#' @param pred_cred Numeric: the desired credibility for the credibility interval.
#' @param parms list: a list contaning extra arguments. In this model, extra parameters are not used.
#'
#' @return A list containing the following values:
#' \itemize{
#'    \item mt Matrix: The filtered mean of the latent variables for each time. Dimensions are n x T.
#'    \item Ct Array: A 3D-array containing the filtered covariance matrix of the latent variable for each time. Dimensions are n x n x T.
#'    \item ft Matrix: The one-step-ahead linear predictor for each time. Dimensions are m x T.
#'    \item qt Array: A 3D-array containing the one-step-ahead covariance matrix for the linear predictor for each time. Dimensions are m x T.
#'    \item a Matrix: The alpha parameter for the gamma prior for each time. Dimensions are m x T.
#'    \item b Matrix: The beta parameter for the gamma prior for each time. Dimensions are m x T.
#'    \item a.post Matrix: The alpha parameter for the gamma posteior for each time. Dimensions are m x T.
#'    \item b.post Matrix: The beta parameter for the gamma posteior for each time. Dimensions are m x T.
#'    \item FF Array: The same as the argument (same values).
#'    \item G Matrix: The same as the argument (same values).
#'    \item D Array: The same as the argument (same values).
#'    \item W Array: The same as the argument (same values).
#'    \item pred Matrix: The one-step-ahead predictions for each time. Dimensions are m x T.
#'    \item var.pred Matrix: The variance for the one-step-ahead predictions for each time. Dimensions are m x T. Note that, in the multivariate Poisson case, the series are supossed independent, so, in particular, they are uncorrelated.
#'    \item icl.pred Matrix: The lower credibility interval for the prediction at each time. Dimensions are m x T.
#'    \item icu.pred Matrix: The upper credibility interval for the prediction at each time. Dimensions are m x T.
#'    \item pred_cred Numeric: Deprecated
#'    \item offset Vector: The same as the argument (same values).
#'    \item log.vero Matrix: log likelyhood for the filtered parameters based on the one-step ahead prediction.
#'    \item outcome Matrix: The same as the argument outcome (same values).
#'    \item parms: The same as the argument.
#' }
#' @export
#'
#' @importFrom MASS mvrnorm
#' @importFrom MASS ginv
#' @importFrom MCMCpack riwish
#'
#' @examples
#' # Not ideal way: should use fit_model function.
#' T <- 200
#' mean <- 0
#' var <- 1
#' outcome <- matrix(rnorm(T) * sqrt(var) + mean)
#' m0 <- c(0, 0)
#' C0 <- diag(c(1, 1))
#' G <- matrix(c(1, 0, 1, 1), 2, 2)
#' FF <- array(matrix(c(1, 0), 2, 1), c(2, 1, T))
#' D <- array(diag(c(0, 0)), c(2, 2, T)) + 1
#' W <- array(diag(c(0, 0)), c(2, 2, T))
#'
#' fitted_data <- GDLM::normal_mcmc_fit(outcome = outcome, m0 = m0, C0 = C0, FF = FF, G = G, D = D, W = W, pred_cred = 0.95, parms = list("sigma2" = 1))
#'
#' plot(outcome)
#' lines(fitted_data$pred[1, ])
normal_mcmc_fit <- function(outcome, m0 = 0, C0 = 1, FF, G, D, W, offset = outcome * 0, pred_cred = 0.95,
                            parms = list(
                              "sample_size" = 500,
                              "sigma2_prior_V" = dim(FF)[2] * diag(dim(FF)[2]),
                              "sigma2_prior_nu" = dim(FF)[2],
                              "skip_size" = 10
                            )) {
  T <- nrow(outcome)
  n <- dim(FF)[1]
  r <- dim(FF)[2]

  sigma2 <- parms$sigma2_prior_V / parms$sigma2_prior_nu
  sample_size <- parms$sample_size
  skip_size <- parms$skip_size

  raw_theta_t_sample <- array(NA, c(n, T, sample_size))
  noise_sample <- array(NA, c(r, T, sample_size))
  raw_sigma2_sample <- array(NA, c(r, r, sample_size))
  for (i in 1:sample_size) {
    cat(paste0(i, "            \r"))
    mu_fit <- normal_fit(outcome, m0 = m0, C0 = C0, FF = FF, G = G, D = D, W = W, offset = offset, pred_cred = pred_cred, parms = list("sigma2" = sigma2))

    theta_t <- mvrnorm(n = 1, mu = mu_fit$mt[, T], Sigma = mu_fit$Ct[, , T])
    raw_theta_t_sample[, T, i] <- theta_t

    noise_sample[, T, i] <- outcome[T, ] - t(FF[, , T]) %*% theta_t
    for (t in (T - 1):1) {
      Rt <- mu_fit$Rt[, , t + 1]
      Ct <- mu_fit$Ct[, , t]
      simple_Rt_inv <- Ct %*% t(G) %*% ginv(Rt)

      mu <- mu_fit$mt[, t] + simple_Rt_inv %*% (theta_t - mu_fit$at[, t + 1])
      Sigma <- Ct - simple_Rt_inv %*% Rt %*% t(simple_Rt_inv)

      theta_t <- mvrnorm(n = 1, mu = mu, Sigma = Sigma)
      raw_theta_t_sample[, t, i] <- theta_t

      noise_sample[, t, i] <- outcome[t, ] - t(FF[, , t]) %*% theta_t
    }
    S <- t(noise_sample[, , i]) %*% noise_sample[, , i]
    sigma2 <- riwish(parms$sigma2_prior_nu + T, parms$sigma2_prior_V + S)
    raw_sigma2_sample[, , i] <- sigma2
  }

  theta_t_sample <- raw_theta_t_sample[, , (1:sample_size) %% skip_size == 0]
  sigma2_sample <- raw_sigma2_sample[, , (1:sample_size) %% skip_size == 0]

  mt <- array(NA, c(n, T))
  Ct <- array(NA, c(n, n, T))

  at <- array(NA, c(n, T))
  Rt <- array(NA, c(n, n, T))

  ft <- array(NA, c(n, T))
  Qt <- array(NA, c(n, n, T))

  for (t in 1:T) {
    mt_sample <- t(theta_t_sample[, t, ])
    mt[, t] <- colMeans(mt_sample)
    Ct[, , t] <- cov(mt_sample)

    at_sample <- t(G %*% t(mt_sample))
    at[, t] <- colMeans(at_sample)
    Rt[, , t] <- cov(at_sample)

    ft_sample <- t(t(FF[, , t]) %*% t(at_sample))
    ft[, t] <- colMeans(ft_sample)
    Qt[, , t] <- cov(ft_sample)
  }

  sigma2_mean <- mean(sigma2_sample)
  sigma2_var <- var(sigma2_sample)

  result <- list(
    "mt" = mt, "Ct" = Ct,
    "at" = at, "Rt" = Rt,
    "ft" = ft, "Qt" = Qt,
    "sigma2_mean" = sigma2_mean, "sigma2_var" = sigma2_var,
    "FF" = FF, "G" = G, "D" = D, "W" = W,
    # "pred"=pred, "var.pred"=var.pred, "icl.pred"=icl.pred, "icu.pred"=icu.pred,
    #' log.vero'=log.vero,
    "sigma2_sample" = sigma2_sample, "raw_sigma2_sample" = raw_sigma2_sample,
    "theta_t_sample" = theta_t_sample, "rawtheta_t_sample" = raw_theta_t_sample,
    "outcome" = outcome, "parms" = parms
  )
  return(result)
}

#' @export
normal_kernel <- list(
  "fit" = normal_fit,
  "filter" = normal_filter,
  "smoother" = generic_smoother,
  "pred" = normal_pred,
  "multi_var" = TRUE
)
