#' poisson_filter
#'
#' Filtering function for the Poisson model.
#'
#' @param outcome Vector: The observed value at time t.
#' @param m0 Vector: The prior mean for the latent vector at time t.
#' @param C0 Matrix: The prior covariance matrix for the latent vector at time t.
#' @param FF Matrix: The regression matrix at time t. Should be compatible with the dimensions of m0 and outcome, i.e., if m0 has dimension n and outcome has dimension m, FF should be n x m.
#' @param G Matrix: The state evolution matrix.
#' @param D Matrix: The discount factor matrix at time t.
#' @param W Matrix: The noise matrix at time t.
#' @param offset Vector: Same dimension as outcome. A vector contaning the offset at time t.
#' @param parms list: a list contaning extra arguments. In this model, extra parameters are not used.
#'
#' @return A list containing:
#' \itemize{
#'  \item at: One-step-ahead mean for the latent vectors.
#'  \item Rt: One-step-ahead covariance matrix for the latent vectors.
#'  \item ft: One-step-ahead linear predictor.
#'  \item qt: One-step-ahead covariance matrix for the linear predictior.
#'  \item a: The alpha parameter of the compatibilized gamma prior (see Ref. Raíra).
#'  \item b: The beta parameter of the compatibilized gamma prior (see Ref. Raíra).
#'  \item a.post: The alpha parameter of the gamma posteirior (see Ref. Raíra).
#'  \item b.post: The beta parameter of the gamma posteirior (see Ref. Raíra).
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
#' filtered_data <- poisson_filter(outcome, m0, C0, FF, G, D, W, offset)
poisson_filter <- function(outcome, m0, C0, FF, G, D, W, offset = 1, parms = list()) {
  na.flag=is.na(outcome) | offset==0
  na.mult=if(na.flag){0}else{1}
  at <- G %*% m0
  Rt <- G %*% C0 %*% (t(G)) * D**na.mult + W*na.mult

  reduc_RFF <- Rt %*% FF

  # One-step-ahead prediction
  ft <- t(FF) %*% at + if(na.flag){0}else{log(offset)}
  qt <- t(FF) %*% reduc_RFF

  # Compatibilizing priors

  h=-3+3*sqrt(1+2*qt/3)
  a <- (1 / h)
  b <- a*exp(-ft - 0.5 * qt)

  # Calculating posterior
  if(!na.flag){
    a.post <- a + outcome
    b.post <- b + 1
    # Compatibilizing posterior

    gt <- digamma(a.post) - log(b.post)
    pt <- trigamma(a.post)

    mt <- at + reduc_RFF * as.vector((gt - ft) * (1 / (qt)))
    if (length(qt) > 1) {
      Ct <- Rt - (reduc_RFF %*% t(reduc_RFF)) %*% diag((1 - pt / qt) * (1 / qt))
    } else {
      Ct <- Rt - (reduc_RFF %*% t(reduc_RFF)) * as.vector((1 - pt / qt) * (1 / qt))
    }
  }else{
    a.post <- a
    b.post <- b
    mt <- at
    Ct <- Rt
  }

  return(list(
    "at" = at, "Rt" = Rt,
    "ft" = ft, "Qt" = qt,
    "a" = a, "b" = b,
    "a.post" = a.post, "b.post" = b.post,
    "mt" = mt, "Ct" = Ct,
    "outcome" = outcome, "parms" = parms
  ))
}

#' poisson_lb_filter
#'
#' Filtering function for the Poisson model using Linear Bayes.
#'
#' @param outcome Vector: The observed value at time t.
#' @param m0 Vector: The prior mean for the latent vector at time t.
#' @param C0 Matrix: The prior covariance matrix for the latent vector at time t.
#' @param FF Matrix: The regression matrix at time t. Should be compatible with the dimensions of m0 and outcome, i.e., if m0 has dimension n and outcome has dimension m, FF should be n x m.
#' @param G Matrix: The state evolution matrix.
#' @param D Matrix: The discount factor matrix at time t.
#' @param W Matrix: The noise matrix at time t.
#' @param offset Vector: Same dimension as outcome. A vector contaning the offset at time t.
#' @param parms list: a list contaning extra arguments. In this model, extra parameters are not used.
#'
#' @return A list containing:
#' \itemize{
#'  \item at: One-step-ahead mean for the latent vectors.
#'  \item Rt: One-step-ahead covariance matrix for the latent vectors.
#'  \item ft: One-step-ahead linear predictor.
#'  \item qt: One-step-ahead covariance matrix for the linear predictior.
#'  \item a: The alpha parameter of the compatibilized gamma prior (see Ref. Raíra).
#'  \item b: The beta parameter of the compatibilized gamma prior (see Ref. Raíra).
#'  \item a.post: The alpha parameter of the gamma posteirior (see Ref. Raíra).
#'  \item b.post: The beta parameter of the gamma posteirior (see Ref. Raíra).
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
#' filtered_data <- poisson_filter(outcome, m0, C0, FF, G, D, W, offset)
poisson_lb_filter <- function(outcome, m0, C0, FF, G, D, W, offset = 1, parms = list()) {
  at <- G %*% m0
  Rt <- G %*% C0 %*% (t(G)) * D + W

  reduc_RFF <- Rt %*% FF

  # One-step-ahead prediction
  ft <- t(FF) %*% at + log(offset)
  qt <- t(FF) %*% reduc_RFF

  # Compatibilizing priors

  a <- (1 / qt)
  b <- (exp(-ft + 0.5 * qt) / (qt))

  # Calculating posterior

  a.post <- a + outcome
  b.post <- b + 1

  # Compatibilizing posterior
  gt <- digamma(a.post) - log(b.post)
  pt <- trigamma(a.post)

  mt <- at + reduc_RFF * as.vector((gt - ft) * (1 / (qt)))
  if (length(qt) > 1) {
    Ct <- Rt - (reduc_RFF %*% t(reduc_RFF)) %*% diag((1 - pt / qt) * (1 / qt))
  } else {
    Ct <- Rt - (reduc_RFF %*% t(reduc_RFF)) * as.vector((1 - pt / qt) * (1 / qt))
  }

  return(list(
    "at" = at, "Rt" = Rt,
    "ft" = ft, "Qt" = qt,
    "a" = a, "b" = b,
    "a.post" = a.post, "b.post" = b.post,
    "mt" = mt, "Ct" = Ct,
    "outcome" = outcome, "parms" = parms
  ))
}

#' poisson_pred
#'
#' Calculate the values for the predictive distribuition given the values of the parameter of the conjugated distribuition of the linear predictor.
#' It's worth noting that, since the conjugated distribuition of the linear predictior is Gamma (see Ref. Raíra), then the predictive distribuition is Negative Binomial.
#'
#' @param model List: contaning the parameters of the conjugated distribuition for the linear predictor (may be based on the prior, filtered or smoothed distribuition).
#' @param pred_cred Numeric: the desired credibility for the credibility interval
#'
#' @return A list containing the following values:
#' \itemize{
#'    \item pred vector/matrix: the mean of the predictive distribuition of a next observation. Same type and shape as the parameter in model.
#'    \item var.pred vector/matrix: the variance of the predictive distribuition of a next observation. Same type and shape as the parameter in model.
#'    \item icl.pred vector/matrix: the percentile of 100*((1-pred_cred)/2)% of the predictive distribuition of a next observation. Same type and shape as the parameter in model.
#'    \item icu.pred vector/matrix: the percentile of 100*(1-(1-pred_cred)/2)% of the predictive distribuition of a next observation. Same type and shape as the parameter in model.
#'    \item log.like vector/matrix: the log likelyhood of the observed data based on the predictive distribuition.
#' }
#'
#' @export
#'
#' @examples
#' # A fitted model shoulb be used as argument, but you can also pass only the parameter themselves.
#'
#' model <- list(
#'   "a" = c(1:3),
#'   "b" = c(3:1)
#' )
#'
#' poisson_pred(model)
poisson_pred <- function(model, pred_cred = 0.95) {
  a <- model$a
  b <- model$b
  T <- length(a)

  if(any(b < 10**-40)){
    b=ifelse(b >= 10**-40,b,NA)
    warning('The beta parameter for the predictive distribution is very low (<1e-40) at some times. Predicition for those times are unviable.')
  }

  pred <- a / b
  var.pred <- a * (b + 1) / (b)^2

  icl.pred <- qnbinom((1 - pred_cred) / 2, a, (b / (b + 1)))
  icu.pred <- qnbinom(1 - (1 - pred_cred) / 2, a, (b / (b + 1)))

  if (all(!is.null(model$outcome))) {
    log.like <- dnbinom(model$outcome, a, (b / (b + 1)), log = TRUE)
  } else {
    log.like <- NA
  }
  return(list(
    "pred" = pred %>% matrix(1, T),
    "var.pred" = var.pred %>% matrix(1, T),
    "icl.pred" = icl.pred %>% matrix(1, T),
    "icu.pred" = icu.pred %>% matrix(1, T),
    "log.like" = log.like
  ))
}

#' poisson_fit
#'
#' Fit the  poisson model giver the observed value and the model parameters.
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
#' @param kernel_filter function: the method for filtering.
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
#'    \item outcome Matrix: The same as the argument outcome (same values).
#'    \item parms: The same as the argument.
#' }
#' @export
#'
#' @examples
#' # Not ideal way: should use fit_model function.
#' T <- 200
#' w <- (200 / 40) * 2 * pi
#' outcome <- matrix(rpois(T, 20 * (sin(w * 1:T / T) + 2)), T, 1)
#' m0 <- c(0, 0, 0)
#' C0 <- diag(c(1, 1, 1))
#' G <- as.matrix(Matrix::bdiag(1, matrix(c(cos(w), sin(w), -sin(w), cos(w)), 2, 2)))
#' FF <- array(matrix(c(1, 1, 0), 3, 1), c(3, 1, T))
#' D <- array(diag(c(0.1, 0, 0)), c(3, 3, T)) + 1
#' W <- array(diag(c(0, 0, 0)), c(3, 3, T))
#' offset <- matrix(1, T, 1)
#'
#' fitted_data <- GDLM::poisson_fit(outcome = outcome, m0 = m0, C0 = C0, FF = FF, G = G, D = D, W = W, offset = offset, pred_cred = 0.95)
#'
#' plot(outcome)
#' lines(fitted_data$pred[1, ])
poisson_fit <- function(outcome, m0 = 0, C0 = 1, FF, G, D, W, offset, pred_cred = 0.95, parms = list(), kernel_filter = poisson_filter) {
  T <- dim(outcome)[1]
  n <- dim(FF)[1]
  r <- dim(FF)[2]

  D.aux <- D
  D <- ifelse(D.aux == 0, 1, D.aux)

  m0 <- matrix(m0, n, 1)
  C0 <- C0
  mt <- matrix(0, nrow = n, ncol = T)
  Ct <- array(rep(diag(n), T), dim = c(n, n, T))
  Rt <- array(rep(0, T), dim = c(n, n, T))
  ft <- matrix(0, nrow = r, ncol = T)
  at <- matrix(0, nrow = n, ncol = T)
  qt <- array(0, dim = c(r, r, T))

  # r=1
  at <- matrix(0, nrow = n, ncol = T)
  mt <- matrix(0, nrow = n, ncol = T)
  ft <- matrix(0, nrow = r, ncol = T)
  qt <- matrix(0, nrow = r, ncol = T)
  Ct <- array(rep(diag(n), T), dim = c(n, n, T))
  Rt <- array(rep(diag(n), T), dim = c(n, n, T))
  a <- b <- matrix(0, r, T)
  a.post <- b.post <- matrix(0, r, T)

  norm_ic <- qnorm(1 - (1 - pred_cred) / 2)

  # Prior

  last_m <- m0
  last_C <- C0

  # mean2_FF=rep(1,n)
  #
  # # for(i in 1:n){
  # #   var=FF[i,,]
  # #   print(any(var!=0 & var!=1))
  # #   if(any(var!=0 & var!=1)){
  # #     mean2_FF[i]=sqrt(mean(var**2))
  # #   }
  # # }
  #
  # mean2_FF=diag(1/mean2_FF)

  for (t in 1:T) {
    filter <- kernel_filter(outcome[t, ], last_m, last_C, FF[, , t], G, D[, , t], W[, , t], offset[t, ])

    at[, t] <- filter$at
    Rt[, , t] <- filter$Rt
    ft[, t] <- filter$ft
    qt[, t] <- filter$Qt
    a[, t] <- filter$a
    b[, t] <- filter$b
    a.post[, t] <- filter$a.post
    b.post[, t] <- filter$b.post
    mt[, t] <- filter$mt
    Ct[, , t] <- filter$Ct

    last_m <- filter$mt
    last_C <- filter$Ct
  }
  result <- list(
    "mt" = mt, "Ct" = Ct,
    "ft" = ft, "Qt" = qt,
    "at" = at, "Rt" = Rt,
    "a" = a, "b" = b,
    "a.post" = a.post, "b.post" = b.post,
    "FF" = FF, "G" = G, "D" = D, "W" = W,
    "pred_cred" = pred_cred, "offset" = offset,
    "outcome" = outcome, "parms" = parms
  )
  return(result)
}



#' @export
poisson_kernel <- list(
  "fit" = poisson_fit,
  "filter" = poisson_filter,
  "smoother" = generic_smoother,
  "pred" = poisson_pred,
  "multi_var" = FALSE
)

#' @export
poisson_lb_kernel <- list(
  "fit" = poisson_fit,
  "filter" = poisson_lb_filter,
  "smoother" = generic_smoother,
  "pred" = poisson_pred,
  "multi_var" = FALSE
)
