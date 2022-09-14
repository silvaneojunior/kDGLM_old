#' system_multinom
#'
#' Evaluate the compatibilizing equation for the multinomial model (see Ref. Raíra).
#'
#' @param x vector: current tau values.
#' @param parms list: auxiliary values for the system.
#'
#' @return a vector with the values of the system (see Ref. Raíra).
system_multinom <- function(x, parms) {
  sub_last <- digamma(x[length(x)] - sum(x[-length(x)]))
  digamma_vec <- digamma(x)

  f_all <- parms$ft - digamma_vec[-length(x)] + sub_last
  last_guy <- parms$media.log + digamma_vec[length(x)] - sub_last

  f_all <- c(f_all, last_guy)

  return(f_all)
}

#' multnom_filter
#'
#' Filtering function for the Multinomial model.
#'
#' @param outcome Vector: The observed values at time t. Must have length greater than 1.
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
#'  \item Qt: One-step-ahead covariance matrix for the linear predictior.
#'  \item tau: BLANK (see Ref. Raíra).
#'  \item tau_star: BLANK (see Ref. Raíra).
#'  \item alpha: BLANK (see Ref. Raíra).
#'  \item alpha_Star: BLANK (see Ref. Raíra).
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
#' outcome <- c(3, 2, 5)
#' m0 <- log(c(0.5, 0.5))
#' C0 <- diag(c(1, 1))
#' G <- matrix(c(1, 0, 0, 1), 2, 2)
#' FF <- matrix(c(1, 0, 0, 1), 2, 2)
#' D <- diag(c(0.1, 0.1), 2, 2) + 1
#' W <- diag(c(0, 0), 2, 2)
#' offset <- c(1, 1, 1)
#'
#' filtered_data <- GDLM::multnom_filter(outcome, m0, C0, FF, G, D, W, offset)
multnom_filter <- function(outcome, m0, C0, FF, G, D, W, offset = c(1, 1, 1), parms = list()) {
  r <- dim(FF)[2]

  na.flag=any(is.na(outcome)) | offset[r + 1]==0
  na.mult=if(na.flag){0}else{1}

  at <- (G %*% m0)[, 1]
  Pt <- G %*% C0 %*% (t(G))
  Rt <- as.matrix(D**na.mult * Pt) + W*na.mult

  # One-step ahead prediction
  ft <- (t(FF) %*% at)[, 1] + if(na.flag){0}else{log(offset[1:r] / offset[r + 1])}
  Qt <- as.matrix(t(FF) %*% Rt %*% FF)

  # Compatibilizing priors
  calc_helper <- 1 + sum(exp(ft))

  H <- exp(ft) %*% t(exp(ft)) / (calc_helper**2)
  diag(H) <- -(exp(ft) * calc_helper - (exp(ft)**2)) / (calc_helper**2)

  media.log <-
    -log(calc_helper) + sum(diag(0.5 * (H %*% Qt)))

  parms <- list("ft" = ft, "media.log" = media.log)

  ss1 <- multiroot(f = system_multinom, start = c(rep(0.01, r), 0.01 * (r + 1)), parms = parms)

  tau <- as.numeric(ss1$root)

  alpha <- tau
  alpha[r + 1] <- tau[r + 1] - sum(tau[-r - 1])

  # Calculating posterior
  if(na.flag){
    alpha_star <- alpha
    tau_star <- tau

    f_star <- ft
    Q_star <- Qt

    mt <- at
    Ct <- Rt
  }else{
    alpha_star <- alpha + outcome

    tau_star <- alpha_star
    tau_star[r + 1] <- sum(alpha_star)
    f_star <- digamma(alpha_star[-r - 1]) - digamma(alpha_star[r + 1])
    Q_star <- matrix(trigamma(alpha_star[r + 1]), r, r)
    diag(Q_star) <- trigamma(alpha_star[-r - 1]) + trigamma(alpha_star[r + 1])

    At <- as.matrix(Rt %*% FF %*% ginv(Qt))
    mt <- at + At %*% (f_star - ft)
    Ct <- Rt + At %*% (Q_star - Qt) %*% t(At)
  }

  return(list(
    "at" = at, "Rt" = Rt,
    "ft" = ft, "Qt" = Qt,
    "tau" = tau, "tau_star" = tau_star,
    "alpha" = alpha, "alpha_star" = alpha_star,
    "mt" = mt, "Ct" = Ct,
    "outcome" = outcome, "parms" = parms
  ))
}

#' multnom_pred
#'
#' This function is a temporary placeholder.
#' Calculate the values for the predictive distribuition given the values of the parameter of the conjugated distribuition of the linear predictor.
#' PLACE PREDICTIVE DISTRIBUITION.
#'
#' @param model List: contaning the parameters of the conjugated distribuition for the linear predictor (may be based on the prior, filtered or smoothed distribuition).
#' @param pred_cred Numeric: the desired credibility for the credibility interval
#'
#' @return A list containing the following values:
#' \itemize{
#'    \item pred vector: the mean of the predictive distribuition of a next observation. Same type and shape as the parameter in model.
#'    \item var.pred vector: the variance of the predictive distribuition of a next observation. Same type and shape as the parameter in model.
#'    \item icl.pred vector: the percentile of 100*((1-pred_cred)/2)% of the predictive distribuition of a next observation. Same type and shape as the parameter in model.
#'    \item icu.pred vector: the percentile of 100*(1-(1-pred_cred)/2)% of the predictive distribuition of a next observation. Same type and shape as the parameter in model.
#'    \item log.like vector/matrix: the log likelyhood of the observed data based on the predictive distribuition.
#' }
#' @export
#'
#' @examples
#' # MUST BE MADE.
#' # A fitted model shoulb be used as argument, but you can also pass only the parameter themselves.
#'
#' model <- list(
#'   "alpha" = c(10, 2, 6),
#'   "outcome" = c(3, 2, 5)
#' )
#'
#' multnom_pred(model)
multnom_pred <- function(model, pred_cred = 0.95) {
  if(is.null(dim(model$outcome))){
    model$outcome=matrix(model$outcome,1,length(model$outcome))
    model$alpha=matrix(model$alpha,length(model$alpha),1)
  }
  T <- nrow(model$outcome)
  n <- dim(model$FF)[1]
  r <- ncol(model$outcome) - 1

  pred <- matrix(NA, r+1, T)
  var.pred <- array(NA, c(r+1, r+1, T))
  icl.pred <- matrix(NA, r+1, T)
  icu.pred <- matrix(NA, r+1, T)
  log.like <- rep(NA, T)
  for (t in 1:T) {
    outcome <- model$outcome[t, ]
    N <- sum(outcome)

    alpha <- model$alpha[, t]
    alpha0 <- sum(alpha)

    p <- alpha / alpha0
    p_var <- p * (1 - p) / (alpha0 + 1)

    pred[,t] <- N * p
    var.pred[, , t] <- diag(N * p * (1 - p) * (N + alpha0) / (alpha0 + 1))
    for (i in 2:(r+1)) {
      for (j in 1:(i - 1)) {
        var.pred[i, j, t] <- var.pred[j, i, t] <- -N * p[i] * p[j] * (N + alpha0) / (alpha0 + 1)
      }
    }

    const <- lgamma(alpha0) + lgamma(N + 1) - lgamma(N + alpha0)

    x_mat <- matrix(0:N, N + 1, r+1)
    alpha_mat <- matrix(alpha, N + 1, r+1, byrow = TRUE)
    x_alpha_mat <- x_mat + alpha_mat

    prob_mat <- lgamma(x_alpha_mat) - lgamma(x_mat + 1) - lgamma(alpha_mat) + lgamma(N + alpha0 - x_alpha_mat) - lgamma(N - x_mat + 1) - lgamma(alpha0 - alpha_mat)
    prob_mat <- exp(const + prob_mat)
    for (i in 1:(r+1)) {
      probs_acum <- cumsum(prob_mat[, i])

      icl.pred[i, t] <- sum(probs_acum <= ((1 - pred_cred) / 2)) - 1
      icu.pred[i, t] <- sum(probs_acum <= (1 - (1 - pred_cred) / 2))

      icl.pred[i, t] <- max(0, icl.pred[i])
      icu.pred[i, t] <- min(N, icu.pred[i])
    }
    log.like[t] <- const + sum(lgamma(outcome + alpha) - lgamma(outcome + 1) - lgamma(alpha) + lgamma(N + alpha0 - outcome) - lgamma(N - outcome + 1) - lgamma(alpha0 - alpha))
  }

  list(
    "pred"     = pred,
    "var.pred" = var.pred,
    "icl.pred" = icl.pred,
    "icu.pred" = icu.pred,
    "log.like" = log.like
  )
}

# multnom_pred=function(model,pred_cred=0.95){
#   r=length(model$ft)
#   pre_ps=exp(model$ft)/(sum(exp(model$ft))+1)
#
#   p=pre_ps
#   var=model$Qt
#
#   diag_mult=diag(p*(1-p))
#   cov_mult=diag_mult%*%var%*%diag_mult
#
#   vec_rest=1-sum(p)
#
#   n_total=sum(model$outcome)
#
#   p=c(p,vec_rest)*n_total
#   trans_mat=rbind(diag(r),rep(-1,r))
#   var=(trans_mat%*%cov_mult%*%t(trans_mat))%*%(diag(r+1)*(n_total**2))
#
#   p_i=p-2*sqrt(diag(var))
#   p_s=p+2*sqrt(diag(var))
#
#   list(
#     'pred'     = p,
#     'var.pred' = var,
#     'icl.pred' = p_i,
#     'icu.pred' = p_s
#   )
# }

#' multnom_fit
#'
#'  Fit the multinomial model given the observed value and the model parameters.
#'
#' @param outcome Matrix: The observed data. It's dimension shoulb be T x m, where T is the length of the time series and m is the number of outcomes at each time.
#' @param m0 Vector: The prior mean for the latent vector.
#' @param C0 Matrix: The prior covariance matrix for the latent vector.
#' @param FF Array: A 3D-array containing the regression matrix for each time. It's dimension should be n x m x T, where n is the number of latent variables, m is the number of outcomes in the model and T is the time series length.
#' @param G Matrix: The state evolution matrix.
#' @param D Array: A 3D-array containing the discount factor matrix for each time. It's dimension should be n x n x T, where n is the number of latent variables and T is the time series length.
#' @param W Array: A 3D-array containing the covariance matrix of the noise for each time. It's dimension should be the same as D.
#' @param offset Matrix: The offset of the model. It's dimension should be the same as outcome.
#' @param parms list: a list contaning extra arguments. In this model, extra parameters are not used.
#' @param kernel_filter function: the method for filtering.
#'
#' @return A list containing the following values:
#' \itemize{
#'    \item mt Matrix: The filtered mean of the latent variables for each time. Dimensions are n x T.
#'    \item Ct Array: A 3D-array containing the filtered covariance matrix of the latent variable for each time. Dimensions are n x n x T.
#'    \item ft Matrix: The one-step-ahead linear predictor for each time. Dimensions are m x T.
#'    \item Qt Array: A 3D-array containing the one-step-ahead covariance matrix for the linear predictor for each time. Dimensions are m x T.
#'    \item alpha Matrix: BLANK (see Ref. Raíra)
#'    \item alpha_star Matrix: BLANK (see Ref. Raíra)
#'    \item tau Matrix: BLANK (see Ref. Raíra)
#'    \item tau_star Matrix: BLANK (see Ref. Raíra)
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
#'    \item log_offset Vector: The log offset.
#'    \item outcome Matrix: The same as the argument outcome (same values).
#'    \item parms: The same as the argument.
#' }
#' @export
#' @importFrom MASS ginv
#' @importFrom rootSolve multiroot
#'
#'
#' @examples
#' # Not ideal way: should use fit_model function.
#' T <- 200
#' w <- (200 / 40) * 2 * pi
#' y1 <- matrix(rpois(T, 20 * (sin(w * 1:T / T) + 2)), T, 1)
#' y2 <- matrix(rpois(T, 1:200 / 200 + 1), T, 1)
#' y3 <- matrix(rpois(T, 6), T, 1)
#' outcome <- cbind(y1, y2, y3)
#' m0 <- c(0, 0, 0, 0, 0)
#' C0 <- diag(5)
#' G1 <- as.matrix(Matrix::bdiag(1, matrix(c(cos(w), sin(w), -sin(w), cos(w)), 2, 2)))
#' G2 <- matrix(c(1, 0, 1, 1), 2, 2)
#' G <- as.matrix(Matrix::bdiag(G1, G2))
#' FF <- array(matrix(c(1, 1, 0, 0, 0, 0, 0, 0, 1, 0), 5, 2), c(5, 2, T))
#' D <- array(diag(c(0.1, 0, 0, 0.1, 0)), c(5, 5, T)) + 1
#' W <- array(diag(c(0, 0, 0, 0, 0)), c(5, 5, T))
#' offset <- matrix(1, T, 3)
#'
#' fitted_model <- GDLM::multnom_fit(outcome = outcome, m0 = m0, C0 = C0, FF = FF, G = G, D = D, W = W, offset = offset, pred_cred = 0.95)
#'
#' plot(y1, col = "red", ylim = c(0, max(outcome) * 1.2))
#' points(y2, col = "green")
#' points(y3, col = "blue")
#' lines(fitted_model$pred[1, ], col = "red")
#' lines(fitted_model$pred[2, ], col = "green")
#' lines(fitted_model$pred[3, ], col = "blue")
multnom_fit <- function(outcome, m0 = 0, C0 = 1, FF, G, D, W, offset, parms = list(), kernel_filter = multnom_filter) {
  T <- nrow(outcome)
  n <- dim(FF)[1]

  D.aux <- D
  D <- ifelse(D.aux == 0, 1, D.aux)

  r <- ncol(outcome) - 1
  m0 <- matrix(m0, n, 1)
  C0 <- C0
  mt <- matrix(0, nrow = n, ncol = T)
  Ct <- array(rep(diag(n), T), dim = c(n, n, T))
  Rt <- array(rep(0, T), dim = c(n, n, T))
  ft <- matrix(0, nrow = r, ncol = T)
  at <- matrix(0, nrow = n, ncol = T)
  Qt <- array(0, dim = c(r, r, T))

  # f_star <- matrix(0,nrow=T,ncol=r)
  # Q_star <- array(0,c(r,r,T))

  mt <- matrix(0, nrow = n, ncol = T)
  Ct <- array(rep(diag(n), T), dim = c(n, n, T))
  tau <- matrix(NA, nrow = r + 1, ncol = T)
  alpha <- matrix(NA, nrow = r + 1, ncol = T)
  alpha_star <- matrix(NA, nrow = r + 1, ncol = T)
  tau_star <- matrix(NA, nrow = r + 1, ncol = T)

  D <- ifelse(D == 0, 1, D)

  last_m <- m0
  last_C <- C0

  for (t in 1:T) {
    filter <- kernel_filter(outcome[t, ], last_m, last_C, matrix(FF[, , t], n, r), G, D[, , t], W[, , t], offset = offset[t, ])

    at[, t] <- filter$at
    Rt[, , t] <- filter$Rt
    ft[, t] <- filter$ft
    Qt[, , t] <- filter$Qt
    tau[, t] <- filter$tau
    alpha[, t] <- filter$alpha
    alpha_star[, t] <- filter$alpha_star
    tau[, t] <- filter$tau
    tau_star[, t] <- filter$tau_star
    # f_star[t,]     <- filter$f_star
    # Q_star[,,t]    <- filter$Q_star
    mt[, t] <- filter$mt
    Ct[, , t] <- filter$Ct

    last_m <- mt[, t]
    last_C <- Ct[, , t]
  }

  result <- list(
    "mt" = mt, "Ct" = Ct,
    "ft" = ft, "Qt" = qt,
    "at" = at, "Rt" = Rt,
    "alpha" = alpha, "alpha_star" = alpha_star,
    "tau" = tau, "tau_star" = tau_star,
    "FF" = FF, "G" = G, "D" = D, "W" = W,
    "offset" = offset,
    "outcome" = outcome, "parms" = parms
  )
  return(result)
}

#' @export
multnom_kernel <- list(
  "fit" = multnom_fit,
  "filter" = multnom_filter,
  "smoother" = generic_smoother,
  "pred" = multnom_pred,
  "multi_var" = TRUE
)
