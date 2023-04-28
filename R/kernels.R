#' generic_smoother
#'
#' Generic smoother for all models.
#'
#' @param mt Matrix: A matrix containing the filtered mean of the latent variables at each time. Each row should represent one variable.
#' @param Ct Array: A 3D-array representing the filtered covariance matrix of the latent variables at each time. The third dimension should represent the time index.
#' @param at Matrix: A matrix containing the one-step-ahead mean of the latent variables at each time based upon the filtered mean. Each row should represent one variable.
#' @param Rt Array: A 3D-array representing the one-step-ahead covariance matrix of the latent variables at each time based upon the filtered covariance matrix. The third dimension should represent the time index.
#' @param G  Array: A 3D-array representing the transition matrix of the model at each time.
#'
#' @return List: The smoothed mean (mts) and covariance (Cts) of the latent variables at each time. Their dimension follows, respectivelly, the dimensions of mt and Ct.
#' @export
#'
#' @importFrom Rfast transpose
#'
#' @examples
#' T <- 20
#'
#' mt <- matrix(c(cumsum(rnorm(T) + 1), rep(1, T)), 2, T, byrow = TRUE)
#' Ct <- array(diag(c(1, 1)), c(2, 2, T))
#' G <- array(c(1, 0, 1, 1), c(2, 2, T))
#' at <- G[, , 1] %*% mt
#' Rt <- array(G[, , 1] %*% t(G[, , 1]) + diag(c(0.1, 0.1)), c(2, 2, T))
#'
#' smoothed_values <- generic_smoother(mt, Ct, at, Rt, G)
generic_smoother <- function(mt, Ct, at, Rt, G) {
  T <- dim(mt)[2]
  n <- dim(mt)[1]
  mts <- mt
  Cts <- Ct

  for (t in (T - 1):1) {
    mt_now <- mt[, t]
    G_sep <- G[, , t + 1]
    G_now <- G_sep
    G_diff <- rep(0, n)
    if (any(is.na(G_sep))) {
      for (index_col in (1:n)[colSums(is.na(G_sep)) > 0]) {
        index_row <- (1:n)[is.na(G_now[, index_col])]
        G_now[index_row, index_col] <- mt_now[index_col + 1]
        G_now[index_row, index_col + 1] <- mt_now[index_col]
        G_diff[index_row] <- G_diff[index_row] - mt_now[index_col] * mt_now[index_col + 1]
      }
    }
    restricted_Rt <- Rt[, , t + 1]
    restricted_Ct <- Ct[, , t]

    simple_Rt_inv <- restricted_Ct %*% transpose(G_now) %*% ginv(restricted_Rt)
    simple_Rt_inv_t <- transpose(simple_Rt_inv)
    # simple_Rt_inv_t=solve(restricted_Rt,G_now %*% restricted_Ct)
    # simple_Rt_inv <- transpose(simple_Rt_inv_t)

    mts[, t] <- mt_now + simple_Rt_inv %*% (mts[, t + 1] - at[, t + 1])
    # Cts[, , t] <- restricted_Ct + simple_Rt_inv %*% (Cts[, , t + 1] - restricted_Rt) %*% simple_Rt_inv_t
    Cts[, , t] <- restricted_Ct + crossprod(simple_Rt_inv_t, (Cts[, , t + 1] - restricted_Rt)) %*% simple_Rt_inv_t
  }
  return(list("mts" = mts, "Cts" = Cts))
}

#' analytic_filter
#'
#' Fit the model given the observed value and the model parameters.
#'
#' @param outcomes List: The observed data. It should contain objects of the class dlm_distr.
#' @param m0 Vector: The prior mean for the latent vector.
#' @param C0 Matrix: The prior covariance matrix for the latent vector.
#' @param FF Array: A 3D-array containing the regression matrix for each time. It's dimension should be n x m x T, where n is the number of latent variables, m is the number of outcomes in the model and T is the time series length.
#' @param G Array: A 3D-array containing the state evolution matrix at each time.
#' @param D Array: A 3D-array containing the discount factor matrix for each time. It's dimension should be n x n x T, where n is the number of latent variables and T is the time series length.
#' @param W Array: A 3D-array containing the covariance matrix of the noise for each time. It's dimension should be the same as D.
#' @param p_monit numeric (optional): The prior probability of changes in the latent space variables that are not part of it's dinamic.
#' @param c_monit numeric (optional, if p_monit is not specified): The relative cost of false alarm in the monitoring compared to the cost of not detecting anormalities.
#'
#' @importFrom Rfast transpose
#'
#' @return A list containing the following values:
#' \itemize{
#'    \item mt Matrix: The filtered mean of the latent variables for each time. Dimensions are n x T.
#'    \item Ct Array: A 3D-array containing the filtered covariance matrix of the latent variable for each time. Dimensions are n x n x T.
#'    \item at Matrix: The one-step-ahead mean of the latent variables at each time. Dimensions are n x T.
#'    \item Rt Array: A 3D-array containing the one-step-ahead covariance matrix for latent variables at each time. Dimensions are n x n x T.
#'    \item ft Matrix: The one-step-ahead linear predictor for each time. Dimensions are m x T.
#'    \item Qt Array: A 3D-array containing the one-step-ahead covariance matrix for the linear predictor for each time. Dimensions are m x m x T.
#'    \item FF Array: The same as the argument (same values).
#'    \item G Matrix: The same as the argument (same values).
#'    \item D Array: The same as the argument (same values) when there is no monitoring. When monitoring for anormalities, the value in times where anormalities were detected is increased.
#'    \item W Array: The same as the argument (same values).
#'    \item outcome List: The same as the argument outcome (same values).
#'    \item var_names Vector: The names of the linear predictors.
#' }
analytic_filter <- function(outcomes, m0 = 0, C0 = 1, FF, G, G_labs, D, W, p_monit = NA, c_monit = 1) {
  # Defining quantities
  T <- dim(FF)[3]
  r <- sum(sapply(outcomes, function(x) {
    x$r
  }))
  n <- dim(FF)[1]
  k <- dim(FF)[2]

  var_names <- colnames(FF)
  D_flags <- (D == 0)
  D <- ifelse(D == 0, 1, D)
  m0 <- matrix(m0, n, 1)
  C0 <- C0
  mt <- matrix(NA, nrow = n, ncol = T)
  Ct <- array(NA, dim = c(n, n, T))
  Rt <- array(NA, dim = c(n, n, T))
  ft <- matrix(NA, nrow = k, ncol = T)
  at <- matrix(NA, nrow = n, ncol = T)
  Qt <- array(NA, dim = c(k, k, T))


  last_m <- m0
  last_C <- C0
  models <- list()
  D_mult <- list("null_model" = 1, "alt_model" = 100)
  W_add <- list("null_model" = 0, "alt_model" = 0.001)
  monit_win <- 1

  for (outcome_name in names(outcomes)) {
    if (class(outcomes[[outcome_name]]) != "dlm_distr") {
      stop(paste0("Error: Outcome contains is not of the right class Expected a dlm_distr object, got a ", class(outcomes[[outcome_name]]), " object."))
    }
    param_names <- outcomes[[outcome_name]]$param_names(outcomes[[outcome_name]]$outcome)
    outcomes[[outcome_name]]$conj_prior_param <- matrix(NA, T, length(param_names)) %>% as.data.frame()
    names(outcomes[[outcome_name]]$conj_prior_param) <- param_names

    outcomes[[outcome_name]]$log.like.null <- rep(NA, T)
    outcomes[[outcome_name]]$log.like.alt <- rep(NA, T)
    outcomes[[outcome_name]]$alt.flags <- rep(0, T)


    pred_index <- match(outcomes[[outcome_name]]$var_names, var_names)
    k_i <- length(pred_index)

    outcomes[[outcome_name]]$ft <- matrix(NA, nrow = k_i, ncol = T)
    outcomes[[outcome_name]]$Qt <- array(NA, dim = c(k_i, k_i, T))
  }

  c <- c_monit
  p <- if (is.na(p_monit)) {
    0
  } else {
    p_monit
  }
  threshold <- log(c_monit) + log(p) - log(1 - p)
  last_C_D <- last_C

  for (t in 1:T) {
    model_list <- if (is.na(p_monit)) {
      c("null_model")
    } else {
      c("alt_model", "null_model")
    }
    for (model in model_list) {
      D_p <- D[, , t]
      D_p[!D_flags[, , t]] <- D_p[!D_flags[, , t]] * D_mult[[model]]
      next_step <- one_step_evolve(last_m, last_C, G[, , t] %>% matrix(n, n), G_labs, D_p**0, W[, , t] + diag(n) * W_add[[model]] + last_C_D * (D_p - 1))


      models[[model]] <- list(
        "at" = next_step$at,
        "Rt" = next_step$Rt,
        "at_step" = next_step$at,
        "Rt_step" = next_step$Rt
      )
    }
    for (outcome_name in names(outcomes)) {
      outcome <- outcomes[[outcome_name]]
      pred_index <- match(outcome$var_names, var_names)
      k_i <- length(pred_index)
      FF_step <- matrix(FF[, pred_index, t], n, k_i)

      offset_step <- outcome$offset[t, ]
      na.flag <- any(is.null(offset_step) | any(offset_step == 0) | any(is.na(offset_step)) | any(is.na(outcome$outcome[t, ])))
      for (model in model_list) {
        at_step <- models[[model]]$at_step
        Rt_step <- models[[model]]$Rt_step

        lin_pred <- calc_lin_pred(at_step, Rt_step, FF_step)

        if (outcome$convert_canom_flag) {
          ft_canom <- outcome$convert_mat_canom %*% lin_pred$ft
          Qt_canom <- outcome$convert_mat_canom %*% lin_pred$Qt %*% transpose(outcome$convert_mat_canom)
        } else {
          ft_canom <- lin_pred$ft
          Qt_canom <- lin_pred$Qt
        }

        next_step <- outcome$apply_offset(ft_canom, Qt_canom, if (na.flag) {
          1
        } else {
          offset_step
        })

        conj_prior <- outcome$conj_prior(next_step$ft, next_step$Qt, parms = outcome$parms)

        models[[model]] <- list(
          "at" = models[[model]]$at,
          "Rt" = models[[model]]$Rt,
          "at_step" = at_step,
          "Rt_step" = Rt_step,
          "ft_step" = next_step$ft,
          "Qt_step" = next_step$Qt,
          "FF_step" = lin_pred$FF,
          "conj_prior" = conj_prior,
          "log_like" = outcome$calc_pred(conj_prior, if (is.na(p_monit)) {
            NULL
          } else {
            outcome$outcome[t, ]
          }, parms = outcome$parms, pred_cred = NA)$log.like
        )
      }
      outcomes[[outcome_name]]$log.like.null[t] <- if (is.na(p_monit)) {
        NA
      } else {
        models$null_model$log_like
      }
      outcomes[[outcome_name]]$log.like.alt[t] <- if (is.na(p_monit)) {
        -Inf
      } else {
        models$alt_model$log_like
      }
      if (monit_win > 0 & !is.na(p_monit)) {
        bayes_factor <- sum(outcomes[[outcome_name]]$log.like.null[t:(t - monit_win + 1)] - outcomes[[outcome_name]]$log.like.alt[t:(t - monit_win + 1)])
      } else {
        bayes_factor <- -1e-10
      }
      bayes_factor <- ifelse(is.nan(bayes_factor), 0, bayes_factor)
      if (bayes_factor < threshold) {
        model <- models$alt_model
        conj_prior <- models$alt_model$conj_prior
        D[, , t] <- D[, , t] * D_mult$alt_model
        W[, , t] <- W[, , t] * W_add$alt_model
        monit_win <- -5
        outcomes[[outcome_name]]$alt.flags[t] <- 1
      } else {
        outcomes[[outcome_name]]$alt.flags[t] <- 0
        model <- models$null_model
        if (bayes_factor < 0) {
          monit_win <- monit_win + 1
        } else {
          monit_win <- 1
        }
      }

      mt_step <- model$at_step
      Ct_step <- model$Rt_step
      ft_step <- model$ft_step
      Qt_step <- model$Qt_step
      if (na.flag) {
        norm_post <- list(ft = ft_step, Qt = Qt_step)
      } else {
        if (outcome$convert_canom_flag) {
          ft_canom <- outcome$convert_mat_canom %*% ft_step
          Qt_canom <- outcome$convert_mat_canom %*% Qt_step %*% transpose(outcome$convert_mat_canom)
        } else {
          ft_canom <- ft_step
          Qt_canom <- Qt_step
        }
        conj_prior <- outcome$conj_prior(ft_canom, Qt_canom, parms = outcome$parms)
        conj_post <- outcome$update(conj_prior, ft = ft_canom, Qt = Qt_canom, y = outcome$outcome[t, ], parms = outcome$parms)

        if (outcome$alt_method) {
          norm_post <- conj_post
        } else {
          norm_post <- outcome$conj_post(conj_post, parms = outcome$parms)
        }

        if (outcome$convert_canom_flag) {
          ft_star <- norm_post$ft <- outcome$convert_mat_default %*% norm_post$ft
          Qt_star <- norm_post$Qt <- outcome$convert_mat_default %*% norm_post$Qt %*% transpose(outcome$convert_mat_default)
        } else {
          ft_star <- norm_post$ft <- norm_post$ft
          Qt_star <- norm_post$Qt <- norm_post$Qt
        }
        At <- Ct_step %*% model$FF %*% ginv(Qt_step)
        At_t <- t(At)
        # At_t <- solve(Qt_step,t(model$FF) %*% Ct_step)
        # At <- t(At_t)
        models[["null_model"]]$at_step <- mt_step <- mt_step + At %*% (ft_star - ft_step)
        models[["null_model"]]$Rt_step <- Ct_step <- Ct_step + crossprod(At_t, (Qt_star - Qt_step)) %*% At_t

        last_C_D <- Ct_step
      }

      at[, t] <- model$at
      Rt[, , t] <- model$Rt

      lin_pred_ref <- calc_lin_pred(model$at, model$Rt, FF[, , t])
      ft[, t] <- lin_pred_ref$ft
      Qt[, , t] <- lin_pred_ref$Qt

      for (outcome_name in names(outcomes)) {
        outcome <- outcomes[[outcome_name]]
        pred_index <- match(outcome$var_names, var_names)

        offset_step <- outcome$offset[t, ]
        na.flag <- any(is.null(offset_step) | any(offset_step == 0) | any(is.na(offset_step)))

        lin_pred <- list(
          ft = lin_pred_ref$ft[pred_index, , drop = FALSE],
          Qt = lin_pred_ref$Qt[pred_index, pred_index, drop = FALSE]
        )

        if (outcome$convert_canom_flag) {
          ft_canom <- outcome$convert_mat_canom %*% lin_pred$ft
          Qt_canom <- outcome$convert_mat_canom %*% lin_pred$Qt %*% transpose(outcome$convert_mat_canom)
        } else {
          ft_canom <- lin_pred$ft
          Qt_canom <- lin_pred$Qt
        }

        next_step <- outcome$apply_offset(ft_canom, Qt_canom, if (na.flag) {
          1
        } else {
          offset_step
        })

        conj_prior <- outcome$conj_prior(next_step$ft, next_step$Qt, parms = outcome$parms)
        outcomes[[outcome_name]]$conj_prior_param[t, ] <- conj_prior
      }

      mt[, t] <- last_m <- mt_step
      Ct[, , t] <- last_C <- Ct_step
      eigen_Ct=eigen(last_C)$values
      pos_eigen_Ct=eigen_Ct[eigen_Ct>0]
      if(log10(max(pos_eigen_Ct))-log10(min(pos_eigen_Ct))>15){
        warning('Covariance for the latent variables is numerical instable. Results may be unreliable.')
      }
    }
  }

  result <- list(
    "mt" = mt, "Ct" = Ct,
    "at" = at, "Rt" = Rt,
    # Consultar Mariane sobre oq fazer com o preditor linear.
    "ft" = ft, "Qt" = Qt,
    "FF" = FF, "G" = G, "D" = D, "W" = W,
    "outcomes" = outcomes, "var_names" = var_names
  )
  return(result)
}

#' one_step_evolve
#'
#' @importFrom Rfast transpose
#'
one_step_evolve <- function(m0, C0, G, G_labs, D, W) {
  n <- dim(G)[1]
  G_now <- G
  G_diff <- rep(0, n)
  if (any(G_labs!='const')) {
    for (index_col in (1:n)[colSums(is.na(G)) > 0]) {
      index_row <- (1:n)[is.na(G_now[, index_col])]
      if(G_labs[index_row, index_col]=='constrained'){
        p=1/(1+exp(-m0[index_col + 1]))
        G_now[index_row, index_col] <- (2*p-1)
        G_now[index_row, index_col + 1] <- 2*p*(1-p)*m0[index_col]
        G_diff[index_row] <- G_diff[index_row] -2*p*(1-p)*m0[index_col]*m0[index_col+1]
      }else{
        G_now[index_row, index_col] <- m0[index_col+1]
        G_now[index_row, index_col + 1] <- m0[index_col+1]*m0[index_col]
        G_diff[index_row] <- G_diff[index_row] - m0[index_col+1]*m0[index_col]
      }
    }
  }
  at <- (G_now %*% m0) + G_diff
  G_now_t <- transpose(G_now)
  Pt <- G_now %*% C0 %*% G_now_t
  Rt <- as.matrix(D * Pt) + W

  list("at" = at, "Rt" = Rt)
}

calc_lin_pred <- function(at, Rt, FF) {
  if (is.null(dim(FF))) {
    FF <- matrix(FF, length(FF), 1)
  }
  n <- dim(FF)[1]
  k <- dim(FF)[2]

  charge <- matrix(0, k, 1)
  at_mod <- at[, 1]
  for (i in 1:k) {
    FF_step <- FF[, i]
    FF_flags <- is.na(FF_step)
    vals_na <- (1:n)[FF_flags]
    n_na <- length(vals_na)
    if (n_na > 1) {
      vals_coef <- vals_na[((1:(n_na / 2) * 2) - 1)]
      vals_noise <- vals_na[((1:(n_na / 2)) * 2)]
      flags_na <- diag(Rt)[vals_noise] == W_correl & diag(Rt)[vals_coef] > 0

      FF[vals_coef, i] <- ifelse(flags_na,
        (at_mod[vals_noise]) * exp(at_mod[vals_coef]),
        (at_mod[vals_noise])
      )
      FF[vals_noise, i] <- ifelse(flags_na,
        exp(at_mod[vals_coef]),
        at_mod[vals_coef]
      )
      charge[i, 1] <- charge[i, 1] + sum(ifelse(flags_na,
        (at_mod[vals_noise]) * exp(at_mod[vals_coef]),
        (at_mod[vals_noise]) * at_mod[vals_coef]
      ))
    }
  }

  ft <- crossprod(FF, at)
  Qt <- as.matrix(crossprod(FF, Rt) %*% FF)
  list("ft" = ft - charge, "Qt" = Qt, "FF" = FF)
}

format_ft <- function(ft, Qt, parms) {
  return(do.call(c, list(ft, Qt)))
}

format_param <- function(conj_param, parms) {
  if (is.null(dim(conj_param))) {
    r_star <- length(conj_param)
    r <- (-1 + sqrt(1 + 4 * r_star)) / 2
    t <- 1
    ft <- conj_param[1:r] %>% matrix(r, t)
    Qt <- conj_param[(r + 1):(r * r + r)] %>% array(c(r, r, t))
  } else {
    r_star <- dim(conj_param)[2]
    t <- dim(conj_param)[1]
    r <- (-1 + sqrt(1 + 4 * r_star)) / 2
    ft <- conj_param[, 1:r] %>%
      data.frame.to_matrix() %>%
      t()
    Qt <- conj_param[, (r + 1):(r * r + r)] %>%
      data.frame.to_matrix() %>%
      t() %>%
      array(c(r, r, t))
  }
  return(list("ft" = ft, "Qt" = Qt))
}
