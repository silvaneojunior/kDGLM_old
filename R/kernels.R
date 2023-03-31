#' generic_smoother
#'
#' Generic smoother for all models.
#'
#' @param mt Matrix: A matrix containing the filtered mean of the latent variables at each time. Each row should represent one variable.
#' @param Ct Array: A 3D-array representing the filtered covariance matrix of the latent variables at each time. The third dimension should represent the time index.
#' @param at Matrix: A matrix containing the one-step-ahead mean of the latent variables at each time based upon the filtered mean. Each row should represent one variable.
#' @param Rt Array: A 3D-array representing the one-step-ahead covariance matrix of the latent variables at each time based upon the filtered covariance matrix. The third dimension should represent the time index.
#' @param G  Matrix: A matrix representing the transition matrix of the model.
#'
#' @return List: The smoothed mean (mts) and covariance (Cts) of the latent variables at each time. Their dimension follows, respectivelly, the dimensions of mt and Ct.
#' @export
#' @importFrom MASS ginv
#'
#' @examples
#' T <- 20
#'
#' mt <- matrix(c(cumsum(rnorm(T) + 1), rep(1, T)), 2, T, byrow = TRUE)
#' Ct <- array(diag(c(1, 1)), c(2, 2, T))
#' G <- matrix(c(1, 0, 1, 1), 2, 2)
#' at <- G %*% mt
#' Rt <- array(G %*% t(G) + diag(c(0.1, 0.1)), c(2, 2, T))
#'
#' smoothed_values <- generic_smoother(mt, Ct, at, Rt, G)
generic_smoother <- function(mt, Ct, at, Rt, G) {
  T <- dim(mt)[2]
  n <- dim(mt)[1]
  mts <- mt
  Cts <- Ct

  # var_index <- matrix(apply(Ct, 3, diag), n, T) != 0

  for (t in (T - 1):1) {
    mt_now <- mt[, t]
    G_now <- G
    G_diff <- rep(0, n)
    if (any(is.na(G))) {
      for (index_col in (1:n)[colSums(is.na(G)) > 0]) {
        index_row <- (1:n)[is.na(G_now[, index_col])]
        G_now[index_row, index_col] <- mt_now[index_col + 1]
        G_now[index_row, index_col + 1] <- mt_now[index_col]
        G_diff[index_col] <- -mt_now[index_col] * mt_now[index_col + 1]
      }
    }

    # var_ref <- var_index[, t]
    # restricted_Rt <- Rt[var_ref, var_ref, t + 1]
    # restricted_Ct <- Ct[var_ref, var_ref, t]
    # print(Ct[, , t])
    restricted_Rt <- Rt[, , t + 1]
    restricted_Ct <- Ct[, , t]

    # simple_Rt_inv <- restricted_Ct %*% t(G_now[var_ref, var_ref]) %*% ginv(restricted_Rt)
    simple_Rt_inv <- restricted_Ct %*% t(G_now[, ]) %*% ginv(restricted_Rt)

    # mts[var_ref, t] <- mt_now[var_ref] + simple_Rt_inv %*% (mts[var_ref, t + 1] - at[var_ref, t + 1])
    # Cts[var_ref, var_ref, t] <- restricted_Ct - simple_Rt_inv %*% (restricted_Rt - Cts[var_ref, var_ref, t + 1]) %*% t(simple_Rt_inv)


    mts[, t] <- mt_now[] + simple_Rt_inv %*% (mts[, t + 1] - at[, t + 1])
    Cts[, , t] <- restricted_Ct - simple_Rt_inv %*% (restricted_Rt - Cts[, , t + 1]) %*% t(simple_Rt_inv)
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
#' @param G Matrix: The state evolution matrix.
#' @param D Array: A 3D-array containing the discount factor matrix for each time. It's dimension should be n x n x T, where n is the number of latent variables and T is the time series length.
#' @param W Array: A 3D-array containing the covariance matrix of the noise for each time. It's dimension should be the same as D.
#' @param p_monit numeric (optional): The prior probability of changes in the latent space variables that are not part of it's dinamic.
#' @param c_monit numeric (optional, if p_monit is not specified): The relative cost of false alarm in the monitoring compared to the cost of not detecting anormalities.
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
analytic_filter <- function(outcomes, m0 = 0, C0 = 1, FF, G, D, W, p_monit = NA, c_monit = 1) {
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
    param_names <- outcomes[[outcome_name]]$family$param_names(outcomes[[outcome_name]]$outcome)
    outcomes[[outcome_name]]$conj_prior_param <- matrix(NA, T, length(param_names)) %>% as.data.frame()
    names(outcomes[[outcome_name]]$conj_prior_param) <- param_names
    outcomes[[outcome_name]]$conj_post_param <- outcomes[[outcome_name]]$conj_prior_param

    outcomes[[outcome_name]]$log.like.null <- rep(NA, T)
    outcomes[[outcome_name]]$log.like.alt <- rep(NA, T)
    outcomes[[outcome_name]]$alt.flags <- rep(0, T)
  }

  c <- c_monit
  p <- if (is.na(p_monit)) {
    0
  } else {
    p_monit
  }
  threshold <- log(c_monit) + log(p) - log(1 - p)

  for (t in 1:T) {
    model_list <- if (is.na(p_monit)) {
      c("null_model")
    } else {
      c("alt_model", "null_model")
    }
    for (model in model_list) {
      D_p <- D[, , t]
      D_p[!D_flags[, , t]] <- D_p[!D_flags[, , t]] * D_mult[[model]]

      next_step <- one_step_evolve(last_m, last_C, G, D_p, W[, , t] + diag(n) * W_add[[model]])
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
      family <- outcome$family
      FF_step <- matrix(FF[, pred_index, t], n, k_i)

      offset_step <- outcome$offset[t, ]
      na.flag <- any(is.null(offset_step) | any(offset_step == 0) | any(is.na(offset_step)))

      for (model in model_list) {
        at_step <- models[[model]]$at_step
        Rt_step <- models[[model]]$Rt_step

        lin_pred <- calc_lin_pred(at_step, Rt_step, FF_step)
        next_step <- family$offset(lin_pred$ft, lin_pred$Qt, if (na.flag) {
          1
        } else {
          offset_step
        })
        ft_canom <- outcome$convert_mat_canom %*% next_step$ft
        Qt_canom <- outcome$convert_mat_canom %*% next_step$Qt %*% t(outcome$convert_mat_canom)

        conj_prior <- family$conj_prior(ft_canom, Qt_canom, parms = outcome$parms)

        models[[model]] <- list(
          "at" = models[[model]]$at,
          "Rt" = models[[model]]$Rt,
          "at_step" = at_step,
          "Rt_step" = Rt_step,
          "ft_step" = next_step$ft,
          "Qt_step" = next_step$Qt,
          "FF_step" = lin_pred$FF,
          "conj_prior" = conj_prior,
          "log_like" = family$log.like(conj_prior, outcome$outcome[t, ], parms = outcome$parms)
        )
      }
      outcomes[[outcome_name]]$log.like.null[t] <- models$null_model$log_like
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
      # ft[, t] <- ft_step <- model$ft_step
      # Qt[, , t] <- Qt_step <- model$Qt_step
      ft_step <- model$ft_step
      Qt_step <- model$Qt_step
      outcomes[[outcome_name]]$conj_prior_param[t, ] <- conj_prior
      if (na.flag) {
        norm_post <- list(ft = ft_step, Qt = Qt_step)
        outcomes[[outcome_name]]$conj_post_param[t, ] <- conj_prior
      } else {
        ft_canom <- outcome$convert_mat_canom %*% ft_step
        Qt_canom <- outcome$convert_mat_canom %*% Qt_step %*% t(outcome$convert_mat_canom)
        conj_prior <- family$conj_prior(ft_canom, Qt_canom, parms = outcome$parms)
        conj_post <- family$update(conj_prior, outcome$outcome[t, ], parms = outcome$parms)

        outcomes[[outcome_name]]$conj_post_param[t, ] <- conj_post
        norm_post <- family$conj_post(conj_post, parms = outcome$parms)

        ft_star <- norm_post$ft <- outcome$convert_mat_default %*% norm_post$ft
        Qt_star <- norm_post$Qt <- outcome$convert_mat_default %*% norm_post$Qt %*% t(outcome$convert_mat_default)

        At <- Ct_step %*% model$FF %*% ginv(Qt_step)
        models[["null_model"]]$at_step <- mt_step <- mt_step + At %*% (ft_star - ft_step)
        models[["null_model"]]$Rt_step <- Ct_step <- Ct_step + At %*% (Qt_star - Qt_step) %*% t(At)
      }

      at[, t] <- model$at
      Rt[, , t] <- model$Rt
      if (length(model_list) == 1) {
        ft[, t] <- ft_step
        Qt[, , t] <- Qt_step
      }
      mt[, t] <- last_m <- mt_step
      Ct[, , t] <- last_C <- Ct_step
    }
  }

  result <- list(
    "mt" = mt, "Ct" = Ct,
    "at" = at, "Rt" = Rt,
    "ft" = ft, "Qt" = Qt,
    "FF" = FF, "G" = G, "D" = D, "W" = W,
    "outcomes" = outcomes, "var_names" = var_names
  )
  return(result)
}

one_step_evolve <- function(m0, C0, G, D, W) {
  n <- dim(G)[1]
  G_now <- G
  G_diff <- rep(0, n)
  if (any(is.na(G))) {
    for (index_col in (1:n)[colSums(is.na(G)) > 0]) {
      index_row <- (1:n)[is.na(G_now[, index_col])]
      G_now[index_row, index_col] <- m0[index_col + 1]
      G_now[index_row, index_col + 1] <- m0[index_col]
      G_diff[index_row] <- G_diff[index_row] - m0[index_col] * m0[index_col + 1]
    }
  }
  at <- (G_now %*% m0) + G_diff
  Pt <- G_now %*% C0 %*% (t(G_now))
  Rt <- as.matrix(D * Pt) + W

  list("at" = at, "Rt" = Rt)
}

calc_lin_pred <- function(at, Rt, FF) {
  n <- dim(FF)[1]
  k <- dim(FF)[2]
  FF_step <- FF
  FF_flags <- is.na(FF)
  at_mod <- at
  at_mod[at_mod == 0] <- 0
  at_matrix <- matrix(at_mod, n, k)
  at_matrix[!FF_flags] <- 1
  charge <- colProd(at_matrix) - colProd(!FF_flags)

  index <- which(FF_flags)
  aux_k <- length(index)
  index_alt <- index[(1:aux_k) + rep(c(1, -1), aux_k / 2)]
  FF_step[index] <- at_matrix[index_alt]

  ft <- (t(FF_step) %*% at)
  Qt <- as.matrix(t(FF_step) %*% Rt %*% FF_step)
  list("ft" = ft - charge, "Qt" = Qt, "FF" = FF_step)
}

ident_offset <- function(ft, Qt, offset) {
  r <- dim(ft)[1]
  t <- dim(ft)[2]
  if (t == 1) {
    return(
      list("ft" = ft * offset, "Qt" = Qt * offset %*% t(offset))
    )
  } else {
    offset_ft <- matrix(offset, r, t)
    offset_Qt <- array(offset %*% t(offset), c(r, r, t))
  }
}

ident_log_offset <- function(ft, Qt, offset) {
  ft[1, ] <- ft[1, ] * offset
  ft[2, ] <- ft[2, ] - log(offset)

  Qt[1, ] <- Qt[1, ] * offset
  Qt[, 1] <- Qt[, 1] * offset
  return(
    list("ft" = ft, "Qt" = Qt)
  )
}

log_offset <- function(ft, Qt, offset) {
  r <- dim(ft)[1]
  t <- dim(ft)[2]
  if (t > 1) {
    offset <- matrix(offset, length(offset), t)
  }
  list("ft" = ft + log(offset), "Qt" = Qt)
}

logit_offset <- function(ft, Qt, offset) {
  r <- dim(ft)[1]
  t <- dim(ft)[2]
  if (t > 1) {
    offset <- matrix(offset, length(offset), t)
    return(list("ft" = ft + log(offset[1:r, ] / offset[r + 1, ]), "Qt" = Qt))
  } else {
    return(list("ft" = ft + log(offset[1:r] / offset[r + 1]), "Qt" = Qt))
  }
}

log_offset_half <- function(ft, Qt, offset) {
  list("ft" = ft + matrix(c(0, log(offset)), 2, dim(ft)[2]), "Qt" = Qt)
}
