#' fit_model
#'
#' Fit a model given it's structure and the observed data. This function can be used for any supported family (see vignette).
#'
#' @param ... <undefined class> or list: The structural block of the model.
#' @param outcome Matrix: A matrix containing the observed data. Dimensions should be T x k, where T is the time series length and k is the number of outcomes.
#' @param family String or <undefined class>: The list of functions (or it's name) used to fit the data. Should be choosed based on the distribution of the outcomes.
#' @param offset Matrix: A matrix containing the offset value for the data. Dimesions should be the same as outcome.
#' @param pred_cred Numeric: A number between 0 and 1 (not included) indicanting the credibility interval for predictions. If not within the valid range of values, predicitions are not made.
#' @param smooth_flag Bool: A flag indicating if the smoothed distribuition for the latent variables should be calculated.
#'
#' @return The fitted model (an object of the <undefined class> or list). Contains the values of the estimated parameter and some extra info regarding the quality of the fit.
#' @export
#'
#' @examples
#' library(GDLM)
#'
#' # Normal with unkown variance case
#' T <- 200
#' mu <- rnorm(T, 0, 0.1)
#' outcome <- rnorm(T, cumsum(mu))
#'
#' level <- polynomial_block(
#'   order = 1,
#'   values = c(1, 0),
#'   D = 1 / 0.98,
#'   by_time = FALSE
#' )
#' variance <- polynomial_block(
#'   order = 1,
#'   values = c(0, 1),
#'   D = 1 / 1,
#'   by_time = FALSE
#' )
#'
#' fitted_data <- fit_model(level, variance, outcome = outcome, family = "normal_gamma")
#' show_fit(fitted_data, smooth = TRUE)$plot
#'
#' # Poisson case
#' T <- 200
#' w <- (200 / 40) * 2 * pi
#' outcome <- rpois(T, 20 * (sin(w * 1:T / T) + 2))
#'
#' level <- polynomial_block(order = 1, values = 1, D = 1 / 0.95)
#' season <- harmonic_block(period = 40, values = 1, D = 1 / 0.98)
#'
#' fitted_data <- fit_model(level, season, outcome = outcome, family = "Poisson")
#' show_fit(fitted_data, smooth = TRUE)$plot
#'
#'
#' # Multinomial case
#' T <- 200
#' y1 <- rpois(T, exp(5 + (-T:T / T) * 5))
#' y2 <- rpois(T, exp(6 + (-T:T / T) * 5 + sin((-T:T) * (2 * pi / 12))))
#' y3 <- rpois(T, exp(5))
#'
#' y <- cbind(y1, y2, y3)
#'
#' level <- polynomial_block(2, k = 2)
#' season <- harmonic_block(12, values = c(0, 1), by_time = FALSE)
#'
#' fitted_data <- fit_model(level, season, outcome = y, family = "Multinomial", pred_cred = 0.95)
#'
#' show_fit(fitted_data, smooth = TRUE)$plot
#'
#' # Gamma case
#' T <- 200
#' w <- (200 / 40) * 2 * pi
#' phi <- 2.5
#' outcome <- matrix(rgamma(T, phi, phi / (20 * (sin(w * 1:T / T) + 2))), T, 1)
#'
#' level <- polynomial_block(order = 1, values = 1, D = 1 / 0.95)
#' season <- harmonic_block(period = 40, values = 1, D = 1 / 0.98)
#'
#' fitted_data <- fit_model(level, season, outcome = outcome, family = "Gamma", parms = list("phi" = phi))
#'
#' show_fit(fitted_data, smooth = TRUE)$plot
fit_model <- function(..., outcome, family, offset = outcome * 0 + 1, parms = list(), pred_cred = 0.95, smooth_flag = TRUE, p_monit = NA, c_monit = 1) {
  if (typeof(family) == typeof("family")) {
    if (tolower(family) == "fgamma") {
      warning("This family is numerically unstable. Be careful with it's usage.")
    }
    family <- kernel_list[[tolower(family)]]
  }
  if (is.null(dim(outcome))) {
    outcome <- matrix(outcome, length(outcome), 1)
  }
  if (is.null(dim(offset))) {
    offset <- matrix(offset, length(offset), 1)
  }
  if (typeof(outcome) == typeof(list())) {
    outcome <- as.matrix(outcome)
  }

  structure <- block_merge(...)
  if (any(dim(outcome) != dim(offset))) {
    stop("Erro: offset does not have the same dim as outcome.")
  }

  if (structure$t == 1) {
    structure$t <- dim(outcome)[1]
    structure$D <- array(structure$D, c(structure$n, structure$n, structure$t))
    structure$W <- array(structure$W, c(structure$n, structure$n, structure$t))
    structure$FF <- array(structure$FF, c(structure$n, structure$k, structure$t))
  }
  if (dim(outcome)[1] != structure$t) {
    stop(paste0("Error: outcome does not have the same time length as structure: got ", dim(outcome)[1], " from outcome, expected ", structure$t))
  }
  # if (dim(outcome)[2] != structure$k) {
  #   stop(paste0("Error: outcome does not have the same number of outcomes as structure: got ", dim(outcome)[2], " from outcome, expected ", structure$k))
  # }
  model <- analytic_filter(
    outcome = outcome,
    m0 = structure$m0,
    C0 = structure$C0,
    FF = structure$FF,
    G = structure$G,
    D = structure$D,
    W = structure$W,
    offset = offset,
    parms = parms,
    family = family,
    p_monit = p_monit,
    c_monit = c_monit
  )
  if (smooth_flag) {
    smoothed <- family$smoother(model$mt, model$Ct, model$at, model$Rt, model$G)
    model$mts <- smoothed$mts
    model$Cts <- smoothed$Cts
  }
  if (is.numeric(pred_cred)) {
    if (0 < pred_cred & 1 > pred_cred) {
      prediction <- family$pred(model$conj_prior_param, model$outcome, parms = parms, pred_cred = pred_cred)

      model$pred <- prediction$pred
      model$var.pred <- prediction$var.pred
      model$icl.pred <- prediction$icl.pred
      model$icu.pred <- prediction$icu.pred
    }
  }

  model$log.like <- sum(family$log.like(model$conj_prior_param, model$outcome, parms = parms))

  model$m0 <- structure$m0
  model$C0 <- structure$C0
  model$names <- structure$names
  model$family <- family

  return(model)
}

#' forecast
#'
#' Makes predictions for t times ahead using the fitted model.
#'
#' @param model <undefined class> or list: The fitted model to be use for predictions.
#' @param t Numeric: Time window for prediction.
#' @param offset Matrix or scalar: offset for predictions. Should have dimensions k x t, where k is the number of outcomes of the model. If offset is not specified, the last value observed by the model will be used.
#' @param FF Array: A 3D-array containing the regression matrix for each time. It's dimension should be n x k x t, where n is the number of latent variables, k is the number of outcomes in the model. If not specified, the last value given to the model will be used.
#' @param D Array: A 3D-array containing the discount factor matrix for each time. It's dimension should be n x n x t, where n is the number of latent variables and T is the time series length. If not specified, the last value given to the model will be used in the first step, and 1 will be use thereafter.
#' @param W Array: A 3D-array containing the covariance matrix of the noise for each time. It's dimension should be n x n x t, where n is the number of latent variables and T is the time series length. If not specified, 0 will be used.
#' @param plot Bool: A flag indicating if a plot should be produced.
#' @param pred_cred Numeric: The credibility level for the I.C. intervals.
#' @param labels Vector: A string vector with the names for the series of prediction. If none is given, will use generic names.
#'
#'
#' @return A list containing:
#' \itemize{
#'    \item pred Matrix: A matrix with the predictive mean at each time. Dimensions are k x t, where k is the number of outcomes.
#'    \item var.pred Array: A 3D-array with the predictive covariance matrix at each time. Dimensions are k x k x t, where k is the number of outcomes.
#'    \item icl.pred Matrix: A matrix with the lower bound of the I.C. based on the credibility given in the arguments. Dimensions are k x t, where k is the number of outcomes.
#'    \item icu.pred Matrix: A matrix with the upper bound of the I.C. based on the credibility given in the arguments. Dimensions are k x t, where k is the number of outcomes.
#'    \item at Matrix: A matrix with the values for the linear predictor at each time. Dimensions are k x t, where k is the number of outcomes.
#'    \item Rt Array: A 3D-array with the covariance of the linear predictor matrix at each time. Dimensions are k x k x t, where k is the number of outcomes.
#'    \item plot (if so choosed): A plotly object.
#' }
#' @import ggplot2
#' @importFrom plotly ggplotly
#' @import dplyr
#' @import tidyr
#' @export
#'
#' @examples
#' T <- 200
#' w <- ((T + 20) / 40) * 2 * pi
#' y1 <- matrix(rpois((T + 20), 20 * (sin(w * 1:(T + 20) / (T + 20)) + 2)), (T + 20), 1)
#' y2 <- matrix(rpois((T + 20), 1:(T + 20) / (T + 20) + 1), (T + 20), 1)
#' y3 <- matrix(rpois((T + 20), 6), (T + 20), 1)
#' outcome <- cbind(y1, y2, y3)
#' y_pred <- outcome[T:(T + 20), ]
#'
#' outcome <- outcome[1:T, ]
#'
#' level_1 <- polynomial_block(order = 1, values = c(1, 0), by_time = FALSE)
#' level_2 <- polynomial_block(order = 2, values = c(0, 1), by_time = FALSE)
#' season_2 <- harmonic_block(period = 20, values = c(0, 1), by_time = FALSE)
#'
#'
#' fitted_data <- fit_model(level_1, level_2, season_2, outcome = outcome, family = "Multinomial")
#' show_fit(fitted_data, smooth = TRUE)$plot
#'
#' forecast(fitted_data, 20, outcome = y_pred)$plot
forecast <- function(model, t = 1, outcome = NULL, offset = NULL, FF = NULL, D = NULL, W = NULL, plot = TRUE, pred_cred = 0.95, labels = NULL) {
  n <- dim(model$mt)[1]
  t_last <- dim(model$mt)[2]
  r <- dim(model$FF)[2]
  r_out <- dim(model$outcome)[2]

  show_y <- TRUE
  if (is.null(outcome)) {
    show_y <- FALSE
    outcome <- model$outcome[t_last, ]
  }
  if (length(outcome) == 1) {
    outcome <- c(outcome, rep(0, r_out - 1))
  }
  if (length(dim(outcome)) < 2) {
    outcome <- matrix(outcome, t, r_out, byrow = TRUE)
  }

  if (length(dim(FF)) > 3) {
    stop(paste0("Error: FF should have at most 3 dimensions. Got ", length(dim(FF)), "."))
  }
  if (length(dim(D)) > 3) {
    stop(paste0("Error: D should have at most 3 dimensions. Got ", length(dim(D)), "."))
  }
  if (length(dim(W)) > 3) {
    stop(paste0("Error: W should have at most 3 dimensions. Got ", length(dim(W)), "."))
  }
  if (length(dim(offset)) > 2) {
    stop(paste0("Error: D should have at most 2 dimensions. Got ", length(dim(offset)), "."))
  }

  if (t > 10) {
    warning("Warning: Prediction window is big, results will probabily be unreliable.")
  }

  #### Consistency check ####
  if (is.null(FF)) {
    FF <- array(model$FF[, , t_last], c(n, r, t))
  }
  if (is.null(D)) {
    D <- array(model$D[, , t_last], c(n, n, t))
    D[, , -1] <- 1
  } else {
    if (all(dim(D) == 1)) {
      D <- array(D, c(n, n, t))
    } else {
      if (length(dim(D)) == 2 | (length(dim(D)) == 3 & dim(D)[3] == 1)) {
        D <- array(D, c(n, n, t))
        D[, , -1] <- 1
      }
    }
  }
  if (is.null(W)) {
    W <- array(model$W[, , t_last], c(n, n, t))
    W[, , -1] <- 0
  } else {
    if (all(dim(W) == 1)) {
      W <- array(diag(n) * W, c(n, n, t))
    } else {
      if (length(dim(W)) == 2 | (length(dim(W)) == 3 & dim(W)[3] == 1)) {
        W <- array(W, c(n, n, t))
        W[, , -1] <- 0
      }
    }
    # if(t>1){
    #   D[,,2:t]=0
    # }
    # if(t>1){
    #   W[,,2:t]=0
    # }
  }
  if (is.null(offset)) {
    offset <- model$offset[t_last, ]
  }
  if (length(dim(offset)) < 2) {
    offset <- matrix(offset, t, r_out, byrow = TRUE)
  }

  if (any(dim(FF) != c(n, r, t))) {
    stop(paste0("Error: FF has wrong dimesions. Expected ", paste(n, r, t, sep = "x"), ". Got ", paste(dim(FF), colapse = "x")))
  }
  if (any(dim(D) != c(n, n, t))) {
    stop(paste0("Error: D has wrong dimesions. Expected ", paste(n, n, t, sep = "x"), ". Got ", paste(dim(D), colapse = "x")))
  }
  if (any(dim(W) != c(n, n, t))) {
    stop(paste0("Error: W has wrong dimesions. Expected ", paste(n, n, t, sep = "x"), ". Got ", paste(dim(W), colapse = "x")))
  }
  if (any(dim(offset) != c(t, r_out))) {
    stop(paste0("Error: W has wrong dimesions. Expected ", paste(t, r_out, sep = "x"), " or a scalar. Got ", paste(dim(offset), colapse = "x")))
  }
  #####

  G <- model$G

  m0 <- model$mt[, t_last]
  C0 <- model$Ct[, , t_last]

  D <- ifelse(D == 0, 1, D)

  last_m <- m0
  last_C <- C0

  # Definindo objetos
  pred <- matrix(NA, dim(model$pred)[1], t)
  var.pred <- array(NA, c(dim(model$var.pred)[1], dim(model$var.pred)[2], t))
  icl.pred <- matrix(NA, dim(model$icl.pred)[1], t)
  icu.pred <- matrix(NA, dim(model$icu.pred)[1], t)


  for (t_i in c(1:t)) {
    next_step <- one_step_evolve(last_m, last_C, FF[, , t_i] %>% matrix(n, r), G, D[, , t_i]**0, W[, , t_i]+C0*D[, , t_i])
    last_m <- next_step$at
    last_C <- next_step$Rt
    next_step <- model$family$offset(next_step$ft, next_step$Qt, model$offset[t_i, ])
    conj_distr <- model$family$conj_prior(next_step$ft, next_step$Qt)
    prediction <- model$family$pred(conj_distr, outcome[t_i, , drop = FALSE], model$parms, pred_cred)

    pred[, t_i] <- prediction$pred
    var.pred[, , t_i] <- prediction$var.pred
    icl.pred[, t_i] <- prediction$icl.pred
    icu.pred[, t_i] <- prediction$icu.pred
  }

  return_list <- list("pred" = pred, "var.pred" = var.pred, "icl.pred" = icl.pred, "icu.pred" = icu.pred)
  if (plot) {
    r <- dim(pred)[1]
    if (is.null(labels)) {
      labels <- c("Serie_" %>% paste0(1:r))
    }

    pred <- cbind(t_last + c(1:t), t(pred) %>% as.data.frame())
    names(pred) <- c("Time", labels)
    pred <- pred %>% pivot_longer(2:(r + 1))
    names(pred) <- c("Time", "Serie", "Prediction")

    obs <- pred
    obs$Prediction <- NULL
    obs$Observation <- NA

    icl.pred <- cbind(t_last + c(1:t), t(icl.pred) %>% as.data.frame())
    names(icl.pred) <- c("Time", labels)
    icl.pred <- icl.pred %>% pivot_longer(2:(r + 1))
    names(icl.pred) <- c("Time", "Serie", "C.I.lower")

    icu.pred <- cbind(t_last + c(1:t), t(icu.pred) %>% as.data.frame())
    names(icu.pred) <- c("Time", labels)
    icu.pred <- icu.pred %>% pivot_longer(2:(r + 1))
    names(icu.pred) <- c("Time", "Serie", "C.I.upper")

    plot_data <- obs %>%
      inner_join(pred, c("Time", "Serie")) %>%
      inner_join(icl.pred, c("Time", "Serie")) %>%
      inner_join(icu.pred, c("Time", "Serie")) %>%
      mutate(Serie = as.factor(Serie))

    plot_data2 <- eval_past(model, smooth = TRUE, pred_cred = pred_cred)
    levels(plot_data2$Serie) <- labels


    max_value <- calcula_max(plot_data2$Observation - min(plot_data2$Observation))[[3]] + min(plot_data2$Observation)
    min_value <- -calcula_max(-(plot_data2$Observation - max(plot_data2$Observation)))[[3]] + max(plot_data2$Observation)

    plot_data <- rbind(plot_data2, plot_data)

    plot <- ggplotly(
      ggplot(plot_data, aes(x = Time, color = Serie, fill = Serie, linetype = ifelse(Time > t_last, "Forecast", "One-step ahead prediction"))) +
        geom_line(aes(y = Prediction)) +
        geom_ribbon(aes(ymin = C.I.lower, ymax = C.I.upper, linetype = "C.I."), alpha = 0.25) +
        geom_point(aes(y = Observation)) +
        scale_fill_hue("", na.value = NA) +
        scale_color_hue("", na.value = NA) +
        scale_linetype_manual("", values = c("solid", "dashed", "solid")) +
        scale_x_continuous("Time") +
        scale_y_continuous("$y_t$") +
        theme_bw() +
        coord_cartesian(ylim = c(min_value, max_value))
    )
    return_list$plot <- plot
  }

  return(return_list)
}

#' eval_past
#'
#' Evaluates the predictive values for the observed values used to fit the model. Predictions can be made with smoothed values or with filtered values with a time offset.
#'
#' @param model <undefined class> or list: The fitted model to be use for evaluation.
#' @param smooth bool: The flag indicating if smoothed values should be used. If TRUE, t_offset will not be used.
#' @param t_offset positive integer: The relative offset for forecast. Values for time t will be calculated based on the filtered values of time t-t_offset. Will be ignored if smooth is TRUE.
#' @param pred_cred Numeric: The credibility level for the I.C. intervals.
#'
#' @return A list containg:
#' \itemize{
#'    \item pred Matrix: A matrix with the predictive mean at each time. Dimensions are k x t, where k is the number of outcomes.
#'    \item var.pred Array: A 3D-array with the predictive covariance matrix at each time. Dimensions are k x k x t, where k is the number of outcomes.
#'    \item icl.pred Matrix: A matrix with the lower bound of the I.C. based on the credibility given in the arguments. Dimensions are k x t, where k is the number of outcomes.
#'    \item icu.pred Matrix: A matrix with the upper bound of the I.C. based on the credibility given in the arguments. Dimensions are k x t, where k is the number of outcomes.
#' }
#' @export
#'
#' @examples
#' T <- 200
#' w <- ((T + 20) / 40) * 2 * pi
#' y1 <- matrix(rpois((T + 20), 20 * (sin(w * 1:(T + 20) / (T + 20)) + 2)), (T + 20), 1)
#' y2 <- matrix(rpois((T + 20), 1:(T + 20) / (T + 20) + 1), (T + 20), 1)
#' y3 <- matrix(rpois((T + 20), 6), (T + 20), 1)
#' y <- cbind(y1, y2, y3)
#' y_pred <- y[T:(T + 20), ]
#'
#' y <- y[1:T, ]
#'
#' level <- polynomial_block(order = 1, values = 1, k = 2)
#' season_2 <- harmonic_block(period = 20, values = c(0, 1), by_time = FALSE)
#'
#'
#' fitted_data <- fit_model(level, season_2, outcome = y, family = "Multinomial")
#'
#' past <- eval_past(fitted_data, smooth = TRUE)
eval_past <- function(model, smooth = FALSE, t_offset = 0, pred_cred = 0.95, labels = NULL) {
  if (smooth & t_offset > 0) {
    t_offset <- 0
    warning("t_offset is only used if smooth is set to TRUE.")
  }
  if (t_offset < 0 | round(t_offset) != t_offset) {
    stop(paste0("ERROR: t_offset should be a positive integer. Got ", t_offset, "."))
  }
  n <- dim(model$mt)[1]
  t_last <- dim(model$mt)[2]
  r <- dim(model$FF)[2]
  k <- dim(model$outcome)[2]

  FF <- model$FF
  G <- diag(n)
  pred <- matrix(NA, k, t_last)
  var.pred <- array(NA, c(k, k, t_last))
  icl.pred <- matrix(NA, k, t_last)
  icu.pred <- matrix(NA, k, t_last)
  log.like <- rep(NA, t_last)

  if (smooth) {
    ref_mt <- model$mts
    ref_Ct <- model$Cts
    D <- model$D**0
    W <- model$W * 0
  } else {
    ref_mt <- model$mt
    ref_Ct <- model$Ct
    D <- model$D
    W <- model$W
    if (t_offset > 0) {
      # for (i in c(1:t_offset)) {
      #   G <- G %*% model$G
      # }
      G <- model$G
    }
  }

  for (i in c(1:t_last)) {
    mt <- if (i <= t_offset) {
      model$m0
    } else {
      ref_mt[, (i - t_offset):(i - t_offset)]
    }
    Ct <- if (i <= t_offset) {
      model$C0
    } else {
      ref_Ct[, , i - t_offset]
    }
    # model$family$filter(model$outcome[i, ], mt, Ct, FF[, , i] %>% matrix(n, r), G, D[, , i], W[, , i], model$offset[i, ], parms = model$parms)
    next_step <- list("at" = mt, "Rt" = Ct)
    for (t in c(1:t_offset)) {
      next_step <- one_step_evolve(next_step$at, next_step$Rt, FF[, , i] %>% matrix(n, r), G, D[, , i]**(t == 1), W[, , i])
    }
    # next_step <- one_step_evolve(mt, Ct, FF[, , i] %>% matrix(n, r), G, D[, , i], W[, , i])
    next_step <- model$family$offset(next_step$ft, next_step$Qt, model$offset[i, ])
    conj_distr <- model$family$conj_prior(next_step$ft, next_step$Qt)
    prediction <- model$family$pred(conj_distr, model$outcome[i, , drop = FALSE], model$parms, pred_cred)

    pred[, i] <- prediction$pred
    var.pred[, , i] <- prediction$var.pred
    icl.pred[, i] <- prediction$icl.pred
    icu.pred[, i] <- prediction$icu.pred
    log.like[i] <- sum(model$family$log.like(conj_distr, model$outcome[i, , drop = FALSE], model$parms))
  }
  data_list <- list("obs" = matrix(model$outcome, t_last, k), "pred" = pred %>% as.matrix() %>% t(), "icl.pred" = icl.pred %>% as.matrix() %>% t(), "icu.pred" = icu.pred %>% as.matrix() %>% t())
  data_name <- c("Observation", "Prediction", "C.I.lower", "C.I.upper")
  data_raw <- lapply(1:4, function(i) {
    data <- cbind(as.character(1:t_last) %>% as.data.frame(), data_list[[i]]) %>%
      as.data.frame() %>%
      pivot_longer(1:k + 1) %>%
      mutate(name = as.factor(name))
    names(data) <- c("Time", "Serie", data_name[i])
    levels(data$Serie) <- paste0("Serie_", 1:k)
    data
  })
  data <- do.call(function(...) {
    data_list <- list(...)
    data <- data_list[[1]]
    for (i in 2:4) {
      data <- data %>% left_join(data_list[[i]], by = c("Time", "Serie"))
    }
    data
  }, data_raw) %>% mutate(Time = as.numeric(Time), Serie = as.factor(Serie))
  if (!is.null(labels)) {
    levels(data$Serie) <- labels
  }
  return(data)
}

#' FFBS_sampling
#'
#' @param model <undefined class>: A fitted model from which to sample.
#' @param sample_size integer: The number of samples to draw.
#'
#' @return A list containing the following values:
#' \itemize{
#'    \item mt Array: An array containing the samples of latent states. Dimensions are n x T x sample_size, where n is the number of latent variable in the model and T is the number of observed values.
#'    \item ft Array: An array containing the samples of linear predictors. Dimensions are m x T x sample_size, where m is the number of linear predictors in the model and T is the number of observed values.
#'    \item param Array: An array containing the samples of the parameters of the observational model. Dimensions are k x T x sample_size, where k is the number of parameters in the observational model and T is the number of observed values.
#' }
#' @export
#'
#' @examples
#'
#' T <- 200
#' w <- (200 / 50) * 2 * pi
#' S <- exp(2 * (sin(w * 1:T / T))) * exp(-20 * (1:T / T) * ((1:T / T) - 1)) / 20
#' mu <- 20 * 1:T / T
#' outcome <- rnorm(T, mu, sqrt(1 / S))
#'
#' level <- polynomial_block(
#'   order = 2,
#'   values = c(1, 0),
#'   D = 1 / 1,
#'   by_time = FALSE
#' )
#' variance1 <- polynomial_block(
#'   order = 3,
#'   values = c(0, 1),
#'   D = 1 / 1,
#'   C0 = diag(c(1, 0.1, 0.01)),
#'   by_time = FALSE
#' )
#' variance2 <- harmonic_block(
#'   period = 50,
#'   values = c(0, 1),
#'   D = 1 / 1,
#'   by_time = FALSE
#' )
#'
#' fitted_data <- fit_model(level, variance1, variance2,
#'   outcome = outcome,
#'   family = "normal_gamma"
#' )
#'
#' show_fit(fitted_data, smooth = TRUE)$plot
#'
#' sample <- FFBS_sampling(fitted_data, 2000)
FFBS_sampling <- function(model, sample_size) {
  G <- model$G
  at <- model$at
  mt <- model$mt
  FF <- model$FF
  T_len <- dim(mt)[2]
  n <- dim(mt)[1]
  k <- dim(FF)[2]
  inv_link <- model$family$inv_link_function
  alt_chol_local <- if (n == 1) {
    sqrt
  } else {
    alt_chol
  }

  mt_sample <- rnorm(n * T_len * sample_size) %>% array(c(n, T_len, sample_size))
  ft_sample <- array(NA, c(k, T_len, sample_size))

  offset_step <- model$offset[T_len, ]
  na.flag <- any(is.null(offset_step) | any(offset_step == 0) | any(is.na(offset_step)))

  Ct_chol <- tryCatch(
    {
      chol(model$Ct[, , T_len])
    },
    error = function(e) {
      alt_chol_local(model$Ct[, , T_len])
    }
  )
  mt_sample_i <- t(Ct_chol) %*% mt_sample[, T_len, ] + mt[, T_len]
  ft_sample_i <- t(FF[, , T_len]) %*% mt_sample_i

  ft_sample_i <- model$family$offset(ft_sample_i, diag(k) * 0, if (na.flag) {
    1
  } else {
    offset_step
  })$ft
  param_sample_i <- inv_link(ft_sample_i)

  l <- dim(param_sample_i)[1]
  param_sample <- array(NA, c(l, T_len, sample_size))
  param_sample[, T_len, ] <- param_sample_i
  mt_sample[, T_len, ] <- mt_sample_i
  ft_sample[, T_len, ] <- ft_sample_i

  for (t in (T_len - 1):1) {
    offset_step <- model$offset[t, ]
    na.flag <- any(is.null(offset_step) | any(offset_step == 0) | any(is.na(offset_step)))

    Rt <- model$Rt[, , t + 1]
    Ct <- model$Ct[, , t]

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
    simple_Rt_inv <- Ct %*% t(G_now) %*% ginv(Rt)

    mts <- mt[, t] + simple_Rt_inv %*% (mt_sample_i - at[, t + 1])
    Cts <- Ct - simple_Rt_inv %*% Rt %*% t(simple_Rt_inv)
    Ct_chol <- tryCatch(
      {
        chol(Cts)
      },
      error = function(e) {
        alt_chol_local(Cts)
      }
    )
    mt_sample_i <- t(Ct_chol) %*% mt_sample[, t, ] + mts
    ft_sample_i <- t(FF[, , t]) %*% mt_sample_i

    ft_sample_i <- model$family$offset(ft_sample_i, diag(k) * 0, if (na.flag) {
      1
    } else {
      offset_step
    })$ft
    param_sample_i <- inv_link(ft_sample_i)

    param_sample[, t, ] <- param_sample_i
    mt_sample[, t, ] <- mt_sample_i
    ft_sample[, t, ] <- ft_sample_i
  }
  return(list(
    "mt" = mt_sample,
    "ft" = ft_sample,
    "param" = param_sample
  ))
}

#' @export
kernel_list <- list(
  "poisson" = poisson_kernel,
  "poisson_ord1" = poisson_kernel_ord1,
  "poisson_rn_default" = poisson_kernel_RN_default,
  "poisson_true" = poisson_kernel_true,
  "poisson_lb" = poisson_lb_kernel,
  "multinomial" = multnom_kernel,
  "normal_gamma" = normal_gamma_kernel,
  "normal_gamma_cor" = normal_gamma_cor_kernel,
  "normal" = normal_kernel,
  "gamma" = gamma_kernel,
  "gamma_lb" = gamma_lb_kernel,
  "fgamma" = Fgamma_kernel,
  "fgamma2" = Fgamma2_kernel,
  "fgamma3" = Fgamma3_kernel
)
