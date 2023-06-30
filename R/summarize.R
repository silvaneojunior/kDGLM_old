#' Summary for a fitted kDGLM model
#'
#' Prints a report for a fitted_dlm object.
#'
#' @param fitted_dlm A fitted_dlm object.
#' @param t Integer: The time index for the latent states.
#' @param smooth Bool: A flag indicating if the smoothed distribution should be used for the latent states.
#' @param metric_cutoff Integer: The cutoff time index for the metric calculation. Values before that time will be ignored.
#'
#' @export
#' @importFrom stats pnorm
#'
#' @examples
#'
#' T <- 200
#' w <- (200 / 40) * 2 * pi
#' data <- rpois(T, 20 * (sin(w * 1:T / T) + 2))
#'
#' level <- polynomial_block(rate = 1, D = 0.95)
#' season <- harmonic_block(rate = 1, period = 40, D = 0.98)
#'
#' outcome <- Poisson(lambda = "rate", outcome = data)
#'
#' fitted_data <- fit_model(level, season, outcomes = outcome)
#' report_dlm(fitted_data)
#'
#' @family {auxiliary visualization functions for the fitted_dlm class}
report_dlm <- function(fitted_dlm, t = fitted_dlm$t, smooth = fitted_dlm$smooth, metric_cutoff = round(fitted_dlm$t / 10)) {
  r <- length(fitted_dlm$outcomes)
  k <- dim(fitted_dlm$mt)[1]
  T <- dim(fitted_dlm$mt)[2]
  distr_names <- list()
  distr_like <- rep(NA, r)
  distr_rae <- rep(NA, r)
  for (outcome_index in 1:r) {
    outcome <- fitted_dlm$outcomes[[outcome_index]]
    distr_names[names(fitted_dlm$outcomes)[outcome_index]] <- outcome$name

    in_range <- (1:T) > metric_cutoff
    prediction <- outcome$calc_pred(outcome$conj_prior_param[in_range, ], outcome$outcome[in_range, ], parms = outcome$parms, pred_cred = 0.95)
    distr_like[outcome_index] <- sum(prediction$log.like, na.rm = TRUE)

    pred <- t(prediction$pred)
    out <- outcome$outcome[in_range, ]
    out_div <- ifelse(out == 0, 1, out)
    distr_rae[outcome_index] <- mean(abs((pred - out) / out_div), na.rm = TRUE)
  }
  distr_names_len <- max(sapply(names(distr_names), nchar))

  coef_label <- if (fitted_dlm$smooth & smooth) {
    "smoothed"
  } else {
    "filtered"
  }
  coef_mean_name <- if (fitted_dlm$smooth & smooth) {
    "mts"
  } else {
    "mt"
  }
  coef_var_name <- if (fitted_dlm$smooth & smooth) {
    "Cts"
  } else {
    "Ct"
  }
  coef_names <- rep(NA, k)
  for (name in names(fitted_dlm$var_names)) {
    name_len <- length(fitted_dlm$var_names[[name]])
    name_i <- name
    if (name_len > 1) {
      name_i <- paste0(name, "_", 1:length(fitted_dlm$var_names[[name]]))
    }
    coef_names[fitted_dlm$var_names[[name]]] <- name_i
  }
  len_names <- max(sapply(as.character(coef_names), function(x) {
    nchar(x)
  }))
  coef_names <- format(coef_names, width = len_names, justify = "r")

  mean_coef <- fitted_dlm[[coef_mean_name]][, t]
  var_mat <- fitted_dlm[[coef_var_name]][, , t]
  if (length(var_mat) == 1) {
    std_coef <- sqrt(abs(var_mat))
  } else {
    std_coef <- sqrt(abs(diag(var_mat)))
  }
  t_coef <- mean_coef / std_coef
  p_val <- 2 * (1 - pnorm(abs(mean_coef) / std_coef))
  status <- rep(" ", length(coef_names))
  status[p_val <= 0.01] <- "."
  status[p_val <= 0.05] <- "*"
  status[p_val <= 0.01] <- "**"
  status[p_val <= 0.001] <- "***"

  mean_coef <- format(round(mean_coef, 5), width = 8, justify = "l")
  std_coef <- format(round(std_coef, 5), width = 10, justify = "l")
  t_coef <- format(round(t_coef, 5), width = 7, justify = "l")
  p_val_str <- ifelse(p_val < 0.001,
    format(p_val, digits = 3, justify = "l", scientific = TRUE),
    format(round(p_val, 3), width = 8, justify = "l", scientific = FALSE)
  )
  p_val_str <- ifelse(p_val < 1e-12,
    "  <1e-12",
    p_val_str
  )

  distr_like <- ifelse(abs(distr_like) < 0.00001,
    format(distr_like, digits = 4, width = 14, justify = "l", scientific = TRUE),
    format(round(distr_like, 5), width = 14, justify = "l", scientific = FALSE)
  )
  distr_rae <- ifelse(abs(distr_rae) < 0.00001,
    format(distr_rae, digits = 4, width = 21, justify = "l", scientific = TRUE),
    format(round(distr_rae, 5), width = 21, justify = "l", scientific = FALSE)
  )

  cat(paste0(
    "Fitted DGLM with ", length(fitted_dlm$outcomes), " outcomes.\n\n",
    "distributions:\n",
    paste0(names(distr_names), ": ", distr_names, "\n", collapse = ""), "\n",
    "Coeficients (", coef_label, ") at time ", t, ":\n",
    paste(format(" ", width = len_names, justify = "l"), "Estimate", "Std. Error", " t value", "Pr(>|t|)"), "\n",
    paste(coef_names, mean_coef, std_coef, t_coef, p_val_str, status, "\n", collapse = ""),
    "---\n",
    "Signif. codes:  0 \xe2\x80\x98***\xe2\x80\x99 0.001 \xe2\x80\x98**\xe2\x80\x99 0.01 \xe2\x80\x98*\xe2\x80\x99 0.05 \xe2\x80\x98.\xe2\x80\x99 0.1 \xe2\x80\x98 \xe2\x80\x99 1\n\n",
    "---\n",
    format(" ", width = distr_names_len, justify = "l"), "  Pred. log-like  Relative abs. Error\n",
    paste0(format(names(distr_names), width = distr_names_len, justify = "l"), ": ", distr_like, distr_rae, "\n", collapse = ""),
    "---"
  ))
}

#' Summary for a kDGLM outcome
#'
#' Prints a report for a dlm_distr object.
#'
#' @param dlm_distr A fitted_dlm object.
#'
#' @export
#' @keywords internal
#' @family {auxiliary visualization functions for the fitted_dlm class}
report_distr <- function(dlm_distr) {
  cat(paste0(
    dlm_distr$name, " distribution.\n\nUnkown parameters:\n",
    paste0(names(dlm_distr$pred_names), ": ", dlm_distr$pred_names, "\n", collapse = ""),
    if (length(dlm_distr$parms) > 0) {
      paste0(names(dlm_distr$parms), ": ", dlm_distr$parms, "\n", collapse = "\n")
    } else {
      "\n"
    },
    paste0("Serie length: ", dlm_distr$t, "\n"),
    paste0("Number of outcomes: ", dlm_distr$r)
  ))
}

#' Summary for a searched_dlm object
#'
#' Prints a report for a searched_dlm object.
#'
#' @param searched_dlm A searched_dlm object.
#'
#' @export
#' @keywords internal
#' @family {auxiliary visualization functions for the fitted_dlm class}
report_searched_dlm <- function(searched_dlm) {
  print(searched_dlm$search.data[1:5, ])
  report_dlm(searched_dlm$model)
}
