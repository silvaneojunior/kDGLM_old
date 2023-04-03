#' summary.dlm_distr
#'
#' summary method for class dlm_distr.
#'
#' @export
summary.dlm_distr <- function(dlm_distr) {
  cat(paste0(
    dlm_distr$name, " distribuition.\n\nUnkown parameters:\n",
    paste0(names(dlm_distr$var_names), ": ", dlm_distr$var_names, "\n", collapse = ""),
    if (length(dlm_distr$parms) > 0) {
      paste0(names(dlm_distr$parms), ": ", dlm_distr$parms, "\n", collapse = "\n")
    } else {
      "\n"
    },
    paste0("Serie length: ", dlm_distr$t, "\n"),
    paste0("Number of outcomes: ", dlm_distr$r)
  ))
}

#' summary.fitted_dlm
#'
#' summary method for class fitted_dlm
#'
#' @export
summary.fitted_dlm <- function(fitted_dlm, t = fitted_dlm$t, smooth = TRUE, metric_cutoff = round(fitted_dlm$t / 10)) {
  distr_names <- c(sapply(fitted_dlm$outcomes, function(x) {
    x$name
  }))
  distr_names_len <- max(nchar(names(distr_names)))
  pred_flag <- FALSE
  distr_like <- NULL
  distr_rae <- NULL
  if (is.numeric(fitted_dlm$pred_cred)) {
    if (0 < fitted_dlm$pred_cred & 1 > fitted_dlm$pred_cred) {
      distr_like <- c(sapply(fitted_dlm$outcomes, function(x) {
        # x$log.like
        sum(x$family$log.like(x$conj_prior_param, x$outcome, parms = x$parms)[-(1:metric_cutoff)], na.rm = TRUE)
      }))
      distr_rae <- c(sapply(fitted_dlm$outcomes, function(x) {
        pred <- t(x$pred)
        out <- x$outcome
        out <- ifelse(out == 0, 1, out)
        mean(abs((pred - out) / out)[-(1:metric_cutoff), ], na.rm = TRUE)
      }))
      pred_flag <- TRUE
    }
  }

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
  coef_names <- c()
  for (name in names(fitted_dlm$names)) {
    name_len <- length(fitted_dlm$names[[name]])
    if (name_len > 1) {
      coef_names <- c(coef_names, paste0(name, "_", 1:length(fitted_dlm$names[[name]])))
    } else {
      coef_names <- c(coef_names, name)
    }
  }
  len_names <- max(sapply(as.character(coef_names), function(x) {
    nchar(x)
  }))
  coef_names <- format(coef_names, width = len_names, justify = "r")

  mean_coef <- fitted_dlm[[coef_mean_name]][, t]
  std_coef <- sqrt(diag(fitted_dlm[[coef_var_name]][, , t]))
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

  if (pred_flag) {
    distr_like <- ifelse(abs(distr_like) < 0.00001,
      format(distr_like, digits = 4, width = 14, justify = "l", scientific = TRUE),
      format(round(distr_like, 5), width = 14, justify = "l", scientific = FALSE)
    )
    distr_rae <- ifelse(abs(distr_rae) < 0.00001,
      format(distr_rae, digits = 4, width = 21, justify = "l", scientific = TRUE),
      format(round(distr_rae, 5), width = 21, justify = "l", scientific = FALSE)
    )
  }

  cat(paste0(
    "Fitted DGLM with ", length(fitted_dlm$outcomes), " outcomes.\n\n",
    "Distribuitions:\n",
    paste0(names(distr_names), ": ", distr_names, "\n", collapse = ""), "\n",
    "Coeficients (", coef_label, ") at time ", t, ":\n",
    paste(format(" ", width = len_names, justify = "l"), "Estimate", "Std. Error", " t value", "Pr(>|t|)"), "\n",
    paste(coef_names, mean_coef, std_coef, t_coef, p_val_str, status, "\n", collapse = ""),
    "---\n",
    "Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1\n\n",
    "---\n",
    if (pred_flag) {
      paste0(
        format(" ", width = distr_names_len, justify = "l"), "  Pred. log-like  Relative abs. Error\n",
        paste0(format(names(distr_names), width = distr_names_len, justify = "l"), ": ", distr_like, distr_rae, "\n", collapse = ""),
        "---"
      )
    } else {
      ""
    }
  ))
}

#' plot.fitted_dlm
#'
#' plot method for class fitted_dlm
#'
#' @export
plot.fitted_dlm <- function(model, pred_cred = 0.95, smooth = TRUE, dynamic = TRUE, h = 0) {
  show_fit(model, pred_cred, smooth, dynamic, h)$plot
}

#' print.fitted_dlm
#'
#' print method for class fitted_dlm
#'
#' @export
print.fitted_dlm <- summary.fitted_dlm

#' effects.fitted_dlm
#'
#' effects method for class fitted_dlm
#'
#' @export
effects.fitted_dlm <- function(fitted_dlm, t = fitted_dlm$t, smooth = TRUE) {
  if (fitted_dlm$smooth & smooth) {
    fitted_dlm$mts[, t]
  } else {
    fitted_dlm$mt[, t]
  }
}

#' fitted.values.fitted_dlm
#'
#' fitted.values method for class fitted_dlm
#'
#' @export
fitted.values.fitted_dlm <- function(fitted_dlm, smooth = TRUE, h = 0, pred_cred = 0.95) {
  eval_past(fitted_dlm, smooth, h, pred_cred)
}

#' +.fitted_dlm
#'
#' Define add operator for class dlm_block
#'
#' @export
`+.dlm_block` <- function(e1, e2) {
  block_merge(e1, e2)
}

#' *.fitted_dlm
#'
#' Define product operator for class dlm_block
#'
#' @export
`*.dlm_block` <- function(block, k) {
  if (is.numeric(k)) {
    return(block_mult(block, k))
  } else {
    return(block_mult(k, block))
  }
}
