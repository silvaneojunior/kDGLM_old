#' calcula_max
#'
#' Auxiliary function to calculate the axis limits and gradation for plots.
#'
#' @param pre_max Numeric: A vector/matrix from which to calculate the axis limits and gradation.
#'
#' @return A list contaning the gradation for the axis, the number of ticks in the axis and the maximum value.
calcula_max <- function(pre_max) {
  if (length(pre_max) == 0 | sum(pre_max**2) < 10**-20) {
    pre_max <- 1
  } else {
    pre_max <- max(pre_max)
  }
  scaled_max <- log10(pre_max)
  category <- scaled_max %% 1
  value <- 10**(floor(log10(max(pre_max))))
  if (category < 0.1) {
    value <- value / 10
  } else {
    if (category < 0.25) {
      value <- value / 5
    } else {
      if (category < 0.5) {
        value <- value / 2
      }
    }
  }
  interval_size <- (pre_max %/% value) + 2
  max_value <- value * interval_size

  return(list(value, interval_size, max_value))
}

#' show_fit
#'
#' Calculate the preditive mean and some quantiles for the observed data and show a plot.
#'
#' @param model <undefined class> or list: A fitted GDLM model.
#' @param pred_cred Numeric: The credibility value for the credibility interval.
#' @param smooth Bool: A flag indicating if the smoothed should be used. If false, the filtered distribuition will be used.
#' @param dynamic Bool: A flag indicating if the created plot should be dynamic.
#' @param t_offset Integer: A integer with the amount of steps ahead should be used for prediciton. Only used if smooth is false.
#' @param labels Character: A vector containg the names for each time series.
#'
#' @return A list containg:
#' \itemize{
#'    \item data tibble object: A data frame containing the observations, predictions and credibility intervals at each time.
#'    \item plot ggplot or plotly object: A plot showing the predictive mean and credibility interval with the observed data.
#' }
#' @export
#'
#' @import ggplot2
#' @importFrom plotly ggplotly
#' @import dplyr
#' @import tidyr
#'
#' @examples
#' T <- 200
#' w <- (200 / 40) * 2 * pi
#' outcome <- rpois(T, 20 * (sin(w * 1:T / T) + 2))
#'
#' level <- polynomial_block(order = 1, values = 1, D = 1 / 0.95)
#' season <- harmonic_block(period = 40, values = 1, D = 1 / 0.98)
#'
#' fitted_data <- fit_model(level, season, outcome = outcome, family = "Poisson")
#' show_fit(fitted_data, smooth = TRUE)$plot
show_fit <- function(model, pred_cred = 0.95, smooth = TRUE, dynamic_plot = TRUE, t_offset = 0, labels = NULL) {
  n <- dim(model$mt)[1]
  t_last <- dim(model$mt)[2]
  eval <- eval_past(model, smooth = smooth, t_offset = t_offset, pred_cred=pred_cred)

  max_value <- calcula_max(model$outcome - min(model$outcome))[[3]] + min(model$outcome)
  min_value <- -calcula_max(-(model$outcome - max(model$outcome)))[[3]] + max(model$outcome)

  plt <- ggplot(eval, aes(x = Time, fill = Serie, color = Serie)) +
    geom_ribbon(aes(ymin = C.I.lower, ymax = C.I.upper, linetype = "Fitted values"), alpha = 0.25) +
    geom_line(aes(y = Prediction, linetype = "Fitted values")) +
    geom_point(aes(y = Observation), alpha = 0.5) +
    scale_fill_hue("", na.value = NA) +
    scale_color_hue("", na.value = NA) +
    scale_y_continuous(name = "$y_t$") +
    scale_x_continuous("Time") +
    theme_bw() +
    coord_cartesian(ylim = c(min_value, max_value))
  if (dynamic_plot) {
    plt <- ggplotly(plt)
  }
  return(list('data'=eval,'plot'=plt))
}

#' plot_lat_var
#'
#' @param model <undefined class> or list: A fitted GDLM model.
#' @param var Character: The name of the variables to plot (same value passed while creating the structure). Any variable whose name partially match this variable will be ploted.
#' @param smooth Bool: A flag indicating if the smoothed distribuition should be used. If false, the filtered distribution shall be used.
#' @param cut_off Integer: The number of initial steps that should be skipped in the plot. Usually, the model is still learning in the initial steps, so the estimated values are not realiable.
#' @param pred_cred Numeric: The credibility value for the credibility interval.
#' @param dynamic Bool: A flag indicating if the created plot should be dynamic.
#'
#' @return A list containg:
#' \itemize{
#'    \item plot ggplot or plotly object: A plot showing the predictive mean and credibility interval with the observed data.
#'    \item data tibble: A data frame containg the data used in the plot.
#' }
#' @export
#'
#' @import ggplot2
#' @importFrom plotly ggplotly
#' @import dplyr
#' @import tidyr
#'
#' @examples
#' T <- 200
#' w <- (200 / 40) * 2 * pi
#' outcome <- rpois(T, 20 * (sin(w * 1:T / T) + 2))
#'
#' level <- polynomial_block(order = 1, values = 1, D = 1 / 0.95,name='level_effect')
#' season <- harmonic_block(period = 40, values = 1, D = 1 / 0.98,name='season_effect')
#'
#' fitted_data <- fit_model(level, season, outcome = outcome, family = "Poisson")
#' plot_lat_var(fitted_data,'effect', smooth = TRUE)$plot
plot_lat_var <- function(model, var, smooth = TRUE, cut_off = 10, pred_cred = 0.95, dynamic = TRUE) {
  if (!any(grepl(var, names(model$names)))) {
    stop(paste0("Error: Invalid selected variable. Got ", var, ", expected one of the following:\n", names(model$names)))
  }
  if (pred_cred >= 1 | pred_cred <= 0) {
    stop(paste0("Error: Invalid value for I.C. width. Must be between 0 and 1, got ", pred_cred))
  }

  indice <- c()
  names_var <- c()
  for (i in names(model$names)) {
    if (grepl(var, i)) {
      indice <- c(indice, model$names[[i]])
      count <- 0
      for (index in model$names[[i]]) {
        count <- count + 1
        names_var <- c(names_var, paste0(i, ":latent value ", count))
      }
    }
  }
  size <- length(indice)
  t <- dim(model$mts)[2]
  m1 <- if (smooth) {
    model$mts[indice, ]
  } else {
    model$mt[indice, ]
  }
  m1 <- m1 %>%
    matrix(size, t) %>%
    t()
  std_mat <- if (smooth) {
    model$Cts[indice, indice, ]
  } else {
    model$Ct[indice, indice, ]
  }
  if (size > 1) {
    std_mat <- std_mat %>% apply(3, diag)
  }
  std_mat <- std_mat %>%
    sqrt() %>%
    matrix(size, t) %>%
    t()

  lim_i <- m1 + qnorm((1 - pred_cred) / 2) * std_mat
  lim_s <- m1 + qnorm(1 - (1 - pred_cred) / 2) * std_mat

  m1 <- as.data.frame(m1)
  lim_i <- as.data.frame(lim_i)
  lim_s <- as.data.frame(lim_s)

  names(m1) <- names_var
  names(lim_i) <- names_var
  names(lim_s) <- names_var

  max_value <- calcula_max(m1 - min(m1))[[3]] + min(m1)
  min_value <- -calcula_max(-(m1 - max(m1)))[[3]] + max(m1)

  m1$time <- c(1:dim(m1)[1])
  lim_i$time <- c(1:dim(lim_i)[1])
  lim_s$time <- c(1:dim(lim_s)[1])

  m1 <- m1[-c(1:cut_off), ] %>%
    pivot_longer(1:size) %>%
    rename(media = value)
  lim_i <- lim_i[-c(1:cut_off), ] %>%
    pivot_longer(1:size) %>%
    rename(lim_i = value)
  lim_s <- lim_s[-c(1:cut_off), ] %>%
    pivot_longer(1:size) %>%
    rename(lim_s = value)

  IC_label <- paste0("I.C. (", pred_cred * 100 %>% round(), "%)")

  plot_data <- m1 %>%
    inner_join(lim_i, by = c("time", "name")) %>%
    inner_join(lim_s, by = c("time", "name"))

  n_var <- length(unique(plot_data$name))
  color_list <- rainbow(n_var, s = 0.5)
  names(color_list) <- paste(unique(plot_data$name), "point estimate")

  fill_list <- rainbow(n_var, s = 0.5)
  names(fill_list) <- paste(unique(plot_data$name), IC_label)

  plt <- ggplot(plot_data, aes(x = time, fill = paste(name, "point estimate"), color = paste(name, "point estimate"))) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_x_continuous("Time") +
    scale_color_manual("", values = color_list, na.value = NA) +
    scale_fill_manual("", values = fill_list, na.value = NA) +
    labs(title = paste0(var, " (", ifelse(smooth, "smoothed", "only filtered"), ")")) +
    scale_y_continuous("Parameter value") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90)) +
    geom_ribbon(aes(ymin = lim_i, ymax = lim_s, fill = paste(name, IC_label), color = paste(name, IC_label)), alpha = 0.25) +
    geom_line(aes(y = media)) +
    coord_cartesian(ylim = c(min_value, max_value))
  if (dynamic) {
    plt <- ggplotly(plt)
  }
  return(list("plot" = plt, "data" = plot_data))
}
