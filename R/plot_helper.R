#' show_fit
#'
#' Calculate the preditive mean and some quantiles for the observed data and show a plot.
#'
#' @param model fitted_dlm: A fitted DGLM.
#' @param pred_cred Numeric: The credibility value for the credibility interval.
#' @param smooth Bool: A flag indicating if the smoothed should be used. If false, the filtered distribuition will be used.
#' @param dynamic Bool: A flag indicating if the created plot should be dynamic.
#' @param h Integer: A integer with the amount of steps ahead should be used for prediciton. Only used if smooth is false.
#'
#' @return A list containg:
#' \itemize{
#'    \item data tibble object: A data frame containing the observations, predictions and credibility intervals at each time.
#'    \item plot ggplot or plotly object: A plot showing the predictive mean and credibility interval with the observed data.
#' }
#' @export
#'
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#'
#' @examples
#' T <- 200
#' w <- (200 / 40) * 2 * pi
#' data <- rpois(T, 20 * (sin(w * 1:T / T) + 2))
#'
#' level <- polynomial_block(rate = 1, D = 1 / 0.95)
#' season <- harmonic_block(rate = 1, period = 40, D = 1 / 0.98)
#'
#' outcome <- Poisson(lambda = "rate", outcome = data)
#'
#' fitted_data <- fit_model(level, season, outcomes = outcome)
#' summary(fitted_data)
#'
#' show_fit(fitted_data, smooth = TRUE)$plot
show_fit <- function(model, pred_cred = 0.95, smooth = TRUE, dynamic = TRUE, h = 0) {
  n <- dim(model$mt)[1]
  t_last <- dim(model$mt)[2]
  eval <- eval_past(model, smooth = smooth, h = h, pred_cred = pred_cred)

  obs_na_rm <- eval$Observation[!is.na(eval$Observation)]
  max_value <- calcula_max(obs_na_rm - min(obs_na_rm))[[3]] + min(obs_na_rm)
  min_value <- -calcula_max(-(obs_na_rm - max(obs_na_rm)))[[3]] + max(obs_na_rm)

  n_colors <- length(unique(eval$Serie))
  colors <- rainbow(n_colors, s = 0.6)
  names(colors) <- unique(eval$Serie)

  plt <- ggplot(eval, aes(x = Time, fill = Serie, color = Serie)) +
    geom_ribbon(aes(ymin = C.I.lower, ymax = C.I.upper, linetype = "Fitted values"), alpha = 0.25) +
    geom_line(aes(y = Prediction, linetype = "Fitted values")) +
    geom_point(aes(y = Observation, linetype = "Observations"), alpha = 0.5) +
    scale_fill_manual("", na.value = NA, values = colors) +
    scale_color_manual("", na.value = NA, values = colors) +
    scale_y_continuous(name = "$y_t$") +
    scale_x_continuous("Time") +
    theme_bw() +
    coord_cartesian(ylim = c(min_value, max_value))

  if (any(model$outcomes[[1]]$alt.flags == 1)) {
    colors[["Detected changes"]] <- "black"
    plt <- plt +
      geom_vline(data = data.frame(xintercept = (1:t_last)[model$outcomes[[1]]$alt.flags == 1], linetype = "Detected changes"), aes(xintercept = xintercept, linetype = linetype, fill = linetype, color = linetype)) +
      scale_linetype_manual("", values = c("dashed", "solid", "solid")) +
      scale_color_manual("", na.value = NA, values = colors)
  }
  if (dynamic) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("The plotly package is required for dynamic plots")
    } else {
      plt <- plotly::ggplotly(plt)
    }
  }
  return(list("data" = eval, "plot" = plt))
}

#' plot_lat_var
#'
#' @param model <undefined class> or list: A fitted DGLM model.
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
#' @import dplyr
#' @import tidyr
#'
#' @examples
#'
#' T <- 200
#' w <- (200 / 40) * 2 * pi
#' data <- rpois(T, 20 * (sin(w * 1:T / T) + 2))
#'
#' level <- polynomial_block(rate = 1, D = 1 / 0.95, name = "level_effect")
#' season <- harmonic_block(rate = 1, period = 40, D = 1 / 0.98, name = "season_effect")
#'
#' outcome <- Poisson(lambda = "rate", outcome = data)
#'
#' fitted_data <- fit_model(level, season, outcomes = outcome)
#' plot_lat_var(fitted_data, "effect", smooth = TRUE)$plot
plot_lat_var <- function(model, var = "", smooth = TRUE, cut_off = 10, pred_cred = 0.95, dynamic = TRUE) {
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
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("The plotly package is required for dynamic plots.")
    } else {
      plt <- plotly::ggplotly(plt)
    }
  }
  return(list("plot" = plt, "data" = plot_data))
}
