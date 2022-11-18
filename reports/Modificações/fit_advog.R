devtools::load_all()

# usethis::use_mit_license()
devtools::document()
devtools::install('.', upgrade='never')

# library(GDLM)
library(tidyverse)
library(plotly)

T <- 200
mu <- rnorm(T, 0, 0.1)
outcome <- rnorm(T, cumsum(mu))

level <- polynomial_block(
  order = 1,
  values = c(1, 0),
  D = 1 / 0.98,
  by_time=FALSE
)
variance <- polynomial_block(
  order = 1,
  values = c(0, 1),
  D = 1 / 1,
  by_time=FALSE
)

fitted_data <- fit_model(level, variance, outcome = outcome, family = "normal_gamma")
show_fit(fitted_data, smooth = TRUE)$plot

# Poisson case
T <- 200
w <- (200 / 40) * 2 * pi
outcome <- rpois(T, exp((sin(w * 1:T / T))))

level <- polynomial_block(order = 1, values = 1, D = 1 / 0.95)
season <- harmonic_block(period = 40, values = 1, D = 1 / 0.98)

fitted_data <- fit_model(level, season, outcome = outcome, family = "Poisson")
show_fit(fitted_data, smooth = TRUE)$plot

plot(fitted_data$conj_post_param[,1]/(fitted_data$conj_post_param[,2]))
lines(exp((sin(w * 1:T / T))))

# Normal with unkown variance case
T <- 200
w <- (200 / 30) * 2 * pi
S=exp((sin(w * 1:T / T)))*10
mu <- 0
outcome <- rnorm(T, mu,sqrt(S))

level <- polynomial_block(
  order = 1,
  values = c(1, 0),
  D = 1 / 1,
  by_time=FALSE
)
variance1 <- polynomial_block(
  order = 1,
  values = c(0, 1),
  D = 1 / 1,
  by_time=FALSE
)
variance2 <- harmonic_block(
  period = 30,
  values = c(0, 1),
  D = 1 / 1,
  by_time=FALSE
)

fitted_data <- fit_model(level, variance1, variance2, outcome = outcome, family = "normal_gamma")
show_fit(fitted_data, smooth = TRUE)$plot

alpha=rep(NA,T)
beta=rep(NA,T)
for(i in 1:T){
  ft=t(fitted_data$FF[,,i]) %*% fitted_data$mts[,i]
  Qt=t(fitted_data$FF[,,i]) %*% fitted_data$Cts[,,i] %*% fitted_data$FF[,,i]

  alpha[i]=1/(-3+3*sqrt(1+2*Qt[2,2]/3))
  beta[i]=alpha[i]*exp(-ft[2]-Qt[2,2]/2)
}

plot(beta/(alpha-1),type='l',ylim=c(0,50))
lines(1/qgamma(0.025,alpha,beta),lty=2)
lines(1/qgamma(0.975,alpha,beta),lty=2)
points(S)

plot(fitted_data$conj_post_param[,4]/(fitted_data$conj_post_param[,3]-1),type='l')
lines(1/qgamma(0.025,fitted_data$conj_post_param[,3],fitted_data$conj_post_param[,4]),lty=2)
lines(1/qgamma(0.975,fitted_data$conj_post_param[,3],fitted_data$conj_post_param[,4]),lty=2)
points(S)
