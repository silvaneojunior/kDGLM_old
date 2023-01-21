# usethis::use_mit_license()

# styler::style_pkg()
# devtools::document()
devtools::install('.', upgrade='never')
devtools::install_github('silvaneojunior/GDLM')

install.packages('vctrs',clean = TRUE)

remove.packages('vctrs')

#################################################################################
############################## Teste monitoramento ##############################
#################################################################################


devtools::load_all()

# library(GDLM)
library(tidyverse)
library(plotly)

T=200

# set.seed(13031998)
t_inter=150

T <- 200
w <- (200 / 50) * 2 * pi
S=10
#S=exp(-10*(1:T/T)*((1:T/T)-1))
set.seed(13031998)
outcome <- rpois(T, S+5*(1:T>t_inter))
plot(outcome)

level <- polynomial_block(
  order = 1,
  D = 1 /1,
  by_time=FALSE
)

fitted_data <- fit_model(level, outcome = outcome, family = "Poisson", p_monit=0.05)

(fitted_data$log.like.null-fitted_data$log.like.alt)[t_inter:(t_inter+4)]

show_fit(fitted_data,smooth=TRUE)$plot

plot(outcome)
lines(fitted_data$pred.model$null_model,lty=1)
lines(fitted_data$pred.model$alt_model,lty=2)

#################################################################################

devtools::load_all()
dados=read.csv('data/dados_wania_gravidade.csv',row.names = 1) %>% t %>% as.data.frame
dados=dados[1:(101-3),][1:86,]
T=dim(dados)[1]
indice_inter=54

nivel=polynomial_block(1,values=1,D=1/1,name='nivel')
# inter=polynomial_block(1,values=1:T>=indice_inter,D=1/1,C0=0,W=1:T==indice_inter,k=1,name='intervenção')

fitted_data <- fit_model(nivel,
                         outcome=dados$N_OBT %>% as.matrix,
                         offset=dados$Pac_dia %>% as.matrix, family = "Poisson",p_monit=0.05)

(fitted_data$log.like.null-fitted_data$log.like.alt)[indice_inter:(indice_inter+4)]

show_fit(fitted_data,smooth=TRUE,dynamic_plot = TRUE)$plot
plot_lat_var(fitted_data,'nivel',smooth=TRUE)$plot

plot(dados$N_OBT %>% as.matrix)
lines(fitted_data$pred.model$null_model,lty=1)
lines(fitted_data$pred.model$alt_model,lty=2)

#################################################################################

devtools::load_all()
set.seed(13031998)
T=200
y1=rpois(T,exp(5+(-T:T/T)*5))+10*((1:T)>100)
y2=rpois(T,exp(6+(-T:T/T)*5+sin((-T:T)*(2*pi/12))))
y3=rpois(T,exp(3))

y=cbind(y1,y2,y3)

level=polynomial_block(2,k=2)
season=harmonic_block(12,values=c(0,1),by_time = FALSE)

fitted_data=fit_model(level,season,outcome=y,family='Multinomial',pred_cred = 0.95, p_monit=0.05)

show_fit(fitted_data,smooth=TRUE)$plot

plot(outcome)
lines(fitted_data$pred.model$null_model,lty=1)
lines(fitted_data$pred.model$alt_model,lty=2)

#################################################################################

devtools::load_all()
T <- 200
w <- (200 / 40) * 2 * pi
phi <- 5
mu=(2 * (sin(w * 1:T / T) + 2))+5*((1:T)>100)
plot(mu)
set.seed(1331)
outcome <- matrix(rgamma(T, phi, phi / mu), T, 1)

level <- polynomial_block(order = 1, values = 1, D = 1 / 1)
season <- harmonic_block(period = 40, values = 1, D = 1 / 1)

fitted_data <- fit_model(level, season, outcome = outcome, family = "Gamma", parms = list("phi" = phi), p_monit=0.05)

show_fit(fitted_data,smooth=TRUE)$plot

plot(outcome)
lines(fitted_data$pred.model$null_model,lty=1)
lines(fitted_data$pred.model$alt_model,lty=2)

#################################################################################

devtools::load_all()
T <- 200
w <- (200 / 40) * 2 * pi
phi <- 2.5
mu=(2 * (sin(w * 1:T / T) + 2))+5*((1:T)>100)
plot(mu)
set.seed(1331)
outcome <- matrix(rgamma(T, phi, phi / mu), T, 1)

level <- polynomial_block(order = 1, values = 1, D = 1 / 1)
season <- harmonic_block(period = 40, values = 1, D = 1 / 1)

fitted_data <- fit_model(level, season, outcome = outcome, family = "Gamma", parms = list("phi" = phi), p_monit=0.05)

show_fit(fitted_data,smooth=TRUE)$plot

plot(outcome)
lines(fitted_data$pred.model$null_model,lty=1)
lines(fitted_data$pred.model$alt_model,lty=2)

#################################################################################

devtools::load_all()
T <- 200
w <- (200 / 40) * 2 * pi
phi <- 1.5
mu=(2 * (sin(w * 1:T / T) + 2))+5*((1:T)>100)
plot(mu)
set.seed(1331)
outcome <- matrix(rgamma(T, phi, phi / mu), T, 1)

level <- polynomial_block(order = 1, values = 1, D = 1 / 1)
season <- harmonic_block(period = 40, values = 1, D = 1 / 1)

fitted_data <- fit_model(level, season, outcome = outcome, family = "Gamma", parms = list("phi" = phi), p_monit=0.05)

show_fit(fitted_data,smooth=TRUE)$plot

plot(outcome)
lines(fitted_data$pred.model$null_model,lty=1)
lines(fitted_data$pred.model$alt_model,lty=2)

#################################################################################

devtools::load_all()
T <- 200
w <- (200 / 40) * 2 * pi
phi <- 0.75
mu=(2 * (sin(w * 1:T / T) + 2))+5*((1:T)>100)
plot(mu)
set.seed(1331)
outcome <- matrix(rgamma(T, phi, phi / mu), T, 1)

level <- polynomial_block(order = 1, values = 1, D = 1 / 1)
season <- harmonic_block(period = 40, values = 1, D = 1 / 1)

fitted_data <- fit_model(level, season, outcome = outcome, family = "Gamma", parms = list("phi" = phi), p_monit=0.05)

show_fit(fitted_data,smooth=TRUE)$plot

plot(outcome)
lines(fitted_data$pred.model$null_model,lty=1)
lines(fitted_data$pred.model$alt_model,lty=2)

#################################################################################

devtools::load_all()
T <- 200
w <- (200 / 40) * 2 * pi
phi <- 0.5
mu=(2 * (sin(w * 1:T / T) + 2))+5*((1:T)>100)
plot(mu)
set.seed(1331)
outcome <- matrix(rgamma(T, phi, phi / mu), T, 1)

level <- polynomial_block(order = 1, values = 1, D = 1 / 1)
season <- harmonic_block(period = 40, values = 1, D = 1 / 1)

fitted_data <- fit_model(level, season, outcome = outcome, family = "Gamma", parms = list("phi" = phi), p_monit=0.05)

show_fit(fitted_data,smooth=TRUE)$plot

plot(outcome)
lines(fitted_data$pred.model$null_model,lty=1)
lines(fitted_data$pred.model$alt_model,lty=2)

#################################################################################


devtools::load_all()

#usethis::use_mit_license()
# devtools::document()
# devtools::install('.', upgrade='never')

# library(GDLM)
library(tidyverse)
library(plotly)

T <- 200
w <- (200 / 50) * 2 * pi
S=1
#S=exp(-10*(1:T/T)*((1:T/T)-1))
mu <- 20*1:T / T + 2*((1:T)>100)
set.seed(13031998)
outcome <- rnorm(T, mu,sqrt(S))

level <- polynomial_block(
  order = 2,
  values=c(1,0),
  D = 1 / 1,
  by_time=FALSE
)
# level$D[,,]=0
# level$D[1,1,]=1
variance1 <- polynomial_block(
  order = 1,
  values=c(0,1),
  D = 1 / 1,
  by_time=FALSE
)
variance1$D[,,]=1
# variance1$D[1,1,]=1

fitted_data <- fit_model(level, variance1, outcome = outcome, family = "normal_gamma_cor",p_monit=0.05)
(fitted_data$log.like.null-fitted_data$log.like.alt)[t_inter:(t_inter+4)]

show_fit(fitted_data,smooth=TRUE)$plot

plot(outcome)
lines(fitted_data$pred.model$null_model,lty=1)
lines(fitted_data$pred.model$alt_model,lty=2)

#################################################################################


devtools::load_all()

#usethis::use_mit_license()
# devtools::document()
# devtools::install('.', upgrade='never')

# library(GDLM)
library(tidyverse)
library(plotly)

T <- 200
w <- (200 / 50) * 2 * pi
S=exp(3*(1:T/T))
#S=exp(-10*(1:T/T)*((1:T/T)-1))
mu <- 20*1:T / T + 10*((1:T)>100)
set.seed(13031998)
outcome <- rnorm(T, mu,sqrt(S))

level <- polynomial_block(
  order = 2,
  values=c(1,0),
  D = 1 / 1,
  by_time=FALSE
)
# level$D[,,]=0
# level$D[1,1,]=1
variance1 <- polynomial_block(
  order = 2,
  values=c(0,1),
  D = 1 / 1,
  C0= diag(c(1,0.1)),
  by_time=FALSE
)
# variance1$D[,,]=0
# variance1$D[1,1,]=1

fitted_data <- fit_model(level, variance1, outcome = outcome, family = "normal_gamma_cor",p_monit=0.05)
(fitted_data$log.like.null-fitted_data$log.like.alt)[t_inter:(t_inter+4)]

show_fit(fitted_data,smooth=TRUE)$plot

plot(outcome)
lines(fitted_data$pred.model$null_model,lty=1)
lines(fitted_data$pred.model$alt_model,lty=2)

#################################################################################


devtools::load_all()

#usethis::use_mit_license()
# devtools::document()
# devtools::install('.', upgrade='never')

# library(GDLM)
library(tidyverse)
library(plotly)

T <- 200
w <- (200 / 50) * 2 * pi
S=exp(1 + 1*((1:T)>100))
#S=exp(-10*(1:T/T)*((1:T/T)-1))
mu <- 20*1:T / T
set.seed(13031998)
outcome <- rnorm(T, mu,sqrt(S))

level <- polynomial_block(
  order = 2,
  values=c(1,0),
  D = 1 / 1,
  by_time=FALSE
)
# level$D[,,]=0
# level$D[1,1,]=1
variance1 <- polynomial_block(
  order = 1,
  values=c(0,1),
  D = 1 / 1,
  by_time=FALSE
)
# variance1$D[,,]=0
# variance1$D[1,1,]=1

fitted_data <- fit_model(level, variance1, outcome = outcome, family = "normal_gamma_cor",p_monit=0.05)
(fitted_data$log.like.null-fitted_data$log.like.alt)[t_inter:(t_inter+4)]

show_fit(fitted_data,smooth=TRUE)$plot

plot(outcome)
lines(fitted_data$pred.model$null_model,lty=1)
lines(fitted_data$pred.model$alt_model,lty=2)


#################################################################################


devtools::load_all()

#usethis::use_mit_license()
# devtools::document()
# devtools::install('.', upgrade='never')

# library(GDLM)
library(tidyverse)
library(plotly)

T <- 200
w <- (200 / 50) * 2 * pi
S=exp(2*(1:T/T) + 2*((1:T)>100))
#S=exp(-10*(1:T/T)*((1:T/T)-1))
mu <- 20*1:T / T
set.seed(13031998)
outcome <- rnorm(T, mu,sqrt(S))

level <- polynomial_block(
  order = 2,
  values=c(1,0),
  D = 1 / 1,
  by_time=FALSE
)
# level$D[,,]=0
# level$D[1,1,]=1
variance1 <- polynomial_block(
  order = 2,
  values=c(0,1),
  D = 1 / 1,
  C0= diag(c(1,0.1)),
  by_time=FALSE
)
# variance1$D[,,]=0
# variance1$D[1,1,]=1

fitted_data <- fit_model(level, variance1, outcome = outcome, family = "normal_gamma_cor",p_monit=0.05)
(fitted_data$log.like.null-fitted_data$log.like.alt)[t_inter:(t_inter+4)]

show_fit(fitted_data,smooth=TRUE)$plot

plot(outcome)
lines(fitted_data$pred.model$null_model,lty=1)
lines(fitted_data$pred.model$alt_model,lty=2)


#################################################################################

devtools::load_all()
library(tidyverse)
library(plotly)

T <- 200
w <- (200 / 50) * 2 * pi
S=exp(1)
#S=exp(-10*(1:T/T)*((1:T/T)-1))
mu <- (20*1:T / T)
mu <- mu+(10*1:T / T-5)*(1:T>100)
set.seed(13031998)
outcome <- rnorm(T, mu,sqrt(S))

level <- polynomial_block(
  order = 2,
  values=c(1,0),
  D = 1 / 1,
  by_time=FALSE
)
variance1 <- polynomial_block(
  order = 1,
  values=c(0,1),
  D = 1 / 1,
  by_time=FALSE
)

fitted_data <- fit_model(level, variance1, outcome = outcome, family = "normal_gamma_cor",p_monit=0.05)
(fitted_data$log.like.null-fitted_data$log.like.alt)[t_inter:(t_inter+4)]

show_fit(fitted_data)$plot

#################################################################################


devtools::load_all()

#usethis::use_mit_license()
# devtools::document()
# devtools::install('.', upgrade='never')

# library(GDLM)
library(tidyverse)
library(plotly)

T <- 200
w <- (200 / 50) * 2 * pi
S=exp(2*(sin(w * 1:T / T)))*exp(-20*(1:T/T)*((1:T/T)-1))/20
#S=exp(-10*(1:T/T)*((1:T/T)-1))
mu <- 20*1:T / T
set.seed(13031998)
outcome <- rnorm(T, mu,sqrt(S))

level <- polynomial_block(
  order = 2,
  values=c(1,0),
  D = 1 / 1,
  by_time=FALSE
)
variance1 <- polynomial_block(
  order = 3,
  values=c(0,1),
  D = 1 / 1,
  C0= diag(c(1,0.1,0.01)),
  by_time=FALSE
)
variance2 <- harmonic_block(
  period = 50,
  values = c(0, 1),
  D = 1 / 1,
  by_time=FALSE
)

fitted_data <- fit_model(level, variance1, variance2, outcome = outcome, family = "normal_gamma_cor",p_monit=0.05)
(fitted_data$log.like.null-fitted_data$log.like.alt)[t_inter:(t_inter+4)]

show_fit(fitted_data,smooth=TRUE)$plot

plot(outcome)
lines(fitted_data$pred.model$null_model,lty=1)
lines(fitted_data$pred.model$alt_model,lty=2)

#################################################################################


devtools::load_all()

#usethis::use_mit_license()
# devtools::document()
# devtools::install('.', upgrade='never')

# library(GDLM)
library(tidyverse)
library(plotly)

T <- 200
w <- (200 / 50) * 2 * pi
S=exp(2*(sin(w * 1:T / T)))*exp(-20*(1:T/T)*((1:T/T)-1))/20
#S=exp(-10*(1:T/T)*((1:T/T)-1))
t_inter=T/2
mu <- 20*1:T / T + 20*((1:T)>t_inter)
set.seed(13031998)
outcome <- rnorm(T, mu,sqrt(S))

level <- polynomial_block(
  order = 2,
  values=c(1,0),
  D = 1 / 1,
  by_time=FALSE
)
variance1 <- polynomial_block(
  order = 3,
  values=c(0,1),
  D = 1 / 1,
  C0= diag(c(1,0.1,0.01)),
  by_time=FALSE
)
variance2 <- harmonic_block(
  period = 50,
  values = c(0, 1),
  D = 1 / 1,
  by_time=FALSE
)

fitted_data <- fit_model(level, variance1, variance2, outcome = outcome, family = "normal_gamma_cor",p_monit=0.05)
(fitted_data$log.like.null-fitted_data$log.like.alt)[t_inter:(t_inter+4)]

show_fit(fitted_data,smooth=TRUE)$plot

plot(outcome)
lines(fitted_data$pred.model$null_model,lty=1)
lines(fitted_data$pred.model$alt_model,lty=2)

#################################################################################


devtools::load_all()

#usethis::use_mit_license()
# devtools::document()
# devtools::install('.', upgrade='never')

# library(GDLM)
library(tidyverse)
library(plotly)

T <- 200
t_inter=T/2
w <- (200 / 50) * 2 * pi
S=exp(2*(sin(w * 1:T / T)))*exp(-20*(1:T/T)*((1:T/T)-1))*exp(5*((1:T)>t_inter))/20
#S=exp(-10*(1:T/T)*((1:T/T)-1))
mu <- 20*1:T / T
set.seed(13031998)
outcome <- rnorm(T, mu,sqrt(S))

level <- polynomial_block(
  order = 2,
  values=c(1,0),
  D = 1 / 1,
  by_time=FALSE
)
variance1 <- polynomial_block(
  order = 3,
  values=c(0,1),
  D = 1 / 1,
  C0= diag(c(1,0.1,0.01)),
  by_time=FALSE
)
variance2 <- harmonic_block(
  period = 50,
  values = c(0, 1),
  D = 1 / 1,
  by_time=FALSE
)

fitted_data <- fit_model(level, variance1, variance2, outcome = outcome, family = "normal_gamma_cor",p_monit=0.05)
(fitted_data$log.like.null-fitted_data$log.like.alt)[t_inter:(t_inter+4)]

show_fit(fitted_data,smooth=TRUE)$plot

plot(outcome)
lines(fitted_data$pred.model$null_model,lty=1)
lines(fitted_data$pred.model$alt_model,lty=2)

