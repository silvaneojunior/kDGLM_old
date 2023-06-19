usethis::use_gpl3_license()
devtools::document()
styler::style_pkg()
devtools::build_manual()

usethis::use_vignette('fit_normal','Fitting a model for Normal outcomes')

usethis::use_vignette('fit_normal','Fitting a model for Normal outcomes')

devtools::document()
devtools::load_all()
help(Gamma)

# devtools::install('C:\\Jupyter\\Mestrado\\Pacote\\kDGLM', upgrade='never')
devtools::install(upgrade='never')
devtools::install_github('silvaneojunior/kDGLM', upgrade='never')

remove.packages('kDGLM')


devtools::load_all()

# Poisson case
T <- 200
w <- (200 / 40) * 2 * pi
data <- rpois(T, 20 * (sin(w * 1:T / T) + 2))

level <- polynomial_block(rate = 1, D = 0.95)*2
season <- harmonic_block(rate = 1, period = 40, D = 1 / 0.98)*2

outcome1 <- Poisson(lambda = "rate_1", outcome = data)
outcome2 <- Poisson_alt(lambda = "rate_2", outcome = data)

fitted_data <- fit_model(level, season,
                         outcomes = list('Original'=outcome1,
                                         'Alternative'=outcome2
                                         )
                         )
summary(fitted_data)

show_fit(fitted_data, smooth = TRUE)$plot
forecast(fitted_data, 20, outcome = list('Original'=rpois(20,50),'Alternative'=rpois(20,50)))$plot
dlm_sampling(fitted_data, 2000)

##################################################################

devtools::load_all()
# Multinomial case
T <- 200
y1 <- rpois(T, exp(5 + (-T:T / T) * 5))
y2 <- rpois(T, exp(6 + (-T:T / T) * 5 + sin((-T:T) * (2 * pi / 12))))
y3 <- rpois(T, exp(5))

y <- cbind(y1, y2, y3)

level <- (polynomial_block(p1 = 1,order=2) + polynomial_block(p2 = 1,order=2))*2
season <- harmonic_block(p2 = 1, period = 12)*2
outcome1 <- Multinom(p = c("p1_1", "p2_1"), outcome = y)
outcome2 <- Multinom_alt(p = c("p1_2", "p2_2"), outcome = y)

fitted_data <- fit_model(level, season, outcomes = list('Original'=outcome1,'Alternative'=outcome2))
summary(fitted_data)

show_fit(fitted_data, smooth = TRUE, plotly=TRUE)$plot

forecast(fitted_data, 20,
         outcome = list(
           'Original'=cbind(rpois(20,50),rpois(20,50),rpois(20,50)),
           'Alternative'=cbind(rpois(20,50),rpois(20,50),rpois(20,50)))
         # plot=TRUE
)$plot
dlm_sampling(fitted_data, 2000)

##################################################################

devtools::load_all()
# Multinomial case
# T <- 200
# y1 <- rbinom(T,1,0.8)
#
# y <- cbind(y1, 1-y1)

level <- polynomial_block(p1 = 1,order=1)*2
outcome1 <- Multinom(p = c("p1_1"), outcome = y)
outcome2 <- Multinom_alt(p = c("p1_2"), outcome = y)

fitted_data <- fit_model(level, outcomes = list('Original'=outcome, 'Alternative'=outcome2),pred_cred=NA)
summary(fitted_data)

#
# show_fit(fitted_data, smooth = TRUE, plotly=TRUE)$plot
#
# forecast(fitted_data, 20,
#          outcome = list(
#            'Original'=cbind(rpois(20,50),rpois(20,50),rpois(20,50)),
#            'Alternative'=cbind(rpois(20,50),rpois(20,50),rpois(20,50)))
#          # plot=TRUE
# )$plot
# dlm_sampling(fitted_data, 2000)
#
#
# smp=dlm_sampling(fitted_data, 2000)

index=1
plot(cumsum(y[,index])+fitted_data$outcomes$Original$conj_prior_param[1,index],type='l',lty=2)
lines(fitted_data$outcomes$Original$conj_prior_param[,index],lty=1,col='red')
lines(fitted_data$outcomes$Alternative$conj_prior_param[,index],lty=1,col='blue')


##################################################################

devtools::load_all()

# Normal case
T <- 200
mu <- rnorm(T, 0, 0.1)
data <- rnorm(T, cumsum(mu))

level <- polynomial_block(
  mu = 1,
  D = 1 / 0.95
)*2
variance <- polynomial_block(
  sigma2 = 1
)*2

# Known variance
outcome1 <- Normal(mu = "mu_1", sigma2 = 1, outcome = data)
outcome2 <- Normal_alt(mu = "mu_2", sigma2 = 1, outcome = data)

fitted_data <- fit_model(level, outcomes = list('Original'=outcome1,'Alternative'=outcome2))
summary(fitted_data)

show_fit(fitted_data, smooth = TRUE)$plot
forecast(fitted_data, 20)$plot
dlm_sampling(fitted_data, 2000)

# Unknown variance
outcome1 <- Normal(mu = "mu_1", sigma2 = "sigma2_1", outcome = data)
outcome2 <- Normal_alt(mu = "mu_2", sigma2 = "sigma2_2", outcome = data)

fitted_data <- fit_model(level, variance, outcomes = list('Original'=outcome1,'Alternative'=outcome2))
summary(fitted_data)

show_fit(fitted_data, smooth = TRUE)$plot
forecast(fitted_data, 20)$plot
dlm_sampling(fitted_data, 2000)


##################################################################

devtools::load_all()

# Gamma case
set.seed(13031998)
T <- 200
w <- (200 / 40) * 2 * pi
phi <- 2.5
data <- matrix(rgamma(T, phi, phi / (20 * (sin(w * 1:T / T) + 2))), T, 1)

level <- polynomial_block(mu = 1, D = 1 / 0.95)*2
season <- harmonic_block(mu = 1, period = 40, D = 1 / 0.98)*2
scale <- polynomial_block(phi = 1, D = 1 / 1)*2

# Known shape
outcome1 <- Gamma(phi = phi, mu = "mu_1", outcome = data)
outcome2 <- Gamma_alt(phi = phi, mu = "mu_2", outcome = data)

fitted_data <- fit_model(level, season, outcomes = list('Original'=outcome1,'Alternative'=outcome2))
summary(fitted_data)

show_fit(fitted_data, smooth = TRUE)$plot
forecast(fitted_data, 20)$plot
dlm_sampling(fitted_data, 2000)

# Unknown shape
outcome1 <- Gamma(phi = "phi_1", mu = "mu_1", outcome = data)
outcome2 <- Gamma_alt(phi = "phi_2", mu = "mu_2", outcome = data)

fitted_data <- fit_model(level, season, scale,
                         outcomes = list(
                           # 'Original'=outcome1,
                           'Alternative'=outcome2
                           )
                         )
summary(fitted_data)

show_fit(fitted_data, smooth = TRUE)$plot
forecast(fitted_data, 20)$plot

##################################################################
library(tidyverse)
devtools::load_all()

r=5
k=r+r*(r+1)/2
T=500

# set.seed(1531531)
# # A=cbind(rnorm(r*5,0,1) %>% matrix(r,5),rnorm(r*(r-5),0,0.01) %>% matrix(r,r-5))
A=rnorm(r*r,0,1) %>% matrix(r,r)
C=A%*%t(A)
ein_C=eigen(C)
var=diag(r)
diag(var)=1/rgamma(r,5,5)
# diag(var)=c(10,1,0.1)
C=eigen(C)$vectors%*%var%*%t(eigen(C)$vectors)

# rho=0.9
# C=matrix(c(1,rho,-0.2,rho,1,-0.2*rho,-0.2,-0.2*rho,1),3,3)
# D=diag(c(1,1,1))
# C=sqrt(D)%*%C%*%sqrt(D)
# C=C#[1:r,1:r]


# A=cbind(rnorm(r*5,0,1) %>% matrix(r,5),rnorm(r*(r-5),0,0.01) %>% matrix(r,r-5))
# A=rnorm(r*r,0,1) %>% matrix(r,r)
# C=A%*%t(A)
# D=diag(c(10,1,0.1))
# Q=eigen(C)$vector
# C=t(Q)%*%D%*%Q
# eigen(C)
# diag(1/sqrt(diag(C)))%*%C%*%diag(1/sqrt(diag(C)))

mu=rpois(r,1)
data=rnorm(r*T) %>% matrix(T,r)
data=data%*%chol(C)+matrix(mu,T,r,byrow=TRUE)

# y=data

devtools::load_all()
level=polynomial_block(mu=1,name='media')*r
var=polynomial_block(sigma=1,name='var')*(r*(r+1)/2)
Sigma_var=matrix(NA,r,r)
Sigma_var[lower.tri(Sigma_var,diag=TRUE)]=var$var_names

y=data[1:500,]
y=data[,order(diag(var(y)))]
t_i=dim(y)[1]

outcome=Normal(mu=level$var_names,
               Tau=Sigma_var,
               outcome=y,
               alt_method = TRUE)

A=fit_model(level,var,outcomes=outcome,pred_cred=-1)

view=show_fit(A)
view$plot
data_view=view$data
mean(data_view$Observation>data_view$C.I.upper | data_view$Observation<data_view$C.I.lower)

smp=dlm_sampling(A,2000)

mu_index=outcome$parms$mu_index
var_index=outcome$parms$var_index
cor_index=outcome$parms$cor_index
x=smp$param$Serie_1[,t_i,] %>% rowMeans
rho=matrix(1,r,r)
rho[upper.index] <- rho[lower.index] <- x[cor_index]
var=diag(sqrt(1/x[var_index]))
Sigma=var%*%rho%*%var

var_Data=var(y)

diag(var_Data)
diag(Sigma)

var_Data
Sigma

image(var_Data,zlim=c(min(var_Data,Sigma),max(var_Data,Sigma)),col=hcl.colors(200, "RdYlGn", rev = TRUE))
image(Sigma,zlim=c(min(var_Data,Sigma),max(var_Data,Sigma)),col=hcl.colors(200, "RdYlGn", rev = TRUE))


cor_Data=cor(data)
cor_Data
rho

image(cor_Data,zlim=c(-1,1),col=hcl.colors(200, "RdYlGn", rev = TRUE))
image(rho,zlim=c(-1,1),col=hcl.colors(200, "RdYlGn", rev = TRUE))

  #######################
