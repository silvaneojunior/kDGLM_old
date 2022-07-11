source('R/main.R')
source('R/kernels.R')
source('R/plot_helper.R')
source('R/structure_helper.R')
library(dplyr)
library(ggplot2)
library(plotly)
library(tidyr)
library(MASS)
library(rootSolve)
library(extraDistr)

T=200
phi=0.5
y= matrix(rgamma(T,phi,phi/5),T,1)

level=polynomial_block(order=1,values=matrix(1,1,T),D=1/1)

fitted_data=fit_model(level,data_out=y,kernel='Gamma',parms=list('phi'=phi))

true_alpha=1.16188+phi*(1:T)
true_beta=0.5083452+phi*cumsum(y)
mt=true_beta/(true_alpha-1)

#plot(y[,1])
plot(mt)
lines(exp(fitted_data$mt[1,]))

