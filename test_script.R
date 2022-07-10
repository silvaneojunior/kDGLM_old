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

#library(GDLM)

# Multinomial case
T=200
w=(200/40)*2*pi
phi=2.5
y= matrix(rgamma(T,phi,phi/20*(sin(w*1:T/T)+2)),T,1)

level=polynomial_block(order=1,values=matrix(1,1,T),D=1/0.95)
season=harmonic_block(period=40,values=matrix(1,1,T),D=1/0.98)

final_block=block_join(level,season)

fitted_data=fit_model(final_block,data_out=y,kernel='Gamma',parms=list('phi'=phi))

#show_fit(fitted_data,smooth = TRUE)$plot

predict(fitted_data,20)$plot

##################################
T=200
w=((T+20)/40)*2*pi
y1= matrix(rpois((T+20),20*(sin(w*1:(T+20)/(T+20))+2)),(T+20),1)
y2= matrix(rpois((T+20),1:(T+20)/(T+20)+1),(T+20),1)
y3= matrix(rpois((T+20),6),(T+20),1)
y=cbind(y1,y2,y3)
y_pred=y[T:(T+20),]

y=y[1:T,]

level_1=polynomial_block(order=1,values=matrix(c(rep(1,T),rep(0,T)),2,T,byrow=TRUE))
level_2=polynomial_block(order=2,values=matrix(c(rep(0,T),rep(1,T)),2,T,byrow=TRUE))
season_2=harmonic_block(period=20,values=matrix(c(rep(0,T),rep(1,T)),2,T,byrow=TRUE))


fitted_data=fit_model(level_1,level_2,season_2,data_out=y,kernel='Multinomial')
show_fit(fitted_data,smooth = TRUE)$plot

predict(fitted_data,20,y=y_pred)$plot

#usethis::use_mit_license()
devtools::document()
devtools::install('.', upgrade='never')


devtools::install_github('silvaneojunior/GDLM')
