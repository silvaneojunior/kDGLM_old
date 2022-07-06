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

library(GDLM)
T=200
w=(200/40)*2*pi
phi=2.5
y= matrix(rgamma(T,phi,phi/20*(sin(w*1:T/T)+2)),T,1)

level=polynomial_block(order=1,values=matrix(1,1,T),D=1/0.95)
season=harmonic_block(period=40,values=matrix(1,1,T),D=1/0.98)

final_block=block_join(level,season)

fitted_data=fit_model(final_block,data_out=y,kernel='Gamma',parms=list('phi'=phi))
show_fit(fitted_data,smooth = TRUE)$plot

#usethis::use_mit_license()
devtools::document()
devtools::install('.', upgrade='never')


devtools::install_github('silvaneojunior/GDLM')
