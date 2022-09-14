devtools::load_all()

library(tidyverse)
library(plotly)


T=200
phi=2.5
mu=5
set.seed(13031998)
y= matrix(rgamma(T,phi,phi/mu),T,1)

level=polynomial_block(order=1,values=1,D=1/1)

fitted_data=fit_model(level,outcome=y,kernel='Gamma',parms=list('phi'=phi))

true_alpha=phi*(fitted_data$tau0[1]+(1:T))
true_beta=phi*(fitted_data$tau1[1]+cumsum(y))
mt=true_beta/(true_alpha-1)

#plot(y[,1])
plot(mt)
lines(exp(fitted_data$mt+fitted_data$Qt/2) %>% t)

plot(true_beta/(true_alpha-1))
lines((fitted_data$tau1_star/(fitted_data$tau0_star-1)) %>% t)

plot((true_beta**2)/((true_alpha**2)*(true_alpha-2)))
lines(((fitted_data$tau1_star**2)/((fitted_data$tau0_star**2)*(fitted_data$tau0_star-2))) %>% t)

plot(true_beta)
lines(fitted_data$tau1_star %>% t)

plot((fitted_data$tau0_star %>% t)/phi)

ft=0
qt=1

tau0=1/qt
tau0*exp(ft-qt/2) -> tau1

ft_star <- log(tau1_star)-digamma(tau0_star)
# qt_star <- 2*(log(tau0_star/tau1_star)+ft_star)


0=-digamma(tau0)+log(tau0)-qt/2
log(tau0)+ft-qt/2=log(tau1)
