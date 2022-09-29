devtools::load_all()

#usethis::use_mit_license()
# devtools::document()
# devtools::install('.', upgrade='never')

# library(GDLM)
library(tidyverse)
library(plotly)

T=200

set.seed(13031998)

mu=0
s=1

y=rnorm(T,mu,s)

level=polynomial_block(values=1,order=1,D=1/1,k=2,C0=1)

resultado=fit_model(level,
                    outcome = y,
                    family='normal_gamma',
                    pred_cred = 0.95,
                    smooth_flag =TRUE)

show_fit(resultado,smooth=TRUE)$plot

mu0=resultado$conj_prior_param$mu0[1]
c0=resultado$conj_prior_param$c0[1]
alpha=resultado$conj_prior_param$alpha[1]
beta=resultado$conj_prior_param$beta[1]

n=1:T
obs_y=y
y_mean=cumsum(obs_y)/n
y_mean2=cumsum(obs_y**2)/n

s=(cumsum(obs_y**2)-(n*(y_mean)**2))/n

mu0_star <- (c0 * mu0 + n*y_mean) / (c0 + n)
c0_star <- c0 + n
alpha_star <- alpha + n/2
beta_star <- beta + 0.5 * (n*s+((c0*n*(mu0 - y_mean)**2) / (c0 + n)))

plot(obs_y)
lines(mu0_star,col='blue',main='m0')
lines(resultado$conj_post_param$mu0,col='red')
legend(x=0,y=max(y),c('observações','ajuste','valor correto'),pch=c(1,NA,NA),lty=c(0,1,1),col=c('black','red','blue'))

plot(c0_star,col='blue',type='l',main='c0')
lines(resultado$conj_post_param$c0)
legend(x=0,y=max(c0_star),c('ajuste','valor correto'),lty=c(1,1),col=c('black','blue'))

plot(alpha_star,col='blue',type='l',main='alpha')
lines(resultado$conj_post_param$alpha)
legend(x=0,y=max(alpha_star),c('ajuste','valor correto'),lty=c(1,1),col=c('black','blue'))


plot(beta_star,col='blue',type='l',main='beta')
lines(resultado$conj_post_param$beta)
legend(x=0,y=max(beta_star),c('ajuste','valor correto'),lty=c(1,1),col=c('black','blue'))


plot(alpha_star/beta_star,type='l',col='blue')
lines(resultado$conj_post_param$alpha/resultado$conj_post_param$beta)
legend(x=0,y=max(c0_star),c('ajuste','valor correto'),lty=c(1,1),col=c('black','blue'))


plot(beta_star/((alpha_star-1))*(1+1/c0_star),type='l',col='blue')
lines(resultado$conj_post_param$beta/((resultado$conj_post_param$alpha-1))*(1+1/resultado$conj_post_param$c0))

show_fit(resultado)$plot

#################################################################################
############################## Teste de composição ##############################
#################################################################################

a=function(ft,qt){1/(-3+3*sqrt(1+2*qt/3))}
b=function(ft,qt){a(ft,qt)*exp(-ft - 0.5 * qt)}


f=function(a.post,b.post){digamma(a.post) - log(b.post)}
q=function(a.post,b.post){trigamma(a.post)}

inter_f=function(ft,qt){f(a(ft,qt),b(ft,qt))}
inter_q=function(ft,qt){q(a(ft,qt),b(ft,qt))}

inter_a=function(at,bt){a(f(at,bt),q(at,bt))}
inter_b=function(at,bt){b(f(at,bt),q(at,bt))}

ft_init=1
qt_init=2

at_init=10
bt_init=1

ft=ft_init
qt=qt_init

at=at_init
bt=bt_init

plot_ft=c(ft_init)
plot_qt=c(qt_init)

plot_at=c(at_init)
plot_bt=c(bt_init)

t_final=1

for(i in 1:t_final){
  fts=inter_f(ft,qt)
  qts=inter_q(ft,qt)
  ft=fts
  qt=qts

  plot_ft=c(plot_ft,fts)
  plot_qt=c(plot_qt,qts)

  ats=inter_a(at,bt)
  bts=inter_b(at,bt)
  at=ats
  bt=bts

  plot_at=c(plot_at,ats)
  plot_bt=c(plot_bt,bts)
}

plot_at[2]/plot_at[1]

plot(0:t_final,plot_ft,type='l')
plot(0:t_final,plot_qt,type='l')

plot(0:t_final,plot_at,type='l')
plot(0:t_final,plot_bt,type='l')



ft_init=0
qt_init=1

at_init=1
bt_init=1

ft=ft_init
qt=qt_init

at=at_init
bt=bt_init

plot_ft_in=seq(-5,5,l=200)
plot_qt_in=seq(-5,100,l=200)

plot_at_in=seq(0,5,l=200)
plot_bt_in=seq(-5,5,l=200)

t_final=length(plot_ft_in)

plot_ft_out=rep(NA,t_final)
plot_qt_out=rep(NA,t_final)

plot_at_out=rep(NA,t_final)
plot_bt_out=rep(NA,t_final)

for(i in 1:t_final){
  plot_ft_out[i]=inter_f(plot_ft_in[i],qt_init)
  plot_qt_out[i]=inter_q(ft_init,plot_qt_in[i] %>% exp)
  plot_at_out[i]=inter_f(plot_at_in[i] %>% exp,bt_init)
  plot_bt_out[i]=inter_b(at_init,plot_bt_in[i] %>% exp)
}
plot(plot_ft_in,plot_ft_out-plot_ft_in,type='l')
plot(plot_qt_in,plot_qt_out/exp(plot_qt_in),type='l')

plot(plot_at_in,plot_at_out/exp(plot_at_in),type='l')
plot(plot_bt_in,plot_bt_out/exp(plot_bt_in),type='l')

################################################################################
################################# Caso Poisson #################################
################################################################################

devtools::load_all()

library(tidyverse)
library(plotly)

T=200

lambda=1

set.seed(13031998)
y=rpois(T,lambda)


level=polynomial_block(order=1,values=1,D=1/1)

resultado=fit_model(level,
                    outcome = y,
                    family='poisson')

show_fit(resultado,smooth=TRUE)$plot

a=resultado$conj_prior_param[1,1]
b=resultado$conj_prior_param[1,2]

a.post=a+cumsum(y)
b.post=b+1:T

plot(a.post,col='blue',type='l',main='alpha')
lines(resultado$conj_post_param[,1])

plot(b.post,col='blue',type='l',main='beta')
lines(resultado$conj_post_param[,2])

################################################################################
################################## Caso Gamma ##################################
################################################################################

devtools::load_all()

library(tidyverse)
library(plotly)

T=200

phi=1
mu=1

set.seed(13031998)
y=rgamma(T,phi,phi/mu)


level=polynomial_block(order=1,values=1,D=1/1)

resultado=fit_model(level,
                    outcome = y,
                    family='gamma',
                    parms=list('phi'=phi))

show_fit(resultado,smooth=TRUE)$plot

a=resultado$conj_prior_param[1,1]
b=resultado$conj_prior_param[1,2]

a.post=a+1:T*phi
b.post=b+cumsum(y)*phi

plot(a.post,col='blue',type='l',main='alpha')
lines(resultado$conj_post_param[,1])

plot(b.post,col='blue',type='l',main='beta')
lines(resultado$conj_post_param[,2])

################################################################################
############################### Caso Multinomial ###############################
################################################################################

devtools::load_all()

library(tidyverse)
library(plotly)

T=10

lambda1=10
lambda2=10
lambda3=10

set.seed(13031998)
y1=rpois(T,lambda1)
y2=rpois(T,lambda2)
y3=rpois(T,lambda3)
y=cbind(y1,y2,y3)


level=polynomial_block(order=1,values=1,D=1/1,k=2)

resultado=fit_model(level,
                    outcome = y,
                    family='multinomial')

show_fit(resultado,smooth=TRUE)$plot

l1=resultado$conj_prior_param[1,1]
l2=resultado$conj_prior_param[1,2]
l3=resultado$conj_prior_param[1,3]

l1.post=l1+cumsum(y1)
l2.post=l2+cumsum(y2)
l3.post=l3+cumsum(y3)

plot(l1.post,col='blue',type='l',main='lambda1')
lines(resultado$conj_post_param[,1])

plot(l2.post,col='blue',type='l',main='lambda2')
lines(resultado$conj_post_param[,2])

plot(l3.post,col='blue',type='l',main='lambda3')
lines(resultado$conj_post_param[,3])

#################################################################################
############################## Teste de composição ##############################
#################################################################################

alpha=rep(1,3)*10000

t_final=1

for(i in 1:t_final){
  alpha=do.call(convert_Dir_Normal,convert_Normal_Dir(alpha))
}
print(alpha)

################################################################################
############################## Ajuste Multinomial ##############################
################################################################################
devtools::load_all()

library(tidyverse)
library(plotly)

T=200
y1=rpois(T,exp(5+(-T:T/T)*5))
y2=rpois(T,exp(6+(-T:T/T)*5+sin((-T:T)*(2*pi/12))))
y3=rpois(T,exp(5))

y=cbind(y1,y2,y3)

level=polynomial_block(2,k=2)
season=harmonic_block(12,values=c(0,1),by_time = FALSE)

fitted_data=fit_model(level,season,outcome=y,family='Multinomial',pred_cred = 0.95)

show_fit(fitted_data,smooth = TRUE)$plot

#################################################################################
###################### Caso Normal com variância conhecida ######################
#################################################################################

devtools::load_all()

#usethis::use_mit_license()
devtools::document()
devtools::install('.', upgrade='never')

# library(GDLM)
library(tidyverse)
library(plotly)

T=200

set.seed(13031998)

s=1/sqrt(rgamma(1,1,1/5))
mu=rnorm(1,0,2*s)

y1=rnorm(T,mu,s)
y2=rnorm(T,mu,2*s)
y=cbind(y1,y2)

level=polynomial_block(values=1,order=1,D=1/1,k=2,C0=1)

resultado=fit_model(level,
                    outcome = y,
                    family='normal',
                    pred_cred =  0.95,
                    smooth_flag =TRUE,
                    parms=list('Sigma'=diag(c(1,1))*(s**2)))

eval_past(resultado)
show_fit(resultado,smooth=TRUE)$plot

ft=resultado$conj_prior_param$ft[1]
Qt=resultado$conj_prior_param$Qt[1]

n=1:T
y_mean=cumsum(y)/n

ft_star <- (ft*Qt + n*y_mean/(s**2)) / (1/Qt + n/((s**2)))
Qt_star <- 1/ (1/Qt + n/((s**2)))

plot(y)
lines(ft_star,col='blue',main='ft')
lines(resultado$conj_post_param$ft,col='red')
legend(x=0,y=max(y),c('observações','ajuste','valor correto'),pch=c(1,NA,NA),lty=c(0,1,1),col=c('black','red','blue'))

plot(Qt_star,col='blue',type='l',main='Qt')
lines(resultado$conj_post_param$Qt)
legend(x=0,y=max(Qt_star),c('ajuste','valor correto'),lty=c(1,1),col=c('black','blue'))

