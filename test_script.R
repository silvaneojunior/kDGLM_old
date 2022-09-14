devtools::load_all()

library(tidyverse)
library(plotly)

T=200

mu=1:T/20
s=0.2*(1:T)

mu0=0
C0=1
# offset da média dos dados observado.
offset=1

set.seed(13031998)
y=rnorm(T,mu,sqrt(1/s))+offset
y=cbind(y,0)

level=polynomial_block(order=2,values=1,D=1/0.95,k=2)
# level$C0[2,2]=0.001
level$C0[3,3]=0.0001
level$C0[4,4]=0.0001

resultado=fit_model(level,
                    outcome = y,
                    kernel='normal_gamma')

show_fit(resultado,smooth=TRUE)$plot


set.seed(13031998)

mu=0
s=1

y=rnorm(T,mu,s)
y=cbind(y,0)

level=polynomial_block(order=1,values=1,D=1/1,k=2,C0=1)

resultado=fit_model(level,
                    outcome = y,
                    kernel='normal_gamma')

show_fit(resultado,smooth=TRUE)$plot

mu0=resultado$mu0[1]
c0=resultado$c0[1]
alpha=resultado$alpha[1]
beta=resultado$beta[1]

n=1:T
obs_y=y[,1]
y_mean=cumsum(obs_y)/n
y_mean2=cumsum(obs_y**2)/n

s=(cumsum(obs_y**2)-(n*(y_mean)**2))/n
# s=(y_mean2-y_mean**2)/n

mu0_star <- (c0 * mu0 + n*y_mean) / (c0 + n)
c0_star <- c0 + n
alpha_star <- alpha + n/2
beta_star <- beta + 0.5 * (n*s+((c0*n*(mu0 - y_mean)**2) / (c0 + n)))

plot(obs_y)
lines(mu0_star,col='blue',main='m0')
lines(resultado$mu0_star,col='red')

plot(c0_star,col='blue',type='l',main='c0')
lines(resultado$c0_star)

plot(alpha_star,col='blue',type='l',main='alpha')
lines(resultado$alpha_star)

plot(beta_star,col='blue',type='l',main='beta')
lines(resultado$beta_star)

plot(alpha_star/beta_star,type='l',col='blue')
lines(resultado$alpha_star/resultado$beta_star)

plot(beta_star/((alpha_star-1))*(1+1/c0_star),type='l',col='blue')
lines(resultado$beta_star/((resultado$alpha_star-1))*(1+1/resultado$c0_star))

show_fit(resultado)$plot

#usethis::use_mit_license()
devtools::document()
devtools::install('.', upgrade='never')


#devtools::install_github('silvaneojunior/GDLM')

library(rootSolve)

f1 <- 0
f2 <- 0
q1 <- 1
q2 <- 0.0001
q12 <- 0

1/tau1+(tau0**2)*tau2/tau3
exp(f2+q2/2)*(q1+(f1+q12)**2)

system <- function(x, parms) {
    y=exp(x)
    out <- digamma(y)-x + parms$q2 / 2
  return(out)
}

s <- multiroot(
  f = system,
  start = c(0),
  parms = list("q2" = q2)
)


tau0 <- f1+q12
tau1 <- 1/(exp(f2 + q2 / 2)*q1)
helper=-3+3*sqrt(1+2*q2/3)
tau2=exp(s$root) #1/helper
tau3 <- tau2/exp(f2 + q2 / 2)

#tau1=tau1*tau2/(tau2-1)
tau2=tau2

tau2/tau3-exp(f2+q2/2)
digamma(tau2)-log(tau3)-f2

print(tau0)
print(tau1)
print(tau2)
print(tau3)

tau0->tau0_star
tau1->tau1_star
tau2->tau2_star
tau3->tau3_star

f1star <- tau0_star
f2star <- digamma(tau2_star)-log(tau3_star)
q1star <- tau3_star / (tau1_star * (tau2_star - 1))
q2star <- trigamma(tau2_star)

print(f1star)
print(f2star)
print(q1star)
print(q2star)

f1 <- f1star
f2 <- f2star
q1 <- q1star
q2 <- q2star
q12 <- 0

# 1/tau1+(tau0**2)*tau2/tau3
# exp(f2+q2/2)*(q1+(f1+q12)**2)

system <- function(x, parms) {
  y=exp(x)
  out <- digamma(y)-x + parms$q2 / 2
  return(out)
}

s <- multiroot(
  f = system,
  start = c(0),
  parms = list("q2" = q2)
)


tau0 <- f1+q12
tau1 <- 1/(exp(f2 + q2 / 2)*q1)
helper=-3+3*sqrt(1+2*q2/3)
tau2=exp(s$root) #1/helper
tau3 <- tau2/exp(f2 + q2 / 2)

#tau1=tau1*(tau2-1)/tau2

#digamma(tau2)-log(tau3)
print(tau0_star)
print(tau1_star)
print(tau2_star)
print(tau3_star)
print(tau0)
print(tau1)
print(tau2)
print(tau3)


T <- 200
w <- ((T + 20) / 40) * 2 * pi
y1 <- matrix(rpois((T + 20), 20 * (sin(w * 1:(T + 20) / (T + 20)) + 2)), (T + 20), 1)
y2 <- matrix(rpois((T + 20), 1:(T + 20) / (T + 20) + 1), (T + 20), 1)
y3 <- matrix(rpois((T + 20), 6), (T + 20), 1)
y <- cbind(y1, y2, y3)
y_pred <- y[T:(T + 20), ]

y <- y[1:T, ]

level <- polynomial_block(order = 1, values = 1,k=2)
season_2 <- harmonic_block(period = 20, values = c(0, 1),by_time = FALSE)


fitted_data <- fit_model(level, season_2, outcome = y, kernel = "Multinomial")

show_fit(fitted_data,smooth=TRUE)$plot


###############################################

a=function(ft,qt){1/(-3+3*sqrt(1+2*qt/3))}
b=function(ft,qt){a(ft,qt)*exp(-ft - 0.5 * qt)}


f=function(a.post,b.post){digamma(a.post) - log(b.post)}
q=function(a.post,b.post){trigamma(a.post)}

inter_f=function(ft,qt){f(a(ft,qt),b(ft,qt))}
inter_q=function(ft,qt){q(a(ft,qt),b(ft,qt))}

inter_a=function(at,bt){a(f(at,bt),q(at,bt))}
inter_b=function(at,bt){b(f(at,bt),q(at,bt))}

ft_init=0
qt_init=1

at_init=5
bt_init=5

ft=ft_init
qt=qt_init

at=at_init
bt=bt_init

plot_ft=c(ft_init)
plot_qt=c(qt_init)

plot_at=c(at_init)
plot_bt=c(bt_init)

t_final=5

for(i in 1:t_final){
  fts=inter_f(ft,qt)
  qts=inter_q(ft,qt)
  ft=fts
  qt=qts

  plot_ft=c(plot_ft,fts)
  plot_qt=c(plot_qt,qts)

  ats=inter_f(at,bt)
  bts=inter_q(at,bt)
  at=ats
  bt=bts

  plot_at=c(plot_at,ats)
  plot_bt=c(plot_bt,bts)
}
plot(0:t_final,plot_ft,type='l')
plot(0:t_final,plot_qt,type='l')

plot(0:t_final,plot_at,type='l')
plot(0:t_final,plot_bt,type='l')

###############################################

devtools::load_all()

library(tidyverse)
library(plotly)

T=200

lambda=200

set.seed(13031998)
y=rpois(T,lambda)


level=polynomial_block(order=1,values=1,D=1/1)

resultado=fit_model(level,
                    outcome = y,
                    kernel='poisson')

show_fit(resultado,smooth=TRUE)$plot

a=resultado$a[1]
b=resultado$b[1]

a.post=a+cumsum(y)
b.post=b+1:T

plot(a.post,col='blue',type='l',main='alpha')
lines(resultado$a.post[1,])

plot(b.post,col='blue',type='l',main='beta')
lines(resultado$b.post[1,])

