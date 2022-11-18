library(tidyverse)
library(rootSolve)
library(lamW)

source('R\\full_gamma_kernel.R')

f1=0
f2=0
Q1=1
Q2=1
Q12=-0.2

f1=f1
f2=f2
q1=Q1
q2=Q2
q12=Q12

# f1=ft[1,]
# f2=ft[1,]-ft[2,]
# q1=Qt[1,1]
# q2=Qt[1,1]+Qt[2,2]+2*Qt[1,2]
# q12=Qt[1,1]-Qt[1,2]

ft=matrix(c(f1,f2),2,1)
Qt=matrix(c(Q1,Q12,Q12,Q2),2,2)

conj_prior=convert_FGamma_Normal(ft,Qt)

conj_norm=convert_Normal_FGamma(conj_prior)
ft=conj_norm$ft
Qt=conj_norm$Qt
conj_norm

x=seq(0,10,l=1000)
plot(x,dlnorm(x,f1,sqrt(q1)),type='l')
# lines(x,dlnorm(x,ft[1],sqrt(Qt[1,1])),type='l')
lines(x,dgamma(x,(conj_prior$n+1)/2,conj_prior$k-conj_prior$n+conj_prior$k*log(conj_prior$tau/conj_prior$k)-conj_prior$theta),lty=2)
# lines(x,dgamma(x,(conj_prior$n+2)/2,conj_prior$n*log(conj_prior$tau/conj_prior$n)-conj_prior$theta),lty=2)


plot(x,dlnorm(x,f2,sqrt(q2)),type='l')
lines(x,dlnorm(x,ft[2],sqrt(Qt[2,2])),type='l')

# n0=1
# tau0=log(2)
# theta0=0
#
# conj_prior=list(n0=n0,tau0=tau0,theta0=theta0)

conj_prior=update_FGamma(conj_prior,rgamma(1,1))
conj_prior=update_FGamma(conj_prior,rgamma(1,1))
conj_prior=update_FGamma(conj_prior,rgamma(1,1))
conj_prior=update_FGamma(conj_prior,rgamma(1,1))
conj_prior=update_FGamma(conj_prior,rgamma(1,1))

x=seq(0,10,l=1000)
print(conj_prior)
n=conj_prior$n
k=conj_prior$k
tau=conj_prior$tau
theta=conj_prior$theta
a=(n+1)/2
b=(k-n+k*log(tau/k)-theta)
# print(trigamma(a))
plot(x,dgamma(x,a,b),type='l',col=rainbow(3,s=0.5)[1])

for(i in 2:2){
  norm_prior=convert_Normal_FGamma(conj_prior)
  conj_prior=convert_FGamma_Normal(norm_prior$ft,norm_prior$Qt)
  print(conj_prior)
  n=conj_prior$n
  tau=conj_prior$tau
  theta=conj_prior$theta
  k=conj_prior$k
  a=(n+1)/2
  b=(k-n+k*log(tau/k)-theta)
  # print(trigamma(a))
  lines(x,dgamma(x,a,b),col=rainbow(3,s=0.5)[i])
}
# n=conj_prior$n
# tau=conj_prior$tau
# theta=conj_prior$theta
# k=n

y=rgamma(5,10,10/3)
y=rgamma(20,1,1)
x=rgamma(20000,1,1)

n=length(y)
tau=sum(y)
theta=sum(log(y))
k=n

# n=1.2774475
# tau=0.8251163
# theta=-1.951927
# k=n

f_densi=function(x){exp(-k*lgamma(x)+theta*x+lgamma(n*x+1)-(n*x)*log(tau))}
c_val2=mean(f_densi(x)/dgamma(x,1,1))
# c_val2
c_val=integrate(f_densi,0,Inf)$value

x=seq(0,100,l=1000)
plot(x,f_densi(x)/c_val,type='l')
lines(x,dgamma(x,(n+3)/2,n*log(tau/n)-theta))

f_densi=function(x){x*exp(-k*lgamma(x)+theta*x+lgamma(n*x+1)-(n*x)*log(tau))}
mean(f_densi(x)/dgamma(x,1,1))/c_val2
integrate(f_densi,0,Inf)$value/c_val
a/b

f_densi=function(x){((n*x+1)/(tau))*exp(-k*lgamma(x)+theta*x+lgamma(n*x+1)-(n*x+1)*log(tau))}
mean(f_densi(x)/dgamma(x,1,1))/c_val2
integrate(f_densi,0,Inf)$value/c_val

f_densi=function(x){(x*(digamma(n*x+1)-log(tau))-lgamma(x))*exp(-k*lgamma(x)+theta*x+lgamma(n*x+1)-(n*x+1)*log(tau))}
mean(f_densi(x)/dgamma(x,1,1))/c_val2
integrate(f_densi,0,Inf)$value/c_val

x=a/b
# x*(digamma(n*x+1)-log(tau))-lgamma(x)
# (digamma(n*x+1)-log(tau))+x*n*trigamma(n*x+1)-digamma(x)
d2=n*trigamma(n*x+1)+n*trigamma(n*x+1)+x*(n^2)*psigamma(n*x+1,2)-trigamma(x)

(x*(digamma(n*x+1)-log(tau))-lgamma(x))+d2*(a/(b**2))

f_densi=function(x){log(x)*exp(-k*lgamma(x)+theta*x+lgamma(n*x+1)-(n*x+1)*log(tau))}
f1_star=mean(f_densi(x)/dgamma(x,1,1))/c_val2
f1=integrate(f_densi,0,Inf)$value/c_val

f_densi=function(x){(digamma(n*x+1)-log(tau))*exp(-k*lgamma(x)+theta*x+lgamma(n*x+1)-(n*x+1)*log(tau))}
f2_star=mean(f_densi(x)/dgamma(x,1,1))/c_val2
f2=integrate(f_densi,0,Inf)$value/c_val

f_densi=function(x){((digamma(n*x+1)-log(tau))**2)*exp(-k*lgamma(x)+theta*x+lgamma(n*x+1)-(n*x+1)*log(tau))}
mean(f_densi(x)/dgamma(x,1,1))/c_val2-f2_star**2
integrate(f_densi,0,Inf)$value/c_val-f2**2

f1=ft[1,]
f2=ft[2,]
q1=Qt[1,1]
q2=Qt[2,2]
q12=Qt[1,2]

Hq1=exp(f1+q1/2)
Hq2=exp(f2+q2/2)
Hq3=(f2+q12)*Hq1
# x=rnorm(20000,f1,sqrt(q1))
# Hq4=mean(lgamma(exp(x)))

# Hq4
Hq4=lgamma(exp(f1+q1/2))+(trigamma(exp(f1))*exp(2*f1)+digamma(exp(f1))*exp(f1))*q1

# digamma(exp(x))*exp(x)
# trigamma(exp(x))*exp(2*x)+digamma(exp(x))*exp(x)

parms=list(
  'Hq1'=Hq1,
  'Hq2'=Hq2,
  'Hq3'=Hq3,
  'Hq4'=Hq4)
Hq1
Hq2
(parms$Hq3-parms$Hq4)
