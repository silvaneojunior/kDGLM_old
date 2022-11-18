library(mvtnorm)
library(rootSolve)
library(tidyverse)
library(plotly)

y=0
f1=0
s=100

#### Laplace ####
f=function(x){y-exp(x)-(x-f1)/(s**2)}
ss1=multiroot(f,start=0)

m_star=ss1$root
v_star=1/(exp(m_star)+1/(s**2))

#### Taylor ####

# x*y-exp(x)-((x-f1)**2)/(s**2)

f=function(x){dpois(y,exp(x))}
f_star=ss1$root
epsilon=1e-20
d0=f(f_star)
d1=(f(f_star+epsilon)-f(f_star-epsilon))/(2*epsilon)
d2=(f(f_star+epsilon)-2*f(f_star)+f(f_star-epsilon))/(epsilon**2)

c=d0+d1*(f1-f_star)+d2*(s**2+(f1-f_star)**2)

f=function(x){x*dpois(y,exp(x))*dnorm(x,f1,s)}
ss1=multiroot(function(x){(f(x+epsilon)-f(x-epsilon))/(2*epsilon)},start=-2)
f=function(x){x*dpois(y,exp(x))}
f_star=ss1$root
epsilon=1e-20
d0=f(f_star)
d1=(f(f_star+epsilon)-f(f_star-epsilon))/(2*epsilon)
d2=(f(f_star+epsilon)-2*f(f_star)+f(f_star-epsilon))/(epsilon**2)

m2_star=(d0+d1*(f1-f_star)+d2*(s**2+(f1-f_star)**2))/c

f=function(x){(x**2)*dpois(y,exp(x))*dnorm(x,f1,s)}
ss1=multiroot(function(x){(f(x+epsilon)-f(x-epsilon))/(2*epsilon)},start=-2)
f=function(x){(x**2)*dpois(y,exp(x))}
f_star=ss1$root
epsilon=1e-20
d0=f(f_star)
d1=(f(f_star+epsilon)-f(f_star-epsilon))/(2*epsilon)
d2=(f(f_star+epsilon)-2*f(f_star)+f(f_star-epsilon))/(epsilon**2)

v2_star=(d0+d1*(f1-f_star)+d2*(s**2+(f1-f_star)**2))/c-m2_star**2

#### Taylor na verossimilhança ####

f=function(x){dpois(y,exp(x))}
f_star=ss1$root
epsilon=1e-20
d0=f(f_star)
d1=(f(f_star+epsilon)-f(f_star-epsilon))/(2*epsilon)
d2=(f(f_star+epsilon)-2*f(f_star)+f(f_star-epsilon))/(epsilon**2)

c=d0+d1*(f1-f_star)+d2*(s**2+(f1-f_star)**2)

f=function(x){x*dpois(y,exp(x))*dnorm(x,f1,s)}
ss1=multiroot(function(x){(f(x+epsilon)-f(x-epsilon))/(2*epsilon)},start=-2)
f=function(x){x*dpois(y,exp(x))}
f_star=ss1$root
epsilon=1e-20
d0=f(f_star)
d1=(f(f_star+epsilon)-f(f_star-epsilon))/(2*epsilon)
d2=(f(f_star+epsilon)-2*f(f_star)+f(f_star-epsilon))/(epsilon**2)

m2_star=(d0+d1*(f1-f_star)+d2*(s**2+(f1-f_star)**2))/c

f=function(x){(x**2)*dpois(y,exp(x))*dnorm(x,f1,s)}
ss1=multiroot(function(x){(f(x+epsilon)-f(x-epsilon))/(2*epsilon)},start=-2)
f=function(x){(x**2)*dpois(y,exp(x))}
f_star=ss1$root
epsilon=1e-20
d0=f(f_star)
d1=(f(f_star+epsilon)-f(f_star-epsilon))/(2*epsilon)
d2=(f(f_star+epsilon)-2*f(f_star)+f(f_star-epsilon))/(epsilon**2)

v2_star=(d0+d1*(f1-f_star)+d2*(s**2+(f1-f_star)**2))/c-m2_star**2

#### Monte Carlo ####

f_densi=function(x){dpois(y,exp(x))}


sample=rnorm(20000,f1,s)

c=mean(f_densi(sample))
m1=mean(sample*f_densi(sample))/c
m2=mean((sample**2)*f_densi(sample))/c
v=m2-m1**2

l=1000
x=seq(-10,10,l=l)
c_true=integrate(function(x){f_densi(x)*dnorm(x,f1,s)},-Inf,Inf)
fx=f_densi(x)*dnorm(x,f1,s)/c

plot(x,fx,type='l')
lines(x,dnorm(x,m1,sqrt(v)),lty=2,col='red')
lines(x,dnorm(x,m_star,sqrt(v_star)),lty=2,col='blue')

ggplotly(ggplot()+
  geom_line(aes(x=x,y=fx))+
  geom_line(aes(x=x,y=dnorm(x,m1,sqrt(v)),linetype='MC',color='MC'))+
    geom_line(aes(x=x,y=dnorm(x,m_star,sqrt(v_star)),linetype='laplace',color='laplace'))+
    geom_line(aes(x=x,y=dnorm(x,m2_star,sqrt(v2_star)),linetype='taylor',color='taylor'))+
  theme_bw())

print(m1)
print(m_star)
print(m2_star)
print(v)
print(v_star)
print(v2_star)
