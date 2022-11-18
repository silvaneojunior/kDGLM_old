devtools::load_all()

# Gamma case
set.seed(13031998)
T <- 200
a=20
b=20/5
outcome <- matrix(rgamma(T,a, b), T, 1)

level <- polynomial_block(order = 1, values = 1, D = 1 / 1,k=2)

fitted_data <- fit_model(level, outcome = outcome, family = "FGamma")

sample=FFBS_sampling(fitted_data,10000)
params=rowMeans(sample$param[,T,])
mu=mean(sample$param[1,T,]/sample$param[2,T,])
ic_up=mean(qgamma(0.975,sample$param[1,T,]/sample$param[2,T,]))
ic_down=mean(qgamma(0.025,sample$param[1,T,]/sample$param[2,T,]))
print(mu)

plot(outcome)
lines(1:T,rep(mu,T))
lines(1:T,rep(ic_up,T),lty=2)
lines(1:T,rep(ic_down,T),lty=2)

mean(sample$param[2,T,])
quantile(sample$param[2,T,],0.975)
quantile(sample$param[2,T,],0.025)


mean(sample$param[1,T,])
quantile(sample$param[1,T,],0.975)
quantile(sample$param[1,T,],0.025)


show_fit(fitted_data, smooth = TRUE)$plot



x=seq(0.5,1000,l=1000)
fx1=rep(NA,1000)
fx2=rep(NA,1000)
for(i in 1:1000){
  tau=x[i]
  # tau=1
  n=(parms$Hq2*tau-1)/parms$Hq1#+(1.5/parms$Hq1)
  theta=n*log(tau/n)-(n+3)/(2*parms$Hq1)
  k=n


  # f_densi=function(x){exp(-k*lgamma(x)+theta*x+lgamma(n*x+1)-(n*x+1)*log(tau))}
  # c_val=integrate(f_densi,0,Inf)$value
  #
  # f_m=function(x){(x*digamma(n*x+1)-lgamma(x))*f_densi(x)}
  # Hp3=integrate(f_m,0,Inf)$value/c_val-log(tau)*parms$Hq1
  a=(n+3)/2
  b=(n*log(tau/n)-theta)
  y=rgamma(20000,a,b)
  f_densi=function(x){(x*digamma(n*x+1)-lgamma(x))}
  Hp3=mean(f_densi(y))
  fx1[i]=Hp3
  # f_densi=function(x){
  #   x*log(n)-(1/(2*n)+1/(12*(n**2)))/x+x+log(x)/2-1/(12*x)-11/12
  #   }
  # Hp3=mean(f_densi(y))
  Hp3=(a/b)*log(n)-(1/(2*n)+1/(12*(n**2))+1/12)*b/(a-1)+a/b+0.5*(digamma(a)-log(b))-11/12
  fx2[i]=Hp3
}
plot(x,fx1,type='l')
lines(x,fx2,lty=2)

(a/b)*log(n)-(1/(2*n)+1/(12*(n**2)+1/12)*b/(a-1)-a/b+0.5*(digamma(a)-log(b))-11/12

-1/(6*x)-1/2-log(n*x+1)/n+x-1/n+(1/2)*log(x)
(x*digamma(n*x+1)-lgamma(x))

x*log(n+1/x)-x/(2*(n*x+1))-x/(12*((n*x+1)**2))
+x+log(x)/2-1/(12*x)-11/12
