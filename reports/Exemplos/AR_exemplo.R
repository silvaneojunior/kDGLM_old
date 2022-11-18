T=200
x=rep(NA,T)
x[1]=rnorm(1)
x[2]=rnorm(1)
phi1=0.99
phi2=-0.1
for(t in 3:T){
  x[t]=phi1*x[t-1]+phi2*x[t-2]+rnorm(1,0,sqrt(0.1))
}
S=0.001
x=x+rnorm(T,0,sqrt(S))
plot(x)

devtools::load_all()

block1=AR_block(2,values=c(1,0),by_time=FALSE,W=diag(c(0.1,0,0,0)))#,D=diag(c(1/0.99,1,1,1)),m0=c(0,0.1,0,0.1))
block2=polynomial_block(1,values=c(0,1),by_time=FALSE)
model=fit_model(block1,block2,outcome = x,family='normal_gamma')
# show_fit(model,smooth=FALSE,t_offset=10)$plot
show_fit(model)$plot

model$mts[c(2,4),T]
exp(-model$mts[5,T])

arima(x,order=c(2,0,0),include.mean = FALSE, transform.pars = FALSE)


sample_size=2000
phi_sample=matrix(NA,2,sample_size)
mt_sample=array(NA,c(2,T,sample_size))
phi=c(phi1,phi2)
for(i in 1:sample_size){
  cat(i,'                            \r')
  block1=polynomial_block(2,values=1,by_time=TRUE,W=diag(c(0.1,0)),D=1)
  block1$G[1,1]=phi[1]
  block1$G[1,2]=phi[2]
  block1$G[2,1]=1
  block1$G[2,2]=0
  model=fit_model(block1,outcome = x,family='normal',parms=list('Sigma'=S))
  mt_sample_i=FFBS_sampling(model,1)$mt[,,1]

  block1=polynomial_block(1,values=mt_sample_i[1,-T],by_time=TRUE,W=0,D=1)
  block2=polynomial_block(1,values=mt_sample_i[2,-T],by_time=TRUE,W=0,D=1)
  model=fit_model(block1,block2,outcome = mt_sample_i[1,-1],family='normal',parms=list('Sigma'=0.1))
  phi_sample_i=FFBS_sampling(model,1)$mt[,T-1,1]

  phi=phi_sample_i
  mt_sample[,,i]=mt_sample_i
  phi_sample[,i]=phi
}

plot(phi_sample[1,])
acf(phi_sample[1,])
rowMeans(phi_sample[,seq(1,sample_size,5)])
save.image(file = "AR_exemplo.RData")

