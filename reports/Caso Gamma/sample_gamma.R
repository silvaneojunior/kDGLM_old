devtools::load_all()
library(tidyverse)
library(extraDistr)
library(latex2exp)

densi_conj=function(x,tau,k){
  p1=k*(x*log(x)-lgamma(x))
  return(-tau*x+p1)
}

# densi_conj=function(x,alpha,beta){
#   return(dgamma(x,alpha,beta,log=TRUE))
# }

# phi=2
# mu=exp(1+1:T/T)
#
# T=200
#
# outcome=rgamma(T,phi,phi/mu)

# Dados IBM
outcome=read.csv('data/m-ibmln.txt',header=FALSE) %>% as.matrix %>% as.numeric
outcome=(outcome-mean(outcome))**2
# outcome=outcome[1:10]
# outcome=rgamma(100,2,2/1)
T=length(outcome)

phi=1/2
mu=mean(outcome)

sample_size=2000
tau_0=1
k_0=1

phi_sample=rep(NA,sample_size)
mu_sample=matrix(NA,T,sample_size)
mt_sample=matrix(NA,T,sample_size)
pred_sample=matrix(NA,T,sample_size)

mu_i=mean(outcome)
phi_i=1/2
gamma_size=2000
gamma_sample=rgamma(gamma_size,(T/2),(T/2)/phi_i)
gamma_densi=dgamma(gamma_sample,(T/2),(T/2)/phi_i,log=TRUE)

alpha_prop=(T/2)
beta_prop=(T/2)/phi_i

level=polynomial_block(1,D=1/0.95)

for(i in 1:sample_size){
  cat(paste0(i,'                          \r'))
  fitted_model=fit_model(level,outcome=outcome,parms=list('phi'=phi_i),family='gamma',pred_cred=-1,smooth_flag = FALSE)
  sample=FFBS_sampling(fitted_model,1)
  mu_i=sample$param

  tau_i=tau_0+sum(outcome/mu_i-log(outcome/mu_i))
  k_i=k_0+T
  # tau_i=tau_0
  # k_i=k_0

  weights_unscaled=densi_conj(gamma_sample,tau_i,k_i)-gamma_densi
  # outcome_mat=matrix(outcome,T,gamma_size)
  # gamma_mat=matrix(gamma_sample,T,gamma_size,byrow=TRUE)
  # mu_i_mat=matrix(mu_i,T,gamma_size)
  #
  # weights_mat=colSums(dgamma(outcome_mat,gamma_mat,gamma_mat/mu_i_mat,log=TRUE))

  # weights_unscaled=densi_conj(gamma_sample,tau_i,k_i)+weights_mat-gamma_densi
  weights_raw=exp(weights_unscaled-max(weights_unscaled))
  weights=weights_raw/sum(weights_raw)
  phi_i=sample(gamma_sample,1,prob=weights)

  alpha_prop=(T/2)
  beta_prop=(T/2)/sum(weights*gamma_sample)

  gamma_sample=rgamma(gamma_size,alpha_prop,beta_prop)
  gamma_densi=dgamma(gamma_sample,alpha_prop,beta_prop,log=TRUE)
  pred_i=rgamma(T, phi_i, phi_i/mu_i[1,,1])

  mu_sample[,i]=mu_i[1,,]
  mt_sample[,i]=sample$mt[1,,]
  phi_sample[i]=phi_i
  pred_sample[,i]=pred_i
}
sample_index=(1:1000)*2
ts.plot(phi_sample)
acf(phi_sample[sample_index])
print(mean(phi_sample[sample_index]))
print(quantile(phi_sample[sample_index],0.975))
print(quantile(phi_sample[sample_index],0.025))

acf(mu_sample[1,sample_index])
print(mean(mu_sample[1,sample_index]))
print(quantile(mu_sample[1,sample_index],0.975))
print(quantile(mu_sample[1,sample_index],0.025))

acf(pred_sample[1,sample_index])
print(mean(pred_sample[1,sample_index]))
print(quantile(pred_sample[1,sample_index],0.975))
print(quantile(pred_sample[1,sample_index],0.025))

year <- seq(as.Date("1926/10/1"), by = "month", length.out = length(outcome)+12)
year_label <- seq(as.Date("1930/1/1"), by = "10 years", length.out = 8)

font_size=16
(ggplot()+
    geom_line(aes(x=1:T,y=mt_sample[,sample_index] %>% rowMeans,linetype='Point estimation'))+
    geom_ribbon(aes(x=1:T,
                    ymax=mt_sample[,sample_index] %>% apply(1,function(x){quantile(x,0.975)}),
                    ymin=mt_sample[,sample_index] %>% apply(1,function(x){quantile(x,0.025)})),
                alpha=0.25,fill='black')+
    scale_x_continuous('Date',breaks=c(0:14)*12*10+1+3+36,labels=function(x){substr(year[x],1,4)},expand=c(0,0),limits=c(39,NA))+
    scale_y_continuous('Estimated value',expand=c(0,0,0.01,0),limits=c(NA,0.03))+
    scale_linetype_manual('',values=c('solid'))+
    # coord_cartesian(ylim=c(0,0.1),xlim=c(700,800))+
    theme_bw()+
    theme(text=element_text(size=font_size)))+
  labs(title='Level')

ggsave(
  'C:\\Jupyter\\Mestrado\\Pacote\\GDLM\\reports\\Graficos artigo raira\\gamma_mt.png',
  scale = 1,
  units='px',
  width=800*4,
  height=600*4,
  dpi = 300
)


(ggplot()+
    geom_line(aes(x=1:T,y=mu_sample[,sample_index] %>% rowMeans,linetype='Point estimation'))+
    geom_ribbon(aes(x=1:T,
                    ymax=mu_sample[,sample_index] %>% apply(1,function(x){quantile(x,0.975)}),
                    ymin=mu_sample[,sample_index] %>% apply(1,function(x){quantile(x,0.025)})),
                alpha=0.25,fill='black')+
    geom_line(aes(x=2:T,y=exp(-fitted_data$mts[2,]+fitted_data$Cts[2,2,]/2),linetype='Normal-gamma'))+
    geom_ribbon(aes(x=2:T,
                    ymin=qlnorm(0.025,-fitted_data$mts[2,],sqrt(fitted_data$Cts[2,2,])),
                    ymax=qlnorm(0.975,-fitted_data$mts[2,],sqrt(fitted_data$Cts[2,2,])),
                    linetype='Normal-gamma'),alpha=0.25,fill='blue')+
    scale_x_continuous('Date',breaks=c(0:14)*12*10+1+3+36,labels=function(x){substr(year[x],1,4)},expand=c(0,0),limits=c(39,NA))+
    scale_y_continuous('Estimated value',expand=c(0,0,0.01,0),limits=c(NA,0.03))+
    scale_linetype_manual('',values=c('solid','dashed'))+
    # coord_cartesian(ylim=c(0,0.1),xlim=c(700,800))+
    theme_bw()+
    theme(text=element_text(size=font_size)))+
  labs(title='Estimation of \\mu' %>% TeX)

ggsave(
  'C:\\Jupyter\\Mestrado\\Pacote\\GDLM\\reports\\Graficos artigo raira\\gamma_mu.png',
  scale = 1,
  units='px',
  width=800*4,
  height=600*4,
  dpi = 300
)

zoom_factor=2

# Tamanho da fonte
font_size=16

(ggplot()+
  geom_line(aes(x=1:T,y=pred_sample[,sample_index] %>% rowMeans,linetype='Predictive\nmean'))+
  geom_ribbon(aes(x=1:T,
                  ymax=pred_sample[,sample_index] %>% apply(1,function(x){quantile(x,0.975)}),
                  ymin=pred_sample[,sample_index] %>% apply(1,function(x){quantile(x,0.025)})),
              alpha=0.25,color=NA,fill='black')+
  geom_point(aes(x=1:T,y=outcome,shape='Observation'))+
  scale_x_continuous('Date',breaks=c(0:6)*12*10+1+3+36,labels=function(x){substr(year[x],1,4)},expand=c(0,0),limits=c(39,NA))+
  scale_y_continuous('Squared returns',expand=c(0,0,0.01,0))+
    scale_linetype_manual('',values='dotdash')+
    scale_shape('')+
  # coord_cartesian(ylim=c(0,0.1),xlim=c(700,800))+
  theme_bw()+
  theme(text=element_text(size=font_size))) %>% plot_zoom(mini_range=list('x'=c(700,820),'y'=c(-0.001,0.05)),
                            mini_place=list('x'=100+c(0,zoom_factor*320),'y'=0.1+c(0,zoom_factor*0.05)))+
  labs(title='IBM monthly returns')


ggsave(
  'C:\\Jupyter\\Mestrado\\Pacote\\GDLM\\reports\\Graficos artigo raira\\gamma_pred.png',
  scale = 1,
  units='px',
  width=800*4,
  height=600*4,
  dpi = 300
)

ggplot(data.frame(phi=phi_sample[sample_index]),aes(x=phi))+
  geom_histogram(binwidth = 0.01,fill='white',color='black')+
  # geom_density(aes(),color='black')+
  scale_x_continuous('\\phi' %>% TeX,expand=c(0.02,0))+
  scale_y_continuous('Frequency',expand=c(0,0,0,10))+
  labs(title='Distribuition of the shape parameter (\\phi)' %>% TeX)+
  theme_bw()

ggsave(
  'C:\\Jupyter\\Mestrado\\Pacote\\GDLM\\reports\\Graficos artigo raira\\gamma_phi.png',
  scale = 1,
  units='px',
  width=800*4,
  height=600*4,
  dpi = 300
)


#########################################################

devtools::load_all()

library(tidyverse)
library(plotly)

T <- 200
w <- (200 / 50) * 2 * pi
S=exp(2*(sin(w * 1:T / T)))*exp(-20*(1:T/T)*((1:T/T)-1))/20
#S=exp(-10*(1:T/T)*((1:T/T)-1))
mu <- 20*1:T / T
set.seed(13031998)
outcome <- rnorm(T, mu,sqrt(1/S))

level <- polynomial_block(
  order = 2,
  values=c(1,0),
  D = 1 / 1,
  by_time=FALSE
)
variance1 <- polynomial_block(
  order = 3,
  values=c(0,1),
  D = 1 / 1,
  C0= diag(c(1,0.1,0.01)),
  by_time=FALSE
)
variance2 <- harmonic_block(
  period = 50,
  values = c(0, 1),
  D = 1 / 1,
  by_time=FALSE
)

fitted_data <- fit_model(level, variance1, variance2,
                         outcome = outcome,
                         family = "normal_gamma_cor")

show_fit(fitted_data,smooth = TRUE)$plot

sample=FFBS_sampling(fitted_data,2000)

sample$param[2,,] %>% rowMeans %>% plot(type='l')
points(1/S)
lines(1/S)
points(1/S)

ggplot()+
  geom_line(aes(x=1:T,y=sample$param[2,,] %>% rowMeans),linetype='dashed')+
  geom_ribbon(aes(x=1:T,
                  ymax=sample$param[2,,] %>% apply(1,function(x){quantile(x,0.975)}),
                  ymin=sample$param[2,,] %>% apply(1,function(x){quantile(x,0.025)})),
              alpha=0,linetype='dashed',color='black')+
  geom_line(aes(x=1:T,y=S))+
  theme_bw()

#########################################################

devtools::load_all()

library(tidyverse)
library(plotly)

T <- 200
y1 <- rpois(T, exp(5 + (-T:T / T) * 5))
y2 <- rpois(T, exp(6 + (-T:T / T) * 5 + sin((-T:T) * (2 * pi / 12))))
y3 <- rpois(T, exp(5))

y <- cbind(y1, y2, y3)

level <- polynomial_block(2, k = 2)
season <- harmonic_block(12, values = c(0, 1), by_time = FALSE)

fitted_data <- fit_model(level, season, outcome = y, family = "Multinomial", pred_cred = 0.95)

show_fit(fitted_data, smooth = TRUE)$plot

sample=FFBS_sampling(fitted_data,2000)

sample$param[2,,] %>% rowMeans %>% plot(type='l')
points(1/S)
lines(1/S)
points(1/S)

ggplot()+
  geom_line(aes(x=1:T,y=sample$param[2,,] %>% rowMeans),linetype='dashed')+
  geom_ribbon(aes(x=1:T,
                  ymax=sample$param[2,,] %>% apply(1,function(x){quantile(x,0.975)}),
                  ymin=sample$param[2,,] %>% apply(1,function(x){quantile(x,0.025)})),
              alpha=0,linetype='dashed',color='black')+
  geom_line(aes(x=1:T,y=S))+
  theme_bw()
