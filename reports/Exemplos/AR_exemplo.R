# Gerando o dado
set.seed(1331)
T=200
x=rep(NA,T)
x[1]=rnorm(1)
x[2]=rnorm(1)
phi1=0.99
phi2=-0.1
phi2_list=rep(-0.1,T)
for(t in 3:T){
  phi1=phi1
  phi2=phi2+rnorm(1,0,sqrt(0.0001))
  phi2_list[t]=phi2
  x[t]=phi1*x[t-1]+phi2*x[t-2]+rnorm(1,0,sqrt(0.1))
}
S=0.001
x=x+rnorm(T,0,sqrt(S))
plot(x)
plot(phi2_list,main='phi2')

# Lendo o dado
# dados=read.csv('data/dados_wania_completo.csv')
# x=as.numeric(dados[1,-1])

# Instalando o pacote
# devtools::install_github('silvaneojunior/GDLM')
# Importando o pacote
library(GDLM)

### Modelo AR ###
# theta_i=phi1*theta_{i-1}+phi2*theta_{i-2}

# Fator de desconto para phi1
d_phi1=1
# Fator de desconto para phi2
d_phi2=1
# Fator de desconto para theta_{i-1}
d_theta1=1
# Fator de desconto para theta_{i-2}
d_theta2=1

# Variância do ruído para phi1
w_phi1=0
# Variância do ruído para phi2
w_phi2=1e-4
# FVariância do ruído para theta_{i-1}
w_theta1=0.1
# Variância do ruído para theta_{i-2}
w_theta2=0

# Especificação do choque aleatório:
#  - D especifica o fator de desconto.
#  - W especifica a matriz de convarância do ruído.
# Ordem das variáveis:
#  Para um modelo AR de ordem k, temos 2*k parametros.
#  A ordem dos parâmetros é  theta_{i-1}, phi_1, theta_{i-2}, phi_2, ... , theta_{i-k}, phi_k.
#  Ao especificar o fator de desconto ou a variância do ruído, essa ordem das variáveis deve ser levada em consideração.
block1=AR_block(order=2,
                values=c(1,0),
                by_time=FALSE,
                D=diag(c(1/d_theta1,1/d_phi1,1/d_theta2,1/d_phi2)),
                W=diag(c(w_theta1,w_phi1,w_theta2,w_phi2)),
                name='AR')
               #D=diag(c(1/0.99,1,1,1)),m0=c(0,0.1,0,0.1))

block2=polynomial_block(order=1,
                        values=c(0,1),
                        by_time=FALSE,
                        name='var')

model=fit_model(block1,block2,
                outcome = x,
                family='normal_gamma')

# Visualizando o ajuste
show_fit(model,t_offset=1,smooth=FALSE)$plot

plot(phi2_list,main='phi2',ylim=c(-0.2,0.2))
lines(model$mts[c(4),])
lines(model$mts[c(4),]+1.96*sqrt(model$Cts[4,4,]),lty=2)
lines(model$mts[c(4),]-1.96*sqrt(model$Cts[4,4,]),lty=2)
lines(1:T,rep(0,T),lty=2)
legend('topleft',legend=c('verdadeiro','estimado','i.c. 95%'),lty=c(NA,1,2),pch=c(1,NA,NA))

# Verificando as variáveis disponíveis em model
names(model)
# ver help(analytic_filter)
# verificando phi1 e phi2 no último tempo
model$mts[c(2,4),length(x)]
# verificando variancia no último tempo
exp(-model$mts[5,length(x)])

# Graficos para as variáveis latentes
# Bloco da variância
plot_lat_var(model,'var')$plot
# Bloco da AR
# Suavizado
plot_lat_var(model,'AR')$plot
# Filtrado
plot_lat_var(model,'AR',smooth=FALSE)$plot


# Apenas testes de consistência, não precisa se preocupar com isso.

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

