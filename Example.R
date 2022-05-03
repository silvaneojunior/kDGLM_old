library(Matrix)   # Pode ser substituido pela bilbioteca Matrix.
library(tidyr)
library(dplyr)

library(ggplot2)
library(plotly)
library(latex2exp)

source('R/main.R')

### Exemplos ####

### Importando dados ####
{dados=read.csv('data\\varicela internacoes.csv')[,c(1,7:162)]
dados[1:2,1]='00 a 04 anos'
dados[5:8,1]='15 a 49 anos'
dados[9:12,1]='50 anos e mais'
dados=aggregate(.~FaixaEtaria,dados,sum)[,-1]

pre_exp=read.csv2('data\\populacao 2000-2020.csv')[-12,c(1,10:22)]
pre_exp[4:7,1]='15 a 49 anos'
pre_exp[8:11,1]='50 anos e mais'
pre_exp=aggregate(.~FaixaEtaria,pre_exp,sum)[,-1]

dummy=matrix(0,dim(pre_exp)[1],0)
nomes=c()
for(ano in c(2008:2020)){
  for(mes in c(1:12)){
    nomes=c(nomes,paste0('X',ano,'.',mes))
    dummy=cbind(dummy,pre_exp[,ano-2007])
  }
}
pre_exp=dummy

idade_indice=1
inter=as.numeric(dados[idade_indice,])
exp=as.data.frame(pre_exp)[idade_indice,]
names(exp)=nomes
exp=as.numeric(exp)

N <- dim(dados)[2]
indice_inter=69
data_lab=names(dados)

indic_inter=c(rep(0,indice_inter-1),rep(1,N-indice_inter+1))
covid=c(rep(0,146),rep(1,N-146))
var_covid=ifelse(is.na(covid),0,covid)
}
### Criando modelo ####

nivel_bloc=gera_bloco_poly(2,D=1/0.9,C0=10**2,name='Nível')
inter_bloc=gera_bloco_poly(1,value=indic_inter,D=1/0.9,C0=0.2**2,name='Vacina')
covid_bloc=gera_bloco_poly(1,value=var_covid,D=1/1,C0=0.2**2,name='Covid')
sazo_bloc=gera_bloco_sazo(12,D=1/0.98,name='Sazonalidade')

estrutura=concat_bloco(nivel_bloc,inter_bloc,covid_bloc,sazo_bloc)

teste <- ajusta_modelo(estrutura,data_out=inter,offset=exp)

show_fit(teste,smooth = T)$plot
plot_lat_var(teste,'Nível',smooth = T)

predicao=predict(teste,t=10,offset=exp[N],log_offset=NULL,FF=NULL,D=NULL,W=NULL,plot=T)

# # exp=rpois(100,10)+2
# # y=rpois(100,10*exp*c(rep(1,50),rep(2,50)))
# # D_value=log(y)/rnorm(100,2)
# # indic=c(rep(0,50),rep(1,50))
# #
# #
# # A=gera_bloco_poly(1,value=log(exp),m0=1,C0=0)
# # B=gera_bloco_sazo(12,D=1/0.95)
# # C=gera_bloco_poly(3,D=1/0.9)
# # D=gera_bloco_poly(1,value=D_value,D=1/0.9)
# # E=gera_bloco_poly(1,value=indic,m0=0,C0=0,D=1/0.9)
# # E$W[,,51]=1
# #
# # estrutura=concat_bloco(A,E,D,B,C)
# #
# # teste1=ajusta_modelo(y,structure=estrutura)
# #
# #
# # B=gera_bloco_sazo(12,D=1/0.95)
# # C=gera_bloco_poly(3,D=1/0.9)
# # D=gera_bloco_poly(1,value=D_value,D=1/0.9)
# # E=gera_bloco_poly(1,value=indic,m0=0,C0=0,D=1/0.9)
# # E$W[,,51]=1
# #
# # estrutura=concat_bloco(E,D,B,C)
# #
# # teste2=ajusta_modelo(y,structure=estrutura,offset=exp)
# #
# # plot(teste1$pred)
# # points(teste2$pred)
# #
# # print(mean(abs(teste1$pred-teste2$pred)/teste2$pred))
# #
# # plot(teste1$mts[2,])
#
#

# # Exemplo de uso --------
# y <- c(131.7, 322.6, 285.6, 105.7,
#        80.4, 285.1, 347.8, 68.9,
#        203.3, 375.9, 415.9, 65.8,
#        177.0, 438.3, 463.2, 136.0,
#        192.2, 442.8, 509.6, 201.2,
#        196.0, 478.6, 688.6, 259.8,
#        352.5, 508.1, 701.5, 325.6,
#        305.9, 422.2, 771.0, 329.3,
#        384.0, 472.0, 852.0)
# y <- trunc(y/10)
#
# plot(y, type="l", xlab = "")
# n <- 6
# N <- length(y)
# w  <- 2*pi/4
#
# FF <- matrix(c(rep(1,N),rep(0,N),rep(1,N),rep(0,N),rep(1,N),rep(0,N)), ncol=N,nrow=n, byrow=T)
#
#
# m0 <- matrix(0,ncol=1, nrow=n)
# C0 <- diag(n)
#
# G1 <- matrix(c(1,1,0,1), nrow=2, byrow =T)
# G2 <- matrix(c(cos(w),-sin(w),sin(w),cos(w)),2,2)
# G3 <- matrix(c(cos(2*w),-sin(2*w),sin(2*w),cos(2*w)),2,2)
# G  <- as.matrix(bdiag(G1, G2, G3))
#
# # Matriz de desconto
#
# D=array(0,c(n,n,N))
#
# # Tendência
# delta.1 <- 0.9
# D[1:2,1:2,] <- 1/delta.1
#
# # Sazonalidade
# delta.2 <- 0.9
# D[3:4,3:4,] <- 1/delta.2
#
# # Sazonalidade
# delta.3 <- 0.9
# D[5:6,5:6,] <- 1/delta.3
#
#
# teste <- poisson_gi_exp(y,m0, C0, F1 = FF ,G1 = G, D1 = D, pop=log(1))
# aux <- rep(NA, 9)
# plot(c(aux,y[10:N]), pch=20, ylim=c(0,140), xlab = "", ylab = "")
# lines(c(aux,teste$pred[10:N]), col="#d8790d", lwd = 2)
