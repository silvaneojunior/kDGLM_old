library(dlm)   # Pode ser substituido pela bilbioteca Matrix.
library(kableExtra)
library(DT)
library(scales)

library(plotly)
library(ggplot2)

library(shiny)
library(shinythemes)

library(MASS)
library(rootSolve)
library(tidyr)
library(dplyr)
source('R\\main.R')
source('R\\kernels.R')
source('R\\plot_helper.R')
source('R\\structure_helper.R')

data=read.csv('data/m-ibmln.txt')[,1]*10

# GDLM

T=length(data)
data=data.frame(mu=data,phi=rep(0,T)) %>% as.matrix

# Ainda n fiz documentação das funções usadas nesse script, por isso estou incluindo comentários explicando o básico.

# gera_bloco_poly é uma função auxiliar q cria um bloco polinomial com a ordem desejada.
# o argumento value deve receber uma matriz com os valores da matriz F.
# os valores da primeira linha da matriz serão associados ao nível mu.
# os valores da segunda linha da matriz serão associados a log precisão.
nivel_mu_bloc=gera_bloco_poly(order=1,value=c(rep(1,T),rep(0,T)) %>% matrix(2,T,byrow=TRUE),
                              D=1/0.95,m0=0,C0=1e0,
                              name='Nivel media')

nivel_sig_bloc=gera_bloco_poly(order=1,value=c(rep(0,T),rep(1,T)) %>% matrix(2,T,byrow=TRUE),
                               D=1/1,m0=0,C0=1e0,
                               name='Nivel var')
# Escolhendo o gamma.
nivel_sig_bloc$G=1

# a função ajusta modelo faz o ajuste do modelo (surprise!), os argumentos sem nome são os blocos de variáveis.
# Internamente a função fará a concatenação dos blocos.
resultado=ajusta_modelo(nivel_mu_bloc,nivel_sig_bloc,
                        data_out = data,
                        kernel='Normal')

# Essas funções auxiliam na vizualização dos dados
# show_fit mostra o ajuste (pode ser suavizado, filtrado ou a previsão k passos a frente)
# para a previsão k passos a frente coloque smooth=FALSE e offset=k
show_fit(resultado,smooth=FALSE)$plot
# plot_lat_var exibe o valor estimado das variáveis latentes.
# o argumento var deve ser igual ao nome do bloco (argumento name).
# smooth indica se o valor exibido deve ser filtrado ou suavizado
plot_lat_var(resultado,var='Nivel media',smooth=TRUE)
plot_lat_var(resultado,var='Nivel var',smooth=TRUE)

# Dados simulados

T=200

mu=0
s=1

mu0=0
C0=1
# offset da média dos dados observado.
offset=0

set.seed(13031998)
y=rnorm(T,mu,s)+offset

nivel_mu_bloc=gera_bloco_poly(order=1,value=c(rep(1,T),rep(0,T)) %>% matrix(2,T,byrow=TRUE),
                              D=1/1,m0=mu0,C0=C0,
                              name='Nivel media')

nivel_sig_bloc=gera_bloco_poly(order=1,value=c(rep(0,T),rep(1,T)) %>% matrix(2,T,byrow=TRUE),
                               D=1/1,m0=0,C0=1,
                               name='Nivel var')

resultado=ajusta_modelo(nivel_mu_bloc,nivel_sig_bloc,
                        data_out = y,
                        kernel='Normal')


show_fit(resultado)$plot

# priori normal gamma (segundo a compatilização feita no código)
mu0=0.0000000
c0=0.6065307
n0=2.2754474
d0=1.3801286

n_t=1:T
mean_t=cumsum(y)/n_t
s_t=rep(0,T)
for(i in 1:T){
  sub_x=y[1:i]
  s_t[i]=(sub_x-mean(sub_x))/i
}

post.c=c0+n_t
post.mu=((c0*mu0)+n_t*mean_t)/post.c
post.n=n0+n_t/2
post.d=d0+0.5*(n_t*s_t+c0*n_t*((mean_t-mu0)**2)/post.c)

post.mean=post.mu

t=3
resultado$mt[1,t]
mean(y[1:t])

# linha = média filtrada de mu no modelo
# pontos = média a posteriori mu no modelo normal-gamma a cada tempo t
plot(resultado$mt[1,],
     type='l',ylim=c(0.8*mu+offset,1.2*mu+offset),
     main='media do parametro mu')
points(post.mean,col='red')
legend(150,1+offset,legend=c('GDLM','sol. analítica'),col=c('black','red'),lty=c(1,1))

outcome=resultado$Ct[1,1,]
plot(outcome)
