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
T=length(data)
data=data.frame(mu=data,phi=rep(0,T)) %>% as.matrix

# GDLM
values=c(700:1000)/1000
erros=rep(0,length(values))

# Ainda n fiz documentação das funções usadas nesse script, por isso estou incluindo comentários explicando o básico.

# gera_bloco_poly é uma função auxiliar q cria um bloco polinomial com a ordem desejada.
# o argumento value deve receber uma matriz com os valores da matriz F.
# os valores da primeira linha da matriz serão associados ao nível mu.
# os valores da segunda linha da matriz serão associados a log precisão.
nivel_mu_bloc=gera_bloco_poly(order=1,value=c(rep(1,T),rep(0,T)) %>% matrix(2,T,byrow=TRUE),
                              D=1/0.95,m0=0,C0=1e0,
                              name='Nivel media')

nivel_sig_bloc=gera_bloco_poly(order=1,value=c(rep(0,T),rep(1,T)) %>% matrix(2,T,byrow=TRUE),
                               D=1/0.98,m0=0,C0=1e0,
                               name='Nivel var')
for(i in 1:length(values)){
  print(i)

  # Escolhendo o gamma.
  nivel_sig_bloc$G=values[i]

  # a função ajusta modelo faz o ajuste do modelo (surprise!), os argumentos sem nome são os blocos de variáveis.
  # Internamente a função fará a concatenação dos blocos.
  resultado=ajusta_modelo(nivel_mu_bloc,nivel_sig_bloc,
                          data_out = data,
                          kernel='Normal')

  resultado$log.ver=ifelse(is.na(resultado$log.vero),Inf,resultado$log.vero)
  if(any(is.na(resultado$log.vero))){

  }
  print(resultado$log.vero[1])
  erros[i]=sum(resultado$log.vero)
}

plot(values,erros,type='l')
points(values,erros,type='p',pch='O',cex=0.5)

# Testa ajuste para um valor de gamma
nivel_sig_bloc$G=0.9

# a função ajusta modelo faz o ajuste do modelo (surprise!), os argumentos sem nome são os blocos de variáveis.
# Internamente a função fará a concatenação dos blocos.
resultado=ajusta_modelo(nivel_mu_bloc,nivel_sig_bloc,
                        data_out = data,
                        kernel='Normal')

show_fit(resultado,smooth=FALSE, t_offset=1)$plot

plot(resultado$pred[1,],type='l')
points(data[,1])


plot_lat_var(resultado,'Nivel var')
