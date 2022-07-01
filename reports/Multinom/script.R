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

# library(GDLM)
source('reports/Multinom/R/main.R')
source('reports/Multinom/R/kernels.R')
source('reports/Multinom/R/plot_helper.R')
source('reports/Multinom/R/structure_helper.R')

dados=read.csv('reports/Multinom/data/varicela internacoes.csv')[,c(1,7:162)]
dados[1:2,1]='00 a 04 anos'
dados[5:8,1]='15 a 49 anos'
dados[9:12,1]='50 anos e mais'
dados=aggregate(.~FaixaEtaria,dados,sum)
labels=dados[,1]
dados=dados[,-1]

pre_exp=read.csv2('reports/Multinom/data/populacao 2000-2020.csv')[-12,c(1,10:22)]
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
#idade_indice=3

proxy_vac=read.csv2('reports/Multinom/data/proxy_vac.csv')

vac_corb=dados[1,]*0
for(ano in c(2013:2019)){
  indices=substr(names(vac_corb),2,5)==as.character(ano)
  vac_corb[indices]=proxy_vac[ano-2012,2]
}
indices=substr(names(vac_corb),2,5)==as.character(2020)
vac_corb[indices]=proxy_vac[2019-2012,2]

t_offset=1
indice_inter=69
true_indice_inter=indice_inter+9

T_final=dim(dados)[2]

FaixaEtaria=4

out_var=4
ref_var=FaixaEtaria %>% as.numeric
indices=c(c(1:5)[c(1:5)!=ref_var],ref_var)

y=t(dados)
y=y[,indices]

ord_labels=labels[indices]

offset=pre_exp
offset=offset[indices,]
offset[1,]=offset[1,]/offset[5,]
offset[2,]=offset[2,]/offset[5,]
offset[3,]=offset[3,]/offset[5,]
offset[4,]=offset[4,]/offset[5,]

offset=log(offset[-5,])
#### offset ####

bloc_final=gera_bloco_poly(order=1,
                           value=offset,
                           name='offset',
                           D=1/1,
                           m0=1,
                           C0=0,
                           W=0)


for(i in c(1:out_var)){
  #### nível ####
  FF_static=matrix(0,out_var,T_final)
  FF_static[i,]=1
  bloc_nivel=gera_bloco_poly(order=2,value=FF_static,name='nivel_serie_' %>% paste0(i),
                             D=1/0.95,m0=0,C0=1,W=0)
  #### sazonalidade ####
  bloc_sazo=gera_bloco_sazo(period=12,value=FF_static,name='sazo_serie_' %>% paste0(i),
                            D=1/0.98,m0=0,C0=1,W=0)
  #### vacina ####
  FF_vac=matrix(0,out_var,T_final)
  # Caso seja usado a indicadora
  #FF_vac[i,true_indice_inter:T_final]=1
  # Caso seja usado a cobertura vacinal
  FF_vac[i,]=as.numeric(vac_corb)
  W=array(0,c(1,1,T_final))
  # Caso seja usado a indicadora
  #W[,,true_indice_inter]=1
  # Caso seja usado a cobertura vacinal
  W[,,true_indice_inter]=1
  bloc_vac=gera_bloco_poly(order=1,value=FF_vac,name='vac_serie_' %>% paste0(i),
                           D=1/0.99,m0=0,C0=0,W=W)
  #### covid ####
  FF_cov=matrix(0,out_var,T_final)
  FF_cov[i,146:T_final]=1
  W=array(0,c(1,1,T_final))
  W[,,146]=1
  bloc_cov=gera_bloco_poly(order=1,value=FF_cov,name='cov_serie_' %>% paste0(i),
                           D=1/1,m0=0,C0=0,W=W)

  bloc_final=concat_bloco(bloc_final,bloc_nivel,bloc_sazo,bloc_vac,bloc_cov)
}

resultado=ajusta_modelo(bloc_final,
                        data_out=y,
                        kernel='Multinomial')

indices=c(c(1:5)[c(1:5)!=ref_var],ref_var)

ord_labels=labels[indices]

(show_fit(resultado,labels=ord_labels,dinamic=FALSE,t_offset=1,smooth=FALSE)$plot+
    geom_vline(xintercept=true_indice_inter,linetype='dashed')+
    labs(title='Previsão um passo à frente')+
    scale_y_continuous('Internações')+
    scale_x_continuous('Data', breaks=c(0:12)*12+1,labels=c(2008:2020))+
    theme(axis.text = element_text(angle=90))+
    coord_cartesian(ylim=c(0,1000))) %>% ggplotly

(plot_lat_var(resultado,'vac_serie',dinamic=FALSE)+
  geom_vline(xintercept=true_indice_inter,linetype='dashed')+
  labs(title='Previsão um passo à frente')+
  scale_y_continuous('Internações')+
  scale_x_continuous('Data', breaks=c(0:12)*12+1,labels=c(2008:2020))+
  theme(axis.text = element_text(angle=90))) %>% ggplotly

(plot_lat_var(resultado,'vac_serie',dinamic=FALSE)+
    geom_vline(xintercept=true_indice_inter,linetype='dashed')+
    labs(title='Previsão um passo à frente')+
    scale_y_continuous('Internações')+
    scale_x_continuous('Data', breaks=c(0:12)*12+1,labels=c(2008:2020))+
    theme(axis.text = element_text(angle=90))) %>% ggplotly
