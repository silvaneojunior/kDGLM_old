library(dlm)   # Pode ser substituido pela bilbioteca Matrix.
library(tidyr) # Opcional: ajuda na clareza do código da barra de loading.
library(hms)   # Opcional: serve para formatar o ETA na barra de loading.
library(beepr) # Opcional: serve para chamar a função beep.
library(feather)
library(plotly)
library(kableExtra)

library(readxl)
library(janitor)
library(ggplot2)
library(dplyr)
library(Matrix)
library(MASS)
library(rootSolve)

source(getwd() %>% paste0('/R/main.R'))

# Ajustando os dados

dados=read.csv('data/varicela internacoes.csv')[,c(1,7:162)]
dados[1:2,1]='00 a 04 anos'
dados[5:8,1]='15 a 49 anos'
dados[9:12,1]='50 anos e mais'
dados=aggregate(.~FaixaEtaria,dados,sum)[,-1]

pre_exp=read.csv2('data/populacao 2000-2020.csv')[-12,c(1,10:22)]
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

t_offset=1
indice_inter=69
true_indice_inter=indice_inter+9

T_final=dim(dados)[2]

out_var=4

#### nível ####

FF_nivel=matrix(0,out_var,T_final)
FF_nivel[1,]=1
bloc_final=gera_bloco_poly(order=2,
                           value=FF_nivel,
                           name='nivel_serie_1',
                           D=1/0.95,
                           m0=0,
                           C0=500,
                           W=0)


for(i in c(2:out_var)){
  FF_nivel=matrix(0,out_var,T_final)
  FF_nivel[i,]=1
  bloc_final=concat_bloco(bloc_final,
                          gera_bloco_poly(order=2,
                                          value=FF_nivel,
                                          name='nivel_serie_' %>% paste0(i),
                                          D=1/0.95,
                                          m0=0,
                                          C0=500,
                                          W=0)
  )
}
for(i in c(1:out_var)){
  FF_sazo=matrix(0,out_var,T_final)
  FF_sazo[i,]=1
  bloc_final=concat_bloco(bloc_final,
                          gera_bloco_sazo(period=2*pi/12,
                                          value=FF_sazo,
                                          name='sazo_serie_' %>% paste0(i),
                                          D=1/0.98,
                                          m0=0,
                                          C0=500,
                                          W=0)
  )
}
for(i in c(1:out_var)){
  FF_vac=matrix(0,out_var,T_final)
  FF_vac[i,true_indice_inter:T_final]=1
  W=array(0,c(1,1,T_final))
  W[,,true_indice_inter]=1
  bloc_final=concat_bloco(bloc_final,
                          gera_bloco_poly(order=1,
                                          value=FF_vac,
                                          name='vac_serie_' %>% paste0(i),
                                          D=1/1,
                                          m0=0,
                                          C0=0,
                                          W=W)
  )
}
for(i in c(1:out_var)){
  FF_cov=matrix(0,out_var,T_final)
  FF_cov[i,146:T_final]=1
  W=array(0,c(1,1,T_final))
  W[,,146]=1
  bloc_final=concat_bloco(bloc_final,
                          gera_bloco_poly(order=1,
                                          value=FF_cov,
                                          name='cov_serie_' %>% paste0(i),
                                          D=1/1,
                                          m0=0,
                                          C0=0,
                                          W=W)
  )
}


y=t(dados)
y=y[,c(1,2,3,5,4)]

resultado=ajusta_modelo(bloc_final,
                        data_out=y,
                        kernel=multinom_gi)

pre_ps=exp(resultado$ft)/(rowSums(exp(resultado$ft[,1:out_var]))+1)

ps=matrix(NA,T_final,5)
ps_i=matrix(NA,T_final,5)
ps_s=matrix(NA,T_final,5)

for(i in c(1:T_final)){
  p=pre_ps[i,]
  var=resultado$Qt[,,i]
  
  diag_mult=diag(p*(1-p))
  cov_mult=diag_mult%*%var%*%diag_mult
  
  p_i=p-2*sqrt(diag(cov_mult))
  p_s=p+2*sqrt(diag(cov_mult))
  
  vec=matrix(1,4,1)
  vec_rest=1-sum(p)
  var_rest=t(vec)%*%cov_mult%*%vec
  
  p=c(p,vec_rest)
  p_i=c(p_i,vec_rest-2*sqrt(var_rest))
  p_s=c(p_s,vec_rest+2*sqrt(var_rest))
  
  ps[i,]=p
  ps_i[i,]=p_i
  ps_s[i,]=p_s
}

plot_data=dados[c(1,2,3,5,4),]
for(i in c(1:T_final)){
  plot_data[,i]=plot_data[,i]/sum(plot_data[,i])
}
row.names(plot_data)=c(1:5)

plot_data=plot_data %>%
  t %>%
  as.data.frame %>%
  mutate(t=1:T_final,) %>%
  pivot_longer(1:5) %>%
  mutate(name='V'%>%paste0(name),
         value=value)

plot_ps=ps %>%
  as.data.frame %>%
  mutate(t=1:T_final) %>%
  pivot_longer(1:5)

plot_ps_i=ps_i %>%
  as.data.frame %>%
  mutate(t=1:T_final) %>%
  pivot_longer(1:5)

plot_ps_s=ps_s %>%
  as.data.frame %>%
  mutate(t=1:T_final) %>%
  pivot_longer(1:5)

plot_data=plot_data %>%
  inner_join(plot_ps,c('t','name')) %>%
  inner_join(plot_ps_i,c('t','name')) %>%
  inner_join(plot_ps_s,c('t','name'))

plot_data$name=as.factor(plot_data$name)
levels(plot_data$name)=c('00 a 04 anos',
                         '05 a 09 anos',
                         '10 a 14 anos',
                         '50 anos e mais',
                         '15 a 49 anos')

names(plot_data)=c('Data','Faixa.Etaria','Valor.Observado','Valor.esperado','I.C.inf','I.C.sup')

ggplotly(
  ggplot(plot_data)+
    geom_point(aes(x=Data,y=Valor.Observado,color=Faixa.Etaria,linetype='Observado',shape='Observado'))+
    geom_line(aes(x=Data,y=Valor.esperado,color=Faixa.Etaria,linetype='Estimado',shape='Estimado'))+
    scale_x_continuous('Data',breaks=c(0:13)*12,labels=c(0:13)+2008)+
    scale_y_continuous('Proporção de internações', breaks=,labels=function(x){round(100*x) %>% paste('%')})+
    scale_linetype('')+
    scale_shape('')+
    geom_ribbon(aes(x=Data,ymin=I.C.inf,ymax=I.C.sup,color=Faixa.Etaria,fill=Faixa.Etaria,shape='I.C.'),alpha=0.5)+
    geom_vline(xintercept=indice_inter+9,linetype='dashed')+
    theme_bw()+
    coord_cartesian(xlim=c(12,156),ylim=c(0,1))
)