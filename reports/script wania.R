library(dlm)   # Pode ser substituido pela bilbioteca Matrix.
library(tidyr) # Opcional: ajuda na clareza do código da barra de loading.
library(hms)   # Opcional: serve para formatar o ETA na barra de loading.
library(beepr) # Opcional: serve para chamar a função beep.
library(mgcv)
library(plsmselect)
source('C:\\Jupyter\\TCC\\GDLM\\R\\main.R')

# Ajustando os dados
raw_dados=read.csv('data/dados_wania_completo.CSV',row.names=1)
raw_dados=t(raw_dados) %>% as.data.frame
dados=raw_dados
dados$Tempo=c(1:dim(dados)[1])
dados$N_OBITOS=dados$N_OBITOS %>% as.integer
dados$log_pop=dados$Pac_dia %>% log
dados$intervencao=c(rep(0,54-1),rep(1,101-54+1))
dados$prorp_covid=dados$Entradas_COVID/dados$Entradas

dados$outcome=dados$Soma_Antimicro

suv_t=gam(outcome~intervencao+prorp_covid+s(Tempo),
          family=Gamma ,data=dados)
plot(dados$Tempo,dados$outcome)
lines(dados$Tempo,suv_t$fitted.values)

summary(suv_t)


dados$outcome=dados$Pipetazo

suv_t=gam(outcome~intervencao+prorp_covid+s(Tempo),
          family=Gamma(link='log') ,data=dados)
plot(dados$Tempo,dados$outcome)
lines(dados$Tempo,suv_t$fitted.values)

summary(suv_t)

vars=c("Mero","PoliB","Zerbaxa","Torgena",
       "Amica","Aztreo","Tige","Cefepime",
       "Ceftaroline","Pipetazo","Erta","Cipro",
       "Soma_Antimicro")

vars_index=c()
for(name_index in 1:length(names(dados))){
  if(names(dados)[name_index] %in% vars){
    vars_index=c(vars_index,name_index)
  }
}

plot_data=dados[,c(22,vars_index)] %>% pivot_longer(2:14)

ggplot(plot_data)+
  geom_point(aes(x=Tempo,y=value))+
  geom_vline(xintercept=54,linetype='dashed')+
  geom_hline(yintercept=0,linetype='dashed')+
  facet_wrap(~name,scale='free_y')+
  theme_bw()


tabela=matrix(0,length(vars),9) %>% as.data.frame
count=0

for(var in vars){
  dados$outcome=dados[[var]]
  count=count+1

  suv_t=gam(outcome~intervencao+prorp_covid+s(Tempo),
            family=gaussian ,data=dados)
  plot(dados$Tempo,dados$outcome,main=var)
  lines(dados$Tempo,suv_t$fitted.values)

  mu=suv_t$coefficients[2:3]
  std=sqrt(diag(suv_t$Vp)[2:3])
  cut_off=qt(0.975,df=suv_t$df.residual)

  print(var)
  print(suv_t$df.residual)
  print(abs(std/mu))
  tabela[count,1]=var
  tabela[count,c(2,6)]=mu
  tabela[count,c(3,7)]=mu-cut_off*std
  tabela[count,c(4,8)]=mu+cut_off*std
  tabela[count,c(5,9)]=2*pt(-abs(mu/std),df=suv_t$df.residual)
}


dados$outcome=dados$Mero
suv_t_norm=gam(outcome~intervencao+prorp_covid+s(Tempo),
          family=gaussian ,data=dados)
suv_t_gamm=gam(outcome~intervencao+prorp_covid+s(Tempo),
          family=Gamma, link='log' ,data=dados)
suv_t_lnor=gam(log(outcome)~intervencao+prorp_covid+s(Tempo),
          family=gaussian ,data=dados)
plot(dados$Tempo,dados$outcome)
lines(dados$Tempo,suv_t_norm$fitted.values,col='red')
lines(dados$Tempo,suv_t_gamm$fitted.values,col='green')
lines(dados$Tempo,exp(suv_t_lnor$fitted.values),col='blue')

metric=function(x,y){
  return(abs(x-y)/x)
}

mean(metric(dados$outcome,suv_t_norm$fitted.values))
mean(metric(dados$outcome,suv_t_gamm$fitted.values))
mean(metric(dados$outcome,exp(suv_t_lnor$fitted.values)))

suv_t_norm$aic+672.8049
suv_t_gamm$aic+686.4652
suv_t_lnor$aic+38.62672

# GDLM

T=dim(dados)[1]

dados$outcome=(dados$Mero-mean(dados$Mero))/sd(dados$Mero)

nivel_mu_bloc=gera_bloco_poly(order=1,value=c(rep(1,T),rep(0,T)) %>% matrix(2,T,byrow=TRUE),
                          D=1/0.95,m0=0,C0=1e3,
                          name='Nivel media')

nivel_sig_bloc=gera_bloco_poly(order=1,value=c(rep(0,T),rep(1,T)) %>% matrix(2,T,byrow=TRUE),
                           D=1/1,m0=0,C0=1e3,
                           name='Nivel var')

inter_bloc=gera_bloco_poly(order=1,value=c(dados$intervencao,rep(0,T)) %>% matrix(2,T,byrow=TRUE),
                           D=1/1,m0=0,C0=1,
                           name='Intervenção')

covid_bloc=gera_bloco_poly(order=1,value=c(dados$prorp_covid,rep(0,T)) %>% matrix(2,T,byrow=TRUE),
                           D=1/1,m0=0,C0=1,
                           name='COVID')

resultado=ajusta_modelo(nivel_mu_bloc,nivel_sig_bloc,inter_bloc,covid_bloc,
                        data_out = dados$outcome,
                        kernel='Normal')

show_fit(resultado)$plot
plot_lat_var(resultado,'Intervenção')
plot_lat_var(resultado,'Nivel media')
plot_lat_var(resultado,'Nivel var')

# Dados simulados

T=200

mu=0
s=1

mu0=0
C0=1

y=rnorm(T,mu,s)

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

post.mean=(mu0/C0+cumsum(y)/((s**2)*1:T))/c(1/((s**2)*1:T)+1/C0)

t=3
resultado$mt[1,t]
mean(y[1:t])

plot(resultado$mt[1,],type='l',ylim=c(0.8*mu,1.2*mu))
points(post.mean,col='red')

outcome=resultado$Ct[1,1,]
plot(outcome)
