library(tidyverse)
library(plotly)
library(hms)

library(GDLM)
dados=read.csv('data/dados_wania_gravidade.csv',row.names = 1) %>% t %>% as.data.frame
dados=dados[1:(101-3),]
T=dim(dados)[1]
indice_inter=54

calcula_erro=function(outcome,multi_flag,d_nivel,d_inter,d_covid,d_gravi,ord_nivel){
  nivel=polynomial_block(ord_nivel,values=1,D=1/d_nivel,name='nivel',k=1)
  inter=polynomial_block(1,values=1:T>=indice_inter,D=1/d_inter,C0=0,k=1,name='intervenção')
  inter$W[,,indice_inter]=1
  covid=polynomial_block(1,values=dados$Entradas_COVID/dados$Entradas,C0=0,D=1/d_covid,k=1,name='COVID')
  covid$W[,,87]=1
  gravi=polynomial_block(1,values=dados$N_VM/dados$Pac_dia,D=1/d_gravi,k=1,name='gravidade')

  blocos=block_merge(nivel,inter,covid,gravi)

  if(multi_flag=='TRUE'){
    NERC_bloc=polynomial_block(1,values=dados$N_ERC ,name='NERC',D=1/1,k=1)
    PaRCarba_bloc=polynomial_block(1,values=dados$N_PaR_carb,name='PaRCarba',D=1/1,k=1)
    ESBL_bloc=polynomial_block(1,values=dados$N_ESBL,name='ESBL',D=1/1,k=1)

    blocos=block_merge(blocos,NERC_bloc,PaRCarba_bloc,ESBL_bloc)
  }else{
    if(multi_flag=='relativo'){
      NERC_bloc=polynomial_block(1,values=dados$N_ERC/dados$Pac_dia,name='NERC',D=1/1,k=1)
      PaRCarba_bloc=polynomial_block(1,values=dados$N_PaR_carb/dados$Pac_dia,name='PaRCarba',D=1/1,k=1)
      ESBL_bloc=polynomial_block(1,values=dados$N_ESBL/dados$Pac_dia,name='ESBL',D=1/1,k=1)

      blocos=block_merge(blocos,NERC_bloc,PaRCarba_bloc,ESBL_bloc)
    }
  }

  ajuste=fit_model(blocos,outcome=dados[[outcome]] %>% as.matrix,offset=dados$Pac_dia %>% as.matrix,kernel='Poisson',smooth_flag=FALSE)
  erro=abs(t(ajuste$pred)-ajuste$outcome)/t(ajuste$pred)
  return(mean(erro[erro<Inf]))
}

fileConn1<-file(paste0('reports/Wania/grid_data.csv'))
open(fileConn1,'w')
writeLines(paste('Outcome','multi_flag','Delta.M1','Delta.Inter','Delta.COV','Delta.gravi','Ord.M1','Error',sep=','),fileConn1)

count=0 # Contagem do progresso
init=Sys.time() # Armazenando o valor do tempo de início.

total=(3+3)*21*(11**3)*1

for(outcome in c("N_OBT","N_ERC","N_PaR_carb","N_ESBL")){
  multi_flag_list=if(outcome=="N_OBT"){c('TRUE','FALSE','relativo')}else{c('FALSE')}
  for(multi_flag in multi_flag_list){
    for(d_nivel in 0.8+(0:20/100)){
      for(d_inter in 0.9+(0:10/100)){
        for(d_covid in 0.9+(0:10/100)){
          for(d_gravi in 0.9+(0:10/100)){
            for(ord_nivel in c(1)){
              erro=calcula_erro(outcome,multi_flag,d_nivel,d_inter,d_covid,d_gravi,ord_nivel)
              writeLines(paste(outcome,multi_flag,d_nivel,d_inter,d_covid,d_gravi,ord_nivel,erro,sep=','),fileConn1)
              count=count+1

              perc=count/total
              cur_time=Sys.time()

              qtd1=min(49,floor(49*perc))
              qtd2=49-qtd1

              cat(paste0('[',paste0(rep('=',qtd1),collapse = ''),'>',paste0(rep(' ',qtd2),collapse = ''),']  ',
                         paste(count,'/',total) %>% paste0(' ',(100*perc) %>% round(2) %>% format(nsmall=2) %>% paste0('%')),
                         ' - ETA: ',((1-perc)*difftime(cur_time, init, unit="secs")/perc) %>% as.numeric %>% round %>% hms,
                         '\r'))
            }
          }
        }
      }
    }
  }
}

close(fileConn1)


data=read.csv(paste0('reports/Wania/grid_data.csv'))
ord=order(data$Error)
data=data[ord,]

data_sub=data[data$Delta.Inter==1 & data$Delta.COV==1 & data$Ord.M1==1,]

(ggplot(data_sub)+
  geom_tile(aes(x=Delta.gravi,y=Delta.Inter, fill=Error))+
  theme_bw()) %>% ggplotly

show_fit(ajuste,smooth=FALSE,t_offset = 1)$plot
show_fit(ajuste,smooth=TRUE)$plot

plot_lat_var(ajuste,'nivel')$plot
plot_lat_var(ajuste,'intervenção')$plot
plot_lat_var(ajuste,'COVID')$plot
plot_lat_var(ajuste,'gravidade')$plot
