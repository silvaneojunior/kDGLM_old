library(GDLM)
library(plotly)
library(tidyverse)

dados=read.csv('data/dados_wania_gravidade.csv',row.names = 1) %>% t %>% as.data.frame
dados=dados[1:(101-3),]
T=dim(dados)[1]
indice_inter=54

nivel=polynomial_block(1,values=1,D=1/1,name='nivel')
inter=polynomial_block(1,values=1:T>=indice_inter,D=1/1,C0=0,W=1:T==indice_inter,k=1,name='intervenção')
covid=polynomial_block(1,values=dados$Entradas_COVID/dados$Entradas,C0=0,W=1:T==sum(dados$Entradas_COVID==0),D=1/1,k=1,name='COVID')
gravi=polynomial_block(1,values=dados$N_VM/dados$Pac_dia,D=1/0.97,k=1,name='gravidade')

N_ERC=polynomial_block(1,values=dados$N_ERC,D=1/1,k=1,name='N_ERC')
N_PRC=polynomial_block(1,values=dados$N_PaR_carb,D=1/1,k=1,name='N_PaRCarb')
N_ESB=polynomial_block(1,values=dados$N_ESBL,D=1/1,k=1,name='N_ESBL')


ajuste=fit_model(nivel,inter,covid,gravi,
                 N_ERC,N_PRC,N_ESB,
                 outcome=dados$N_OBT %>% as.matrix,
                 offset=dados$Pac_dia %>% as.matrix,
                 kernel='Poisson')

plot_data=show_fit(ajuste,smooth=FALSE,t_offset = 1)

fill_list=c('#2596be','black','#2596be')

(ggplot()+
    geom_point(aes(x=as.numeric(1:T),y=plot_data$pred$Prediction,color='Predição',fill='Predição'))+
    geom_ribbon(aes(x=as.numeric(1:T),
                    ymin=plot_data$icl.pred$I.C.lower,
                    ymax=plot_data$icu.pred$I.C.upper,
                    color='I.C.',fill='I.C.'),alpha=0.25)+
    geom_point(aes(x=as.numeric(1:T),y=dados$N_OBT,
                   color='Observações',fill='Observações'))+
    geom_vline(xintercept=indice_inter,linetype='dashed')+
    scale_x_continuous('Ano',breaks=round(c(0:8)*12),lim=c(1,T),labels=c(2013:2021))+
    scale_y_continuous('Óbitos',expand=c(0,0))+
    scale_color_manual('',values=fill_list)+
    scale_fill_manual('',values=fill_list)+
    coord_cartesian(ylim=c(0,100))+
    theme_bw()+
    theme(axis.text.x = element_text(angle=90))) %>% ggplotly

show_fit(ajuste,smooth=TRUE)$plot

plot_lat_var(ajuste,'nivel')$plot
plot_lat_var(ajuste,'intervenção')$plot
plot_lat_var(ajuste,'COVID')$plot
plot_lat_var(ajuste,'gravidade')$plot

plot_lat_var(ajuste,'N_ERC')$plot
plot_lat_var(ajuste,'N_PaRCarb')$plot
plot_lat_var(ajuste,'N_ESBL')$plot



# anti_micro=c("Mero","PoliB","Zerbaxa",
# "Torgena","Amica","Aztreo","Tige",
# "Cefepime","Ceftaroline","Pipetazo","Erta","Cipro")

anti_micro_dados=dados[,10:21]
anti_micro_dados=anti_micro_dados[,apply(anti_micro_dados==0,2,sum)<5]
anti_micro=names(anti_micro_dados)

anti_micro_diff=anti_micro_dados[2:98,]-anti_micro_dados[1:97,]
diff_mat=as.dist(1-cor(anti_micro_diff))
cluster=hclust(diff_mat,method='average')
plot(cluster)

plot(anti_micro_dados$Pipetazo,anti_micro_dados$Cipro)
plot(anti_micro_diff$Pipetazo,anti_micro_diff$Cipro)


# Segundo ajuste


nivel=polynomial_block(1,values=1,D=1/0.9,name='nivel')
inter=polynomial_block(1,values=1:T>=indice_inter,D=1/1,k=1,name='intervenção')
covid=polynomial_block(1,values=dados$Entradas_COVID/dados$Entradas,D=1/1,k=1,name='COVID')
gravi=polynomial_block(1,values=dados$N_VM/dados$Pac_dia,D=1/1,k=1,name='gravidade')
modelo=block_merge(nivel,inter,covid,gravi)
for(name in anti_micro){
  bloco=polynomial_block(1,values=anti_micro_dados[[name]],D=1/1,name=paste0('anti_micro_',name),k=1,C0=1/mean(anti_micro_dados[[name]]**2))
  modelo=block_merge(modelo,bloco)
}
ajuste=fit_model(modelo,data_out=dados$N_OBT %>% as.matrix,offset=dados$Pac_dia %>% as.matrix,kernel='Poisson')

show_fit(ajuste,smooth=FALSE,t_offset = 1)$plot
show_fit(ajuste,smooth=TRUE)$plot

plot_lat_var(ajuste,'intervenção')$plot
plot_lat_var(ajuste,'COVID')$plot
plot_lat_var(ajuste,'gravidade')$plot


plot_lat_var(ajuste,'anti_micro')$plot



cov_matrix=ajuste$Cts[1:9,1:9,98]

n=dim(cov_matrix)[1]
labels=c('nivel','intervenção','COVID','gravidade',anti_micro)

for(i in c(1:n)){
  var=cov_matrix[i,i]
  cov_matrix[i,]=cov_matrix[i,]/sqrt(var)
  cov_matrix[,i]=cov_matrix[,i]/sqrt(var)
}

plot_tile_data=matrix(0,0,3)
for(i in c(1:n)){
  for(j in c(1:n)){
    if(i==j){
      plot_tile_data=rbind(plot_tile_data,c(i,j,NA))
    }else{
      plot_tile_data=rbind(plot_tile_data,c(i,j,cov_matrix[i,j]))
    }
  }
}


plot_tile_data=as.data.frame(plot_tile_data)
names(plot_tile_data)=c('X','Y','color')

plot_tile_data$X=as.numeric(plot_tile_data$X)
plot_tile_data$Y=as.numeric(plot_tile_data$Y)

ggplot(plot_tile_data)+
  geom_tile(aes(x=X,y=Y,fill=color))+
  geom_text(aes(x=X,y=Y,label=round(color,2)))+
  scale_fill_gradient2(low='red',mid='green',high='blue',limits=c(-1,1))+
  scale_x_continuous(labels=labels,breaks=c(1:n),expand=c(0,0))+
  scale_y_reverse(labels=labels,breaks=c(1:n),expand=c(0,0))+
  coord_equal()+
  theme_bw()
