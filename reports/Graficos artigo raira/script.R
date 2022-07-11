library(ggplot2)
library(plotly)
library(dplyr)
library(tidyr)

# Se o pacote GDLM não estiver instaldo
# devtools::install_github('silvaneojunior/GDLM')
library(GDLM)

# Tamanho da font
font_size=16

# Dados IBM
data=read.csv('data/m-ibmln.txt')

T = dim(data)[1]
scale=10
y <- (cbind(data,rep(0,T))*scale) %>% as.matrix
year <- seq(as.Date("1926/10/1"), by = "month", length.out = dim(y)[1])
year_label <- seq(as.Date("1930/1/1"), by = "10 years", length.out = 8)

level=polynomial_block(order=1,value=matrix(c(1,0),2,T),D=1/0.85,C0=scale**2)
volatility=polynomial_block(order=1,value=matrix(c(0,1),2,T),D=1/0.9,mu=-2*log(scale))

fitted_data=fit_model(level,volatility,data_out=y,kernel='Normal')

show_fit(fitted_data)$plot

ggplot()+
  geom_line(aes(x=year,y=fitted_data$mt[1,],linetype='On-line'))+
  geom_line(aes(x=year,y=fitted_data$mts[1,],linetype='Smoothed'))+
  geom_ribbon(aes(x=year,
                  ymin=fitted_data$mts[1,]-1.96*sqrt(fitted_data$Cts[1,1,]),
                  ymax=fitted_data$mts[1,]+1.96*sqrt(fitted_data$Cts[1,1,]),
                  linetype='Smoothed'),alpha=0.25)+
  scale_linetype_manual('',values=c('dotted','solid'))+
  scale_x_date('Year',breaks=year_label,labels=function(x){substr(x,1,4)})+
  scale_y_continuous('Estimated value')+
  theme_bw()+
  theme(text=element_text(size=font_size))

ggplot()+
  geom_line(aes(x=year,y=exp(fitted_data$mt[2,]),linetype='On-line'))+
  geom_line(aes(x=year,y=exp(fitted_data$mts[2,]),linetype='Smoothed'))+
  geom_ribbon(aes(x=year,
                  ymin=qlnorm(0.025,fitted_data$mts[2,],sqrt(fitted_data$Cts[2,2,])),
                  ymax=qlnorm(0.975,fitted_data$mts[2,],sqrt(fitted_data$Cts[2,2,])),
                  linetype='Smoothed'),alpha=0.25)+
  scale_linetype_manual('',values=c('dotted','solid'))+
  scale_x_date('Year',breaks=year_label,labels=function(x){substr(x,1,4)})+
  scale_y_continuous('Estimated value')+
  theme_bw()+
  theme(text=element_text(size=font_size))


# Dados Diarréia
library(readxl)
library(janitor)

rs_diarreia=read_excel("data/brasil_diarreia-2.xlsx",
                sheet = "regiao_sudeste_v2") %>%  clean_names()
xx = rs_diarreia$data
rs_diarreia <-
  rs_diarreia %>%
  mutate(x10_a_19_anos = x10_a_14_anos + x15_a_19_anos,
         x20_a_49_anos = x20_a_29_anos +  x30_a_39_anos  +  x40_a_49_anos,
         x50_anos_e_mais = x50_a_59_anos + x60_a_69_anos + x70_a_79_anos + x80_anos_e_mais,
         menor_que_4 = x1_a_4_anos + menor_1_ano,
         total_sem = total - x50_anos_e_mais - menor_que_4,
         outros = total - x50_anos_e_mais - menor_que_4)

y1 <- cbind(rs_diarreia$x50_anos_e_mais, rs_diarreia$menor_que_4, rs_diarreia$total_sem)

y <- trunc(y1)
T <- nrow(y)

level=polynomial_block(order=2,value=1,D=1/0.95,k=2)
season=harmonic_block(period=12,value=1,D=1/0.975,k=2)

fitted_data=fit_model(level,season,data_out=y,kernel='Multinomial')

predictions=show_fit(fitted_data,labels=c('Senior Group','Early childhood','Remaining groups'))

plot_data=predictions$pred %>%
  full_join(predictions$icl.pred, by=c('time','Serie')) %>%
  full_join(predictions$icu.pred, by=c('time','Serie'))


ggplot(plot_data[plot_data$Serie!='Remaining groups',])+
  geom_line(aes(x=time,y=Prediction,linetype=Serie))+
  geom_ribbon(aes(x=time,
                  ymin=I.C.lower,
                  ymax=I.C.upper,
                  linetype=Serie),alpha=0.25)+
  scale_linetype_manual('',values=c('dotted','solid'))+
  scale_y_continuous('One-step ahead prediction')+
  scale_x_continuous('Year',breaks=c(0:10)*12*5+1+24,labels =function(x){xx[x] %>% substr(1,4)},limits=c(25,T))+
  theme_bw()+
  theme(text=element_text(size=font_size))

total=rowSums(y)
for(i in 1:T){
  plot_data$Prediction[plot_data$time==i]=plot_data$Prediction[plot_data$time==i]/total[i]
  plot_data$I.C.lower[plot_data$time==i]=plot_data$I.C.lower[plot_data$time==i]/total[i]
  plot_data$I.C.upper[plot_data$time==i]=plot_data$I.C.upper[plot_data$time==i]/total[i]
}

ggplot(plot_data[plot_data$Serie!='Remaining groups',])+
  geom_line(aes(x=time,y=Prediction,linetype=Serie))+
  geom_ribbon(aes(x=time,
                  ymin=I.C.lower,
                  ymax=I.C.upper,
                  linetype=Serie),alpha=0.25)+
  scale_linetype_manual('',values=c('dotted','solid'))+
  scale_y_continuous('One-step ahead prediction',labels=function(x){round(100*x)%>%paste0('%')},breaks=c(0:4)*0.25,limits=c(0,1))+
  scale_x_continuous('Year',breaks=c(0:10)*12*5+1+24,labels =function(x){xx[x] %>% substr(1,4)},limits=c(25,T))+
  theme_bw()+
  theme(text=element_text(size=font_size))


ggplot()+
  geom_line(aes(x=25:T,y=fitted_data$mts[2,25:T]),linetype='solid')+
  geom_ribbon(aes(x=25:T,
                  ymin=fitted_data$mts[2,25:T]-1.96*sqrt(fitted_data$Cts[2,2,25:T]),
                  ymax=fitted_data$mts[2,25:T]+1.96*sqrt(fitted_data$Cts[2,2,25:T])),
              alpha=0.25)+
  geom_hline(linetype='dashed',yintercept = 0)+
  scale_y_continuous('Estimated value',limits=c(-0.02,0.02))+
  scale_x_continuous('Year',breaks=c(0:10)*12*5+1+24,labels =function(x){xx[x] %>% substr(1,4)})+
  theme_bw()+
  theme(text=element_text(size=16))

ggplot()+
  geom_line(aes(x=25:T,y=fitted_data$mts[4,25:T]),linetype='dashed')+
  geom_ribbon(aes(x=25:T,
                  ymin=fitted_data$mts[4,25:T]-1.96*sqrt(fitted_data$Cts[4,4,25:T]),
                  ymax=fitted_data$mts[4,25:T]+1.96*sqrt(fitted_data$Cts[4,4,25:T])),
              alpha=0.25)+
  geom_hline(linetype='dashed',yintercept = 0)+
  scale_y_continuous('Estimated value',limits=c(-0.02,0.02))+
  scale_x_continuous('Year',breaks=c(0:10)*12*5+1+24,labels =function(x){xx[x] %>% substr(1,4)})+
  theme_bw()+
  theme(text=element_text(size=font_size))
