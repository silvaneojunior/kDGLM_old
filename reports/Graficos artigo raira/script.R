library(tidyverse)
devtools::load_all()
source('reports/Graficos artigo raira/plot_zoom.R')

library(extrafont)
loadfonts(device = "win")

# Tamanho da fonte
font_size=12
family_font='serif'
mult_scale=2

base.h=540
base.w=960

# Dados IBM
data=read.csv('data/m-ibmln.txt')

T = dim(data)[1]
y <- data
year <- seq(as.Date("1926/10/1"), by = "month", length.out = dim(y)[1])
year_label <- seq(as.Date("1930/1/1"), by = "10 years", length.out = 8)

level=polynomial_block(order=1,value=c(1,0),D=1/0.98,by_time = F)
volatility1=polynomial_block(order=1,value=c(0,1),D=1/1,by_time = F)
volatility2=AR_block(order=1,value=c(0,1),D=diag(c(1,1/0.98)),W=diag(c(0.1,0)),by_time = F)

fitted_data=fit_model(level,volatility1,volatility2,outcome=y[,1],family='Normal_Gamma_cor')

acf(fitted_data$mts[3,])
fitted_data$mts[4,T]
arima(fitted_data$mts[3,],c(1,0,0),include.mean=FALSE,transform.pars = FALSE)

show_fit(fitted_data,smooth=FALSE)$plot

pred_data=show_fit(fitted_data,smooth=FALSE)$data

zoom_factor=2

(ggplot()+
    geom_point(aes(x=year,y=pred_data$Observation,shape='Observations'))+
    geom_line(aes(x=year,y=pred_data$Prediction,linetype='Predictive mean'))+
    geom_ribbon(aes(x=year,
                    ymin=pred_data$C.I.lower,
                    ymax=pred_data$C.I.upper,
                    linetype='C.I. 95%'),alpha=0.25)+
    scale_linetype_manual('',values=c('dotdash','solid'))+
    scale_shape('')+
    scale_x_date('Year',breaks=year_label,labels=function(x){substr(x,1,4)},expand=c(0,0))+
    scale_y_continuous('Estimated value')+
    coord_cartesian(ylim=c(-0.5,0.5))+
    theme_bw()+
    theme(text=element_text(size=font_size,family=family_font)))+
  labs(title='One-step ahead prediction')

ggsave(
  'C:\\Jupyter\\Mestrado\\Pacote\\GDLM\\reports\\Graficos artigo raira\\NG-1-step-pred.png',
  scale = 1,
  units='px',
  width=base.w*mult_scale,
  height=base.h*mult_scale,
  dpi = 300
)


((ggplot()+
  geom_line(aes(x=year,y=fitted_data$mt[1,],linetype='On-line'))+
  geom_line(aes(x=year,y=fitted_data$mts[1,],linetype='Smoothed'))+
  geom_ribbon(aes(x=year,
                  ymin=fitted_data$mts[1,]-1.96*sqrt(fitted_data$Cts[1,1,]),
                  ymax=fitted_data$mts[1,]+1.96*sqrt(fitted_data$Cts[1,1,]),
                  linetype='Smoothed'),alpha=0.25)+
  scale_linetype_manual('',values=c('dotdash','solid'))+
  scale_x_date('Year',breaks=year_label,labels=function(x){substr(x,1,4)},expand=c(0,0))+
  scale_y_continuous('Estimated value')+
  coord_cartesian(ylim=c(-0.2,0.1))+
  theme_bw()+
  theme(text=element_text(size=font_size,family=family_font))) %>% plot_zoom(mini_range=list('x'=c(as.Date("1990-01-01"),as.Date("1998-12-31")),'y'=c(-0.1,0.1)),
                                                         mini_place=list('x'=c(as.Date("1940-01-01"),as.Date(paste0("1980-12-31"))),'y'=c(-0.2,-0.05))))+
  labs(title='Mean estimation')

ggsave(
  'C:\\Jupyter\\Mestrado\\Pacote\\GDLM\\reports\\Graficos artigo raira\\NG-mean-dual.png',
  scale = 1,
  units='px',
  width=base.w*mult_scale,
  height=base.h*mult_scale,
  dpi = 300
)

(ggplot()+
    geom_line(aes(x=year,y=fitted_data$mt[1,],linetype='On-line'))+
    geom_line(aes(x=year,y=fitted_data$mts[1,],linetype='Smoothed'))+
    geom_ribbon(aes(x=year,
                    ymin=fitted_data$mts[1,]-1.96*sqrt(fitted_data$Cts[1,1,]),
                    ymax=fitted_data$mts[1,]+1.96*sqrt(fitted_data$Cts[1,1,]),
                    linetype='Smoothed'),alpha=0.25)+
    scale_linetype_manual('',values=c('dotdash','solid'))+
    scale_x_date('Year',breaks=year_label,labels=function(x){substr(x,1,4)},expand=c(0,0))+
    scale_y_continuous('Estimated value')+
    coord_cartesian(ylim=c(-0.15,0.15))+
    theme_bw()+
    theme(text=element_text(size=font_size,family=family_font)))+
  labs(title='Mean estimation')

ggsave(
  'C:\\Jupyter\\Mestrado\\Pacote\\GDLM\\reports\\Graficos artigo raira\\NG-mean-full.png',
  scale = 1,
  units='px',
  width=base.w*mult_scale,
  height=base.h*mult_scale,
  dpi = 300
)


ggplot()+
  geom_line(aes(x=year,y=fitted_data$mt[1,],linetype='On-line'))+
  geom_line(aes(x=year,y=fitted_data$mts[1,],linetype='Smoothed'))+
  geom_ribbon(aes(x=year,
                  ymin=fitted_data$mts[1,]-1.96*sqrt(fitted_data$Cts[1,1,]),
                  ymax=fitted_data$mts[1,]+1.96*sqrt(fitted_data$Cts[1,1,]),
                  linetype='Smoothed'),alpha=0.25)+
  scale_linetype_manual('',values=c('dotdash','solid'))+
  coord_cartesian(ylim=c(-0.1,0.1))+
  scale_x_date('Year',breaks=seq(as.Date("1930/1/1"), by = "1 years", length.out = 70),labels=function(x){substr(x,1,4)},expand=c(0,0),limits=c(as.Date("1990-01-01"),NA))+
  scale_y_continuous('Estimated value')+
  labs(title='Mean estimation')+
  theme_bw()+
  theme(text=element_text(size=font_size,family=family_font))

ggsave(
  'C:\\Jupyter\\Mestrado\\Pacote\\GDLM\\reports\\Graficos artigo raira\\NG-mean-zoomed.png',
  scale = 1,
  units='px',
  width=base.w*mult_scale,
  height=base.h*mult_scale,
  dpi = 300
)

ggplot()+
  geom_line(aes(x=year,y=fitted_data$mt[4,],linetype='On-line'))+
  geom_line(aes(x=year,y=fitted_data$mts[4,],linetype='Smoothed'))+
  geom_ribbon(aes(x=year,
                  ymin=fitted_data$mts[4,]-1.96*sqrt(fitted_data$Cts[4,4,]),
                  ymax=fitted_data$mts[4,]+1.96*sqrt(fitted_data$Cts[4,4,]),
                  linetype='Smoothed'),alpha=0.25)+
  scale_linetype_manual('',values=c('dotdash','solid'))+
  coord_cartesian(ylim=c(0,1.5))+
  scale_x_date('Year',breaks=seq(as.Date("1930/1/1"), by = "1 years", length.out = 70),labels=function(x){substr(x,1,4)},expand=c(0,0),limits=c(as.Date("1990-01-01"),NA))+
  scale_y_continuous('Estimated value')+
  labs(title='AR coefficient')+
  theme_bw()+
  theme(text=element_text(size=font_size,family=family_font))

ggsave(
  'C:\\Jupyter\\Mestrado\\Pacote\\GDLM\\reports\\Graficos artigo raira\\NG-AR.png',
  scale = 1,
  units='px',
  width=base.w*mult_scale,
  height=base.h*mult_scale,
  dpi = 300
)

fts=matrix(NA,2,T)
Qts=array(NA,c(2,2,T))
for(i in 1:T){
  fts[,i]=t(fitted_data$FF[,,i])%*%fitted_data$mts[,i]
  Qts[,,i]=t(fitted_data$FF[,,i])%*%fitted_data$Cts[,,i]%*%fitted_data$FF[,,i]
}

((ggplot()+
    geom_line(aes(x=year,y=exp(fitted_data$ft[2,]+fitted_data$Qt[2,2,]/2),linetype='On-line'))+
    geom_line(aes(x=year,y=exp(fts[2,]+Qts[2,2,]/2),linetype='Smoothed'))+
    geom_ribbon(aes(x=year,
                    ymin=qlnorm(0.025,fts[2,],sqrt(Qts[2,2,])),
                    ymax=qlnorm(0.975,fts[2,],sqrt(Qts[2,2,])),
                    linetype='Smoothed'),alpha=0.25)+
    scale_linetype_manual('',values=c('dotdash','solid'))+
    scale_x_date('Year',breaks=year_label,labels=function(x){substr(x,1,4)},expand=c(0,0))+
    coord_cartesian(ylim=c(0,5000))+
    scale_y_continuous('Estimated value',expand=c(0,0),breaks=seq(0,5000,500),labels=function(x){formatC(x,big.mark='.')})+
    theme_bw()+
    theme(text=element_text(size=font_size,family=family_font))) %>% plot_zoom(mini_range=list('x'=c(as.Date("1990-01-01"),as.Date("1998-12-31")),'y'=c(0,1000)),
                                                            mini_place=list('x'=c(as.Date("1940-01-01"),as.Date(paste0("1980-12-31"))),'y'=c(2500,4500))))+
  labs(title='Precision estimation')

# ((ggplot()+
#   geom_line(aes(x=year,y=exp(-fitted_data$ft[2,]+fitted_data$Qt[2,2,]/2),linetype='On-line'))+
#   geom_line(aes(x=year,y=exp(-fts[2,]+Qts[2,2,]/2),linetype='Smoothed'))+
#   geom_ribbon(aes(x=year,
#                   ymin=qlnorm(0.025,-fts[2,],sqrt(Qts[2,2,])),
#                   ymax=qlnorm(0.975,-fts[2,],sqrt(Qts[2,2,])),
#                   linetype='Smoothed'),alpha=0.25)+
#   scale_linetype_manual('',values=c('dotdash','solid'))+
#   scale_x_date('Year',breaks=year_label,labels=function(x){substr(x,1,4)},expand=c(0,0))+
#   coord_cartesian(ylim=c(NA,0.03))+
#   scale_y_continuous('Estimated value',expand=c(0,0),breaks=seq(0,5000,500),labels=function(x){formatC(x,big.mark='.')})+
#   theme_bw()+
#   theme(text=element_text(size=font_size,family=family_font))))+
#   labs(title='Variance estimation')

ggsave(
  'C:\\Jupyter\\Mestrado\\Pacote\\GDLM\\reports\\Graficos artigo raira\\NG-volatility-dual.png',
  scale = 1,
  units='px',
  width=base.w*mult_scale,
  height=base.h*mult_scale,
  dpi = 300
)

((ggplot()+
    geom_line(aes(x=year,y=exp(fitted_data$ft[2,]+fitted_data$Qt[2,2,]/2),linetype='On-line'))+
    geom_line(aes(x=year,y=exp(fts[2,]+Qts[2,2,]/2),linetype='Smoothed'))+
    geom_ribbon(aes(x=year,
                    ymin=qlnorm(0.025,fts[2,],sqrt(Qts[2,2,])),
                    ymax=qlnorm(0.975,fts[2,],sqrt(Qts[2,2,])),
                    linetype='Smoothed'),alpha=0.25)+
    scale_linetype_manual('',values=c('dotdash','solid'))+
    scale_x_date('Year',breaks=year_label,labels=function(x){substr(x,1,4)},expand=c(0,0))+
    coord_cartesian(ylim=c(0,3000))+
    scale_y_continuous('Estimated value',expand=c(0,0),breaks=seq(0,3000,500),labels=function(x){formatC(x,big.mark='.')})+
    theme_bw()+
    theme(text=element_text(size=font_size,family=family_font))))+
  labs(title='Precision estimation')

ggsave(
  'C:\\Jupyter\\Mestrado\\Pacote\\GDLM\\reports\\Graficos artigo raira\\NG-volatility-full.png',
  scale = 1,
  units='px',
  width=base.w*mult_scale,
  height=base.h*mult_scale,
  dpi = 300
)

((ggplot()+
    geom_line(aes(x=year,y=exp(fitted_data$ft[2,]+fitted_data$Qt[2,2,]/2),linetype='On-line'))+
    geom_line(aes(x=year,y=exp(fts[2,]+Qts[2,2,]/2),linetype='Smoothed'))+
    geom_ribbon(aes(x=year,
                    ymin=qlnorm(0.025,fts[2,],sqrt(Qts[2,2,])),
                    ymax=qlnorm(0.975,fts[2,],sqrt(Qts[2,2,])),
                    linetype='Smoothed'),alpha=0.25)+
    scale_linetype_manual('',values=c('dotdash','solid'))+
    scale_x_date('Year',breaks=seq(as.Date("1930/1/1"), by = "1 years", length.out = 70),labels=function(x){substr(x,1,4)},expand=c(0,0),limits=c(as.Date("1990-01-01"),NA))+
    coord_cartesian(ylim=c(0,1000))+
    scale_y_continuous('Estimated value',expand=c(0,0),breaks=seq(0,5000,250),labels=function(x){formatC(x,big.mark='.')})+
    theme_bw()+
    theme(text=element_text(size=font_size,family=family_font))))+
  labs(title='Precision estimation')

ggsave(
  'C:\\Jupyter\\Mestrado\\Pacote\\GDLM\\reports\\Graficos artigo raira\\NG-volatility-zoomed.png',
  scale = 1,
  units='px',
  width=base.w*mult_scale,
  height=base.h*mult_scale,
  dpi = 300
)


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

fitted_data=fit_model(level,season,outcome=y,family='Multinomial')
show_fit(fitted_data,labels=c('Senior Group','Early childhood','Remaining groups'))$plot

predictions=show_fit(fitted_data,labels=c('Senior Group','Early childhood','Remaining groups'),smooth=FALSE)$data

plot_data=predictions %>% filter(Serie!='Remaining groups')


((ggplot()+
  geom_line(aes(x=Time,linetype=Serie,y=Prediction),data=plot_data)+
  geom_point(aes(x=Time,y=Observation,shape=Serie),data=plot_data)+
  geom_ribbon(aes(x=Time,linetype=Serie,
                  ymin=C.I.lower,
                  ymax=C.I.upper),alpha=0.25,data=plot_data)+
  scale_linetype_manual('Prediction',values=c('dotdash','solid'))+
  scale_shape('Observation')+
  scale_y_continuous('Hospital admissions',limits=c(0,4500),expand=c(0,0),labels=function(x){formatC(x,big.mark='.')})+
  scale_x_continuous('Year',breaks=c(0:10)*12*5+1+24,labels =function(x){xx[x] %>% substr(1,4)},limits=c(25,T),expand=c(0,0))+
  theme_bw()+
  theme(text=element_text(size=font_size,family=family_font))) %>% plot_zoom(mini_range=list('x'=c(17*12+1,22*12+1),'y'=c(250,1250)),
                                                          mini_place=list('x'=c(100,250),'y'=c(2250,4250))))+
  labs(title='One-step ahead prediction')

ggsave(
  'C:\\Jupyter\\Mestrado\\Pacote\\GDLM\\reports\\Graficos artigo raira\\Multinom-pred-dual-obs.png',
  scale = 1,
  units='px',
  width=base.w*mult_scale,
  height=base.h*mult_scale,
  dpi = 300
)

(ggplot()+
    geom_line(aes(x=Time,linetype=Serie,y=Prediction),data=plot_data)+
    geom_point(aes(x=Time,y=Observation,shape=Serie),data=plot_data)+
    geom_ribbon(aes(x=Time,linetype=Serie,
                    ymin=C.I.lower,
                    ymax=C.I.upper),alpha=0.25,data=plot_data)+
    scale_linetype_manual('Prediction',values=c('dotdash','solid'))+
    scale_shape('Observation')+
    scale_y_continuous('Hospital admissions',limits=c(250,1250),expand=c(0,0),labels=function(x){formatC(x,big.mark='.')})+
    scale_x_continuous('Year',breaks=c(0:10)*12*5+1+24,labels =function(x){xx[x] %>% substr(1,4)},limits=c(17*12+1,22*12+1),expand=c(0,0))+
    theme_bw()+
    theme(text=element_text(size=font_size,family=family_font)))+
  labs(title='One-step ahead prediction')

ggsave(
  'C:\\Jupyter\\Mestrado\\Pacote\\GDLM\\reports\\Graficos artigo raira\\Multinom-pred-zoomed-obs.png',
  scale = 1,
  dpi = 300,
  units='px',
  width=base.w*mult_scale,
  height=base.h*mult_scale
)

(ggplot()+
    geom_line(aes(x=Time,linetype=Serie,y=Prediction),data=plot_data)+
    geom_point(aes(x=Time,y=Observation,shape=Serie),data=plot_data)+
    geom_ribbon(aes(x=Time,linetype=Serie,
                    ymin=C.I.lower,
                    ymax=C.I.upper),alpha=0.25,data=plot_data)+
    scale_linetype_manual('Prediction',values=c('dotdash','solid'))+
    scale_shape('Observation')+
    scale_y_continuous('Hospital admissions',limits=c(0,2500),expand=c(0,0),labels=function(x){formatC(x,big.mark='.')})+
    scale_x_continuous('Year',breaks=c(0:10)*12*5+1+24,labels =function(x){xx[x] %>% substr(1,4)},limits=c(25,T),expand=c(0,0))+
    theme_bw()+
    theme(text=element_text(size=font_size,family=family_font)))+
  labs(title='One-step ahead prediction')

ggsave(
  'C:\\Jupyter\\Mestrado\\Pacote\\GDLM\\reports\\Graficos artigo raira\\Multinom-pred-full-obs.png',
  scale = 1,
  dpi = 300,
  units='px',
  width=base.w*mult_scale,
  height=base.h*mult_scale
)

#################################### Sem observações ############################################

((ggplot()+
    geom_line(aes(x=Time,linetype=Serie,y=Prediction),data=plot_data)+
    geom_ribbon(aes(x=Time,linetype=Serie,
                    ymin=C.I.lower,
                    ymax=C.I.upper),alpha=0.25,data=plot_data)+
    scale_linetype_manual('Prediction',values=c('dotdash','solid'))+
    scale_shape('Observation')+
    scale_y_continuous('Hospital admissions',limits=c(0,4500),expand=c(0,0),labels=function(x){formatC(x,big.mark='.')})+
    scale_x_continuous('Year',breaks=c(0:10)*12*5+1+24,labels =function(x){xx[x] %>% substr(1,4)},limits=c(25,T),expand=c(0,0))+
    theme_bw()+
    theme(text=element_text(size=font_size,family=family_font))) %>% plot_zoom(mini_range=list('x'=c(17*12+1,22*12+1),'y'=c(250,1250)),
                                                            mini_place=list('x'=c(100,250),'y'=c(2250,4250))))+
  labs(title='One-step ahead prediction')

ggsave(
  'C:\\Jupyter\\Mestrado\\Pacote\\GDLM\\reports\\Graficos artigo raira\\Multinom-pred-dual.png',
  scale = 1,
  units='px',
  width=base.w*mult_scale,
  height=base.h*mult_scale,
  dpi = 300
)

(ggplot()+
    geom_line(aes(x=Time,linetype=Serie,y=Prediction),data=plot_data)+
    geom_ribbon(aes(x=Time,linetype=Serie,
                    ymin=C.I.lower,
                    ymax=C.I.upper),alpha=0.25,data=plot_data)+
    scale_linetype_manual('Prediction',values=c('dotdash','solid'))+
    scale_shape('Observation')+
    scale_y_continuous('Hospital admissions',limits=c(250,1250),expand=c(0,0),labels=function(x){formatC(x,big.mark='.')})+
    scale_x_continuous('Year',breaks=c(0:10)*12*5+1+24,labels =function(x){xx[x] %>% substr(1,4)},limits=c(17*12+1,22*12+1),expand=c(0,0))+
    theme_bw()+
    theme(text=element_text(size=font_size,family=family_font)))+
  labs(title='One-step ahead prediction')

ggsave(
  'C:\\Jupyter\\Mestrado\\Pacote\\GDLM\\reports\\Graficos artigo raira\\Multinom-pred-zoomed.png',
  scale = 1,
  dpi = 300,
  units='px',
  width=base.w*mult_scale,
  height=base.h*mult_scale
)

(ggplot()+
    geom_line(aes(x=Time,linetype=Serie,y=Prediction),data=plot_data)+
    geom_ribbon(aes(x=Time,linetype=Serie,
                    ymin=C.I.lower,
                    ymax=C.I.upper),alpha=0.25,data=plot_data)+
    scale_linetype_manual('Prediction',values=c('dotdash','solid'))+
    scale_shape('Observation')+
    scale_y_continuous('Hospital admissions',limits=c(0,2500),expand=c(0,0),labels=function(x){formatC(x,big.mark='.')})+
    scale_x_continuous('Year',breaks=c(0:10)*12*5+1+24,labels =function(x){xx[x] %>% substr(1,4)},limits=c(25,T),expand=c(0,0))+
    theme_bw()+
    theme(text=element_text(size=font_size,family=family_font)))+
  labs(title='One-step ahead prediction')

ggsave(
  'C:\\Jupyter\\Mestrado\\Pacote\\GDLM\\reports\\Graficos artigo raira\\Multinom-pred-full.png',
  scale = 1,
  dpi = 300,
  units='px',
  width=base.w*mult_scale,
  height=base.h*mult_scale
)

total=rowSums(y)
for(i in 1:T){
  plot_data$Observation[plot_data$Time==i]=plot_data$Observation[plot_data$Time==i]/total[i]
  plot_data$Prediction[plot_data$Time==i]=plot_data$Prediction[plot_data$Time==i]/total[i]
  plot_data$C.I.lower[plot_data$Time==i]=plot_data$C.I.lower[plot_data$Time==i]/total[i]
  plot_data$C.I.upper[plot_data$Time==i]=plot_data$C.I.upper[plot_data$Time==i]/total[i]
}

ggplot()+
  geom_line(aes(x=Time,linetype=Serie,y=Prediction),data=plot_data)+
  geom_point(aes(x=Time,y=Observation,shape=Serie),data=plot_data)+
  geom_ribbon(aes(x=Time,linetype=Serie,
                  ymin=C.I.lower,
                  ymax=C.I.upper),alpha=0.25,data=plot_data)+
  scale_linetype_manual('Prediction',values=c('dotdash','solid'))+
  scale_shape('Observation')+
  scale_y_continuous('Proportion of total admissions',labels=function(x){round(100*x)%>%paste0('%')},breaks=c(0:4)*0.2,limits=c(0,0.8),expand=c(0,0))+
  scale_x_continuous('Year',breaks=c(0:10)*12*5+1+24,labels =function(x){xx[x] %>% substr(1,4)},limits=c(25,T),expand=c(0,0))+
  theme_bw()+
  theme(text=element_text(size=font_size,family=family_font))+
  labs(title='One-step ahead prediction')

ggsave(
  'C:\\Jupyter\\Mestrado\\Pacote\\GDLM\\reports\\Graficos artigo raira\\Multinom-pred-comp-obs.png',
  scale = 1,
  dpi = 300,
  units='px',
  width=base.w*mult_scale,
  height=base.h*mult_scale
)

ggplot()+
  geom_line(aes(x=Time,linetype=Serie,y=Prediction),data=plot_data)+
  geom_ribbon(aes(x=Time,linetype=Serie,
                  ymin=C.I.lower,
                  ymax=C.I.upper),alpha=0.25,data=plot_data)+
  scale_linetype_manual('Prediction',values=c('dotdash','solid'))+
  scale_shape('Observation')+
  scale_y_continuous('Proportion of total admissions',labels=function(x){round(100*x)%>%paste0('%')},breaks=c(0:4)*0.2,limits=c(0,0.8),expand=c(0,0))+
  scale_x_continuous('Year',breaks=c(0:10)*12*5+1+24,labels =function(x){xx[x] %>% substr(1,4)},limits=c(25,T),expand=c(0,0))+
  theme_bw()+
  theme(text=element_text(size=font_size,family=family_font))+
  labs(title='One-step ahead prediction')

ggsave(
  'C:\\Jupyter\\Mestrado\\Pacote\\GDLM\\reports\\Graficos artigo raira\\Multinom-pred-comp.png',
  scale = 1,
  dpi = 300,
  units='px',
  width=base.w*mult_scale,
  height=base.h*mult_scale
)

var_index=1

data_plot=rbind(
                data.frame(Time=25:T,
                           var='Level',
                           values=fitted_data$mts[var_index,25:T],
                           c.i.lower=fitted_data$mts[var_index,25:T]-1.96*sqrt(fitted_data$Cts[var_index,var_index,25:T]),
                           c.i.upper=fitted_data$mts[var_index,25:T]+1.96*sqrt(fitted_data$Cts[var_index,var_index,25:T]),
                           age='Senior Group'),
                data.frame(Time=25:T,
                           var='Level',
                           values=fitted_data$mts[var_index+2,25:T],
                           c.i.lower=fitted_data$mts[var_index+2,25:T]-1.96*sqrt(fitted_data$Cts[var_index+2,var_index+2,25:T]),
                           c.i.upper=fitted_data$mts[var_index+2,25:T]+1.96*sqrt(fitted_data$Cts[var_index+2,var_index+2,25:T]),
                           age='Early childhood'))

ggplot(data_plot)+
  geom_line(aes(x=Time,y=values,linetype=age))+
  geom_ribbon(aes(x=Time,
                  ymin=c.i.lower,
                  ymax=c.i.upper,
                  group=age),
              alpha=0.25)+
  geom_hline(linetype='dashed',yintercept = 0)+
  scale_linetype('Age group')+
  scale_y_continuous('Estimated value',expand=c(0,0))+
  scale_x_continuous('Year',breaks=c(0:10)*12*5+1+24,labels =function(x){xx[x] %>% substr(1,4)},expand=c(0,0))+
  theme_bw()+
  theme(text=element_text(size=font_size,family=family_font))

ggsave(
  'C:\\Jupyter\\Mestrado\\Pacote\\GDLM\\reports\\Graficos artigo raira\\Multinom-vars1A.png',
  scale = 1,
  dpi = 300,
  units='px',
  width=base.w*mult_scale,
  height=base.h*mult_scale
)

var_index=2

data_plot=rbind(
                data.frame(Time=25:T,
                           var='Trend',
                           values=fitted_data$mts[var_index,25:T],
                           c.i.lower=fitted_data$mts[var_index,25:T]-1.96*sqrt(fitted_data$Cts[var_index,var_index,25:T]),
                           c.i.upper=fitted_data$mts[var_index,25:T]+1.96*sqrt(fitted_data$Cts[var_index,var_index,25:T]),
                           age='Senior Group'),
                data.frame(Time=25:T,
                           var='Trend',
                           c.i.lower=fitted_data$mts[var_index+2,25:T]-1.96*sqrt(fitted_data$Cts[var_index+2,var_index+2,25:T]),
                           c.i.upper=fitted_data$mts[var_index+2,25:T]+1.96*sqrt(fitted_data$Cts[var_index+2,var_index+2,25:T]),
                           values=fitted_data$mts[var_index+2,25:T],
                           age='Early childhood'))

ggplot(data_plot)+
  geom_line(aes(x=Time,y=values,linetype=age))+
  geom_ribbon(aes(x=Time,
                  ymin=c.i.lower,
                  ymax=c.i.upper,
                  group=age),
              alpha=0.25)+
  geom_hline(linetype='dashed',yintercept = 0)+
  scale_linetype('Age group')+
  scale_y_continuous('Estimated value',expand=c(0,0))+
  scale_x_continuous('Year',breaks=c(0:10)*12*5+1+24,labels =function(x){xx[x] %>% substr(1,4)},expand=c(0,0))+
  theme_bw()+
  theme(text=element_text(size=font_size,family=family_font))

ggsave(
  'C:\\Jupyter\\Mestrado\\Pacote\\GDLM\\reports\\Graficos artigo raira\\Multinom-vars1B.png',
  scale = 1,
  dpi = 300,
  units='px',
  width=base.w*mult_scale,
  height=base.h*mult_scale
)

var_index=5

data_plot=rbind(
                data.frame(Time=25:T,
                           var='Seasonality',
                           values=fitted_data$mts[var_index,25:T],
                           c.i.lower=fitted_data$mts[var_index,25:T]-1.96*sqrt(fitted_data$Cts[var_index,var_index,25:T]),
                           c.i.upper=fitted_data$mts[var_index,25:T]+1.96*sqrt(fitted_data$Cts[var_index,var_index,25:T]),
                           age='Senior Group'),
                data.frame(Time=25:T,
                           var='Seasonality',
                           values=fitted_data$mts[var_index+2,25:T],
                           c.i.lower=fitted_data$mts[var_index+2,25:T]-1.96*sqrt(fitted_data$Cts[var_index+2,var_index+2,25:T]),
                           c.i.upper=fitted_data$mts[var_index+2,25:T]+1.96*sqrt(fitted_data$Cts[var_index+2,var_index+2,25:T]),
                           age='Early childhood'))

ggplot(data_plot)+
  geom_line(aes(x=Time,y=values,linetype=age))+
  geom_ribbon(aes(x=Time,
                  ymin=c.i.lower,
                  ymax=c.i.upper),
              alpha=0.25)+
  geom_hline(linetype='dashed',yintercept = 0)+
  scale_linetype('Age group')+
  scale_y_continuous('Estimated value',expand=c(0,0),limits=c(-0.8,0.8))+
  scale_x_continuous('Year',breaks=c(0:10)*12*5+1+24,labels =function(x){xx[x] %>% substr(1,4)},expand=c(0,0))+
  theme_bw()+
  theme(text=element_text(size=font_size,family=family_font))+
  facet_wrap(.~age,ncol=1)

ggsave(
  'C:\\Jupyter\\Mestrado\\Pacote\\GDLM\\reports\\Graficos artigo raira\\Multinom-vars2.png',
  scale = 1,
  dpi = 300,
  units='px',
  width=base.w*mult_scale,
  height=base.h*mult_scale*1.5
)


