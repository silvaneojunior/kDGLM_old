# library(ggplot2)
# library(plotly)
# library(latex2exp)

calcula_max=function(pre_max){
  if(length(pre_max)==0 | sum(pre_max**2)<10**-20){
    pre_max=1
  }else{
    pre_max=max(pre_max)
  }
  scaled_max=log10(pre_max)
  category=scaled_max%%1
  value=10**(floor(log10(max(pre_max))))
  if(category<0.1){
    value=value/10
  }else{
    if(category<0.25){
      value=value/5
    }else{
      if(category<0.5){
        value=value/2
      }
    }
  }
  interval_size=(pre_max%/%value)+2
  max_value=value*interval_size

  return(list(value,interval_size,max_value))
}

show_fit=function(model,IC_prob=0.95,smooth=TRUE,dinamic=TRUE,t_offset=0){
  n=dim(model$mt)[1]
  t_last=dim(model$mt)[2]
  eval=eval_past(model,smooth,t_offset)
  pred=eval$pred
  a=eval$a
  b=eval$b
  icl.pred<-qnbinom((1-IC_prob)/2, a, (b/(b +1)))
  icu.pred<-qnbinom(1-(1-IC_prob)/2, a, (b/(b +1)))

  fill_list=c('#2596be','#2596be','black')
  names(fill_list)=c(paste0('Prediction I.C. (',(100*IC_prob) %>% round(),'%)'),'Prediction','Observed values')
  color_list=c('#2596be','black')
  names(color_list)=c('Prediction','Observed values')

  max_value=calcula_max(model$data_out-min(model$data_out))[[3]]+min(model$data_out)
  min_value=-calcula_max(-(model$data_out-max(model$data_out)))[[3]]+max(model$data_out)

  plt=ggplot()+
    geom_point(aes(x=c(1:t_last),y=pred,color='Prediction',fill='Prediction'))+
    geom_ribbon(aes(x=c(1:t_last),ymin=icl.pred,ymax=icu.pred,fill=paste0('Prediction I.C. (',(100*IC_prob) %>% round(),'%)'),color=paste0('Prediction I.C. (',(100*IC_prob) %>% round(),'%)')),alpha=0.25)+
    geom_point(aes(x=c(1:t_last),y=model$data_out,color='Observed values',fill='Observed values'))+
    scale_fill_manual('',na.value=NA,values=fill_list)+
    scale_color_manual('',na.value=NA,values=color_list)+
    scale_y_continuous(name='$y_t$')+
    scale_x_continuous('Time')+
    theme_bw()+
    coord_cartesian(ylim=c(min_value,max_value))
  if(dinamic){
    plt=ggplotly(plt)
  }
  return(list('plot'=plt,'pred'=pred,'icl.pred'=icl.pred,'icu.pred'=icu.pred))
}

plot_lat_var=function(model,var,smooth=TRUE,cut_off=10,IC_prob=0.95,dinamic=TRUE,tranform_y=function(y){y}){
  if(!(var %in% names(model$names))){
    stop(paste0('Error: Invalid selected variable. Got ',var,', expected one of the following:\n',names(model$names)))
  }
  if(IC_prob>=1 | IC_prob<=0){
    stop(paste0('Error: Invalid value for I.C. width. Must be between 0 and 1, got ',IC_prob))
  }

  indice=model$names[[var]]
  size=length(indice)
  t=dim(model$mts)[2]
  m1=if(smooth){model$mts[indice,]}else{model$mt[indice,]}
  m1=m1 %>% matrix(size,t) %>% t
  std_mat=if(smooth){model$Cts[indice,indice,]}else{model$Ct[indice,indice,]}
  if(size>1){
    std_mat=std_mat %>% apply(3,diag)
  }
  std_mat=std_mat %>% sqrt %>% matrix(size,t) %>% t

  lim_i=m1+qnorm((1-IC_prob)/2)*std_mat
  lim_s=m1+qnorm(1-(1-IC_prob)/2)*std_mat

  m1=m1 %>% tranform_y
  lim_i=lim_i %>% tranform_y
  lim_s=lim_s %>% tranform_y

  m1=as.data.frame(m1)
  lim_i=as.data.frame(lim_i)
  lim_s=as.data.frame(lim_s)

  names(m1)=paste0('Variable ',c(1:dim(m1)[2]))
  names(lim_i)=paste0('Variable ',c(1:dim(lim_i)[2]))
  names(lim_s)=paste0('Variable ',c(1:dim(lim_s)[2]))

  max_value=calcula_max(m1-min(m1))[[3]]+min(m1)
  min_value=-calcula_max(-(m1-max(m1)))[[3]]+max(m1)

  m1$time=c(1:dim(m1)[1])
  lim_i$time=c(1:dim(lim_i)[1])
  lim_s$time=c(1:dim(lim_s)[1])

  m1=m1[-c(1:cut_off),] %>% pivot_longer(1:size) %>% rename(media=value)
  lim_i=lim_i[-c(1:cut_off),] %>% pivot_longer(1:size) %>% rename(lim_i=value)
  lim_s=lim_s[-c(1:cut_off),] %>% pivot_longer(1:size) %>% rename(lim_s=value)

  plot_data= m1 %>%
    inner_join(lim_i,by=c('time','name')) %>%
    inner_join(lim_s,by=c('time','name'))

  n_var=length(unique(plot_data$name))
  color_list=rainbow(n_var,s=0.5)
  names(color_list)=paste(unique(plot_data$name),'point estimate')

  fill_list=rainbow(n_var,s=0.5)
  names(fill_list)=paste0(unique(plot_data$name),' I.C. (',IC_prob*100 %>% round(),'%)')

  plt=ggplot(plot_data)+
    geom_hline(yintercept=0,linetype='dashed')+
    scale_x_continuous('Time')+
    scale_color_manual('',values=color_list,na.value=NA)+
    scale_fill_manual('',values=fill_list,na.value=NA)+
    labs(title=paste0(var,' (',ifelse(smooth,'smoothed','only filtered'),')'))+
    scale_y_continuous('Parameter value')+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90))+
    geom_ribbon(aes(x=time,ymin=lim_i,ymax=lim_s,fill=paste0(name,' I.C. (',IC_prob*100 %>% round(),'%)'),color=paste0(name,' I.C. (',IC_prob*100 %>% round(),'%)')),alpha=0.25)+
    geom_line(aes(x=time,y=media,color=paste(name,'point estimate'),fill=paste(name,'point estimate')))+
    coord_cartesian(ylim=c(min_value,max_value))
  if(dinamic){
    plt=ggplotly(plt)
  }
  plt
}
