# library(Matrix)   # Pode ser substituido pela bilbioteca Matrix.
# library(tidyr)
# library(dplyr)

source('R/plot_helper.R')
source('R/structure_helper.R')
source('R/kernels.R')

ajusta_modelo <- function(...,data_out,kernel='Poisson (univariada)',offset=NULL,log_offset=NULL){
  if(typeof(kernel)==typeof('kernel')){
    kernel=kernel_list[[kernel]]
  }

  structure=concat_bloco(...)
  if(is.null(dim(data_out))){
    data_out=array(data_out,c(length(data_out),1))
  }
  if(!is.null(offset) & !is.null(log_offset)){
    stop('Erro: Cannot set both offset and log_offset. Choose only one.')
  }else{if(!is.null(offset)){
    log_offset=log(offset)
  }else{if(!is.null(log_offset)){
    offset=exp(log_offset)
  }else{
    offset=1
    log_offset=0
  }}}
  if(1==length(log_offset)){
    log_offset=rep(log_offset,dim(data_out)[1])
    offset=rep(offset,dim(data_out)[1])
  }
  if(dim(data_out)[1]!=length(log_offset)){
    stop('Erro: offset/log_offset does not have the same length as data_out.')
  }

  if(structure$t==1){
    structure$t=dim(data_out)[1]
    structure$D=array(structure$D,c(structure$n,structure$n,structure$t))
    structure$W=array(structure$W,c(structure$n,structure$n,structure$t))
    structure$FF=array(structure$FF,c(structure$n,structure$k,structure$t))
  }
  if(dim(data_out)[1]!=structure$t){
    stop(paste0('Erro: data_out does not have the same time length as structure: got ',dim(data_out)[1],' from data_out, expected ',structure$t))
  }
  #if(dim(data_out)[2]!=structure$k){
  #  stop(paste0('Erro: data_out does not have the same dimensions as structure: got ',dim(data_out)[2],' from data_out, expected ',structure$k))
  #}
  model=kernel$fit(y=data_out,
               m0=structure$m0,
               C0=structure$C0,
               FF=structure$FF,
               G=structure$G,
               D=structure$D,
               W=structure$W,
               pop=offset)

  model$m0=structure$m0
  model$C0=structure$C0
  model$names=structure$names
  model$kernel=kernel

  return(model)

}

predict=function(model,t=1,offset=NULL,log_offset=NULL,FF=NULL,D=NULL,W=NULL,plot=TRUE,IC_prob=0.95,labels=NULL){
  n=dim(model$mt)[1]
  t_last=dim(model$mt)[2]
  r=dim(model$FF)[2]

  if(t>10){
    warning('Warning: Prediction window is big, results will probabily be unreliable.')
  }

  #### Consistency check ####
  if(is.null(FF)){
    FF=array(model$FF[,,t_last],c(n,r,t))
  }
  if(is.null(D)){
    D=array(model$D[,,t_last],c(n,n,t))
    D[,,-1]=1
  }else{if(all(dim(D)==1)){
    D=array(D,c(n,n,t))
  }else{if(length(dim(D))==2 | (length(dim(D))==3 & dim(D)[3]==1)){
    D=array(D,c(n,n,t))
    D[,,-1]=1
  }}
  }
  if(is.null(W)){
    W=array(model$W[,,t_last],c(n,n,t))
    W[,,-1]=0
  }else{if(all(dim(W)==1)){
    W=array(diag(n)*W,c(n,n,t))
  }else{if(length(dim(W))==2 | (length(dim(W))==3 & dim(W)[3]==1)){
    W=array(W,c(n,n,t))
    W[,,-1]=0
  }}
  # if(t>1){
  #   D[,,2:t]=0
  # }
  # if(t>1){
  #   W[,,2:t]=0
  # }
  }
  if(!is.null(offset) & !is.null(log_offset)){
    stop('Erro: Cannot set both offset and log_offset. Choose only one.')
  }else{if(!is.null(offset)){
    log_offset=log(offset)
  }else{if(!is.null(log_offset)){
    offset=exp(log_offset)
  }else{
    offset=1
    log_offset=0
  }}}
  if(1==length(log_offset)){
    log_offset=rep(log_offset,t)
    offset=rep(offset,t)
  }
  print('prediction')

  if(dim(FF)[3]!=t){
    stop(paste0('Error: FF should have one matrix for each time or exactly 1 matrix, got ',dim(FF)[3],'!=',t,'.'))
  }
  if(dim(FF)[2]!=r){
    stop(paste0('Error: FF should have one column for each serie or exactly 1 column, got ',dim(FF)[2],'!=',r,'.'))
  }
  if(dim(FF)[1]!=n){
    stop(paste0('Error: FF should have one line for each latent variable in the model, got ',dim(FF)[1],'!=',n,'.'))
  }
  if(dim(D)[3]!=t){
    stop(paste0('Error: D should have 3º dimention equal to t or 1, got ',dim(D)[3],'.'))
  }
  if(dim(D)[1]!=n | dim(D)[2]!=n){
    stop(paste0('Error: D should have 1º and 2º dimentions equal the number of latent variables in the model, got (',dim(D)[1],',',dim(D)[2],')!=(',n,',',n,').'))
  }
  if(dim(W)[1]!=n | dim(W)[2]!=n){
    stop(paste0('Error: W should have 1º and 2º dimentions equal the number of latent variables in the model, got (',dim(W)[1],',',dim(W)[2],')!=(',n,',',n,').'))
  }
  if(dim(W)[3]!=t){
    stop(paste0('Error: W should have 3º dimention equal to t or 1, got ',dim(W)[3],'.'))
  }
  if(length(offset)!=t){
    stop(paste0('Error: Offset should have length 1 or equal to t, got ',dim(offset)[1],'!=',n,'.'))
  }
  #####

  G=model$G

  m0=model$mt[,t_last]
  C0=model$Ct[,,t_last]

  r=dim(model$FF)[2]

  D <- ifelse(D == 0, 1, D)

  last_m=m0
  last_C=C0

  # Definindo objetos
  at <- matrix(0, nrow=n, ncol=t)
  Rt <- array(0,dim=c(n,n,t))
  ft <- matrix(0, nrow=t, ncol=r)
  Qt <- array(0,dim=c(r,r,t))

  pred <- matrix(NA,dim(model$pred)[1],t)
  var.pred <- array(NA,c(dim(model$var.pred)[1],dim(model$var.pred)[2],t))
  icl.pred <- matrix(NA,dim(model$icl.pred)[1],t)
  icu.pred <- matrix(NA,dim(model$icu.pred)[1],t)

  for(t_i in c(1:t)){

    filter=model$kernel$filter(0,last_m,last_C,FF[,,t_i],G,D[,,t_i],W[,,t_i],offset[t_i])

    print(filter$Rt)

    last_m=filter$at
    last_C=filter$Rt
    Rt[,,t_i] <- filter$Rt
    at[,t_i] <- filter$at
    ft[t_i,] <- filter$ft
    Qt[,,t_i] <-  filter$Qt

    prediction=model$kernel$pred(filter,IC_prob)

    pred[,t_i]     <- prediction$pred
    var.pred[,,t_i] <- prediction$var.pred
    icl.pred[,t_i] <- prediction$icl.pred
    icu.pred[,t_i] <- prediction$icu.pred
  }

  return_list=list('pred'=pred,'var.pred'=var.pred,'icl.pred'=icl.pred,'icu.pred'=icu.pred,'at'=at,'Rt'=Rt)
  if(plot){
    r=dim(pred)[1]
    if(is.null(labels)){
      labels=c('Serie_' %>% paste0(1:r))
    }

    pred=cbind(t_last+c(1:t),t(pred) %>% as.data.frame)
    names(pred)=c('time',labels)
    pred=pred %>% pivot_longer(2:(r+1))
    names(pred)=c('time','Serie','Prediction')

    obs=pred
    obs$Prediction=NULL
    obs$Observation=NA

    icl.pred=cbind(t_last+c(1:t),t(icl.pred) %>% as.data.frame)
    names(icl.pred)=c('time',labels)
    icl.pred=icl.pred %>% pivot_longer(2:(r+1))
    names(icl.pred)=c('time','Serie','I.C.lower')

    icu.pred=cbind(t_last+c(1:t),t(icu.pred) %>% as.data.frame)
    names(icu.pred)=c('time',labels)
    icu.pred=icu.pred %>% pivot_longer(2:(r+1))
    names(icu.pred)=c('time','Serie','I.C.upper')

    plot_data=obs %>%
      inner_join(pred,c('time','Serie')) %>%
      inner_join(icl.pred,c('time','Serie')) %>%
      inner_join(icu.pred,c('time','Serie'))

    obs=cbind(c(1:t_last),model$data_out %>% as.data.frame)
    names(obs)=c('time',labels)
    obs=obs %>% pivot_longer(2:(r+1))
    names(obs)=c('time','Serie','Observation')

    pred=cbind(c(1:t_last),t(model$pred) %>% as.data.frame)
    names(pred)=c('time',labels)
    pred=pred %>% pivot_longer(2:(r+1))
    names(pred)=c('time','Serie','Prediction')

    icl.pred=cbind(c(1:t_last),t(model$icl.pred) %>% as.data.frame)
    names(icl.pred)=c('time',labels)
    icl.pred=icl.pred %>% pivot_longer(2:(r+1))
    names(icl.pred)=c('time','Serie','I.C.lower')

    icu.pred=cbind(c(1:t_last),t(model$icu.pred) %>% as.data.frame)
    names(icu.pred)=c('time',labels)
    icu.pred=icu.pred %>% pivot_longer(2:(r+1))
    names(icu.pred)=c('time','Serie','I.C.upper')

    plot_data2=obs %>%
      inner_join(pred,c('time','Serie')) %>%
      inner_join(icl.pred,c('time','Serie')) %>%
      inner_join(icu.pred,c('time','Serie'))

    max_value=calcula_max(plot_data2$Observation-min(plot_data2$Observation))[[3]]+min(plot_data2$Observation)
    min_value=-calcula_max(-(plot_data2$Observation-max(plot_data2$Observation)))[[3]]+max(plot_data2$Observation)

    plot_data=rbind(plot_data2,plot_data)

    plot=ggplotly(
        ggplot(plot_data)+
          geom_point(aes(x=time,y=Prediction,color=Serie,fill=Serie))+
          geom_ribbon(aes(x=time,ymin=I.C.lower,ymax=I.C.upper,fill=Serie,color=Serie),alpha=0.25)+
          geom_point(aes(x=time,y=Observation,color=Serie,fill=Serie))+
          scale_fill_hue('',na.value=NA)+
          scale_color_hue('',na.value=NA)+
          scale_x_continuous('Time')+
          scale_y_continuous('$y_t$')+
          theme_bw()+
          coord_cartesian(ylim=c(min_value,max_value))
      )
    return_list$plot=plot
  }

  return(return_list)
}
eval_past=function(model,smooth=FALSE,t_offset=0,IC_prob=0.95){
  if(smooth & t_offset>0){
    t_offset=0
    warning('t_offset is only used if smooth is set to TRUE.')
  }
  n=dim(model$mt)[1]
  t_last=dim(model$mt)[2]
  r=dim(model$FF)[2]
  k=dim(model$pred)[1]

  FF=model$FF
  G=diag(n)
  pred=matrix(NA,k,t_last)
  var.pred=array(NA,c(k,k,t_last))
  icl.pred=matrix(NA,k,t_last)
  icu.pred=matrix(NA,k,t_last)

  if(smooth){
    ref_mt=model$mts
    ref_Ct=model$Cts
    D=model$D**0
    W=model$W*0
  }else{
    ref_mt=model$mt
    ref_Ct=model$Ct
    D=model$D
    W=model$W
    for(i in t_offset){
      G=G%*%model$G
    }
  }

  at=array(0,c(n,t_last))
  Rt=array(0,c(n,n,t_last))
  for(i in c(1:t_last)){
    mt=if(i<=t_offset){model$m0}else{ref_mt[,i-t_offset]}
    Ct=if(i<=t_offset){model$C0}else{ref_Ct[,,i-t_offset]}

    filter=model$kernel$filter(model$data_out[i,],mt,Ct,FF[,,i]  %>% matrix(n,r),G,D[,,i],W[,,i],model$log_offset[i])
    prediction=model$kernel$pred(filter,IC_prob)

    pred[,i]      <- prediction$pred
    var.pred[,,i] <- prediction$var.pred
    icl.pred[,i]  <- prediction$icl.pred
    icu.pred[,i]  <- prediction$icu.pred
  }


  return(list('pred'=pred,'icl.pred'=icl.pred,'icu.pred'=icu.pred))
}
