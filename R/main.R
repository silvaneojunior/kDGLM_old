# library(Matrix)   # Pode ser substituido pela bilbioteca Matrix.
# library(tidyr)
# library(dplyr)

source('R/plot_helper.R')
source('R/structure_helper.R')
source('R/kernels.R')

ajusta_modelo <- function(...,data_out,kernel=poisson_gi_exp,offset=NULL,log_offset=NULL){
  structure=concat_bloco(...)
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
    log_offset=rep(log_offset,length(data_out))
    offset=rep(offset,length(data_out))
  }
  if(length(data_out)!=length(log_offset)){
    stop('Erro: offset/log_offset does not have the same length as data_out')
  }

  if(structure$t==1){
    structure$D=array(structure$D,c(structure$n,structure$n,length(data_out)))
    structure$W=array(structure$W,c(structure$n,structure$n,length(data_out)))
    structure$FF=matrix(structure$FF,structure$n,length(data_out))
    structure$t=length(data_out)
  }
  if(length(data_out)!=structure$t){
    stop('Erro: data_out does not have the same length as structure.')
  }

  model=kernel(y=data_out,
               m0=structure$m0,
               C0=structure$C0,
               F1=structure$FF,
               G1=structure$G,
               D1=structure$D,
               W1=structure$W,
               pop=log_offset)

  model$names=structure$names
  model$data_out=data_out

  return(model)

}

predict=function(model,t=1,offset=NULL,log_offset=NULL,FF=NULL,D=NULL,W=NULL,plot=TRUE,IC_prob=0.95){
  n=dim(model$mt)[1]
  t_last=dim(model$mt)[2]

  if(t>10){
    warning('Warning: Prediction window is big, results will probabily be unreliable.')
  }

  #### Consistency check ####
  if(is.null(FF)){
    FF=array(model$F[,t_last],c(n,t))
  }
  if(is.null(D)){
    D=array(model$D[,,t_last],c(n,n,t))
  }else{if(all(dim(D)==1)){
    D=array(D,c(n,n,t))
  }else{if(length(dim(D))==2 | (length(dim(D))==3 & dim(D)[3]==1)){
    D=array(D,c(n,n,t))
  }}
  }
  if(is.null(W)){
    W=array(model$W[,,t_last],c(n,n,t))
  }else{if(all(dim(W)==1)){
    W=array(diag(n)*W,c(n,n,t))
  }else{if(length(dim(W))==2 | (length(dim(W))==3 & dim(W)[3]==1)){
    W=array(W,c(n,n,t))
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

  if(dim(FF)[2]!=t){
    stop(paste0('Error: FF should have one column for each time or exactly 1 column, got ',dim(FF)[2],'!=',t,'.'))
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

  D <- ifelse(D == 0, 1, D)

  # Definindo objetos
  at <- matrix(0, ncol=t, nrow=n)
  mt <- matrix(0, ncol=t, nrow=n)
  ft <- matrix(0, ncol=1, nrow=t)
  qt <- matrix(0, ncol=1, nrow=t)
  Ct <- array(rep(diag(n),t),dim=c(n,n,t))
  Rt <- array(rep(diag(n),t),dim=c(n,n,t))
  a = b= 0
  pred = var.pred = icl.pred = icu.pred = matrix(0, ncol=1, nrow=t)

  ## Algoritmo

  # Priori

  at[,1] <- G%*%m0
  Rt[,,1] <-G%*%C0%*%(t(G))*D[,,1]+W[,,1]

  reduc_RFF=Rt[,,1]%*%FF[,1]

  # Previsão 1 passo a frente
  ft[1,] <- t(FF[,1])%*%at[,1] + log_offset[1]
  qt[1,] <- t(FF[,1])%*%reduc_RFF

  a[1] <- (1/qt[1,])
  b[1] <- (exp(-ft[1,] -0.5*qt[1,])/(qt[1,]))

  # Preditiva em t = 1

  pred[1] <- a[1]/ b[1]
  var.pred <- a[1]*(b[1]+1)/(b[1])^2
  icl.pred[1]<-qnbinom((1-IC_prob)/2, a[1], (b[1]/(b[1] +1)))
  icu.pred[1]<-qnbinom(1-(1-IC_prob)/2, a[1], (b[1]/(b[1] +1)))
  if(t>1){
  for(i in c(2:t)){
    # Priori

    at[,i] <- G%*%at[,i-1]
    Rt[,,i] <-G%*%Rt[,,i-1]%*%(t(G))*D[,,i]+W[,,i]

    reduc_RFF=Rt[,,i]%*%FF[,i]

    # Previsão 1 passo a frente
    ft[i,] <- t(FF[,i])%*%at[,i] + log_offset[i]
    qt[i,] <- t(FF[,i])%*%reduc_RFF

    a[i] <- (1/qt[i,])
    b[i] <- (exp(-ft[i,] -0.5*qt[i,])/(qt[i,]))

    pred[i] <- a[i]/ b[i]
    var.pred <- a[i]*(b[i]+1)/(b[i])^2
    icl.pred[i]<-qnbinom((1-IC_prob)/2, a[i], (b[i]/(b[i] +1)))
    icu.pred[i]<-qnbinom(1-(1-IC_prob)/2, a[i], (b[i]/(b[i] +1)))
  }}
  return_list=list('pred'=pred,'var.pred'=var.pred,'icl.pred'=icl.pred,'icu.pred'=icu.pred,'at'=at,'Rt'=Rt)
  if(plot){
    fill_list=c('#2596be','#2596be','black')
    names(fill_list)=c(paste0('Prediction I.C. (',(100*IC_prob) %>% round(),'%)'),'Prediction','Observed values')
    color_list=c('#2596be','black')
    names(color_list)=c('Prediction','Observed values')

    plot=ggplotly(
        ggplot()+
          geom_point(aes(x=c(1:t)+t_last,y=pred,color='Prediction',fill='Prediction'))+
          geom_ribbon(aes(x=c(1:t)+t_last,ymin=icl.pred,ymax=icu.pred,fill=paste0('Prediction I.C. (',(100*IC_prob) %>% round(),'%)'),color=paste0('Prediction I.C. (',(100*IC_prob) %>% round(),'%)')),alpha=0.25)+
          geom_point(aes(x=c(1:t_last),y=model$data_out,color='Observed values',fill='Observed values'))+
          scale_fill_manual('',na.value=NA,values=fill_list)+
          scale_color_manual('',na.value=NA,values=color_list)+
          scale_x_continuous('Time')+
          scale_y_continuous('$y_t$')+
          theme_bw()
      )
    return_list$plot=plot
  }

  return(return_list)
}
eval_past=function(model,smooth=FALSE,t_offset=0){
  if(smooth & t_offset>0){
    t_offset=0
    warning('t_offset is only used if smooth is set to TRUE.')
  }
  n=dim(model$mt)[1]
  t_last=dim(model$mt)[2]

  at=array(0,c(n,t_last))
  Rt=array(0,c(n,n,t_last))

  ft=array(0,c(t_last))
  qt=array(0,c(t_last))
  a=array(0,c(t_last))
  b=array(0,c(t_last))
  FF=model$F
  G=array(model$G,c(n,n,t_last))
  D=model$D
  W=model$W

  pred=c(1:t_last)*0

  for(i in c(1:(t_last-t_offset))+t_offset){
    if(smooth){
      at[,i]=model$mts[,i]
      Rt[,,i]=model$Cts[,,i]
    }else{
      at[,i]=model$mt[,i-t_offset]
      Rt[,,i]=model$Ct[,,i-t_offset]
      if(t_offset>0){
        at[,i]=G[,,i-t_offset+1]%*%at[,i]
        Rt[,,i]=G[,,i-t_offset+1]%*%Rt[,,i]%*%t(G[,,i-t_offset+1])*D[,,1]+W[,,1]
        if(t_offset>1){
          multi_G=diag(n)
          for(j in c(2:t_offset)){
            multi_G=G[,,i-t_offset+j]%*%multi_G
          }
          at[,i]=multi_G%*%at[,i]
          Rt[,,i]=multi_G%*%Rt[,,i]%*%t(multi_G)
        }
      }
    }

    ft[i] <- t(FF[,i])%*%at[,i] + model$log_offset[i]
    qt[i] <- t(FF[,i])%*%Rt[,,i]%*%FF[,i]

    a[i] <- (1/qt[i])
    b[i] <- (exp(-ft[i] -0.5*qt[i])/(qt[i]))

    pred[i] <- a[i]/ b[i]
    #icl.pred[i]<-qnbinom((1-IC_prob)/2, a[i], (b[i]/(b[i] +1)))
    #icu.pred[i]<-qnbinom(1-(1-IC_prob)/2, a[i], (b[i]/(b[i] +1)))
  }

  return(list('pred'=pred,'a'=a,'b'=b))
}
