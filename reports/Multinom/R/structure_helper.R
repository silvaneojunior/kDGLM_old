gera_bloco_poly <- function(order,value=1,name='Var_Poly',D=1,m0=0,C0=1,W=0){
  G=diag(order)
  t=ifelse(is.null(dim(value)),length(value),dim(value)[2])
  k=ifelse(is.null(dim(value)),1,dim(value)[1])
  #FF=matrix(c(ifelse(is.na(value),0,value),rep(0,(order-1)*t)),order,t,byrow = TRUE)
  FF=array(0,c(order,k,t))
  FF[1,,]=value

  if(order==2){
    G[1,2]=1
  }else{if(order>2){
    diag(G[1:(order-1),2:order])=1
  }}
  if(length(m0)<order){
    m0=rep(m0,order)
  }

  if(length(D)==1){
    D=array(1,c(order,order,t))*D
    D[,,apply(is.na(FF),2,any)]=1
  }
  if(length(W)==1){
    W=array(diag(order),c(order,order,t))*W
    W[,,apply(is.na(FF),2,any)]=0
  }
  FF=ifelse(is.na(FF),0,FF)

  names=list()
  names[[name]]=c(1:order)
  return(list('FF'=FF,
              'G'=G,
              'D'=D,
              'W'=W,
              'm0'=m0,
              'C0'=diag(order)*C0,
              'names'=names,
              'order'=order,
              'n'=order,
              't'=t,
              'k'=k))
}

gera_bloco_sazo <- function(period,value=1,name='Var_Sazo',D=1,m0=0,C0=1,W=0){
  w=2*pi/period
  order=2
  t=ifelse(is.null(dim(value)),length(value),dim(value)[2])
  k=ifelse(is.null(dim(value)),1,dim(value)[1])
  G <- matrix(c(cos(w),-sin(w),sin(w),cos(w)),order,order)
  if(length(m0)<2){
    m0=rep(m0,order)
  }
  #FF=matrix(c(ifelse(is.na(value),0,value),rep(0,(order-1)*t)),order,t,byrow = TRUE)
  FF=array(0,c(order,k,t))
  FF[1,,]=value
  if(length(D)==1){
    D=array(1,c(order,order,t))*D
    D[,,apply(is.na(FF),2,any)]=1
  }
  if(length(W)==1){
    W=array(diag(order),c(order,order,t))*W
    W[,,apply(is.na(FF),2,any)]=0
  }
  FF=ifelse(is.na(FF),0,FF)
  names=list()
  names[[name]]=c(1:order)
  return(list('FF'=FF,
              'G'=G,
              'D'=D,
              'W'=W,
              'm0'=m0,
              'C0'=diag(order)*C0,
              'names'=names,
              'period'=period,
              'n'=order,
              't'=t,
              'k'=k))
}

transf_block <- function(lag,value,name='Var_Poly_transf',D=1,m0=0,C0=1,W=0){
  G=diag(order)
  t=dim(value)[2]
  k=dim(value)[1]

  x=c(1:lag)
  mat=c()
  for(j in c(1:lag)){
    mat=c(mat,x**j)
  }

  M=matrix(mat,lag,lag,byrow=TRUE)
  pre_time=matrix(0,out_var,lag)
  extended_values=cbind(pre_time,value)
  pre_FF=matrix(0,(lag+1)*out_var,T_final)
  for(t in c(1:T_final)){
    for(out in c(1:out_var)){
      pre_FF[out_var*(k-1)+1:out_var,t]=M%*%t(extended_values[,(t):(t+k-1)])
    }
  }

  for(i in c(1:out_var)){
    placeholder=matrix(0,out_var,T_final)
    placeholder[i,]=vac_flag
    W=array(0,c(1,1,T_final))
    W[,,true_indice_inter]=1
    bloc_final=concat_bloco(bloc_final,
                            gera_bloco_poly(order=1,
                                            value=placeholder,
                                            name='vac_serie_' %>% paste0(i,'_',0),
                                            D=1/1,
                                            m0=0,
                                            C0=0,
                                            W=W))
                            }


  FF=array(0,c(order,k,t))
  FF[1,,]=value
  if(order==2){
    G[1,2]=1
  }else{if(order>2){
    diag(G[1:(order-1),2:order])=1
  }}
  if(length(m0)<order){
    m0=rep(m0,order)
  }

  if(length(D)==1){
    D=array(1,c(order,order,t))*D
    D[,,apply(is.na(value),2,any)]=1
  }
  if(length(W)==1){
    W=array(diag(order),c(order,order,t))*W
    W[,,apply(is.na(value),2,any)]=0
  }

  names=list()
  names[[name]]=c(1:order)
  return(list('FF'=FF,
              'G'=G,
              'D'=D,
              'W'=W,
              'm0'=m0,
              'C0'=diag(order)*C0,
              'names'=names,
              'order'=order,
              'n'=order,
              't'=t,
              'k'=k))
}

concat_bloco <- function(...){
  blocks=list(...)

  n=0
  t=1
  k=1
  names=list()
  for(block in blocks){
    ref_names=block$names
    for(name in names(ref_names)){
      ref_names[[name]]=ref_names[[name]]+n
    }
    names=c(names,ref_names)
    if(block$t>1){
      if(block$t!=t & t>1){
        stop(paste('Error: Blocks should have same length or length equal 1. Got',block$t,'and',t))
      }
      t=block$t
    }
    n=n+block$n
    k=max(block$k,k)
  }
  for(name in names(names)){
    ref_idx=which(names(names)==name)
    n_names=length(ref_idx)
    if(n_names>1){
      names(names)[ref_idx]=paste0(names(names)[ref_idx],'_',c(1:n_names))
    }
  }

  FF=array(0,c(n,k,t))
  G=matrix(0,n,n)
  D=array(0,c(n,n,t))
  W=array(0,c(n,n,t))
  m0=c()
  C0=matrix(0,n,n)
  position=1
  for(block in blocks){
    current_range=position:(position+block$n-1)
    FF[current_range,,]=block$FF
    G[current_range,current_range]=block$G
    D[current_range,current_range,]=block$D
    W[current_range,current_range,]=block$W
    m0=c(m0,block$m0)
    C0[current_range,current_range]=block$C0
    position=position+block$n
  }
  return(list('FF'=FF,'G'=G,'D'=D,'W'=W,'m0'=m0,'C0'=C0,'n'=n,'t'=t,'k'=k,'names'=names))
}
