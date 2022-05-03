gera_bloco_poly <- function(order,value=1,name='Var_Poly',D=1,m0=0,C0=1,W=0){
  G=diag(order)
  t=length(value)
  FF=matrix(c(ifelse(is.na(value),0,value),rep(0,(order-1)*t)),order,t,byrow = TRUE)
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
    D[,,is.na(value)]=1
  }
  if(length(W)==1){
    W=array(diag(order),c(order,order,t))*W
    W[,,is.na(value)]=0
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
              't'=t))
}
gera_bloco_sazo <- function(period,value=1,name='Var_Sazo',D=1,m0=0,C0=1,W=0){
  w=2*pi/period
  order=2
  t=length(value)
  G <- matrix(c(cos(w),-sin(w),sin(w),cos(w)),order,order)
  if(length(m0)<2){
    m0=rep(m0,order)
  }
  FF=c(ifelse(is.na(value),0,value),rep(0,t)) %>% matrix(order,t,byrow = TRUE)
  if(length(D)==1){
    D=array(1,c(order,order,t))*D
    D[,,is.na(value)]=1
  }
  if(length(W)==1){
    W=array(diag(order),c(order,order,t))*W
    W[,,is.na(value)]=0
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
              'period'=period,
              'n'=order,
              't'=t))
}
concat_bloco <- function(...){
  blocks=list(...)

  n=0
  t=1
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
  }
  for(name in names(names)){
    ref_idx=which(names(names)==name)
    n_names=length(ref_idx)
    if(n_names>1){
      names(names)[ref_idx]=paste0(names(names)[ref_idx],'_',c(1:n_names))
    }
  }

  FF=matrix(0,n,t)
  G=matrix(0,n,n)
  D=array(0,c(n,n,t))
  W=array(0,c(n,n,t))
  m0=c()
  C0=matrix(0,n,n)
  position=1
  for(block in blocks){
    current_range=position:(position+block$n-1)
    FF[current_range,]=block$FF
    G[current_range,current_range]=block$G
    D[current_range,current_range,]=block$D
    W[current_range,current_range,]=block$W
    m0=c(m0,block$m0)
    C0[current_range,current_range]=block$C0
    position=position+block$n
  }
  return(list('FF'=FF,'G'=G,'D'=D,'W'=W,'m0'=m0,'C0'=C0,'n'=n,'t'=t,'names'=names))
}
