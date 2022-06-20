normal.analise1 <- function(y,m01, C01,m02,C02, F1,F2,G1,G2,D1,D2){
  
  # Função a ser otimizada
  
  model_tau <- function(x, parms){
    output=c(
      F2 =  -digamma(exp(x[1])/2)+x[1]-q2/2-log(2)
    )
    return(output)
  }
  
  jacfunc=function(x,parms){
    jac=matrix(c(-0.5*exp(x[1])*trigamma(exp(x[1])/2)+1),
               1,1)
    return(jac)
  }
  
  
  # Definindo quantidades
  n1 <- nrow(F1)
  n2 <- nrow(F2)
  T <- length(y)
  n <- n1 + n2
  F <- as.matrix(bdiag(F1[,1], F2[,1]))
  G <- as.matrix(bdiag(G1,G2))
  
  D.aux <- as.matrix(bdiag(D1, D2))
  D <- ifelse(D.aux == 0, 1, D.aux)
  
  m0 <- as.matrix(c(m01, m02), nrow = 1)
  C0 <- as.matrix(bdiag(C01, C02))
  
  r = 2
  mt <- matrix(0,nrow=n,ncol=T)
  Ct <- array(rep(diag(n),T),dim=c(n,n,T))
  Rt <- array(rep(0,T),dim=c(n,n,T))
  Pt <- array(rep(0,T),dim=c(n,n,T))
  Wt <- array(rep(0,T),dim=c(n,n,T))
  ft <- matrix(0,nrow=T,ncol=r)
  at <- matrix(0,nrow=n,ncol=T)
  Qt <-  array(rep(0,T),dim=c(r,r,T))
  et <- matrix(0,nrow=T,ncol=r)
  At <- array(rep(0,T),dim=c(n,r,T))
  dt <- matrix(0,nrow=T,ncol=r)
  nt <- matrix(0,nrow=T,ncol=r)
  
  f1star <- matrix(0,nrow=T)
  f2star <- matrix(0,nrow=T)
  Q1star <- matrix(0,nrow=T)
  Q2star <- matrix(0,nrow=T)
  Q12star <- matrix(0,nrow=T)
  fstar <-  matrix(0,nrow=T,ncol=r)
  Qstar <-  matrix(0,nrow=T,ncol=r)
  
  m1 <- matrix(0,nrow=n1,ncol=T)
  C1 <- array(rep(diag(n1),T),dim=c(n1,n1,T))
  
  tau0 <- rep(NA,l=T)
  tau1 <- rep(NA,l=T)
  tau2 <- rep(NA,l=T)
  tau3 <- rep(NA,l=T)
  
  tau0_star = rep(NA, l = T)
  tau1_star = rep(NA, l = T)
  tau2_star = rep(NA, l = T)
  tau3_star = rep(NA, l = T)
  
  # #Auxiliar do desconto
  if(is.null(dim(G)) == FALSE){
    matrixaux1 <- G == 0
    if(G[1,1] == 1 & G[1,2] == 1 & G[2,2] == 1 ) {matrixaux1[2] <- 0}
    tira1 <- which(matrixaux1 == 1)
    mantem1 <- which(matrixaux1 == 0)}
  # 
  # if(is.null(dim(G2)) == FALSE){
  #   matrixaux2 <- G2 == 0 
  #   if(G2[1,1] == 1 & G2[1,2] == 1 & G2[2,2] == 1 ) {matrixaux2[2] <- 0}
  #   tira2 <- which(matrixaux2 == 1)
  #   mantem2 <- which(matrixaux2 == 0)}
  # 
  
  # Priori em t = 1
  
  D=ifelse(D==0,1,D)
  
  at[,1]          = G%*%m0
  Pt            <-    G%*%C0%*%(t(G))
  Rt[,,1]       <- D*Pt
  
  
  # Previsao em t = 1
  ft[1,]        <- t(F)%*%at[,1]
  Qt[,,1]       <- t(F)%*%Rt[,,1]%*%F 
  
  # minimizando...
  
  f1=ft[1,1]
  f2=ft[1,2]
  q1 = Qt[1,1,1]
  q2 = Qt[2,2,1]
  q12 = Qt[1,2,1]
  
  c0=1/(exp(f2 + q2/2)*q1)
  mu0=f1
  
  # parms = c(f1,f2, q1, q2,c0,m0)
  # 
  # ss1 <- multiroot(f = model_tau , start = c(0), parms = parms,jacfunc = jacfunc, jactype='fullusr')
  # 
  # 
  # n0=exp(ss1$root[1])
  
  n0=2/q2
  d0=n0/exp(f2+q2/2)
  
  tau0[1]=n0/2-1/2
  tau1[1]=-c0/2
  tau2[1]=c0*mu0
  tau3[1]=-(c0*(mu0**2)+d0)/2
  
  tau0_star[1]  <-  tau0[1] + 1/2
  tau1_star[1]  <-  tau1[1] - 1/2
  tau2_star[1]  <-  tau2[1] + y[1]
  tau3_star[1]  <-  tau3[1] - y[1]^2
  
  trigamma(n0/2+1/2)+(digamma(n0/2+1/2)-log(d0/2))
  
  aux_1=digamma(tau0_star[1] + 0.5) -log(((tau2_star[1]^2)/(4*tau1_star[1])) - tau3_star[1])
  
  # Posteriori
  f1star[1,] <-   -tau2_star[1]/(2*tau1_star[1])
  f2star[1,] <-   aux_1
  Q1star[1,] <-   (tau2_star[1]^2)/(4*tau1_star[1]^2) - 8*(tau1_star[1]^2)*(tau0_star[1] + 1/2)/(tau2_star[1]^2 - 4*tau1_star[1]*(tau3_star[1]))
  Q2star[1,] <-   trigamma(tau0_star[1] + 0.5) + aux_1^2
  Q12star[1,] <-   (-tau2_star[1]/(2*tau1_star[1]))*aux_1
  
  fstar <- c(f1star[1,],   f2star[1,])
  Qstar <- matrix(c( Q1star[1,],  Q12star[1,],  Q12star[1,],  Q2star[1,]), byrow =F, ncol = 2)
  
  At[,,1] <- Rt[,,1]%*%F%*%ginv(Qt[,,1])
  mt[,1] <- at[,1] + At[,,1]%*%(fstar -ft[1,])
  Ct[,,1] <- Rt[,,1] +  At[,,1]%*%(Qstar - Qt[,,1])%*%t(At[,,1])
  
  for(t in 2:T){
    at[,t] = G%*%mt[,t-1] 
    Pt <- G%*%Ct[,,t-1]%*%(t(G))
    Rt[,,t] <- D*Pt
    
    # Previsao em t = 1
    ft[t,] <-  t(F)%*%at[,t]
    Qt[,,t] <- t(F)%*%Rt[,,t]%*%F
    
    # minimizando...
    
    f1=ft[t,1]
    f2=ft[t,2]
    q1 = Qt[1,1,t]
    q2 = Qt[2,2,t]
    q12 = Qt[1,2,t]
    
    c0=1/(exp(f2 + q2/2)*q1)
    mu0=f1
    
    # parms = c(f1,f2, q1, q2,c0,mu0)
    # 
    # ss1 <- multiroot(f = model_tau , start = c(0), parms = parms,jacfunc = jacfunc, jactype='fullusr')
    # 
    # n0=exp(ss1$root[1])
    
    
    n0=2/q2
    d0=n0/exp(f2+q2/2)
    
    tau0[t]=n0/2-1/2
    tau1[t]=-c0/2
    tau2[t]=c0*mu0
    tau3[t]=-(c0*(mu0**2)+d0)/2
    
    tau0_star[t]  <-  tau0[t] + 1/2
    tau1_star[t]  <-  tau1[t] - 1/2
    tau2_star[t]  <-  tau2[t] + y[t]
    tau3_star[t]  <-  tau3[t] - y[t]^2
    
    aux_1=digamma(tau0_star[t] + 0.5) -log(((tau2_star[t]^2)/(4*tau1_star[t])) - tau3_star[t])
    
    # Posteriori
    f1star[t,] <-   -tau2_star[t]/(2*tau1_star[t])
    f2star[t,] <-   aux_1
    Q1star[t,] <-   (tau2_star[t]^2)/(4*tau1_star[t]^2) - 8*(tau1_star[t]^2)*(tau0_star[t] + 1/2)/(tau2_star[t]^2 - 4*tau1_star[t]*(tau3_star[t]))
    Q2star[t,] <-   trigamma(tau0_star[t] + 0.5) + aux_1^2
    Q12star[t,] <-   (-tau2_star[t]/(2*tau1_star[t]))*aux_1
    
    
    fstar <- c(f1star[t,],   f2star[t,])
    Qstar <- matrix(c( Q1star[t,],  Q12star[t,],  Q12star[t,],  Q2star[t,]), byrow =F, ncol = 2)
    
    At[,,t] <- Rt[,,t]%*%F%*%ginv(Qt[,,t])
    mt[,t] <- at[,t] + At[,,t]%*%(fstar -ft[t,])
    Ct[,,t] <- Rt[,,t] +  At[,,t]%*%(Qstar - Qt[,,t])%*%t(At[,,t])
    
  }
  
  
  mts <- matrix(0, ncol=T, nrow=n)
  Cts <- array(rep(diag(n),T),dim=c(n,n,T))
  
  mts[,T] <- mt[, T]
  Cts[,,T] <- Ct[,,T]
  for(t in (T-1):1){
    mts[,t] <- mt[,t] + Ct[,,t]%*%t(G)%*%solve(Rt[,,t+1])%*%(mts[,t+1] - at[,t+1])
    Cts[,,t] <- Ct[,,t] - Ct[,,t]%*%t(G)%*%solve(Rt[,,t+1])%*%(Rt[,,t+1] - Cts[,,t+1])%*%t(Ct[,,t]%*%t(G)%*%solve(Rt[,,t+1]))
  }
  
  result <- list(mt,Ct,at,At,Rt,f1star, Q1star,ft, Qt,
                 mts, Cts, tau0_star, tau1_star,tau2_star ,tau3_star, F, G, D, tau0, tau1, tau2, tau3)
  names(result) <- c("mt",  "Ct", "at","At","Rt",  "f1star", "Q1star","ft", "qt",
                     "mts","Cts", "tau0_star", "tau1_star",  "tau_2_star", "tau3_star",
                     "F", "G", "D", "tau0", "tau1", "tau2", "tau4")
  return(result)
}


normal.analise2 <- function(y,m01, C01,m02,C02, F1,F2,G1,G2,D1,D2){
  
  # Função a ser otimizada
  
  model_tau <- function(x, parms){
    output=c(
      F2 =  -digamma(exp(x[1])/2)+x[1]-q2/2-log(2)
    )
    return(output)
  }
  
  jacfunc=function(x,parms){
    jac=matrix(c(-0.5*exp(x[1])*trigamma(exp(x[1])/2)+1),
               1,1)
    return(jac)
  }
  
  
  # Definindo quantidades
  n1 <- nrow(F1)
  n2 <- nrow(F2)
  T <- length(y)
  n <- n1 + n2
  F <- as.matrix(bdiag(F1[,1], F2[,1]))
  G <- as.matrix(bdiag(G1,G2))
  
  D.aux <- as.matrix(bdiag(D1, D2))
  D <- ifelse(D.aux == 0, 1, D.aux)
  
  m0 <- as.matrix(c(m01, m02), nrow = 1)
  C0 <- as.matrix(bdiag(C01, C02))
  
  r = 2
  mt <- matrix(0,nrow=n,ncol=T)
  Ct <- array(rep(diag(n),T),dim=c(n,n,T))
  Rt <- array(rep(0,T),dim=c(n,n,T))
  Pt <- array(rep(0,T),dim=c(n,n,T))
  Wt <- array(rep(0,T),dim=c(n,n,T))
  ft <- matrix(0,nrow=T,ncol=r)
  at <- matrix(0,nrow=n,ncol=T)
  Qt <-  array(rep(0,T),dim=c(r,r,T))
  et <- matrix(0,nrow=T,ncol=r)
  At <- array(rep(0,T),dim=c(n,r,T))
  dt <- matrix(0,nrow=T,ncol=r)
  nt <- matrix(0,nrow=T,ncol=r)
  
  f1star <- matrix(0,nrow=T)
  f2star <- matrix(0,nrow=T)
  Q1star <- matrix(0,nrow=T)
  Q2star <- matrix(0,nrow=T)
  Q12star <- matrix(0,nrow=T)
  fstar <-  matrix(0,nrow=T,ncol=r)
  Qstar <-  matrix(0,nrow=T,ncol=r)
  
  m1 <- matrix(0,nrow=n1,ncol=T)
  C1 <- array(rep(diag(n1),T),dim=c(n1,n1,T))
  
  tau0 <- rep(NA,l=T)
  tau1 <- rep(NA,l=T)
  tau2 <- rep(NA,l=T)
  tau3 <- rep(NA,l=T)
  
  tau0_star = rep(NA, l = T)
  tau1_star = rep(NA, l = T)
  tau2_star = rep(NA, l = T)
  tau3_star = rep(NA, l = T)
  
  # #Auxiliar do desconto
  if(is.null(dim(G)) == FALSE){
    matrixaux1 <- G == 0
    if(G[1,1] == 1 & G[1,2] == 1 & G[2,2] == 1 ) {matrixaux1[2] <- 0}
    tira1 <- which(matrixaux1 == 1)
    mantem1 <- which(matrixaux1 == 0)}
  # 
  # if(is.null(dim(G2)) == FALSE){
  #   matrixaux2 <- G2 == 0 
  #   if(G2[1,1] == 1 & G2[1,2] == 1 & G2[2,2] == 1 ) {matrixaux2[2] <- 0}
  #   tira2 <- which(matrixaux2 == 1)
  #   mantem2 <- which(matrixaux2 == 0)}
  # 
  
  # Priori em t = 1
  
  D=ifelse(D==0,1,D)
  
  at[,1]          = G%*%m0
  Pt            <-    G%*%C0%*%(t(G))
  Rt[,,1]       <- D*Pt
  
  
  # Previsao em t = 1
  ft[1,]        <- t(F)%*%at[,1]
  Qt[,,1]       <- t(F)%*%Rt[,,1]%*%F 
  
  # minimizando...
  
  f1=ft[1,1]
  f2=ft[1,2]
  q1 = Qt[1,1,1]
  q2 = Qt[2,2,1]
  q12 = Qt[1,2,1]
  
  c0=1/(exp(f2 + q2/2)*q1)
  mu0=f1
  
  parms = c(f1,f2, q1, q2,c0,m0)

  ss1 <- multiroot(f = model_tau , start = c(0), parms = parms,jacfunc = jacfunc, jactype='fullusr')


  n0=exp(ss1$root[1])
  d0=n0/exp(f2+q2/2)
  
  tau0[1]=n0/2-1/2
  tau1[1]=-c0/2
  tau2[1]=c0*mu0
  tau3[1]=-(c0*(mu0**2)+d0)/2
  
  tau0_star[1]  <-  tau0[1] + 1/2
  tau1_star[1]  <-  tau1[1] - 1/2
  tau2_star[1]  <-  tau2[1] + y[1]
  tau3_star[1]  <-  tau3[1] - y[1]^2
  
  trigamma(n0/2+1/2)+(digamma(n0/2+1/2)-log(d0/2))
  
  aux_1=digamma(tau0_star[1] + 0.5) -log(((tau2_star[1]^2)/(4*tau1_star[1])) - tau3_star[1])
  
  # Posteriori
  f1star[1,] <-   -tau2_star[1]/(2*tau1_star[1])
  f2star[1,] <-   aux_1
  Q1star[1,] <-   (tau2_star[1]^2)/(4*tau1_star[1]^2) - 8*(tau1_star[1]^2)*(tau0_star[1] + 1/2)/(tau2_star[1]^2 - 4*tau1_star[1]*(tau3_star[1]))
  Q2star[1,] <-   trigamma(tau0_star[1] + 0.5) + aux_1^2
  Q12star[1,] <-   (-tau2_star[1]/(2*tau1_star[1]))*aux_1
  
  fstar <- c(f1star[1,],   f2star[1,])
  Qstar <- matrix(c( Q1star[1,],  Q12star[1,],  Q12star[1,],  Q2star[1,]), byrow =F, ncol = 2)
  
  At[,,1] <- Rt[,,1]%*%F%*%ginv(Qt[,,1])
  mt[,1] <- at[,1] + At[,,1]%*%(fstar -ft[1,])
  Ct[,,1] <- Rt[,,1] +  At[,,1]%*%(Qstar - Qt[,,1])%*%t(At[,,1])
  
  for(t in 2:T){
    at[,t] = G%*%mt[,t-1] 
    Pt <- G%*%Ct[,,t-1]%*%(t(G))
    Rt[,,t] <- D*Pt
    
    # Previsao em t = 1
    ft[t,] <-  t(F)%*%at[,t]
    Qt[,,t] <- t(F)%*%Rt[,,t]%*%F
    
    # minimizando...
    
    f1=ft[t,1]
    f2=ft[t,2]
    q1 = Qt[1,1,t]
    q2 = Qt[2,2,t]
    q12 = Qt[1,2,t]
    
    c0=1/(exp(f2 + q2/2)*q1)
    mu0=f1
    
    parms = c(f1,f2, q1, q2,c0,mu0)

    ss1 <- multiroot(f = model_tau , start = c(0), parms = parms,jacfunc = jacfunc, jactype='fullusr')

    n0=exp(ss1$root[1])
    d0=n0/exp(f2+q2/2)
    
    tau0[t]=n0/2-1/2
    tau1[t]=-c0/2
    tau2[t]=c0*mu0
    tau3[t]=-(c0*(mu0**2)+d0)/2
    
    tau0_star[t]  <-  tau0[t] + 1/2
    tau1_star[t]  <-  tau1[t] - 1/2
    tau2_star[t]  <-  tau2[t] + y[t]
    tau3_star[t]  <-  tau3[t] - y[t]^2
    
    aux_1=digamma(tau0_star[t] + 0.5) -log(((tau2_star[t]^2)/(4*tau1_star[t])) - tau3_star[t])
    
    # Posteriori
    f1star[t,] <-   -tau2_star[t]/(2*tau1_star[t])
    f2star[t,] <-   aux_1
    Q1star[t,] <-   (tau2_star[t]^2)/(4*tau1_star[t]^2) - 8*(tau1_star[t]^2)*(tau0_star[t] + 1/2)/(tau2_star[t]^2 - 4*tau1_star[t]*(tau3_star[t]))
    Q2star[t,] <-   trigamma(tau0_star[t] + 0.5) + aux_1^2
    Q12star[t,] <-   (-tau2_star[t]/(2*tau1_star[t]))*aux_1
    
    
    fstar <- c(f1star[t,],   f2star[t,])
    Qstar <- matrix(c( Q1star[t,],  Q12star[t,],  Q12star[t,],  Q2star[t,]), byrow =F, ncol = 2)
    
    At[,,t] <- Rt[,,t]%*%F%*%ginv(Qt[,,t])
    mt[,t] <- at[,t] + At[,,t]%*%(fstar -ft[t,])
    Ct[,,t] <- Rt[,,t] +  At[,,t]%*%(Qstar - Qt[,,t])%*%t(At[,,t])
    
  }
  
  
  mts <- matrix(0, ncol=T, nrow=n)
  Cts <- array(rep(diag(n),T),dim=c(n,n,T))
  
  mts[,T] <- mt[, T]
  Cts[,,T] <- Ct[,,T]
  for(t in (T-1):1){
    mts[,t] <- mt[,t] + Ct[,,t]%*%t(G)%*%solve(Rt[,,t+1])%*%(mts[,t+1] - at[,t+1])
    Cts[,,t] <- Ct[,,t] - Ct[,,t]%*%t(G)%*%solve(Rt[,,t+1])%*%(Rt[,,t+1] - Cts[,,t+1])%*%t(Ct[,,t]%*%t(G)%*%solve(Rt[,,t+1]))
  }
  
  result <- list(mt,Ct,f1star, Q1star,ft, Qt,
                 mts, Cts, tau0_star, tau1_star,tau2_star ,tau3_star, F, G, D, tau0, tau1, tau2, tau3)
  names(result) <- c("mt",  "Ct",  "f1star", "Q1star","ft", "qt",
                     "mts","Cts", "tau0_star", "tau1_star",  "tau_2_star", "tau3_star",
                     "F", "G", "D", "tau0", "tau1", "tau2", "tau4")
  return(result)
}
