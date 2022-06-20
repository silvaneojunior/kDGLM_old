mult.analise <- function(y,m01, C01,m02,C02, F1,F2,G1,G2,D1,D2, N){
  
  # Função a ser otimizada
  # 
  # otim1 <- function(x, f1,f2,q1,q2,q12, Omega){
  #   alpha1 <- x[1]
  #   alpha2 <- x[2]
  #   alpha3 <- x[3]
  #   
  #   eqs = rbind(
  #     f1 - digamma(alpha1) + digamma(alpha1 + alpha2 + alpha3),
  #     f2 - digamma(alpha2) + digamma(alpha1 + alpha2 + alpha3),
  #     q1 - trigamma(alpha1)  - trigamma(alpha1 + alpha2 + alpha3),
  #     q2 - trigamma(alpha2) - trigamma(alpha1 + alpha2 + alpha3),
  #     q12 - trigamma(alpha1 + alpha2 + alpha3))
  #   return(t(eqs)%*%Omega%*%eqs)
  # }
  # 
  
  model_tau0_e_tau1 <- function(x, parms) c(F1 =  f1 -  digamma(x[1]) + digamma(x[3] - x[2] - x[1]),
                                            F2 =  f2 -  digamma(x[2]) + digamma(x[3] - x[2] - x[1]),
                                            F3 =  media.log +  digamma(x[3]) - digamma(x[3] - x[2] - x[1]))
  
  
  
  f <- function(m) {
    m[upper.tri(m)] <- t(m)[upper.tri(m)]
    m
  }
  
  otim2 <- function(x, f1,f2,q1,q2,q12, Omega){
    tau1 <- x[1]
    tau2 <- x[2]
    tau0 <- x[3]
    media.log = log(1/(1+exp(f1) + exp(f2))) - 
      (q1/2)*(exp(f1)*exp(f2) + exp(f1))/(exp(f1)+exp(f2)+ 1)^2 - 
      (q2/2)*(exp(f2)*exp(f1) + exp(f2))/(exp(f1)+exp(f2)+1)^2 + 
      q12*(exp(f1 + f2)/(1+exp(f1)+exp(f2))) 
    
    eqs = rbind(
      f1  - digamma(tau1) + digamma(tau0 - tau1 - tau2),
      f2  - digamma(tau2) + digamma(tau0 - tau1 - tau2),
      digamma(tau0) - digamma(tau0 - tau1 - tau2) + media.log)
    
    return(t(eqs)%*%Omega%*%eqs)
  }
  
  
  
  # Definindo quantidades
  n1 <- nrow(F1)
  n2 <- nrow(F2)
  T <- nrow(y)
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
  
  alpha_1 <- rep(NA,l=T)
  alpha_2 <- rep(NA,l=T)
  alpha_3 <- rep(NA,l=T)
  
  alpha_1_star <- rep(NA,l=T)
  alpha_2_star <- rep(NA,l=T)
  alpha_3_star <- rep(NA,l=T)
  
  tau0_star = rep(NA, l = T)
  tau1_star = rep(NA, l = T)
  tau2_star = rep(NA, l = T)
  
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
  
  media.log = 
    log(1/(1 + exp(f1) + exp(f2))) - 
    (q1/2)*(exp(f1)*exp(f2) + exp(f1))/((exp(f1)+exp(f2)+ 1)^2) - 
    (q2/2)*(exp(f2)*exp(f1) + exp(f2))/((exp(f1)+exp(f2)+ 1)^2) + 
    q12*(exp(f1 + f2)/((exp(f1)+exp(f2)+ 1)^2))
  
  parms = c(f1,f2, media.log)
  
  ss1 <- multiroot(f = model_tau0_e_tau1 , start = c(0.01,0.01,0.03), parms = parms)
  
  
  
  tau1[1] <- ss1$root[1]
  tau2[1] <- ss1$root[2]
  tau0[1] <- ss1$root[3]
  
  
  alpha_1[1]      <- tau1[1] 
  alpha_2[1]      <- tau2[1]
  alpha_3[1]      <- tau0[1] - tau1[1] - tau2[1]
  
  alpha_1_star[1] <-    alpha_1[1]  +  y[1,1]
  alpha_2_star[1] <-    alpha_2[1]  +  y[1,2]
  alpha_3_star[1]  <-   alpha_3[1]  +  N[1] -  y[1,1] -  y[1,2]
  
  tau0_star[1]  <-  alpha_3_star[1] + alpha_2_star[1] + alpha_1_star[1]
  tau1_star[1]  <-  alpha_1_star[1] 
  tau2_star[1]  <-  alpha_2_star[1]
  
  # Posteriori
  f1star[1,] <-   digamma(alpha_1_star[1]) -  digamma(alpha_3_star[1])
  f2star[1,] <-   digamma(alpha_2_star[1]) -  digamma(alpha_3_star[1])
  Q1star[1,] <-   trigamma(alpha_1_star[1]) + trigamma(alpha_3_star[1])
  Q2star[1,] <-   trigamma(alpha_2_star[1]) + trigamma(alpha_3_star[1])
  Q12star[1,] <-  trigamma(alpha_3_star[1] )
  
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
    
    media.log = 
      log(1/(1 + exp(f1) + exp(f2))) - 
      (q1/2)*(exp(f1)*exp(f2) + exp(f1))/((exp(f1)+exp(f2)+ 1)^2) - 
      (q2/2)*(exp(f2)*exp(f1) + exp(f2))/((exp(f1)+exp(f2)+ 1)^2) + 
      q12*(exp(f1 + f2)/((exp(f1)+exp(f2)+ 1)^2))
    
    parms = c(f1,f2, media.log)
    
    
    ss1 <- multiroot(f = model_tau0_e_tau1 , start = c(0.01,0.01,0.04), parms = parms)
    
    tau1[t] <- ss1$root[1]
    tau2[t] <- ss1$root[2]
    tau0[t] <- ss1$root[3]
    
    
    alpha_1[t]      <- tau1[t]  
    alpha_2[t]      <- tau2[t] 
    alpha_3[t]      <- tau0[t] - tau1[t] - tau2[t] 
    
    alpha_1_star[t] <-    alpha_1[t]  +  y[t,1]
    alpha_2_star[t] <-    alpha_2[t]  +  y[t,2]
    alpha_3_star[t]  <-   alpha_3[t]  +  N[t] -  y[t,1] -  y[t,2]
    
    tau0_star[t]  <-  alpha_3_star[t] + alpha_2_star[t] + alpha_1_star[t]
    tau1_star[t]  <-  alpha_1_star[t] 
    tau2_star[t]  <-  alpha_2_star[t]
    
    # Posteriori
    f1star[t,] <-   digamma(alpha_1_star[t]) -  digamma(alpha_3_star[t])
    f2star[t,] <-   digamma(alpha_2_star[t]) -  digamma(alpha_3_star[t])
    Q1star[t,] <-   trigamma(alpha_1_star[t]) + trigamma(alpha_3_star[t])
    Q2star[t,] <-   trigamma(alpha_2_star[t]) + trigamma(alpha_3_star[t])
    Q12star[t,] <-  trigamma(alpha_3_star[t] )
    
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
                 mts, Cts, tau0_star, tau1_star,tau2_star , alpha_1_star,alpha_2_star,alpha_3_star, F, G, D, tau0, tau1, tau2)
  names(result) <- c("mt",  "Ct",  "f1star", "Q1star","ft", "qt",
                     "mts","Cts", "tau0_star", "tau1_star",  "tau_2_star",
                     "alpha_1_star", "alpha_2_star","alpha_3_star", "F", "G", "D", "tau0", "tau1", "tau2")
  return(result)
}
