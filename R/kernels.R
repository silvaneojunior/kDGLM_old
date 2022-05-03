poisson_gi_exp <- function(y,m0 = 0, C0 = 1, F1,G1,D1,W1, pop, IC_prob=0.95){

  # Definindo quantidades
  n1 <- nrow(F1)

  r <- 1
  T <- length(y)
  n <- n1
  FF <- F1
  G <- G1

  # D1.aux <- matrix(rep(D1,n1^2), ncol = n1 )
  # D2.aux <- matrix(rep(D2,n2^2), ncol = n2 )
  # D.aux <- as.matrix(bdiag(D1.aux, D2.aux))
  D.aux <- D1
  D <- ifelse(D.aux == 0, 1, D.aux)

  W <- W1

  # Definindo objetos
  at <- matrix(0, ncol=T, nrow=n)
  mt <- matrix(0, ncol=T, nrow=n)
  ft <- matrix(0, ncol=1, nrow=T)
  qt <- matrix(0, ncol=1, nrow=T)
  alphat <- matrix(0, ncol=1, nrow=T)
  betat <- matrix(0, ncol=1, nrow=T)
  Ct <- array(rep(diag(n),T),dim=c(n,n,T))
  Rt <- array(rep(diag(n),T),dim=c(n,n,T))
  gt = pt = a = b= 0
  rep <- 1
  pred = var.pred = icl.pred = icu.pred = matrix(0, ncol=rep, nrow=T)
  media.post = var.post = icu.post = icl.post = matrix(0, ncol=rep, nrow=T)
  a.post = b.post = eqm =0
  E.th3 = E.th4 = raiz2 = raiz1= matrix(0, ncol=rep, nrow=T)
  tau0_star <- rep(NA, l = T)
  tau1_star <- rep(NA, l = T)
  tau0 <- rep(NA, l = T)
  tau1 <- rep(NA, l = T)
  fstar <- rep(NA, l = T)
  qstar <- rep(NA, l = T)

  norm_ic=qnorm(1-(1-IC_prob)/2)

  ## Algoritmo

  # Priori

  at[,1] <- G%*%m0
  Rt[,,1] <-G%*%C0%*%(t(G))*D[,,1]+W[,,1]

  reduc_RFF=Rt[,,1]%*%FF[,1]

  # Previsão 1 passo a frente
  ft[1,] <- t(FF[,1])%*%at[,1] + pop[1]
  qt[1,] <- t(FF[,1])%*%reduc_RFF

  # Compatibilizando prioris

  a[1] <- (1/qt[1,])
  b[1] <- (exp(-ft[1,] -0.5*qt[1,])/(qt[1,]))

  # Posteriori

  a.post[1] <- a[1] + y[1]
  b.post[1] <- b[1] + 1

  gt[1] <- log(a.post[1]/b.post[1]) + 1/(2*a.post[1])
  pt[1] <- (2*a.post[1]-1)/(2*a.post[1]^2)

  mt[,1] <- at[,1]+reduc_RFF*(gt[1]-ft[1,])*(1/(qt[1,]))
  Ct[,,1] <- Rt[,,1] - (reduc_RFF%*%t(reduc_RFF))*(1 - pt[1]/qt[1,])*(1/qt[1,])

  media.post[1] <- a.post[1]/(b.post[1])
  var.post[1] <- a.post[1]/(b.post[1]^2)
  icu.post[1] <- media.post[1] + norm_ic*sqrt(var.post[1])
  icl.post[1] <- media.post[1] - norm_ic*sqrt(var.post[1])

  # Preditiva em t = 1

  pred[1] <- a[1]/ b[1]
  var.pred <- a[1]*(b[1]+1)/(b[1])^2
  icl.pred[1]<-qnbinom((1-IC_prob)/2, a[1], (b[1]/(b[1] +1)))
  icu.pred[1]<-qnbinom(1-(1-IC_prob)/2, a[1], (b[1]/(b[1] +1)))

  # Passo 2 até t

  start<- proc.time()
  for(t in 2:T){

    # Priori

    at[,t] <-  G%*%mt[,t-1]
    Rt[,,t] <- G%*%Ct[,,t-1]%*%(t(G))*D[,,t]+W[,,t]

    reduc_RFF=Rt[,,t]%*%FF[,t]

    # Previsão 1 passo a frente

    ft[t,] <- t(FF[,t])%*%at[,t] + pop[t]
    qt[t,] <- t(FF[,t])%*%reduc_RFF

    # Compatibilizando prioris

    a[t] <- (1/qt[t,])
    b[t] <- (exp(-ft[t,] )/(qt[t,]))

    # Posteriori

    a.post[t] <- a[t] + y[t]
    b.post[t] <- b[t] + 1

    gt[t] <- log(a.post[t]/b.post[t])+1/(2*a.post[t])
    pt[t] <- (2*a.post[t]-1)/(2*a.post[t]^2)

    mt[,t] <- at[,t]+reduc_RFF*(gt[t]-ft[t,])*(1/(qt[t,]))
    Ct[,,t] <- Rt[,,t] - (reduc_RFF%*%t(reduc_RFF))*(1 - pt[t]/qt[t,])*(1/qt[t,])

    media.post[t] <- a.post[t]/b.post[t]
    var.post[t] <- a.post[t]/(b.post[t]^2)
    icu.post[t] <- media.post[t] + norm_ic*sqrt(var.post[t])
    icl.post[t] <- media.post[t] - norm_ic*sqrt(var.post[t])

    # Preditiva em t = 1

    pred[t] <- a[t]/b[t]
    var.pred <- a[t]*(b[t]+1)/(b[t])^2
    icl.pred[t] <- qnbinom((1-IC_prob)/2, a[t], (b[t]/(b[t] +1)))
    icu.pred[t] <- qnbinom(1-(1-IC_prob)/2, a[t], (b[t]/(b[t] +1)))
  }

  mts <- mt
  Cts <- Ct

  var_index=matrix(apply(Ct,3,diag),n,T)!=0

  for(t in (T-1):1){
    var_ref=var_index[,t]
    restricted_Rt=Rt[var_ref,var_ref,t+1]
    restricted_Ct=Ct[var_ref,var_ref,t]
    simple_Rt_inv=restricted_Ct%*%t(G[var_ref,var_ref])%*%solve(restricted_Rt)

    mts[var_ref,t] <- mt[var_ref,t] + simple_Rt_inv%*%(mts[var_ref,t+1] - at[var_ref,t+1])
    Cts[var_ref,var_ref,t] <- restricted_Ct - simple_Rt_inv%*%(restricted_Rt - Cts[var_ref,var_ref,t+1])%*%t(simple_Rt_inv)
  }

  result <- list(mt,Ct,
                 ft, qt,
                 a,b,
                 a.post, b.post,
                 FF, G, D,W,
                 pred, icl.pred, icu.pred,
                 mts, Cts ,
                 IC_prob,var.pred,exp(pop),pop)
  names(result) <- c("mt",  "Ct",
                     "ft", "qt",
                     "alpha", "beta",
                     "alpha_star", "beta_star",
                     "F", "G", "D","W",
                     "pred", "icl.pred", "icu.pred",
                     "mts", "Cts",
                     "IC_prob","var.pred",'offset','log_offset')
  return(result)

}

poisson_testa_par <- function(y,m0 = 0, C0 = 1, F1,G1,D1,W1, pop, IC_prob=0.95){

  # Definindo quantidades
  n1 <- nrow(F1)

  r <- 1
  T <- length(y)
  n <- n1
  FF <- F1
  G <- G1

  # D1.aux <- matrix(rep(D1,n1^2), ncol = n1 )
  # D2.aux <- matrix(rep(D2,n2^2), ncol = n2 )
  # D.aux <- as.matrix(bdiag(D1.aux, D2.aux))
  D.aux <- D1
  D <- ifelse(D.aux == 0, 1, D.aux)

  W <- W1

  # Definindo objetos
  at <- matrix(0, ncol=T, nrow=n)
  mt <- matrix(0, ncol=T, nrow=n)
  ft <- matrix(0, ncol=1, nrow=T)
  qt <- matrix(0, ncol=1, nrow=T)
  alphat <- matrix(0, ncol=1, nrow=T)
  betat <- matrix(0, ncol=1, nrow=T)
  Ct <- array(rep(diag(n),T),dim=c(n,n,T))
  Rt <- array(rep(diag(n),T),dim=c(n,n,T))
  gt = pt = a = b= 0
  rep <- 1
  pred = var.pred = icl.pred = icu.pred = matrix(0, ncol=rep, nrow=T)
  media.post = var.post = icu.post = icl.post = matrix(0, ncol=rep, nrow=T)
  a.post = b.post = eqm =0
  E.th3 = E.th4 = raiz2 = raiz1= matrix(0, ncol=rep, nrow=T)
  tau0_star <- rep(NA, l = T)
  tau1_star <- rep(NA, l = T)
  tau0 <- rep(NA, l = T)
  tau1 <- rep(NA, l = T)
  fstar <- rep(NA, l = T)
  qstar <- rep(NA, l = T)

  norm_ic=qnorm(1-(1-IC_prob)/2)

  ## Algoritmo

  # Priori

  at[,1] <- G%*%m0
  Rt[,,1] <-G%*%C0%*%(t(G))*D[,,1]+W[,,1]

  reduc_RFF=Rt[,,1]%*%FF[,1]

  # Previsão 1 passo a frente
  ft[1,] <- t(FF[,1])%*%at[,1] + pop[1]
  qt[1,] <- t(FF[,1])%*%reduc_RFF

  # Compatibilizando prioris

  a[1] <- (1/qt[1,])
  b[1] <- (exp(-ft[1,] -0.5*qt[1,])/(qt[1,]))

  # Posteriori

  a.post[1] <- a[1] + y[1]
  b.post[1] <- b[1] + 1

  gt[1] <- log(a.post[1]/b.post[1]) + 1/(2*a.post[1])
  pt[1] <- (2*a.post[1]-1)/(2*a.post[1]^2)

  mt[,1] <- at[,1]+reduc_RFF*(gt[1]-ft[1,])*(1/(qt[1,]))
  Ct[,,1] <- Rt[,,1] - (reduc_RFF%*%t(reduc_RFF))*(1 - pt[1]/qt[1,])*(1/qt[1,])

  # Preditiva em t = 1

  pred[1] <- a[1]/ b[1]
  var.pred <- a[1]*(b[1]+1)/(b[1])^2

  # Passo 2 até t

  start<- proc.time()
  for(t in 2:T){

    # Priori

    at[,t] <-  G%*%mt[,t-1]
    Rt[,,t] <- G%*%Ct[,,t-1]%*%(t(G))*D[,,t]+W[,,t]

    reduc_RFF=Rt[,,t]%*%FF[,t]

    # Previsão 1 passo a frente

    ft[t,] <- t(FF[,t])%*%at[,t] + pop[t]
    qt[t,] <- t(FF[,t])%*%reduc_RFF

    # Compatibilizando prioris

    a[t] <- (1/qt[t,])
    b[t] <- (exp(-ft[t,] )/(qt[t,]))

    # Posteriori

    a.post[t] <- a[t] + y[t]
    b.post[t] <- b[t] + 1
    gt[t] <- log(a.post[t]/b.post[t])+1/(2*a.post[t])
    pt[t] <- (2*a.post[t]-1)/(2*a.post[t]^2)

    mt[,t] <- at[,t]+reduc_RFF*(gt[t]-ft[t,])*(1/(qt[t,]))
    Ct[,,t] <- Rt[,,t] - (reduc_RFF%*%t(reduc_RFF))*(1 - pt[t]/qt[t,])*(1/qt[t,])

    # Preditiva em t = 1

    pred[t] <- a[t]/b[t]
  }

  result <- list(mt,Ct,ft, qt, a,b, FF, G, D,W, pred,var.pred,exp(pop),pop)
  names(result) <- c("mt",  "Ct","ft", "qt",
                     "alpha", "beta","F", "G", "D","W","pred","var.pred",'offset','log_offset')
  return(result)

}
