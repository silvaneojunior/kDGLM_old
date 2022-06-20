library(MASS)
library(Matrix)
library(rootSolve)

normal.analise <- function(y,m01, C01,m02,C02, F1,F2,G1,G2,D1,D2){
  
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
  
  model_tau <- function(x, parms) c(
    F1 =  (q1 + f1^2)*exp(f2 + q2/2) - ((2*x[1] + 1)*x[3]^2)/(2*x[2]*(x[3]^2) - 8*(x[1]^2)*x[4]) - 1/(2*x[2])  ,
    F2 =  f1*exp(f2 + q2/2) + (2*x[1] + 1)*x[3]/(x[3]^2 - 4*x[2]*x[3]),
    F3 =  exp(f1 + 0.5*q2) - 4*x[2]*(x[1] + 1/2)/(x[3]^2 - 4*x[2]*x[4]),
    F4 =  f2 -digamma(x[1]+ 1/2) + log(x[3]^2/(4*x[2]) -x[4]))
  
  
  
  f <- function(m) {
    m[upper.tri(m)] <- t(m)[upper.tri(m)]
    m
  }
  
  # otim2 <- function(x, f1,f2,q1,q2,q12, Omega){
  #   tau1 <- x[1]
  #   tau2 <- x[2]
  #   tau0 <- x[3]
  #   media.log = log(1/(1+exp(f1) + exp(f2))) - 
  #     (q1/2)*(exp(f1)*exp(f2) + exp(f1))/(exp(f1)+exp(f2)+ 1)^2 - 
  #     (q2/2)*(exp(f2)*exp(f1) + exp(f2))/(exp(f1)+exp(f2)+1)^2 + 
  #     q12*(exp(f1 + f2)/(1+exp(f1)+exp(f2))) 
  #   
  #   eqs = rbind(
  #     f1  - digamma(tau1) + digamma(tau0 - tau1 - tau2),
  #     f2  - digamma(tau2) + digamma(tau0 - tau1 - tau2),
  #     digamma(tau0) - digamma(tau0 - tau1 - tau2) + media.log)
  #   
  #   return(t(eqs)%*%Omega%*%eqs)
  # }
  # 
  
  
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
  
  parms = c(f1,f2, q1, q2)
  
  ss1 <- multiroot(f = model_tau , start = c(0.3,0.3,0.6,0.3), parms = parms)

  
  tau0[1] <- ss1$root[1]
  tau1[1] <- ss1$root[2]
  tau2[1] <- ss1$root[3]
  tau3[1] <- ss1$root[4] 
  
  tau0_star[1]  <-  tau0[1] + 1/2
  tau1_star[1]  <-  tau1[1] - 1/2
  tau2_star[1]  <-  tau2[1] + y[1]
  tau3_star[1]  <-  tau3[1] - y[1]^2
  
  
  # Posteriori
  f1star[1,] <-   -2*tau2_star[1]/(2*tau1_star[1])
  f2star[1,] <-   digamma(tau0_star[1] + 1/2) - log(tau2_star[1]^2/(4*tau1_star[1]) - tau3_star[1])
  Q1star[1,] <-   (tau2_star[1]^2)/(4*tau1_star[1]^2) - 8*(tau1_star[1]^2)*(tau0_star[1] + 1/2)/(tau2_star[1]^2 - 4*tau1_star[1]*tau3_star[1]^2)
  Q2star[1,] <-   trigamma(tau0_star[1] + 0.5) + (digamma(tau0_star[1] + 0.5) - log(tau2_star[1]^2/(4*tau1_star[1]) - tau3_star[1]))^2
  Q12star[1,] <-   (-tau2_star[1]/(2*tau1_star[1]))*(digamma(tau0_star + 0.5) - log(tau2_star[1]^2/(4*tau1_star[1]) - tau3_star[1]))
  
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
  
    parms = c(f1,f2, q1, q2)
    
    ss1 <- multiroot(f = model_tau , start = c(0.01,0.01,0.01,0.01), parms = parms)
    
    
    tau0[t] <- ss1$root[1]
    tau1[t] <- ss1$root[2]
    tau2[t] <- ss1$root[3]
    tau3[t] <- ss1$root[4] 
    
    tau0_star[t]  <-  tau0[t] + 1/2
    tau1_star[t]  <-  tau1[t] - 1/2
    tau2_star[t]  <-  tau2[t] + y[t]
    tau3_star[t]  <-  tau3[t] - y[t]^2
    
    
    # Posteriori
    f1star[t,] <-   -2*tau2_star[t]/(2*tau1_star[t])
    f2star[t,] <-   digamma(tau0_star[t] + 1/2) - log(tau2_star[t]^2/(4*tau1_star[t]) - tau3_star[t])
    Q1star[t,] <-   (tau2_star[t]^2)/(4*tau1_star[t]^2) - 8*(tau1_star[t]^2)*(tau0_star[t] + 1/2)/(tau2_star[t]^2 - 4*tau1_star[t]*tau3_star[t]^2)
    Q2star[t,] <-   trigamma(tau0_star[t] + 0.5) + (digamma(tau0_star[t] + 0.5) - log(tau2_star[t]^2/(4*tau1_star[t]) - tau3_star[t]))^2
    Q12star[t,] <-   (-tau2_star[t]/(2*tau1_star[t]))*(digamma(tau0_star + 0.5) - log(tau2_star[t]^2/(4*tau1_star[t]) - tau3_star[t]))
    
    
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




# ibm ---------------------
m.ibmln <-read.table("~/Dropbox/Figuras Dissertacao/Antigos/GamaOmega/m-ibmln.txt", quote="\"", comment.char="")
y <- m.ibmln$V1

year <- seq(as.Date("1926/10/1"), by = "month", length.out = length(y))
par(mar=c(4,3,3,1))
plot(y,  pch = 20, xlab= "", ylab = "", cex = 0.6)
T <- length(y)
y <- y*100


n1 <- 1
n2 <- 1
T <- length(y)
r <- 1

# DEFININDO O VETOR F

F2 <- matrix(rep(1,n2),nrow=n2,ncol=T)
F1 <- matrix(rep(c(1),n1),nrow=n1,ncol=T)


# DEFININDO A MATRIZ G

# Bloco de tendęncia

G1 <- 1


G2 <- 1


## Passo t=0

m01 <- matrix(0,nrow=n1,ncol=r)
C01 <- matrix(c(100,0,0,10),n1,n1)

m02 <- matrix(0,nrow=n2,ncol=r)
C02 <- diag(1,n2,n2)


## Utilizando fator de desconto
D1 <- 1/0.85
D2 <- 1/0.90
resultados_phi_dinamico<- normal.analise(y, m01, C01,m02,C02,F1, F2, G1, G2, D1, D2)




xx <- seq(1, length(y))
m1 = resultados_phi_dinamico$mt1[1,]
ms1 = resultados_phi_dinamico$mts1[1,]
icu.mu <- resultados_phi_dinamico$mt1[1,] + 2*sqrt(resultados_phi_dinamico$C1[1,1,])
icl.mu <- resultados_phi_dinamico$mt1[1,] - 2*sqrt(resultados_phi_dinamico$C1[1,1,])




xx <- seq(1, length(y))
m2 = resultados_phi_dinamico$mt2[1,]
ms2 = resultados_phi_dinamico$mts2[1,]
icu.phi <- resultados_phi_dinamico$mt2[1,] + 2*sqrt(resultados_phi_dinamico$Cs2[1,1,])
icl.phi <- resultados_phi_dinamico$mt2[1,] - 2*sqrt(resultados_phi_dinamico$Cs2[1,1,])


f1 = resultados_phi_dinamico$f1
Q1 = resultados_phi_dinamico$Q1

xx <- seq(1,T,1)
iclpred <- f1 - 2*sqrt(Q1)
icupred <- f1 + 2*sqrt(Q1)

par(mfrow=c(1,1))
plot(y, type="l", ylab=" ",col=1, ylim=c(-30,40), xlab = " ") 
polygon(c(xx, rev(xx)),c(iclpred,rev(icupred)),col="grey80",border=NA)
points(y, lwd=2, pch = 20, cex = 0.6)
lines(m1, col=2, lwd =2)
legend("topleft", legend = c("Valor Verdadeiro", "Pred. 1 passo a frente"),
       lty=c(1,1), pch=c(NA, NA), lwd = c(2,2), col= c(1,2), bty = "n", cex =0.8)



plot(y, pch = 20, cex = 0.6, ylab=" ", ylim=c(-30,40), xlab = " ", col = "#949799") 
lines(f1, col="#e09690", lwd =2)
legend("topleft", legend = c("y", "one step forecast"),
       lty=c(NA,1), pch=c(20, NA), lwd = c(NA,2), col= c("#949799","#e09690"), bty = "n", cex =0.8)




library(plotly)

plot_ly(y = y, x = year) %>%
  add_markers(y = ~ y, name = "true value", marker = list(size = 3,color = "rgb(10, 10, 10)" )) %>%   
  add_trace(y = iclpred, type = 'scatter', mode = 'lines',fill = 'tonexty', fillcolor='rgba(232, 232, 232,0.8)',line = list(color = 'rgba(255,255,255,1)'),
            showlegend = FALSE, name = 'Low 2014')  %>% 
  add_trace(y = icupred, type = 'scatter', mode = 'lines',
            fill = 'tonexty', fillcolor='rgba(232, 232, 232,0.8)',line = list(color = 'rgba(255,255,255,1)'),
            showlegend = FALSE, name = 'Low 2014')  %>% 
  add_lines(y = f1, name = "Our proposal", line = list(color = 'rgb(133, 21, 11)', width = 1.2))  %>%
  layout(showlegend = FALSE,
         xaxis = list(title = " ", zeroline = TRUE,showline = TRUE),
         yaxis = list (title = " ", zeroline = FALSE,showline = TRUE))


plot_ly(y = m1, x = year, type = 'scatter', mode = 'lines',line = list(color = 'rgb(205, 12, 24)', width = 1.0)) %>%
  add_trace(y = icu.mu, type = 'scatter', mode = 'lines',
            fill = 'tonexty', fillcolor='rgba(232, 232, 232,0.8)',line = list(color = 'rgba(255,255,255,1)'),
            showlegend = FALSE, name = 'Low 2014')  %>% 
  add_trace(y = icl.mu, type = 'scatter', mode = 'lines',
            fill = 'tonexty', fillcolor='rgba(232, 232, 232,0.8)',line = list(color = 'rgba(255,255,255,1)'),
            showlegend = FALSE, name = 'Low 2014')  %>% 
  add_lines(y = m1, line = list(color = 'rgb(133, 21, 11)', width = 1.0))  %>% 
  add_lines(y = ms1, line = list(color = 'rgb(137, 120, 245)', width = 1.3, dash = 'dash'))  %>%
  layout(showlegend = FALSE,
         xaxis = list(title = " ", zeroline = FALSE,showline = TRUE),
         yaxis = list (title = " ", zeroline = FALSE,showline = TRUE))


plot_ly(y = m2, x = year, type = 'scatter', mode = 'lines',line = list(color = 'rgb(205, 12, 24)', width = 1.0)) %>%
  add_trace(y = icu.phi, type = 'scatter', mode = 'lines',
            fill = 'tonexty', fillcolor='rgba(232, 232, 232,0.8)',line = list(color = 'rgba(255,255,255,1)'),
            showlegend = FALSE, name = 'Low 2014')  %>% 
  add_trace(y = icl.phi, type = 'scatter', mode = 'lines',
            fill = 'tonexty', fillcolor='rgba(232, 232, 232,0.9)',line = list(color = 'rgba(255,255,255,1)'),
            showlegend = FALSE, name = 'Low 2014')  %>% 
  add_lines(y = m2, line = list(color = 'rgb(224, 150, 144)', width = 1))  %>% 
  add_lines(y = ms2, line = list(color = 'rgb(137, 120, 245)', width = 2.0))  %>%
  layout(showlegend = FALSE,
         xaxis = list(title = " ", showline = TRUE),
         yaxis = list (title = " "))

aux = (1/exp(2*resultados_phi_dinamico$mt2[1,]))
icu.phi <- 1/exp(resultados_phi_dinamico$mt2[1,]) + 2*sqrt(aux*resultados_phi_dinamico$C2[1,1,])
icl.phi <- 1/exp(resultados_phi_dinamico$mt2[1,]) - 2*sqrt(aux*resultados_phi_dinamico$C2[1,1,])
icl.phi <- ifelse(icl.phi <0, 0, icl.phi)

plot_ly(y = 1/exp(m2), x = year, type = 'scatter', mode = 'lines',line = list(color = 'rgb(205, 12, 24)', width = 1.0)) %>%
  add_trace(y = icu.phi, type = 'scatter', mode = 'lines',
            fill = 'tonexty', fillcolor='rgba(232, 232, 232,0.8)',line = list(color = 'rgba(255,255,255,1)'),
            showlegend = FALSE, name = 'Low 2014')  %>% 
  add_trace(y = icl.phi, type = 'scatter', mode = 'lines',
            fill = 'tonexty', fillcolor='rgba(232, 232, 232,0.9)',line = list(color = 'rgba(255,255,255,1)'),
            showlegend = FALSE, name = 'Low 2014')  %>% 
  add_lines(y = 1/exp(m2), line = list(color = 'rgb(133, 21, 11)', width = 1))  %>% 
  add_lines(y = 1/exp(ms2), line = list(color = 'rgb(68, 38, 255)', width = 1.3, dash = 'dash'))  %>%
  layout(showlegend = FALSE,
         xaxis = list(title = " ", zeroline = TRUE, showline = TRUE),
         yaxis = list (title = " ", zeroline = TRUE, showline = TRUE))


plot(1/exp(m2),lwd=2,type="l",xlab=" ",ylab=" ", col = "#E8E8E8")
polygon(c(xx, rev(xx)),c(icl.phi,rev(icu.phi)),col="#E8E8E8",border=NA)
lines(m2,lwd =1, col = "#e09690")
lines(ms2,lty=2,col="#85150b",lwd=2)
legend("topleft",cex =0.8,lwd =c(2,2),col = c("#e09690","#85150b"), lty = c(1,2), legend = c("Online mean", "Smoothed mean"), bty="n")


