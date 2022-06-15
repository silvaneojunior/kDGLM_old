library(ggplot2)
library(plotly)
library(Matrix)
library(rootSolve)
library(MASS)

source('norm_script.R')


# ibm ---------------------
m.ibmln <-read.table("m-ibmln.txt", quote="\"", comment.char="")
y <- m.ibmln$V1

year <- seq(as.Date("1926/10/1"), by = "month", length.out = length(y))
par(mar=c(4,3,3,1))
plot(y,  pch = 20, xlab= "", ylab = "", cex = 0.6)
T <- length(y)
y <- y


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
C01 <- matrix(c(1,0,0,1),n1,n1)

m02 <- matrix(0,nrow=n2,ncol=r)
C02 <- diag(1,n2,n2)


## Utilizando fator de desconto
D1 <- 1/0.85
D2 <- 1/0.90
resultados_phi_dinamico<- normal.analise(y, m01, C01,m02,C02,F1, F2, G1, G2, D1, D2)



plot(y,  pch = 20, cex = .4)
lines(resultados_phi_dinamico$mts1[1,], col = 4)
lines(resultados_phi_dinamico$f1 + 2*sqrt(resultados_phi_dinamico$Q1), col = 4, lty = 2)
lines(resultados_phi_dinamico$f1 - 2*sqrt(resultados_phi_dinamico$Q1), col = 4, lty = 2)
legend("topright", legend =c(expression(paste(y[t]), "f_pred")), 
       pch = c(20, NA), 
       col = c(1,4), 
       lty = c(NA,1), bty = "n", cex = 0.6)


plot(resultados_phi_dinamico$mts1[1,], col =4, type ="l", lwd = 1.5)
legend("topright", legend =c(expression(paste(mu[t]), m1[t])), 
       pch = c(20, NA), 
       col = c(1,4), 
       lty = c(NA,1), bty = "n", cex = 0.6)

plot(exp(resultados_phi_dinamico$mts2[1,]), ylim =c(0,4),  type ="l", col = 2)
lines(qgamma(0.975, rate = resultados_phi_dinamico$b, 
             shape= resultados_phi_dinamico$a), lty = 2, col = 2)
lines(qgamma(0.025, rate = resultados_phi_dinamico$b, 
             shape= resultados_phi_dinamico$a), lty = 2, col = 2)
legend("topright", legend =c(expression(paste(alpha[t]/beta[t]), phi[t])), 
       pch = c(NA, 20), 
       col = c(2,1), 
       lty = c(1,NA), bty = "n", cex = 0.6)





xx <- seq(1, length(y))
m1 = resultados_phi_dinamico$mt1[1,]
ms1 = resultados_phi_dinamico$mts1[1,]
icu.mu <- resultados_phi_dinamico$mt1[1,] + 2*sqrt(resultados_phi_dinamico$C1[1,1,])
icl.mu <- resultados_phi_dinamico$mt1[1,] - 2*sqrt(resultados_phi_dinamico$C1[1,1,])
par(mar=c(4,3,3,1)) 
plot(m1,lwd=2,type="l",xlab=" ",ylab=" ", ylim=c(-15,15), col = "#E8E8E8")
polygon(c(xx, rev(xx)),c(icl.mu,rev(icu.mu)),col="#E8E8E8",border=NA)
lines(m1,lwd =1, col = "#85150b")
lines(ms1,lty=2,col=2,lwd=1.4)
legend("topleft",cex =0.8,lwd =c(2,2),col = c("#e09690","#85150b"), lty = c(1,2), legend = c("Online mean", "Smoothed mean"), bty="n")




xx <- seq(1, length(y))
m2 = resultados_phi_dinamico$mt2[1,]
ms2 = resultados_phi_dinamico$mts2[1,]
icu.phi <- resultados_phi_dinamico$mt2[1,] + 2*sqrt(resultados_phi_dinamico$Cs2[1,1,])
icl.phi <- resultados_phi_dinamico$mt2[1,] - 2*sqrt(resultados_phi_dinamico$Cs2[1,1,])
par(mar=c(4,3,3,1))
plot(m2,lwd=2,type="l",xlab=" ",ylab=" ", ylim=c(-7,12), col = "#E8E8E8")
polygon(c(xx, rev(xx)),c(icl.phi,rev(icu.phi)),col="#E8E8E8",border=NA)
lines(m2,lwd =1, col = "#e09690")
lines(ms2,lty=2,col="#85150b",lwd=2)
legend("topleft",cex =0.8,lwd =c(2,2),col = c("#e09690","#85150b"), lty = c(1,2), legend = c("Online mean", "Smoothed mean"), bty="n")


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


