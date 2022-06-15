library(readxl)
library(janitor)
library(ggplot2)
library(dplyr)
library(Matrix)
library(MASS)
library(rootSolve)

source("mult_maio_2021.R",encoding ='UTF-8')

rs_diarreia <- read_excel("brasil_diarreia-2.xlsx", 
                          sheet = "regiao_sudeste_v2") %>%  clean_names()
xx = rs_diarreia$data
rs_diarreia <-
  rs_diarreia %>% 
  mutate(x10_a_19_anos = x10_a_14_anos + x15_a_19_anos,
         x20_a_49_anos = x20_a_29_anos +  x30_a_39_anos  +  x40_a_49_anos,
         x50_anos_e_mais = x50_a_59_anos + x60_a_69_anos + x70_a_79_anos + x80_anos_e_mais,
         menor_que_4 = x1_a_4_anos + menor_1_ano,
         total_sem = total - x50_anos_e_mais - menor_que_4,
         outros = total - x50_anos_e_mais - menor_que_4)


ggplot(rs_diarreia) +
  geom_line(aes(x = data, y = x50_anos_e_mais, colour = "#52a19d"), size = 0.5) +
  geom_line(aes(x = data, y = menor_que_4, colour = "#db2757"), size = 0.5) +
  geom_line(aes(x = data, y = total_sem, colour = "#aac92a"), size = 0.5) +
  theme_minimal() + theme(legend.text = element_text(size=10),
                          legend.position="bottom") +
  ylim(0L, 3000L) + labs(y="Admissions", x = " ") + 
  scale_color_identity(name = "Age",
                       breaks = c( "#52a19d","#db2757", "#aac92a"),
                       labels = c( "More than 50 years" , " Less than 4 years", " Others"),
                       guide = "legend") 




ggplot(rs_diarreia) +
  geom_line(aes(x = data, y = menor_que_4, colour = "#db2757"), size = 0.5) +
  geom_line(aes(x = data, y = x50_anos_e_mais, colour = "#52a19d"), size = 0.5) +
  geom_vline(xintercept = as.numeric(rs_diarreia$data[109]),linetype=4) +
  theme_minimal() + theme(legend.text = element_text(size=10),
                          legend.position="bottom") +
  ylim(0L, 3000L) + labs(y="Internações", x = " ") + 
  scale_color_identity(name = "Age",
                       breaks = c( "#db2757", "#52a19d"),
                       labels = c( " <= 4 years" , " >= 50 years"),
                       guide = "legend") 



ggplot(rs_diarreia) +
  geom_line(aes(x = data, y = menor_que_4/total, colour = "#db2757"), size = 0.5) +
  geom_line(aes(x = data, y = x50_anos_e_mais/total, colour = "#52a19d"), size = 0.5) +
  geom_vline(xintercept = as.numeric(rs_diarreia$data[109]),linetype=4) +
  theme_minimal() + theme(legend.text = element_text(size=10),
                          legend.position="bottom") +
  ylim(0L, 1L) + labs(y="Internações", x = " ") + 
  scale_color_identity(name = "Age",
                       breaks = c( "#db2757", "#52a19d"),
                       labels = c(  " <4 years", " maior que 50" ),
                       guide = "legend") 


# Incluindo tendência e 1 harmonico ----------


y1 <- cbind(rs_diarreia$x50_anos_e_mais, rs_diarreia$menor_que_4, rs_diarreia$total_sem)

y <- trunc(y1)
T <- nrow(y)
n <- 4
m01 <- rep(0,n)
m02 = m01

C01 <- diag(50,n)
C02 = C01

F1 <- matrix(c(1,0,1,0),nrow=n,ncol=T)
F2 = F1

w  <- 2*pi/12
G.nivel <- 1
G.trend <- matrix(c(1,1,0,1), nrow=2, byrow =T)
G.saz <- matrix(c(cos(w),-sin(w),sin(w),cos(w)),2,2)
#G.saz2 <- matrix(c(cos(2*w),-sin(2*w),sin(2*w),cos(2*w)),2,2)
#G.saz3 <- matrix(c(cos(3*w),-sin(3*w),sin(3*w),cos(3*w)),2,2)

G1 <- as.matrix(bdiag( G.trend, G.saz))
G2 = G1

D.nivel <- matrix(1/0.9999,1,1)
D.trend <- matrix(1/0.95,2,2)
D.saz   <- matrix(1/0.975,2,2)
D1 <- as.matrix(bdiag(D.trend, D.saz))
D2 = D1

#total = rs_diarreia$total
total = rowSums((y))
resultados_mult <- mult.analise(y, m01, C01, m02, C02, F1,F2,G1,G2,D1,D2, N = total)

f1 <- resultados_mult$ft[,1]
f2 <- resultados_mult$ft[,2]
af = exp(f1)/(exp(f1) + exp(f2) + 1)
bf = exp(f2)/(exp(f1) + exp(f2) + 1)


par(mfrow = c(1,1))

par(mfrow = c(1,1))
year <- rs_diarreia$data
# Predição um passo a frente
plot(xx, rs_diarreia$x50_anos_e_mais/rs_diarreia$total, ylim =c(0,0.8), lty = 2, type = "p",
     ylab = "number of cases", xlab = "year", col = "#cccbca", pch = 20, cex = 0.4)
lines(xx,af, col = "#52a19d",  cex = 1,  lty = 1, lwd = 1.3)
points(xx,rs_diarreia$menor_que_4/rs_diarreia$total, cex = 0.4, pch = 20, col = "#cccbca")
lines(xx,bf, col =  "#db2757",  cex = 1,  lty = 1, lwd = 1.3)
# points(xx,rs_diarreia$total_sem/rs_diarreia$total, cex = 0.4, pch = 20, col = "#cccbca")
# lines(xx, 1-af - bf, col = "#aac92a",  cex = 1,  lty = 1, lwd = 1.3)




plot(xx, rs_diarreia$x50_anos_e_mais/rs_diarreia$total, ylim =c(0.1,0.5), lty = 2, type = "p", main = " >= 50", 
     ylab = "number of cases", xlab = "number of days", col = "#cccbca", pch = 20, cex = 0.4)
lines(xx,af, col = "#52a19d",  cex = 1,  lty = 1, lwd = 1.3)

plot(xx, rs_diarreia$menor_que_4/rs_diarreia$total, pch = 20, cex = 0.4, ylim =c(0,0.8), main = "< 4",  col = "#cccbca")
lines(xx,bf, col =  "#db2757",  cex = 1,  lty = 1, lwd = 1.3)

plot(xx,rs_diarreia$total_sem/rs_diarreia$total, pch = 20, cex = 0.4, ylim =c(0.,.6), main = "outros",  col = "#cccbca",)
lines(xx, 1-af - bf, col = "#aac92a",  cex = 1,  lty = 1, lwd = 1.3)


# Predição um passo a frente
plot(xx, rs_diarreia$x50_anos_e_mais, ylim =c(0,3000), lty = 2, type = "p",
     ylab = "number of cases", xlab = "year", col = "#cccbca", pch = 20, cex = 0.4)
lines(xx,af*rs_diarreia$total, col = "#52a19d",  cex = 1,  lty = 1, lwd = 1.3)
points(xx,rs_diarreia$menor_que_4, cex = 0.4, pch = 20, col = "#cccbca")
lines(xx,bf*rs_diarreia$total, col =  "#db2757",  cex = 1,  lty = 1, lwd = 1.3)
# points(xx,rs_diarreia$total_sem/rs_diarreia$total, cex = 0.4, pch = 20, col = "#cccbca")
# lines(xx, 1-af - bf, col = "#aac92a",  cex = 1,  lty = 1, lwd = 1.3)

plot_ly(y = rs_diarreia$x50_anos_e_mais, x = year) %>%
  add_lines(y =  rs_diarreia$x50_anos_e_mais, name = "true value", line = list(width = 1.2,color = "#52a19d" )) %>%   
  add_lines(y =  rs_diarreia$menor_que_4, name = "true value", line = list(width = 1.2,color = "#db2757" )) %>%   
  layout(showlegend = FALSE,
         xaxis = list(title = "year", zeroline = TRUE,showline = TRUE),
         yaxis = list (title = "# of admissions ", zeroline = FALSE,showline = TRUE))

plot_ly(y = rs_diarreia$x50_anos_e_mais, x = year) %>%
  add_markers(y =  rs_diarreia$x50_anos_e_mais, name = "true value", marker = list(size = 3,color = "#adadad" )) %>%   
  add_lines(y = af*rs_diarreia$total, name = "Our proposal", line = list(color = '#52a19d', width = 1.2))  %>%
  add_markers(y =  rs_diarreia$menor_que_4, name = "true value", marker = list(size = 3,color = "#adadad" )) %>%   
  add_lines(y = bf*rs_diarreia$total, name = "Our proposal", line = list(color = '#db2757', width = 1.2))  %>%
  layout(showlegend = FALSE,
         xaxis = list(title = "year", zeroline = TRUE,showline = TRUE),
         yaxis = list (title = "# of admissions ", zeroline = FALSE,showline = TRUE))

plot_ly(y = rs_diarreia$x50_anos_e_mais, x = year) %>%
  add_markers(y =  rs_diarreia$x50_anos_e_mais/total, name = "true value", marker = list(size = 3,color = "#adadad" )) %>%   
  add_lines(y = af, name = "Our proposal", line = list(color = '#52a19d', width = 1.2))  %>%
  add_markers(y =  rs_diarreia$menor_que_4/total, name = "true value", marker = list(size = 3,color = "#adadad" )) %>%   
  add_lines(y = bf, name = "Our proposal", line = list(color = '#db2757', width = 1.2))  %>%
  layout(showlegend = FALSE,
         xaxis = list(title = "year", zeroline = TRUE,showline = TRUE),
         yaxis = list (title = "proportion of admissions ", zeroline = FALSE,showline = TRUE))



    alpha1 <- resultados_mult$alpha_1_star
alpha2 <- resultados_mult$alpha_2_star
alpha3 <- resultados_mult$alpha_3_star
soma = alpha1 + alpha2+ alpha3

plot_ly(y = rs_diarreia$x50_anos_e_mais, x = year) %>%
  add_markers(y =  rs_diarreia$x50_anos_e_mais, name = "true value", marker = list(size = 3,color = "black" )) %>%   
  add_lines(y = (alpha1/soma)*rs_diarreia$total, name = "Our proposal", line = list(color = '#52a19d', width = 1.2))  %>%
  add_trace(y = qbinom(0.025, rs_diarreia$total, (alpha1/soma)) , type = 'scatter', mode = 'lines',fill = 'tonexty', fillcolor='rgba(232, 232, 232,0.8)',line = list(color = 'rgba(255,255,255,1)'),
            showlegend = FALSE, name = 'Low 2014')  %>% 
  add_trace(y = qbinom(0.975, rs_diarreia$total, (alpha1/soma)), type = 'scatter', mode = 'lines',
            fill = 'tonexty', fillcolor='rgba(232, 232, 232,0.8)',line = list(color = 'rgba(255,255,255,1)'),
            showlegend = FALSE, name = 'Low 2014')  %>% 
  add_lines(y = (alpha1/soma)*rs_diarreia$total, name = "Our proposal", line = list(color = '#52a19d', width = 1.2))  %>%
  # add_markers(y = rs_diarreia$x50_anos_e_mais, name = "Our proposal", line = list(color = '#52a19d', width = 1.2))  %>%
    layout(showlegend = FALSE,
         xaxis = list(title = "year", zeroline = TRUE,showline = TRUE),
         yaxis = list (title = "# of admissions ", zeroline = FALSE,showline = TRUE))





# resultados_mult$mts[,1:5] <- NA
# ms1 <- as.vector(resultados_mult$F[,1]%*%resultados_mult$mts)
# ms2 <- as.vector(resultados_mult$F[,2]%*%resultados_mult$mts)
# asf = exp(ms1)/(exp(ms1) + exp(ms2) + 1)
# bsf = exp(ms2)/(exp(ms1) + exp(ms2) + 1)
# 
# 
# plot(xx, rs_diarreia$x50_anos_e_mais/rs_diarreia$total, ylim =c(0.1,0.5), lty = 2, type = "p", main = " >= 50", 
#      ylab = "number of cases", xlab = "number of days", col = "#cccbca", pch = 20, cex = 0.4)
# lines(xx,asf, col = "#52a19d",  cex = 1,  lty = 1, lwd = 1.3)
# 
# plot(xx, rs_diarreia$menor_que_4/rs_diarreia$total, pch = 20, cex = 0.4, ylim =c(0,0.8), main = "< 4",  col = "#cccbca")
# lines(xx,bsf, col =  "#db2757",  cex = 1,  lty = 1, lwd = 1.3)
# 
# plot(xx,rs_diarreia$total_sem/rs_diarreia$total, pch = 20, cex = 0.4, ylim =c(0.,.6), main = "outros",  col = "#cccbca",)
# lines(xx, 1-asf - bsf, col = "#aac92a",  cex = 1,  lty = 1, lwd = 1.3)
# 
# 

resultados_mult$mts[,1:25] <- NA
ms1 <- as.vector(resultados_mult$F[1:2,1]%*%resultados_mult$mts[1:2,])
ms2 <- as.vector(resultados_mult$F[5:6,2]%*%resultados_mult$mts[5:6,])
asf = exp(ms1)/(exp(ms1) + exp(ms2) + 1)
bsf = exp(ms2)/(exp(ms1) + exp(ms2) + 1)


plot(xx, rs_diarreia$x50_anos_e_mais/rs_diarreia$total, ylim =c(0.1,0.5), lty = 2, type = "p", main = " >= 50", 
     ylab = "number of cases", xlab = "number of days", col = "#cccbca", pch = 20, cex = 0.4)
lines(xx,asf, col = "#52a19d",  cex = 1,  lty = 1, lwd = 1.3)

plot(xx, rs_diarreia$menor_que_4/rs_diarreia$total, pch = 20, cex = 0.4, ylim =c(0,0.8), main = "< 4",  col = "#cccbca")
lines(xx,bsf, col =  "#db2757",  cex = 1,  lty = 1, lwd = 1.3)

plot(xx,rs_diarreia$total_sem/rs_diarreia$total, pch = 20, cex = 0.4, ylim =c(0.,.6), main = "outros",  col = "#cccbca",)
lines(xx, 1-asf - bsf, col = "#aac92a",  cex = 1,  lty = 1, lwd = 1.3)


plot_ly(y = rs_diarreia$x50_anos_e_mais, x = year) %>%
  add_markers(y =  rs_diarreia$x50_anos_e_mais/rs_diarreia$total, name = "true value", marker = list(size = 3,color = "#adadad" )) %>%   
  add_lines(y = asf, name = "Our proposal", line = list(color = '#52a19d', width = 1.2))  %>%
  add_markers(y =  rs_diarreia$menor_que_4/rs_diarreia$total, name = "true value", marker = list(size = 3,color = "#adadad" )) %>%   
  add_lines(y = bsf, name = "Our proposal", line = list(color = '#db2757', width = 1.2))  %>%
  layout(showlegend = FALSE,
         xaxis = list(title = "year", zeroline = TRUE,showline = TRUE),
         yaxis = list (title = "proportion of admissions ", zeroline = FALSE,showline = TRUE))




resultados_mult$mts[,1:25] <- NA
mm1 <- as.vector(resultados_mult$F[,1]%*%resultados_mult$mt)
mm2 <- as.vector(resultados_mult$F[,2]%*%resultados_mult$mt)
amf = exp(mm1)/(exp(mm1) + exp(mm2) + 1)
bmf = exp(mm2)/(exp(mm1) + exp(mm2) + 1)

plot_ly(y = rs_diarreia$x50_anos_e_mais, x = year) %>%
  add_markers(y =  rs_diarreia$x50_anos_e_mais/rs_diarreia$total, name = "true value", marker = list(size = 3,color = "#adadad" )) %>%   
  add_lines(y = amf, name = "Our proposal", line = list(color = '#52a19d', width = 1.2))  %>%
  add_markers(y =  rs_diarreia$menor_que_4/rs_diarreia$total, name = "true value", marker = list(size = 3,color = "#adadad" )) %>%   
  add_lines(y = bmf, name = "Our proposal", line = list(color = '#db2757', width = 1.2))  %>%
  layout(showlegend = FALSE,
         xaxis = list(title = "year", zeroline = TRUE,showline = TRUE),
         yaxis = list (title = "proportion of admissions ", zeroline = FALSE,showline = TRUE))


plot(xx, amf, type = "o", lty = 1, col = "#52a19d",
     main = "online mean", ylab = " ", xlab = " ", ylim = c(0, 1.0), pch = 20, cex = 0.3)
lines(xx, bmf, type = "o", lty = 1, col = "#db2757",
      main = "nivel 50", ylab = " ", xlab = " ", ylim = c(0, 1.0), pch = 20, cex = 0.3)
# lines(xx, 1- amf - bmf, type = "o", lty = 1, col = "#aac92a",
#       main = "nivel 50", ylab = " ", xlab = " ", ylim = c(0, 1.0), pch = 20, cex = 0.3)

sup <-  resultados_mult$mt[1,] + 2*sqrt(resultados_mult$Ct[1,1,])
inf <-  resultados_mult$mt[1,] - 2*sqrt(resultados_mult$Ct[1,1,])
plot(xx, resultados_mult$mt[1,], type = "o", lty = 1, col = "#52a19d",
     main = "nivel 50", ylab = " ", xlab = " ", ylim = c(-1.5, 0.5), pch = 20, cex = 0.3)
# lines(xx, resultados_mult$mt[1,] + 2*sqrt(resultados_mult$Ct[1,1,]), col = "#52a19d",  cex = 1,  lty = 2, lwd = 1.3)
# lines(xx, resultados_mult$mt[1,] - 2*sqrt(resultados_mult$Ct[1,1,]), col = "#52a19d",  cex = 1,  lty = 2, lwd = 1.3)
polygon(c(xx, rev(xx)),c(inf,rev(sup)),col="#cccbca",border=NA)
lines(xx, resultados_mult$mt[1,], type = "o", lty = 1, col = "#52a19d",pch = 20, cex = 0.3 )




sup <-  resultados_mult$mt[2,] + 2*sqrt(resultados_mult$Ct[2,2,])
inf <-  resultados_mult$mt[2,] - 2*sqrt(resultados_mult$Ct[2,2,])
plot(xx, resultados_mult$mt[2,], type = "o", lty = 2, col = "#52a19d",
     main = "fator crescimento 50", ylab = " ", xlab = " ", ylim = c(-0.05, 0.05), pch = 20, cex = 0.3)
# lines(xx, resultados_mult$mt[2,] + 2*sqrt(resultados_mult$Ct[2,2,]), col = "#52a19d",  cex = 2,  lty = 2, lwd = 1.3)
# lines(xx, resultados_mult$mt[2,] - 2*sqrt(resultados_mult$Ct[2,2,]), col = "#52a19d",  cex = 2,  lty = 2, lwd = 1.3)
polygon(c(xx, rev(xx)),c(inf,rev(sup)),col="#cccbca",border=NA)
lines(xx, resultados_mult$mt[2,], type = "o", lty = 2, col = "#52a19d",pch = 20, cex = 0.3 )
abline(h = 0)


# plot(xx, resultados_mult$mt[2,], type = "l", lty = 1, col = "#52a19d", 
#      main = "fator crescimento 50", ylab = " ", xlab = " ", ylim = c(-0.5, 1.5))
# lines(xx, resultados_mult$mt[2,] + 2*sqrt(resultados_mult$Ct[2,2,]), col = "#52a19d",  cex = 1,  lty = 2, lwd = 1.3)
# lines(xx, resultados_mult$mt[2,] - 2*sqrt(resultados_mult$Ct[2,2,]), col = "#52a19d",  cex = 1,  lty = 2, lwd = 1.3)
# abline(h=0)




sup <-  resultados_mult$mt[3,] + 2*sqrt(resultados_mult$Ct[3,3,])
inf <-  resultados_mult$mt[3,] - 2*sqrt(resultados_mult$Ct[3,3,])
plot(xx, resultados_mult$mt[3,], type = "o", lty = 2, col = "#52a19d",
     main = "sazonalidade 50", ylab = " ", xlab = " ", ylim = c(-0.5, 0.5), pch = 20, cex = 0.3)
polygon(c(xx, rev(xx)),c(inf,rev(sup)),col="#cccbca",border=NA)
lines(xx, resultados_mult$mt[3,], type = "o", lty = 3, col = "#52a19d",pch = 20, cex = 0.3 )
abline(h = 0)


sup <-  resultados_mult$mt[5,] + 2*sqrt(resultados_mult$Ct[5,5,])
inf <-  resultados_mult$mt[5,] - 2*sqrt(resultados_mult$Ct[5,5,])
plot(xx, resultados_mult$mt[5,], type = "o", lty = 2, col = "#db2757",
     main = "nivel < 4", ylab = " ", xlab = " ", ylim = c(-0.5, 1.5), pch = 20, cex = 0.3)
polygon(c(xx, rev(xx)),c(inf,rev(sup)),col="#cccbca",border=NA)
lines(xx, resultados_mult$mt[5,], type = "o", lty = 5, col = "#db2757",pch = 20, cex = 0.3 )
abline(h = 0)


# plot(xx, resultados_mult$mt[5,], type = "o", lty = 1, col = "#db2757",
#      main = "nivel <4", ylab = " ", xlab = " ", ylim = c(-0.5, 1.5), pch = 20, cex = 0.3)
# lines(xx, resultados_mult$mt[5,] + 2*sqrt(resultados_mult$Ct[5,5,]), col = "#db2757",  cex = 1,  lty = 2, lwd = 1.3)
# lines(xx, resultados_mult$mt[5,] - 2*sqrt(resultados_mult$Ct[5,5,]), col = "#db2757",  cex = 1,  lty = 2, lwd = 1.3)
# abline(h=0)


sup <-  resultados_mult$mt[6,] + 2*sqrt(resultados_mult$Ct[6,6,])
inf <-  resultados_mult$mt[6,] - 2*sqrt(resultados_mult$Ct[6,6,])
plot(xx, resultados_mult$mt[6,], type = "o", lty = 2, col = "#db2757",
     main = "fator crescimento < 4", ylab = " ", xlab = " ", ylim = c(-0.05, 0.05), pch = 20, cex = 0.3)
polygon(c(xx, rev(xx)),c(inf,rev(sup)),col="#cccbca",border=NA)
lines(xx, resultados_mult$mt[6,], type = "o", lty = 6, col = "#db2757",pch = 20, cex = 0.3 )
abline(h = 0)



sup <-  resultados_mult$mt[7,] + 2*sqrt(resultados_mult$Ct[7,7,])
inf <-  resultados_mult$mt[7,] - 2*sqrt(resultados_mult$Ct[7,7,])
plot(xx, resultados_mult$mt[7,], type = "o", lty = 2, col = "#db2757",
     main = "sazonalidade < 4", ylab = " ", xlab = " ", ylim = c(-0.9, 0.9), pch = 20, cex = 0.3)
polygon(c(xx, rev(xx)),c(inf,rev(sup)),col="#cccbca",border=NA)
lines(xx, resultados_mult$mt[7,], type = "o", lty = 7, col = "#db2757",pch = 20, cex = 0.3 )
abline(h = 0)



# Suavizado

aux <- rep(0,length(xx))

sup <-  resultados_mult$mts[1,] + 2*sqrt(resultados_mult$Cts[1,1,])
inf <-  resultados_mult$mts[1,] - 2*sqrt(resultados_mult$Cts[1,1,])
plot(xx, resultados_mult$mts[1,], type = "o", lty = 1, col = "#52a19d",
     main = "nivel 50", ylab = " ", xlab = " ", ylim = c(-1.5, 0.5), pch = 20, cex = 0.3)
# lines(xx, resultados_mult$mts[1,] + 2*sqrt(resultados_mult$Cts[1,1,]), col = "#52a19d",  cex = 1,  lty = 2, lwd = 1.3)
# lines(xx, resultados_mult$mts[1,] - 2*sqrt(resultados_mult$Cts[1,1,]), col = "#52a19d",  cex = 1,  lty = 2, lwd = 1.3)
polygon(c(xx, rev(xx)),c(inf,rev(sup)),col="#cccbca",border=NA)
lines(xx, resultados_mult$mts[1,], type = "o", lty = 1, col = "#52a19d",pch = 20, cex = 0.3 )



plot_ly(y = resultados_mult$mts[1,], x = year) %>%
  add_markers(y = ~ resultados_mult$mts[1,], name = "true value", marker = list(size = 3,color = "#52a19d" )) %>%   
  add_trace(y = inf, type = 'scatter', mode = 'lines',fill = 'tonexty', fillcolor='rgba(232, 232, 232,0.8)',line = list(color = 'rgba(255,255,255,1)'),
            showlegend = FALSE, name = 'Low 2014')  %>% 
  add_trace(y = sup, type = 'scatter', mode = 'lines',
            fill = 'tonexty', fillcolor='rgba(232, 232, 232,0.8)',line = list(color = 'rgba(255,255,255,1)'),
            showlegend = FALSE, name = 'Low 2014')  %>% 
  add_lines(y = resultados_mult$mts[1,], name = "Our proposal", line = list(color = '#52a19d', width = 1.2))  %>%
  add_trace(y = 0,type = 'scatter', mode = 'lines',
            line = list(color = 'black', width = 0.8, dash="dot"),name = '') %>%
  layout(showlegend = FALSE,
         xaxis = list(title = " ", zeroline = TRUE,showline = TRUE),
         yaxis = list (title = " ", zeroline = FALSE,showline = TRUE))


sup <-  resultados_mult$mts[2,] + 2*sqrt(resultados_mult$Cts[2,2,])
inf <-  resultados_mult$mts[2,] - 2*sqrt(resultados_mult$Cts[2,2,])
plot(xx, resultados_mult$mts[2,], type = "o", lty = 2, col = "#52a19d",
     main = "fator crescimento 50", ylab = " ", xlab = " ", ylim = c(-0.05, 0.05), pch = 20, cex = 0.3)
# lines(xx, resultados_mult$mts[2,] + 2*sqrt(resultados_mult$Cts[2,2,]), col = "#52a19d",  cex = 2,  lty = 2, lwd = 1.3)
# lines(xx, resultados_mult$mts[2,] - 2*sqrt(resultados_mult$Cts[2,2,]), col = "#52a19d",  cex = 2,  lty = 2, lwd = 1.3)
polygon(c(xx, rev(xx)),c(inf,rev(sup)),col="#cccbca",border=NA)
lines(xx, resultados_mult$mts[2,], type = "o", lty = 2, col = "#52a19d",pch = 20, cex = 0.3 )
abline(h = 0)

plot_ly(y = resultados_mult$mts[2,], x = year) %>%
  add_markers(y = ~ resultados_mult$mts[2,], name = "true value", marker = list(size = 3,color = "#52a19d" )) %>%   
  add_trace(y = inf, type = 'scatter', mode = 'lines',fill = 'tonexty', fillcolor='rgba(232, 232, 232,0.8)',line = list(color = 'rgba(255,255,255,1)'),
            showlegend = FALSE, name = 'Low 2014')  %>% 
  add_trace(y = sup, type = 'scatter', mode = 'lines',
            fill = 'tonexty', fillcolor='rgba(232, 232, 232,0.8)',line = list(color = 'rgba(255,255,255,1)'),
            showlegend = FALSE, name = 'Low 2014')  %>% 
  add_lines(y = resultados_mult$mts[2,], name = "Our proposal", line = list(color = '#52a19d', width = 1.2))  %>%
  add_trace(y = 0,type = 'scatter', mode = 'lines',
            line = list(color = 'black', width = 0.8, dash="dot"),name = '') %>%
  layout(showlegend = FALSE,
         xaxis = list(title = " ", zeroline = TRUE,showline = TRUE),
         yaxis = list (title = " ", zeroline = FALSE,showline = TRUE))

# plot(xx, resultados_mult$mts[2,], type = "l", lty = 1, col = "#52a19d", 
#      main = "fator crescimento 50", ylab = " ", xlab = " ", ylim = c(-0.5, 1.5))
# lines(xx, resultados_mult$mts[2,] + 2*sqrt(resultados_mult$Cts[2,2,]), col = "#52a19d",  cex = 1,  lty = 2, lwd = 1.3)
# lines(xx, resultados_mult$mts[2,] - 2*sqrt(resultados_mult$Cts[2,2,]), col = "#52a19d",  cex = 1,  lty = 2, lwd = 1.3)
# abline(h=0)




sup <-  resultados_mult$mts[3,] + 2*sqrt(resultados_mult$Cts[3,3,])
inf <-  resultados_mult$mts[3,] - 2*sqrt(resultados_mult$Cts[3,3,])
plot(xx, resultados_mult$mts[3,], type = "o", lty = 2, col = "#52a19d",
     main = "sazonalidade 50", ylab = " ", xlab = " ", ylim = c(-0.5, 0.5), pch = 20, cex = 0.3)
polygon(c(xx, rev(xx)),c(inf,rev(sup)),col="#cccbca",border=NA)
lines(xx, resultados_mult$mts[3,], type = "o", lty = 3, col = "#52a19d",pch = 20, cex = 0.3 )
abline(h = 0)


plot_ly(y = resultados_mult$mts[3,], x = year) %>%
  add_markers(y = ~ resultados_mult$mts[3,], name = "true value", marker = list(size = 3,color = "#52a19d" )) %>%   
  add_trace(y = inf, type = 'scatter', mode = 'lines',fill = 'tonexty', fillcolor='rgba(232, 232, 232,0.8)',line = list(color = 'rgba(255,255,255,1)'),
            showlegend = FALSE, name = 'Low 2014')  %>% 
  add_trace(y = sup, type = 'scatter', mode = 'lines',
            fill = 'tonexty', fillcolor='rgba(232, 232, 232,0.8)',line = list(color = 'rgba(255,255,255,1)'),
            showlegend = FALSE, name = 'Low 2014')  %>% 
  add_lines(y = resultados_mult$mts[3,], name = "Our proposal", line = list(color = '#52a19d', width = 1.2))  %>%
  add_trace(y = 0,type = 'scatter', mode = 'lines',
            line = list(color = 'black', width = 0.8, dash="dot"),name = '') %>%
  layout(showlegend = FALSE,
         xaxis = list(title = " ", zeroline = TRUE,showline = TRUE),
         yaxis = list (title = " ", zeroline = FALSE,showline = TRUE, range = c(-0.3,0.3)))

sup <-  resultados_mult$mts[5,] + 2*sqrt(resultados_mult$Cts[5,5,])
inf <-  resultados_mult$mts[5,] - 2*sqrt(resultados_mult$Cts[5,5,])
plot(xx, resultados_mult$mts[5,], type = "o", lty = 2, col = "#db2757",
     main = "nivel < 4", ylab = " ", xlab = " ", ylim = c(-0.5, 1.5), pch = 20, cex = 0.3)
polygon(c(xx, rev(xx)),c(inf,rev(sup)),col="#cccbca",border=NA)
lines(xx, resultados_mult$mts[5,], type = "o", lty = 5, col = "#db2757",pch = 20, cex = 0.3 )
abline(h = 0)

plot_ly(y = resultados_mult$mts[5,], x = year) %>%
  add_markers(y = ~ resultados_mult$mts[5,], name = "true value", marker = list(size = 3,color = "#db2757" )) %>%   
  add_trace(y = inf, type = 'scatter', mode = 'lines',fill = 'tonexty', fillcolor='rgba(232, 232, 232,0.8)',line = list(color = 'rgba(255,255,255,1)'),
            showlegend = FALSE, name = 'Low 2014')  %>% 
  add_trace(y = sup, type = 'scatter', mode = 'lines',
            fill = 'tonexty', fillcolor='rgba(232, 232, 232,0.8)',line = list(color = 'rgba(255,255,255,1)'),
            showlegend = FALSE, name = 'Low 2014')  %>% 
  add_lines(y = resultados_mult$mts[5,], name = "Our proposal", line = list(color = '#db2757', width = 1.2))  %>%
  add_trace(y = 0,type = 'scatter', mode = 'lines',
            line = list(color = 'black', width = 0.8, dash="dot"),name = '') %>%
  layout(showlegend = FALSE,
         xaxis = list(title = " ", zeroline = TRUE,showline = TRUE),
         yaxis = list (title = " ", zeroline = FALSE,showline = TRUE))



# plot(xx, resultados_mult$mts[5,], type = "o", lty = 1, col = "#db2757",
#      main = "nivel <4", ylab = " ", xlab = " ", ylim = c(-0.5, 1.5), pch = 20, cex = 0.3)
# lines(xx, resultados_mult$mts[5,] + 2*sqrt(resultados_mult$Cts[5,5,]), col = "#db2757",  cex = 1,  lty = 2, lwd = 1.3)
# lines(xx, resultados_mult$mts[5,] - 2*sqrt(resultados_mult$Cts[5,5,]), col = "#db2757",  cex = 1,  lty = 2, lwd = 1.3)
# abline(h=0)


sup <-  resultados_mult$mts[6,] + 2*sqrt(resultados_mult$Cts[6,6,])
inf <-  resultados_mult$mts[6,] - 2*sqrt(resultados_mult$Cts[6,6,])
plot(xx, resultados_mult$mts[6,], type = "o", lty = 2, col = "#db2757",
     main = "fator crescimento < 4", ylab = " ", xlab = " ", ylim = c(-0.05, 0.05), pch = 20, cex = 0.3)
polygon(c(xx, rev(xx)),c(inf,rev(sup)),col="#cccbca",border=NA)
lines(xx, resultados_mult$mts[6,], type = "o", lty = 6, col = "#db2757",pch = 20, cex = 0.3 )
abline(h = 0)


plot_ly(y = resultados_mult$mts[6,], x = year) %>%
  add_markers(y = ~ resultados_mult$mts[6,], name = "true value", marker = list(size = 3,color = "#db2757" )) %>%   
  add_trace(y = inf, type = 'scatter', mode = 'lines',fill = 'tonexty', fillcolor='rgba(232, 232, 232,0.8)',line = list(color = 'rgba(255,255,255,1)'),
            showlegend = FALSE, name = 'Low 2014')  %>% 
  add_trace(y = sup, type = 'scatter', mode = 'lines',
            fill = 'tonexty', fillcolor='rgba(232, 232, 232,0.8)',line = list(color = 'rgba(255,255,255,1)'),
            showlegend = FALSE, name = 'Low 2014')  %>% 
  add_lines(y = resultados_mult$mts[6,], name = "Our proposal", line = list(color = '#db2757', width = 1.2))  %>%
  add_trace(y = 0,type = 'scatter', mode = 'lines',
            line = list(color = 'black', width = 0.8, dash="dot"),name = '') %>%
  layout(showlegend = FALSE,
         xaxis = list(title = " ", zeroline = TRUE,showline = TRUE),
         yaxis = list (title = " ", zeroline = FALSE,showline = TRUE))




sup <-  resultados_mult$mts[7,] + 2*sqrt(resultados_mult$Cts[7,7,])
inf <-  resultados_mult$mts[7,] - 2*sqrt(resultados_mult$Cts[7,7,])
plot(xx, resultados_mult$mts[7,], type = "o", lty = 2, col = "#db2757",
     main = "sazonalidade < 4", ylab = " ", xlab = " ", ylim = c(-0.9, 0.9), pch = 20, cex = 0.3)
polygon(c(xx, rev(xx)),c(inf,rev(sup)),col="#cccbca",border=NA)
lines(xx, resultados_mult$mts[7,], type = "o", lty = 7, col = "#db2757",pch = 20, cex = 0.3 )
abline(h = 0)

plot_ly(y = resultados_mult$mts[7,], x = year) %>%
  add_markers(y = ~ resultados_mult$mts[7,], name = "true value", marker = list(size = 3,color = "#db2757" )) %>%   
  add_trace(y = inf, type = 'scatter', mode = 'lines',fill = 'tonexty', fillcolor='rgba(232, 232, 232,0.8)',line = list(color = 'rgba(255,255,255,1)'),
            showlegend = FALSE, name = 'Low 2014')  %>% 
  add_trace(y = sup, type = 'scatter', mode = 'lines',
            fill = 'tonexty', fillcolor='rgba(232, 232, 232,0.8)',line = list(color = 'rgba(255,255,255,1)'),
            showlegend = FALSE, name = 'Low 2014')  %>% 
  add_lines(y = resultados_mult$mts[7,], name = "Our proposal", line = list(color = '#db2757', width = 1.2))  %>%
  add_trace(y = 0,type = 'scatter', mode = 'lines',
            line = list(color = 'black', width = 0.8, dash="dot"),name = '') %>%
  layout(showlegend = FALSE,
         xaxis = list(title = " ", zeroline = TRUE,showline = TRUE),
         yaxis = list (title = " ", zeroline = FALSE,showline = TRUE, range= c(-0.9, 0.9)))






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



