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

y
m0=matrix(c(m01,m02),2*n,1)
C0=bdiag(C01,C02)
F=array(bdiag(F1[,1],F2[,1]),c(2*n,2,length(y)))
G=bdiag(G1,G2)
D=array(bdiag(D1,D2),c(2*n,2*n,length(y)))
N=total

resultados=multinom_gi(
  y=y,
  m0=m0,
  C0=C0,
  FF=F,
  G=G,
  D=D,
  W=0,
  pop=N)


f1 <- resultados_mult$ft[,1]
f2 <- resultados_mult$ft[,2]
af_1 = exp(f1)/(exp(f1) + exp(f2) + 1)
bf_1 = exp(f2)/(exp(f1) + exp(f2) + 1)

f1 <- resultados$ft[,1]
f2 <- resultados$ft[,2]
af_2 = exp(f1)/(exp(f1) + exp(f2) + 1)
bf_2 = exp(f2)/(exp(f1) + exp(f2) + 1)

par(mfrow = c(1,1))
year <- rs_diarreia$data
# Predição um passo a frente
plot(xx, rs_diarreia$x50_anos_e_mais/rs_diarreia$total, ylim =c(0,0.8), lty = 2, type = "p",
     ylab = "number of cases", xlab = "year", col = "#cccbca", pch = 20, cex = 0.4)
lines(xx,af_1, col = "#52a19d",  cex = 1,  lty = 1, lwd = 1.3)
points(xx,rs_diarreia$menor_que_4/rs_diarreia$total, cex = 0.4, pch = 20, col = "#cccbca")
lines(xx,bf_1, col =  "#db2757",  cex = 1,  lty = 1, lwd = 1.3)
# points(xx,rs_diarreia$total_sem/rs_diarreia$total, cex = 0.4, pch = 20, col = "#cccbca")
# lines(xx, 1-af - bf, col = "#aac92a",  cex = 1,  lty = 1, lwd = 1.3)



