library(GDLM)

dados=read.csv('reports/Multinom/data/varicela internacoes.csv')[,c(1,7:162)]
dados[1:2,1]='00 a 04 anos'
dados[5:8,1]='15 a 49 anos'
dados[9:12,1]='50 anos e mais'
dados=aggregate(.~FaixaEtaria,dados,sum)
labels=dados[,1]
dados=dados[,-1]

pre_exp=read.csv2('reports/Multinom/data/populacao 2000-2020.csv')[-12,c(1,10:22)]
pre_exp[4:7,1]='15 a 49 anos'
pre_exp[8:11,1]='50 anos e mais'
pre_exp=aggregate(.~FaixaEtaria,pre_exp,sum)[,-1]

dummy=matrix(0,dim(pre_exp)[1],0)
nomes=c()
for(ano in c(2008:2020)){
  for(mes in c(1:12)){
    nomes=c(nomes,paste0('X',ano,'.',mes))
    dummy=cbind(dummy,pre_exp[,ano-2007])
  }
}
pre_exp=dummy
#idade_indice=3

proxy_vac=read.csv2('reports/Multinom/data/proxy_vac.csv')

vac_corb=dados[1,]*0
for(ano in c(2013:2019)){
  indices=substr(names(vac_corb),2,5)==as.character(ano)
  vac_corb[indices]=proxy_vac[ano-2012,2]
}
indices=substr(names(vac_corb),2,5)==as.character(2020)
vac_corb[indices]=proxy_vac[2019-2012,2]

t_offset=1
indice_inter=69
true_indice_inter=indice_inter+9

T_final=dim(dados)[2]


FaixaEtaria=4

out_var=4
ref_var=FaixaEtaria %>% as.numeric
indices=c(c(1:5)[c(1:5)!=ref_var],ref_var)

y=t(dados)
y=y[,indices]

ord_labels=labels[indices]

offset=pre_exp
offset=offset[indices,]

level=polynomial_block(2,1,k=4,D=1/0.9)
season=harmonic_block(12,1,k=4,D=1/0.98)
vac=polynomial_block(2,values=1:T>=true_indice_inter,k=4,D=1/0.9)
covid=polynomial_block(2,1,k=1:T>=146,D=1/0.9)

resultado=fit_model(level,season,
                        outcome=y,
                    offset=t(offset),
                        family='Multinomial')
show_fit(resultado)$plot
