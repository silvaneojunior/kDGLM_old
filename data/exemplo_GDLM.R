# ## Se o pacote não estiver instalado
# devtools::install_github('silvaneojunior/GDLM')
#
# library(GDLM)
#
# ### Carregando dados ###
# dados=read.csv('dados_wania_exemplo.CSV')
# obitos=as.numeric(dados[1,-1])
# pac.dia=as.numeric(dados[2,-1])
# entradas=as.numeric(dados[6,-1])
# entradas.covid=as.numeric(dados[7,-1])
# T=length(obitos)
# intervencao=1:T>=54
#
# ### Criando estrutura para o modelo ###
#
# # Nível
# level=polynomial_block(order=1,D=1/0.9,name='nível')
# # Intervenção
# inter=polynomial_block(order=1,values=intervencao,C0=0,W=as.numeric(1:T==54),name='intervenção')
# # covid
# covid=polynomial_block(order=1,values=entradas.covid/entradas,name='Covid')
#
# ### Ajustando modelo ###
#
# fitted.model=fit_model(level,inter,covid,outcome=obitos,offset=pac.dia,family='Poisson')
#
# ### Visualizando ajuste ###
#
# show_fit(fitted.model)$plot
# plot_lat_var(fitted.model,'intervenção')$plot
# plot_lat_var(fitted.model,'')$plot
#
#
# ### Acessando algumas informações úteis ###
#
# fitted.model$mt # Média filtrada dos parâmetros latentes
# fitted.model$Ct # Variância filtrada dos parâmetros latentes
# fitted.model$ft # Média filtrada do preditor linear
# fitted.model$Qt # Variância filtrada do preditor linear
# fitted.model$conj_prior_param # Parâmetros da priori conjugada
# fitted.model$conj_post_param # Parâmetros da posteriori conjugada
# fitted.model$pred # Média da previsão um passo à frente
# fitted.model$var.pred # Variância da previsão um passo à frente
# fitted.model$icl.pred # Intervalo de credibilidade inferior para a previsão um passo à frente (sujeito a mudança de nome na próxima atualização do pacote).
# fitted.model$icu.pred # Intervalo de credibilidade superior para a previsão um passo à frente (sujeito a mudança de nome na próxima atualização do pacote).
# fitted.model$log.like # Verossimilhança da previsão um passo à frente
# fitted.model$mts # Média suavizada dos parâmetros latentes
# fitted.model$Cts # Variância suavizada dos parâmetros latentes
#
# ### Previsões para o futuro ###
#
# forecast(fitted.model,t=20)$plot
#
# ### Todas as funções foram documentadas. Mais detalhes sobre os argumentos de cada função estão presentes no help ###
