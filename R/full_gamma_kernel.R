#' system_full_gamma
#'
#' Evaluate the compatibilizing equation for the full gamma model (see Ref. Raíra).
#'
#' @param x vector: current tau values.
#' @param parms list: auxiliary values for the system.
#'
#' @return a vector with the values of the system (see Ref. Raíra).
system_full_gamma <- function(x, parms) {
  # A função que resolve o sistema tem um parâmetro que restringe x a valores positivos.
  # Restringir n,k a valores positivos através desse argumento pode ser numericamente mais estável em alguns casos.
  # Quando esta restrição estiver sendo feita, usar as linhas abaixo.
  # n=x[1]
  # k=x[2]
  # Quando a restrição não estiver sendo feita, usar as linhas abaixo.
  n <- exp(x[1])
  k <- exp(x[2])
  # tau=exp(x[3])
  # n=2
  # k=2

  # Calculando tau e theta dado n e k
  tau <- (n * parms$Hq1 + 1) / parms$Hq2
  theta <- n - k + n * log(tau / n) - (k + 1) / (2 * parms$Hq1)

  # (n * a / b + 1) * tau

  # Parâmetros da distribuição aproximada de phi.
  a <- (k + 1) / 2
  b <- (n - k + n * log(tau / n) - theta)

  # Se a<=5, calculamos a integral por quadratura gaussiana.
  # Se a>5 usamos aproximações de segunda ordem para obter uma expressão analítica para os valores esperados.
  # Em geral, a<=5, então usamos quadratura gaussiana praticamente em todos os casos.
  # Também testamos usar apenas quadratura gaussiana, mas os problemas persistiram.
  # print(n)
  # print(k)
  # print(a)
  # print(b)

  # n= 0.2322536
  # k= 1.865727
  # tau= 0.5087484
  # theta= -2.320433



  if (a <= 5) {
    # Densidade marginal aproximada de alpha (uso opcional).
    # f_densi=function(x){dgamma(a,b))}
    # c_val=1
    # Densidade marginal exata de phi.
    f_densi_raw <- function(x) {
      exp(k * (x + 1) * log(x) + lgamma(n * x + 1) + theta * x - k * lgamma(x + 1) - (n * x + 1) * log(x * tau))
    }
    c_val <- integrate(f_densi_raw, 0, Inf)$value
    f_densi <- function(x) {
      f_densi_raw(x) / c_val
    }

    # f=function(x){x*f_densi(x)}
    # Hp1=integrate(f,0,Inf)$value

    f <- function(x) {
      x * (digamma(x * n + 1) - log(x) - log(tau)) * f_densi(x)
    }
    Hp3 <- integrate(f, 0, Inf)$value

    f <- function(x) {
      (x * log(x) - lgamma(x)) * f_densi(x)
    }
    Hp4 <- integrate(f, 0, Inf)$value
  } else {
    c_val <- 1
    Hp3 <- log(tau / n) * a / b - 1 / n + b / (12 * (n**2) * (a - 1))
    Hp4 <- a / b + 0.5 * (digamma(a) - log(b)) - b / (12 * (a - 1)) - 11 / 12
  }

  f_all <- c(
    (Hp3 - parms$Hq3),
    (Hp4 - parms$Hq4)
  )
  # print(f_all)
  return(f_all)
}

#' system_full_gamma
#'
#' Evaluate the compatibilizing equation for the full gamma model (see Ref. Raíra).
#'
#' @param x vector: current tau values.
#' @param parms list: auxiliary values for the system.
#'
#' @return a vector with the values of the system (see Ref. Raíra).
system_full_gamma2 <- function(x, parms) {
  # A função que resolve o sistema tem um parâmetro que restringe x a valores positivos.
  # Restringir n,k a valores positivos através desse argumento pode ser numericamente mais estável em alguns casos.
  # Quando esta restrição estiver sendo feita, usar as linhas abaixo.
  # n=x[1]
  # k=x[2]
  # Quando a restrição não estiver sendo feita, usar as linhas abaixo.
  n <- exp(x)#exp(x[1])
  k <- n #exp(x[2])
  # tau=exp(x[3])
  # n=2
  # k=2

  # Calculando tau e theta dado n e k
  tau <- (n * parms$Hq1 + 1) / parms$Hq2
  theta <- n - k + n * log(tau / n) - (k + 1) / (2 * parms$Hq1)

  # (n * a / b + 1) * tau

  # Parâmetros da distribuição aproximada de phi.
  a <- (k + 1) / 2
  b <- (n - k + n * log(tau / n) - theta)

  # Se a<=5, calculamos a integral por quadratura gaussiana.
  # Se a>5 usamos aproximações de segunda ordem para obter uma expressão analítica para os valores esperados.
  # Em geral, a<=5, então usamos quadratura gaussiana praticamente em todos os casos.
  # Também testamos usar apenas quadratura gaussiana, mas os problemas persistiram.
  # print(n)
  # print(k)
  # print(a)
  # print(b)

  # n= 0.2322536
  # k= 1.865727
  # tau= 0.5087484
  # theta= -2.320433

  # print(a)
  # print(b)

  if (a <= 5) {
    # Densidade marginal aproximada de alpha (uso opcional).
    # f_densi=function(x){dgamma(a,b))}
    # c_val=1
    # Densidade marginal exata de phi.
    f_densi_raw <- function(x) {
      exp(k * (x + 1) * log(x) + lgamma(n * x + 1) + theta * x - k * lgamma(x + 1) - (n * x + 1) * log(x * tau))
    }
    c_val <- integrate(f_densi_raw, 0, Inf)$value
    f_densi <- function(x) {
      f_densi_raw(x) / c_val
    }

    f <- function(x) {
      x * (digamma(x * n + 1) - log(x) - log(tau)) * f_densi(x)
    }
    Hp3 <- integrate(f, 0, Inf)$value

    f <- function(x) {
      (x * log(x) - lgamma(x)) * f_densi(x)
    }
    Hp4 <- integrate(f, 0, Inf)$value
    Hp5=Hp3+Hp4


    # f <- function(x) {
    #   (lgamma(x)-x*digamma(n*x+1)) * f_densi(x)
    # }
    # Hp5 <- integrate(f, 0, Inf)$value+parms$Hq1*log(tau)
  } else {
    c_val <- 1
    Hp3 <- log(tau / n) * a / b - 1 / n + b / (12 * (n**2) * (a - 1))
    Hp4 <- a / b + 0.5 * (digamma(a) - log(b)) - b / (12 * (a - 1)) - 11 / 12
    Hp5=Hp3+Hp4
    # Hp5=parms$Hq1*(log(tau/n)-1)-0.5*digamma((n+1)/2)+0.5*log(n*log(tau/n)-theta)+(log(tau/n)-theta/n)/6+11/12-1/(2*n)
  }

  f_all <- c(
    (Hp5 - (parms$Hq3+parms$Hq4))
  )
  # print(f_all)
  return(f_all)
}

#' system_full_gamma
#'
#' Evaluate the compatibilizing equation for the full gamma model (see Ref. Raíra).
#'
#' @param x vector: current tau values.
#' @param parms list: auxiliary values for the system.
#'
#' @return a vector with the values of the system (see Ref. Raíra).
system_full_gamma3 <- function(x, parms) {
  # A função que resolve o sistema tem um parâmetro que restringe x a valores positivos.
  # Restringir n,k a valores positivos através desse argumento pode ser numericamente mais estável em alguns casos.
  # Quando esta restrição estiver sendo feita, usar as linhas abaixo.
  # n=x[1]
  # k=x[2]
  # Quando a restrição não estiver sendo feita, usar as linhas abaixo.
  n <- exp(x)#exp(x[1])
  k <- n #exp(x[2])
  # tau=exp(x[3])
  # n=2
  # k=2

  # Calculando tau e theta dado n e k
  tau <- (n * parms$Hq1 + 1) / parms$Hq2
  theta <- n - k + n * log(tau / n) - (k + 1) / (2 * parms$Hq1)

  # (n * a / b + 1) * tau

  # Parâmetros da distribuição aproximada de phi.
  a <- (k + 1) / 2
  b <- (n - k + n * log(tau / n) - theta)

  # Se a<=5, calculamos a integral por quadratura gaussiana.
  # Se a>5 usamos aproximações de segunda ordem para obter uma expressão analítica para os valores esperados.
  # Em geral, a<=5, então usamos quadratura gaussiana praticamente em todos os casos.
  # Também testamos usar apenas quadratura gaussiana, mas os problemas persistiram.
  # print(n)
  # print(k)
  # print(a)
  # print(b)

  # n= 0.2322536
  # k= 1.865727
  # tau= 0.5087484
  # theta= -2.320433

  # print(a)
  # print(b)

  if (a <= 5) {
    # Densidade marginal aproximada de alpha (uso opcional).
    # f_densi=function(x){dgamma(a,b))}
    # c_val=1
    # Densidade marginal exata de phi.
    f_densi_raw <- function(x) {
      exp(k * (x + 1) * log(x) + lgamma(n * x + 1) + theta * x - k * lgamma(x + 1) - (n * x + 1) * log(x * tau))
    }
    c_val <- integrate(f_densi_raw, 0, Inf)$value
    f_densi <- function(x) {
      f_densi_raw(x) / c_val
    }
    # print('a')
    f <- function(x) {
      (x * digamma(x * n + 1)  -x*log(x)- x*log(tau)) * f_densi(x)
    }
    Hp3 <- integrate(f, 0, Inf)$value

    # print('b')
    f <- function(x) {
      (x*log(x)- lgamma(x)) * f_densi(x)
    }
    Hp4 <- integrate(f, 0, Inf)$value

    # print('c')
    # f <- function(x) {
    #   (x * digamma(x * n + 1)  - lgamma(x)- x*log(tau)) * f_densi(x)
    # }
    # Hp5 <- integrate(f, 0, Inf)$value
    # print('sd')
    Hp5=Hp3+Hp4


    # f <- function(x) {
    #   (lgamma(x)-x*digamma(n*x+1)) * f_densi(x)
    # }
    # Hp5 <- integrate(f, 0, Inf)$value+parms$Hq1*log(tau)
  } else {
    c_val <- 1
    # Hp3 <- log(tau / n) * a / b - 1 / n + b / (12 * (n**2) * (a - 1))
    # Hp4 <- a / b + 0.5 * (digamma(a) - log(b)) - b / (12 * (a - 1)) - 11 / 12
    # Hp5=Hp3+Hp4
    Hp5=parms$Hq1*(log(tau/n)-1)-0.5*digamma((n+1)/2)+0.5*log(n*log(tau/n)-theta)+(log(tau/n)-theta/n)/6+11/12-1/(2*n)
  }

  f_all <- c(
    (Hp5 - (parms$Hq3+parms$Hq4))
  )
  # print(f_all)
  return(f_all)
}

#' convert_FGamma_Normal
#'
#' DESCRIPTION
#'
#' @param ft vector: A vector representing the means from the normal distribution.
#' @param Qt matrix: A matrix representing the covariance matrix of the normal distribution.
#' @param parms list: A list of extra known parameters of the distribuition. Not used in this function.
#'
#' @return The parameters of the conjugated distribuition of the linear predictor.
#' @export
convert_FGamma_Normal <- function(ft, Qt, parms) {

  # As contas foram feitas
  s <- exp(ft[2, ] + 3)
  f1 <- ft[1, ]
  f2 <- ft[2, ] - log(s)
  q1 <- Qt[1, 1]
  q2 <- Qt[2, 2]
  q12 <- Qt[1, 2]

  # Obs.: No relatório eu escrevi Hq3 com sinal trocado. Precisamos que f2+q12>0.
  Hq1 <- exp(f1 + q1 / 2)
  Hq2 <- exp(f1 - f2 + (q1 + q2 - 2 * q12) / 2)
  Hq3 <- -(f2 + q12) * Hq1

  # Valor de Hq4 usando aproximações de segunda ordem para a função log gamma.
  # Hq4=Hq1+0.5*f1-exp(-f1+q1/2)/12-11/12

  # Valor de Hq4 calculado por Monte Carlo.
  # Esse trecho é usado apenas para verificar se a aproximação acima está boa.
  # x=rnorm(20000,f1,sqrt(q1))
  # Hq4=mean(exp(x)*x-lgamma(exp(x)))

  # Valor de Hq4 calculado por quadratura gaussiana.
  # É preferível fazer a conta usando a densidade da log normal do que da normal.
  Hq4 <- integrate(function(x) {
    (x * log(x) - lgamma(x)) * dlnorm(x, f1, sqrt(q1))
  }, 0, Inf)$value

  parms <- list(
    "Hq1" = Hq1,
    "Hq2" = Hq2,
    "Hq3" = Hq3,
    "Hq4" = Hq4
  )

  # print(parms)
  # print(f1)
  # print(f2)
  # print(q1)
  # print(q2)
  # print(q12)
  ss1 <- multiroot(f = system_full_gamma, start = c(0,0), parms = parms) # , maxiter = 2000, atol = 10**-20)
  # A função que resolve o sistema tem um parâmetro que restringe x a valores positivos.
  # Restringir n,k a valores positivos através desse argumento pode ser numericamente mais estável em alguns casos.
  # Quando esta restrição estiver sendo feita, usar as linhas abaixo.
  # n=x[1]
  # k=x[2]
  # Quando a restrição não estiver sendo feita, usar as linhas abaixo.
  x <- as.numeric(ss1$root)
  n <- exp(x[1])
  k <- exp(x[2])

  # Calculando tau e theta dado n e k
  tau <- ((n * parms$Hq1 + 1) / parms$Hq2)
  # tau=exp(x[3])
  theta <- (n - k + n * log(tau / n) - (k + 1) / (2 * parms$Hq1))
  tau <- tau * s
  theta <- theta + n * log(s)
  return(list("n" = n, "k" = k, "tau" = tau, "theta" = theta))
}

#' convert_FGamma_Normal2
#'
#' DESCRIPTION
#'
#' @param ft vector: A vector representing the means from the normal distribution.
#' @param Qt matrix: A matrix representing the covariance matrix of the normal distribution.
#' @param parms list: A list of extra known parameters of the distribuition. Not used in this function.
#'
#' @return The parameters of the conjugated distribuition of the linear predictor.
#' @export
convert_FGamma_Normal2 <- function(ft, Qt, parms) {

  # As contas foram feitas
  s <- exp(ft[2, ] + 0)
  f1 <- ft[1, ]
  f2 <- ft[2, ] - log(s)
  q1 <- Qt[1, 1]
  q2 <- Qt[2, 2]
  q12 <- Qt[1, 2]

  # Obs.: No relatório eu escrevi Hq3 com sinal trocado. Precisamos que f2+q12>0.
  Hq1 <- exp(f1 + q1 / 2)
  Hq2 <- exp(f1 - f2 + (q1 + q2 - 2 * q12) / 2)
  Hq3 <- -(f2 + q12) * Hq1

  # Valor de Hq4 usando aproximações de segunda ordem para a função log gamma.
  # Hq4=Hq1+0.5*f1-exp(-f1+q1/2)/12-11/12

  # Valor de Hq4 calculado por Monte Carlo.
  # Esse trecho é usado apenas para verificar se a aproximação acima está boa.
  # x=rnorm(20000,f1,sqrt(q1))
  # Hq4=mean(exp(x)*x-lgamma(exp(x)))

  # Valor de Hq4 calculado por quadratura gaussiana.
  # É preferível fazer a conta usando a densidade da log normal do que da normal.
  Hq4 <- integrate(function(x) {
    (x * log(x) - lgamma(x)) * dlnorm(x, f1, sqrt(q1))
  }, 0, Inf)$value

  parms <- list(
    "Hq1" = Hq1,
    "Hq2" = Hq2,
    "Hq3" = Hq3,
    "Hq4" = Hq4
  )

  # print(parms)
  # print(f1)
  # print(f2)
  # print(q1)
  # print(q2)
  # print(q12)
  ss1 <- multiroot(f = system_full_gamma2, start = c(0), parms = parms) # , maxiter = 2000, atol = 10**-20)
  # A função que resolve o sistema tem um parâmetro que restringe x a valores positivos.
  # Restringir n,k a valores positivos através desse argumento pode ser numericamente mais estável em alguns casos.
  # Quando esta restrição estiver sendo feita, usar as linhas abaixo.
  # n=x[1]
  # k=x[2]
  # Quando a restrição não estiver sendo feita, usar as linhas abaixo.
  x <- as.numeric(ss1$root)
  n <- exp(x)#exp(x[1])
  k <- n#exp(x[2])

  # Calculando tau e theta dado n e k
  tau <- ((n * parms$Hq1 + 1) / parms$Hq2)
  # tau=exp(x[3])
  theta <- (n - k + n * log(tau / n) - (k + 1) / (2 * parms$Hq1))
  tau <- tau * s
  theta <- theta + n * log(s)
  return(list("n" = n, "k" = k, "tau" = tau, "theta" = theta))
}

#' convert_FGamma_Normal3
#'
#' DESCRIPTION
#'
#' @param ft vector: A vector representing the means from the normal distribution.
#' @param Qt matrix: A matrix representing the covariance matrix of the normal distribution.
#' @param parms list: A list of extra known parameters of the distribuition. Not used in this function.
#'
#' @return The parameters of the conjugated distribuition of the linear predictor.
#' @export
convert_FGamma_Normal3 <- function(ft, Qt, parms) {

  # As contas foram feitas
  s <- exp(ft[2, ] + 0)
  f1 <- ft[1, ]
  f2 <- ft[2, ] - log(s)
  q1 <- Qt[1, 1]
  q2 <- Qt[2, 2]
  q12 <- Qt[1, 2]

  # Obs.: No relatório eu escrevi Hq3 com sinal trocado. Precisamos que f2+q12>0.
  Hq1 <- exp(f1 + q1 / 2)
  Hq2 <- exp(f1 - f2 + (q1 + q2 - 2 * q12) / 2)
  Hq3 <- -(f2 + q12) * Hq1

  # Valor de Hq4 usando aproximações de segunda ordem para a função log gamma.
  # Hq4=Hq1+0.5*f1-exp(-f1+q1/2)/12-11/12

  # Valor de Hq4 calculado por Monte Carlo.
  # Esse trecho é usado apenas para verificar se a aproximação acima está boa.
  # x=rnorm(20000,f1,sqrt(q1))
  # Hq4=mean(exp(x)*x-lgamma(exp(x)))

  # Valor de Hq4 calculado por quadratura gaussiana.
  # É preferível fazer a conta usando a densidade da log normal do que da normal.
  Hq4 <- integrate(function(x) {
    (x * log(x) - lgamma(x)) * dlnorm(x, f1, sqrt(q1))
  }, 0, Inf)$value

  parms <- list(
    "Hq1" = Hq1,
    "Hq2" = Hq2,
    "Hq3" = Hq3,
    "Hq4" = Hq4
  )

  # print(parms)
  # print(f1)
  # print(f2)
  # print(q1)
  # print(q2)
  # print(q12)
  ss1 <- multiroot(f = system_full_gamma3, start = c(0), parms = parms) # , maxiter = 2000, atol = 10**-20)
  # A função que resolve o sistema tem um parâmetro que restringe x a valores positivos.
  # Restringir n,k a valores positivos através desse argumento pode ser numericamente mais estável em alguns casos.
  # Quando esta restrição estiver sendo feita, usar as linhas abaixo.
  # n=x[1]
  # k=x[2]
  # Quando a restrição não estiver sendo feita, usar as linhas abaixo.
  x <- as.numeric(ss1$root)
  n <- exp(x)#exp(x[1])
  k <- n#exp(x[2])

  # Calculando tau e theta dado n e k
  tau <- ((n * parms$Hq1 + 1) / parms$Hq2)
  # tau=exp(x[3])
  theta <- (n - k + n * log(tau / n) - (k + 1) / (2 * parms$Hq1))
  tau <- tau * s
  theta <- theta + n * log(s)
  return(list("n" = n, "k" = k, "tau" = tau, "theta" = theta))
}

#' convert_Normal_FGamma
#'
#' Calculate the parameters of the log-Normal that best approximates the given Inverse-Gamma distribuition.
#' The approximation is the best in the sense that it minimizes the KL divergence from the Inverse-Gamma to the log-Normal
#'
#' @param conj_prior list: A vector containing the parameters of the Inverse-Gamma (alpha,beta).
#' @param parms list: A list of extra known parameters of the distribuition. Not used in this function.
#'
#' @return The parameters of the Normal distribuition of the linear predictor.
#' @export
convert_Normal_FGamma <- function(conj_prior, parms) {
  n <- conj_prior$n
  tau <- conj_prior$tau
  theta <- conj_prior$theta
  k <- conj_prior$k

  s <- (10 * tau / n)

  tau <- tau / s
  theta <- theta - n * log(s)


  # print(n)
  # print(tau)
  # print(theta)
  # print(k)

  # Parâmetros da densidade aproximada de alpha
  a <- (k + 1) / 2
  b <- (n - k + n * log(tau / n) - theta)

  # print(a)
  # print(b)

  # print("a")
  # Comentar essas linhas caso a densidade aproximada seja usada.
  f_densi <- function(x) {
    exp(k * (x + 1) * log(x) + lgamma(n * x + 1) + theta * x - k * lgamma(x + 1) - (n * x + 1) * log(x * tau))
  }
  c_val <- integrate(f_densi, 0, Inf)$value
  # print("b")

  # Média 1 calculada com a densidade exata.
  f <- function(x) {
    log(x) * f_densi(x)
  }
  f1 <- integrate(f, 0, Inf)$value / c_val
  # Média 1 calculada com a densidade aproximada para phi.
  # f1=digamma(a)-log(b)

  # Média 2 calculada com a densidade exata.
  # Lembremos que mu|phi ~ IG(n*phi+1,phi*tau), logo E[log(mu)]=E[E[log(mu)|phi]]=E[digamma(n*phi+1)-log(tau*phi)]
  f <- function(x) {
    (-digamma(n * x + 1) + log(x * tau)) * f_densi(x)
  }
  f2 <- integrate(f, 0, Inf)$value / c_val
  # Média 2 calculada com a densidade aproximada
  # f2=log(tau/n)-(1/(2*n))*(b/(a-1))+(1/(12*n**2))*(b**2)/((a-1)*(a-2))

  # Variância 1 calculada com a densidade exata.
  f <- function(x) {
    (log(x)**2) * f_densi(x)
  }
  Q1 <- integrate(f, 0, Inf)$value / c_val - f1**2
  # Variância 1 calculada com a densidade aproximada.
  # Q1=trigamma(a)

  # Variância 2 calculada com a densidade exata.
  # O mesmo argumento para a média foi usado para o segundo momento.
  f <- function(x) {
    ((-digamma(n * x + 1) + log(x * tau))**2) * f_densi(x)
  }
  Q2 <- integrate(f, 0, Inf)$value / c_val - f2**2
  # Variância 2 calculada com a densidade aproximada.

  # Covariância calculada com a densidade exata.
  # O mesmo argumento para a média e para o segundo momento foi usado para a covariância.
  f <- function(x) {
    (log(x) - f1) * (-digamma(n * x + 1) + log(x * tau) - f2) * f_densi(x)
  }
  Q12 <- integrate(f, 0, Inf)$value / c_val
  # Covariância calculada com a densidade aproximada

  ft <- matrix(c(f1, f2 + log(s)), 2, 1)
  Qt <- matrix(c(Q1, Q12, Q12, Q2), 2, 2)
  # print("c")
  return(list("ft" = ft, "Qt" = Qt))
}

#' update_FGamma
#'
#' Calculate posterior parameter for the Inverser-Gamma, assuming that the observed values came from a Gamma model from which the shape parameter (phi) is known and the mean (mu) have prior distribuition Inverse-Gamma.
#'
#' @param conj_prior list: A vector containing the parameters of the Inverse-Gamma (alpha,beta).
#' @param y vector: A vector containing the observations.
#' @param parms list: A list of extra known parameters of the distribuition. For this kernel, parms should containg the shape parameter (phi) for the observational gamma model.
#'
#' @return The parameters of the posterior distribution.
#' @export
update_FGamma <- function(conj_prior, y, parms) {
  n0 <- conj_prior$n
  k0 <- conj_prior$k
  tau0 <- conj_prior$tau
  theta0 <- conj_prior$theta

  a <- (k0 + 1) / 2
  b <- (n0 - k0 + n0 * log(tau0 / n0) - theta0)

  # print(a)
  # print(b)

  n1 <- n0 + 1
  k1 <- k0 + 1
  tau1 <- tau0 + y
  theta1 <- theta0 + log(y)

  a <- (k1 + 1) / 2
  b <- (n1 - k1 + n1 * log(tau1 / n1) - theta1)

  # print(a)
  # print(b)

  return(list("n" = n1, "k" = k1, "tau" = tau1, "theta" = theta1))
}

#' Fgamma_pred
#'
#' Calculate the values for the predictive distribuition given the values of the parameter of the conjugated distribuition of the linear predictor.
#' The data is assumed to have Gamma distribuition with known shape parameter phi and it's mean having distribuition Inverse-Gamma with shape parameter a e rate parameter b.
#' In this scenario, the marginal distribuition of the data is Beta prime with parameters phi, alpha, beta / phi.
#'
#' @param conj_param List or data.frame: The parameters of the conjugated distribuitions of the linear predictor.
#' @param outcome Vector or matrix (optional): The observed values at the current time. Not used in this function.
#' @param parms List: A list of extra parameters for the model. For this function, it must contain the shape parameter phi of the observational model.
#' @param pred_cred Numeric: the desired credibility for the credibility interval.
#'
#' @return A list containing the following values:
#' \itemize{
#'    \item pred vector/matrix: the mean of the predictive distribuition of a next observation. Same type and shape as the parameter in model.
#'    \item var.pred vector/matrix: the variance of the predictive distribuition of a next observation. Same type and shape as the parameter in model.
#'    \item icl.pred vector/matrix: the percentile of 100*((1-pred_cred)/2)% of the predictive distribuition of a next observation. Same type and shape as the parameter in model.
#'    \item icu.pred vector/matrix: the percentile of 100*(1-(1-pred_cred)/2)% of the predictive distribuition of a next observation. Same type and shape as the parameter in model.
#' }
#'
#' @importFrom extraDistr qbetapr
#' @export
#'
#' @examples
#'
#' conj_param <- data.frame(
#'   "alpha" = c(1:3),
#'   "beta" = c(3:1)
#' )
#'
#' gamma_pred(conj_param, parms = list("phi" = 1))
Fgamma_pred <- function(conj_param, outcome = NULL, parms = list(), pred_cred = 0.95) {
  # DUMMY
  # Essa função não está pronta, é apenas um placeholder.
  alpha <- parms$alpha

  outcome_list <- list(
    "pred"     = alpha,
    "var.pred" = alpha,
    "icl.pred" = alpha,
    "icu.pred" = alpha
  )
  return(outcome_list)
}

#' Fgamma_log_like
#'
#' Calculate the predictive log-likelyhood of the data given the values of the parameter of the conjugated distribuition of the linear predictor.
#' The data is assumed to have Gamma distribuition with known shape parameter phi and it's mean having distribuition Inverse-Gamma with shape parameter a e rate parameter b.
#' In this scenario, the marginal distribuition of the data is Beta prime with parameters phi, alpha, beta / phi.
#'
#' @param conj_param List or data.frame: The parameters of the conjugated distribuition (Gamma) of the linear predictor.
#' @param outcome Vector or matrix: The observed values at the current time.
#' @param parms List: A list of extra parameters for the model. For this function, it must contain the shape parameter phi of the observational model.
#'
#' @return The predictive log-likelyhood of the observed data.
#' @export
#'
#' @importFrom extraDistr dbetapr
#' @examples
#'
#' conj_param <- data.frame(
#'   "alpha" = c(1:3),
#'   "beta" = c(3:1)
#' )
#'
#' gamma_log_like(conj_param, rgamma(3, c(1:3), c(3:1)), parms = list("phi" = 1))
Fgamma_log_like <- function(conj_param, outcome, parms = list()) {
  # DUMMY
  # Essa função não está pronta, é apenas um placeholder.
  return(outcome * 0)
}

#' @export
Fgamma_kernel <- list(
  "conj_prior" = convert_FGamma_Normal,
  "conj_post" = convert_Normal_FGamma,
  "update" = update_FGamma,
  "smoother" = generic_smoother,
  "pred" = Fgamma_pred,
  "log.like" = Fgamma_log_like,
  "offset" = log_offset_half,
  "link_function" = log,
  "inv_link_function" = exp,
  "param_names" = function(y) {
    c("n", "k", "tau", "theta")
  },
  "multi_var" = FALSE
)

#' @export
Fgamma2_kernel <- list(
  "conj_prior" = convert_FGamma_Normal2,
  "conj_post" = convert_Normal_FGamma,
  "update" = update_FGamma,
  "smoother" = generic_smoother,
  "pred" = Fgamma_pred,
  "log.like" = Fgamma_log_like,
  "offset" = log_offset_half,
  "link_function" = log,
  "inv_link_function" = exp,
  "param_names" = function(y) {
    c("n", "k", "tau", "theta")
  },
  "multi_var" = FALSE
)

#' @export
Fgamma3_kernel <- list(
  "conj_prior" = convert_FGamma_Normal3,
  "conj_post" = convert_Normal_FGamma,
  "update" = update_FGamma,
  "smoother" = generic_smoother,
  "pred" = Fgamma_pred,
  "log.like" = Fgamma_log_like,
  "offset" = log_offset_half,
  "link_function" = log,
  "inv_link_function" = exp,
  "param_names" = function(y) {
    c("n", "k", "tau", "theta")
  },
  "multi_var" = FALSE
)
