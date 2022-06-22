# Plotting stuff
library(ggplot2)
library(plotly)

# Manipulation  stuff
library(Matrix)
library(tidyr)
library(dplyr)
library(purrr)
library(stringr)
library(abind)

# Misc stuff
library(hms)
library(beepr)
library(latex2exp)

# Math stuff
library(rootSolve)
library(MASS)

#' generic_smoother
#'
#' Generic smoother for all models.
#'
#' @param mt Matrix: A matrix containing the filtered mean of the latent variables at each time. Each row should represent one variable.
#' @param Ct Array: A 3D-array representing the filtered covariance matrix of the latent variables at each time. The third dimension should represent the time index.
#' @param at Matrix: A matrix containing the one-step-ahead mean of the latent variables at each time based upon the filtered mean. Each row should represent one variable.
#' @param Rt Array: A 3D-array representing the one-step-ahead covariance matrix of the latent variables at each time based upon the filtered covariance matrix. The third dimension should represent the time index.
#' @param G  Matrix: A matrix representing the transition matrix of the model.
#'
#' @return List: The smoothed mean (mts) and covariance (Cts) of the latent variables at each time. Their dimension follows, respectivelly, the dimensions of mt and Ct.
#' @export
#'
#' @examples
#' T=20
#'
#' mt=matrix(c(cumsum(rnorm(T)+1),rep(1,T)),2,T,byrow=TRUE)
#' Ct=array(diag(c(1,1)),c(2,2,T))
#' G=matrix(c(1,0,1,1),2,2)
#' at=G%*%mt
#' Rt=array(G%*%t(G)+diag(c(0.1,0.1)),c(2,2,T))
#'
#' smoothed_values=generic_smoother(mt,Ct,at,Rt,G)
generic_smoother= function(mt,Ct,at,Rt,G){
  T=dim(mt)[2]
  n=dim(mt)[1]
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
  return(list('mts'=mts,'Cts'=Cts))
}

#' poisson_filter
#'
#' Filtering function for the Poisson model.
#'
#' @param y Vector: The observed value at time t.
#' @param m0 Vector: The prior mean for the latent vector at time t.
#' @param C0 Matrix: The prior covariance matrix for the latent vector at time t.
#' @param FF Matrix: The regression matrix at time t. Should be compatible with the dimensions of m0 and y, i.e., if m0 has dimension n and y has dimension m, FF should be n x m.
#' @param G Matrix: The state evolution matrix.
#' @param D Matrix: The discount factor matrix at time t.
#' @param W Matrix: The noise matrix at time t.
#' @param pop Vector: Same dimension as y. A vector contaning the offset at time t.
#'
#' @return A list containing:
#' \itemize{
#'  \item at: One-step-ahead mean for the latent vectors.
#'  \item Rt: One-step-ahead covariance matrix for the latent vectors.
#'  \item ft: One-step-ahead linear predictor.
#'  \item qt: One-step-ahead covariance matrix for the linear predictior.
#'  \item a: The alpha parameter of the compatibilized gamma prior (see Ref. Raíra).
#'  \item b: The beta parameter of the compatibilized gamma prior (see Ref. Raíra).
#'  \item a.post: The alpha parameter of the gamma posteirior (see Ref. Raíra).
#'  \item b.post: The beta parameter of the gamma posteirior (see Ref. Raíra).
#'  \item mt: The filtered mean of the latent vector at time t.
#'  \item Ct: The filtered covariance matrix of the latent vector at time t.
#'  \item y: The observed value at time t.
#'  }
#' @export
#'
#' @examples
#' y=5
#' m0=c(4,1)
#' C0=diag(c(1,1))
#' G=matrix(c(1,0,1,1),2,2)
#' FF=matrix(c(1,0),2,1)
#' D=diag(c(0.1,0.1),2,2)+1
#' W=diag(c(0,0),2,2)
#' offset=1
#'
#' filtered_data=poisson_filter(y,m0,C0,FF,G,D,W,offset)
poisson_filter = function(y,m0,C0,FF,G,D,W,offset){
  at <- G%*%m0
  Rt <-G%*%C0%*%(t(G))*D+W

  reduc_RFF=Rt%*%FF

  # One-step-ahead prediction
  ft <- t(FF)%*%at + log(offset)
  qt <- t(FF)%*%reduc_RFF

  # Compatibilizing priors

  a <- (1/qt)
  b <- (exp(-ft -0.5*qt)/(qt))

  # Calculating posterior

  a.post <- a + y
  b.post <- b + 1

  # Compatibilizing posterior

  gt <- digamma(a.post)-log(b.post)
  pt <- trigamma(a.post)

  mt <- at+reduc_RFF*as.vector((gt-ft)*(1/(qt)))
  if(length(qt)>1){
    Ct <- Rt - (reduc_RFF%*%t(reduc_RFF))%*%diag((1 - pt/qt)*(1/qt))
  }else{
    Ct <- Rt - (reduc_RFF%*%t(reduc_RFF))*as.vector((1 - pt/qt)*(1/qt))
  }

  return(list('at'=at,   'Rt'=Rt,
              'ft'=ft,   'Qt'=qt,
              'a'=a,     'b'=b,
              'a.post'=a.post,'b.post'=b.post,
              'mt'=mt,   'Ct'=Ct,
              'y'=y))
}

#' Title
#'
#' @param y
#' @param m0
#' @param C0
#' @param FF
#' @param G
#' @param D
#' @param W
#' @param pop
#'
#' @return
#'
#' @examples
poisson_LB_filter = function(y,m0,C0,FF,G,D,W,pop){
  at <- G%*%m0
  Rt <-G%*%C0%*%(t(G))*D+W

  reduc_RFF=Rt%*%FF

  # Previsão 1 passo a frente
  ft <- t(FF)%*%at + log(pop)
  qt <- t(FF)%*%reduc_RFF

  # Compatibilizando prioris

  a <- (1/qt)
  b <- (exp(-ft -0.5*qt)/(qt))

  # Posteriori

  a.post <- a + y
  b.post <- b + 1

  gt <- log(a.post/b.post) + 1/(2*a.post)
  pt <- (2*a.post-1)/(2*a.post^2)

  mt <- at+reduc_RFF*as.vector((gt-ft)*(1/(qt)))
  Ct <- Rt - (reduc_RFF%*%t(reduc_RFF))*as.vector((1 - pt/qt)*(1/qt))

  return(list('at'=at,   'Rt'=Rt,
              'ft'=ft,   'Qt'=qt,
              'a'=a,     'b'=b,
              'a.post'=a,'b.post'=b,
              'mt'=mt,   'Ct'=Ct,
              'y'=y))
}

poisson_pred=function(filter,IC_prob){
  a=filter$a
  b=filter$b
  list(
    'pred'     = a/ b,
    'var.pred' = a*(b+1)/(b)^2,
    'icl.pred' = qnbinom((1-IC_prob)/2, a, (b/(b +1))),
    'icu.pred' = qnbinom(1-(1-IC_prob)/2, a, (b/(b +1)))
  )
}

#' poisson_fit
#'
#' @param y Matrix: The observed data. It's dimension shoulb be T x m, where T is the length of the time series and m is the number of outcomes at each time.
#' @param m0 Vector: The prior mean for the latent vector.
#' @param C0 Matrix: The prior covariance matrix for the latent vector.
#' @param FF Array: A 3D-array containing the regression matrix for each time. It's dimension should be n x m x T, where n is the number of latent variables, m is the number of outcomes in the model and T is the time series length.
#' @param G Matrix: The state evolution matrix.
#' @param D Array: A 3D-array containing the discount factor matrix for each time. It's dimension should be n x n x T, where n is the number of latent variables and T is the time series length.
#' @param W Array: A 3D-array containing the covariance matrix of the noise for each time. It's dimension should be the same as D.
#' @param offset Matrix: The offset of the model. It's dimension should be the same as y.
#' @param IC_prob Numeric: Deprecated.
#'
#' @return A list containing the following values:
#' \itemize{
#'    \item mt Matrix: The filtered mean of the latent variables for each time. Dimensions are n x T.
#'    \item Ct Array: A 3D-array containing the filtered covariance matrix of the latent variable for each time. Dimensions are n x n x T.
#'    \item ft Matrix: The one-step-ahead linear predictor for each time. Dimensions are m x T.
#'    \item qt Array: A 3D-array containing the one-step-ahead covariance matrix for the linear predictor for each time. Dimensions are m x T.
#'    \item a Matrix: The alpha parameter for the gamma prior for each time. Dimensions are m x T.
#'    \item b Matrix: The beta parameter for the gamma prior for each time. Dimensions are m x T.
#'    \item a.post Matrix: The alpha parameter for the gamma posteior for each time. Dimensions are m x T.
#'    \item b.post Matrix: The beta parameter for the gamma posteior for each time. Dimensions are m x T.
#'    \item FF Array: The same as the argument (same values).
#'    \item G Matrix: The same as the argument (same values).
#'    \item D Array: The same as the argument (same values).
#'    \item W Array: The same as the argument (same values).
#'    \item pred Matrix: The one-step-ahead predictions for each time. Dimensions are m x T.
#'    \item var.pred Matrix: The variance for the one-step-ahead predictions for each time. Dimensions are m x T. Note that, in the multivariate Poisson case, the series are supossed independent, so, in particular, they are uncorrelated.
#'    \item icl.pred Matrix: The lower credibility interval for the prediction at each time. Dimensions are m x T.
#'    \item icu.pred Matrix: The upper credibility interval for the prediction at each time. Dimensions are m x T.
#'    \item mts Matrix: The smoothed mean of the latent variables for each time. Dimensions are n x T.
#'    \item Cts Array: A 3D-array containing the smoothed covariance matrix of the latent variable for each time. Dimensions are n x n x T.
#'    \item IC_prob Numeric: Deprecated
#'    \item offset Vector: The same as the argument (same values).
#'    \item log_offset Vector: The log offset.
#'    \item data_out Matrix: The same as the argument y (same values).
#' }
#' @export
#'
#' @examples
#' # Not ideal way: should use fit_model function.
#' T=200
#' w=(200/40)*2*pi
#' y= matrix(rpois(T,20*(sin(w*1:T/T)+2)),T,1)
#' m0=c(0,0,0)
#' C0=diag(c(1,1,1))
#' G=as.matrix(Matrix::bdiag(1,matrix(c(cos(w),sin(w),-sin(w),cos(w)),2,2)))
#' FF=array(matrix(c(1,1,0),3,1),c(3,1,T))
#' D=array(diag(c(0.1,0,0)),c(3,3,T))+1
#' W=array(diag(c(0,0,0)),c(3,3,T))
#' offset=matrix(1,T,1)
#'
#' filtered_data=GDLM::poisson_fit(y=y,m0 = m0, C0 = C0, FF=FF,G=G,D=D,W=W, offset=offset, IC_prob=0.95)
#'
#' plot(y)
#' lines(filtered_data$pred[1,])
poisson_fit <- function(y,m0 = 0, C0 = 1, FF,G,D,W, offset, IC_prob=0.95){
  T <- dim(y)[1]
  n <- dim(FF)[1]
  r <- dim(FF)[2]

  D.aux <- D
  D <- ifelse(D.aux == 0, 1, D.aux)

  m0 <- matrix(m0,n,1)
  C0 <- C0
  mt <- matrix(0,nrow=n,ncol=T)
  Ct <- array(rep(diag(n),T),dim=c(n,n,T))
  Rt <- array(rep(0,T),dim=c(n,n,T))
  ft <- matrix(0,nrow=T,ncol=r)
  at <- matrix(0,nrow=n,ncol=T)
  qt <-  array(0,dim=c(r,r,T))

  #r=1
  at <- matrix(0, nrow=n, ncol=T)
  mt <- matrix(0, nrow=n, ncol=T)
  ft <- matrix(0, nrow=r, ncol=T)
  qt <- matrix(0, nrow=r, ncol=T)
  Ct <- array(rep(diag(n),T),dim=c(n,n,T))
  Rt <- array(rep(diag(n),T),dim=c(n,n,T))
  pred = var.pred = icl.pred = icu.pred = matrix(0, nrow=r, ncol=T)
  media.post = var.post = icu.post = icl.post = matrix(0, nrow=r, ncol=T)
  a = b = matrix(0,r,T)
  a.post = b.post = matrix(0,r,T)

  norm_ic=qnorm(1-(1-IC_prob)/2)

  # Prior

  last_m=m0
  last_C=C0

  start<- proc.time()
  for(t in 1:T){
    filter=poisson_filter(y[t,],last_m,last_C,FF[,,t],G,D[,,t],W[,,t],offset[t,])

    at[,t]  <- filter$at
    Rt[,,t] <- filter$Rt
    ft[,t]  <- filter$ft
    qt[,t]  <- filter$Qt
    a[,t]  <- filter$a
    b[,t] <- filter$b
    a.post[,t]  <- filter$a.post
    b.post[,t] <- filter$b.post
    mt[,t]  <- filter$mt
    Ct[,,t] <- filter$Ct

    last_m=mt[,t]
    last_C=Ct[,,t]

    media.post[t] <- a.post[t]/(b.post[t])
    var.post[t] <- a.post[t]/(b.post[t]^2)
    icu.post[t] <- media.post[t] + norm_ic*sqrt(var.post[t])
    icl.post[t] <- media.post[t] - norm_ic*sqrt(var.post[t])

    prediction=poisson_pred(filter,IC_prob)

    pred[,t]     <- prediction$pred
    var.pred[,t] <- prediction$var.pred
    icl.pred[,t] <- prediction$icl.pred
    icu.pred[,t] <- prediction$icu.pred
  }

  smoothed=generic_smoother(mt,Ct,at,Rt,G)
  mts <- smoothed$mts
  Cts <- smoothed$Cts

  result <- list(mt,Ct,
                 ft, qt,
                 a,b,
                 a.post, b.post,
                 FF, G, D,W,
                 pred, var.pred, icl.pred, icu.pred,
                 mts, Cts ,
                 IC_prob,exp(offset),offset,
                 y)
  names(result) <- c("mt",  "Ct",
                     "ft", "qt",
                     "a", "b",
                     "a.post", "b.post",
                     "FF", "G", "D","W",
                     "pred", "var.pred", "icl.pred", "icu.pred",
                     "mts", "Cts",
                     "IC_prob",'offset','log_offset',
                     "data_out")
  return(result)

}

#' multnom_filter
#'
#' @param y
#' @param m0
#' @param C0
#' @param FF
#' @param G
#' @param D
#' @param W
#' @param pop
#'
#' @return
#' @importFrom(rootSolve,multiroot)
#' @export
#'
#' @examples
multnom_filter = function(y,m0,C0,FF,G,D,W,pop=NULL){
  r=dim(FF)[2]

  at = (G%*%m0)[,1]
  Pt <-  G%*%C0%*%(t(G))
  Rt <- as.matrix(D*Pt)+W

  # Previsao em t = 1
  ft <- (t(FF)%*%at)[,1]
  Qt <- as.matrix(t(FF)%*%Rt%*%FF)

  # minimizando...

  calc_helper=1 + sum(exp(ft))

  H=exp(ft)%*%t(exp(ft))/(calc_helper**2)
  diag(H)=-(exp(ft)*calc_helper-(exp(ft)**2))/(calc_helper**2)

  media.log =
    -log(calc_helper) + 0.5*(H%*%Qt) %>% diag %>% sum

  parms = list('ft'=ft, 'media.log'=media.log)

  ss1 <- multiroot(f = system_multinom , start = c(rep(0.01,r),0.01*(r+1)), parms = parms)

  tau <- as.numeric(ss1$root)

  alpha      <- tau
  alpha[r+1]      <- tau[r+1] - sum(tau[-r-1])

  alpha_star <-    alpha  +  y

  tau_star <-  alpha_star
  tau_star[r+1]  <-  sum(alpha_star)

  # Posteriori
  f_star <-   digamma(alpha_star[-r-1]) -  digamma(alpha_star[r+1])
  Q_star <-   matrix(trigamma(alpha_star[r+1]),r,r)
  diag(Q_star) <- trigamma(alpha_star[-r-1])+trigamma(alpha_star[r+1])

  At <- as.matrix(Rt%*%FF%*%ginv(Qt))
  mt <- at + At%*%(f_star -ft)
  Ct <- Rt +  At%*%(Q_star - Qt)%*%t(At)

  return(list('at'=at,   'Rt'=Rt,
              'ft'=ft,   'Qt'=Qt,
              'tau'=tau,     'tau_star'=tau_star,
              'alpha'=alpha,     'alpha_star'=alpha_star,
              'f_star'=f_star,'Q_star'=Q_star,
              'At'=At,
              'mt'=mt,   'Ct'=Ct,
              'y'=y))
}

multnom_pred=function(filter,IC_prob){
  r=length(filter$ft)
  pre_ps=exp(filter$ft)/(sum(exp(filter$ft))+1)

  p=pre_ps
  var=filter$Qt

  diag_mult=diag(p*(1-p))
  cov_mult=diag_mult%*%var%*%diag_mult

  vec_rest=1-sum(p)

  n_total=sum(filter$y)

  p=c(p,vec_rest)*n_total
  trans_mat=rbind(diag(r),rep(-1,r))
  var=(trans_mat%*%cov_mult%*%t(trans_mat))%*%(diag(r+1)*(n_total**2))

  p_i=p-2*sqrt(diag(var))
  p_s=p+2*sqrt(diag(var))

  list(
    'pred'     = p,
    'var.pred' = var,
    'icl.pred' = p_i,
    'icu.pred' = p_s
  )
}

system_multinom <- function(x, parms){
  sub_last=digamma(x[length(x)] - sum(x[1:(length(x)-1)]))
  digamma_vec=digamma(x)

  f_all=parms$ft-digamma_vec[-length(x)]+sub_last
  last_guy=parms$media.log+digamma_vec[length(x)]-sub_last

  f_all=c(f_all,last_guy)

  return(f_all)
}

multnom_fit <- function(y,m0=0, C0=1, FF,G,D,W, pop, IC_prob=0.95){
  # Definindo quantidades
  T <- nrow(y)
  n <- dim(FF)[1]

  D.aux <- D
  D <- ifelse(D.aux == 0, 1, D.aux)

  r = ncol(y)-1
  m0 <- matrix(m0,n,1)
  C0 <- C0
  mt <- matrix(0,nrow=n,ncol=T)
  Ct <- array(rep(diag(n),T),dim=c(n,n,T))
  Rt <- array(rep(0,T),dim=c(n,n,T))
  ft <- matrix(0,nrow=T,ncol=r)
  at <- matrix(0,nrow=n,ncol=T)
  Qt <-  array(0,dim=c(r,r,T))
  At <- array(0,dim=c(n,r,T))

  pred <-  matrix(0,r+1,T)
  var.pred <- array(0,c(r+1,r+1,T))
  icl.pred <- matrix(0,r+1,T)
  icu.pred <- matrix(0,r+1,T)

  f_star <- matrix(0,nrow=T,ncol=r)
  Q_star <- array(0,c(r,r,T))

  mt <- matrix(0,nrow=n,ncol=T)
  Ct <- array(rep(diag(n),T),dim=c(n,n,T))
  tau <- matrix(NA,nrow=r+1,ncol=T)
  alpha <- matrix(NA,nrow=r+1,ncol=T)
  alpha_star <- matrix(NA,nrow=r+1,ncol=T)
  tau_star = matrix(NA,nrow=r+1,ncol=T)

  D=ifelse(D==0,1,D)

  last_m=m0
  last_C=C0

  for(t in 1:T){
    filter=multnom_filter(y[t,],last_m,last_C,FF[,,t] %>% matrix(n,r),G,D[,,t],W[,,t])

    at[,t]         <- filter$at
    Rt[,,t]        <- filter$Rt
    ft[t,]         <- filter$ft
    Qt[,,t]        <- filter$Qt
    tau[,t]        <- filter$tau
    alpha[,t]      <- filter$alpha
    alpha_star[,t] <- filter$alpha_star
    tau[,t]        <- filter$tau
    tau_star[,t]   <- filter$tau_star
    f_star[t,]     <- filter$f_star
    Q_star[,,t]    <- filter$Q_star
    At[,,t]        <- filter$At
    mt[,t]         <- filter$mt
    Ct[,,t]        <- filter$Ct

    last_m=mt[,t]
    last_C=Ct[,,t]

    prediction=multnom_pred(filter,IC_prob)

    pred[,t]      <- prediction$pred
    var.pred[,,t] <- prediction$var.pred
    icl.pred[,t]  <- prediction$icl.pred
    icu.pred[,t]  <- prediction$icu.pred
  }

  smoothed=generic_smoother(mt,Ct,at,Rt,G)
  mts <- smoothed$mts
  Cts <- smoothed$Cts

  result <- list(mt,Ct,
                 ft, Qt,
                 f_star, Q_star,
                 alpha,alpha_star,
                 tau,tau_star,
                 FF, G, D,W,
                 mts, Cts ,
                 pred, var.pred, icl.pred, icu.pred,
                 exp(pop),pop,
                 y)
  names(result) <- c("mt",  "Ct",
                     "ft", "Qt",
                     'f_star', 'Q_star',
                     'alpha','alpha_star',
                     'tau','tau_star',
                     "FF", "G", "D","W",
                     "mts", "Cts",
                     "pred", "var.pred", "icl.pred", "icu.pred",
                     'offset','log_offset',
                     "data_out")
  return(result)
}

normal_filter = function(y,m0,C0,FF,G,D,W,pop=NULL){
  r=2

  at = (G%*%m0)
  Pt <-  G%*%C0%*%(t(G))
  Rt <- as.matrix(D*Pt)+W

  # Previsao em t = 1
  ft <- (t(FF)%*%at)
  Qt <- as.matrix(t(FF)%*%Rt%*%FF)

  # minimizando...

  f1=ft[1]
  f2=ft[2]
  q1 = Qt[1,1]
  q2 = Qt[2,2]
  q12 = Qt[1,2]

  c0=1/(exp(f2 + q2/2)*q1)
  mu0=f1

  # n0=2/q2
  parms = c(q2)
  model_tau <- function(x, parms){
    output=c(
      F4 =  -digamma(exp(x[1])/2)+x[1]-q2/2-log(2)
    )
    return(output)
  }
  jacfunc=function(x,parms){
    jac=matrix(c(-0.5*exp(x[1])*trigamma(exp(x[1])/2)+1),
               1,1)
    return(jac)
  }

  ss1 <- multiroot(f = model_tau , start = c(0), parms = parms,jacfunc = jacfunc, jactype='fullusr', rtol = 1e-40)
  n0=exp(ss1$root[1])

  #print(exp(-f2-q2/2))
  d0=n0/exp(f2+q2/2)

  tau0=n0/2-1/2
  tau1=-c0/2
  tau2=c0*mu0
  tau3=-(c0*(mu0**2)+d0)/2

  # parms = c(f1,f2, q1, q2)
  #
  # model_tau <- function(x, parms) c(
  #   F1 =  (q1 + f1^2)*exp(f2 + q2/2) - ((2*x[1] + 1)*x[3]^2)/(2*x[2]*(x[3]^2) - 8*(x[2]^2)*x[4]) + 1/(2*x[2])  ,
  #   F2 =  f1*exp(f2 + q2/2) + (2*x[1] + 1)*x[3]/(x[3]^2 - 4*x[2]*x[4]),
  #   F3 =  exp(f2 + 0.5*q2) - 4*x[2]*(x[1] + 1/2)/(x[3]^2 - 4*x[2]*x[4]),
  #   F4 =  f2 -digamma(x[1]+ 1/2) + log(x[3]^2/(4*x[2]) -x[4]))
  #
  # print(model_tau(x=c(tau0,tau1,tau2,tau3),parms=parms))

  tau0_star  <-  tau0 + 1/2
  tau1_star  <-  tau1 - 1/2
  tau2_star  <-  tau2 + y[1]
  tau3_star  <-  tau3 - y[1]^2

  # Posteriori
  f1star <-   -tau2_star/(2*tau1_star)
  f2star <-   digamma(tau0_star + 0.5) -log(((tau2_star^2)/(4*tau1_star)) - tau3_star)
  Q1star <-   f1star**2 - 8*(tau1_star^2)*(tau0_star + 1/2)/(tau2_star^2 - 4*tau1_star*(tau3_star))
  Q2star <-   trigamma(tau0_star + 0.5) + f2star**2
  Q12star<-   f1star*f2star

  fstar <- c(f1star,   f2star)
  Qstar <- matrix(c( Q1star,  Q12star,  Q12star,  Q2star), byrow =F, ncol = 2)

  At <- Rt%*%FF%*%ginv(Qt)
  mt <- at + At%*%(fstar -ft)
  Ct <- Rt +  At%*%(Qstar - Qt)%*%t(At)

  return(list('at'=at,   'Rt'=Rt,
              'ft'=ft,   'Qt'=Qt,
              'tau'=c(tau0,tau1,tau2,tau3),     'tau_star'=c(tau0_star,tau1_star,tau2_star,tau3_star),
              'f_star'=fstar,'Q_star'=Qstar,
              'At'=At,
              'mt'=mt,   'Ct'=Ct,
              'y'=y))
}

normal_pred=function(filter,IC_prob){
  c0=-2*filter$tau[2]
  mu0=filter$tau[3]/c0
  alpha=(filter$tau[3]**2)/(4*filter$tau[2])-filter$tau[4]
  beta=filter$tau[1]+0.5

  mu=mu0
  nu=2*alpha
  sigma2=beta/(c0*alpha)
  # print(nu)

  pred=mu
  var.pred=sigma2*nu/(nu-2)

  icl.pred=qt((1-IC_prob)/2,nu)*sqrt(sigma2)+mu
  icu.pred=qt(1-(1-IC_prob)/2,nu)*sqrt(sigma2)+mu

  list(
    'pred'     = pred,
    'var.pred' = var.pred,
    'icl.pred' = icl.pred,
    'icu.pred' = icu.pred
  )
}

normal_fit <- function(y,m0=0, C0=1, FF,G,D,W, pop, IC_prob=0.95){

  # Definindo quantidades
  T <- nrow(y)
  n <- dim(FF)[1]

  D.aux <- D
  D <- ifelse(D.aux == 0, 1, D.aux)
  r = 2
  m0 <- matrix(m0,n,1)
  C0 <- C0
  mt <- matrix(0,nrow=n,ncol=T)
  Ct <- array(rep(diag(n),T),dim=c(n,n,T))
  Rt <- array(rep(0,T),dim=c(n,n,T))
  ft <- matrix(0,nrow=T,ncol=r)
  at <- matrix(0,nrow=n,ncol=T)
  Qt <-  array(0,dim=c(r,r,T))
  At <- array(0,dim=c(n,r,T))

  pred <-  matrix(0,r-1,T)
  var.pred <- array(0,c(r-1,r-1,T))
  icl.pred <- matrix(0,r-1,T)
  icu.pred <- matrix(0,r-1,T)

  f_star <- matrix(0,nrow=T,ncol=r)
  Q_star <- array(0,c(r,r,T))

  mt <- matrix(0,nrow=n,ncol=T)
  Ct <- array(rep(diag(n),T),dim=c(n,n,T))
  tau <- matrix(NA,nrow=4,ncol=T)
  tau_star <- matrix(NA,nrow=4,ncol=T)

  D=ifelse(D==0,1,D)

  last_m=m0
  last_C=C0

  for(t in 1:T){
    filter=normal_filter(y[t,],last_m,last_C,FF[,,t] %>% matrix(n,r),G,D[,,t],W[,,t])

    at[,t]         <- filter$at
    Rt[,,t]        <- filter$Rt
    ft[t,]         <- filter$ft
    Qt[,,t]        <- filter$Qt
    tau[,t]        <- filter$tau
    tau_star[,t]   <- filter$tau_star
    f_star[t,]     <- filter$f_star
    Q_star[,,t]    <- filter$Q_star
    At[,,t]        <- filter$At
    mt[,t]         <- filter$mt
    Ct[,,t]        <- filter$Ct

    last_m=mt[,t]
    last_C=Ct[,,t]

    prediction=normal_pred(filter,IC_prob)

    pred[,t]      <- prediction$pred
    var.pred[,,t] <- prediction$var.pred
    icl.pred[,t]  <- prediction$icl.pred
    icu.pred[,t]  <- prediction$icu.pred
  }

  smoothed=generic_smoother(mt,Ct,at,Rt,G)
  mts <- smoothed$mts
  Cts <- smoothed$Cts

  result <- list("mt"=mt,  "Ct"=Ct,
                 "at"=at,  "At"=At, "Rt"=Rt,
                 "ft"=ft, "Qt"=Qt,
                 'f_star'=f_star, 'Q_star'=Q_star,
                 'tau'=tau,'tau_star'=tau_star,
                 "FF"=FF, "G"=G, "D"=D,"W"=W,
                 "mts"=mts, "Cts"=Cts,
                 "pred"=pred, "var.pred"=var.pred, "icl.pred"=icl.pred, "icu.pred"=icu.pred,
                 "data_out"=y)
  return(result)
}

poisson_kernel=list('fit'=poisson_fit,
                    'filter'=poisson_filter,
                    'smoother'=generic_smoother,
                    'pred'=poisson_pred,
                    'multi_var'=FALSE)

multnom_kernel=list('fit'=multnom_fit,
                    'filter'=multnom_filter,
                    'smoother'=generic_smoother,
                    'pred'=multnom_pred,
                    'multi_var'=TRUE)

normal_kernel=list('fit'=normal_fit,
                   'filter'=normal_filter,
                   'smoother'=generic_smoother,
                   'pred'=normal_pred,
                   'multi_var'=FALSE)

kernel_list=list('Poisson (univariada)'=poisson_kernel,
                 'Multinomial'=multnom_kernel,
                 'Normal'=normal_kernel)
