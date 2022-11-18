#' #' system_full_gamma
#' #'
#' #' Evaluate the compatibilizing equation for the full gamma model (see Ref. Raíra).
#' #'
#' #' @param x vector: current tau values.
#' #' @param parms list: auxiliary values for the system.
#' #'
#' #' @return a vector with the values of the system (see Ref. Raíra).
#' system_full_gamma <- function(x, parms) {
#'   b=exp(x)#+(1.5/parms$Hq1)
#'
#'   a <- parms$Hq1*b
#'   n0 <- 2*a-3
#'   tau0=(n0*parms$Hq1+1)/parms$Hq2
#'
#'   # a=2.049694
#'   # b=0.717326
#'   # tau0=0.717326
#'   #
#'   # x=rgamma(20000,a,b)
#'   # Hp3=mean(lgamma(x)+x*(digamma(x)-log(tau0)))
#'   #
#'   # Hp3
#'   # f=function(x){x*(digamma(n0*x+1)-log(tau0))}
#'   # Hp3=integrate(function(x){f(x)*dgamma(x,a,b)},
#'   #               0,
#'   #               Inf)$value
#'
#'   # x=rgamma(20000,a,b)
#'   mu=a/b
#'   v=a/(b**2)
#'   alnb=mu*digamma(n0*mu+1)+v*(n0*trigamma(n0*mu+1)+n0*trigamma(n0*mu+1)+mu*(n0**2)*psigamma(n0*mu+1,2))
#'   Hp3=digamma(a)-log(b)+2*(alnb-log(tau0))
#'
#'   # x=rgamma(20000,a,b)
#'   # f=function(x){log(x)+2*x*(digamma(x)-log(tau0))}
#'   # mean(f(x))
#'
#'   f_all=c(Hp3-parms$Hq3)
#'   return(f_all)
#' }
#'
#' #' convert_FGamma_Normal
#' #'
#' #' DESCRIPTION
#' #'
#' #' @param ft vector: A vector representing the means from the normal distribution.
#' #' @param Qt matrix: A matrix representing the covariance matrix of the normal distribution.
#' #' @param parms list: A list of extra known parameters of the distribuition. Not used in this function.
#' #'
#' #' @return The parameters of the conjugated distribuition of the linear predictor.
#' #' @export
#' convert_FGamma_Normal <- function(ft, Qt, parms) {
#'
#'   f1=ft[1,]
#'   f2=ft[1,]-ft[2,]
#'   q1=Qt[1,1]
#'   q2=Qt[1,1]+Qt[2,2]+2*Qt[1,2]
#'   q12=Qt[1,1]-Qt[1,2]
#'   Hq1=exp(f1+q1/2)
#'   Hq2=exp(f2+q2/2)
#'   Hq3=f1+2*(f2+q12)*Hq1
#'
#'   parms=list(
#'     'Hq1'=Hq1,
#'     'Hq2'=Hq2,
#'     'Hq3'=Hq3)
#'
#'   ss1 <- multiroot(f = system_full_gamma, start = c(0), parms = parms, maxiter = 1000, rtol = 10**-20)
#'
#'
#'   beta=exp(as.numeric(ss1$root))#+(1.5/Hq1)
#'   alpha=Hq1*b
#'   n0 <- 2*alpha-3
#'   tau=(n0*Hq1+1)/Hq2
#'
#'   if(alpha<1.5){
#'     alpha=1.51
#'     n0 <- 2*alpha-3
#'     beta=alpha/Hq1
#'     tau=(n0*Hq1+1)/Hq2
#'   }
#'   return(list("alpha" = alpha,"beta" = beta,"tau" = tau))
#' }
#'
#' #' convert_Normal_FGamma
#' #'
#' #' Calculate the parameters of the log-Normal that best approximates the given Inverse-Gamma distribuition.
#' #' The approximation is the best in the sense that it minimizes the KL divergence from the Inverse-Gamma to the log-Normal
#' #'
#' #' @param conj_prior list: A vector containing the parameters of the Inverse-Gamma (alpha,beta).
#' #' @param parms list: A list of extra known parameters of the distribuition. Not used in this function.
#' #'
#' #' @return The parameters of the Normal distribuition of the linear predictor.
#' #' @export
#' convert_Normal_FGamma <- function(conj_prior, parms) {
#'
#'   n0=conj_prior$n
#'   tau0=conj_prior$tau
#'   theta0=conj_prior$theta
#'
#'   a=conj_prior$alpha
#'   b=conj_prior$beta
#'   tau=conj_prior$tau
#'   print(a)
#'   print(b)
#'   print(tau)
#'
#'   # a=25.689675
#'   # b=0.05121638
#'   # tau=369.50988
#'   #
#'   # a=6.594885
#'   # b=1.03109885
#'   # tau=25.08951
#'
#'   n <- 2*a-3
#'
#'
#'
#'   f1 <- digamma(a) - log(b)
#'   f1
#'
#'   sample_gamma=rgamma(20000,a,b)
#'
#'   f=function(x){digamma(n*x+1)}
#'   # f2 <- integrate(function(x){f(x)*dgamma(x,a,b)},
#'   #                 0,
#'   #                 Inf)$value
#'   f2=mean(f(sample_gamma))-log(tau)
#'   Q1 <- trigamma(a)
#'   f=function(x){(digamma(n*x+1)-log(tau)-f2)**2}
#'   # Q2 <- integrate(function(x){f(x)*dgamma(x,a,b)},
#'   #                 0,
#'   #                 Inf)$value
#'   Q2=mean(f(sample_gamma))
#'   f=function(x){(log(x)-f1)*(digamma(n*x+1)-log(tau)-f2)}
#'   # Q12 <- integrate(function(x){f(x)*dgamma(x,a,b)},
#'   #                 0,
#'   #                 Inf)$value
#'   Q12=mean(f(sample_gamma))
#'
#'   ft=matrix(c(f1,f1-f2),2,1)
#'   Qt=matrix(c(Q1,Q1-Q12,Q1-Q12,Q1+Q2+2*Q12),2,2)
#'   return(list("ft" = ft, "Qt" = Qt))
#' }
#'
#' #' update_FGamma
#' #'
#' #' Calculate posterior parameter for the Inverser-Gamma, assuming that the observed values came from a Gamma model from which the shape parameter (phi) is known and the mean (mu) have prior distribuition Inverse-Gamma.
#' #'
#' #' @param conj_prior list: A vector containing the parameters of the Inverse-Gamma (alpha,beta).
#' #' @param y vector: A vector containing the observations.
#' #' @param parms list: A list of extra known parameters of the distribuition. For this kernel, parms should containg the shape parameter (phi) for the observational gamma model.
#' #'
#' #' @return The parameters of the posterior distribution.
#' #' @export
#' update_FGamma <- function(conj_prior, y, parms) {
#'   n0=2*conj_prior$alpha-3
#'   tau0=conj_prior$tau
#'   theta0=log(tau0/n0)*n0-conj_prior$beta
#'
#'   n1=n0+1
#'   tau1=tau0+y
#'   theta1=theta0+log(y)
#'
#'   alpha=(n1+3)/2
#'   tau=tau1
#'   beta=log(tau1/n1)*n1-theta1
#'   return(list("alpha" = alpha, "beta" = beta, "tau" = tau))
#' }
#'
#' #' Fgamma_pred
#' #'
#' #' Calculate the values for the predictive distribuition given the values of the parameter of the conjugated distribuition of the linear predictor.
#' #' The data is assumed to have Gamma distribuition with known shape parameter phi and it's mean having distribuition Inverse-Gamma with shape parameter a e rate parameter b.
#' #' In this scenario, the marginal distribuition of the data is Beta prime with parameters phi, alpha, beta / phi.
#' #'
#' #' @param conj_param List or data.frame: The parameters of the conjugated distribuitions of the linear predictor.
#' #' @param outcome Vector or matrix (optional): The observed values at the current time. Not used in this function.
#' #' @param parms List: A list of extra parameters for the model. For this function, it must contain the shape parameter phi of the observational model.
#' #' @param pred_cred Numeric: the desired credibility for the credibility interval.
#' #'
#' #' @return A list containing the following values:
#' #' \itemize{
#' #'    \item pred vector/matrix: the mean of the predictive distribuition of a next observation. Same type and shape as the parameter in model.
#' #'    \item var.pred vector/matrix: the variance of the predictive distribuition of a next observation. Same type and shape as the parameter in model.
#' #'    \item icl.pred vector/matrix: the percentile of 100*((1-pred_cred)/2)% of the predictive distribuition of a next observation. Same type and shape as the parameter in model.
#' #'    \item icu.pred vector/matrix: the percentile of 100*(1-(1-pred_cred)/2)% of the predictive distribuition of a next observation. Same type and shape as the parameter in model.
#' #' }
#' #'
#' #' @importFrom extraDistr qbetapr
#' #' @export
#' #'
#' #' @examples
#' #'
#' #' conj_param <- data.frame(
#' #'   "alpha" = c(1:3),
#' #'   "beta" = c(3:1)
#' #' )
#' #'
#' #' gamma_pred(conj_param, parms = list("phi" = 1))
#' Fgamma_pred <- function(conj_param, outcome = NULL, parms = list(), pred_cred = 0.95) {
#'   alpha <- parms$alpha
#'
#'   outcome_list <- list(
#'     "pred"     = alpha,
#'     "var.pred" = alpha,
#'     "icl.pred" = alpha,
#'     "icu.pred" = alpha
#'   )
#'   return(outcome_list)
#' }
#'
#' #' Fgamma_log_like
#' #'
#' #' Calculate the predictive log-likelyhood of the data given the values of the parameter of the conjugated distribuition of the linear predictor.
#' #' The data is assumed to have Gamma distribuition with known shape parameter phi and it's mean having distribuition Inverse-Gamma with shape parameter a e rate parameter b.
#' #' In this scenario, the marginal distribuition of the data is Beta prime with parameters phi, alpha, beta / phi.
#' #'
#' #' @param conj_param List or data.frame: The parameters of the conjugated distribuition (Gamma) of the linear predictor.
#' #' @param outcome Vector or matrix: The observed values at the current time.
#' #' @param parms List: A list of extra parameters for the model. For this function, it must contain the shape parameter phi of the observational model.
#' #'
#' #' @return The predictive log-likelyhood of the observed data.
#' #' @export
#' #'
#' #' @importFrom extraDistr dbetapr
#' #' @examples
#' #'
#' #' conj_param <- data.frame(
#' #'   "alpha" = c(1:3),
#' #'   "beta" = c(3:1)
#' #' )
#' #'
#' #' gamma_log_like(conj_param, rgamma(3, c(1:3), c(3:1)), parms = list("phi" = 1))
#' Fgamma_log_like <- function(conj_param, outcome, parms = list()) {
#'   return(outcome*0)
#' }
#'
#' #' @export
#' Fgamma_kernel <- list(
#'   "conj_prior" = convert_FGamma_Normal,
#'   "conj_post" = convert_Normal_FGamma,
#'   "update" = update_FGamma,
#'   "smoother" = generic_smoother,
#'   "pred" = Fgamma_pred,
#'   "log.like" = Fgamma_log_like,
#'   "offset" = log_offset_half,
#'   "link_function" = log,
#'   "inv_link_function" = exp,
#'   "param_names" = function(y) {
#'     c("alpha", "beta",'tau')
#'   },
#'   "multi_var" = FALSE
#' )
