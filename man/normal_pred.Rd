% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kernel_normal.R
\name{normal_pred}
\alias{normal_pred}
\title{normal_pred}
\usage{
normal_pred(conj_param, outcome = NULL, parms = list(), pred_cred = 0.95)
}
\arguments{
\item{conj_param}{List or data.frame: The parameters of the conjugated distributions of the linear predictor.}

\item{outcome}{Vector or matrix (optional): The observed values at the current time. Not used in this function.}

\item{parms}{List: A list of extra parameters for the model. For this function, it must contain the observational covariance matrix, Sigma.}

\item{pred_cred}{Numeric: the desired credibility for the credibility interval.}
}
\value{
A list containing the following values:
\itemize{
   \item pred vector/matrix: the mean of the predictive distribution of a next observation. Same type and shape as the parameter in model.
   \item var.pred vector/matrix: the variance of the predictive distribution of a next observation. Same type and shape as the parameter in model.
   \item icl.pred vector/matrix: the percentile of 100*((1-pred_cred)/2)% of the predictive distribution of a next observation. Same type and shape as the parameter in model.
   \item icu.pred vector/matrix: the percentile of 100*(1-(1-pred_cred)/2)% of the predictive distribution of a next observation. Same type and shape as the parameter in model.
   \item log.like vector: the The log likelihood for the outcome given the conjugated parameters.
}
}
\description{
Calculate the values for the predictive distribution given the values of the parameter of the conjugated distribution of the linear predictor.
The data is assumed to have Normal distribution with known variance and it's mean having distribution Normal.
In this scenario, the marginal distribution of the data is also Normal.
}
\seealso{
Other {auxiliary functions for a Normal outcome}: 
\code{\link{convert_multi_NG_Normal}()},
\code{\link{multi_normal_gamma_pred}()},
\code{\link{update_NG_alt}()},
\code{\link{update_Normal}()},
\code{\link{update_multi_NG_correl}()}
}
\concept{{auxiliary functions for a Normal outcome}}
\keyword{internal}
