% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kernel_poisson.R
\name{update_Poisson_alt}
\alias{update_Poisson_alt}
\title{update_Poisson_alt}
\usage{
update_Poisson_alt(conj_prior, ft, Qt, y, parms)
}
\arguments{
\item{conj_prior}{list: A vector containing the parameters of the Gamma (alpha,beta). Not used in the alternative method.}

\item{ft}{vector: A vector representing the means from the normal distribution.}

\item{Qt}{matrix: A matrix representing the covariance matrix of the normal distribution.}

\item{y}{vector: A vector containing the observations.}

\item{parms}{list: A list of extra known parameters of the distribution. Not used in this kernel.}
}
\value{
The parameters of the posterior distribution.
}
\description{
Calculate the (approximated) posterior parameter for the linear predictors, assuming that the observed values came from a Poisson model from which the rate parameter have prior distribution in the log-Normal family.
}
\details{
For evaluating the posterior parameters, we use a modified version of the method proposed in \insertCite{ArtigokParametrico;textual}{kDGLM}.

For computational efficiency, we also use a Laplace approximations to obtain the first and second moments of the posterior \insertCite{@see @TierneyKadane1 and @TierneyKadane2 }{kDGLM}.

For the details about the implementation see  \insertCite{ArtigoPacote;textual}{kDGLM}.

For the detail about the modification of the method proposed in \insertCite{ArtigokParametrico;textual}{kDGLM}, see \insertCite{ArtigoAltMethod;textual}{kDGLM}.
}
\references{
\insertAllCited{}
}
\seealso{
Other {auxiliary functions for a Poisson outcome}: 
\code{\link{convert_Normal_Poisson}()},
\code{\link{convert_Poisson_Normal}()},
\code{\link{poisson_pred}()},
\code{\link{update_Poisson}()}
}
\concept{{auxiliary functions for a Poisson outcome}}
\keyword{internal}
