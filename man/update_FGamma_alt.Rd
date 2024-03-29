% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kernel_gamma.R
\name{update_FGamma_alt}
\alias{update_FGamma_alt}
\title{update_FGamma_alt}
\usage{
update_FGamma_alt(conj_prior, ft, Qt, y, parms)
}
\arguments{
\item{conj_prior}{list: A vector containing the parameters of the Inverse-Gamma (alpha,beta). Not used in the alternative method.}

\item{ft}{vector: A vector representing the means from the normal distribution.}

\item{Qt}{matrix: A matrix representing the covariance matrix of the normal distribution.}

\item{y}{vector: A vector containing the observations.}

\item{parms}{list: A list of extra known parameters of the distribution. For this kernel, parms should containg the shape parameter (phi) for the observational gamma model.}
}
\value{
The parameters of the posterior distribution.
}
\description{
Calculate the (approximated) posterior parameter for the linear predictors, assuming that the observed values came from a Gamma model from which the shape and mean parameters have prior distribution in the log-Normal family.
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
Other {auxiliary functions for a Gamma outcome with unknowned shape}: 
\code{\link{Fgamma_pred_alt}()},
\code{\link{Fgamma_pred}()},
\code{\link{convert_FGamma_Normal}()},
\code{\link{convert_Normal_FGamma}()},
\code{\link{update_FGamma}()}
}
\concept{{auxiliary functions for a Gamma outcome with unknowned shape}}
\keyword{internal}
