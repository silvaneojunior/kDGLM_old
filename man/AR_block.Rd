% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/structure_helper.R
\name{AR_block}
\alias{AR_block}
\title{Structural blocks for auto regressive trends and regressions}
\usage{
AR_block(
  ...,
  order,
  noise_var,
  pulse = 0,
  name = "Var_AR",
  AR_support = "free",
  D = 1,
  W = 0,
  m0 = c(1, rep(0, order - 1)),
  C0 = 1,
  drift = 0,
  m0_states = 0,
  C0_states = diag(order),
  D_states = 1,
  drift_states = 0,
  m0_pulse = 0,
  C0_pulse = 1,
  D_pulse = 1,
  W_pulse = 0,
  drift_pulse = 0
)
}
\arguments{
\item{...}{Named values for the planing matrix.}

\item{order}{Positive integer: The order of the AR block.}

\item{noise_var}{Non negative scalar: The variance of the white noise added to the latent state.}

\item{pulse}{Vector or scalar: An optional argument providing the values for the pulse for a Transfer Function. Default is 0.}

\item{name}{String: An optional argument providing the name for this block. Can be useful to identify the models with meaningful labels, also, the name used will be used in some auxiliary functions.}

\item{AR_support}{String: Either "constrained" or "free". If AR_support is "constrained", then the AR coefficients will be forced to be on the interval (-1,1), otherwise, the coefficients will be unrestricted. Beware that, under no restriction on the coefficients, there is no guarantee that the estimated coefficients will imply in a stationary process, furthermore, if the order of the AR block is greater than 1, then the restriction imposed when AR_support is equal to "constrained" does NOT guarantee that the process will be stationary (although it may help).}

\item{D}{Array, Matrix, vector or  scalar: The values for the discount factors associated with the AR coefficients at each time. If D is a array, it's dimensions should be n x n x t, where n is the order of the AR block and t is the length of the outcomes. If D is a matrix, it's dimensions should be n x n and it's values will be used for each time. If D is a vector or scalar, a discount factor matrix will be created as a diagonal matrix with the values of D in the diagonal.}

\item{W}{Array, Matrix, vector or  scalar: The values for the covariance matrix for the noise factor associated with the AR coefficients at each time. If W is a array, it's dimensions should be n x n x t, where n is the order of the AR block and t is the length of the outcomes. If W is a matrix, it's dimensions should be n x n and it's values will be used for each time. If W is a vector or scalar, a discount factor matrix will be created as a diagonal matrix with the values of W in the diagonal.}

\item{m0}{Vector or scalar: The prior mean for the coefficients associated with this block. If m0 is a vector, it's dimension should be equal to the order of the AR block. If m0 is a scalar, it's value will be used for all coefficients. If the coefficients are restricted to the interval (-1,1), the m0 is interpreted as the mean for logit((rho+1)/2), where rho is the AR coefficient.}

\item{C0}{Matrix, vector or scalar: The prior covariance matrix for the coefficients associated with this block. If C0 is a matrix, it's dimensions should be n x n, where n is the order of the AR block. If C0 is a vector or scalar, a covariance matrix will be created as a diagonal matrix with the values of C0 in the diagonal. If the coefficients are restricted to the interval (-1,1), the C0 is interpreted as the covariance matrix for logit((rho+1)/2), where rho is the AR coefficient.}

\item{drift}{Matrix, vector or scalar: A drift to be add in the AR coefficients after the temporal evoltion (can be interpreted as the mean of the random noise at each time). If a matrix, it's dimension should be n x T, where n is the order of the AR block and T is the length of the serie. If a scalar, the passed value will be used for all coefficients at each time.}

\item{m0_states}{Vector or scalar: The prior mean for the states associated with this block. If m0_states is a vector, it's dimension should be equal to the order of the AR block. If m0_states is a scalar, it's value will be used for all coefficients.}

\item{C0_states}{Matrix, vector or scalar: The prior covariance matrix for the states associated with this block. If C0_states is a matrix, it's dimensions should be n x n, where n is the order of the AR block. If C0_state is a vector or scalar, a covariance matrix will be created as a diagonal matrix with the values of C0_state in the diagonal.}

\item{D_states}{Array, Matrix, vector or  scalar: The values for the discount factors for the states associated with this block. If D_states is a array, it's dimensions should be n x n x t, where n is the order of the AR block and t is the length of the outcomes. If D_states is a matrix, it's dimensions should be n x n and it's values will be used for each time. If D_states is a vector or scalar, a discount factor matrix will be created as a diagonal matrix with the values of D_states in the diagonal.}

\item{drift_states}{Matrix, vector or scalar: A drift to be add in the states after the temporal evoltion (can be interpreted as the mean of the random noise at each time). If a matrix, it's dimension should be n x T, where n is the order of the AR block and T is the length of the serie. If a vector, it should have size T, and each value will be applied to the first latent variable (the one which affects the linear predictors) in their respective time. If a scalar, the passed value will be used for the first latent variable at each time.}

\item{m0_pulse}{Vector or scalar: The prior mean for the coefficients associated with the pulses. If m0_pulse is a vector, it's dimension should be equal to the number of pulses. If m0_pulse is a scalar, it's value will be used for all coefficients.}

\item{C0_pulse}{Matrix, vector or scalar: The prior covariance matrix for the coefficients associated with the pulses. If C0_pulse is a matrix, it's dimensions should be n x n, where n is the number of pulses. If C0_pulse is a vector or scalar, a covariance matrix will be created as a diagonal matrix with the values of C0_pulse in the diagonal.}

\item{D_pulse}{Array, Matrix, vector or  scalar: The values for the discount factors associated with pulse coefficients at each time. If D_pulse is a array, it's dimensions should be n x n x t, where n is the number of pulses and t is the length of the outcomes. If D_pulse is a matrix, it's dimensions should be n x n and it's values will be used for each time. If D_pulse is a vector or scalar, a discount factor matrix will be created as a diagonal matrix with the values of D_pulse in the diagonal.}

\item{W_pulse}{Array, Matrix, vector or  scalar: The values for the covariance matrix for the noise factor associated with pulse coefficients at each time. If W_pulse is a array, it's dimensions should be n x n x t, where n is the number of pulses and t is the length of the outcomes. If W_pulse is a matrix, it's dimensions should be n x n and it's values will be used for each time. If W_pulse is a vector or scalar, a covariance matrix will be created as a diagonal matrix with the values of W_pulse in the diagonal.}

\item{drift_pulse}{Matrix, vector or scalar: A drift to be add in the pulse effect after the temporal evoltion (can be interpreted as the mean of the random noise at each time). If a matrix, it's dimension should be n x T, where n is the number of pulses and T is the length of the serie. If a scalar, the passed value will be used for all latent variable at each time.}
}
\value{
An object of the class dlm_block containing the following values:
\itemize{
   \item FF Array: A 3D-array containing the regression matrix for each time. It's dimension should be n x m x T, where n is the number of latent variables, m is the number of outcomes in the model and T is the time series length.
   \item FF_labs Array: DESCRIPTION
   \item G Matrix: The state evolution matrix.
   \item D Array: A 3D-array containing the discount factor matrix for each time. It's dimension should be n x n x T, where n is the number of latent variables and T is the time series length.
   \item W Array: A 3D-array containing the covariance matrix of the noise for each time. It's dimension should be the same as D.
   \item m0 Vector: The prior mean for the latent vector.
   \item C0 Matrix: The prior covariance matrix for the latent vector.
   \item var_names list: A list containing the variables indexes by their name.
   \item order Positive integer: Same as argument.
   \item n Positive integer: The number of latent variables associated with this block (2).
   \item t Positive integer: The number of time steps associated with this block. If 1, the block is compatible with blocks of any time length, but if t is greater than 1, this block can only be used with blocks of the same time length.
   \item k Positive integer: The number of outcomes associated with this block. This block can only be used with blocks with the same outcome length.
   \item pred_names Vector: The name of the linear predictors associated with this block,
   \item type Character: The type of block (AR).
   \item AR_support Character: Same as argument.
}
}
\description{
Creates the structure for a Auto Regressive (AR) block (see West & Harrison (1997), chapter 9) with desired order.
As the package suppose that the structure of the model is linear, a linearization is applied to the evolution equation, as described in West & Harrison (1997), chapter 13.
This block also supports Transfer Functions, being necessary to specify the associated pulse when calling the AR_block function (see arg.).
}
\details{
For the ..., noise_var, D, W, m0, C0, m0_states, C0_states, D_states, m0_pulse, C0_pulse, D_pulse, W_pulse arguments, the user may set one or more of it's values as a string.
By doing so, the user will leave the block partially undefined and it can no longer be used in the \code{\link{fit_model}} function.
Instead, the user must use the \code{\link{search_model}} function to search the best hyper parameters among a defined range of possible values.
See the \code{\link{search_model}} function for details on it's usage.

For the details about the implementation see \insertCite{ArtigoPacote;textual}{kDGLM}.

For the details about Auto regressive models in the context of DLM's, see \insertCite{WestHarr-DLM;textual}{kDGLM}, chapter 9.

For the details about the linearization of non-linear evolution equations in the context of DLM's, see \insertCite{WestHarr-DLM;textual}{kDGLM}, chapter 13.

For the details about dynamic regression models in the context of DLM's, see \insertCite{WestHarr-DLM;textual}{kDGLM}, chapters 6 and 9.
}
\examples{

AR_block(mu = 1, pulse = rnorm(200), order = 3, noise_var = 0.1)

}
\references{
\insertAllCited{}
}
\seealso{
\code{\link{fit_model}}

Other {auxiliary functions for structural blocks}: 
\code{\link{block_merge}()},
\code{\link{block_mult}()},
\code{\link{block_rename}()},
\code{\link{harmonic_block}()},
\code{\link{polynomial_block}()},
\code{\link{report_block}()}
}
\concept{{auxiliary functions for structural blocks}}
