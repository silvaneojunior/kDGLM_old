% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/func_helper.R
\name{rowQuantile}
\alias{rowQuantile}
\title{rowQuantile}
\usage{
rowQuantile(X, q)
}
\arguments{
\item{X}{Matrix.}

\item{q}{Numeric: A number between 0 and 1.}
}
\value{
Vector: The chosen quantile for each row of X.
}
\description{
A function that calculates the row-wise quantile of a matrix.
}
\keyword{internal}
