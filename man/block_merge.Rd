% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/structure_helper.R
\name{block_merge}
\alias{block_merge}
\title{Auxiliary function to merge blocks}
\usage{
block_merge(...)
}
\arguments{
\item{...}{dlm_block: A sequence of block to be merged.}
}
\value{
The merged block as a dlm_block.
}
\description{
An auxiliary function to merge blocks.
}
\examples{

# Long way
level_1 <- polynomial_block(alpha1 = 1, order = 1)
level_2 <- polynomial_block(alpha2 = 1, order = 2)
season_2 <- harmonic_block(alpha2 = 1, period = 20)

final_block <- block_merge(level_1, level_2, season_2)

# Short way
final_block <- polynomial_block(alpha1 = 1, order = 1) +
  polynomial_block(alpha2 = 1, order = 2) +
  harmonic_block(alpha2 = 1, period = 20)

}
\seealso{
Other {auxiliary functions for structural blocks}: 
\code{\link{AR_block}()},
\code{\link{block_mult}()},
\code{\link{block_rename}()},
\code{\link{harmonic_block}()},
\code{\link{polynomial_block}()},
\code{\link{report_block}()}
}
\concept{{auxiliary functions for structural blocks}}
