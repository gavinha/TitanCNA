\name{viterbiClonalCN}
\alias{viterbiClonalCN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Function to run the Viterbi algorithm for \pkg{TitanCNA}.
}
\description{
  Function to run the Viterbi algorithm to find the optimal state path in the \pkg{TitanCNA} hidden Markov model (HMM). The states returned will indicate the optimal copy number and LOH state as well as the most likely clonal cluster for each data point. After running EM, use the converge parameters and the input data to infer the optimal state for each position. This function makes calls to a C subroutine that allows the algorithm to be run more efficiently.
}
\usage{
  viterbiClonalCN(data, convergeParams, genotypeParams = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
    \code{\link[base:list]{list}} object that contains the components for the data to be analyzed. \code{chr}, \code{posn}, \code{ref}, and \code{tumDepth} that can be obtained using \code{\link{loadAlleleCounts}}, and \code{logR} that can be obtained using \code{\link{correctReadDepth}} and \code{\link{getPositionOverlap}} (see Example).
  }
  \item{convergeParams}{
    \code{\link[base:list]{list}} object that is returned from the function \code{\link{runEMclonalCN}} in \pkg{TitanCNA}.
  }
  \item{genotypeParams}{
    If \code{convergeParams} does not contain a \code{genotypeParams} element, then the user can pass this as an argument.
  }
}
\details{
  It is difficult to interpret the output of this function directly. The user should use the function \code{\link{outputTitanResults}} after.
}

\value{
  \code{numeric \link[base:array]{array}} containing the integer states corresponding to each data point in \code{data}.
}

\references{
  G. Ha, A. Roth, N. Melnyk, A. McPherson, L. M. Prentice, A. Bashashati, J. Ding, A. Le, J. Rosner, K. Shumansky, M. A. Marra, D. G. Huntsman, J. N. McAlpine, S. A. J. R. Aparicio, S. P. Shah. Modeling the evolution of copy number architecture in tumour whole genome sequencing data. Under Review 2013.
}

\author{
  Gavin Ha <gavinha@gmail.com>
}

\seealso{
  \code{\link{outputTitanResults}}, \code{\link{loadAlleleCounts}}
}

\examples{
data(EMresults)

#### COMPUTE OPTIMAL STATE PATH USING VITERBI ####
optimalPath <- viterbiClonalCN(data, convergeParams)
}

\keyword{manip}
\keyword{iteration}
\keyword{htest}
\keyword{models}
\keyword{multivariate}

