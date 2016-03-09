\name{removeEmptyClusters}
\alias{removeEmptyClusters}

\title{
  Post-process TitanCNA results by removing clusters with proportion of data points altered lower than a threshold.  The number of clonal clusters, cellular prevalence, and normal contamination will be adjusted to reflect the remaining clonal clusters.   
}
\description{
  Filters all vectors in list based on specified chromosome(s) of interest, minimum and maximum read depths, missing data, mappability score threshold
}

\usage{
removeEmptyClusters(convergeParams, results, proportionThreshold = 0.001, proportionThresholdClonal = 0.3)
}


\arguments{
  \item{convergeParams}{
    \code{\link[base:list]{list}} object that is returned from the function \code{\link{runEMclonalCN}} in \pkg{TitanCNA}.      
  }
  \item{results}{
   \code{\link[base:data.frame]{data.frame}} returned by \code{\link{outputTitanResults}}. Each row corresponds to a genomic SNP position in the analysis.
  }
  \item{proportionThreshold}{
  	Minimum proportion of the genome altered (by SNPs) for a cluster to be retained.  Clonal clusters having lower proportion of alteration are removed. 
  }
  \item{proportionThresholdClonal}{
  	Minimum proportion of genome altered by clonal events (by SNPs) for the highest cellular prevalence cluster.
  }
}

\details{
  All vectors in the input \code{data} list object, and \code{map}, must all have the same number of rows.
}

\value{
	\code{\link[base:list]{list}} with components:
	\item{convergeParams}{The same data structure and format as the input \code{convergeParams} with removed clusters and adjusted parameters.}
	\item{results}{The same data structure and format as the input \code{results} with removed clusters and adjusted clonal cluster and cellular prevalence values.}
}

\references{
Ha, G., Roth, A., Khattra, J., Ho, J., Yap, D., Prentice, L. M., Melnyk, N., McPherson, A., Bashashati, A., Laks, E., Biele, J., Ding, J., Le, A., Rosner, J., Shumansky, K., Marra, M. A., Huntsman, D. G., McAlpine, J. N., Aparicio, S. A. J. R., and Shah, S. P. (2014). TITAN: Inference of copy number architectures in clonal cell populations from tumour whole genome sequence data. Genome Research, 24: 1881-1893. (PMID: 25060187)
}

\author{
  Gavin Ha <gavinha@gmail.com>
}

\seealso{
  \code{\link{outputTitanResults}}
}

\examples{
data(EMresults)

#### COMPUTE OPTIMAL STATE PATH USING VITERBI ####
optimalPath <- viterbiClonalCN(data, convergeParams)

#### FORMAT RESULTS ####
results <- outputTitanResults(data, convergeParams, optimalPath, 
                              filename = NULL, posteriorProbs = FALSE,
                              subcloneProfiles = TRUE)

#### REMOVE EMPTY CLONAL CLUSTERS ####
corrResults <- removeEmptyClusters(convergeParams, results, proportionThreshold = 0.001,
																		proportionThresholdClonal = 0.3)
convergeParams <- corrResults$convergeParams
results <- corrResults$results

}


\keyword{manip}
