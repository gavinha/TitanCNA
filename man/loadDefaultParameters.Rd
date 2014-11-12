\name{loadDefaultParameters}
\alias{loadDefaultParameters}

\title{
  Load TITAN parameters
}

\description{
  Load TITAN model parameters based on maximum copy number and number of clonal clusters.
}

\usage{
  loadDefaultParameters(copyNumber = 5, numberClonalClusters = 1, skew = 0,
      symmetric = TRUE)
}

\arguments{
  \item{copyNumber}{
    Maximum number of absolute copies to account for in the model. Default (and recommended) is 5.
  }
  \item{numberClonalClusters}{
    Number of clonal clusters to use. Using \sQuote{1} represents more tumour subclonality. \sQuote{2} or higher treats the tumour data as being subclonal. 
  }
  \item{skew}{
    \code{numeric} float indicating the heterozygous baseline shift for the allelic ratios towards 1. This is may be required for SOLiD data, but for most cases, this argument can be omitted.
  }
  \item{symmetric}{
  	\code{logical}; if \code{TRUE}, then treat genotypes as symmetric. See Details.
  }
}

\details{
  Generally, \pkg{TitanCNA} should be run once for each number of clonal clusters in the range of 1 to 5. Then, use model selection to choose the run with the optimal number of clusters.
  
  If the allelic ratio data is skewed towards one allele, then use \code{skew} to help define the baseline. For example, if the data is skewed towards the reference, then use 0.1 so that the heterozygous baseline is at 0.6. The allelic ratio baseline is normally at 0.5.
  
  \code{sParams}, which represents the parameters for estimation of subclonality, always contains values for one cluster that represents the clonally dominant cluster (events present in nearly all tumour cells) with an initial value of \code{sParams$s_0[1] = 0.001}.
  
  Setting \code{symmetric} to \code{TRUE} will treat reference and non-reference alleles the same. For example, genotypes \code{AA} (homozygous for reference allele) and \code{BB} (homozygous for non-reference allele) as being equivalent. This will reduce the state space substantially. 
}

\value{
  \code{\link[base:list]{list}} containing 4 sets of parameters, each as a component: 
  \item{genotypeParams}{Parameters for copy number and allelic ratios geneotype states}
  \item{normalParams}{Parameters for normal contamination}
  \item{ploidyParams}{Parameters for average tumour ploidy}
  \item{sParams}{Parameters for modeling subclonality: clonal clusters and cellular prevalence}
}

\references{
Ha, G., Roth, A., Khattra, J., Ho, J., Yap, D., Prentice, L. M., Melnyk, N., McPherson, A., Bashashati, A., Laks, E., Biele, J., Ding, J., Le, A., Rosner, J., Shumansky, K., Marra, M. A., Huntsman, D. G., McAlpine, J. N., Aparicio, S. A. J. R., and Shah, S. P. (2014). TITAN: Inference of copy number architectures in clonal cell populations from tumour whole genome sequence data. Genome Research, 24: 1881-1893. (PMID: 25060187)
}
\author{
  Gavin Ha <gavinha@gmail.com>
}

\seealso{
  \code{\link{loadAlleleCounts}}
}

\examples{
  #### LOAD PARAMETERS ####
  numClusters <- 2
  params <- loadDefaultParameters(copyNumber = 5, 
                                  numberClonalClusters = numClusters)
}


\keyword{attribute}
