\name{runEMclonalCN}
\alias{runEMclonalCN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Function to run the Expectation Maximization Algorithm in \pkg{TitanCNA}.
}
\description{
  Function to run the Expectation Maximization Algorithm for inference of model parameters: cellular prevalence, normal proportion, tumour ploidy. This is a key function in the \pkg{TitanCNA} package and is the most computationally intense. This function makes calls to a C subroutine that allows the algorithm to be run more efficiently.
}
\usage{
runEMclonalCN(data, gParams, nParams, pParams, sParams,
              txnExpLen = 1e15, txnZstrength = 5e05, maxiter = 15,
              maxiterUpdate = 1500, pseudoCounts = 1e-300,
              normalEstimateMethod = "map", estimateS = TRUE, 
              estimatePloidy = TRUE, useOutlierState = FALSE, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
    \code{\link[base:list]{list}} object that contains the components for the data to be analyzed. \code{chr}, \code{posn}, \code{ref}, and \code{tumDepth} that can be obtained using \code{\link{loadAlleleCounts}}, and \code{logR} that can be obtained using \code{\link{correctReadDepth}} and \code{\link{getPositionOverlap}} (see Example).
  }
  \item{gParams}{
    \code{\link[base:list]{list}} object that contains the copy number and allelic ratio genotype parameters. Can be obtained from \code{\link{loadDefaultParameters}}.
  }
  \item{nParams}{
    \code{\link[base:list]{list}} object that contains the normal contamination parameters. Can be obtained from \code{\link{loadDefaultParameters}}.
  }
  \item{pParams}{
    \code{\link[base:list]{list}} object that contains the tumour ploidy parameters. Can be obtained from \code{\link{loadDefaultParameters}}.
  }
  \item{sParams}{
    \code{\link[base:list]{list}} object that contains the subclonality (cellular prevalence and clonal cluster) parameters. Can be obtained from \code{\link{loadDefaultParameters}}.
  }
  \item{txnExpLen}{
    Influences prior probability of genotype transitions in the HMM. Smaller value have lower tendency to change state; however, too small and it produces underflow problems. \code{1e-9} works well for up to 3 million total positions.
  }
  \item{txnZstrength}{
    Influences prior probability of clonal cluster transitions in the HMM. Smaller value have lower tendency to change clonal cluster state. \code{1e-9} works well for up to 3 million total positions.
  }
  \item{pseudoCounts}{
    Small, machine precision values to add to probabilities to avoid underflow. For example, \code{.Machine$double.eps}.
  }
  \item{maxiter}{
    Maximum number of expectation-maximization iterations allowed. In practice, for \pkg{TitanCNA}, it will usually not exceed 20.
  }
  \item{maxiterUpdate}{
    Maximum number of coordinate descent iterations during the M-step (of EM algorithm) when parameters are estimated.
  }
  \item{normalEstimateMethod}{
    Specifies how to handle normal proportion estimation. Using \code{map} will use the maximum a posteriori estimation. Using \code{fixed} will not estimate the normal proportion; the normal proportion will be fixed to whatever is specified in \code{params$normalParams$n_0}. See Details.
  }
  \item{estimateS}{
    Logical indicating whether to account for clonality and estimate subclonal events.  See Details.
  }
  \item{estimatePloidy}{
    Logical indicating whether to estimate and account for tumour ploidy.
  }
  \item{useOutlierState}{
    Logical indicating whether an additional outlier state should be used. In practice, this is usually not necessary.
  }
  \item{verbose}{
    Set to FALSE to suppress program messages.
  }
}

\details{
  This function is implemented with the \code{"\link{foreach}"} package and therefore supports parallelization. See \code{"doMC"} or \code{"doMPI"} for some parallelization packages.
  
  The forwards-backwards algorithm is used for the E-step in the EM algorithm. This is done using a call to a C subroutine for each chromosome. The maximization step uses maximum a posteriori (MAP) for estimation of parameters.

  If the sample has absolutely no normal contamination, then assign \code{nParams$n_0 <- 0} and use argument \code{normalEstimateMethod="fixed"}.
  
  \code{estimateS} should always be set to \code{TRUE}. If no subclonality is expected, then use \code{\link{loadDefaultParameters}(numberClonalClusters=1)}. Using \code{estimateS=FALSE} and \code{\link{loadDefaultParameters}(numberClonalClusters=0)} is gives more or less the same results.
}

\value{
  \code{\link[base:list]{list}} with components for results returned from the EM algorithm, including converged parameters, posterior marginal responsibilities, log likelihood, and original parameter settings.
  \item{n}{Converged estimate for normal contamination parameter. \code{numeric array} containing estimates at each EM iteration.}
  \item{s}{Converged estimate(s) for cellular prevalence parameter(s). This value is defined as the proportion of tumour sample that does \emph{not} contain the aberrant genotype. This will contrast what is output in \code{\link{outputTitanResults}}. \code{numeric array} containing estimates at each EM iteration. If more than one cluster is specified, then \code{s} is a \code{numeric matrix}.}
  \item{var}{Converged estimates for variance parameter of the Gaussian mixtures used to model the log ratio data. \code{numeric matrix} containing estimates at each EM iteration.}
  \item{phi}{Converged estimate for tumour ploidy parameter. \code{numeric array} containing estimates at each EM iteration.}
  \item{piG}{Converged estimate for initial genotype state distribution. \code{numeric matrix} containing estimates at each EM iteration.}
  \item{piZ}{Converged estimate for initial clonal cluster state distribution. \code{numeric matrix} containing estimates at each EM iteration.}
  \item{muR}{Mean of binomial mixtures computed as a function of \code{s} and \code{n}. \code{numeric matrix} containing estimates at each EM iteration. See References for mathematical details.}
  \item{muC}{Mean of Gaussian mixtures computed as a function of \code{s}, \code{n}, and \code{phi}. \code{numeric matrix} containing estimates at each EM iteration. See References for mathematical details.}
  \item{loglik}{Posterior Log-likelihood that includes data likelihood and the priors. \code{numeric array} containing estimates at each EM iteration.}
  \item{rhoG}{Posterior marginal probabilities for the genotype states computed during the E-step. Only the final iteration is returned as a \code{numeric matrix}.}
  \item{rhoZ}{Posterior marginal probabilities for the clonal cluster states computed during the E-step. Only the final iteration is returned as a \code{numeric matrix}.}
  \item{genotypeParams}{Original genotype parameters. See \code{\link{loadDefaultParameters}}.}
  \item{ploidyParams}{Original tumour ploidy parameters. See \code{\link{loadDefaultParameters}}.}
  \item{normalParams}{Original normal contamination parameters. See \code{\link{loadDefaultParameters}}.}
  \item{clonalParams}{Original subclonal parameters. See \code{\link{loadDefaultParameters}}.}
  \item{txnExpLen}{Original genotype transition expected length. See \code{\link{loadDefaultParameters}}.}
  \item{txnZstrength}{Original clonal cluster transition expected length. See \code{\link{loadDefaultParameters}}.}
  \item{useOutlierState}{Original setting indicating usage of outlier state. See \code{\link{loadDefaultParameters}}.}
  
}

\references{
Ha, G., Roth, A., Khattra, J., Ho, J., Yap, D., Prentice, L. M., Melnyk, N., McPherson, A., Bashashati, A., Laks, E., Biele, J., Ding, J., Le, A., Rosner, J., Shumansky, K., Marra, M. A., Huntsman, D. G., McAlpine, J. N., Aparicio, S. A. J. R., and Shah, S. P. (2014). TITAN: Inference of copy number architectures in clonal cell populations from tumour whole genome sequence data. Genome Research, Published in advance July 24 2014. (PMID: 25060187)
}

\author{
  Gavin Ha <gavinha@gmail.com>
}

\seealso{
  \code{"\link{foreach}"}, \code{"doMC"}, \code{"doMPI"}, \code{\link{loadAlleleCounts}}, \code{\link{loadDefaultParameters}}, \code{\link{viterbiClonalCN}}
}

\examples{
message('Running TITAN ...')
#### LOAD DATA ####
infile <- system.file("extdata", "test_alleleCounts_chr2.txt", 
              package = "TitanCNA")
data <- loadAlleleCounts(infile)

#### LOAD PARAMETERS ####
message('titan: Loading default parameters')
numClusters <- 2
params <- loadDefaultParameters(copyNumber = 5, 
              numberClonalClusters = numClusters, skew = 0.1)

#### READ COPY NUMBER FROM HMMCOPY FILE ####
message('titan: Correcting GC content and mappability biases...')
tumWig <- system.file("extdata", "test_tum_chr2.wig", package = "TitanCNA")
normWig <- system.file("extdata", "test_norm_chr2.wig", package = "TitanCNA")
gc <- system.file("extdata", "gc_chr2.wig", package = "TitanCNA")
map <- system.file("extdata", "map_chr2.wig", package = "TitanCNA")
cnData <- correctReadDepth(tumWig, normWig, gc, map)
logR <- getPositionOverlap(data$chr, data$posn, cnData)
data$logR <- log(2^logR) #transform to natural log

#### FILTER DATA FOR DEPTH, MAPPABILITY, NA, etc ####
data <- filterData(data, 1:24, minDepth = 10, maxDepth = 200, map = NULL)

#### EM (FWD-BACK) TO TRAIN PARAMETERS ####
#### Can use parallelization packages ####
K <- length(params$genotypeParams$alphaKHyper)
params$genotypeParams$alphaKHyper <- rep(500, K)
params$ploidyParams$phi_0 <- 1.5
convergeParams <- runEMclonalCN(data, gParams = params$genotypeParams, 
                                nParams = params$normalParams, 
                                pParams = params$ploidyParams, 
                                sParams = params$cellPrevParams, 
                                maxiter = 3, maxiterUpdate = 500, 
                                txnExpLen = 1e9, txnZstrength = 1e9, 
                                useOutlierState = FALSE, 
                                normalEstimateMethod = "map", 
                                estimateS = TRUE, estimatePloidy = TRUE)
}

\keyword{manip}
\keyword{iteration}
\keyword{htest}
\keyword{models}
\keyword{multivariate}
