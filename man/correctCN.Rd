\name{correctIntegerCopyNumber}
\alias{correctIntegerCN}

\title{
  Compute purity and ploidy corrected log ratios; recompute integer CN for high-level amplifications.
}
\description{
  TitanCNA uses a finite state space that defines a maximum number of copies to model. High-level amplifications that exceed this defined maximum need to be corrected and reported as the likely copy number based on the observed data. \code{correctIntegerCN} performs two tasks: (1) correct log ratio based on purity and ploidy, and then convert to decimal CN value; (2) Correct bins (from \code{cn}) and segments (from \code{segs}) in which the original predicted integer copy number was assigned the maximum CN state; bins and segments for all of chromosome X are also corrected, if provided in the input.
}

\usage{
correctIntegerCN(cn, segs, purity, ploidy, maxCNtoCorrect.autosomes = NULL, 
		maxCNtoCorrect.X = NULL, correctHOMD = TRUE, minPurityToCorrect = 0.2, gender = "male", 
		chrs = c(1:22, "X"))
}


\arguments{
  \item{cn}{
    \link[data.table]{data.table} object output from the function \link[TitanCNA]{outputTitanResults}
  }
  \item{segs}{
    \link[data.table]{data.table} object output from the function \link[TitanCNA]{outputTitanSegments}
  }
  \item{purity}{
    Float type of the 1 minus the normal contamination estimate from TitanCNA
  }
  \item{ploidy}{
    Float type of the average tumor ploidy estimate from TitanCNA
  }
  \item{maxCNtoCorrect.autosomes}{
    Bins and segments in autosomes with this copy number value or higher will be corrected. If \code{NULL}, then it will use the original copy number value from the input data.
  }
  \item{maxCNtoCorrect.X}{
    Bins and segments in chromosome X, if provided, with this copy number value or higher will be corrected. If \code{NULL}, then it will use the original copy number value from the input data.
  }
  \item{minPurityToCorrect}{
    If \code{purity} is less than \code{minPurityToCorrect}, then \code{Corrected_Copy_Number} will retain the same copy number values as the input copy number.
  }
  \item{correctHOMD}{
  	If \code{TRUE}, then will correct the copy number of homozygous deletion bins and segments based on purity and ploidy corrected log ratios.
  }
  \item{gender}{
  	data.frame containing list of centromere regions. This should contain 3 columns: chr, start, and end. If this argument is used, then data at and flanking the centromeres will be removed.
  }
  \item{chrs}{
  	Chromosomes to consider for copy number correction.
  }
}

\details{
  
}

\value{
  \code{cn}: \code{\link[data.table]{data.table}} object that contains the same columns as the input object but also includes new columns \code{logR_Copy_Number}, \code{Corrected_Copy_Number}, \code{Corrected_Call}.
  \code{segs}: \code{\link[data.table]{data.table}} object that contains the same columns as the input object but also includes new columns \code{logR_Copy_Number}, \code{Corrected_Copy_Number}, \code{Corrected_Call}, \code{Corrected_MajorCN}, \code{Corrected_MinorCN}.
  Column definitions: 
  \item{logR_Copy_Number}{Purity and ploidy corrected log ratios that have been converted to a decimal-based copy number value.}
  \item{Corrected_Copy_Number}{\code{round(logR_Copy_Number)}}
  \item{Corrected_Call}{String representation of \code{Corrected_Copy_Number}; \code{HLAMP}=high-level amplification is assigned to bins/segments that have been corrected.}
  \item{Corrected_MajorCN}{Purity and ploidy corrected integer (rounded) major copy number value.}
  \item{Corrected_MinorCN}{Purity and ploidy corrected integer (rounded) minor copy number value.}
}

\references{
Ha, G., Roth, A., Khattra, J., Ho, J., Yap, D., Prentice, L. M., Melnyk, N., McPherson, A., Bashashati, A., Laks, E., Biele, J., Ding, J., Le, A., Rosner, J., Shumansky, K., Marra, M. A., Huntsman, D. G., McAlpine, J. N., Aparicio, S. A. J. R., and Shah, S. P. (2014). TITAN: Inference of copy number architectures in clonal cell populations from tumour whole genome sequence data. Genome Research, 24: 1881-1893. (PMID: 25060187)
}

\author{
  Gavin Ha <gavinha@gmail.com>
}

\seealso{
  \code{\link{outputTitanResults}}, \code{\link{outputTitanSegments}}
}

\examples{
data(EMresults)

#### COMPUTE OPTIMAL STATE PATH USING VITERBI ####
optimalPath <- viterbiClonalCN(data, convergeParams)

#### FORMAT RESULTS ####
results <- outputTitanResults(data, convergeParams, optimalPath,
                              filename = NULL, posteriorProbs = FALSE, 
                              subcloneProfiles = TRUE, correctResults = TRUE, 
                              proportionThreshold = 0.05, recomputeLogLik = FALSE,
                              proportionThresholdClonal = 0.05,
                              is.haplotypeData = FALSE)
## use corrected parameters
convergeParams <- results$convergeParam 
## use corrected results
results <- results$corrResults 
## get normal contamination and ploidy estimates
norm <- tail(convergeParams$n,1)
ploidy <- tail(convergeParams$phi,1)

#### OUTPUT SEGMENTS ####
segs <- outputTitanSegments(results, id = "test", convergeParams, 
  filename = NULL, igvfilename = NULL)
corrIntCN.results <- correctIntegerCN(results, segs, 1 - norm, ploidy, maxCNtoCorrect.autosomes = NULL, 
		maxCNtoCorrect.X = NULL, correctHOMD = TRUE, minPurityToCorrect = 0.2, gender = "female", chrs = 2)
}


\keyword{manip}
