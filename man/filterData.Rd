\name{filterData}
\alias{filterData}

\title{
  Filter list object based on read depth and missing data 
}
\description{
  Filters all vectors in list based on specified chromosome(s) of interest, minimum and maximum read depths, missing data, mappability score threshold
}

\usage{
filterData(data ,chrs = NULL, minDepth = 10, maxDepth = 200, 
    positionList = NULL, map = NULL, mapThres = 0.9)
}


\arguments{
  \item{data}{
    \code{\link[base:list]{list}} object that contains an arbitrary number of components. Should include \sQuote{chr}, \sQuote{tumDepth}. All vector elements must have the same number of rows where each row corresponds to information pertaining to a chromosomal position.
  }
  \item{chrs}{
    \code{character} or vector of \code{character} specifying the chromosomes to keep. Chromosomes not included in this \code{array} will be filtered out. Chromosome style must match the \code{genomeStyle} used when running \code{\link{loadAlleleCounts}}
  }
  \item{minDepth}{
    \code{Numeric integer} specifying the minimum tumour read depth to include. Positions >= \code{minDepth} are kept. 
  }
  \item{maxDepth}{
    \code{Numeric integer} specifying the maximum tumour read depth to include. Positions <= \code{maxDepth} are kept.
  }
  \item{positionList}{
    \code{\link[base:data.frame]{data.frame}} with two columns: \sQuote{chr} and \sQuote{posn}. \code{positionList} lists the chromosomal positions to use in the analysis. All positions not overlapping this list will be excluded. Use \code{NULL} to use all current positions in \code{data}.
  }
  \item{map}{
    \code{Numeric array} containing map scores corresponding to each position in \code{data}. Optional for filtering positions based on mappability scores.
  }
  \item{mapThres}{
    \code{Numeric float} specifying the mappability score threshold. Only applies if \code{map} is specified. \code{map} scores >= \code{mapThres} are kept.
  }
}

\details{
  All vectors in the input \code{data} list object, and \code{map}, must all have the same number of rows.
}

\value{
  The same \code{\link[base:list]{list} data} containing filtered components.
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
infile <- system.file("extdata", "test_alleleCounts_chr2.txt", 
                      package = "TitanCNA")
tumWig <- system.file("extdata", "test_tum_chr2.wig", package = "TitanCNA")
normWig <- system.file("extdata", "test_norm_chr2.wig", package = "TitanCNA")
gc <- system.file("extdata", "gc_chr2.wig", package = "TitanCNA")
map <- system.file("extdata", "map_chr2.wig", package = "TitanCNA")

#### LOAD DATA ####
data <-  loadAlleleCounts(infile, genomeStyle = "NCBI")

#### GC AND MAPPABILITY CORRECTION ####
cnData <- correctReadDepth(tumWig, normWig, gc, map)


#### READ COPY NUMBER FROM HMMCOPY FILE ####
logR <- getPositionOverlap(data$chr, data$posn, cnData)
data$logR <- log(2^logR) #use natural logs

#### FILTER DATA FOR DEPTH, MAPPABILITY, NA, etc ####
filtereData <- filterData(data, as.character(1:24), minDepth = 10, 
				maxDepth = 200, map = NULL, mapThres=0.9)
}


\keyword{manip}
