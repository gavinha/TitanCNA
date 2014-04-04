\name{getPositionOverlap}
\alias{getPositionOverlap}

\title{
  Function to assign values to given chromosome-position that overlaps a list of chromosomal segments 
}

\description{
  Given a list of chromosomes and positions, uses a C-based function that searches a list of segments to find the overlapping segment. Then, takes the value (4th column in segment data.frame) of the overlapping segment and assigns to the given chromosome and position.
}

\usage{
  getPositionOverlap(chr, posn, cnData)
}

\arguments{
  \item{chr}{
    \code{Numeric \link[base:array]{array}} denoting the chromosome for a list of positions. Must have the same number of elements as \code{posn}.
  }
  \item{posn}{
    \code{Numeric \link[base:array]{array}} denoting the position in the chromosome for a list of positions. Must have the same number of elements as \code{chr}.
  }
  \item{cnData}{
    \code{\link[base:data.frame]{data.frame}} containing a list of segments described with 4 columns: chromosome, start coordinate, end coordinate, value of interest (e.g. log ratio). Chromosome can be all numeric or chrX and chrY can use \sQuote{X} and \sQuote{Y}.
  }
}

\value{
  \code{Numeric \link[base:array]{array}} of values from the 4th column of \code{\link[base:data.frame]{data.frame} cnData}. Each element corresponds to a genomic location from \code{chr} and \code{posn} that overlapped the segment in \code{cnData}. 
}

\references{
  G. Ha, A. Roth, N. Melnyk, A. McPherson, L. M. Prentice, A. Bashashati, J. Ding, A. Le, J. Rosner, K. Shumansky, M. A. Marra, D. G. Huntsman, J. N. McAlpine, S. A. J. R. Aparicio, S. P. Shah. Modeling the evolution of copy number architecture in tumour whole genome sequencing data. Under Review 2013.
}

\author{
  Gavin Ha <gavinha@gmail.com>
}

\seealso{
  \code{\link{loadAlleleCountsFromFile}}, \code{\link{correctReadDepth}}
}

\examples{
infile <- system.file("extdata", "test_alleleCounts_chr2.txt", 
                      package = "TitanCNA")
tumWig <- system.file("extdata", "test_tum_chr2.wig", package = "TitanCNA")
normWig <- system.file("extdata", "test_norm_chr2.wig", package = "TitanCNA")
gc <- system.file("extdata", "gc_chr2.wig", package = "TitanCNA")
map <- system.file("extdata", "map_chr2.wig", package = "TitanCNA")

#### LOAD DATA ####
data <- loadAlleleCountsFromFile(infile)

#### GC AND MAPPABILITY CORRECTION ####
cnData <- correctReadDepth(tumWig, normWig, gc, map)

#### READ COPY NUMBER FROM HMMCOPY FILE ####
logR <- getPositionOverlap(data$chr, data$posn, cnData)
}


\keyword{manip}
