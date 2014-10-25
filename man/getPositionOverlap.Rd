\name{getPositionOverlap}
\alias{getPositionOverlap}

\title{
  Function to assign values to given chromosome-position that overlaps a list of chromosomal segments 
}

\description{
  Given a list of chromosomes and positions, uses a C-based function that searches a list of segments to find the overlapping segment. Then, takes the value (4th column in segment data.frame) of the overlapping segment and assigns to the given chromosome and position.
}

\usage{
  getPositionOverlap(chr, posn, dataVal)
}

\arguments{
  \item{chr}{
    \code{Numeric \link[base:array]{array}} denoting the chromosome for a list of positions. Must have the same number of elements as \code{posn}.
  }
  \item{posn}{
    \code{Numeric \link[base:array]{array}} denoting the position in the chromosome for a list of positions. Must have the same number of elements as \code{chr}.
  }
  \item{dataVal}{
    \code{\link[base:data.frame]{data.frame}} containing a list of segments described with 4 columns: chromosome, start coordinate, end coordinate, value of interest (e.g. log ratio). Chromosome can be all numeric or chrX and chrY can use \sQuote{X} and \sQuote{Y}.
  }
}

\value{
  \code{Numeric \link[base:array]{array}} of values from the 4th column of \code{\link[base:data.frame]{data.frame} cnData}. Each element corresponds to a genomic location from \code{chr} and \code{posn} that overlapped the segment in \code{cnData}. 
}

\references{
Ha, G., Roth, A., Khattra, J., Ho, J., Yap, D., Prentice, L. M., Melnyk, N., McPherson, A., Bashashati, A., Laks, E., Biele, J., Ding, J., Le, A., Rosner, J., Shumansky, K., Marra, M. A., Huntsman, D. G., McAlpine, J. N., Aparicio, S. A. J. R., and Shah, S. P. (2014). TITAN: Inference of copy number architectures in clonal cell populations from tumour whole genome sequence data. Genome Research, Published in advance July 24 2014. (PMID: 25060187)
}

\author{
  Gavin Ha <gavinha@gmail.com>
}

\seealso{
  \code{\link{loadAlleleCounts}}, \code{\link{correctReadDepth}}
}

\examples{
infile <- system.file("extdata", "test_alleleCounts_chr2.txt", 
                      package = "TitanCNA")
tumWig <- system.file("extdata", "test_tum_chr2.wig", package = "TitanCNA")
normWig <- system.file("extdata", "test_norm_chr2.wig", package = "TitanCNA")
gc <- system.file("extdata", "gc_chr2.wig", package = "TitanCNA")
map <- system.file("extdata", "map_chr2.wig", package = "TitanCNA")

#### LOAD DATA ####
data <- loadAlleleCounts(infile)

#### GC AND MAPPABILITY CORRECTION ####
cnData <- correctReadDepth(tumWig, normWig, gc, map)

#### READ COPY NUMBER FROM HMMCOPY FILE ####
logR <- getPositionOverlap(data$chr, data$posn, cnData)
}


\keyword{manip}
