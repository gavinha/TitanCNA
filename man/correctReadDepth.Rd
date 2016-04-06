\name{correctReadDepth}
\alias{correctReadDepth}

\title{
  Correct GC content and mappability biases in sequencing data read counts
}
\description{
  Correct GC content and mappability biases in tumour sequence read counts using Loess curve fitting. Wrapper for function in \pkg{HMMcopy}.
}
\usage{
  correctReadDepth(tumWig, normWig, gcWig, mapWig, 
  			genomeStyle = "NCBI", targetedSequence = NULL)
}

\arguments{
  \item{tumWig}{
    File path to fixedStep WIG format file for the tumour sample. See \code{\link[HMMcopy:wigToRangedData]{wigToRangedData}} in the \pkg{HMMcopy} for more details.
  }
  \item{normWig}{
    File path to fixedStep WIG format file for the normal sample.
  }
  \item{gcWig}{
    File path to fixedStep WIG format file for the GC content based on the specific reference genome sequence used. 
  }
  \item{mapWig}{
    File path to fixedStep WIG format file for the mappability scores computed on the specific reference genome used.
  }
  \item{genomeStyle}{
  	The genome style to use for chromosomes by \pkg{TitanCNA}. Use one of \sQuote{NCBI} or \sQuote{UCSC}. It does not matter what style is found in \code{inCounts}, \code{genomeStyle} will be the style returned. 
  }
  \item{targetedSequence}{
    data.frame with 3 columns: chr, start position, stop position.  Use this argument for exome capture sequencing or targeted deep sequencing data. This is experimental and may not work as desired. 
  }
}
\details{
  Wrapper for \code{\link[HMMcopy:correctReadcount]{correctReadcount}} in \pkg{HMMcopy} package. It uses a sampling of 50000 bins to find the Loess fit. Then, the log ratio for every bin is returned as the log base 2 of the ratio between  the corrected tumour read count and the corrected normal read count.
}

\value{
  \code{\link[base:data.frame]{data.frame}} containing columns:
  \item{chr}{Chromosome; uses 'X' and 'Y' for sex chromosomes}
  \item{start}{Start genomic coordinate for bin in which read count is corrected}
  \item{end}{End genomic coordinate for bin in which read count is corrected}
  \item{logR}{Log ratio, log2(tumour:normal), for bin in which read count is corrected}
}

\references{
Ha, G., Roth, A., Lai, D., Bashashati, A., Ding, J., Goya, R., Giuliany, R., Rosner, J., Oloumi, A., Shumansky, K., Chin, S.F., Turashvili, G., Hirst, M., Caldas, C., Marra, M. A., Aparicio, S., and Shah, S. P. (2012). Integrative analysis of genome wide loss of heterozygosity and monoallelic expression at nucleotide resolution reveals disrupted pathways in triple negative breast cancer. Genome Research, 22(10):1995,2007. (PMID: 22637570)

Ha, G., Roth, A., Khattra, J., Ho, J., Yap, D., Prentice, L. M., Melnyk, N., McPherson, A., Bashashati, A., Laks, E., Biele, J., Ding, J., Le, A., Rosner, J., Shumansky, K., Marra, M. A., Huntsman, D. G., McAlpine, J. N., Aparicio, S. A. J. R., and Shah, S. P. (2014). TITAN: Inference of copy number architectures in clonal cell populations from tumour whole genome sequence data. Genome Research, 24: 1881-1893. (PMID: 25060187)
}

\author{
  Gavin Ha <gavinha@gmail.com>, Daniel Lai <jujubix@cs.ubc.ca>, Yikan Wang <ykwang@bccrc.ca>
}

\seealso{
  \code{\link[HMMcopy:correctReadcount]{correctReadcount}} and \code{\link[HMMcopy:wigToRangedData]{wigToRangedData}} in the \pkg{HMMcopy} package. WIG: \url{http://genome.ucsc.edu/goldenPath/help/wiggle.html}
}

\examples{
tumWig <- system.file("extdata", "test_tum_chr2.wig", package = "TitanCNA")
normWig <- system.file("extdata", "test_norm_chr2.wig", package = "TitanCNA")
gc <- system.file("extdata", "gc_chr2.wig", package = "TitanCNA")
map <- system.file("extdata", "map_chr2.wig", package = "TitanCNA")

#### GC AND MAPPABILITY CORRECTION ####
cnData <- correctReadDepth(tumWig, normWig, gc, map)
}

\keyword{IO}
\keyword{manip}
