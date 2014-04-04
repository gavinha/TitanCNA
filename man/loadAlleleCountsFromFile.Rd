\name{loadAlleleCountsFromFile}
\alias{loadAlleleCountsFromFile}

\title{
  Function to load tumour allele counts from a text file
}

\description{
  Function to load in the allele counts from tumour sequencing data from a delimited text file.
}

\usage{
  loadAlleleCountsFromFile(infile, symmetric = TRUE, sep = "\t")
}

\arguments{
  \item{infile}{
    Full file path to text file containing tumour allele count data. \code{infile} must be 6 columns: chromosome, position, reference base, reference read counts, non-reference base, non-reference read counts.  \sQuote{chromosome} column must not contain the letters \dQuote{chr}, only the chromosome number [1-22,X,Y,M or MT]. The reference and non-reference base columns can be any arbitrary character; it is not used by \pkg{TitanCNA}.
  }
  \item{symmetric}{
    \code{logical}; if \code{TRUE}, then the symmetric allelic counts will be used. \code{ref} will equal \code{max(ref,nonRef)}.  
  }
  \item{sep}{
    Character indicating the delimiter used for the columns for \code{infile}. Default is tab-delimited, "\\t".
  }
}

\value{
  \code{\link[base:list]{list}} containing components for 
  \item{chr}{Chromosome; numeric with chrX=23, chrY=24, chrM=25}
  \item{posn}{Position; numeric}
  \item{ref}{Reference counts; numeric}
  \item{nonRef}{Non-reference counts; numeric}
  \item{tumDepth}{Tumour depth; numeric}
  \item{normDepth}{Optionally (see Note) normal depth}
}

\references{
  G. Ha, A. Roth, N. Melnyk, A. McPherson, L. M. Prentice, A. Bashashati, J. Ding, A. Le, J. Rosner, K. Shumansky, M. A. Marra, D. G. Huntsman, J. N. McAlpine, S. A. J. R. Aparicio, S. P. Shah. Modeling the evolution of copy number architecture in tumour whole genome sequencing data. Under Review 2013.
}

\author{
  Gavin Ha <gavinha@gmail.com>
}

\note{
  An additional 7th column can be included to represent the read counts for the normal genome. This is not currently used by \pkg{TitanCNA} but is supported for future extensions.
}

\seealso{
  \code{\link{loadDefaultParameters}}
}


\examples{
  infile <- system.file("extdata", "test_alleleCounts_chr2.txt", 
                        package = "TitanCNA")
  #### LOAD DATA ####
  data <- loadAlleleCountsFromFile(infile)
}


\keyword{IO}
