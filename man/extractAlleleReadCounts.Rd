\name{extractAlleleReadCounts}
\alias{extractAlleleReadCounts}

\title{
  Function to extract allele read counts from a sequence alignment (BAM) file
}

\description{
  Function to extract allele read counts from a sequence alignment (BAM) file at specific positions of interest. The positions are passed in as the file path to a file in variant call format (VCF).
}

\usage{
  extractAlleleReadCounts(bamFile, bamIndex, positions, 
  		outputFilename = NULL, pileupParam = PileupParam())
}

\arguments{
  \item{bamFile}{
    File path location to the sequencing alignment file (BAM format) from which to extract read counts.
  }
  \item{bamIndex}{
     File path location to the BAM index file (usually with extension .bai) corresponding to the sequencing alignment file \code{bamFile}.
  }
  \item{positions}{
  	File path location to the variant call format (VCF) file containing the positions at which read counts are to be extracted.
  }
  \item{outputFilename}{
    If given, will specify the file path to which the result will be output as tab-delimited text. Otherwise, the no output is written to file.
  }
  \item{pileupParam}{
  	\code{\link[Rsamtools:PileupParam]{PileupParam}} object from the \pkg{Rsamtools}. See Details.  
  }
}

\value{
  \code{\link[base:data.frame]{data.frame}} containing columns: 
  \item{chr}{Chromosome; character}
  \item{position}{Position; numeric}
  \item{ref}{Reference counts; character}
  \item{refCount}{Reference counts; numeric}
  \item{Nref}{Non-reference counts; character}
  \item{NrefCount}{Non-reference counts; numeric}
}

\details{
	The \code{pileupParam} object allows users to specify the sequencing parameters to consider when generating the pileup from which read counts are extracted. This includes \sQuote{max_depth}, \sQuote{min_base_quality}, \sQuote{min_mapq}, \sQuote{min_nucleotide_depth}=10 (recommended), \sQuote{min_minor_allele_depth}, \sQuote{distinguish_strands}, \sQuote{distinguish_nucleotides}=TRUE (must be TRUE).
	
}

\references{
 Ha, G., Roth, A., Khattra, J., Ho, J., Yap, D., Prentice, L. M., Melnyk, N., McPherson, A., Bashashati, A., Laks, E., Biele, J., Ding, J., Le, A., Rosner, J., Shumansky, K., Marra, M. A., Huntsman, D. G., McAlpine, J. N., Aparicio, S. A. J. R., and Shah, S. P. (2014). TITAN: Inference of copy number architectures in clonal cell populations from tumour whole genome sequence data. Genome Research, Published in advance July 24 2014. (PMID: 25060187)
}

\author{
  Gavin Ha <gavinha@gmail.com>
}

\seealso{
  \code{\link[Rsamtools:PileupParam]{PileupParam}}; \url{http://samtools.sourceforge.net/}
}


\examples{
  \dontrun{
	countsDF <- extractAlleleReadCounts(bamFile, bamIndex, 
				positions, outputFilename = NULL, 
				pileupParam = PileupParam())
	data <- loadAlleleCounts(countsDF, symmetric = TRUE, 
  				genomeStyle = "NCBI")
  }
}


\keyword{IO}
