\name{loadAlleleCounts}
\alias{loadAlleleCounts}
\alias{setGenomeStyle}

\title{
  Function to load tumour allele counts from a text file or data.frame and returns a \link[data.table]{data.table}.
}

\description{
  Function to load in the allele counts from tumour sequencing data from a delimited text file or data.frame object.
}

\usage{
  loadAlleleCounts(inCounts, symmetric = TRUE, 
  		genomeStyle = "NCBI", sep = "\t", header = TRUE)

	setGenomeStyle(x, genomeStyle = "NCBI", species = "Homo_sapiens", filterExtraChr = TRUE)
}

\arguments{
  \item{inCounts}{
    Full file path to text file or data.frame containing tumour allele count data. \code{inCounts} must be 6 columns: chromosome, position, reference base, reference read counts, non-reference base, non-reference read counts.  \sQuote{chromosome} column can be in \sQuote{NCBI} or \sQuote{UCSC} genome style; only autosomes, sex chromosomes, and mitochondrial chromosome are included (e.g. 1-22,X,Y,MT). The reference and non-reference base columns can be any arbitrary character; it is not used by \pkg{TitanCNA}.
  }
  \item{symmetric}{
    \code{logical}; if \code{TRUE}, then the symmetric allelic counts will be used. \code{ref} will equal \code{max(ref,nonRef)}. This parameter must be the same as the \code{symmetric} parameter for \code{\link{loadDefaultParameters}}. 
  }
  \item{genomeStyle}{
  	The genome style to use for chromosomes by \pkg{TitanCNA}. Use one of \sQuote{NCBI} or \sQuote{UCSC}. It does not matter what style is found in \code{inCounts}, \code{genomeStyle} will be the style returned. 
  }
  \item{sep}{
    Character indicating the delimiter used for the columns for \code{infile}. Default is tab-delimited, "\\t".
  }
  \item{header}{
  	\code{logical} to indicate if the input tumour counts file contains a header line.
  }
  \item{x}{
    \code{character} vector of chromosome names to change.
  }
  \item{species}{
  	\code{character} denoting the species
  }
  \item{filterExtraChr}{
  	\code{logical}; if \code{TRUE}, then will return the list of chromosomes given by \code{extractSeqlevelsByGroup} for the \code{species} and for autosomes and sex chromosomes, which means that only the major chromosomes are returned (i.e. 1:22, X, Y). 
  }
}

\value{
	\code{loadAlleleCounts} returns a \link[data.table]{data.table} containing components for 
  \item{chr}{Chromosome; character, NCBI or UCSC genome style format}
  \item{posn}{Position; integer}
  \item{ref}{Reference counts; numeric}
  \item{nonRef}{Non-reference counts; numeric}
  \item{tumDepth}{Tumour depth; numeric}
}

\references{
Ha, G., Roth, A., Khattra, J., Ho, J., Yap, D., Prentice, L. M., Melnyk, N., McPherson, A., Bashashati, A., Laks, E., Biele, J., Ding, J., Le, A., Rosner, J., Shumansky, K., Marra, M. A., Huntsman, D. G., McAlpine, J. N., Aparicio, S. A. J. R., and Shah, S. P. (2014). TITAN: Inference of copy number architectures in clonal cell populations from tumour whole genome sequence data. Genome Research, 24: 1881-1893. (PMID: 25060187)
}

\author{
  Gavin Ha <gavinha@gmail.com>
}

\seealso{
  \code{\link{loadDefaultParameters}}
}


\examples{
  infile <- system.file("extdata", "test_alleleCounts_chr2.txt", 
                        package = "TitanCNA")
  #### LOAD DATA FROM TEXT FILE ####
  data <- loadAlleleCounts(infile, symmetric = TRUE, 
  		genomeStyle = "NCBI", header = TRUE)
  
  ## use the UCSC chromosome naming convention instead ##
  data$chr <- setGenomeStyle(data$chr, genomeStyle = "UCSC")
  \dontrun{
	data <- loadAlleleCounts(countsDF, symmetric = TRUE, 
  			genomeStyle = "NCBI")
  }
}


\keyword{IO}
