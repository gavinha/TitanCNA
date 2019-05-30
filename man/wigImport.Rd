\name{WIG Import Functions}
\alias{wigToRangedData}
\alias{wigToGRanges}
\title{
  WIG Import Functions. wigToGRanges (new) and wigToRangedData (deprecated)
}
\description{
  Fast fixedStep WIG file reading and parsing
}
\usage{
wigToGRanges(wigfile, verbose = TRUE)
wigToRangedData(wigfile, verbose = TRUE)
}
\arguments{
  \item{wigfile}{
    Filepath to fixedStep WIG format file
}
  \item{verbose}{
    Set to FALSE to suppress messages
}
}
\details{
  Reads the entire file into memory, then processes the file in rapid fashion,
  thus performance will be limited by memory capacity.

  The WIG file is expected to conform to the minimal fixedStep WIG format
  (see References), where each chromsome is started by a \dQuote{fixedStep}
  declaration line. The function assumes only a single track in the WIG file,
  and will ignore any lines before the first line starting with \dQuote{fixedStep}.
}
\value{
  \code{\link{GRanges}} object with chromosome and position information, sorted
  in decreasing chromosal size and increasing position.

}
\references{
  \describe{
    \item{WIG}{\url{http://genome.ucsc.edu/goldenPath/help/wiggle.html}}
  }
}
\author{
  Gavin Ha, Daniel Lai
}

\seealso{
  \code{\link{wigToGRanges}} is a wrapper around these functions for easy
  WIG file loading and structure formatting. It is modified from the \pkg{HMMcopy} package.
}

\examples{
	map <- system.file("extdata", "map_chr2.wig", package = "TitanCNA")
	mScore <- as.data.frame(wigToGRanges(map))
}
\keyword{IO}
