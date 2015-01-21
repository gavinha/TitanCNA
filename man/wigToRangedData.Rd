\name{WIG Import Functions}
\alias{wigToRangedData}
\title{
  WIG Import Functions
}
\description{
  Fast fixedStep WIG file reading and parsing
}
\usage{
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
  \code{\link[IRanges:RangedData-class]{RangedData}} for
  \code{\link{wigToRangedData}} with chromosome and position information, sorted
  in decreasing chromosal size and increasing position.

}
\references{
  \describe{
    \item{WIG}{\url{http://genome.ucsc.edu/goldenPath/help/wiggle.html}}
  }
}
\author{
  Daniel Lai
}

\seealso{
  \code{\link{wigToRangedData}} is a wrapper around these functions for easy
  WIG file loading and structure formatting. It is taken from the \pkg{HMMcopy} package.
}

\examples{
	map <- system.file("extdata", "map_chr2.wig", package = "TitanCNA")
	mScore <- as.data.frame(wigToRangedData(map))
}
\keyword{IO}
