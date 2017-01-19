\name{Formatting and output of Titan results}
\alias{outputTitanResults}
\alias{outputModelParameters}
\alias{outputTitanSegments}

\title{
  Formatting and printing \pkg{TitanCNA} results.
}
\description{
  Function to format \pkg{TitanCNA} results in to a data.frame and output the results to a tab-delimited file.
}

\usage{
  outputTitanResults(data, convergeParams, optimalPath, filename = NULL, 
      posteriorProbs = FALSE, subcloneProfiles = TRUE)
      
  outputModelParameters(convergeParams, results, filename, 
  		S_Dbw.scale = 1, S_Dbw.method = "Tong")
  
  outputTitanSegments(results, id, convergeParams, filename = NULL, 
  		igvfilename = NULL)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{id}{
		Character string identifier for sample 
	}
  \item{data}{
    \code{\link[base:list]{list}} object that contains the components for the data to be analyzed. \code{chr}, \code{posn}, \code{ref}, and \code{tumDepth} that can be obtained using \code{\link{loadAlleleCounts}}, and \code{logR} that can be obtained using \code{\link{correctReadDepth}} and \code{\link{getPositionOverlap}} (see Example).
  }
  \item{convergeParams}{
    \code{\link[base:list]{list}} object that is returned from the function \code{\link{runEMclonalCN}} in \pkg{TitanCNA}.      
  }
  \item{optimalPath}{
    \code{numeric \link[base:array]{array}} containing the optimal \pkg{TitanCNA} genotype and clonal cluster states for each data point in the analysis. \code{optimalPath} is obtained from running \code{\link{viterbiClonalCN}}.  
  }
  \item{results}{
    Formatted \pkg{TitanCNA} results output from \code{\link{outputTitanResults}}.
  }
  \item{filename}{
    Path of the file to write the \pkg{TitanCNA} results.
  }
  \item{igvfilename}{
  	Path of the file to write the IGV seg file.
  }
  \item{posteriorProbs}{
    \code{Logical TRUE} to include the posterior marginal probabilities in printing to \code{filename}. 
  }
  \item{subcloneProfiles}{
  	\code{Logical TRUE} to include the subclone profiles to the output \code{data.frame}. Currently, this only works for 1 or 2 clonal clusters.
  }
  \item{S_Dbw.scale}{
  	The S_Dbw validity index can be adjusted to account for differences between datasets. \code{SDbw.scale} can be used to penalize the S_Dbw \code{dens.bw} component.  The default is 1.
  }
  \item{S_Dbw.method}{
  	Compute S_Dbw validity index using \code{Halkidi} or \code{Tong} method. See \code{\link{computeSDbwIndex}}.
  }
}
\details{
  \code{\link{outputModelParameters}} outputs to a file with the estimated TITAN model parameters and model selection index. Each row contains information regarding different parameters:
  
1) Normal contamination estimate - proportion of normal content in the sample; tumour content is 1 minus this number

2) Average tumour ploidy estimate - average number of estimated copies in the genome; 2 represents diploid

3) Clonal cluster cellular prevalence - Z denotes the number of clonal clusters; each value (space-delimited) following are the cellular prevalence estimates for each cluster. Cellular prevalence here is defined as the proportion of tumour sample that does contain the aberrant genotype.

4) Genotype binomial means for clonal cluster Z - set of 21 binomial estimated parameters for each specified cluster

5) Genotype Gaussian means for clonal cluster Z - set of 21 Gaussian estimated means for each specified cluster

6) Genotype Gaussian variance - set of 21 Gaussian estimated variances; variances are shared for across all clusters

7) Number of iterations - number of EM iterations needed for convergence

8) Log likelihood - complete data log-likelihood for current cluster run

9) S_Dbw dens.bw - density component of S_Dbw index; see \code{\link{computeSDbwIndex}}

10) S_Dbw scat - scatter component of S_Dbw index; see \code{\link{computeSDbwIndex}}

11) S_Dbw validity index - used for model selection where the run with optimal number of clusters based on lowest S_Dbw index. This value is slightly modified from that computed from \code{\link{computeSDbwIndex}}. It is computed as S_Dbw= S_Dbw.scale * dens.bw + scat

12) S_Dbw dens.bw, scat, validity index is computed for \code{LogRatio} and \code{AllelicRatio} datatypes, as well as the combination of \code{Both}. For \code{Both}, the values are summed for both datatypes.

  \code{\link{outputTitanResults}} outputs a file that has the similar format described in \sQuote{Value} section.
}

\value{
  \code{\link{outputTitanResults}} also returns a \code{\link[base:data.frame]{data.frame}}, where each row corresponds to a position in the analysis, and with the following columns:
  \item{Chr}{character denoting chromosome number.  ChrX and ChrY uses \sQuote{X} and \sQuote{Y}.}
  \item{Position}{genomic coordinate}
  \item{RefCount}{number of reads matching the reference base}
  \item{NRefCount}{number of reads matching the non-reference base}
  \item{Depth}{total read depth at the position}
  \item{AllelicRatio}{RefCount/Depth}
  \item{LogRatio}{log2 ratio between normalized tumour and normal read depths}
  \item{CopyNumber}{predicted \pkg{TitanCNA} copy number}
  \item{TITANstate}{internal state number used by \pkg{TitanCNA}; see Reference}
  \item{TITANcall}{interpretable \pkg{TitanCNA} state; string (HOMD,DLOH,HET,NLOH,ALOH,ASCNA,BCNA,UBCNA); See Reference}
  \item{ClonalCluster}{predicted \pkg{TitanCNA} clonal cluster; lower cluster numbers represent clusters with higher cellular prevalence}
  \item{CellularPrevalence}{proportion of tumour cells containing event; not to be mistaken as proportion of sample (including normal)}
  
  If \code{subcloneProfiles} is set to \code{TRUE}, then the subclone profiles will be appended to the output \code{data.frame}. 
  \item{Subclone1.CopyNumber}{Integer copy number for Subclone 1.}
  \item{Subclone1.TITANcall}{States for Subclone 1}
  \item{Subclone1.Prevalence}{The cellular prevalence of Subclone 1, or sometimes referred to as the subclone fraction.}
  
  \code{\link{outputModelParameters}} returns a \code{\link[base:list]{list}} containing the S_Dbw model selection:
  	\item{dens.bw}{}
  	\item{scat}{}
  	\item{S_Dbw}{S_Dbw.scale * dens.bw + scat}
}

\references{
Ha, G., Roth, A., Khattra, J., Ho, J., Yap, D., Prentice, L. M., Melnyk, N., McPherson, A., Bashashati, A., Laks, E., Biele, J., Ding, J., Le, A., Rosner, J., Shumansky, K., Marra, M. A., Huntsman, D. G., McAlpine, J. N., Aparicio, S. A. J. R., and Shah, S. P. (2014). TITAN: Inference of copy number architectures in clonal cell populations from tumour whole genome sequence data. Genome Research, 24: 1881-1893. (PMID: 25060187)
}

\author{
  Gavin Ha <gavinha@gmail.com>
}

\seealso{
  \code{\link{runEMclonalCN}}, \code{\link{viterbiClonalCN}}, \code{\link{computeSDbwIndex}}
}

\examples{
data(EMresults)

#### COMPUTE OPTIMAL STATE PATH USING VITERBI ####
optimalPath <- viterbiClonalCN(data, convergeParams)

#### FORMAT RESULTS ####
results <- outputTitanResults(data, convergeParams, optimalPath, 
                              filename = NULL, posteriorProbs = FALSE,
                              subcloneProfiles = TRUE)

#### OUTPUT RESULTS TO FILE ####
outparam <- paste0("cluster2_params.txt")
outputModelParameters(convergeParams, results, outparam)

#### OUTPUT SEGMENTS TO FILE ####
outseg <- paste0("cluster2_segs.txt")
outigv <- paste0("cluster2.seg")
segs <- outputTitanSegments(results, id = "test", convergeParams, filename = outseg, igvfilename = outigv)
# segment results also stored in data.frame "segs" 
}


\keyword{IO}
\keyword{manip}
