\name{computeSDbwIndex}
\alias{computeSDbwIndex}

\title{
  Compute the S_Dbw Validity Index for \pkg{TitanCNA} model selection
}
\description{
  Compute the S_Dbw Validity Index internal cluster validation from the \pkg{TitanCNA} results to use for model selection.
}
\usage{
  computeSDbwIndex(x, centroid.method = "median", 
  				data.type = "LogRatio", use.corrected.cn =TRUE,
  				S_Dbw.method = "Halkidi", symmetric = TRUE)
}

\arguments{
  \item{x}{
    Formatted \pkg{TitanCNA} results output from \code{\link{outputTitanResults}}. See Example.
  }
  \item{centroid.method}{
    \code{median} or \code{mean} method to compute cluster centroids during internal cluster validation.
  }
  \item{data.type}{
  	Compute S_Dbw validity index based on copy number (use \sQuote{\code{LogRatio}}) or allelic ratio (use \sQuote{\code{AllelicRatio}}). 
  }
  \item{symmetric}{
  	\code{TRUE} if the TITAN analysis was carried out using symmetric genotypes.  See \code{\link{loadAlleleCounts}}.
  }
  \item{S_Dbw.method}{
  	Compute S_Dbw validity index using \code{Halkidi} or \code{Tong} method. See details and references.
  }
  \item{use.corrected.cn}{
  	\code{TRUE}: Will use corrected copy number calls for computing S_Dbw validity index. 
  }
}
\details{
S_Dbw Validity Index is an internal clustering evaluation that is used for model selection (Halkidi et al. 2002). It attempts to choose the model that minimizes within cluster variances (scat) and maximizes density-based cluster separation (Dens). Then, S_Dbw(|c_{T}|x z)=Dens(|c_{T}|x z)+scat(|c_{T}|x z). 

In the context of \pkg{TitanCNA}, if \code{data.type}=\sQuote{\code{LogRatio}}, then the S_Dbw internal data consists of copy number log ratios, and the resulting joint states of copy number (c_{T}, forall c_{T} in \{0 : max.copy.number\}) and clonal cluster (z) make up the clusters in the internal evaluation. If \code{data.type}=\sQuote{\code{AllelicRatio}}, then the S_Dbw internal data consists of the allelic ratios. The optimal \pkg{TitanCNA} run is chosen as the run with the minimum S_Dbw. If \code{data.type}=\sQuote{\code{Both}}, then the sum of the S_Dbw for \sQuote{\code{LogRatio}} and \sQuote{\code{AllelicRatio}} are added together. This helps account for both data types when choosing the optimal solution.

Note that for \code{S_Dbw.method}, the \code{Tong} method has an incorrect formulation of the \code{scat(c)} function. The function should be a weighted sum, but that is not the formulation shown in the publication. \code{computeSDbwIndex} uses \code{(ni/N)} instead of \code{(N-ni)/N} in the \code{scat} formula, where \code{ni} is the number of datapoints in cluster \code{i} and \code{N} is the total number of datapoints.
}
\value{
  \code{\link[base:list]{list}} with components: 
  \item{dens.bw}{density component of S_Dbw index}
  \item{scat}{scatter component of S_Dbw index}
  \item{S_DbwIndex}{Sum of dens.bw and scat.}
}

\references{
	Halkidi, M., Batistakis, Y., and Vazirgiannis, M. (2002). Clustering validity checking methods: part ii. SIGMOD Rec., 31(3):19--27.
  
	Tong, J. and Tan, H. Clustering validity based on the improved S_Dbw* index. (2009). Journal of Electronics (China), Volume 26, Issue 2, pp 258-264.

	Ha, G., Roth, A., Khattra, J., Ho, J., Yap, D., Prentice, L. M., Melnyk, N., McPherson, A., Bashashati, A., Laks, E., Biele, J., Ding, J., Le, A., Rosner, J., Shumansky, K., Marra, M. A., Huntsman, D. G., McAlpine, J. N., Aparicio, S. A. J. R., and Shah, S. P. (2014). TITAN: Inference of copy number architectures in clonal cell populations from tumour whole genome sequence data. Genome Research, 24: 1881-1893. (PMID: 25060187)
}
\author{
  Gavin Ha <gavinha@gmail.com>
}

\seealso{
  \code{\link{outputModelParameters}}, \code{\link{loadAlleleCounts}}
}

\examples{
data(EMresults)

#### COMPUTE OPTIMAL STATE PATH USING VITERBI ####
#options(cores=1)
optimalPath <- viterbiClonalCN(data, convergeParams)

#### FORMAT RESULTS ####
results <- outputTitanResults(data, convergeParams, optimalPath, 
                              filename = NULL, posteriorProbs = FALSE,
                              correctResults = TRUE, 
                              proportionThreshold = 0.05,
															proportionThresholdClonal = 0.05)
															
results <- results$corrResults ## use corrected results
#### COMPUTE S_Dbw Validity Index FOR MODEL SELECTION ####
s_dbw <- computeSDbwIndex(results, data.type = "LogRatio", 
					centroid.method = "median", S_Dbw.method = "Tong")
}

\keyword{manip}
