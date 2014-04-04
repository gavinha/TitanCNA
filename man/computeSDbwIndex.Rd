\name{computeSDbwIndex}
\alias{computeSDbwIndex}

\title{
  Compute the S_Dbw Validity Index for \pkg{TitanCNA} model selection
}
\description{
  Compute the S_Dbw Validity Index internal cluster validation from the \pkg{TitanCNA} results to use for model selection.
}
\usage{
  computeSDbwIndex(x, method = "median")
}

\arguments{
  \item{x}{
    Formatted \pkg{TitanCNA} results output from \code{\link{outputTitanResults}}. See Example.
  }
  \item{method}{
    \code{median} or \code{mean} distance method to use during internal cluster validation.
  }
}
\details{
S_Dbw Validity Index is an internal clustering evaluation that is used for model selection (Halkidi et al. 2002). It attempts to choose the model that minimizes within cluster variances (scat) and maximizes density-based cluster separation (Dens). Then, S_Dbw(\|c_{T}\|x z)=Dens(\|c_{T}\|x z)+scat(\|c_{T}\|x z). 

In the context of \pkg{TitanCNA}, the copy number log ratio represents the internal data and the resulting joint states of copy number (c_{T}, forall c_{T} in \{0 : 5\}) and clonal cluster (z) make up the clusters in the internal evaluation. The optimal \pkg{TitanCNA} run is chosen as the run with the minimum S_Dbw. 
}
\value{
  \code{\link[base:list]{list}} with components: 
  \item{dens.bw}{density component of S_Dbw index}
  \item{scat}{scatter component of S_Dbw index}
  \item{S_DbwIndex}{Sum of dens.bw and scat.}
}

\references{
  Halkidi, M., Batistakis, Y., and Vazirgiannis, M. (2002). Clustering validity checking methods: part ii. SIGMOD Rec., 31(3):19--27.

 G. Ha, A. Roth, N. Melnyk, A. McPherson, L. M. Prentice, A. Bashashati, J. Ding, A. Le, J. Rosner, K. Shumansky, M. A. Marra, D. G. Huntsman, J. N. McAlpine, S. A. J. R. Aparicio, S. P. Shah. Modeling the evolution of copy number architecture in tumour whole genome sequencing data. Under Review 2013.
}
\author{
  Gavin Ha <gavinha@gmail.com>
}

\seealso{
  \code{\link{outputModelParameters}}
}

\examples{
data(EMresults)

#### COMPUTE OPTIMAL STATE PATH USING VITERBI ####
#options(cores=1)
optimalPath <- viterbiClonalCN(data, convergeParams)

#### FORMAT RESULTS ####
results <- outputTitanResults(data, convergeParams, optimalPath, 
                              filename = NULL, posteriorProbs = FALSE)

#### COMPUTE S_Dbw Validity Index FOR MODEL SELECTION ####
s_dbw <- computeSDbwIndex(results, method = "median")
}

\keyword{manip}
