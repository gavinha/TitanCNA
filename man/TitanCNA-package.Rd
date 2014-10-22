\name{TitanCNA-package}
\alias{TitanCNA-package}
\alias{TitanCNA}
\docType{package}
\title{
TITAN: Subclonal copy number and LOH prediction whole genome sequencing of tumours
}
\description{
  TITAN is a software tool for inferring subclonal copy number alterations (CNA) and loss of heterozygosity (LOH). The algorithm also infers clonal group cluster membership for each event and the tumour proportion, or cellular prevalence, for each event. 
}
\details{


\tabular{ll}{
Package: \tab TitanCNA\cr
Type: \tab Package\cr
Version: \tab 1.5.2\cr
Date: \tab 2014-10-10\cr
License: \tab see LICENSE\cr
}

\code{example("TitanCNA-package")} for quick tour of functionality and visualization

\code{vignette("TitanCNA")} for detailed example

}

\author{
Gavin Ha, Sohrab P Shah
Maintainer: Gavin Ha <gavinha@gmail.com>
}

\references{
Ha, G., Roth, A., Khattra, J., Ho, J., Yap, D., Prentice, L. M., Melnyk, N., McPherson, A., Bashashati, A., Laks, E., Biele, J., Ding, J., Le, A., Rosner, J., Shumansky, K., Marra, M. A., Huntsman, D. G., McAlpine, J. N., Aparicio, S. A. J. R., and Shah, S. P. (2014). TITAN: Inference of copy number architectures in clonal cell populations from tumour whole genome sequence data. Genome Research, Published in advance July 24 2014. (PMID: 25060187)
}

\keyword{IO}
\keyword{manip}
\keyword{package}

\examples{
message('Running TITAN ...')
#### LOAD DATA ####
infile <- system.file("extdata", "test_alleleCounts_chr2.txt", package = "TitanCNA")
data <- loadAlleleCounts(infile)

#### LOAD PARAMETERS ####
message('titan: Loading default parameters')
numClusters <- 2
params <- loadDefaultParameters(copyNumber = 5, 
                                numberClonalClusters = numClusters, skew = 0.1)

#### READ COPY NUMBER FROM HMMCOPY FILE ####
message('titan: Correcting GC content and mappability biases...')
tumWig <- system.file("extdata", "test_tum_chr2.wig", package = "TitanCNA")
normWig <- system.file("extdata", "test_norm_chr2.wig", package = "TitanCNA")
gc <- system.file("extdata", "gc_chr2.wig", package = "TitanCNA")
map <- system.file("extdata", "map_chr2.wig", package = "TitanCNA")
cnData <- correctReadDepth(tumWig, normWig, gc, map)
logR <- getPositionOverlap(data$chr, data$posn, cnData)
data$logR <- log(2^logR) #transform to natural log

#### FILTER DATA FOR DEPTH, MAPPABILITY, NA, etc ####
data <- filterData(data, c(1:22,"X"), minDepth = 10, maxDepth = 200, map = NULL)

#### EM (FWD-BACK) TO TRAIN PARAMETERS ####
#### Can use parallelization packages ####
K <- length(params$genotypeParams$alphaKHyper)
params$genotypeParams$alphaKHyper <- rep(500, K)
params$ploidyParams$phi_0 <- 1.5 
convergeParams <- runEMclonalCN(data, gParams = params$genotypeParams, 
                                nParams = params$normalParams, 
                                pParams = params$ploidyParams, 
                                sParams = params$cellPrevParams, 
                                maxiter = 3, maxiterUpdate = 500, 
                                txnExpLen = 1e9, txnZstrength = 1e9, 
                                useOutlierState = FALSE, 
                                normalEstimateMethod = "map", 
                                estimateS = TRUE, estimatePloidy = TRUE)                                
#### COMPUTE OPTIMAL STATE PATH USING VITERBI ####
optimalPath <- viterbiClonalCN(data, convergeParams)

#### FORMAT RESULTS ####
results <- outputTitanResults(data, convergeParams, optimalPath, 
                              filename = NULL, posteriorProbs = FALSE,
                              subcloneProfiles = TRUE)

#### PLOT RESULTS ####
norm <- tail(convergeParams$n, 1)
ploidy <- tail(convergeParams$phi, 1)

par(mfrow=c(4, 1))    
plotCNlogRByChr(results, chr = 2, ploidy = ploidy, geneAnnot = NULL, 
                ylim = c(-2, 2), cex = 0.5, xlab = "", main = "Chr 2")
plotAllelicRatio(results, chr = 2, geneAnnot = NULL, ylim = c(0, 1), cex = 0.5, 
                xlab = "", main = "Chr 2")
plotClonalFrequency(results, chr = 2, normal = norm, geneAnnot = NULL, 
                    ylim = c(0, 1), cex = 0.5, xlab = "", main = "Chr 2")
plotSubcloneProfiles(results, chr = 2, cex = 2, main = "Chr 2")
}
