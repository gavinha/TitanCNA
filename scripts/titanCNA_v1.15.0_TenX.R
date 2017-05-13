library(optparse)

option_list <- list(
	make_option(c("--id"), type = "character", help = "Sample ID"),
	make_option(c("--hetFile"), type = "character", help = "File containing allelic read counts at HET sites. [Required]"),
	make_option(c("--cnFile"), type = "character", help = "File containing normalized coverage as log2 ratios. [Required]"),
	make_option(c("--outDir"), type = "character", help = "Output directory to output the results. [Required]"),
	make_option(c("--numClusters"), type = "integer", default = 1, help = "Number of clonal clusters. [Default: 1]"),
	make_option(c("--numCores"), type = "integer", default = 1, help = "Number of cores to use. [Default: %default]"),
	make_option(c("--ploidy_0"), type = "numeric", default = 2, help = "Initial ploidy value; float [Default: %default]"),
	make_option(c("--estimatePloidy"), type = "logical", default = TRUE, help = "Estimate ploidy; TRUE or FALSE [Default: %default]"),
	make_option(c("--normal_0"), type = "numeric", default = 0.5, help = "Initial normal contamination (1-purity); float [Default: %default]"),
	make_option(c("--estimateNormal"), type = "character", default = "map", help = "Estimate normal contamination method; string {'map', 'fixed'} [Default: %default]"),
	make_option(c("--estimateClonality"), type="logical", default=TRUE, help="Estimate cellular prevalence. [Default: %default]"),
	make_option(c("--maxCN"), type = "integer", default = 8, help = "Maximum number of copies to model; integer [Default: %default]"),
	make_option(c("--alphaK"), type = "numeric", default = 10000, help = "Hyperparameter on Gaussian variance for log ratio data; for WES, use 1000; for WGS, use 10000; float [Default: %default]"),
	make_option(c("--alphaKHigh"), type = "numeric", default = 10000, help = "Hyperparameter on Gaussian variance for extreme copy number states; for WES, use 1000; for WGS, use 10000; float [Default: %default]"),
	make_option(c("--alleleModel"), type = "character", default = "binomial", help = "Emission density to use for allelic input data (binomial or Gaussian). [Default: %default]"),
	make_option(c("--alphaR"), type = "numeric", default = 1000, help = "Hyperparaemter on the Gaussian variance for allelic fraction data; used if --alleleModel=\"Gaussian\". [Defaule: %default]"),
	make_option(c("--txnExpLen"), type = "numeric", default = 1e15, help = "Expected length of segments; higher leads to longer (less sensitive) segments; float [Default: %default]"),
	make_option(c("--txnZStrength"), type = "numeric", default = 1, help = "Expected length of clonal cluster segmentation (factor of txnExpLen); float [Default: %default]"),
	make_option(c("--minDepth"), type = "integer", default = 10, help = "Minimum read depth of a HET site to include in analysis; integer [Default: %default]"),
	make_option(c("--maxDepth"), type = "integer", default = 1000, help = "Maximum read depth of a HET site to include in analysis; integer [Default: %default]"),
	make_option(c("--skew"), type = "numeric", default=0, help = "Allelic reference skew for all states except heterozygous states (e.g. 1:1, 2:2, 3:3). Value is additive to baseline allelic ratios. float [Default: %default]"),
	make_option(c("--hetBaselineSkew"), type="numeric", default=NULL, help="Allelic reference skew for heterozygous states (e.g. 1:1, 2:2, 3:3). Value is the additive to baseline allelic ratios. float [Default: %default]"), 
	make_option(c("--genomeStyle"), type = "character", default = "NCBI", help = "NCBI or UCSC chromosome naming convention; use UCSC if desired output is to have \"chr\" string. [Default: %default]"),
	make_option(c("--chrs"), type = "character", default = "c(1:22, 'X')", help = "Chromosomes to analyze; string [Default: %default"),
	make_option(c("--mapWig"), type = "character", default = NULL, help = "Mappability score file for bin sizes matching cnfile. [Default: %default]"),
	make_option(c("--mapThres"), type = "numeric", default = 0.9, help = "Minimum mappability score threshold to use; float [Default: %default]"),
	make_option(c("--centromere"), type = "character", default=NULL, help = "Centromere gap file. [Default: %default]"),
	make_option(c("--libdir"), type = "character", default=NULL, help = "Directory containing source code. Specify if changes have been made to source code and want to over-ride package code. [Default: %default]"),
	make_option(c("--outFile"), type = "character", default = NULL, help = "Output file to write position-level file. (default uses extension: *.titan.txt]"),
	make_option(c("--outSeg"), type = "character", default = NULL, help = "Output file to write detailed segments. (default uses extension: *.segs.txt]"),
	make_option(c("--outIGV"), type = "character", default = NULL, help = "Output file to write segments for loading into IGV. (default uses extension: *.seg]"),
	make_option(c("--outParam"), type = "character", default = NULL, help = "Output file to write parameters. [Default: %default]"),
	make_option(c("--outPlotDir"), type = "character", default = NULL, help = "Output directory to save plots. [Default: %default]"),
	make_option(c("--plotYlim"), type = "character", default = "c(-2,4)", help = "The Y-axis limits to use for plotting log ratio coverage results. [Default: %default]"),
	make_option(c("--haplotypeBinSize"), type = "integer", default = 1e5, help = "Bin size for summarizing phased haplotypes. [Default: %default]"),
	make_option(c("--phaseSummarizeFun"), type = "character", default = "mean", help = "Strategy to summarize specified haplotype bins: mean, sum, SNP. [Default: %default]")
)

parseobj <- OptionParser(option_list=option_list, usage = "usage: Rscript %prog [options]")
opt <- parse_args(parseobj)
print(opt)

library(TitanCNA)
library(data.table)
library(GenomicRanges)
library(dplyr)
library(doMC)
library(SNPchip)
sessionInfo()
options(bitmapType='cairo', scipen=0)

#setwd("~/Documents/Projects/CRPC_SV/TITAN_10X")
libdir <- opt$libdir
#libdir = "~/home_unix/code/git/TitanCNA/"
if (!is.null(libdir)){
	source(paste0(libdir, "/R/plotting.R"))
	source(paste0(libdir, "/R/utils.R"))
  source(paste0(libdir, "/R/hmmClonal.R"))
  source(paste0(libdir, "/R/paramEstimation.R"))
  source(paste0(libdir, "/R/correction.R"))
  source(paste0(libdir, "/R/haplotype.R"))
}

dyn.load(paste0(libdir, "src/TitanCNA.so"))
#dyn.load("/Users/gavinha/Documents/Code/git/github/TitanCNA/src/fwd_backC_clonalCN.so")

id <- opt$id
hetfile <- opt$hetFile
cnfile <- opt$cnFile
numClusters <- opt$numClusters
numCores <- opt$numCores
ploidy_0 <- opt$ploidy_0
boolEstPloidy <- opt$estimatePloidy
norm_0 <- opt$normal_0
normEstMeth <- opt$estimateNormal
estimateS <- opt$estimateClonality
maxCN <- opt$maxCN
alphaK <- opt$alphaK
alphaHigh <- opt$alphaKHigh
alleleEmissionModel <- opt$alleleModel
alphaR <- opt$alphaR
txn_exp_len <- opt$txnExpLen
txn_z_strength <- opt$txnZStrength
mapThres <- opt$mapThres
minDepth <- opt$minDepth
maxDepth <- opt$maxDepth
skew <- opt$skew
hetBaselineSkew <- opt$hetBaselineSkew
chrs <- eval(parse(text = opt$chrs))
genomeStyle <- opt$genomeStyle
mapWig <- opt$mapWig
centromere <- opt$centromere
haplotypeBinSize <- opt$haplotypeBinSize
phaseSummarizeFun <- opt$phaseSummarizeFun
outdir <- opt$outDir
outfile <- opt$outFile
outparam <- opt$outParam
outseg <- opt$outSeg
outigv <- opt$outIGV
outplot <- opt$outPlotDir
plotYlim <- eval(parse(text = opt$plotYlim))

## check arguments ##
if (!normEstMeth %in% c("map", "fixed")){
	stop("--estimateNormal must be \"map\" or \"fixed\"")
}

### SETUP OUTPUT FILE NAMES ###
numClustersStr <- as.character(numClusters)
if (numClusters < 10) { 
	numClustersStr <- paste0("0", numClusters)
}
if (is.null(outfile)){
	outfile <- paste0(outdir, "/", id, "_cluster", numClustersStr, ".titan.txt")
}
if (is.null(outparam)){
	outparam <- gsub(".titan.txt", ".param.txt", outfile)
}
if (is.null(outseg)){
	outseg <- gsub(".titan.txt", ".segs.txt", outfile)
}
if (is.null(outigv)){
	outigv <- gsub(".titan.txt", ".seg", outfile)
}
if (is.null(outplot)){
	outplot <- paste0(outdir, "/", id, "_cluster", numClustersStr, "/")
	dir.create(outplot)
}
outImage <- gsub(".titan.txt", ".RData", outfile)

## set up chromosome naming convention ##
seqinfo.hg19 <- readRDS("/home/unix/gavinha/software/code/git/scripts/references/Seqinfo_hg19.rds")
chrs <- setGenomeStyle(chrs, genomeStyle = genomeStyle)

pseudo_counts <- 1e-300
centromereFlank <- 100000
maxI <- 50

message('Running TITAN...')
save.image(file=outImage)
#### LOAD DATA ####
#data <- loadAlleleCounts(hetfile, header=T, genomeStyle = genomeStyle)
data <- loadHaplotypeAlleleCounts(hetfile, seqinfo = seqinfo.hg19, haplotypeBinSize = haplotypeBinSize, fun=phaseSummarizeFun) 
data <- data$haplotypeData

#### REMOVE CENTROMERES ####
if (!is.null(centromere)){
	centromere <- read.delim(centromere,header=T,stringsAsFactors=F,sep="\t")
}

#### LOAD GC AND MAPPABILITY CORRECTED COVERAGE LOG RATIO FILE ####
message('titan: Loading GC content and mappability corrected log2 ratios...')
cnData <- fread(cnfile)
cnData$chr <- setGenomeStyle(cnData$chr, genomeStyle = genomeStyle)

#### ADD CORRECTED LOG RATIOS TO DATA OBJECT ####
message('titan: Extracting read depth...')
logR <- getPositionOverlap(data$chr,data$posn, cnData)
data$logR <- log(2^logR)
rm(logR,cnData)

#### FILTER DATA FOR DEPTH, MAPPABILITY, NA, etc ####
if (!is.null(mapWig)){
	mScore <- as.data.frame(wigToRangedData(mapWig))
	mScore <- getPositionOverlap(data$chr,data$posn,mScore[,-4])
	data <- filterData(data,chrs,minDepth=minDepth,maxDepth=maxDepth,map=mScore,mapThres=mapThres, centromeres = centromere)
	rm(mScore)
}else{
	data <- filterData(data,chrs,minDepth=minDepth,maxDepth=maxDepth,centromeres = centromere)
}
## reassign chromosomes ##
chrs <- unique(data$chr)

#### LOAD PARAMETERS ####
message('titan: Loading default parameters')
params <- loadDefaultParameters(copyNumber=maxCN,numberClonalClusters=numClusters, 
																skew=skew, hetBaselineSkew=hetBaselineSkew, 
																alleleEmissionModel = alleleEmissionModel, data=data)

#### MODEL SELECTION USING EM (FWD-BACK) TO SELECT NUMBER OF CLUSTERS ####
registerDoMC()
options(cores=numCores)
message("Using ",getDoParWorkers()," cores.")
K <- length(params$genotypeParams$rt)
params$genotypeParams$alphaKHyper <- rep(alphaK,K)
params$genotypeParams$betaKHyper <- rep(25,K)
#params$genotypeParams$alphaKHyper[params$genotypeParams$ct == 0] <- alphaK / 500 #HOMD
#params$genotypeParams$alphaKHyper[params$genotypeParams$ct == max(params$genotypeParams$ct)] <- alphaK / 50 #maxCN
#params$genotypeParams$var_0[params$genotypeParams$ct %in% c(2, 4, 8)] <- 1/20 / 100
#params$genotypeParams$alphaKHyper[params$genotypeParams$ct %in% c(1,2,3)] <- alphaK / 2
params$genotypeParams$alphaRHyper <- rep(alphaR,K)
params$genotypeParams$betaRHyper <- rep(25,K)
#params$genotypeParams$alphaRHyper[params$genotypeParams$ct == 0] <- alphaR / 500 #HOMD
#params$genotypeParams$alphaRHyper[params$genotypeParams$ct == max(params$genotypeParams$ct)] <- alphaR / 50 #maxCN
#params$genotypeParams$varR_0[c(4, 9, 16, 25)] <- 1/20 / 100
#params$genotypeParams$alphaRHyper[c(4,9,25)] <- alphaR * 5
#params$genotypeParams$alphaKHyper[c(1,7:K)] <- alphaHigh 
params$ploidyParams$phi_0 <- ploidy_0
params$normalParams$n_0 <- norm_0
#params$genotypeParams$rt[c(4, 9)] <- hetAR

convergeParams <- runEMclonalCN(data, params=params,
                                #gParams=params$genotypeParams,nParams=params$normalParams,
                                #pParams=params$ploidyParams,sParams=params$cellPrevParams,
                                maxiter=maxI,maxiterUpdate=15,
                                txnExpLen=txn_exp_len,txnZstrength=txn_z_strength,
                                useOutlierState=FALSE,
                                normalEstimateMethod=normEstMeth,estimateS=estimateS,
                                estimatePloidy=boolEstPloidy, pseudoCounts=pseudo_counts)
    
#### COMPUTE OPTIMAL STATE PATH USING VITERBI ####
message("Using ",getDoParWorkers()," cores.")
optimalPath <- viterbiClonalCN(data,convergeParams)
save.image(file=outImage)
#### PRINT RESULTS TO FILES ####
results <- outputTitanResults(data,convergeParams,optimalPath,is.haplotypeData=TRUE,
			filename=NULL,posteriorProbs=FALSE,subcloneProfiles=TRUE, 
			proportionThreshold = 0.05, proportionThresholdClonal = 0.05, verbose=FALSE)
#corrResults <- removeEmptyClusters(convergeParams, results, proportionThreshold = 0.05, proportionThresholdClonal = 0.05)
convergeParams <- results$convergeParams
results <- results$corrResults
numClustersToPlot <- nrow(convergeParams$s)
write.table(results, file = outfile, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
outputModelParameters(convergeParams, results, outparam)

# save specific objects to a file
convergeParams$rhoG <- NULL; convergeParams$rhoZ <- NULL
#save(convergeParams, file=outImage)
save.image(file=outImage)
#### OUTPUT SEGMENTS ####
message("Writing segments to ", outseg)
segs <- outputTitanSegments(results, id, convergeParams, filename = outseg, igvfilename = outigv)

#### PLOT RESULTS ####
dir.create(outplot)
norm <- tail(convergeParams$n,1)
ploidy <- tail(convergeParams$phi,1)
for (chr in chrs){
	outfig <- paste0(outplot, "/", id, "_cluster", numClustersStr, "_chr", chr, ".png")
	#if (as.numeric(numClusters) <= 2){
	#  png(outfig,width=1200,height=1250,res=100)
	#	par(mfrow=c(5,1))
	#}else{
	  png(outfig,width=1200,height=1000,res=100)
		par(mfrow=c(4,1))  
	#}
	#source("/home/unix/gavinha/software/code/git/scripts/10XGenomics/utils.R")
	#source("/Volumes/software/code/git/scripts/10XGenomics/utils.R")
	source(paste0(libdir, "R/plotting.R"))
	plotCNlogRByChr(results, chr, segs = NULL, ploidy=ploidy, normal = norm, geneAnnot=NULL,  cex.axis=1.5, 
					cex.lab=1.5, ylim=plotYlim, cex=0.5, xlab="", main=paste("Chr ",chr,sep=""))
	plotAllelicRatio(results, chr, geneAnnot=NULL, spacing=4, cex.axis=1.5, cex.lab=1.5, 
					ylim=c(0,1), xlab="", cex=0.5)
	#plotHaplotypeFraction(data, type = "AllelicRatio", colType = "haplotypes", xlab="", chr="6", cex=0.25)
	#plotHaplotypeFraction(results, chr, type = "AllelicRatio", colType = "haplotype", 
	#  xlab="", cex=0.5, cex.axis=1.5, cex.lab=1.5)
	plotHaplotypeFraction(results, chr, resultType = "HaplotypeRatio", colType = "Haplotypes", 
	  xlab="", ylim=c(0,1), cex=0.5, cex.axis=1.5, cex.lab=1.5)
	#plotHaplotypeFraction(results, chr, type = "HaplotypeRatio", colType = "titan", 
	#  xlab="", cex=0.5, cex.axis=1.5, cex.lab=1.5)
	#plotHaplotypeFraction(data, chr, type = "HaplotypeRatio", colType = "haplotype", xlab="", ylim=c(0,1), cex=0.5, cex.axis=1.5, cex.lab=1.5)
	plotClonalFrequency(results, chr, normal=norm, geneAnnot=NULL, spacing=4, 
					cex.axis=1.5, ylim=c(0,1), xlab="", cex=0.5, main=paste("Chr ",chr,sep=""))
  
  par(xpd = NA)
	#if (as.numeric(numClustersToPlot) <= 2 && as.numeric(numClusters) <= 2){
		#plotSubcloneProfiles(results, chr, cex = 2, spacing=6, main=paste("Chr ",chr,sep=""), cex.axis=1.5)
		#pI <- plotIdiogram(chr, build="hg19", unit="bp", label.y=-4.25, new=FALSE, ylim=c(-2,-1))
	#}else{
		pI <- plotIdiogram(chr, build="hg19", unit="bp", label.y=-0.35, new=FALSE, ylim=c(-0.2,-0.1))
	#}
	
	dev.off()
}

################################################
############## GENOME WIDE PLOTS ###############
################################################
outFile <- paste0(outplot, "/", id, "_cluster", numClustersStr, "_CNA.pdf")
#png(outFile,width=1000,height=300)
pdf(outFile,width=20,height=6)
plotCNlogRByChr(dataIn=results, chr=NULL, segs = segs, ploidy=ploidy,  normal = norm, geneAnnot=genes, spacing=4, main=id, xlab="", ylim=plotYlim, cex=0.5, cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
dev.off()

outFile <- paste0(outplot, "/", id, "_cluster", numClustersStr, "_LOH.pdf")
#png(outFile,width=1000,height=300)
pdf(outFile,width=20,height=6)
plotAllelicRatio(dataIn=results, chr=NULL, geneAnnot=genes, spacing=4, main=id, xlab="", ylim=c(0,1), cex=0.5, cex.axis=1.5, cex.lab=1.5, cex.main=1.5)	
dev.off()

outFile <- paste0(outplot, "/", id, "_cluster", numClustersStr, "_PHASE.pdf")
#png(outFile,width=1000,height=300)
pdf(outFile,width=20,height=6)
par(mfrow=c(1,1))
plotHaplotypeFraction(results, chr=NULL, resultType = "HaplotypeRatio", colType = "Haplotypes", xlab="", ylim=c(0,1), cex=0.5, cex.axis=1.5, cex.lab=1.5)
dev.off()

outFile <- paste0(outplot, "/", id, "_cluster", numClustersStr, "_HAPLO.pdf")
#png(outFile,width=1000,height=300)
pdf(outFile,width=20,height=6)
par(mfrow=c(1,1))
plotHaplotypeFraction(results, chr=NULL, resultType = "HaplotypeRatio", colType = "CopyNumber", xlab="", ylim=c(0,1), cex=0.5, cex.axis=1.5, cex.lab=1.5)
dev.off()

outFile <- paste0(outplot, "/", id, "_cluster", numClustersStr, "_CF.pdf")
#png(outFile,width=1000,height=300)
pdf(outFile,width=20,height=6)
plotClonalFrequency(dataIn=results, chr=NULL, norm, geneAnnot=genes, spacing=4, main=id, xlab="", ylim=c(0,1), cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
dev.off()

outFile <- paste0(outplot, "/", id, "_cluster", numClustersStr, "_LOH-SEG.pdf")
pdf(outFile, width=20, height=6)
plotSegmentMedians(dataIn=segs, chr=NULL, resultType = "AllelicRatio", plotType = "CopyNumber", plot.new=T, ylim=c(0,8), cex.axis=1.5, cex.lab=1.5, cex.main=1.5)



if (as.numeric(numClusters) <= 2){
	outFile <- paste0(outplot, "/", id, "_cluster", numClustersStr, "_subclone.pdf")
	#png(outFile,width=1000,height=300)
	pdf(outFile,width=20,height=6)
	plotSubcloneProfiles(dataIn=results, chr=NULL, cex = 0.5, spacing=4, main=id, cex.axis=1.5, xlab="")
	dev.off()
}

