#!/usr/bin/env Rscript

#' titanCNA.R
#' author: Gavin Ha
#' Dana-Farber Cancer Institute
#' Broad Institute
#' contact: <gavinha@gmail.com> or <gavinha@broadinstitute.org>
#' date:	  May 14, 2017
#' Notes: This script is tested for TitanCNA v1.13.1 and higher

suppressPackageStartupMessages(
    require(optparse, quietly=TRUE, warn.conflicts=FALSE)
)
suppressPackageStartupMessages(
    require(TitanCNA, quietly=TRUE, warn.conflicts=FALSE)
)
suppressPackageStartupMessages(
    require(data.table, quietly=TRUE, warn.conflicts=FALSE)
)
suppressPackageStartupMessages(
    require(GenomicRanges, quietly=TRUE, warn.conflicts=FALSE)
)
suppressPackageStartupMessages(
    require(dplyr, quietly=TRUE, warn.conflicts=FALSE)
)
suppressPackageStartupMessages(
    require(doMC, quietly=TRUE, warn.conflicts=FALSE)
)
suppressPackageStartupMessages(
    require(SNPchip, quietly=TRUE, warn.conflicts=FALSE)
)

option_list <- list(
	make_option(c("--id"), type = "character", help = "Sample ID"),
	make_option(c("--hetFile"), type = "character",
            help = "File containing allelic read counts at HET sites. [Required]"),
	make_option(c("--cnFile"), type = "character",
            help = "File containing normalized coverage as log2 ratios. [Required]"),
	make_option(c("--outDir"), type = "character",
            help = "Output directory to output the results. [Required]"),
	make_option(c("--numClusters"), type = "integer", default = 1,
            help = "Number of clonal clusters. [Default: 1]"),
	make_option(c("--numCores"), type = "integer", default = 1,
            help = "Number of cores to use. [Default: %default]"),
	make_option(c("--ploidy_0"), type = "numeric", default = 2,
            help = "Initial ploidy value; float [Default: %default]"),
	make_option(c("--estimatePloidy"), type = "logical", default = TRUE,
            help = "Estimate ploidy; TRUE or FALSE [Default: %default]"),
	make_option(c("--normal_0"), type = "numeric", default = 0.5,
            help = "Initial normal contamination (1-purity); float [Default: %default]"),
	make_option(c("--estimateNormal"), type = "character", default = "map",
            help = "Estimate normal contamination method; string {'map', 'fixed'} [Default: %default]"),
	make_option(c("--estimateClonality"), type="logical", default=TRUE,
            help="Estimate cellular prevalence. [Default: %default]"),
	make_option(c("--maxCN"), type = "integer", default = 8,
            help = "Maximum number of copies to model; integer [Default: %default]"),
	make_option(c("--alphaK"), type = "numeric", default = 10000,
            help = "Hyperparameter on Gaussian variance; for WES, use 2500; for WGS, use 10000; float [Default: %default]"),
	make_option(c("--alphaKHigh"), type = "numeric", default = 10000,
            help = "Hyperparameter on Gaussian variance for extreme copy number states; for WES, use 2500; for WGS, use 10000; float [Default: %default]"),
	make_option(c("--txnExpLen"), type = "numeric", default = 1e15,
            help = "Expected length of segments; higher leads to longer (less sensitive) segments; float [Default: %default]"),
	make_option(c("--txnZStrength"), type = "numeric", default = 1,
            help = "Expected length of clonal cluster segmentation (factor of txnExpLen); float [Default: %default]"),
	make_option(c("--minDepth"), type = "integer", default = 10,
            help = "Minimum read depth of a HET site to include in analysis; integer [Default: %default]"),
	make_option(c("--maxDepth"), type = "integer", default = 1000,
            help = "Maximum read depth of a HET site to include in analysis; integer [Default: %default]"),
	make_option(c("--skew"), type = "numeric", default=0,
            help = "Allelic reference skew for all states except heterozygous states (e.g. 1:1, 2:2, 3:3). Value is additive to baseline allelic ratios. float [Default: %default]"),
	make_option(c("--hetBaselineSkew"), type="numeric", default=NULL,
            help="Allelic reference skew for heterozygous states (e.g. 1:1, 2:2, 3:3). Value is the additive to baseline allelic ratios. float [Default: %default]"),
	make_option(c("--minClustProportion"), type="numeric", default=0.05,
            help="Minimum proportion of the genome altered (by SNPs) for a cluster to be retained.  Clonal clusters having lower proportion of alteration are removed. [Default: %default]"),
	make_option(c("--genomeStyle"), type = "character", default = "NCBI",
            help = "NCBI or UCSC chromosome naming convention; use UCSC if desired output is to have \"chr\" string. [Default: %default]"),
	make_option(c("--genomeBuild"), type = "character", default = "hg38", help="Genome build to use; will load Seqinfo from GenomeInfoDb."),
	make_option(c("--chrs"), type = "character", default = "c(1:22, 'X')",
            help = "Chromosomes to analyze; string [Default: %default"),
    make_option(c("--gender"), type = "character", default = "male", help = "User specified gender: male or female [Default: %default]"),
	make_option(c("--mapWig"), type = "character", default = NULL,
            help = "Mappability score file for bin sizes matching cnfile. [Default: %default]"),
	make_option(c("--mapThres"), type = "numeric", default = 0.9,
            help = "Minimum mappability score threshold to use; float [Default: %default]"),
	make_option(c("--centromere"), type = "character", default=NULL,
            help = "Centromere gap file. [Default: %default]"),
    make_option(c("--cytobandFile"), type = "character", default = "None", help = "Cytoband file should be provided only if reference genome is hg38 and genomeStyle is UCSC."),
	make_option(c("--libdir"), type = "character", default=NULL,
            help = "Directory containing source code. Specify if changes have been made to source code and want to over-ride package code. [Default: %default]"),
	make_option(c("--outFile"), type = "character", default = NULL,
            help = "Output file to write position-level file. (default uses extension: *.titan.txt]"),
	make_option(c("--outSeg"), type = "character", default = NULL,
            help = "Output file to write detailed segments. (default uses extension: *.segs.txt]"),
	make_option(c("--outIGV"), type = "character", default = NULL,
            help = "Output file to write segments for loading into IGV. (default uses extension: *.seg]"),
	make_option(c("--outParam"), type = "character", default = NULL,
            help = "Output file to write parameters. [Default: %default]"),
	make_option(c("--outPlotDir"), type = "character", default = NULL,
            help = "Output directory to save plots. [Default: %default]"),
	make_option(c("--plotYlim"), type = "character", default = "c(-2,4)",
            help = "The Y-axis limits to use for plotting log ratio coverage results. [Default: %default]"),
	make_option(c("--verbose"), type = "logical", default = FALSE,
            help = "Be verbose; TRUE or FALSE [Default: %default]")
)

parseobj <- OptionParser(option_list=option_list, usage = "usage: Rscript %prog [options]")
opt <- parse_args(parseobj)
options(bitmapType='cairo', scipen=0)

libdir <- opt$libdir
if (!is.null(libdir) & libdir != "None"){
  source(paste0(libdir, "/R/plotting.R"))
  source(paste0(libdir, "/R/utils.R"))
  source(paste0(libdir, "/R/hmmClonal.R"))
  source(paste0(libdir, "/R/paramEstimation.R"))
  source(paste0(libdir, "/R/correction.R"))
}

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
txn_exp_len <- opt$txnExpLen
txn_z_strength <- opt$txnZStrength
mapThres <- opt$mapThres
minDepth <- opt$minDepth
maxDepth <- opt$maxDepth
skew <- opt$skew
hetBaselineSkew <- opt$hetBaselineSkew
minClustProportion <- opt$minClustProportion
chrs <- as.character(eval(parse(text = opt$chrs)))
gender <- opt$gender
genomeStyle <- opt$genomeStyle
genomeBuild <- opt$genomeBuild
cytobandFile <- opt$cytobandFile
mapWig <- opt$mapWig
centromere <- opt$centromere
outdir <- opt$outDir
outfile <- opt$outFile
outparam <- opt$outParam
outseg <- opt$outSeg
outigv <- opt$outIGV
outplot <- opt$outPlotDir
plotYlim <- eval(parse(text = opt$plotYlim))
verbose <- opt$verbose

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
	outparam <- gsub(".titan.txt", ".params.txt", outfile)
}
if (is.null(outseg)){
	outseg <- gsub(".titan.txt", ".segs.txt", outfile)
}
if (is.null(outigv)){
	outigv <- gsub(".titan.txt", ".seg", outfile)
}
if (is.null(outplot)){
	outplot <- paste0(outdir, "/", id, "_cluster", numClustersStr, "/")
}
dir.create(outplot, showWarnings=verbose)
outImage <- gsub(".titan.txt", ".RData", outfile)

bsg <- paste0("BSgenome.Hsapiens.UCSC.", genomeBuild)
if (!require(bsg, character.only=TRUE, quietly=TRUE, warn.conflicts=FALSE)) {
	seqinfo <- Seqinfo(genome=genomeBuild)
} else {
	seqinfo <- seqinfo(get(bsg))
}
seqlevelsStyle(chrs) <- genomeStyle
## exclude chrX if gender==male ##
if (gender == "male" || gender == "Male" || gender == "MALE"){
	chrs <- chrs[chrs!=grep("X", chrs, value=TRUE)]
}

pseudo_counts <- 1e-300
centromereFlank <- 100000
maxI <- 50

message('Running TITAN...')

#### LOAD DATA ####
data <- loadAlleleCounts(hetfile, header=T, genomeStyle = genomeStyle)

#### REMOVE CENTROMERES ####
if (!is.null(centromere)){
	centromere <- read.delim(centromere,header=T,stringsAsFactors=F,sep="\t")
}

#### LOAD GC AND MAPPABILITY CORRECTED COVERAGE LOG RATIO FILE ####
message('titan: Loading GC content and mappability corrected log2 ratios...')
cnData <- fread(cnfile) #read.delim(cnfile, header=T, stringsAsFactors=F, sep="\t")
cnData$chr <- setGenomeStyle(cnData$chr, genomeStyle = genomeStyle)

#### ADD CORRECTED LOG RATIOS TO DATA OBJECT ####
message('titan: Extracting read depth...')
logR <- getPositionOverlap(data$chr,data$posn,cnData)
data$logR <- log(2^logR)
rm(logR,cnData)

#### FILTER DATA FOR DEPTH, MAPPABILITY, NA, etc ####
if (!is.null(mapWig)){
	mScore <- as.data.frame(wigToRangedData(mapWig))
	mScore <- getPositionOverlap(data$chr,data$posn,mScore[,-4])
}else{
  mScore <- NULL
}
data <- filterData(data,chrs,minDepth=minDepth,maxDepth=maxDepth,
                   centromeres = centromere, centromere.flankLength = 1e6,
                   map=mScore,mapThres=mapThres)

#### LOAD PARAMETERS ####
message('titan: Loading default parameters')
params <- loadDefaultParameters(copyNumber=maxCN,numberClonalClusters=numClusters,
                                skew=skew, hetBaselineSkew=hetBaselineSkew, data=data)

#### PARAMETER ESTIMATION USING EM (FWD-BACK) TO SELECT NUMBER OF CLUSTERS ####
registerDoMC()
options(cores=numCores)
message("titan: Using ",getDoParWorkers()," cores.")
K <- length(params$genotypeParams$rt)
params$genotypeParams$alphaKHyper <- rep(alphaK,K)
params$genotypeParams$betaKHyper <- rep(25,K)
#params$genotypeParams$alphaKHyper[c(1,7:K)] <- alphaHigh
params$ploidyParams$phi_0 <- ploidy_0
params$normalParams$n_0 <- norm_0
#params$genotypeParams$rt[c(4, 9)] <- hetAR
message("titan: Parameter estimation")
convergeParams <- runEMclonalCN(data, params,
                                maxiter=maxI,maxiterUpdate=1500,
                                txnExpLen=txn_exp_len,txnZstrength=txn_z_strength,
                                useOutlierState=FALSE,
                                normalEstimateMethod=normEstMeth,estimateS=estimateS,
                                estimatePloidy=boolEstPloidy, pseudoCounts=pseudo_counts,
                                verbose=verbose)

#### COMPUTE OPTIMAL STATE PATH USING VITERBI ####
message("Optimal state path computation: Using ",getDoParWorkers()," cores.")
optimalPath <- viterbiClonalCN(data,convergeParams)
#save.image(file=outImage)
#### PRINT RESULTS TO FILES ####
results <- outputTitanResults(data,convergeParams,optimalPath,
			filename=NULL,posteriorProbs=F,subcloneProfiles=TRUE,
			proportionThreshold = minClustProportion, proportionThresholdClonal = 0.05,
			recomputeLogLik = TRUE, rerunViterbi = FALSE, verbose=verbose)
convergeParams <- results$convergeParams
results <- results$corrResults
norm <- tail(convergeParams$n,1)
ploidy <- tail(convergeParams$phi,1)

# save specific objects to a file
# if you don't specify the path, the cwd is assumed
convergeParams$rhoG <- NULL; convergeParams$rhoZ <- NULL
#save(convergeParams, file=outImage)
save.image(file=outImage)
#### OUTPUT SEGMENTS ####
segs <- outputTitanSegments(results, id, convergeParams, filename = NULL, igvfilename = outigv)
corrIntCN.results <- correctIntegerCN(results, segs, 1 - norm, ploidy, maxCNtoCorrect.autosomes = maxCN, 
		maxCNtoCorrect.X = NULL, minPurityToCorrect = 0.2, gender = gender, chrs = chrs)
results <- corrIntCN.results$cn
segs <- corrIntCN.results$segs
message("Writing results to ", outfile, ", ", outseg, ", ", outparam)
write.table(results, file = outfile, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(segs, file = outseg, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
outputModelParameters(convergeParams, results, outparam)
save.image(file=outImage)

#### PLOT RESULTS ####
numClustersToPlot <- nrow(convergeParams$s)
dir.create(outplot, showWarnings=verbose)

if (genomeBuild == "hg38" && file.exists(cytobandFile)){
	cytoband <- fread(cytobandFile)
	names(cytoband) <- c("chrom", "start", "end", "name", "gieStain")
	cytoband <- cytoband[chrom %in% chrs]
	cytoband$chrom <- setGenomeStyle(cytoband$chrom, genomeStyle = genomeStyle)
	cytoband <- as.data.frame(cytoband)
}
for (chr in unique(results$Chr)){
	chrStr <- chr
	if (genomeStyle == "NCBI"){
		chrStr <- paste0("chr", chr)
	}
	outfig <- paste0(outplot, "/", id, "_cluster", numClustersStr, "_", chrStr, ".png")
	png(outfig,width=1200,height=1200,res=100)
	if (as.numeric(numClusters) <= 2){
		par(mfrow=c(5,1))
	}else{
		par(mfrow=c(4,1))
	}
	plotCNlogRByChr(results, chr, segs = segs, ploidy=ploidy,
                        normal = norm, geneAnnot=NULL,  cex.axis=1.5,
                        ylim=plotYlim, cex=0.5, xlab="", main=paste("Chr ",chr,sep=""))
	plotAllelicRatio(results, chr, geneAnnot=NULL, spacing=4, cex.axis=1.5,
                        ylim=c(0,1), xlab="", cex=0.5, main=paste("Chr ",chr,sep=""))
	plotClonalFrequency(results, chr, normal=norm, geneAnnot=NULL, spacing=4,
                            cex.axis=1.5, ylim=c(0,1), xlab="", cex=0.5,
                            main=paste("Chr ",chr,sep=""))
    maxCorCN <- segs[chr==chr, max(Corrected_Copy_Number, na.rm = TRUE)]
	plotSegmentMedians(segs, chr=chr, resultType = "LogRatio", plotType = "CopyNumber", 
				plot.new=TRUE, ylim=c(0,maxCorCN), xlab="", spacing=4, main=paste("Chr ",chr,sep=""))
	if (as.numeric(numClustersToPlot) <= 2 && as.numeric(numClusters) <= 2){
		plotSubcloneProfiles(results, chr, cex = 2, spacing=6, xlab="",
                                     main=paste("Chr ",chr,sep=""), cex.axis=1.5)
		ylim <- c(-2,-1)
		label.y=-4.25
	}else{
		ylim <- c(-0.2,-0.1)
		label.y <- -0.35
	}
	if (genomeBuild == "hg38" && file.exists(cytobandFile)){
		sl <- seqlengths(seqinfo[chr])
  		pI <- plotIdiogram.hg38(chr, cytoband=cytoband, seqinfo=seqinfo, xlim=c(0, max(sl)), unit="bp", label.y=label.y, new=FALSE, ylim=ylim)	
  	}else{
		pI <- plotIdiogram(chr, build="hg19", unit="bp", label.y=label.y, new=FALSE, ylim=ylim)
  	}

	dev.off()
}

################################################
############## GENOME WIDE PLOTS ###############
################################################
outFile <- paste0(outplot, "/", id, "_cluster", numClustersStr, "_CNA.pdf")
#png(outFile,width=1000,height=300)
pdf(outFile,width=20,height=6)
plotCNlogRByChr(dataIn=results, chr=chrs, segs = segs, ploidy=ploidy,
                normal = norm, geneAnnot=NULL, spacing=4, main=id, xlab="",
                ylim=plotYlim, cex=0.5, cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
dev.off()

outFile <- paste0(outplot, "/", id, "_cluster", numClustersStr, "_CNASEG.pdf")
#png(outFile,width=1000,height=300)
pdf(outFile,width=20,height=6)
maxCorCN <- max(segs$Corrected_Copy_Number, na.rm = TRUE)
plotSegmentMedians(dataIn=segs, chr=chrs, resultType = "LogRatio", plotType = "CopyNumber", 
				plot.new=T, ylim=c(0,maxCorCN), cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
dev.off()

outFile <- paste0(outplot, "/", id, "_cluster", numClustersStr, "_LOH.pdf")
#png(outFile,width=1000,height=300)
pdf(outFile,width=20,height=6)
plotAllelicRatio(dataIn=results, chr=chrs, geneAnnot=NULL, spacing=4,
                 main=id, xlab="", ylim=c(0,1), cex=0.5, cex.axis=1.5,
                 cex.lab=1.5, cex.main=1.5)
dev.off()

outFile <- paste0(outplot, "/", id, "_cluster", numClustersStr, "_LOHSEG.pdf")
#png(outFile,width=1000,height=300)
pdf(outFile,width=20,height=6)
maxCorCN <- max(segs$Corrected_Copy_Number, na.rm = TRUE)
plotSegmentMedians(dataIn=segs, chr=chrs, resultType = "AllelicRatio", plotType = "CopyNumber", 
				plot.new=T, ylim=c(0,maxCorCN), cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
dev.off()

outFile <- paste0(outplot, "/", id, "_cluster", numClustersStr, "_CF.pdf")
#png(outFile,width=1000,height=300)
pdf(outFile,width=20,height=6)
plotClonalFrequency(dataIn=results, chr=chrs, norm, geneAnnot=NULL,
                    spacing=4, main=id, xlab="", ylim=c(0,1), cex.axis=1.5,
                    cex.lab=1.5, cex.main=1.5)
dev.off()

if (as.numeric(numClusters) <= 2){
	outFile <- paste0(outplot, "/", id, "_cluster", numClustersStr, "_subclone.pdf")
	#png(outFile,width=1000,height=300)
	pdf(outFile,width=20,height=6)
	plotSubcloneProfiles(dataIn=results, chr=chrs, cex = 0.5, spacing=4,
                             main=id, cex.axis=1.5, xlab="")
	dev.off()
}

