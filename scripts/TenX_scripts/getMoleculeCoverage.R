# file:   getMoleculeCoverage.R
# author: Gavin Ha, Ph.D.
#               Dana-Farber Cancer Institute
#               Broad Institute
# contact: <gavinha@broadinstitute.org>
# ULP-WGS website: http://www.broadinstitute.org/~gavinha/ULP-WGS/
# HMMcopy website: http://compbio.bccrc.ca/software/hmmcopy/ and https://www.bioconductor.org/packages/release/bioc/html/HMMcopy.html
# date:   September 15, 2016
# description: Hidden Markov model (HMM) to analyze molecule-based coverage in non-overlapping bins, genome-wide.
# This script is the main script to run the HMM.
library(optparse)
options(stringsAsFactors=FALSE)
options(bitmapType='cairo')

option_list <- list(
	make_option(c("-l", "--libdir"), type = "character", help = "Script library path"),
	make_option(c("-d", "--datadir"), type = "character", help = "Data library path containing gc/map wig files"),
	make_option(c("-t", "--tumorBXDir"), type = "character", help = "Path to directory containing tumor bed files for each chromosome containing BX tags."),
	make_option(c("-n", "--normalBXDir"), type = "character", default=NULL, help = "Path to directory containing normal bed files for each chromosome containing BX tags."),
	make_option(c("--minReadsPerBX"), type="integer", default=2, help="Minimum number of reads per barcode. [Default: %default]"),
	make_option(c("--normalPanel"), type="character", default=NULL, help="Median corrected depth from panel of normals"),
	make_option(c("-e", "--exons.bed"), type = "character", default=NULL, help = "Path to bed file containing exon regions."),
	make_option(c("-i", "--id"), type = "character", help = "Sample ID."),
	make_option(c("--gcOnly"), type="logical", default=FALSE, help = "TRUE if only correct for GC content bias"),
	make_option(c("-c", "--centromere"), type="character", default=NULL, help = "File containing Centromere locations"),
	make_option(c("--rmCentromereFlankLength"), type="numeric", default=1e5, help="Length of region flanking centromere to remove"),
	make_option(c("--normal"), type="character", default="0.5", help = "Initial normal contamination"),
	make_option(c("--coverage"), type="numeric", default=NULL, help = "PICARD sequencing coverage"),
	make_option(c("--lambda"), type="character", default="100", help="Initial Student's t precision."),
#	make_option(c("--kappa"), type="character", default=50, help="Initial state distribution"),
	make_option(c("-p", "--ploidy"), type="character", default="2", help = "Initial tumour ploidy"),
	make_option(c("-m", "--maxCN"), type="numeric", default=5, help = "Gain CN states above ploidy initialization"),
	make_option(c("--estimateNormal"), type="logical", default=TRUE, help = "Estimate normal."),
	make_option(c("--estimatePloidy"), type="logical", default=FALSE, help = "Estimate tumour ploidy."),
	make_option(c("--chrNormalize"), type="character", default="c(1:22)", help = "Specify chromosomes to normalize GC/mappability biases"),
	make_option(c("--chrTrain"), type="character", default="c(1:22)", help = "Specify chromosomes to estimate params."),
	make_option(c("--chrs"), type="character", default="c(1:22,\"X\")", help = "Specify chromosomes to analyze."),
	make_option(c("--diploidChrX"), type="logical", default=FALSE, help = "Assume 2 copies of chrX"),
	make_option(c("--normalizeMaleX"), type="logical", default=FALSE, help = "If male, then normalize chrX by median"),
	make_option(c("--fracReadsInChrYForMale"), type="numeric", default=0.001, help = "Threshold for fraction of reads in chrY to assign as male."),
	make_option(c("--includeHOMD"), type="logical", default=TRUE, help="If false, then exclude HOMD state"),
	make_option(c("--txnE"), type="numeric", default=0.9999999, help = "Self-transition probability"),
	make_option(c("--txnStrength"), type="numeric", default=10000000, help = "Transition pseudo-counts"),
	make_option(c("--plotFileType"), type="character", default="pdf", help = "File format for output plots"),
	make_option(c("--plotYLim"), type="character", default="c(-2,4)", help = "ylim to use for chromosome plots"),
	make_option(c("-o", "--outDir"), type="character", default="./", help = "Output Directory")
)
parseobj <- OptionParser(option_list=option_list)
opt <- parse_args(parseobj)
print(opt)
options(scipen=0, stringsAsFactors=F)

id <- opt$id
tumour_file <- opt$tumorBXDir
normal_file <- opt$normalBXDir
normal_panel <- opt$normalPanel
minReadsPerBX <- opt$minReadsPerBX
gcOnly <- as.logical(opt$gcOnly)  # {TRUE, FALSE}
exons.bed <- opt$exons.bed  # "0" if none specified
centromere <- opt$centromere
flankLength <- opt$rmCentromereFlankLength
normal <- eval(parse(text = opt$normal))
lambda <- eval(parse(text = opt$lambda))
estimateNormal <- opt$estimateNormal
estimatePloidy <- opt$estimatePloidy
ploidy <- eval(parse(text = opt$ploidy))
coverage <- opt$coverage
maxCN <- opt$maxCN
txnE <- opt$txnE
txnStrength <- opt$txnStrength
diploidChrX <- as.logical(opt$diploidChrX)
normalizeMaleX <- as.logical(opt$normalizeMaleX)
includeHOMD <- as.logical(opt$includeHOMD)
fracReadsInChrYForMale <- opt$fracReadsInChrYForMale
outDir <- opt$outDir
libdir <- opt$libdir
datadir <- opt$datadir
plotFileType <- opt$plotFileType
plotYLim <- eval(parse(text=opt$plotYLim))
gender <- NULL

maxiter <- 50
diploidChrXThreshold <- 0.5
chrs <- eval(parse(text = opt$chrs))
#chrs = paste0("chr", c(1:22, "X"))
chrsAll <- c(1:22, "X")
chrTrain <- eval(parse(text=opt$chrTrain))
#chrTrain <- paste0("chr", chrTrain)
chrNormalize <- eval(parse(text=opt$chrNormalize))
dir.create(paste0(outDir,"/",id))
outRdata <- paste0(outDir,"/",id,".RData")
save.image(file=outRdata)


library(TitanCNA)
library(ichorCNA)
library(data.table)
library(quantsmooth)

libdir <- opt$libdir
if (!is.null(libdir)){
  source(paste0(libdir, "R/haplotype.R"))
}


if (!is.null(centromere)){
	centromere <- read.delim(centromere,header=T,stringsAsFactors=F,sep="\t")
}
#save.image()

### LOAD TUMOUR AND NORMAL FILES ###
tumour_doc <- loadBXcountsFromBEDDir(tumour_file, minReads = minReadsPerBX)
tumour_doc$BX.medianNorm <- log2(tumour_doc$BX.count / median(tumour_doc$BX.count, na.rm=T))

## LOAD GC/MAP WIG FILES ###
# find the bin size and load corresponding wig files #
binSize <- as.data.frame(tumour_doc[1,])$width / 1000
gcWig <- list.files(libdir, pattern = paste0("gc_hg[0-9]+_",binSize,"kb.wig"), full.names=T)
mapWig <- paste0(libdir,"/map_hg19_",binSize,"kb.wig")
#message("Reading GC and mappability files")
gc <- wigToRangedData(gcWig)
map <- wigToRangedData(mapWig)

## FILTER BY EXONS IF PROVIDED ##
## add gc and map to RangedData object ##
if (!is.null(exons.bed)){
	targetedSequences <- read.delim(exons.bed, header=F, sep="\t", skip=86)
}else{
	targetedSequences <- NULL
}
tumour_counts <- loadReadCountsFromWig(tumour_doc, chrs=chrs, gc=gc, map=map, genomeStyle = "NCBI", 
				centromere=NULL, flankLength = NULL, targetedSequences=targetedSequences)
tumour_copy <- correctReadCounts(tumour_counts, chrNormalize=chrNormalize)
tumour_copy <- filterByMappabilityScore(tumour_copy, map=map, mapScoreThres = 0.9)

### CORRECT NORMAL DATA FOR GC CONTENT AND MAPPABILITY BIASES ###
if (!is.null(normal_file)){
	normal_doc <- loadBXcountsFromBEDDir(normal_file, minReads = minReadsPerBX)
	normal_doc$BX.medianNorm <- log2(normal_doc$BX.count / median(normal_doc$BX.count, na.rm=T))
	normal_counts <- loadReadCountsFromWig(normal_doc, chrs=chrs, gc=gc, map=map, genomeStyle = "NCBI", 
				centromere=NULL, flankLength = NULL, targetedSequences=targetedSequences)
	normal_copy <- correctReadCounts(normal_counts, chrNormalize=chrNormalize)
	normal_copy <- filterByMappabilityScore(normal_copy, map=map, mapScoreThres = 0.9)
}

tumour_copy$tumour.copy <- tumour_copy$copy
tumour_copy$tumor.BX.count <- tumour_copy$reads
tumour_copy$normal.copy <- normal_copy$copy
tumour_copy$normal.BX.count <- normal_copy$reads
tumour_copy$normal.BX.medianNorm <- normal_copy$BX.medianNorm
tumour_copy$tumour.copy <- tumour_copy$copy - normal_copy$copy
tumour_copy$copy <- tumour_copy$tumour.copy
save.image(file=outRdata)

### OUTPUT FILE ###
### PUTTING TOGETHER THE COLUMNS IN THE OUTPUT ###
outMat <- as.data.frame(tumour_copy)
#outMat <- outMat[, -which(colnames(outMat) %in% c("width", "copy"))]
outMat <- outMat[, c("space","start","end","tumour.copy")]
colnames(outMat)[1:3] <- c("chr","start","end")
outFile <- paste0(outDir,"/",id,".BXcounts.txt")
message(paste("Outputting to:", outFile))
write.table(outMat, file=outFile, row.names=F, col.names=T, quote=F, sep="\t")

##############################################
### PERFORM SEGMENTATION USING 6-STATE HMM ###
##############################################
chr <- space(tumour_copy)
chrInd <- chr %in% chrTrain
message("Segmenting")
#valid <- tumour_copy$valid# & autosomes
valid <- !is.na(tumour_copy$copy) & is.finite(tumour_copy$copy)
#valid <- valid & tumour_copy$BX.ratio > -4 & tumour_copy$BX.ratio < 5

results <- list()
loglik <- matrix(NA, nrow = length(normal) * length(ploidy), ncol = 4, 
	dimnames = list(c(), c("init", "n_est", "phi_est", "loglik")))
counter <- 1
compNames <- rep(NA, nrow(loglik))
#### restart for purity and ploidy values ####
for (n in normal){
	for (p in ploidy){
		## uncorrected depth ##
		#param <- getDefaultParameters(log(tumour_copy$reads[valid & chrInd] / median(tumour_copy$reads[valid], na.rm = TRUE)))
		#param$n_0 <- n; param$phi_0 <- p; param$lambda <- lambda
		#hmmResults.raw <- HMMsegment(tumour_copy[valid, ], id = id, dataType = "reads", 
		#		param = param, chrTrain = chrTrain, maxiter = 50, 
		#		estimateNormal = estimateNormal, estimatePloidy = estimatePloidy, logTransform = TRUE, verbose = TRUE) 

		## corrected depth ##
		maxCNum <- maxCN + floor(p)
		#maxCNum <- floor(p) + (maxCN * floor(p / 2) + ceiling(p %% 2))
		if (includeHOMD){
			param <- getDefaultParameters(log(2^tumour_copy$copy[valid & chrInd]), 
							ct=0:maxCNum, ploidy = floor(p), e=txnE, strength=txnStrength)
			K <- length(param$ct)
			param$A[1, 2:K] <- param$A[1, 2:K] * 1e-5; param$A[2:K, 1] <- param$A[2:K, 1] * 1e-5;
			param$A[1, 1] <- param$A[1, 1] * 1e-5
			param$A <- normalize(param$A); param$dirPrior <- param$A * param$strength
		}else{
			param <- getDefaultParameters(log(2^tumour_copy$copy[valid & chrInd]), 
							ct=1:maxCNum, ploidy = floor(p), e=txnE, strength=txnStrength)
		}
		param$lambda[param$ct >= 4] <- lambda
		param$lambda[param$ct == max(param$ct)] <- lambda / 10
		param$lambda[param$ct %in% c(1,2,3)] <- lambda
		param$lambda[param$ct %in% c(0)] <- lambda / 2 		
		param$n_0 <- n; param$phi_0 <- p; #param$lambda <- lambda
		hmmResults.cor <- HMMsegment(tumour_copy[valid, ], id = id, dataType = "copy", 
				param = param, chrTrain = chrTrain, maxiter = maxiter,
				estimateNormal = estimateNormal, estimatePloidy = estimatePloidy, verbose = TRUE) 
    
   	 ## convert full diploid solution to have 1.0 normal or 0.0 purity
		altInd <- hmmResults.cor$cna[hmmResults.cor$cna[, "chr"] %in% chrTrain, "event"] == "NEUT"
		if (sum(!altInd, na.rm=TRUE) / sum(altInd, na.rm = TRUE) < 1e-3){
			hmmResults.cor$results$n <- 1.0
		}
    
    ## plot genome wide figures for each solution ##
		ploidyS <- tail(hmmResults.cor$results$phi, 1)
		normEstS <- tail(hmmResults.cor$results$n, 1)
		purityS <- 1 - normEstS 
		ploidyAllS <- (1 - normEstS) * ploidyS + normEstS * 2

		outPlotFile <- paste0(outDir,"/",id,"/",id,"_genomeWide_", "n", n, "-p", p)
		if (plotFileType == "png"){ 
			outPlotFile <- paste0(outPlotFile, ".png")
			png(outPlotFile,width=20,height=6,units="in",res=350)
		}else{
			outPlotFile <- paste0(outPlotFile, ".pdf")
			pdf(outPlotFile,width=20,height=6)
		}		
		plotCNlogRByChr(hmmResults.cor$cna, segs = hmmResults.cor$results$segs, chr=NULL, 
			ploidy = ploidyAllS, cytoBand=T, yrange=plotYLim)
		mtext(line=-1, paste0("Tumor Fraction: ", format(purityS, digits=2), ", Ploidy: ", format(ploidyS, digits=3)), cex=1.5)
		dev.off()
    
  	results[[counter]] <- hmmResults.cor
  	loglik[counter, "loglik"] <- format(tail(hmmResults.cor$results$loglik, 1), digits = 4)
  	loglik[counter, "init"] <- paste0("n", n, "-p", p)
  	loglik[counter, "n_est"] <- format(tail(hmmResults.cor$results$n, 1), digits = 2)
  	loglik[counter, "phi_est"] <- format(tail(hmmResults.cor$results$phi, 1), digits = 4)
    counter <- counter + 1
	}
}

### SAVE R IMAGE ###
save.image(file=outRdata)
#save(tumour_copy, results, loglik, file=paste0(outDir,"/",id,".RData"))
ind <- order(as.numeric(loglik[, "loglik"]), decreasing=TRUE) 
hmmResults.cor <- results[[ind[1]]]
hmmResults.cor$results$loglik <- as.data.frame(loglik)

outputHMM(hmmResults.cor$cna, hmmResults.cor$results, outDir=outDir)

### OUTPUT ESTIMATED PARAMETERS ###
#outFile <- paste0(outDir, "/", id, "_unCorrected_params.txt")
#outputParametersToFile(hmmResults.raw$results, file = outFile)
outFile <- paste0(outDir, "/", id, ".params.txt")
outputParametersToFile(hmmResults.cor, file = outFile)

### PLOT THE LOG RATIO DATA ALONG WITH COLOUR-CODING FOR PREDICTED CALLS ###
for (i in chrs){
    ## PLOT CNA BY CHROMOSOME ##
     outPlot <- paste0(outDir,"/",id,"/",id,"_CNA_chr",i)
    if (plotFileType == "png"){ 
			outPlot <- paste0(outPlot, ".png")
			png(outPlot,width=15,height=4,units="in",res=300)
		}else{
			outPlot <- paste0(outPlot, ".pdf")
			pdf(outPlot,width=15,height=4)
		}			
  	par(mfrow=c(1,1))
  	#plotCNlogRByChr(hmmResults.raw$cna[,1:13], 
  	#		segs=hmmResults.raw$segs$results, chr=i, 
  	#		cytoBand=T, yrange=c(-4,6))	
  	plotCNlogRByChr(hmmResults.cor$cna, cex = 0.25,
  			segs=hmmResults.cor$segs$results, chr=i, 
  			cytoBand=T, yrange=plotYLim)	
    dev.off()

	### PLOT THE CORRECTION COMPARISONS ###
	outPlot <- paste0(outDir,"/",id,"/",id,"_correct_chr",i)
	if (plotFileType == "png"){ 
		outPlot <- paste0(outPlot, ".png")
		png(outPlot,width=20,height=12,units="in",res=300)
	}else{
		outPlot <- paste0(outPlot, ".pdf")
		pdf(outPlot,width=10,height=12)
	}		
	plotCorrection(tumour_copy, chr=i, pch = ".")
	dev.off()
}

### PLOT CNA LANDSCAPE COMPARISON OF RAW VS CORRECTED ##
#outPlotFile <- paste0(outDir,"/",id,"/",id,"_genomeWide_raw-vs-cor")
#if (plotFileType == "png"){ 
#	outPlotFile <- paste0(outPlotFile, ".png")
#	png(outPlotFile,width=20,height=6,units="in",res=300)
#}else{
#	outPlotFile <- paste0(outPlotFile, ".pdf")
#	pdf(outPlotFile,width=20,height=6)
#}
#pdf(outPlot,width=20,height=16)
#par(mfrow=c(2,1))
#plotCNlogRByChr(hmmResults.raw$cna[,1:13], chr=NULL, cytoBand=F, yrange=c(-2,2))
#plotCNlogRByChr(hmmResults.cor$cna[,1:13], chr=NULL, cytoBand=T, yrange=c(-2,2))
#dev.off()

### PLOT FOR WHOLE GENOME ###
ploidy <- tail(hmmResults.cor$results$phi, 1)
normEst <- tail(hmmResults.cor$results$n, 1)
purity <- 1 - normEst 
ploidyAll <- (1 - normEst) * ploidy + normEst * 2

outPlotFile <- paste0(outDir,"/",id,"/",id,"_genomeWide")
if (plotFileType == "png"){ 
	outPlotFile <- paste0(outPlotFile, ".png")
	png(outPlotFile,width=20,height=6,units="in",res=300)
}else{
	outPlotFile <- paste0(outPlotFile, ".pdf")
	pdf(outPlotFile,width=20,height=6)
}
plotCNlogRByChr(hmmResults.cor$cna, segs = hmmResults.cor$results$segs, chr=NULL, cex = 0.25,
	ploidy = ploidyAll, cytoBand=T, yrange=plotYLim)
annotStr <- paste0("Tumor Fraction: ", format(purity,digits=2), 
											", Ploidy: ", format(ploidy,digits=3))
if (!is.null(coverage)){
	annotStr <- paste0(annotStr, ", Coverage: ", signif(coverage,digits=2))
}
mtext(line=-1, annotStr, cex=1.5)
dev.off()

### PLOT THE CORRECTION COMPARISONS ###
outPlotFile <- paste0(outDir,"/",id,"/",id,"_correct")
if (plotFileType == "png"){ 
	outPlotFile <- paste0(outPlotFile, ".png")
	png(outPlotFile,width=20,height=18,units="in",res=300)
}else{
	outPlotFile <- paste0(outPlotFile, ".pdf")
	pdf(outPlotFile,width=20,height=18)
}		
plotCorrectionGenomeWide(tumour_copy, pch = ".")
dev.off()

### PLOT THE BIAS ###
outPlotFile <- paste0(outDir,"/",id,"/",id,"_bias.pdf")
pdf(outPlotFile)
plotBias(tumour_copy, pch = 20, cex = 0.5)
dev.off()

### PLOT TPDF ##
outPlotFile <- paste0(outDir,"/",id,"/",id,"_tpdf.pdf")
pdf(outPlotFile)
plotParam(hmmResults.cor$results, copy.states = 1:maxCNum)
dev.off()


