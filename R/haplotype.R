#' author: Gavin Ha 
#' 		Dana-Farber Cancer Institute
#'		Broad Institute
#' contact: <gavinha@gmail.com> or <gavinha@broadinstitute.org>
#' date:	  January 11, 2018

#' @import data.table
#' @import GenomicRanges

loadBXcountsFromBEDDir <- function(bxDir, chrs = c(1:22, "X", "Y"), minReads = 2){
	files <- list.files(bxDir, pattern=".bed", full.names = TRUE)
	message("Loading BX counts from ", bxDir)
	bxAll <- NULL
	for (i in files){
		message(i)
		awkcmd <- paste0("awk -F \"\t\" \'{if (!$4) {print $1\"\t\"$2\"\t\"$3\"\t\"\"NA\"} else {print $1\"\t\"$2\"\t\"$3\"\t\"$4}}\' ", i)
		bxChr <- fread(awkcmd, header = FALSE, na.strings = c("NA"))
		colnames(bxChr) <- c("chr", "start", "end", "BX")
		bxChr[, BXcounts:=sapply(BX, function(x){
		if (!is.na(x)){
				bxCodes <- unlist(tstrsplit(x,","))
				bxReads <- as.numeric(tstrsplit(bxCodes, "_")[[2]])
				return(sum(bxReads >= minReads))
			}else{
				return(0)
			}
		})]
		bxChr[, BX:=NULL]
		bxAll <- rbind(bxAll, bxChr)
	}
	bxGR <- RangedData(space = bxAll[[1]], ranges = IRanges(start = bxAll[[2]], 
										 end = bxAll[[3]]), BX.count = bxAll[[4]])
	bxGR <- keepChr(bxGR, chr = chrs)
	return(bxGR)
}


loadReadCountsFromBed <- function(counts, chrs = c(1:22, "X", "Y"), gc = NULL, map = NULL, centromere = NULL, flankLength = 100000, targetedSequences = NULL, genomeStyle = "NCBI"){

	names(counts) <- setGenomeStyle(names(counts), genomeStyle)
	counts <- keepChr(counts, chrs)
	if (!is.null(gc)){ 
		names(gc) <- setGenomeStyle(names(gc), genomeStyle)
		gc <- keepChr(gc, chrs)
		counts$gc <- gc$value
	}
	if (!is.null(map)){ 
		names(map) <- setGenomeStyle(names(map), genomeStyle)
		map <- keepChr(map, chrs)
		counts$map <- map$value
	}
	colnames(counts)[1] <- c("reads")
	
	# remove centromeres
	if (!is.null(centromere)){ 
		centromere$Chr <- setGenomeStyle(centromere$Chr, genomeStyle)
		counts <- excludeCentromere(counts, centromere, flankLength = flankLength)
	}
	# keep targeted sequences
	if (!is.null(targetedSequences)){
		targetedSequences[,1] <- setGenomeStyle(targetedSequences[,1], genomeStyle)
		countsExons <- filterByTargetedSequences(counts, targetedSequences)
		counts <- counts[countsExons$ix,]
	}
	## filter 0 read bins ##
	#counts <- counts[counts$reads > 0, ]
	return(counts)
}

loadHaplotypeAlleleCounts <- function(inCounts, cnfile, fun = "sum", haplotypeBinSize = 1e5, 
      minSNPsInBin = 3, chrs = c(1:22, "X"), minNormQual = 200, 
      genomeStyle = "NCBI", sep = "\t", header = TRUE, seqinfo = NULL,
      mapWig = NULL, mapThres = 0.9, centromere = NULL, minDepth = 10, maxDepth = 1000) {
	if (is.character(inCounts)){
    ## LOAD INPUT READ COUNT DATA 
    	message("titan: Loading data and phasing information ", inCounts)
    	data <- read.delim(inCounts, header = header, stringsAsFactors = FALSE, 
        		sep = sep)
      colnames(data) <- c("chr", "posn", "refBase", "ref", "nonRefBase", "nonRef", "normQual", "genotype", "phaseSet")
      if (typeof(data[,"posn"])!="integer" || typeof(data[,"ref"])!="integer" || 
          typeof(data[,"nonRef"])!="integer" || is.null(data$genotype) || is.null(data$phaseSet)){
        stop("loadHaplotypeAlleleCounts: Input counts file format does not match required specifications.")	
      }
    }else if (is.data.frame(inCounts)){  #inCounts is a data.frame
    	data <- inCounts
  	}else{
    	stop("loadHaplotypeAlleleCounts: Must provide a filename or data.frame to inCounts")
  	}
  
	if (is.null(seqinfo)){
	seqinfo <- readRDS(system.file("extdata", "Seqinfo_hg19.rds", package = "TitanCNA"))
	}
	seqinfo <- keepStandardChromosomes(seqinfo)
	seqlevelsStyle(chrs) <- genomeStyle
	# convert to desired genomeStyle and only include autosomes, sex chromosomes
	data[, 1] <- setGenomeStyle(data[, 1], genomeStyle)
   
	## sort chromosomes
	indChr <- orderSeqlevels(as.character(data[, "chr"]), X.is.sexchrom = TRUE)
	data <- data[indChr, ]
	## sort positions within each chr
	for (x in unique(data[, "chr"])){
		ind <- which(data[, "chr"] == x)
		data[ind, ] <- data[ind[order(data[ind, "posn"], decreasing = FALSE)], ]
	}
  
  ## filter data ##
  data <- cbind(data, start = data$posn, end = data$posn, depth = data$ref + data$nonRef)
  data <- data[data$normQual >= minNormQual, ]
  # get max of allele counts
  #data$refOriginal <- data$ref
  #data$ref <- pmax(data$refOriginal, data$nonRef)
  ## use GRanges to find overlap of tiling haplotypeBinSize ##
  phasedAlleles <- getPhasedAlleleFraction(data)
  data$phasedAlleleFraction <- phasedAlleles$allele.fraction
  data$phasedCount <- phasedAlleles$phasedCount
  data$ref.symmetric <- pmax(data$ref, data$nonRef)
  data.gr <- makeGRangesFromDataFrame(data, keep.extra.columns = TRUE, seqinfo = seqinfo, ignore.strand = TRUE)
  tile.gr <- unlist(tileGenome(seqinfo, tilewidth = haplotypeBinSize))
  data.gr <- keepSeqlevels(data.gr, chrs, pruning.mode="coarse")
  tile.gr <- keepSeqlevels(tile.gr, chrs, pruning.mode="coarse")
  hits <- findOverlaps(query = data.gr, subject = tile.gr)
  data.gr$haplotypeBin <- subjectHits(hits)
  ## use data.table to process haplotype counts by blocks
  data <- as(data.gr, "data.frame")
  data.dt <- data.table(as(data.gr, "data.frame"))
  haploBinSummary <- data.dt[, list(HaplotypeFraction = mean(phasedAlleleFraction), 
      HaplotypeDepth.sum = sum(phasedCount), HaplotypeBinDepth.sum = sum(depth),
      HaplotypeDepth.mean = round(mean(phasedCount)), #round?
      #HaplotypeDepth.mean = round(mean(ref.symmetric)), 
      HaplotypeBinDepth.mean = round(mean(depth)), #round?
      SNPs = length(phasedAlleleFraction)), by = c("phaseSet", "haplotypeBin")]
  # filter bins by number of SNPs #
  haploBinSummary <- haploBinSummary[SNPs >= minSNPsInBin]
  # summary bins with multiple phaseset ID such that haplotypeBin is unique
  haploBinSummary.unique <- haploBinSummary[, list(SNPs = sum(SNPs),
    HaplotypeFraction = sum(HaplotypeFraction * SNPs) / sum(SNPs), 
    HaplotypeDepth.sum = sum(HaplotypeDepth.sum),
    #HaplotypeDepth.sum.symmetric = sum(HaplotypeDepth.sum.symmetric), 
    HaplotypeBinDepth.sum = sum(HaplotypeBinDepth.sum),
    HaplotypeDepth.mean = round(sum(HaplotypeDepth.mean * SNPs) / sum(SNPs)), #round?
    #HaplotypeDepth.mean.symmetric = round(sum(HaplotypeDepth.mean.symmetric * SNPs) / sum(SNPs)), 
    HaplotypeBinDepth.mean = round(sum(HaplotypeBinDepth.mean * SNPs) / sum(SNPs)), #round?
    phaseSet = phaseSet[which.max(SNPs)]), by = haplotypeBin]
    # get symmetric haplotype fraction
  haploBinSummary.unique[, HaplotypeFraction.symmetric := pmax(HaplotypeFraction, 1 - HaplotypeFraction)]
  haploBinSummary.unique[, HaplotypeDepth.sum.symmetric := pmax(HaplotypeDepth.sum, HaplotypeBinDepth.sum - HaplotypeDepth.sum)]
  haploBinSummary.unique[, HaplotypeDepth.mean.symmetric := pmax(HaplotypeDepth.mean, HaplotypeBinDepth.mean - HaplotypeDepth.mean)]
   # set bin column as key so that we can map back to original data
  setkey(haploBinSummary.unique, haplotypeBin) 
  # add the bin summarized values back to data.dt
  data.dt <- cbind(data.dt, select(haploBinSummary.unique[.(subjectHits(hits))], 
      haplotypeBin.aggr = haplotypeBin, HaplotypeFraction.symmetric, HaplotypeDepth.sum,
      HaplotypeDepth.sum.symmetric, HaplotypeBinDepth.sum, HaplotypeDepth.mean, 
      HaplotypeDepth.mean.symmetric, HaplotypeBinDepth.mean, 
      SNPs, phaseSet.aggr = phaseSet))
  data.dt <- na.omit(data.dt)
  data.dt[, phasedCount.haploSymmetric := {
    if (HaplotypeDepth.sum != HaplotypeDepth.sum.symmetric){
      depth - phasedCount 
    }else{
      phasedCount
    }
  }, by=1:nrow(data.dt)]
  alleleData <- select(data.dt, chr=seqnames, posn=start, 
      refOriginal=ref, nonRef=nonRef, tumDepth=depth)
  alleleData$chr <- as.character(alleleData$chr)
  alleleData$ref = pmax(alleleData$refOriginal, alleleData$nonRef)
  haplotypeData <- select(data.dt, chr=seqnames, posn=start, 
      phaseSet=phaseSet.aggr, refOriginal=ref, tumDepthOriginal = depth)
  haplotypeData$chr <- as.character(haplotypeData$chr)
  if (fun == "sum"){
    haplotypeData$ref <- data.dt[, HaplotypeDepth.sum.symmetric]
    haplotypeData$tumDepth <- data.dt[, HaplotypeBinDepth.sum]
    haplotypeData[, HaplotypeRatio := ref / tumDepth]
    #haplotypeData$haplotypeCount <- data.dt[, phasedCount]#data.dt[, HaplotypeDepth.sum]
    haplotypeData$haplotypeCount <- data.dt[, phasedCount.haploSymmetric]
  }else if (fun == "mean"){
    haplotypeData$ref <- data.dt[, HaplotypeDepth.mean.symmetric]
    haplotypeData$tumDepth <- data.dt[, HaplotypeBinDepth.mean]
    haplotypeData$HaplotypeRatio <- data.dt[, HaplotypeDepth.sum.symmetric] / data.dt[, HaplotypeBinDepth.sum]
    #haplotypeData[, HaplotypeRatio := ref / tumDepth]
    #haplotypeData$haplotypeCount <- data.dt[, phasedCount]#data.dt[, HaplotypeDepth.mean]
    haplotypeData$haplotypeCount <- data.dt[, phasedCount.haploSymmetric]
  }else if (fun == "SNP"){
    haplotypeData$ref <- data.dt[, phasedCount.haploSymmetric]
    haplotypeData$tumDepth <- data.dt[, depth]
    haplotypeData$HaplotypeRatio <- data.dt[, HaplotypeDepth.sum.symmetric] / data.dt[, HaplotypeBinDepth.sum]
    haplotypeData$haplotypeCount <- data.dt[, phasedCount.haploSymmetric]
  }
  haplotypeData$nonRef <- haplotypeData$tumDepth - haplotypeData$ref
  
  	## filtering ##
	if (!is.null(centromere)){
		centromere <- read.delim(centromere,header=T,stringsAsFactors=F,sep="\t")
	}

	#### LOAD GC AND MAPPABILITY CORRECTED COVERAGE LOG RATIO FILE ####
	message('titan: Loading GC content and mappability corrected log2 ratios...')
	cnData <- fread(cnfile)
	cnData$chr <- setGenomeStyle(cnData$chr, genomeStyle = genomeStyle)

	#### ADD CORRECTED LOG RATIOS TO DATA OBJECT ####
	message('titan: Extracting read depth...')
	logR <- getPositionOverlap(haplotypeData$chr,haplotypeData$posn, cnData)
	haplotypeData$logR <- log(2^logR)
	rm(logR,cnData)

	#### FILTER DATA FOR DEPTH, MAPPABILITY, NA, etc ####
	if (!is.null(mapWig)){
		mScore <- as.data.frame(wigToRangedData(mapWig))
		mScore <- getPositionOverlap(haplotypeData$chr,haplotypeData$posn,mScore[,-4])
		haplotypeData <- filterData(haplotypeData,chrs, minDepth=minDepth, maxDepth=maxDepth, 
			map=mScore,mapThres=mapThres, centromeres = centromere)
		rm(mScore)
	}else{
		haplotypeData <- filterData(haplotypeData,chrs,minDepth=minDepth,maxDepth=maxDepth,centromeres = centromere)
	}
  
  return(list(haplotypeData=haplotypeData, alleleData=alleleData, data=data.dt))
}

getPhasedAlleleFraction <- function(x){
  altInd <- as.logical(as.numeric(tstrsplit(x$genotype, "\\||\\/")[[1]]))
  refInd <- !altInd
  allele.fraction <- rep(NA, length(x))
  allele.fraction[refInd] <- x[refInd, "ref"] / x[refInd, "depth"]
  allele.fraction[altInd] <- x[altInd, "nonRef"] / x[altInd, "depth"]
  phasedCount <- rep(NA, length(x))
  phasedCount[refInd] <- x[refInd, "ref"]
  phasedCount[altInd] <- x[altInd, "nonRef"]
  return(list(allele.fraction = allele.fraction, phasedCount = phasedCount))
}

getHaplotypesFromVCF <- function(vcfFile, chrs = c(1:22, "X"), build = "hg19", genomeStyle = "NCBI",
                                 filterFlags = c("PASS", "10X_RESCUED_MOLECULE_HIGH_DIVERSITY"), 
                                 minQUAL = 100, minDepth = 10, minVAF = 0.25, altCountField = "AD",
                                 keepGenotypes = c("1|0", "0|1", "0/1"), snpDB = NULL){
  #require(data.table)
  message("Loading ", vcfFile)
  vcf <- readVcf(vcfFile, genome = build)
  chrName <- mapSeqlevels(seqlevels(vcf), style = genomeStyle)
  rowRanges(vcf) <- renameSeqlevels(rowRanges(vcf), na.omit(chrName))
  #keepGenotypes = c("1|0", "0|1", "0/1")
  
  ## filter vcf ##
  message("Filtering VCF ...")
  message("  by chromsomes")
  # keep specified chromosomes
  seqlevelsStyle(chrs) <- genomeStyle
  vcf <- keepSeqlevels(vcf, chrs, pruning.mode="coarse")
  # keep by filter flags
  indFILTER <- rowRanges(vcf)$FILTER %in% filterFlags
  # keep SNPs - ref and alt have length of 1 and only a single allele for ref/alt
  indSNP <- nchar(unstrsplit(CharacterList(rowRanges(vcf)$ALT), sep=",")) == 1 & 
    nchar(unstrsplit(rowRanges(vcf)$REF, sep=",")) == 1
  message("  by quality >=", minQUAL, ")")
  indQUAL <- rowRanges(vcf)$QUAL >= minQUAL
  # keep genotypes
  message("  by genotypes: ", paste0(keepGenotypes, collapse=";"))
  indGeno <- geno(vcf)$GT %in% keepGenotypes
  ind <- indFILTER & indSNP & indQUAL & indGeno
  vcf <- vcf[which(ind)]
  
  message("  by depth (>=", minDepth, ")")
  depth <- geno(vcf)$DP
  indDP <- depth >= minDepth
  message("  by VAF (>=", minVAF, ")")
  minVAF <- min(minVAF, 1 - minVAF) # symmetric from 0-0.5
  if (!is.null(geno(vcf)[[altCountField]])){
    altCounts <- unlist(lapply(geno(vcf)[[altCountField]], min))
  }else{
    stop("Alternate read counts not in AD or AO fields.")
  }
  indVAR <- (altCounts / depth) >= minVAF
  vcf <- vcf[which(indDP & indVAR)]
  rm(depth)
  
  if (!is.null(snpDB)){
    message (" by SNP VCF file ", snpDB)
    snp <- readVcf(snpDB, genome = build)
    snpInd <- nchar(unstrsplit(CharacterList(rowRanges(snp)$ALT), sep=",")) == 1 & 
      nchar(unstrsplit(rowRanges(snp)$REF, sep=",")) == 1
    snp <- snp[which(snpInd)]
    hits <- findOverlaps(query = rowRanges(vcf), subject = rowRanges(snp))
    indSNPDB <- queryHits(hits)
    vcf <- vcf[indSNPDB]
  }
  ## 
  #acounts <- do.call(rbind, geno(vcf)$AD)
  # phased snps
  #indPhased <- grepl(pattern = "\\|", geno(vcf)$GT)
  # phase sets
  #ps.rle <- rle(as.integer(geno(vcf)$PS))
  
  geno.gr <- rowRanges(vcf)	
  values(geno.gr) <- cbind(values(geno.gr), 
                           DataFrame(GT = as.character(geno(vcf)$GT), PS = as.character(geno(vcf)$PS)))
  HT <- getPhasedAllele(geno.gr)
  values(geno.gr) <- cbind(values(geno.gr), DataFrame(HT1 = HT$h1, HT2 = HT$h2))
  geno.gr <- keepSeqlevels(geno.gr, chrs, pruning.mode="coarse")
  
  return(list(vcf.filtered = vcf, geno = geno.gr))
}

# x GRanges object with GT, REF, ALT metadata columns
# returns reference or alternate allele in haplotype
getPhasedAllele <- function(x){
	h1 <- rep(NA, length(x))
	h2 <- rep(NA, length(x))
	# haplotype 1 (h1) ref allele for 1|0
	h1[x$GT == "1|0"] <- as.character(x$REF)[x$GT == "1|0"]
	h2[x$GT == "1|0"] <- as.character(unlist(x$ALT))[x$GT == "1|0"]
	# haplotype 1 (h1) alt allele for 0|1
	h1[x$GT == "0|1"] <- as.character(unlist(x$ALT))[x$GT == "0|1"]	
	h2[x$GT == "0|1"] <- as.character(x$REF)[x$GT == "0|1"]
	return(list(h1 = h1, h2 = h2))
}

plotHaplotypeFraction <- function(dataIn, chr = c(1:22), resultType = "HaplotypeRatio", colType = "Haplotypes", 
	phaseBlockCol = c("#9ad0f3", "#CC79A7"), geneAnnot = NULL, spacing = 4,  xlim = NULL,  ...) {
    if (!resultType %in% c("HaplotypeRatio", "AllelicRatio")){
      stop("plotHaplotypeFraction: resultType must be one of 'HaplotypeRatio' or 'AllelicRatio'.")
    }
    if (!colType %in% c("Haplotypes", "CopyNumber")){
      stop("plotHaplotypeFraction: plotType must be one of 'Haplotypes' or 'CopyNumber'")
    }
    
   	# use consistent chromosome naming convention
  	chr <- as.character(chr)
	seqlevelsStyle(chr) <- seqlevelsStyle(as.character(dataIn$Chr))

    lohCol.hap <- c(`0`=phaseBlockCol[1], `1`=phaseBlockCol[2])
    lohCol.titan <- c("#00FF00", "#006400", "#0000FF", "#8B0000", 
        "#006400", "#BEBEBE", "#FF0000", "#BEBEBE", 
        "#FF0000")
    names(lohCol.titan) <- c("HOMD", "DLOH", "NLOH", "GAIN", 
        "ALOH", "HET", "ASCNA", "BCNA", "UBCNA")
    dataIn <- copy(dataIn)
    colnames(dataIn)[1:2] <- c("Chr", "Position")
    #dataIn$AllelicRatio <- pmax(dataIn$refOriginal, dataIn$tumDepthOriginal - dataIn$refOriginal) / dataIn$tumDepthOriginal
    ## set alternating 0 and 1 for phaseSets
    psRuns <- rle(dataIn$PhaseSet)
    ps.id <- 1:length(psRuns$lengths) %% 2
    dataIn$phaseSet.id <- rep(ps.id, psRuns$lengths)
    #dataIn$HaplotypeCount[as.logical(dataIn$phaseSet.id)] <- dataIn$HaplotypeDepth[as.logical(dataIn$phaseSet.id)] - dataIn$HaplotypeCount[as.logical(dataIn$phaseSet.id)]
    #dataIn$HaplotypeRatio[as.logical(dataIn$phaseSet.id)] <- 1 - dataIn$HaplotypeRatio[as.logical(dataIn$phaseSet.id)]
    #dataIn$AllelicRatio <- dataIn$refOriginal / dataIn$tumDepthOriginal
    dataIn[, HaplotypeRatio.1 := HaplotypeRatio]#dataIn$HaplotypeCount / dataIn$HaplotypeDepth
    dataIn[, HaplotypeRatio.2 := 1 - HaplotypeRatio]#(dataIn$HaplotypeDepth - dataIn$HaplotypeCount) / dataIn$HaplotypeDepth
    
    if (!is.null(chr) && length(chr) == 1) {
        for (i in chr) {
            dataByChr <- dataIn[Chr == i, ]
            #dataByChr <- dataByChr[dataByChr[, "TITANcall"] != "OUT", ]
            # plot the data if (outfile!=''){
            # pdf(outfile,width=10,height=6) }
            if (colType == "Haplotypes" || !"TITANcall" %in% colnames(dataByChr)){
              colors.1 <- lohCol.hap[as.character(dataByChr$phaseSet.id)]
              colors.2 <- lohCol.hap[as.character(as.numeric(!dataByChr$phaseSet.id))]
            }else{
              colors.1 <- lohCol.titan[dataByChr[, TITANcall]]
              colors.2 <- colors.1
            }
            
            par(mar = c(spacing, 8, 2, 2))
            # par(xpd=NA)
            if (missing(xlim)) {
                xlim <- as.numeric(c(1, dataByChr[nrow(dataByChr), Position]))
            }
            if (resultType == "HaplotypeRatio"){
              plot(dataByChr[, Position], dataByChr[, HaplotypeRatio.1], 
                  col = colors.1, 
                  pch = 16, xaxt = "n", las = 1, ylab = "Haplotype Fraction", xlim = xlim, 
                  ...)
              points(dataByChr[, Position], dataByChr[, HaplotypeRatio.2], col = colors.1, pch=16, ...)
            }else if (resultType == "AllelicRatio"){
               plot(dataByChr[, Position], dataByChr[, AllelicRatio], 
                  col = colors.1, 
                  pch = 16, xaxt = "n", las = 1, ylab = "Allelic Fraction", xlim = xlim, 
                  ...)
            }else{
              stop("Need to specify \"resultType\": HaplotypeRatio or AllelicRatio")
            }
            lines(as.numeric(c(1, dataByChr[nrow(dataByChr), Position])), rep(0.5, 2), type = "l", 
                  col = "grey", lwd = 3)
            
            if (!is.null(geneAnnot)) {
                plotGeneAnnotation(geneAnnot, i)
            }
        }
    } else {
        if (colType == "Haplotypes" || !"TITANcall" %in% colnames(dataIn)){
              colors.1 <- lohCol.hap[as.character(dataIn$phaseSet.id)]
              colors.2 <- lohCol.hap[as.character(as.numeric(!dataIn$phaseSet.id))]
            }else{
              colors.1 <- lohCol.titan[dataIn[, TITANcall]]
              colors.2 <- colors.1
        }
        
        # plot for all chromosomes specified
        dataIn <- dataIn[Chr %in% chr]
        coord <- getGenomeWidePositions(dataIn[, Chr], dataIn[, Position])
        if (resultType == "HaplotypeRatio"){
          plot(coord$posns, as.numeric(dataIn[, HaplotypeRatio.1]), 
            col = colors.1, pch = 16, 
            xaxt = "n", bty = "n", las = 1, ylab = "Haplotype Fraction", ...)
          points(coord$posns, dataIn[, HaplotypeRatio.2], col = colors.1, pch=16, ...)
        }else if (resultType == "AllelicRatio"){
          plot(coord$posns, as.numeric(dataIn[, AllelicRatio]), 
            col = colors.2, pch = 16, 
            xaxt = "n", bty = "n", las = 1, ylab = "Allelic Fraction", ...)
        }else{
              stop("Need to specify \"resultType\": HaplotypeRatio or AllelicRatio")
        }
        lines(as.numeric(c(1, coord$posns[length(coord$posns)])), 
            rep(0.5, 2), type = "l", col = "grey", 
            lwd = 3)
        plotChrLines(unique(dataIn[, Chr]), coord$chrBkpt, 
            c(-0.1, 1.1))
        
    }
}

keepChr <- function(tumour_reads, chr = c(1:22,"X","Y")){	
	tumour_reads <- tumour_reads[space(tumour_reads) %in% chr, ]
	tumour_reads <- as.data.frame(tumour_reads)
	tumour_reads$space <- droplevels(tumour_reads$space)
	tumour_reads$space <- factor(tumour_reads$space,levels=chr)
	tumour_reads <- as(tumour_reads,"RangedData")
	return(tumour_reads)
}

