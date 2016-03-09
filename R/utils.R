# author: Gavin Ha 
<<<<<<< HEAD
# 		  Dana-Farber Cancer Institute
#		  Broad Institute
# contact: <gavinha@gmail.com> or <gavinha@broadinstitute.org>
# date:	  January 20, 2015
=======
# 		Dana-Farber Cancer Institute
#		  Broad Institute
# contact: <gavinha@gmail.com> or <gavinha@broadinstitute.org>
# date:	  October 19, 2015
>>>>>>> master

loadDefaultParameters <- function(copyNumber = 5, numberClonalClusters = 1, 
    skew = 0, symmetric = TRUE, data = NULL) {
    if (copyNumber < 3 || copyNumber > 8) {
        stop("loadDefaultParameters: Fewer than 3 or more than 8 copies are 
             being specified. Please use minimum 3 or maximum 8 'copyNumber'.")
    }
    ## Data without allelic skew rn is theoretical
    ## normal reference allelic ratio initialize to
    ## theoretical values
    rn <- 0.5
    if (symmetric) {
        rt = c(rn, 1, 1, 1/2, 1, 2/3, 1, 3/4, 2/4, 
            1, 4/5, 3/5, 1, 5/6, 4/6, 3/6, 1, 6/7, 
            5/7, 4/7, 1, 7/8, 6/8, 5/8, 4/8)
        rt = rt + skew
        rt[rt > 1] <- 1
        rt[rt < (rn + skew)] <- rn + skew
        ## shift heterozygous states to account for noise 
        ##   when using symmetric = TRUE
        hetARshift <- 0.55
        if (!is.null(data)){
        	hetARshift <- median(data$ref / data$tumDepth, na.rm = TRUE)
        }
        rt[c(4, 9, 25)] <- hetARshift
        ZS = 0:24
        ct = c(0, 1, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 
            6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 8)
        highStates <- c(1,10:length(rt))
        hetState <- c(4, 9, 16, 25)
    } else {
        rt = c(rn, 1, 1e-05, 1, 1/2, 1e-05, 1, 2/3, 
            1/3, 1e-05, 1, 3/4, 2/4, 1/4, 1e-05, 1, 
            4/5, 3/5, 2/5, 1/5, 1e-05, 1, 5/6, 4/6, 
            3/6, 2/6, 1/6, 1e-05, 1, 6/7, 5/7, 4/7, 
            3/7, 2/7, 1/7, 1e-05, 1, 7/8, 6/8, 5/8, 
            4/8, 3/8, 2/8, 1/8, 1e-05)
        rt = rt + skew
        rt[rt > 1] <- 1
        rt[rt < 0] <- 1e-05
        ZS = c(0, 1, 1, 2, 3, 2, 4, 5, 5, 4, 6, 7, 
            8, 7, 6, 9, 10, 11, 11, 10, 9, 12, 13, 
            14, 15, 14, 13, 12, 16, 17, 18, 19, 19, 
            18, 17, 16, 20, 21, 22, 23, 24, 23, 22, 
            21, 20)
        ct = c(0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 
            4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 
            6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 
            8, 8, 8, 8, 8, 8, 8)
        highStates <- c(1,16:length(rt))
        hetState <- c(5, 13, 25, 41)
    }
    ZS[hetState[1]] <- -1
    rn = rn + skew
    ind <- ct <= copyNumber
    rt <- rt[ind]
    ZS <- ZS[ind]
    ct <- ct[ind]  #reassign rt and ZS based on specified copy number
    K = length(rt)
    ## VARIANCE for Gaussian to model copy number, var
    var_0 = rep(1/20, K)
    ## Dirichlet hyperparameter for initial state
    ## distribution, kappaG
    kappaGHyper = rep(1, K) + 1
    kappaGHyper[hetState[hetState %in% 1:K]] = 5
    ## Gather all genotype related parameters into a
    ## list
    genotypeParams <- vector("list", 0)
    genotypeParams$rt <- rt
    genotypeParams$rn <- rn
    genotypeParams$ZS <- ZS
    genotypeParams$ct <- ct
    genotypeParams$var_0 <- var_0
    genotypeParams$alphaKHyper <- rep(15000, K)
    varHyperHigh <- 15000
    genotypeParams$alphaKHyper[highStates] <- varHyperHigh  #AMP(11-15),HLAMP(16-21) states
    genotypeParams$betaKHyper <- rep(25, K)
    genotypeParams$kappaGHyper <- kappaGHyper
    genotypeParams$outlierVar <- 10000
    genotypeParams$symmetric <- symmetric
    ## NORMAL, n
    normalParams <- vector("list", 0)
    normalParams$n_0 <- 0.5
    normalParams$alphaNHyper <- 2
    normalParams$betaNHyper <- 2
    rm(list = c("rt", "rn", "ZS", "ct", "var_0", "kappaGHyper", 
        "skew"))
    ## PLOIDY, phi
    ploidyParams <- vector("list", 0)
    ploidyParams$phi_0 <- 2
    ploidyParams$alphaPHyper <- 20
    ploidyParams$betaPHyper <- 42
    
    ## CELLULAR PREVALENCE, s
    sParams <- setupClonalParameters(Z = numberClonalClusters)
    
    # return
    output <- vector("list", 0)
    output$genotypeParams <- genotypeParams
    output$ploidyParams <- ploidyParams
    output$normalParams <- normalParams
    output$cellPrevParams <- sParams
    return(output)
}


loadAlleleCounts <- function(inCounts, symmetric = TRUE, 
			genomeStyle = "NCBI", sep = "\t", header = TRUE) {
	if (is.character(inCounts)){
    #### LOAD INPUT READ COUNT DATA ####
    	message("titan: Loading data ", inCounts)
    	data <- read.delim(inCounts, header = header, stringsAsFactors = FALSE, 
        		sep = sep)
        if (typeof(data[,2])!="integer" || typeof(data[,4])!="integer" || 
        		typeof(data[,6])!="integer"){
        	stop("loadAlleleCounts: Input counts file format does not 
        		match required specifications.")		
        }
    }else if (is.data.frame(inCounts)){  #inCounts is a data.frame
    	data <- inCounts
    }else{
    	stop("loadAlleleCounts: Must provide a filename or data.frame 
    		to inCounts")
    }
    ## use GenomeInfoDb
    #require(GenomeInfoDb)
    # convert to desired genomeStyle and only include autosomes, sex chromosomes
    data[, 1] <- setGenomeStyle(data[, 1], genomeStyle)
   
    ## sort chromosomes
	indChr <- orderSeqlevels(as.character(data[, 1]), X.is.sexchrom = TRUE)
	data <- data[indChr, ]
	## sort positions within each chr
	for (x in unique(data[, 1])){
		ind <- which(data[, 1] == x)
		data[ind, ] <- data[ind[sort(data[ind, 2], index.return = TRUE)$ix], ]
	}
    
    refOriginal <- as.numeric(data[, 4])
    nonRef <- as.numeric(data[, 6])
    tumDepth <- refOriginal + nonRef
    if (symmetric) {
        ref <- apply(cbind(refOriginal, nonRef), 1, max, na.rm = TRUE)
    } else {
        ref <- refOriginal
    }
    
    return(list(chr = data[, 1], posn = data[, 2], ref = ref, 
        refOriginal = refOriginal, nonRef = nonRef, 
        tumDepth = tumDepth))
}

extractAlleleReadCounts <- function(bamFile, bamIndex, 
			positions, outputFilename = NULL, 
			pileupParam = PileupParam()){
	#require(Rsamtools)

## read in vcf file of het positions
	vcf <- BcfFile(positions)
	vcfPosns <- scanBcf(vcf)

## setup the positions of interest to generate the pileup for
	which <- GRanges(as.character(vcfPosns$CHROM), 
		IRanges(vcfPosns$POS, width = 1))
## setup addition BAM filters, such as excluding duplicate reads
	sbp <- ScanBamParam(flag = scanBamFlag(isDuplicate = FALSE), which = which)

## generate pileup using function (Rsamtools >= 1.17.11)
## this step can take a while
	tumbamObj <- BamFile(bamFile, index = bamIndex)
	counts <- pileup(tumbamObj, scanBamParam = sbp,  pileupParam = pileupParam)

## set of command to manipulate the "counts" data.frame output
##     by pileup() such that multiple nucleotides are in a single
##     row rather than in multiple rows.
	countsMerge <- xtabs(count ~ which_label + nucleotide, counts)
	
	label <- do.call(rbind, strsplit(rownames(countsMerge), ":"))
	posn <- do.call(rbind, strsplit(label[, 2],"-"))
	countsMerge <- cbind(data.frame(chr = label[, 1]), 
		position = posn[, 1], countsMerge[,1:7])
	
## GET REFERENCE AND NON-REF READ COUNTS
## this block of code is used to match up the reference and 
##   non-reference nucleotide when assigning read counts
##   final output data.frame is "countMat"
## setup output data.frame
	countMat <- data.frame(chr = vcfPosns$CHROM, 
			position = as.numeric(vcfPosns$POS), 
			ref = vcfPosns$REF, refCount = 0, 
			Nref = vcfPosns$ALT, NrefCount = 0, 
			stringsAsFactors = FALSE)

## match rows with vcf positions of interest
	countMat <- merge(countMat, countsMerge, by = c("chr","position"), 
		sort = FALSE, stringsAsFactors = FALSE)

## assign the flattened table of nucleotide counts to ref, Nref
## note that non-reference (Nref) allele is sum of other bases
##    that is not matching the ref.
	NT <- c("A", "T", "C", "G")
	for (n in 1:length(NT)){	
		indRef <- countMat$ref == NT[n]
		countMat[indRef, "refCount"] <- countMat[indRef, NT[n]]
		countMat[indRef, "NrefCount"] <- rowSums(countMat[indRef, NT[-n]])
	}

## remove "chr" string from chromosome
	countMat$chr <- gsub("chr","",countMat$chr)	
## only use autosomes and sex chrs
	countMat <- countMat[countMat$chr %in% c(as.character(1:22),"X","Y"),]
## only use first 6 columns for TitanCNA
	countMat <- countMat[,1:6]

## OUTPUT TO TEXT FILE 
	if (!is.null(outputFilename)){
		## output text file will have the same format required by TitanCNA
		message("extractAlleleReadCounts: writing to ", outputFilename)
		write.table(countMat, file = outputFilename, row.names = FALSE, 
			col.names = TRUE, quote = FALSE, sep = "\t")
	}
	return(countMat)
	#return(loadAlleleCounts(countMat))
}

## filter data by depth and mappability.  data is a
## logical vector data is a list containing all our
## data: ref, depth, logR, etc.  
filterData <- function(data, chrs = NULL, minDepth = 10, 
    maxDepth = 200, positionList = NULL, map = NULL, 
    mapThres = 0.9, centromeres = NULL, centromere.flankLength = 0) {
    if (!is.null(map)) {
        keepMap <- map >= mapThres
    } else {
        keepMap <- !logical(length = length(data$refOriginal))
    }
    if (!is.null(positionList)) {
        chrPosnList <- paste(positionList[, 1], positionList[, 
            2], sep = ":")  #chr:posn
        chrPosnData <- paste(data$chr, data$posn, sep = ":")
        keepPosn <- is.element(chrPosnData, chrPosnList)
    } else {
        keepPosn <- !logical(length = length(data$chr))
    }
    keepTumDepth <- data$tumDepth <= maxDepth & data$tumDepth >= 
        minDepth
    if (is.null(chrs)){
   		keepChrs <- logical(length = length(data$chr))
   	}else{
   		keepChrs <- is.element(data$chr, chrs)
   	}
    cI <- keepChrs & keepTumDepth & !is.na(data$logR) & 
        keepMap & keepPosn
    for (i in 1:length(data)) {
        if (!is.null(data[[i]])) {
            data[[i]] <- data[[i]][cI]
        }
        
    }
    ## remove centromere SNPs ##
    if (!is.null(centromeres)){
    	colnames(centromeres)[1:3] <- c("space", "start", "ends") 
    	data <- removeCentromere(data, centromeres, flankLength = centromere.flankLength)
    }
    return(data)
}

## input: 
# 1) data object output by loadAlleleCounts(); 6 element list: chr, posn, ref, refOriginal, nonRef, tumDepth)
# 2) data.frame containing coordinates of centromeres; 4 columns: Chr, Start, End, arbitrary
##
removeCentromere <- function(data, centromere, flankLength = 0){
	keepInd <- !logical(length = length(data$chr))
	for (c in 1:nrow(centromere)){
		ind <- which((data$chr == centromere[c, "Chr"]) &
				(data$posn >= (centromere[c, "Start"] - flankLength)) &
				(data$posn <= (centromere[c, "End"] + flankLength)))
		keepInd[ind] <- FALSE			
	}	
	message("Removed ", sum(!keepInd), " centromeric positions")
		## remove positions in all elements of list
	for (i in 1:length(data)) {
        if (!is.null(data[[i]])) {
            data[[i]] <- data[[i]][keepInd]
        }        
    }
    return(data)
}


excludeGarbageState <- function(params, K) {
    newParams <- params
    for (i in 1:length(newParams)) {
        if (length(newParams[[i]]) == K) {
            newParams[[i]] <- newParams[[i]][2:K]
        }
    }
    return(newParams)
}

getPositionOverlap <- function(chr, posn, dataVal) {
# use RangedData to perform overlap
    dataIR <- RangedData(space = dataVal[, 1], 
    				IRanges(start = dataVal[, 2], end = dataVal[, 3]),
    				val = as.numeric(dataVal[, 4]))
    				
    ## load chr/posn as data.frame first to use proper chr ordering by factors/levels
    chrDF <- data.frame(space=chr,start=posn,end=posn)
    chrDF$space <- factor(chrDF$space, levels = unique(chr))    
    chrIR <- as(chrDF, "RangedData")
    
    hits <- findOverlaps(query = chrIR, subject = dataIR)
    
    ## create full dataval list ##
    hitVal <- rep(NA, length = length(chr))
    hitVal[queryHits(hits)] <- dataIR$val[subjectHits(hits)]
    #chrIR$hitVal <- hitVal
    ## reorder to match input chr and posn arguments
    #chrDF <- as.data.frame(chrIR)
    #indReorder <- order(match(chrDF[, 1], chr))
    #return(hitVal[indReorder])
	return(hitVal) 
 
    #cnChr <- cnData[, 1]
    #cnStart <- as.numeric(cnData[, 2])
    #cnStop <- as.numeric(cnData[, 3])
    #cnVal <- as.numeric(cnData[, 4])
    
    #N = length(posn)
    #valByPosn = rep(NA, N)
    
    #for (c in unique(chr)) {
    #    indData <- chr == c
    #    indCN <- cnChr == c
    #    cnStartC <- cnStart[indCN]
    #    cnStopC <- cnStop[indCN]
    #    cnValC <- cnVal[indCN]
    #    posnC <- posn[indData]
    #    cnInd <- .Call("getPositionOverlapC", posnC, 
    #        cnStartC, cnStopC)
    #    if (sum(cnInd > 0) > 0) {
    #        cnValToUse <- rep(NA, length(cnInd))
    #        cnValToUse[which(cnInd > 0)] <- cnValC[cnInd]
    #        valByPosn[indData] <- cnValToUse
    #    }
    #}
    
    #return(as.numeric(valByPosn))
}

setGenomeStyle <- function(x, genomeStyle = "NCBI", species = "Homo_sapiens"){
	#chrs <- genomeStyles(species)[c("NCBI","UCSC")]
	if (!genomeStyle %in% seqlevelsStyle(as.character(x))){
    	x <- suppressWarnings(mapSeqlevels(as.character(x), 
    					genomeStyle, drop = FALSE)[1,])
    }
    
    autoSexMChr <- extractSeqlevelsByGroup(species = species, 
    				style = genomeStyle, group = "all")
    x <- x[x %in% autoSexMChr]
    return(x)
}

correctReadDepth <- function(tumWig, normWig, gcWig, mapWig, 
	genomeStyle = "NCBI", targetedSequence = NULL) {
    #require(HMMcopy)
    
    message("Reading GC and mappability files")
    gc <- wigToRangedData(gcWig)
    map <- wigToRangedData(mapWig)
    
    ### LOAD TUMOUR AND NORMAL FILES ###
    message("Loading tumour file:", tumWig)
    tumour_reads <- wigToRangedData(tumWig)
    message("Loading normal file:", normWig)
    normal_reads <- wigToRangedData(normWig)
    
    ### set the genomeStyle: NCBI or UCSC
    #require(GenomeInfoDb)
	names(gc) <- setGenomeStyle(names(gc), genomeStyle)
	names(map) <- setGenomeStyle(names(map), genomeStyle)
	names(tumour_reads) <- setGenomeStyle(names(tumour_reads), genomeStyle)
	names(normal_reads) <- setGenomeStyle(names(normal_reads), genomeStyle)
    
    ### make sure tumour wig and gc/map wigs have same
    ### chromosomes
    gc <- gc[gc$space %in% tumour_reads$space, ]
    map <- map[map$space %in% tumour_reads$space, ]
    samplesize <- 50000
    
    ### for targeted sequencing (e.g.  exome capture),
    ### ignore bins with 0 for both tumour and normal
    ### targetedSequence = RangedData (IRanges object)
    ### containing list of targeted regions to consider;
    ### 3 columns: chr, start, end
    if (!is.null(targetedSequence)) {
        message("Analyzing targeted regions...")
        targetIR <- RangedData(ranges = IRanges(start = targetedSequence[, 2], 
                    end = targetedSequence[, 3]), space = targetedSequence[, 1])
                    
        #keepInd <- unlist(as.list(findOverlaps(tumour_reads, targetIR, select = "first")))
        #keepInd <- !is.na(keepInd)
        hits <- findOverlaps(query = tumour_reads, subject = targetIR)
        keepInd <- unique(queryHits(hits))    
        
        # ind <- tumour_reads$value>10 &
        # normal_reads$value>10 tumThres <-
        # quantile(tumour_reads$value[ind],1/4) normThres
        # <- quantile(normal_reads$value[ind],1/4) keepInd
        # <- which(ind & !is.na(tumour_reads$value) &
        # !is.na(normal_reads$value) &
        # tumour_reads$value>=tumThres &
        # normal_reads$value>=normThres)
        tumour_reads <- tumour_reads[keepInd, ]
        normal_reads <- normal_reads[keepInd, ]
        gc <- gc[keepInd, ]
        map <- map[keepInd, ]
        samplesize <- min(ceiling(nrow(tumour_reads) * 
            0.1), samplesize)
    }
    
    ### add GC and Map data to IRanges objects ###
    tumour_reads$gc <- gc$value
    tumour_reads$map <- map$value
    colnames(tumour_reads) <- c("reads", "gc", "map")
    normal_reads$gc <- gc$value
    normal_reads$map <- map$value
    colnames(normal_reads) <- c("reads", "gc", "map")
    
    ### CORRECT TUMOUR DATA FOR GC CONTENT AND
    ### MAPPABILITY BIASES ###
    message("Correcting Tumour")
    tumour_copy <- correctReadcount(tumour_reads, samplesize = samplesize)
    
    ### CORRECT NORMAL DATA FOR GC CONTENT AND
    ### MAPPABILITY BIASES ###
    message("Correcting Normal")
    normal_copy <- correctReadcount(normal_reads, samplesize = samplesize)
    
    ### COMPUTE LOG RATIO ###
    message("Normalizing Tumour by Normal")
    tumour_copy$copy <- tumour_copy$copy - normal_copy$copy
    rm(normal_copy)
    
    ### PUTTING TOGETHER THE COLUMNS IN THE OUTPUT ###
    temp <- cbind(chr = as.character(space(tumour_copy)), 
        start = start(tumour_copy), end = end(tumour_copy), 
        logR = tumour_copy$copy)
    temp <- as.data.frame(temp, stringsAsFactors = FALSE)
    mode(temp$start) <- "numeric"
    mode(temp$end) <- "numeric"
    mode(temp$logR) <- "numeric"
    return(temp)
}


setupClonalParameters <- function(Z, sPriorStrength = 2) {
    alphaSHyper = rep(sPriorStrength, Z)
    betaSHyper = rep(sPriorStrength, Z)
    kappaZHyper = rep(1, Z) + 1
    
    # use naive initialization
    s_0 <- (1:Z)/10
    
    ## first cluster should be clonally dominant (use
    ## 0.001) ##
    s_0[1] <- 0.001
    
    output <- vector("list", 0)
    output$s_0 <- s_0
    output$alphaSHyper <- alphaSHyper
    output$betaSHyper <- betaSHyper
    output$kappaZHyper <- kappaZHyper
    return(output)
}

computeBIC <- function(maxLoglik, M, N) {
    bic <- -2 * maxLoglik + (M * log(N))
    return(bic)
}

computeSDbwIndex <- function(x, centroid.method = "median", data.type = "LogRatio", 
						S_Dbw.method = "Halkidi", symmetric = TRUE) {
    ## input x: Titan results dataframe from
    ## 'outputTitanResults()' S_Dbw Validity Index
    ## Halkidi and Vazirgiannis (2001). Clustering
    ## Validity Assessment: Finding the Optimal
    ## Partition of a Data Set
    ## AND
    ## Tong and Tan (2009) Cluster validity based on the 
    ## improved S_Dbw index
    
    if (!data.type %in% c("LogRatio", "AllelicRatio")){
    	stop("computeSDbwIndex: data.type must be either 'LogRatio' or 'AllelicRatio'")
    }
    
      if (!S_Dbw.method %in% c("Halkidi", "Tong")){
    	stop("computeSDbwIndex: S_Dbw.method must be either 'Halkidi' or 'Tong'")
    }
    
    ## flatten copynumber-clonalclusters to single vector
    if (data.type=="LogRatio"){
    	cn <- as.numeric(x[, "CopyNumber"]) + 1
		cn[cn == 3] <- NA  ## remove all CN=2 positions
    	flatState <- (as.numeric(x[, "ClonalCluster"]) - 1) * (max(cn, na.rm = TRUE)) + cn
    	flatState[is.na(flatState)] <- 3 ### assign all the CN=2 positions to cluster 3
    	CNdata <- scale(as.numeric(x[, data.type]))
    	x <- as.matrix(cbind(as.numeric(flatState), CNdata))
    }else if (data.type=="AllelicRatio"){
    	st <- as.numeric(x[, "TITANstate"]) + 1
    	st[x[, "TITANcall"] == "HET"] <- NA
    	flatState <- (as.numeric(x[, "ClonalCluster"]) - 1) * (max(st, na.rm = TRUE)) + st
    	if (symmetric){
    		flatState[is.na(flatState)] <- 4
    	}else{
    		flatState[is.na(flatState)] <- 5
    	}
    	## for allelic ratios, compute the symmetric allelic ratio
    	ARdata <- as.numeric(x[, data.type])
    	ARdata <- apply(cbind(ARdata, 1 - ARdata), 1, max, na.rm = TRUE)
    	ARdata <- scale(ARdata)
    	x <- as.matrix(cbind(as.numeric(flatState), ARdata))
    	rm(ARdata)
    }
    
    clust <- sort(unique(x[, 1]))
    K <- length(clust)
    N <- nrow(x)
    
    ## find average standard deviation and scatter (compactness)
    stdev <- rep(NA, K)
    scat.Ci <- rep(NA, K)
    for (i in 1:K) {
        ind.i <- x[, 1] == clust[i]
        ni <- sum(ind.i)
        stdev[i] <- var(x[ind.i, 2], na.rm = TRUE)
        
        ## compute scatter based on variances within objects
        ## of a cluster (compactness)
        var.Ci <- var(x[ind.i, 2], na.rm = TRUE)
        var.D <- var(x[, 2], na.rm = TRUE)
        if (S_Dbw.method == "Halkidi"){
        	scat.Ci[i] <- var.Ci/var.D
        }else if (S_Dbw.method == "Tong"){
        	###### NOTE #######
        	## The authors originally used ((N-ni)/N), but this is incorrect ##
        	## We want a weighted-sum ##
        	scat.Ci[i] <- ((ni) / N) * (var.Ci/var.D)
        }
    }
    avgStdev <- sqrt(sum(stdev, na.rm = TRUE))/K
    
    ## compute density between clusters (separation)
    sumDensityDiff <- matrix(NA, nrow = K, ncol = K)
    for (i in 1:K) {
        # cat('Calculating S_Dbw for cluster # ',clust[i],'\n')
        ind.i <- x[, 1] == clust[i]
        ni <- sum(ind.i)
        xci <- x[ind.i, 2]

        #density of Ci
        sumDiff.xci <- sdbw.density(xci, avgStdev, method = S_Dbw.method, 
        					centroid.method = centroid.method)
        
        for (j in 1:K) {
            if (i == j) {
                next
            }
            ind.j <- x[, 1] == clust[j]
            nj <- sum(ind.j)
            xcj <- x[ind.j, 2]

            #density of Cj
            sumDiff.xcj <- sdbw.density(xcj, avgStdev, method = S_Dbw.method, 
        					centroid.method = centroid.method)
            
            ## union and midpoint of both clusters
            x.ci.cj <- union(xci, xcj)
            ci <- median(xci, na.rm = TRUE)  #centroid of cluster Ci
            cj <- median(xcj, na.rm = TRUE)  #centroid of cluster Cj
            nij <- ni + nj
            stdCiCj <- (sd(xci) + sd(xcj)) / 2
            if (S_Dbw.method == "Halkidi"){
            	cij <- (ci + cj)/2
            	#cij <- median(x.ci.cj, na.rm = TRUE)
            }else if (S_Dbw.method == "Tong"){
            	lambda <- 0.7
            	cij <- lambda * ((nj * ci + ni * cj) / nij) + (1 - lambda) * 
            		(sumDiff.xci * ci + sumDiff.xcj * cj) /
            		(sumDiff.xci + sumDiff.xcj) 
            }
            
            
            #density of union of both clusters using special centroid
            sumDiff.xci.xcj <- sdbw.density(x.ci.cj, avgStdev, stDev = stdCiCj, 
            				method = S_Dbw.method, centroid = cij, 
            				centroid.method = centroid.method)            
            
            maxDiff <- max(sumDiff.xci, sumDiff.xcj)
            if (maxDiff == 0) {
                maxDiff <- 0.1
            }
            sumDensityDiff[i, j] <- sumDiff.xci.xcj/maxDiff
        }
    }
    
    if (S_Dbw.method == "Halkidi"){
        scat <- sum(scat.Ci, na.rm = TRUE)/(K)
    }else if (S_Dbw.method == "Tong"){
   		scat <- sum(scat.Ci, na.rm = TRUE)/(K - 1)
    }    
    dens.bw <- sum(sumDensityDiff, na.rm = TRUE)/(K * (K - 1))
    
    S_DbwIndex <- scat + dens.bw
    # return(S_DbwIndex)
    return(list(S_DbwIndex = S_DbwIndex, dens.bw = dens.bw, scat = scat))
}


sdbw.density <- function(x, avgStdev, stDev = NULL, method = "Halkidi", 
					centroid = NULL, centroid.method = "median"){
	if (is.null(centroid)){
		if (centroid.method == "median") {
        	centroid <- median(x, na.rm = TRUE)  #centroid of cluster Cj
    	} else if (centroid.method == "mean") {
        	centroid <- mean(x, na.rm = TRUE)  #centroid of cluster Cj
    	}
    }
	#density of Ci
    if (method == "Halkidi"){
        sumDiff <- sum(abs(x - centroid) <= avgStdev, na.rm = TRUE) 
    }else if (method == "Tong"){
    	if (is.null(stDev)){
    		stDev <- sd(x, na.rm = TRUE)
    	}
        conf.int <- 1.96 * (stDev / sqrt(length(x)))
        sumDiff <- sum(abs(x - centroid) <= conf.int, na.rm = TRUE)  
    }
    return(sumDiff)
}


# G = sequence of states for mega-variable K =
# number of unit states per cluster excluding
# outlier state precondition: If outlierState is
# included, it must be at G=1, else HOMD is G=1
decoupleMegaVar <- function(G, K, useOutlierState = FALSE) {
    if (useOutlierState) {
        G <- G - 1  #do this to make OUT=0 and HOMD=1
        G[G == 0] <- NA  #assign NA to OUT states
    }
    newG <- G%%K
    newG[newG == 0] <- K
    newG[is.na(newG)] <- 0
    newZ <- ceiling(G/K)
    output <- vector("list", 0)
    output$G <- newG
    output$Z <- newZ
    return(output)
}

# pre-condition: outlier state is -1 if included
decodeLOH <- function(G, symmetric = TRUE) {
    T <- length(G)
    Z <- rep("NA", T)
    CN <- rep(NA, T)
    
    if (symmetric) {
        DLOH <- G == 1
        NLOH <- G == 2
        ALOH <- G == 4 | G == 6 | G == 9 | G == 12 | 
            G == 16 | G == 20
        HET <- G == 3
        GAIN <- G == 5
        ASCNA <- G == 7 | G == 10 | G == 13 | G == 
            17 | G == 21
        BCNA <- G == 8 | G == 15 | G == 24
        UBCNA <- G == 11 | G == 14 | G == 18 | G == 
            19 | G == 22 | G == 23
    } else {
        DLOH <- G == 1 | G == 2
        NLOH <- G == 3 | G == 5
        ALOH <- G == 6 | G == 9 | G == 10 | G == 14 | 
            G == 15 | G == 20 | G == 22 | G == 28 | 
            G == 29 | G == 36 | G == 37 | G == 45
        HET <- G == 4
        GAIN <- G == 7 | G == 8
        ASCNA <- G == 11 | G == 13 | G == 16 | G == 
            19 | G == 23 | G == 27 | G == 30 | G == 
            35 | G == 38 | G == 44
        BCNA <- G == 12 | G == 25 | G == 41
        UBCNA <- G == 17 | G == 18 | G == 24 | G == 
            26 | G == 31 | G == 32 | G == 33 | G == 
            34 | G == 39 | G == 40 | G == 42 | G == 
            43
    }
    HOMD <- G == 0
    OUT <- G == -1
    
    Z[HOMD] <- "HOMD"
    Z[DLOH] <- "DLOH"
    Z[NLOH] <- "NLOH"
    Z[ALOH] <- "ALOH"
    Z[HET] <- "HET"
    Z[GAIN] <- "GAIN"
    Z[ASCNA] <- "ASCNA"
    Z[BCNA] <- "BCNA"
    Z[UBCNA] <- "UBCNA"
    Z[OUT] <- "OUT"
    
    if (symmetric) {
        CN[HOMD] <- 0
        CN[DLOH] <- 1
        CN[G >= 2 & G <= 3] <- 2
        CN[G >= 4 & G <= 5] <- 3
        CN[G >= 6 & G <= 8] <- 4
        CN[G >= 9 & G <= 11] <- 5
        CN[G >= 12 & G <= 15] <- 6
        CN[G >= 16 & G <= 19] <- 7
        CN[G >= 20 & G <= 24] <- 8
    } else {
        CN[HOMD] <- 0
        CN[DLOH] <- 1
        CN[G >= 3 & G <= 5] <- 2
        CN[G >= 6 & G <= 9] <- 3
        CN[G >= 10 & G <= 14] <- 4
        CN[G >= 15 & G <= 20] <- 5
        CN[G >= 21 & G <= 28] <- 6
        CN[G >= 29 & G <= 36] <- 7
        CN[G >= 37 & G <= 45] <- 8
    }
    
    output <- vector("list", 0)
    output$G <- Z
    output$CN <- CN
    return(output)
}


outputTitanResults <- function(data, convergeParams, 
    optimalPath, filename = NULL, posteriorProbs = FALSE, 
    subcloneProfiles = TRUE) {
    
    # check if useOutlierState is in convergeParams
    if (length(convergeParams$useOutlierState) == 0) {
        stop("convergeParams does not contain element: useOutlierState.")
    }
     useOutlierState <- convergeParams$useOutlierState
     # check if symmetric is in convergeParams
    if (length(convergeParams$symmetric) == 0) {
        stop("convergeParams does not contain element: symmetric.")
    }
   
    
    #### PROCESS HMM RESULTS ####
    numClust <- dim(convergeParams$s)[1]
    K <- dim(convergeParams$var)[1]
    if (useOutlierState) {
        K <- K - 1
    }
    Z <- dim(convergeParams$s)[1]
    i <- dim(convergeParams$s)[2]  #iteration of training to use (last iteration)
    partGZ <- decoupleMegaVar(optimalPath, K, useOutlierState)
    G <- partGZ$G - 1  #assign analyzed points, minus 2 so HOMD=0
    sortS <- sort(convergeParams$s[, i], decreasing = FALSE, 
        index.return = TRUE)
    s <- sortS$x
    Zclust <- partGZ$Z  #assign analyzed points
    Zclust <- sortS$ix[Zclust]  #reassign sorted cluster membership
    
    ### OUTPUT RESULTS #### Output Z ##
    Gdecode <- decodeLOH(G, symmetric = convergeParams$symmetric)
    Gcalls <- Gdecode$G
    CN <- Gdecode$CN
    Zclust[Gcalls == "HET" & CN == 2] <- NA  #diploid HET positions do not have clusters
    Sout <- rep(NA, length(Zclust))  #output cluster frequencies
    Sout[Zclust > 0 & !is.na(Zclust)] <- s[Zclust[Zclust > 
        0 & !is.na(Zclust)]]
    Sout[Zclust == 0] <- 0
    clonalHeaderStr <- rep(NA, Z)
    for (j in 1:Z) {
        clonalHeaderStr[j] <- sprintf("pClust%d", j)
    }
    
    outmat <- cbind(Chr = data$chr, Position = data$posn, 
        RefCount = data$refOriginal, NRefCount = data$tumDepth - 
            data$refOriginal, Depth = data$tumDepth, 
        AllelicRatio = sprintf("%0.2f", data$refOriginal/data$tumDepth), 
        LogRatio = sprintf("%0.2f", log2(exp(data$logR))), 
        CopyNumber = CN, TITANstate = G, TITANcall = Gcalls, 
        ClonalCluster = Zclust, CellularPrevalence = sprintf("%0.2f", 
            1 - Sout))
   	
   	## INCLUDE SUBCLONE PROFILES 
   	if (subcloneProfiles & numClust <= 2){
   		outmat <- as.data.frame(outmat, stringsAsFactors = FALSE)
    	outmat <- cbind(outmat, getSubcloneProfiles(outmat))
    }else{
    	message("outputTitanResults: More than 2 clusters or subclone profiles not requested.")
    }
    if (posteriorProbs) {
    	rhoG <- t(convergeParams$rhoG)
    	rhoZ <- t(convergeParams$rhoZ)
    	rhoZ <- rhoZ[, sortS$ix, drop = FALSE]
   		colnames(rhoZ) <- clonalHeaderStr
        outmat <- cbind(outmat, format(round(rhoZ, 
            4), nsmall = 4, scientific = FALSE), format(round(rhoG, 
            4), nsmall = 4, scientific = FALSE))
    }
    if (!is.null(filename)) {
        message("titan: Writing results to ", filename)
        write.table(outmat, file = filename, col.names = TRUE, 
            row.names = FALSE, quote = FALSE, sep = "\t")
    }
    return(as.data.frame(outmat, stringsAsFactors = FALSE))
}

outputModelParameters <- function(convergeParams, results, filename, 
		S_Dbw.scale = 1, S_Dbw.method = "Tong") {
    message("titan: Saving parameters to ", filename)
    Z <- dim(convergeParams$s)[1]
    i <- dim(convergeParams$s)[2]  #iteration of training to use (last iteration)
    sortS <- sort(convergeParams$s[, i], decreasing = FALSE, 
        index.return = TRUE)
    s <- sortS$x
    fc <- file(filename, "w+")
    norm_str <- sprintf("Normal contamination estimate:\t%0.2f", 
        convergeParams$n[i])
    write.table(norm_str, file = fc, col.names = FALSE, 
        row.names = FALSE, quote = FALSE, sep = "", append = TRUE)
    ploid_str <- sprintf("Average tumour ploidy estimate:\t%0.2f", 
        convergeParams$phi[i])
    write.table(ploid_str, file = fc, col.names = FALSE, 
        row.names = FALSE, quote = FALSE, sep = "", append = TRUE)
    s_str <- sprintf("%0.4f ", 1 - s)
    s_str <- gsub(" ", "", s_str)
    outStr <- sprintf("Clonal cluster cellular prevalence Z=%d:\t%s", 
        Z, paste(s_str, collapse = " "))
    write.table(outStr, file = fc, col.names = FALSE, 
        row.names = FALSE, quote = FALSE, sep = "", append = TRUE)
    for (j in 1:Z) {
        musR_str <- sprintf("%0.2f ", convergeParams$muR[, j, i])
        musR_str <- gsub(" ", "", musR_str)
        outStr <- sprintf("Genotype binomial means for clonal cluster Z=%d:\t%s", j, paste(musR_str, collapse = " "))
        write.table(outStr, file = fc, col.names = FALSE, 
            row.names = FALSE, quote = FALSE, sep = "", 
            append = TRUE)
        musC_str <- sprintf("%0.2f ", log2(exp(convergeParams$muC[, j, i])))
        musC_str <- gsub(" ", "", musC_str)
        outStr <- sprintf("Genotype Gaussian means for clonal cluster Z=%d:\t%s", j, paste(musC_str, collapse = " "))
        write.table(outStr, file = fc, col.names = FALSE, 
            row.names = FALSE, quote = FALSE, sep = "", 
            append = TRUE)
    }
    var_str <- sprintf("%0.4f ", convergeParams$var[, i])
    var_str <- gsub(" ", "", var_str)
    outStr <- sprintf("Genotype Gaussian variance:\t%s", 
        paste(var_str, collapse = " "))
    write.table(outStr, file = fc, col.names = FALSE, 
        row.names = FALSE, quote = FALSE, sep = "", 
        append = TRUE)
    iter_str <- sprintf("Number of iterations:\t%d", 
        length(convergeParams$phi))
    write.table(iter_str, file = fc, col.names = FALSE, 
        row.names = FALSE, quote = FALSE, sep = "", 
        append = TRUE)
    loglik_str <- sprintf("%0.4f ", convergeParams$loglik[i])
    loglik_str <- gsub(" ", "", loglik_str)
    outStr <- sprintf("Log likelihood:\t%s", paste(loglik_str, 
        collapse = " "))
    write.table(outStr, file = fc, col.names = FALSE, 
        row.names = FALSE, quote = FALSE, sep = "", 
        append = TRUE)
    
    # compute SDbw_index
    sdbw.LR <- computeSDbwIndex(results, centroid.method = "median", 
    					data.type = "LogRatio", 
    					S_Dbw.method = S_Dbw.method,
    					symmetric = convergeParams$symmetric)
    sdbw.AR <- computeSDbwIndex(results, centroid.method = "median", 
    					data.type = "AllelicRatio", 
    					S_Dbw.method = S_Dbw.method,
    					symmetric = convergeParams$symmetric)
    ## element-wise addition -> returns list
    ## add the values for allelicRatio and logRatio
    sdbw <- mapply('+', sdbw.LR, sdbw.AR, SIMPLIFY = FALSE)        
    
    ## print out combined S_Dbw ##
	printSDbw(sdbw.LR, fc, S_Dbw.scale, "LogRatio")
    printSDbw(sdbw.AR, fc, S_Dbw.scale, "AllelicRatio")
    printSDbw(sdbw, fc, S_Dbw.scale, "Both")
    close(fc)
    
    return(list(dens.bw = sdbw$dens.bw, scat = sdbw$scat, 
    			S_Dbw = S_Dbw.scale * sdbw$dens.bw + sdbw$scat))
    
} 


printSDbw <- function(sdbw, fc, scale, data.type = ""){
	sdbw_str <- sprintf("S_Dbw dens.bw (%s):\t%0.4f ", data.type, sdbw$dens.bw)
    write.table(sdbw_str, file = fc, col.names = FALSE, 
        row.names = FALSE, quote = FALSE, sep = "", 
        append = TRUE)
    sdbw_str <- sprintf("S_Dbw scat (%s):\t%0.4f ", data.type, sdbw$scat)
    write.table(sdbw_str, file = fc, col.names = FALSE, 
        row.names = FALSE, quote = FALSE, sep = "", 
        append = TRUE)
    sdbw_str <- sprintf("S_Dbw validity index (%s):\t%0.4f ", 
        data.type, scale * sdbw$dens.bw + sdbw$scat)
    write.table(sdbw_str, file = fc, col.names = FALSE, 
        row.names = FALSE, quote = FALSE, sep = "", 
        append = TRUE)
}

<<<<<<< HEAD
=======
## TODO: Add documentation
removeEmptyClusters <- function(convergeParams, results, proportionThreshold = 0.001, 
	proportionThresholdClonal = 0.3){
	clust <- 1:nrow(convergeParams$s)
	names(clust) <- clust
	#newClust <- clust #original clusters
	for (cl in clust){
		ind <- which(results$ClonalCluster == cl)
		if (length(ind) / nrow(results) < proportionThreshold || 
				(length(ind) / nrow(results) < proportionThresholdClonal && cl == 1)){
			#newClust <- newClust[-which(names(newClust) == cl)]
			clust[cl] <- NA #assign cluster without sufficient data with NA
		}
	}
	k <- ncol(convergeParams$s)
	# sort the cellular prevalence since they are sorted in "results"
	convergeParams$s <- convergeParams$s[order(convergeParams$s[, k], decreasing = FALSE), , drop = FALSE]
	# if there is at least 1 cluster with sufficient data
	if (length(which(clust > 0)) > 0){
		#set new normal estimate as cluster 1
		#purity <- as.numeric(unique(results$CellularPrevalence[which(results$ClonalCluster == which(!is.na(clust))[1])]))
		purity <- (1 - convergeParams$s[which(!is.na(clust))[1], k]) * (1 -  convergeParams$n[k])
		convergeParams$n[k] <- 1 - purity
		#set new cellular prevalence using new clusters and renormalize to new cluster 1		
		convergeParams$s <- convergeParams$s[which(!is.na(clust)), , drop = FALSE]
		convergeParams$s[, k] <- 1 - (1 - convergeParams$s[, k]) / (1 - convergeParams$s[1, k])
		#names(newClust) <- 1:length(newClust)
		clust[!is.na(clust)] <- 1:sum(!is.na(clust))

		#set new cellular prevalence and clonal cluster in results file	
		for (cl in 1:length(clust)){
			# assign data in removed cluster cl to next non-NA cluster
			if (is.na(clust[cl])){
				# assign to the right (larger cluster number)
				if (length(which(!is.na(clust) & names(clust) > cl)) > 0){
					results[which(results$ClonalCluster == names(clust)[cl]), "CellularPrevalence"] <- 1 - convergeParams$s[clust[which(!is.na(clust) & names(clust) > cl)[1]], k]
					results[which(results$ClonalCluster == names(clust)[cl]), "ClonalCluster"] <- clust[which(!is.na(clust) & names(clust) > cl)][1]
				# assign to the left (smaller cluster number)
				}else if (length(which(!is.na(clust) & names(clust) < cl)) > 0){
					results[which(results$ClonalCluster == names(clust)[cl]), "CellularPrevalence"] <- 1 - convergeParams$s[clust[tail(which(!is.na(clust) & names(clust) < cl), 1)], k]
					results[which(results$ClonalCluster == names(clust)[cl]), "ClonalCluster"] <- clust[tail(which(!is.na(clust) & names(clust) < cl), 1)]
				}
			}else{ # update cluster and cellPrev info for kept clusters
				results[which(results$ClonalCluster == names(clust)[cl]), "CellularPrevalence"] <- 1 - convergeParams$s[clust[cl], k]
				results[which(results$ClonalCluster == names(clust)[cl]), "ClonalCluster"] <- clust[cl]		
			}
		}
	}else{ # no clusters with sufficient data
		
		# set params to only cluster with data or to default cluster01 if no cluster with data
		clustData <- which.max(table(results$ClonalCluster))
		if (length(clustData) >= 0){
			clustData <- 1
		}
		#set normal contamination to 100%
		convergeParams$n[k] <- 1 - convergeParams$s[clustData, k]
		convergeParams$s <- convergeParams$s[clustData, , drop = FALSE]
		convergeParams$s[, k] <- 0.0
		
			# set all clusters to 1 and all cellular prevalence to 1.0; leave HET as NA
		results[which(results$TITANcall != "HET"), "CellularPrevalence"] <- 1
		results[which(results$TITANcall != "HET"), "ClonalCluster"] <- 1
	}	
		
	return(list(convergeParams = convergeParams, results = results))
}
>>>>>>> master

getSubcloneProfiles <- function(titanResults){
	if (is.character(titanResults)){
		titanResults <- read.delim(titanResults, header = TRUE, 
				stringsAsFactors = FALSE, sep = "\t")
	}else if (!is.data.frame(titanResults)){
		stop("getSubcloneProfiles: titanResults is not character or
				data.frame.")
	}
	
	numClones <- as.numeric(max(titanResults$ClonalCluster,
			na.rm = TRUE))
	if (is.na(numClones)){ numClones <- 0 }
	cellPrev <- unique(cbind(Cluster = titanResults$ClonalCluster, 
			Prevalence = titanResults$CellularPrevalence))
	
	if (numClones == 0){
		subc1 <- as.data.frame(cbind(CopyNumber = as.numeric(titanResults$CopyNumber), 
				TITANcall = titanResults$TITANcall,
				Prevalence = "NA"), stringsAsFactors = FALSE)
	}
	if (numClones == 1){
		subc1Prev <- cellPrev[which(cellPrev[, "Cluster"] == "1"), "Prevalence"]
		subc1 <- as.data.frame(cbind(CopyNumber = as.numeric(titanResults$CopyNumber), 
				TITANcall = titanResults$TITANcall,
				Prevalence = as.numeric(subc1Prev)), stringsAsFactors = FALSE)
	}
	if (numClones == 2){
		subc2Prev <- as.numeric(cellPrev[which(cellPrev[, "Cluster"] == "2"),
				"Prevalence"])
		subc1Prev <- as.numeric(cellPrev[which(cellPrev[, "Cluster"] == "1"),
				"Prevalence"])
		subc1Prev <- subc1Prev - subc2Prev
		subc2 <- as.data.frame(cbind(CopyNumber = as.numeric(titanResults$CopyNumber), 
			TITANcall = titanResults$TITANcall,
			Prevalence = as.numeric(subc2Prev)), stringsAsFactors = FALSE)
		mode(subc2[, 1]) <- "numeric"; mode(subc2[, 3]) <- "numeric"
		subc1 <- subc2
		ind <- which(titanResults$ClonalCluster == 2)
		subc1[ind, c("CopyNumber", "TITANcall")] <- t(matrix(cbind(2, "HET"), 
				ncol = length(ind), nrow = 2))
		subc1[, "Prevalence"] <- subc1Prev
	}
	
	## Add subclone 1, 2 and 3 if they are defined
	outMat <- cbind(Subclone1 = subc1)
	if (exists("subc2")){
		outMat <- cbind(outMat, Subclone2 = subc2)
	}
	#if (exists("subc3")){
	#	outMat <- cbind(outMat, Subclone3 = subc3, stringsAsFactors = FALSE)
	#}
	return(as.data.frame(outMat, stringsAsFactors = FALSE))
}
