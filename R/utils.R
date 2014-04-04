# file : utils.R author: Gavin Ha <gha@bccrc.ca>
# Dept of Molecular Oncolgy British Columbia Cancer
# Agency University of British Columbia date :
# March 2, 2014

loadDefaultParameters <- function(copyNumber = 5, numberClonalClusters = 1, 
    skew = 0, symmetric = TRUE) {
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
        ZS = 0:24
        ct = c(0, 1, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 
            6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 8)
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
    }
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
    kappaGHyper[5] = 5
    ## Gather all genotype related parameters into a
    ## list
    genotypeParams <- vector("list", 0)
    genotypeParams$rt <- rt
    genotypeParams$rn <- rn
    genotypeParams$ZS <- ZS
    genotypeParams$ct <- ct
    genotypeParams$var_0 <- var_0
    genotypeParams$alphaKHyper <- rep(5000, K)
    varHyperHigh <- 15000
    highStates <- c(11:K)
    genotypeParams$alphaKHyper[highStates] <- varHyperHigh  #AMP(11-15),HLAMP(16-21) states
    genotypeParams$betaKHyper <- rep(25, K)
    genotypeParams$kappaGHyper <- kappaGHyper
    genotypeParams$outlierVar <- 10000
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


loadAlleleCountsFromFile <- function(infile, symmetric = TRUE, 
    sep = "\t") {
    options(stringsAsFactors = FALSE)
    #### LOAD INPUT READ COUNT DATA ####
    message("titan: Loading data ", infile)
    data <- read.delim(infile, header = TRUE, stringsAsFactors = FALSE, 
        sep = sep)
    chr <- data[, 1]
    chr[chr == "X"] <- "23"
    chr[chr == "Y"] <- "24"
    chr[chr == "MT"] <- "25"
    chr[chr == "M"] <- "25"
    chr <- as.numeric(chr)
    posn <- as.numeric(data[, 2])
    refOriginal <- as.numeric(data[, 4])
    nonRef <- as.numeric(data[, 6])
    tumDepth <- refOriginal + nonRef
    if (symmetric) {
        ref <- apply(cbind(refOriginal, nonRef), 1, 
            max, na.rm = TRUE)
    } else {
        ref <- refOriginal
    }
    if (ncol(data) == 7) {
        normDepth <- as.numeric(data[, 7])
    } else {
        normDepth <- NULL
    }
    return(list(chr = chr, posn = posn, ref = ref, 
        refOriginal = refOriginal, nonRef = nonRef, 
        tumDepth = tumDepth, normDepth = normDepth))
}


## filter data by depth and mappability.  cI is a
## logical vector data is a list containing all our
## data: ref, depth, logR, etc.  assumes that data
## is pass by value assumes chromosome is numeric
## and ranges from 1:25 where 23=X, 24=Y, 25=MT
filterData <- function(data, chrs = 1:24, minDepth = 10, 
    maxDepth = 200, positionList = NULL, map = NULL, 
    mapThres = 0.9) {
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
    keepChrs <- is.element(data$chr, chrs)
    cI <- keepChrs & keepTumDepth & !is.na(data$logR) & 
        keepMap & keepPosn
    for (i in 1:length(data)) {
        if (!is.null(data[[i]])) {
            data[[i]] <- data[[i]][cI]
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

getPositionOverlap <- function(chr, posn, cnData) {
    cnChr <- cnData[, 1]
    cnChr[cnChr == "X"] <- "23"
    cnChr[cnChr == "Y"] <- "24"
    cnChr[cnChr == "MT"] <- "25"
    cnChr[cnChr == "M"] <- "25"
    cnChr <- as.numeric(cnChr)
    cnStart <- as.numeric(cnData[, 2])
    cnStop <- as.numeric(cnData[, 3])
    cnVal <- as.numeric(cnData[, 4])
    
    N = length(posn)
    valByPosn = rep(NA, N)
    
    for (c in unique(chr)) {
        indData <- chr == c
        indCN <- cnChr == c
        cnStartC <- cnStart[indCN]
        cnStopC <- cnStop[indCN]
        cnValC <- cnVal[indCN]
        posnC <- posn[indData]
        cnInd <- .Call("getPositionOverlapC", posnC, 
            cnStartC, cnStopC)
        if (sum(cnInd > 0) > 0) {
            cnValToUse <- rep(NA, length(cnInd))
            cnValToUse[which(cnInd > 0)] <- cnValC[cnInd]
            valByPosn[indData] <- cnValToUse
        }
    }
    
    return(as.numeric(valByPosn))
}

### Subroutine for correcting read counts ###
correctGCMapBias <- function(x, gc, map, mappability = 0.9, 
    samplesize = 50000, verbose = TRUE) {
    
    if (verbose) {
        message("Applying filter on data...")
    }
    valid <- !logical(length = length(x))
    valid[x <= 0 | gc < 0] <- FALSE
    ideal <- !logical(length = length(x))
    routlier <- 0.01
    range <- quantile(x[valid], prob = c(0, 1 - routlier), 
        na.rm = TRUE)
    doutlier <- 0.001
    domain <- quantile(gc[valid], prob = c(doutlier, 
        1 - doutlier), na.rm = TRUE)
    ideal[!valid | map < mappability | x <= range[1] | 
        x > range[2] | gc < domain[1] | gc > domain[2]] <- FALSE
    
    if (verbose) {
        message("Correcting for GC bias...")
    }
    set <- which(ideal)
    select <- sample(set, min(length(set), samplesize))
    rough = loess(x[select] ~ gc[select], span = 0.03)
    i <- seq(0, 1, by = 0.001)
    final = loess(predict(rough, i) ~ i, span = 0.3)
    cor.gc <- x/predict(final, gc)
    
    if (verbose) {
        message("Correcting for mappability bias...")
    }
    coutlier <- 0.01
    range <- quantile(cor.gc[which(valid)], prob = c(0, 
        1 - coutlier), na.rm = TRUE)
    set <- which(cor.gc < range[2])
    select <- sample(set, min(length(set), samplesize))
    final = approxfun(lowess(map[select], cor.gc[select]))
    y <- cor.gc/final(map)
    
    return(y)
}

correctReadDepth <- function(tumWig, normWig, gcWig, 
    mapWig, targetedSequence = NULL) {
    # require(HMMcopy)
    
    message("Reading GC and mappability files")
    gc <- wigToRangedData(gcWig)
    map <- wigToRangedData(mapWig)
    
    ### LOAD TUMOUR AND NORMAL FILES ###
    message("Loading tumour file:", tumWig)
    tumour_reads <- wigToRangedData(tumWig)
    message("Loading normal file:", normWig)
    normal_reads <- wigToRangedData(normWig)
    
    ### remove 'chr' or 'CHR' from chromosome
    
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
        keepInd <- unlist(as.list(findOverlaps(tumour_reads, 
            targetIR, select = "first")))
        keepInd <- !is.na(keepInd)
        
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
    options(stringsAsFactors = FALSE)
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

computeSDbwIndex <- function(x, method = "median") {
    ## input x: Titan results dataframe from
    ## 'outputTitanResults()' S_Dbw Validity Index
    ## Halkidi and Vazirgiannis (2001). Clustering
    ## Validity Assessment: Finding the Optimal
    ## Partitional of a Data Set
    
    ## flatten copynumber-clonalclusters to single
    ## vector
    cn <- as.numeric(x[, "CopyNumber"])
    flatState <- (as.numeric(x[, "ClonalCluster"]) - 
        1) * (max(cn, na.rm = TRUE) + 1) + cn
    flatState[is.na(flatState)] <- 2
    x <- as.matrix(cbind(as.numeric(flatState), as.numeric(x[, "LogRatio"])))
    
    clust <- sort(unique(x[, 1]))
    K <- length(clust)
    N <- nrow(x)
    
    ## find average standard deviation and scatter
    ## (compactness)
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
        scat.Ci[i] <- var.Ci/var.D
    }
    avgStdev <- sqrt(sum(stdev, na.rm = TRUE))/K
    
    ## compute density between clusters (separation)
    sumDensityDiff <- matrix(NA, nrow = K, ncol = K)
    for (i in 1:K) {
        # cat('Calculating S_Dbw for cluster
        # ',clust[i],'\n')
        ind.i <- x[, 1] == clust[i]
        xci <- x[ind.i, 2]
        if (method == "median") {
            ci <- median(xci, na.rm = TRUE)  #centroid of cluster Ci
        } else if (method == "mean") {
            ci <- mean(xci, na.rm = TRUE)  #centroid of cluster Ci
        }
        sumDiff.xci <- sum(abs(xci - ci) <= avgStdev, 
            na.rm = TRUE)  #density of Ci
        # sumDiff.xci <- sum(abs(xci-ci),na.rm=TRUE)
        
        for (j in 1:K) {
            if (i == j) {
                next
            }
            ind.j <- x[, 1] == clust[j]
            xcj <- x[ind.j, 2]
            if (method == "median") {
                cj <- median(xcj, na.rm = TRUE)  #centroid of cluster Cj
            } else if (method == "mean") {
                cj <- mean(xcj, na.rm = TRUE)  #centroid of cluster Cj
            }
            sumDiff.xcj <- sum(abs(xcj - cj) <= avgStdev, 
                na.rm = TRUE)  #density of Cj
            # sumDiff.xcj <- sum(abs(xcj-cj),na.rm=TRUE)
            
            ## union of both clusters
            x.ci.cj <- union(xci, xcj)
            cij <- (ci + cj)/2
            if (method == "median") {
                cij <- median(x.ci.cj, na.rm = TRUE)
            } else if (method == "mean") {
                cij <- mean(x.ci.cj, na.rm = TRUE)
            }
            sumDiff.xci.xcj <- sum(abs(x.ci.cj - cij) <= 
                avgStdev, na.rm = TRUE)  #density of mid-point
            # sumDiff.xci.xcj <-
            # sum(abs(x.ci.cj-cij),na.rm=TRUE)
            maxDiff <- max(sumDiff.xci, sumDiff.xcj)
            if (maxDiff == 0) {
                maxDiff <- 0.1
            }
            sumDensityDiff[i, j] <- sumDiff.xci.xcj/maxDiff
        }
    }
    scat <- sum(scat.Ci, na.rm = TRUE)/(K)
    dens.bw <- sum(sumDensityDiff, na.rm = TRUE)/(K * 
        (K - 1))
    S_DbwIndex <- scat + dens.bw
    # return(S_DbwIndex)
    return(list(S_DbwIndex = S_DbwIndex, dens.bw = dens.bw, 
        scat = scat))
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
        CN[G >= 20 & G <= 23] <- 8
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
    optimalPath, filename = NULL, posteriorProbs = FALSE) {
    
    # check if useOutlierState is in convergeParams
    if (length(convergeParams$useOutlierState) == 0) {
        stop("convergeParams does not contain element: useOutlierState.")
    }
    useOutlierState <- convergeParams$useOutlierState
    
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
    rhoG <- t(convergeParams$rhoG)
    rhoZ <- t(convergeParams$rhoZ)
    
    ### OUTPUT RESULTS #### Output Z ##
    Gdecode <- decodeLOH(G)
    Gcalls <- Gdecode$G
    CN <- Gdecode$CN
    rhoZ <- rhoZ[, sortS$ix, drop = FALSE]
    Zclust[Gcalls == "HET" & CN == 2] <- NA  #diploid HET positions do not have clusters
    Sout <- rep(NA, length(Zclust))  #output cluster frequencies
    Sout[Zclust > 0 & !is.na(Zclust)] <- s[Zclust[Zclust > 
        0 & !is.na(Zclust)]]
    Sout[Zclust == 0] <- 0
    data$chr[data$chr == 23] <- "X"
    data$chr[data$chr == 24] <- "Y"
    data$chr[data$chr == 25] <- "MT"
    clonalHeaderStr <- rep(NA, Z)
    for (j in 1:Z) {
        clonalHeaderStr[j] <- sprintf("pClust%d", j)
    }
    colnames(rhoZ) <- clonalHeaderStr
    outmat <- cbind(Chr = data$chr, Position = data$posn, 
        RefCount = data$refOriginal, NRefCount = data$tumDepth - 
            data$refOriginal, Depth = data$tumDepth, 
        AllelicRatio = sprintf("%0.2f", data$refOriginal/data$tumDepth), 
        LogRatio = sprintf("%0.2f", log2(exp(data$logR))), 
        CopyNumber = CN, TITANstate = G, TITANcall = Gcalls, 
        ClonalCluster = Zclust, CellularPrevalence = sprintf("%0.2f", 
            1 - Sout))
    if (posteriorProbs) {
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

outputModelParameters <- function(convergeParams, results, 
    filename) {
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
        row.names = FALSE, quote = FALSE, sep = "", 
        append = TRUE)
    ploid_str <- sprintf("Average tumour ploidy estimate:\t%0.2f", 
        convergeParams$phi[i])
    write.table(ploid_str, file = fc, col.names = FALSE, 
        row.names = FALSE, quote = FALSE, sep = "", 
        append = TRUE)
    s_str <- sprintf("%0.4f ", 1 - s)
    s_str <- gsub(" ", "", s_str)
    outStr <- sprintf("Clonal cluster cellular prevalence Z=%d:\t%s", 
        Z, paste(s_str, collapse = " "))
    write.table(outStr, file = fc, col.names = FALSE, 
        row.names = FALSE, quote = FALSE, sep = "", 
        append = TRUE)
    for (j in 1:Z) {
        musR_str <- sprintf("%0.2f ", convergeParams$muR[, 
            j, i])
        musR_str <- gsub(" ", "", musR_str)
        outStr <- sprintf("Genotype binomial means for clonal cluster 
                  Z=%d:\t%s", j, paste(musR_str, collapse = " "))
        write.table(outStr, file = fc, col.names = FALSE, 
            row.names = FALSE, quote = FALSE, sep = "", 
            append = TRUE)
        musC_str <- sprintf("%0.2f ", log2(exp(convergeParams$muC[, j, i])))
        musC_str <- gsub(" ", "", musC_str)
        outStr <- sprintf("Genotype Gaussian means for clonal cluster 
            Z=%d:\t%s", j, paste(musC_str, collapse = " "))
        write.table(outStr, file = fc, col.names = FALSE, 
            row.names = FALSE, quote = FALSE, sep = "", 
            append = TRUE)
    }
    var_str <- sprintf("%0.4f ", convergeParams$var[, 
        i])
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
    sdbw <- computeSDbwIndex(results, method = "median")
    sdbw_str <- sprintf("S_Dbw dens.bw:\t%0.4f ", sdbw$dens.bw)
    write.table(sdbw_str, file = fc, col.names = FALSE, 
        row.names = FALSE, quote = FALSE, sep = "", 
        append = TRUE)
    sdbw_str <- sprintf("S_Dbw scat:\t%0.4f ", sdbw$scat)
    write.table(sdbw_str, file = fc, col.names = FALSE, 
        row.names = FALSE, quote = FALSE, sep = "", 
        append = TRUE)
    sdbw_str <- sprintf("S_Dbw validity index:\t%0.4f ", 
        25 * sdbw$dens.bw + sdbw$scat)
    write.table(sdbw_str, file = fc, col.names = FALSE, 
        row.names = FALSE, quote = FALSE, sep = "", 
        append = TRUE)
    
    close(fc)
    
} 
