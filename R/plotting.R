# file : plotting.R author: Gavin Ha <gha@bccrc.ca>
# Dept of Molecular Oncolgy British Columbia Cancer
# Agency University of British Columbia date :
# November 7, 2012

# data is the output format of TITAN cytoBand = {T,
# F} alphaVal = [0,1] geneAnnot is a dataframe with
# 4 columns: geneSymbol, chr, start, stop spacing
# is the distance between each track
plotAllelicRatio <- function(dataIn, chr = NULL, geneAnnot = NULL, 
    spacing = 4, xlim = NULL, ...) {
    # color coding alphaVal <- ceiling(alphaVal * 255);
    # class(alphaVal) = 'hexmode'
    lohCol <- c("#00FF00", "#006400", "#0000FF", "#8B0000", 
        "#006400", "#BEBEBE", "#FF0000", "#BEBEBE", 
        "#FF0000")
    # lohCol <- paste(lohCol,alphaVal,sep='') lohCol <-
    # col2rgb(c('green','darkgreen','blue','darkgreen','grey','red'))
    names(lohCol) <- c("HOMD", "DLOH", "NLOH", "GAIN", 
        "ALOH", "HET", "ASCNA", "BCNA", "UBCNA")
    
    
    if (!is.null(chr)) {
        for (i in chr) {
            dataByChr <- dataIn[dataIn[, "Chr"] == 
                i, ]
            dataByChr <- dataByChr[dataByChr[, "TITANcall"] != 
                "OUT", ]
            # plot the data if (outfile!=''){
            # pdf(outfile,width=10,height=6) }
            par(mar = c(spacing, 8, 2, 2))
            # par(xpd=NA)
            if (is.null(xlim)) {
                xlim <- as.numeric(c(1, dataByChr[nrow(dataByChr), 
                  "Position"]))
            }
            plot(dataByChr[, "Position"], dataByChr[, 
                "AllelicRatio"], col = lohCol[dataByChr[, 
                "TITANcall"]], pch = 16, xaxt = "n", 
                las = 1, ylab = "Allelic Ratio", xlim = xlim, 
                ...)
            lines(as.numeric(c(1, dataByChr[nrow(dataByChr), 
                "Position"])), rep(0.5, 2), type = "l", 
                col = "grey", lwd = 3)
            
            if (!is.null(geneAnnot)) {
                plotGeneAnnotation(geneAnnot, i)
            }
        }
    } else {
        # plot for all chromosomes
        coord <- getGenomeWidePositions(dataIn[, "Chr"], 
            dataIn[, "Position"])
        plot(coord$posns, as.numeric(dataIn[, "AllelicRatio"]), 
            col = lohCol[dataIn[, "TITANcall"]], pch = 16, 
            xaxt = "n", bty = "n", las = 1, ylab = "Allelic Ratio", 
            ...)
        lines(as.numeric(c(1, coord$posns[length(coord$posns)])), 
            rep(0.5, 2), type = "l", col = "grey", 
            lwd = 3)
        plotChrLines(unique(dataIn[, "Chr"]), coord$chrBkpt, 
            c(-0.1, 1.1))
        
    }
}

# data is the output format of TITAN alphaVal =
# [0,1] geneAnnot is a dataframe with 4 columns:
# geneSymbol, chr, start, stop spacing is the
# distance between each track
plotClonalFrequency <- function(dataIn, chr = NULL, 
    normal = NULL, geneAnnot = NULL, spacing = 4, xlim = NULL, 
    ...) {
    # color coding
    lohCol <- c("#00FF00", "#006400", "#0000FF", "#8B0000", 
        "#006400", "#BEBEBE", "#FF0000", "#FF0000", 
        "#FF0000")
    names(lohCol) <- c("HOMD", "DLOH", "NLOH", "GAIN", 
        "ALOH", "HET", "ASCNA", "BCNA", "UBCNA")
    
    # get unique set of cluster and estimates table:
    # 1st column is cluster number, 2nd column is
    # clonal freq
    clusters <- unique(dataIn[, c("ClonalCluster", 
        "CellularPrevalence")])
    clusters <- clusters[!is.na(clusters[, 1]), , drop = FALSE]  #exclude NA
    if (!is.null(normal)) {
        clusters[, 2] <- (as.numeric(clusters[, 2])) * 
            (1 - as.numeric(normal))
    }
    
    dataToUse <- dataIn[dataIn[, "TITANcall"] != "OUT", 
        ]
    dataToUse[dataToUse[, "CellularPrevalence"] == 
        "NA" | is.na(dataToUse[, "CellularPrevalence"]), 
        c("ClonalCluster", "CellularPrevalence")] <- c(NA, 
        NA)
    # extract clonal info
    clonalFreq <- cbind(as.numeric(dataToUse[, "ClonalCluster"]), 
        as.numeric(dataToUse[, "CellularPrevalence"]))
    # mode(clonalFreq) <- 'numeric' clonalFreq[,2] <- 1
    # - clonalFreq[,2]
    if (!is.null(normal)) {
        clonalFreq[, 2] <- clonalFreq[, 2] * (1 - normal)
    }
    clonalFreq[is.na(clonalFreq[, 2]) | clonalFreq[, 
        2] == "0" | clonalFreq[, 2] == "NA", 2] <- 0
    
    # plot per chromosome
    if (!is.null(chr)) {
        for (i in chr) {
            ind <- dataToUse[, "Chr"] == as.character(i)
            dataByChr <- dataToUse[ind, ]
            clonalFreq <- clonalFreq[ind, ]
            # plot the data
            par(mar = c(spacing, 8, 2, 2), xpd = NA)
            # par(xpd=NA)
            
            # PLOT CLONAL FREQUENCIES
            if (is.null(xlim)) {
                xlim <- as.numeric(c(1, dataByChr[nrow(dataByChr), 
                  "Position"]))
            }
            plot(dataByChr[, "Position"], clonalFreq[, 
                2], type = "h", col = lohCol[dataByChr[, 
                "TITANcall"]], las = 1, xaxt = "n", 
                ylab = "Cellular Prevalence", xlim = xlim, 
                ...)
            
            # plot cluster lines and labels
            for (j in 1:length(clusters[, 1])) {
                chrLen <- as.numeric(dataByChr[dim(dataByChr)[1], 
                  "Position"])
                lines(c(1 - chrLen * 0.02, chrLen * 
                  1.02), rep(clusters[j, 2], 2), type = "l", 
                  col = "grey", lwd = 3)
                mtext(side = 4, at = clusters[j, 2], 
                  text = paste("Z", clusters[j, 1], 
                    "", sep = ""), cex = 1, padj = 0.5, 
                  adj = 1, las = 2, outer = FALSE)
                mtext(side = 2, at = clusters[j, 2], 
                  text = paste("Z", clusters[j, 1], 
                    "", sep = ""), cex = 1, padj = 0.5, 
                  adj = 0, las = 2, outer = FALSE)
            }
            
            if (!is.null(normal)) {
                chrLen <- as.numeric(dataByChr[nrow(dataByChr), 
                  "Position"])
                lines(c(1 - chrLen * 0.02, chrLen * 
                  1.02), rep((1 - normal), 2), type = "l", 
                  col = "#000000", lwd = 3)
                mtext(side = 4, at = (1 - normal), 
                  text = paste("-T-", sep = ""), padj = 0.5, 
                  adj = 1, cex = 1, las = 2, outer = FALSE)
                mtext(side = 2, at = (1 - normal), 
                  text = paste("-T-", sep = ""), padj = 0.5, 
                  adj = 0, cex = 1, las = 2, outer = FALSE)
            }
            
            if (!is.null(geneAnnot)) {
                plotGeneAnnotation(geneAnnot, i)
            }
        }
    } else {
        # plot genome-wide
        coord <- getGenomeWidePositions(dataIn[, "Chr"], 
            dataIn[, "Position"])
        plot(coord$posns, clonalFreq[, 2], type = "h", 
            col = lohCol[dataIn[, "TITANcall"]], pch = 16, 
            xaxt = "n", las = 1, bty = "n", ylab = "Cellular Prevalence", 
            ...)
        plotChrLines(unique(dataIn[, "Chr"]), coord$chrBkpt, 
            c(-0.1, 1.1))
        
        # plot cluster lines and labels
        for (j in 1:length(clusters[, 1])) {
            chrLen <- as.numeric(coord$posns[length(coord$posns)])
            lines(c(1 - chrLen * 0.02, chrLen * 1.02), 
                rep(clusters[j, 2], 2), type = "l", 
                col = "grey", lwd = 3)
            mtext(side = 4, at = clusters[j, 2], text = paste("Z", 
                clusters[j, 1], "", sep = ""), cex = 1, 
                padj = 0.5, adj = 1, las = 2, outer = FALSE)
            mtext(side = 2, at = clusters[j, 2], text = paste("Z", 
                clusters[j, 1], "", sep = ""), cex = 1, 
                padj = 0.5, adj = 0, las = 2, outer = FALSE)
        }
        if (!is.null(normal)) {
            chrLen <- as.numeric(coord$posns[length(coord$posns)])
            lines(c(1 - chrLen * 0.02, chrLen * 1.02), 
                rep((1 - normal), 2), type = "l", col = "#000000", 
                lwd = 3)
        }
        
    }
    
}



# data is the output format of TITAN (*loh.txt)
# alphaVal = [0,1] geneAnnot is a dataframe with 4
# columns: geneSymbol, chr, start, stop spacing is
# the distance between each track
plotCNlogRByChr <- function(dataIn, chr = NULL, geneAnnot = NULL, 
    ploidy = NULL, spacing = 4, alphaVal = 1, xlim = NULL, 
    ...) {
    # color coding
    alphaVal <- ceiling(alphaVal * 255)
    class(alphaVal) = "hexmode"
    cnCol <- c("#00FF00", "#006400", "#0000FF", "#880000", 
        "#BB0000", "#CC0000", "#DD0000", "#EE0000", 
        "#FF0000")
    cnCol <- paste(cnCol, alphaVal, sep = "")
    # cnCol <-
    # col2rgb(c('green','darkgreen','blue','darkred','red','brightred'))
    names(cnCol) <- c("0", "1", "2", "3", "4", "5", 
        "6", "7", "8")
    
    ## adjust logR values for ploidy ##
    if (!is.null(ploidy)) {
        dataIn[, "LogRatio"] <- as.numeric(dataIn[, 
            "LogRatio"]) + log2(ploidy/2)
    }
    
    if (!is.null(chr)) {
        for (i in chr) {
            dataByChr <- dataIn[dataIn[, "Chr"] == 
                i, ]
            dataByChr <- dataByChr[dataByChr[, "TITANcall"] != 
                "OUT", ]
            # plot the data if (outfile!=''){
            # pdf(outfile,width=10,height=6) }
            par(mar = c(spacing, 8, 2, 2))
            # par(xpd=NA)
            if (is.null(xlim)) {
                xlim <- as.numeric(c(1, dataByChr[nrow(dataByChr), 
                  "Position"]))
            }
            coord <- as.numeric(dataByChr[, "Position"])
            plot(coord, as.numeric(dataByChr[, "LogRatio"]), 
                col = cnCol[as.character(dataByChr[, 
                  "CopyNumber"])], pch = 16, xaxt = "n", 
                las = 1, ylab = "Copy Number (log ratio)", 
                xlim = xlim, ...)
            lines(c(1, as.numeric(dataByChr[dim(dataByChr)[1], 
                "Position"])), rep(0, 2), type = "l", 
                col = "grey", lwd = 0.75)
            
            if (!is.null(geneAnnot)) {
                plotGeneAnnotation(geneAnnot, i)
            }
        }
    } else {
        # plot for all chromosomes
        coord <- getGenomeWidePositions(dataIn[, "Chr"], 
            dataIn[, "Position"])
        plot(coord$posns, as.numeric(dataIn[, "LogRatio"]), 
            col = cnCol[as.character(dataIn[, "CopyNumber"])], 
            pch = 16, xaxt = "n", las = 1, bty = "n", 
            ylab = "Copy Number (log ratio)", ...)
        lines(as.numeric(c(1, coord$posns[length(coord$posns)])), 
            rep(0, 2), type = "l", col = "grey", lwd = 2)
        plotChrLines(dataIn[, "Chr"], coord$chrBkpt, 
            par("yaxp")[1:2])
    }
    
}


plotGeneAnnotation <- function(geneAnnot, chr = 1, 
    ...) {
    colnames(geneAnnot) <- c("Gene", "Chr", "Start", 
        "Stop")
    geneAnnot <- geneAnnot[geneAnnot[, "Chr"] == as.character(chr), 
        ]
    if (nrow(geneAnnot) != 0) {
        for (g in 1:dim(geneAnnot)[1]) {
            # print(geneAnnot[g,'Gene'])
            abline(v = as.numeric(geneAnnot[g, "Start"]), 
                col = "black", lty = 3, xpd = FALSE)
            abline(v = as.numeric(geneAnnot[g, "Stop"]), 
                col = "black", lty = 3, xpd = FALSE)
            atP <- (as.numeric(geneAnnot[g, "Stop"]) - 
                as.numeric(geneAnnot[g, "Start"]))/2 + 
                as.numeric(geneAnnot[g, "Start"])
            # if (atP < dataByChr[1,2]){ atP <- dataByChr[1,2]
            # }else if (atP > dataByChr[dim(dataByChr)[1],2]){
            # atP <- dataByChr[dim(dataByChr)[1],2] }
            mtext(geneAnnot[g, "Gene"], side = 3, line = 0, 
                at = atP, ...)
        }
    }
}

plotChrLines <- function(chrs, chrBkpt, yrange) {
    # plot vertical chromosome lines
    for (j in 1:length(chrBkpt)) {
        lines(rep(chrBkpt[j], 2), yrange, type = "l", 
            lty = 2, col = "black", lwd = 0.75)
    }
    numLines <- length(chrBkpt)
    mid <- (chrBkpt[1:(numLines - 1)] + chrBkpt[2:numLines])/2
    chrs[chrs == "X"] <- 23
    chrs[chrs == "Y"] <- 24
    chrsToShow <- sort(unique(as.numeric(chrs)))
    chrsToShow[chrsToShow == 23] <- "X"
    chrsToShow[chrsToShow == 24] <- "Y"
    axis(side = 1, at = mid, labels = c(chrsToShow), 
        cex.axis = 1.5, tick = FALSE)
}

getGenomeWidePositions <- function(chrs, posns) {
    # create genome coordinate scaffold
    positions <- as.numeric(posns)
    chrsNum <- unique(chrs)
    chrBkpt <- rep(0, length(chrsNum) + 1)
    for (i in 2:length(chrsNum)) {
        chrInd <- which(chrs == chrsNum[i])
        prevChrPos <- positions[chrInd[1] - 1]
        chrBkpt[i] = prevChrPos
        positions[chrInd] = positions[chrInd] + prevChrPos
    }
    chrBkpt[i + 1] <- positions[length(positions)]
    return(list(posns = positions, chrBkpt = chrBkpt))
} 
