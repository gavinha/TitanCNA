# author: Gavin Ha
# 		  Dana-Farber Cancer Institute
#		  Broad Institute
# contact: <gavinha@gmail.com> or <gavinha@broadinstitute.org>
# date:	  January 11, 2018

# data is the output format of TITAN cytoBand = {T,
# F} alphaVal = [0,1] geneAnnot is a dataframe with
# 4 columns: geneSymbol, chr, start, stop spacing
# is the distance between each track
plotAllelicRatio <- function(dataIn, chr = NULL, geneAnnot = NULL,
    spacing = 4,  xlim = NULL, ...) {
    # color coding alphaVal <- ceiling(alphaVal * 255);
    # class(alphaVal) = 'hexmode'
    lohCol <- c("#00FF00", "#006400", "#0000FF", "#8B0000",
        "#006400", "#BEBEBE", "#FF0000", "#BEBEBE",
        "#FF0000")
    # lohCol <- paste(lohCol,alphaVal,sep='') lohCol <-
    # col2rgb(c('green','darkgreen','blue','darkgreen','grey','red'))
    names(lohCol) <- c("HOMD", "DLOH", "NLOH", "GAIN",
        "ALOH", "HET", "ASCNA", "BCNA", "UBCNA")

    dataIn <- copy(dataIn)
    if (!is.null(chr)) {
        for (i in chr) {
            dataByChr <- dataIn[Chr == i]
            dataByChr <- dataByChr[TITANcall != "OUT"]
            # plot the data if (outfile!=''){
            # pdf(outfile,width=10,height=6) }
            par(mar = c(spacing, 8, 2, 2))
            # par(xpd=NA)
            if (missing(xlim)) {
                xlim <- as.numeric(c(1, dataByChr[nrow(dataByChr), Position]))
            }
            plot(dataByChr[, Position], dataByChr[,
                AllelicRatio], col = lohCol[dataByChr[,
                TITANcall]], pch = 16, xaxt = "n",
                las = 1, ylab = "Allelic Ratio", xlim = xlim,
                ...)
            lines(as.numeric(c(1, dataByChr[nrow(dataByChr),
                Position])), rep(0.5, 2), type = "l",
                col = "grey", lwd = 3)

            if (!is.null(geneAnnot)) {
                plotGeneAnnotation(geneAnnot, i)
            }
        }
    } else {
        # plot for all chromosomes
        coord <- getGenomeWidePositions(dataIn[, Chr],
            dataIn[, Position])
        plot(coord$posns, as.numeric(dataIn[, AllelicRatio]),
            col = lohCol[dataIn[, TITANcall]], pch = 16,
            xaxt = "n", bty = "n", las = 1, ylab = "Allelic Fraction",
            ...)
        lines(as.numeric(c(1, coord$posns[length(coord$posns)])),
            rep(0.5, 2), type = "l", col = "grey",
            lwd = 3)
        plotChrLines(unique(dataIn[, Chr]), coord$chrBkpt,
            c(-0.1, 1.1))

    }
}

# data is the output format of TITAN alphaVal =
# [0,1] geneAnnot is a dataframe with 4 columns:
# geneSymbol, chr, start, stop spacing is the
# distance between each track
plotClonalFrequency <- function(dataIn, chr = NULL,
    normal = NULL, geneAnnot = NULL, spacing = 4, xlim = NULL, ...) {
    # color coding
    lohCol <- c("#00FF00", "#006400", "#0000FF", "#8B0000",
        "#006400", "#BEBEBE", "#FF0000", "#FF0000",
        "#FF0000", "#FF0000", "#FF0000")
    names(lohCol) <- c("HOMD", "DLOH", "NLOH", "GAIN",
        "ALOH", "HET", "ASCNA", "BCNA", "UBCNA", "AMP", "HLAMP")

    # get unique set of cluster and estimates table:
    # 1st column is cluster number, 2nd column is
    # clonal freq
    clusters <- unique(dataIn[, list(ClonalCluster, CellularPrevalence)])
    clusters <- clusters[!is.na(ClonalCluster), ]  #exclude NA
    if (!is.null(normal)) {
        clusters[, CellularPrevalence := CellularPrevalence * (1 - as.numeric(normal))]
    }
    dataToUse <- copy(dataIn)
    dataToUse <- dataToUse[TITANcall != "OUT", ]
    #dataToUse[CellularPrevalence == "NA" | is.na(CellularPrevalence),
    #          list(ClonalCluster, CellularPrevalence) := c(NA, NA)]
    # extract clonal info
    clonalFreq <- dataToUse[, list(ClonalCluster, CellularPrevalence)]
    # mode(clonalFreq) <- 'numeric' clonalFreq[,2] <- 1 - clonalFreq[,2]
    if (!is.null(normal)) {
        clonalFreq[, CellularPrevalence := CellularPrevalence * (1 - normal)]
    }
    clonalFreq[is.na(CellularPrevalence) | CellularPrevalence == "0" |
                 CellularPrevalence == "NA", CellularPrevalence := 0]

    # plot per chromosome
    if (!is.null(chr)) {
        for (i in chr) {
            ind <- dataToUse[, Chr] == as.character(i)
            dataByChr <- dataToUse[ind, ]
            clonalFreq <- clonalFreq[ind, ]
            # plot the data
            par(mar = c(spacing, 8, 2, 2), xpd = NA)
            # par(xpd=NA)

            # PLOT CLONAL FREQUENCIES
            if (missing(xlim)) {
                xlim <- as.numeric(c(1, dataByChr[nrow(dataByChr), Position]))
            }
            plot(dataByChr[, Position], clonalFreq[, CellularPrevalence], type = "h",
                 col = lohCol[dataByChr[, TITANcall]], las = 1, xaxt = "n",
                ylab = "Cellular Prevalence", xlim = xlim, ...)

            # plot cluster lines and labels
            if (nrow(clusters) > 0){
                for (row in 1:nrow(clusters)) {
                    prevalence <- clusters[row, CellularPrevalence]
                    clustnum <- clusters[row, ClonalCluster]
                    chrLen <- as.numeric(dataByChr[dim(dataByChr)[1], Position])
                    lines(c(1 - chrLen * 0.02, chrLen * 1.02),
                          rep(prevalence , 2), type = "l", col = "grey", lwd = 3)
                    mtext(side = 4, at = prevalence,
                          text = paste("Z", clustnum, "", sep = ""),
                          cex = 1, padj = 0.5, adj = 1, las = 2, outer = FALSE)
                    mtext(side = 2, at = prevalence,
                          text = paste("Z", clustnum, "", sep = ""),
                          cex = 1, padj = 0.5,
                          adj = 0, las = 2, outer = FALSE)
                }
            }

            if (!is.null(normal)) {
                chrLen <- as.numeric(dataByChr[nrow(dataByChr), Position])
                lines(c(1 - chrLen * 0.02, chrLen *
                  1.02), rep((1 - normal), 2), type = "l",
                  col = "#000000", lwd = 3)
                #mtext(side = 4, at = (1 - normal),
                  #text = paste("-T-", sep = ""), padj = 0.5,
                  #adj = 1, cex = 1, las = 2, outer = FALSE)
                #mtext(side = 2, at = (1 - normal),
                  #text = paste("-T-", sep = ""), padj = 0.5,
                  #adj = 0, cex = 1, las = 2, outer = FALSE)
            }

            if (!is.null(geneAnnot)) {
                plotGeneAnnotation(geneAnnot, i)
            }
        }
    } else {
        # plot genome-wide
        coord <- getGenomeWidePositions(dataIn[, Chr],
            dataIn[, Position])
        plot(coord$posns, clonalFreq[, CellularPrevalence], type = "h",
            col = lohCol[dataIn[, TITANcall]], pch = 16,
            xaxt = "n", las = 1, bty = "n", ylab = "Cellular Prevalence",
            ...)
        plotChrLines(unique(dataIn[, Chr]), coord$chrBkpt,
            c(-0.1, 1.1))

        # plot cluster lines and labels
      	if (nrow(clusters) > 0){
             for (row in 1:nrow(clusters)) {
                prevalence <- clusters[row, CellularPrevalence]
                clustnum <- clusters[row, ClonalCluster]
                chrLen <- as.numeric(coord$posns[length(coord$posns)])
                lines(c(1 - chrLen * 0.02, chrLen * 1.02),
                      rep(prevalence, 2), type = "l",
                      col = "grey", lwd = 3)
                mtext(side = 4, at = prevalence, text = paste("Z",
                      clustnum, "", sep = ""), cex = 1,
                      padj = 0.5, adj = 1, las = 2, outer = FALSE)
                mtext(side = 2, at = prevalence, text = paste("Z",
                      clustnum, "", sep = ""), cex = 1,
                      padj = 0.5, adj = 0, las = 2, outer = FALSE)
            }
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
plotCNlogRByChr <- function(dataIn, chr = NULL, segs = NULL, 
	plotCorrectedCN = TRUE, geneAnnot = NULL,
    ploidy = NULL, normal = NULL, spacing = 4, alphaVal = 1, xlim = NULL, ...) {
    # color coding
    alphaVal <- ceiling(alphaVal * 255)
    class(alphaVal) = "hexmode"
    cnCol <- c("#00FF00", "#006400", "#0000FF", "#880000",
        "#BB0000", "#CC0000", "#DD0000", "#EE0000", rep("#FF0000",493))
    cnCol <- paste(cnCol, alphaVal, sep = "")
    # cnCol <-
    # col2rgb(c('green','darkgreen','blue','darkred','red','brightred'))
    names(cnCol) <- c(0:500)
	
	if (plotCorrectedCN && "Corrected_Copy_Number" %in% colnames(dataIn)){
		binCN <- "Corrected_Copy_Number"
		segCN <- "Corrected_Copy_Number"
	}else{
		binCN <- "CopyNumber"
		segCN <- "Copy_Number"
	}
	
    dataIn <- copy(dataIn)
    ## adjust logR values for ploidy ##
    if (!is.null(ploidy)) {
    	if (is.null(normal)){
    		stop("plotCNlogRByChr: Please provide \"normal\" contamination estimate.")
    	}
        dataIn[, LogRatio := LogRatio + log2(((1-normal)*ploidy+normal*2)/2)]

      if (!is.null(segs)){
        segs.sample <- copy(segs)
				segs.sample[, Median_logR := Median_logR + log2(((1-normal)*ploidy+normal*2) / 2)]
			}
    }

    if (!is.null(chr)) {
        for (i in chr) {
            dataByChr <- dataIn[dataIn[, Chr] == i, ]
            dataByChr <- dataByChr[dataByChr[, TITANcall] != "OUT", ]
            # plot the data if (outfile!=''){
            # pdf(outfile,width=10,height=6) }
            par(mar = c(spacing, 8, 2, 2))
            # par(xpd=NA)
            if (missing(xlim)) {
                xlim <- as.numeric(c(1, dataByChr[nrow(dataByChr), Position]))
            }
            coord <- as.numeric(dataByChr[, Position])
            plot(coord, as.numeric(dataByChr[, LogRatio]),
                col = cnCol[as.character(dataByChr[, get(binCN)])], pch = 16, xaxt = "n", las = 1, ylab = "Copy Number (log ratio)", xlim = xlim, ...)
            lines(xlim, rep(0, 2), type = "l", col = "grey", lwd = 0.75)
            if (!is.null(segs)){
							segsByChr <- segs.sample[Chromosome == as.character(i), ]
							tmp <- apply(segsByChr, 1, function(x){
								lines(x[c("Start_Position.bp.","End_Position.bp.")], 
										rep(x["Median_logR"], 2), col = cnCol[as.character(x[segCN])], lwd = 3, lend = 1)
							})
						}


            if (!is.null(geneAnnot)) {
                plotGeneAnnotation(geneAnnot, i)
            }
        }
    } else {
        # plot for all chromosomes
        coord <- getGenomeWidePositions(dataIn[, Chr], dataIn[, Position])
        plot(coord$posns, as.numeric(dataIn[, LogRatio]),
            col = cnCol[as.character(dataIn[, get(binCN)])],
            pch = 16, xaxt = "n", las = 1, bty = "n",
            ylab = "Copy Number (log ratio)", ...)
        lines(as.numeric(c(1, coord$posns[length(coord$posns)])),
            rep(0, 2), type = "l", col = "grey", lwd = 2)
        plotChrLines(dataIn[, Chr], coord$chrBkpt, par("yaxp")[1:2])
        #plot segments
				if (!is.null(segs)){
					coordEnd <- getGenomeWidePositions(segs.sample[, Chromosome], segs.sample[, End_Position.bp.])
					coordStart <- coordEnd$posns - (segs.sample[, End_Position.bp.] - segs.sample[, Start_Position.bp.] + 1)
					xlim <- as.numeric(c(1, coordEnd$posns[length(coordEnd$posns)]))
					col <- cnCol[as.character(segs.sample[, get(segCN)])]
					value <- as.numeric(segs.sample[, Median_logR])
					mat <- as.data.frame(cbind(coordStart, coordEnd$posns, value, col))
					rownames(mat) <- 1:nrow(mat)
					tmp <- apply(mat, 1, function(x){
						lines(x[1:2], rep(x[3], 2), col = x[4], lwd = 3, lend = 1)
					})
				}

    }

}

plotSubcloneProfiles <- function(dataIn, chr = NULL, geneAnnot = NULL,
	spacing = 4, xlim = NULL, ...){
	args <- list(...)
	lohCol <- c("#00FF00", "#006400", "#0000FF", "#8B0000",
        "#006400", "#BEBEBE", "#FF0000", "#FF0000",
        "#FF0000")
    names(lohCol) <- c("HOMD", "DLOH", "NLOH", "GAIN",
        "ALOH", "HET", "ASCNA", "BCNA", "UBCNA")

    ## pull out params from dots ##
    if (!is.null(args$cex.axis)) cex.axis <- args$cex.axis else cex.axis <- 0.75
    if (!is.null(args$cex.lab)) cex.lab <- args$cex.lab else cex.lab <- 0.75
    dataIn <- copy(dataIn)
    numClones <- sum(!is.na(unique(as.numeric(dataIn$ClonalCluster))))
    if (numClones == 0){ numClones <- 1 }
         # plot per chromosome
    if (!is.null(chr)) {
        for (i in chr) {
            ind <- dataIn[, Chr == as.character(i)]
            dataByChr <- dataIn[ind, ]

            ## find x domain #
            if (missing(xlim)) {
    			xlim <- c(1, dataByChr[.N, Position])
    		}

            # plot the data
            par(mar = c(spacing, 8, 2, 2), xpd = NA)

            # PLOT SUBCLONE PROFILES
            # setup plot to include X number of clones (numClones)
            maxCN <- dataByChr[, max(CopyNumber)] + 1
            ylim <- c(0, numClones * (maxCN + 2) - 1)
            plot(0, type = "n", xaxt = "n", ylab = "", xlab = "",
            	xlim = xlim, ylim = ylim, yaxt = "n", ...)
            axis(2, at = seq(ylim[1], ylim[2], 1), las = 1,
            	labels = rep(c(0:maxCN, "---"), numClones), cex.axis=cex.axis)
            for (i in 1:numClones){
            	val <- dataByChr[, get(paste0("Subclone", i, ".CopyNumber"))]
            	cellPrev <- suppressWarnings(unique(dataByChr[, get(paste0("Subclone", i, ".Prevalence"))]))
            	cellPrev <- cellPrev[!is.na(cellPrev)] ## remove NA prevalence, leave subclonal prev
            	if (length(cellPrev) == 0){ cellPrev <- 0.0 } ## if only NA, then assign 0 prev
            	if (i > 1){
            		# shift values up for each subclone
            		val <- val + (numClones - 1) * (maxCN + 2)
            	}
            	call <- dataByChr[, get(paste0("Subclone", i, ".TITANcall"))]
            	points(dataByChr[, Position], val, col = lohCol[call], pch = 15, ...)
            	#lines(dataIn[, Position], val, col = lohCol[call], type = "l", lwd = 3, ...)
               	mtext(text = paste0("Subclone", i, "\n", format(cellPrev, digits = 2)),
               		side = 2, las = 0, line = 3,
               		at = i * (maxCN + 2) - (maxCN + 2) / 2 - 1, cex = cex.lab)
               	chrLen <- as.numeric(dataByChr[.N, Position])
                lines(c(1 - chrLen * 0.035, chrLen *
                  1.035), rep(i * (maxCN + 2) - 1, 2), type = "l",
                  col = "black", lwd = 1.5)
            }

            if (!is.null(geneAnnot)) {
                plotGeneAnnotation(geneAnnot, i)
            }
        }
    } else {
        # plot genome-wide
        coord <- getGenomeWidePositions(dataIn[, Chr], dataIn[, Position])
        # setup plot to include X number of clones (numClones)
		maxCN <- dataIn[, max(CopyNumber)] + 1
		ylim <- c(0, numClones * (maxCN + 2) - 1)
		xlim <- as.numeric(c(1, coord$posns[length(coord$posns)]))
		plot(0, type = "n", xaxt = "n", bty = "n", ylab = "", xlim = xlim,
			ylim = ylim, yaxt = "n", ...)
		axis(2, at = seq(ylim[1], ylim[2], 1), las = 1,
			labels = rep(c(0:maxCN, "---"), numClones))
		for (i in 1:numClones){
			val <- dataIn[, get(paste0("Subclone", i, ".CopyNumber"))]
			if (i > 1){
				# shift values up for each subclone
				val <- val + (numClones - 1) * (maxCN + 2)
			}
			call <- dataIn[, get(paste0("Subclone", i, ".TITANcall"))]
			points(coord$posns, val, col = lohCol[call], pch = 15, ...)
			mtext(text = paste0("Subclone", i), side = 2, las = 0,
					line = 2, at = i * (maxCN + 2) - (maxCN + 2) / 2 - 1, cex = 0.75)
				chrLen <- xlim[2]
			lines(c(1 - chrLen * 0.035, chrLen *
			  1.035), rep(i * (maxCN + 2) - 1, 2), type = "l",
			  col = "black", lwd = 1.5)
		}
        plotChrLines(unique(dataIn[, Chr]), coord$chrBkpt, ylim)
    }

}

## TODO: Not completed ##
plotAllelicCN <- function(dataIn, resultType = "AllelicRatio",
                          chr = NULL, geneAnnot = NULL,
    ploidy = 2, spacing = 4, alphaVal = 1, xlim = NULL, ...) {
    # color coding
    alphaVal <- ceiling(alphaVal * 255)
    class(alphaVal) = "hexmode"
    cnCol <- c("#00FF00", "#006400", "#0000FF", "#880000",
        "#BB0000", "#CC0000", "#DD0000", "#EE0000",
        "#FF0000")
    cnCol <- paste(cnCol, alphaVal, sep = "")
    # cnCol <-
    # col2rgb(c('green','darkgreen','blue','darkred','red','brightred'))
    names(cnCol) <- c("0", "1", "2", "3", "4", "5", "6", "7", "8")
    dataIn <- copy(dataIn)
    ## compute allelic copy number for each
    dataIn[, Allele.1 := get(resultType) * 2^LogRatio*ploidy]
    dataIn[, Allele.2 := (1 - get(resultType)) * 2^LogRatio*ploidy]

    if (!is.null(chr)) {
        for (i in chr) {
            dataByChr <- dataIn[Chr == i, ]
            dataByChr <- dataByChr[dataByChr[, TITANcall] != "OUT", ]
            par(mar = c(spacing, 8, 2, 2))
            # par(xpd=NA)
            if (missing(xlim)) {
                xlim <- as.numeric(c(1, dataByChr[.N, Position]))
            }
            coord <- dataByChr[, Position]
            plot(coord, dataByChr[, Allele.1],
                col = cnCol[as.character(dataByChr[, CopyNumber])], pch = 16,
                xaxt = "n", las = 1, ylab = "Copy Number",
                xlim = xlim, ...)
            points(coord, dataByChr[, Allele.2],
                   col = cnCol[as.character(dataByChr[, CopyNumber])],
                   pch = 16)
            lines(xlim, rep(0, 2), type = "l", col = "grey", lwd = 0.75)

            if (!is.null(geneAnnot)) {
                plotGeneAnnotation(geneAnnot, i)
            }
        }
	}
}



plotSegmentMedians <- function(dataIn, resultType = "LogRatio",
                               plotType = "CopyNumber", plotCorrectedCN = TRUE, 
                               chr = NULL, geneAnnot = NULL, ploidy = NULL, spacing = 4, alphaVal = 1, xlim = NULL, plot.new = FALSE, lwd = 8, ...){

	## check for the possible resultType to plot ##
	if (!resultType %in% c("LogRatio", "AllelicRatio", "HaplotypeRatio")){
		stop("plotSegmentMedians: resultType must be 'LogRatio', 'AllelicRatio', or 'HaplotypeRatio'")
	}
  if (!plotType %in% c("CopyNumber", "Ratio")){
    stop("plotSegmentMedians: plotType must be 'CopyNumber' or 'Ratio'")
  }
	dataType <- c("Median_logR", "Median_Ratio", "Median_HaplotypeRatio")
	names(dataType) <- c("LogRatio", "AllelicRatio", "HaplotypeRatio")
	axisName <- c("Copy Number (log ratio)", "Allelic Ratio", "Haplotype Fraction")
	names(axisName) <- c("LogRatio", "AllelicRatio", "HaplotypeRatio")
	axisNameCN <- c("Copy Number", "Allelic Copy Number")
	names(axisNameCN) <- c("LogRatio", "AllelicRatio")
	colName <- c("Copy_Number","TITAN_call", "TITAN_call")
	names(colName) <- c("LogRatio", "AllelicRatio", "HaplotypeRatio")

	if (plotCorrectedCN && "Corrected_Copy_Number" %in% colnames(dataIn)){
		colName[1] <- "Corrected_Copy_Number"
	}

	dataIn <- copy(dataIn)
	# color coding
    alphaVal <- ceiling(alphaVal * 255)
    class(alphaVal) = "hexmode"

    if (resultType == "LogRatio"){
		cnCol <- c("#00FF00", "#006400", "#0000FF", "#880000",
        "#BB0000", "#CC0000", "#DD0000", "#EE0000", rep("#FF0000",493))
		cnCol <- paste(cnCol, alphaVal, sep = "")
		# cnCol <-
		# col2rgb(c('green','darkgreen','blue','darkred','red','brightred'))
		names(cnCol) <- c(0:500)
	}else if (resultType %in% c("AllelicRatio", "HaplotypeRatio")){
		cnCol <- c("#00FF00", "#006400", "#0000FF", "#8B0000",
        	"#006400", "#BEBEBE", "#FF0000", "#BEBEBE", "#FF0000", "#FF0000", "#FF0000")
    # lohCol <- paste(lohCol,alphaVal,sep='') lohCol <-
    # col2rgb(c('green','darkgreen','blue','darkgreen','grey','red'))
    names(cnCol) <- c("HOMD", "DLOH", "NLOH", "GAIN",
        "ALOH", "HET", "ASCNA", "BCNA", "UBCNA", "AMP", "HLAMP")
	}
    if (plotType == "CopyNumber"){
      axisName <- axisNameCN
    }
    if (is.null(ploidy)){
      ploidy <- 2
    }
    ## adjust logR values for ploidy ##
    if (resultType == "LogRatio") {
      if (plotType == "CopyNumber"){
        #dataIn[, (dataType[resultType]) := 2 ^ get(dataType[resultType]) * 2]
      }else{
        dataIn[, (dataType[resultType]) := get(dataType[resultType]) + log2(ploidy/2)]
      }
    }
    ## allelic or haplotype copy number
    if (resultType %in% c("AllelicRatio") && plotType == "CopyNumber"){
      ## compute allelic copy number for each
      #dataIn[, Allele.1 := get(dataType[resultType]) * 2 ^ Median_logR * ploidy]
      #dataIn[, Allele.2 := (1 - get(dataType[resultType])) * 2 ^ Median_logR * ploidy]
    }

    # plot for specified chromosomes #
	if (!is.null(chr)) {
    	for (i in chr) {
    		dataByChr <- dataIn[Chromosome == i, ]
        dataByChr <- dataByChr[TITAN_call != "OUT", ]
        # plot the data
        par(mar = c(spacing, 8, 2, 2))
        if (missing(xlim)) {
            xlim <- as.numeric(c(1, dataByChr[.N, End_Position.bp.]))
        }
        col <- cnCol[as.character(dataByChr[, get(colName[resultType])])]
        coord <- dataByChr[, .(Start_Position.bp., End_Position.bp.)]
        if (plotType == "CopyNumber"){
          if (resultType %in% c("AllelicRatio")){
            value <- dataByChr[, MajorCN]
            value2 <- dataByChr[, MinorCN]
          }else{
            value <- dataByChr[, Copy_Number]
          }
        }else{
          value <- dataByChr[, get(dataType[resultType])]
        }
        if (plot.new){
        	plot(0, type = "n", col = col, xaxt = "n", las = 1,
        		ylab = axisName[resultType], xlim = xlim, ...)
        }
        tmp <- apply(cbind(coord, value, col), 1, function(x){
        	lines(x[1:2], rep(x[3], 2), col = x[4], lend = 1, lwd = lwd)
        	})
        if (plotType == "CopyNumber" && resultType %in% c("AllelicRatio")){
           tmp <- apply(cbind(coord, value2, col), 1, function(x){
             lines(x[1:2], rep(x[3], 2), col = x[4], lend = 1, lwd = lwd)
             })
        }
        if (plotType == "Ratio"){
          lines(xlim, rep(0, 2), type = "l", col = "grey", lwd = 0.75)
        }
        if (!is.null(geneAnnot)) {
            plotGeneAnnotation(geneAnnot, i)
        }
    	}
    } else {
        # plot for all chromosomes
        coordEnd <- getGenomeWidePositions(dataIn[, Chromosome], dataIn[, End_Position.bp.])
    	  coordStart <- coordEnd$posns - (dataIn[, End_Position.bp.] - dataIn[, Start_Position.bp.] + 1)
        xlim <- as.numeric(c(1, coordEnd$posns[length(coordEnd$posns)]))
    	  col <- cnCol[as.character(dataIn[, get(colName[resultType])])]
    	  if (plotType == "CopyNumber"){
          if (resultType %in% c("AllelicRatio")){
            value <- dataIn[, MajorCN]
            value2 <- dataIn[, MinorCN]
          }else{
            value <- dataIn[, Copy_Number]
          }
    	  }else{
          value <- dataIn[, get(dataType[resultType])]
        }
        #mat <- data.table(cbind(coordStart, coordEnd$posns, value, col))
        #rownames(mat) <- 1:nrow(mat)
        if (plot.new){
        	plot(0, type = "n", col = col, xaxt = "n", las = 1,
           		ylab = axisName[resultType], xlim = xlim, ...)
        }
        tmp <- apply(data.table(coordStart, coordEnd$posns, value, col), 1, function(x){
        	lines(x[1:2], rep(x[3], 2), col = x[4], lend = 1, lwd = lwd)
        	})
        if (plotType == "CopyNumber" && resultType %in% c("AllelicRatio")){
           tmp <- apply(data.table(coordStart, coordEnd$posns, value2, col), 1, function(x){
             lines(x[1:2], rep(x[3], 2), col = x[4], lend = 1, lwd = lwd)
             })
        }
        if (plotType == "Ratio"){
          lines(xlim, rep(0, 2), type = "l", col = "grey", lwd = 2)
        }
        plotChrLines(dataIn[, Chromosome], coordEnd$chrBkpt, par("yaxp")[1:2])
    }
}

plotGeneAnnotation <- function(geneAnnot, chr = 1, ...) {
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
