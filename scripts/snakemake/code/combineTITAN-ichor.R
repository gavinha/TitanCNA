#' combineTITAN-ichor.R
#' author: Gavin Ha 
#' institution: Fred Hutchinson Cancer Research Center
#' contact: <gha@fredhutch.org>
#' date: July 23, 2018

#' requires R-3.3+
#' @import data.table
#' @import GenomicRanges
#' @import stringr
#' @import optparse

library(optparse)

option_list <- list(
	make_option(c("--titanSeg"), type="character", help="TitanCNA segs.txt file. Required."),
	make_option(c("--titanBin"), type="character", help="TitanCNA titan.txt file. Required."),
	make_option(c("--titanParams"), type="character", help="TitanCNA params.txt file. Required."),
  	make_option(c("--ichorSeg"), type="character", help="ichorCNA segs.txt file. Required."),
	make_option(c("--ichorBin"), type="character", help="ichorCNA cna.seg file. Required."),
	make_option(c("--ichorParams"), type="character", help="ichorCNA params.txt file. Required."),
	make_option(c("--ichorNormPanel"), type="character", help="Panel of normals; bin-level black list."),
	make_option(c("--sex"), type="character", default="female", help="female or male. Default [%default]."),
	make_option(c("--libdir"), type="character", help="TitanCNA directory path to source R files if custom changes made."),
  	make_option(c("--outSegFile"), type="character", help="New combined segment file. Required"),
  	make_option(c("--outBinFile"), type="character", help="New combined bin-level file. Required"),
  	make_option(c("--centromere"), type="character", default=NULL, help="Centromere table.")  
  )

parseobj <- OptionParser(option_list=option_list, usage = "usage: Rscript %prog [options]")
opt <- parse_args(parseobj)
print(opt)

titanSeg <- opt$titanSeg
titanBin <- opt$titanBin
titanParams <- opt$titanParams
ichorSeg <- opt$ichorSeg
ichorBin <- opt$ichorBin
ichorParams <- opt$ichorParams
ichorNormPanel <- opt$ichorNormPanel
gender <- opt$sex
outSegFile <- opt$outSegFile
outBinFile <- opt$outBinFile
centromere <- opt$centromere
libdir <- opt$libdir
outImageFile <- gsub(".seg.txt", ".RData", outSegFile)

library(TitanCNA)
library(stringr)
library(data.table)
library(GenomicRanges)

if (!is.null(libdir) && libdir != "None"){
	source(paste0(libdir, "/R/utils.R"))
}

options(stringsAsFactors=F, width=150, scipen=999)
save.image(outImageFile)
## copy number state mappings ##
ichorCNmap <- c("HOMD","HETD","NEUT","GAIN","AMP","HLAMP", rep("HLAMP", 1000))
#ichorCNmap <- list("0"="HOMD", "1"="DLOH", "2"="NEUT", "3"="GAIN", "4"="AMP", "5"="AMP")
maxichorcn <- 5

## get chromosome style
titan$Chromosome <- as.character(titan$Chromosome)
titan.cn$Chr <- as.character(titan.cn$Chr)
genomeStyle <- seqlevelsStyle(titan.cn$Chr)[1]
chrs <- c(1:22, "X")
chrXStr <- grep("X", chrs, value=TRUE)
seqlevelsStyle(chrs) <- genomeStyle

## load parameters ##
params <- read.delim(titanParams, header=F, as.is=T)
purity <- 1 - as.numeric(params[1,2])
ploidyT <- as.numeric(params[2,2])
ploidy <- purity * ploidyT + (1-purity) * 2
params.ichor <- read.delim(ichorParams, header=T, as.is=T)
homd.var <- as.numeric(strsplit(params[params[,1]=="logRatio Gaussian variance:",2], " ")[[1]][1])
homd.sd <- sqrt(homd.var)

## get gender
if (is.null(gender) || gender == "None"){
	gender <- params.ichor[3, 2]
}

## load segments 
titan <- fread(titanSeg)
if (length(titan[["Cellular_Frequency"]] == 0)){
	setnames(titan, "Cellular_Frequency", "Cellular_Prevalence")
}
ichor.segs <- fread(ichorSeg)
setnames(ichor.segs, c("ID", "chrom", "start", "end", "num.mark", "seg.median.logR", "copy.number", "call"), 
		c("Sample", "Chromosome", "Start_Position.bp.", "End_Position.bp.", 
		  "Length.snp.", "Median_logR", "Copy_Number", "TITAN_call"))

## load data points ##
titan.cn <- fread(titanBin)
titan.cn <- cbind(Sample=titan[1,Sample], titan.cn)
#titan.cn[, chr := as.character(Chr)]
id <- titan[1, Sample]
ichor.cn <- fread(ichorBin)
ichor.cn <- cbind(Sample = id, ichor.cn)
#ichor.cn[, CopyNumber := state - 1]
ichor.cn[, Position := start]
setnames(ichor.cn, c("chr", "start", paste0(id,".copy.number"), paste0(id,".event"), paste0(id,".logR"), "end"), 
		c("Chr", "Start", "CopyNumber", "TITANcall", "LogRatio", "End"))

## get bin overlap with SNPs - include ichor bins even if no SNPs overlap 
titan.gr <- titan.cn[, .(Chr, Position)]
titan.gr[, Start := Position]; titan.gr[, End := Position]
titan.gr <- as(titan.gr, "GRanges")
ichor.gr <- as(ichor.cn, "GRanges")
hits <- findOverlaps(query = titan.gr, subject = ichor.gr)
titan.cn[queryHits(hits), Start := ichor.cn[subjectHits(hits), Start]]
titan.cn[queryHits(hits), End := ichor.cn[subjectHits(hits), End]]
titan.ichor.cn <- merge(titan.cn, ichor.cn, by=c("Sample", "Chr", "Start", "End"), all=T, suffix=c("",".ichor"))
titan.ichor.cn[is.na(LogRatio), LogRatio := LogRatio.ichor] # assign ichor log ratio to missing titan SNPs
titan.ichor.cn <- titan.ichor.cn[, -c(grep("ichor", colnames(titan.ichor.cn),value=T)), with=F] 

## include HOMD from ichorCNA if missing in TitanCNA (usually when no SNPs are there)
homdLogRThres.auto <- log2((2*(1-purity)) / (2*(1-purity) + ploidyT*purity)) + homd.sd
ichor.segs[, logR_Copy_Number := logRbasedCN(Median_logR, purity, ploidy, cn=2)]
ichor.segs[, Corrected_Copy_Number := as.integer(round(logR_Copy_Number))]
ichor.segs[, Corrected_Call := ichorCNmap[Corrected_Copy_Number + 1]]
ichor.segs[, Corrected_logR := log2(logR_Copy_Number / ploidy)]
ichor.segs.homd.ind <- ichor.segs[Chromosome != chrXStr & Corrected_Copy_Number == 0 & Median_logR < homdLogRThres.auto, which = TRUE]
titan.gr <- copy(titan)
titan.gr[, Start := Start_Position.bp.]; titan.gr[, End := End_Position.bp.]
titan.gr <- as(titan.gr, "GRanges")
ichor.segs.gr <- copy(ichor.segs)
ichor.segs.gr[, Start := Start_Position.bp.]; ichor.segs.gr[, End := End_Position.bp.]
ichor.homd.gr <- as(ichor.segs.gr[ichor.segs.homd.ind], "GRanges")
#hits <- findOverlaps(query = titan.gr, subject = ichor.homd.gr)
titan.gr.combinedHomd <- sort(c(titan.gr, ichor.homd.gr))
titan.combinedHomd <- as.data.table(disjoin(titan.gr.combinedHomd, with.revmap = TRUE))
mcolData <- data.table()
for (i in 1:nrow(titan.combinedHomd)){
	ind <- titan.combinedHomd[i, revmap][[1]]
	if (length(ind) == 1){
		mcolData <- rbind(mcolData, mcols(titan.gr.combinedHomd)[ind, c(4:19,21)])
	}else{ # ind has more than one mapping index
		ind.homd.gr <- ind[which(titan.gr.combinedHomd[ind]$Corrected_Copy_Number == 0)]
		mcolData <- rbind(mcolData, mcols(titan.gr.combinedHomd)[ind.homd.gr, c(4:19,21)])
	}
}
titan.combinedHomd <- cbind(titan.combinedHomd, as.data.table(mcolData))
titan.combinedHomd <- cbind(Sample = id, titan.combinedHomd[, -c(4:6)])
setnames(titan.combinedHomd, c("seqnames", "start", "end"), c("Chromosome", "Start_Position.bp.", "End_Position.bp."))


## combine TITAN (chr1-22) and ichorCNA (chrX) segments and bin/SNP level data ##
## if male only ##
if (gender == "male"){
	cn <- rbind(titan.ichor.cn[Chr %in% chrs[1:22]], ichor.cn[Chr == chrs[grep("X", chrs)]], fill = TRUE)
	segs <- rbind(titan.combinedHomd[Chromosome %in% chrs[1:22]], ichor.segs[Chromosome == chrs[grep("X", chrs)]], fill = TRUE)
	segs[, subclone.status := NULL]
}else{
	cn <- titan.ichor.cn
	segs <- titan.combinedHomd
}

## sort column order
setnames(segs, c("Start_Position.bp.", "End_Position.bp."), c("Start", "End"))
cols <- c("Sample", "Chr", "Position", "Start", "End")
setcolorder(cn, c(cols, colnames(cn)[!colnames(cn) %in% cols]))

## get major/minor CN from segs and place in SNP/level data ##
#cn.gr <- cn[, .(Chr, Start, End)]
#cn.gr <- as(na.omit(cn.gr), "GRanges")
#segs.gr <- as(segs, "GRanges")
#hits <- findOverlaps(query = cn.gr, subject = segs.gr)
#cn[queryHits(hits), MajorCN := segs[subjectHits(hits), MajorCN]]
#cn[queryHits(hits), MinorCN := segs[subjectHits(hits), MinorCN]]
#cn[is.na(CopyNumber), CopyNumber := MajorCN + MinorCN]

## remove LogRatio for bins that do not overlap segments (i.e. remaining rows with CopyNumber == NA
# these regions that are not part of segments are usually noisy
#cn[is.na(CopyNumber), LogRatio := NA]

## filter outlier HOMD data points
# message("Filtering bins with outlier negative log ratios...")
# homdLenThres <- 10000
# homdNumSNPThres <- 40
# homdLogRThres.auto <- log2((2*(1-purity)) / (2*(1-purity) + ploidyT*purity)) - 1*homd.sd
# ## OUTLIERS IN CHROMOSOME X SHOULD BE AGNOISTIC OF SEX 
# #if (gender == "male"){
# #	homdLogRThres.X <- round(log2((1*(1-purity)) / (1*(1-purity) + ploidyT*purity)), digits = 1) - 0.1
# #}else{
# #	homdLogRThres.X <- homdLogRThres.auto
# #}
# ind.homd.cn <- cn[LogRatio < homdLogRThres.auto, which = TRUE]
# message("Removing ", length(ind.homd.cn), " negative log ratio outlier bins ...")
# ind.segs.remove <- segs[((End-Start < homdLenThres | Length.snp. < homdNumSNPThres) & Corrected_Call == "HOMD") | Median_logR < homdLogRThres.auto, which=T]
# message("Removing ", length(ind.segs.remove), " negative log ratio outlier segments ...")
# if (length(ind.segs.remove) > 0){
# 	segsToRemove <- segs[ind.segs.remove]
# 	hits <- findOverlaps(query = as(as.data.frame(segsToRemove), "GRanges"), subject = as(as.data.frame(cn), "GRanges"))
# 	ind.homd.segs <- union(subjectHits(hits), ind.homd.cn)
# }else{
# 	ind.homd.segs <- ind.homd.cn
# }

#message("Loading panel of normals: ", ichorNormPanel)
#panel <- readRDS(ichorNormPanel)
#panel.dt <- as.data.table(as.data.frame(panel))
#panel.colNames <- colnames(panel.dt[, 11:(ncol(panel.dt)-1)])
#panel.dt[, variance := apply(.SD, 1, var, na.rm=TRUE), .SDcols = panel.colNames]
#panel.sd <- 2 * sd(panel.dt$variance, na.rm=T)
#panel.dt[, outlier := variance < -panel.sd | variance > panel.sd]
#ind.panel <- which(panel$Median < -2*sd(x$panel,na.rm=T))
#ind.panel <- which(panel.dt$outlier == TRUE)
#hits <- findOverlaps(query = as(as.data.frame(cn), "GRanges"), subject = as(panel[ind.panel,], "GRanges"))
#ind.homd.panel <- queryHits(hits)
#ind.homd.all <- union(ind.homd.panel, ind.homd.segs)
# remove the elements 
# if (length(ind.homd.segs) > 0){
# 	cn <- cn[-ind.homd.segs]
# }
# #cn <- cn[ind.homd.segs, FILTER := "EXCLUDE"]
# if (length(ind.segs.remove) > 0){
# 	segs <- segs[-ind.segs.remove]
# }
#segs <- segs[ind.segs.remove, FILTER := "EXCLUDE"]

## correct copy number beyond maximum CN state based on purity and logR
correctCN <- correctIntegerCN(cn, segs, purity, ploidyT, maxCNtoCorrect.autosomes = NULL, 
		maxCNtoCorrect.X = NULL, correctHOMD = TRUE, minPurityToCorrect = 0.05, gender = gender, chrs = chrs)
segs <- correctCN$segs
cn <- correctCN$cn
## extend segments to remove gaps
centromeres <- fread(centromere)
segs <- extendSegments(segs, removeCentromeres = TRUE, centromeres = centromeres, extendToTelomeres = FALSE,
	chrs = chrs, genomeStyle = genomeStyle)

## write segments to file ##
write.table(segs, file = outSegFile, col.names=T, row.names=F, quote=F, sep="\t")
write.table(cn, file = outBinFile, col.names=T, row.names=F, quote=F, sep="\t")
## write segments without germline SNPs
outSegNoSNPFile <- gsub(".txt", ".noSNPs.txt", outSegFile)
write.table(segs[, -c("Start.snp", "End.snp")], file = outSegNoSNPFile, col.names=T, row.names=F, quote=F, sep="\t")

## write segments in IGV / GISTIC format ##
igv <- segs[, .(Sample, Chromosome, Start.snp, End.snp, Length.snp., logR_Copy_Number)]
igv[Chromosome %in% chrs[1:22], Corrected.logR := log2(logR_Copy_Number / 2)]
igv[Chromosome == chrs[grep("X", chrs)], Corrected.logR := log2(logR_Copy_Number / 1)]
igv[, logR_Copy_Number := NULL]
outIGVFile <- gsub("seg.txt", "segIGV.txt", outSegFile)
write.table(igv, file = outIGVFile, col.names=T, row.names=F, quote=F, sep="\t")

save.image(outImageFile)



