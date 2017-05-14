library(stringr)
library(data.table)
args <- commandArgs(TRUE)

phi2Dir <- args[1]
phi3Dir <- args[2]
phi4Dir <- args[3]
threshold <- as.numeric(args[4])
outDir <- args[5]
#outFile <- args[5]
#outLink <- args[6]

homd.snp.thres <- 1500
homd.snp.prop <- 0.02

phi2Files <- list.files(phi2Dir, pattern="params.txt", full.names=T)
if (phi3Dir != "0"){
  phi3Files <- list.files(phi3Dir, pattern="params.txt", full.names=T)
}else{
  phi3Files <- NULL
}
if (phi4Dir != "0"){
  phi4Files <- list.files(phi4Dir, pattern="params.txt", full.names=T)
}else{
  phi4Files <- NULL
}
if (is.null(phi3Files) && is.null(phi4Files)){
  stop("phi3Dir or phi4Dir must be provided.")
}
samples <-  gsub(".params.txt", "", basename(phi2Files))
numSamples <- length(samples)
patients <- unique(str_match(samples, "(.+)_cluster")[,2])
numPatients <- length(patients)

########################################################
##### FUNCTION TO FORMAT PARAMS FOR OUTPUT #############
########################################################
formatParams <- function(params){
	id <- colnames(params)
	barcode <- strsplit(id, "_cluster")[[1]][1]
	cellPrev <- strsplit(params[grepl("Clonal cluster cellular prevalence", 
	  rownames(params)), 1], " ")[[1]]
	numClust <- length(cellPrev)
	cellPrev <- paste0(format(cellPrev, digits=4), collapse=",")
	norm <- as.numeric(params[grepl("Normal contamination estimate", rownames(params)), 1])
	purity <- 1 - norm
	ploidy <- as.numeric(params[grepl("Average tumour ploidy estimate", rownames(params)), 1])
	loglik <- as.numeric(params[grepl("likelihood", rownames(params)), 1])
	sdbw <- as.numeric(params[grepl("S_Dbw validity index \\(Both\\)", rownames(params)), 1])
	return(list(id=id, barcode=barcode, numClust=numClust, cellPrev=cellPrev, 
	  purity=purity, norm=norm, ploidy=ploidy, loglik=loglik, sdbw=sdbw))
}
save.image()

getParamAllClusters <- function(phiSamples, phiStr = "2"){
	phiParams <- data.frame(id=NA, barcode=NA, numClust=NA, 
	  cellPrev=NA, purity=NA, norm=NA, ploidy=NA, loglik=NA, sdbw=NA)
	for (j in 1:length(phiSamples)){
	  phi <- read.delim(phiSamples[j], header=F, row.names=1, stringsAsFactors=F, sep="\t")
	  colnames(phi) <- gsub(".params.txt", "", basename(phiSamples[j]))	
	  param.df <- data.frame(formatParams(phi), stringsAsFactors=F)  
	  phisegs <- fread(gsub(".params.txt", ".segs.txt", phiSamples[j]))
	  stateDist <- phisegs[, sum(Length.snp.), by = TITAN_call]
	  # HOMD total proportion greater than homd.snp.prop
	  if ("HOMD" %in% stateDist$TITAN_call){
      ind <- stateDist[TITAN_call == "HOMD", V1] / sum(stateDist[, V1]) > homd.snp.prop
      # HOMD seg larger than homd.snp.thres
      ind <- ind | sum(phisegs[TITAN_call == "HOMD", Length.snp.] > homd.snp.thres) > 0
      if (ind){
        param.df[, "loglik"] <- -Inf
        param.df[, "sdbw"] <- Inf
      }
    }
    phiParams <- rbind(phiParams, param.df )
  }
  phiParams <- na.omit(phiParams)
  phiParams <- cbind(Phi=phiStr, phiParams, path = gsub(".params.txt", "", phiSamples))
  return(phiParams)
}

########################################################
################### MAIN FUNCTION ######################
########################################################
#optSolution <- matrix(NA, ncol = 9, nrow = numSamples, dimnames = list(c(), c("SampleID_cluster", "SampleID_Barcode", "NumClonalClusters", "CellularPrevalence", "Purity", "NormalContam", "Ploidy", "LogLik", "Selected.Solution")))
outLink <- paste0(basename(outDir), ".sh")
fc <- file(outLink, "w+")
optSolutionAll <- NULL
for (i in 1:numPatients){
	id <- patients[i]
	phi2Samples <- grep(id, phi2Files, value=T)
	phi2Params <- getParamAllClusters(phi2Samples, "2")

  phi3Params <- NULL
  phi3Params$loglik <- NA
  if (!is.null(phi3Files)){
    phi3Samples <- grep(id, phi3Files, value=T)	
    phi3Params <- getParamAllClusters(phi3Samples, "3")
  }
	
	phi4Params <- NULL
	phi4Params$loglik <- NA
  if (!is.null(phi4Files)){
    phi4Samples <- grep(id, phi4Files, value=T)
	  phi4Params <- getParamAllClusters(phi4Samples, "4")
  }

  ## select ploidy based on which set of ploidy solutions has consistently lower loglik
  maxInd <- apply(cbind(phi2Params$loglik * (threshold + 1), phi3Params$loglik, phi4Params$loglik), 1, which.max)
  ploidySolInd <- as.numeric(names(which.max(table(maxInd))))
  optPloidy <- switch(ploidySolInd, phi2Params, phi3Params, phi4Params)
  optSolution <- optPloidy[which.min(optPloidy$sdbw), ]
  optSolutionAll <- rbind(optSolutionAll, optSolution)
  
  lnCmd <- paste0("cp -r ", optSolution$path, "* ", outDir, "/")
	## apply link command ##
	write.table(lnCmd, file = fc, col.names=F, row.names=F, quote=F, sep="\t", append=T)
}
outFile <- paste0(basename(outDir), ".txt")
write.table(optSolutionAll, file = outFile, col.names = T, row.names = F, quote = F, sep = "\t")
close(fc)


