#' getPhasedHETsitesFromLRVCF.R
#' author: Gavin Ha 
#' Dana-Farber Cancer Institute
#' Broad Institute
#' contact: <gavinha@gmail.com> or <gavinha@broadinstitute.org>
#' date:	  May 14, 2017

#' requires R-3.3+
#' @import VariantAnnotation

library(optparse)

option_list <- list(
  make_option(c("-i", "--inVCF"), type = "character", help = "LongRanger 2.1 phased variant result VCF file; typically has filename suffix \"*phased_variants.vcf.gz\". [Required]"),
  make_option(c("-q", "--minQuality"), type = "numeric", default = 100, help = "Heterozygous variants with QUAL greater than or equal to this value are considered. [Default: %default]"),
  make_option(c("-d", "--minDepth"), type = "numeric", default=10, help = "Heterozygous variants with read depth greater than or equal to this value are considered. [Default: %default]"),
  make_option(c("-v", "--minVAF"), type = "numeric", default=0.25, help = "Heterozygous variants with variant allele fraction or reference allele fraction greater than this value are considered. [Default: %default]"),
  make_option(c("-s","--snpDB"), type="character", default=NULL, help= "VCF file of list of known SNPs (e.g. HapMap, 1000 Genomes, dbSNP)"),
  make_option(c("--altCountField"), type="character", default="AD", help="Alternate allele count field name in genotype. [Default: %default]"),
  make_option(c("-o","--outVCF"), type = "character", help = "Output VCF file with suffix \"_phasedHets.vcf\" [Required]"),
  make_option(c("--libdir"), type="character", default=NULL, help="Path to TitanCNA package directory. If provided, will source scripts after loading installed TitanCNA package.")
)

parseobj <- OptionParser(option_list=option_list, usage = "usage: Rscript %prog [options]")
opt <- parse_args(parseobj)
print(opt)

library(TitanCNA)
library(VariantAnnotation)

vcfFile <- opt$inVCF
snpFile <- opt$snpDB
minQUAL <- opt$minQuality
minDepth <- opt$minDepth
minVAF <- opt$minVAF
outVCF <- opt$outVCF
libdir <- opt$libdir
altCountField <- opt$altCountField

if (!is.null(libdir)){
  source(paste0(libdir, "R/haplotype.R"))
}

#dir.create(outDir, showWarnings = FALSE)
chrs <- c(1:22, "X")
build <- "hg19"
filterFlags <- c("PASS", "10X_RESCUED_MOLECULE_HIGH_DIVERSITY")
#minQUAL <- 100
keepGenotypes <- c("1|0", "0|1", "0/1")

#vcfFiles <- list.files(LRdir, pattern = "_phased_variants.vcf.gz$", full.name = TRUE)

#for (i in 1:length(vcfFiles)){
#id <- gsub("_phased_variants.vcf.gz", "", basename(vcfFile))
hap <- getHaplotypesFromVCF(vcfFile, chrs = chrs, build = build,
			filterFlags = filterFlags, minQUAL = minQUAL, minDepth = minDepth,
			minVAF = minVAF, keepGenotypes = keepGenotypes, 
			altCountField = altCountField, snpDB = snpFile)

#outFile <- paste0(outDir, "/", id, "_phasedHets.vcf")
## remove BX genotype field to make things faster
geno(hap$vcf)$BX <- NULL

writeVcf(hap$vcf, filename = outVCF)
bgzipFile <- bgzip(outVCF, dest = paste0(outVCF, ".gz"), 
                   overwrite = TRUE)
indexTabix(bgzipFile, format = "vcf")
	
#	outFile <- paste0(outDir, "/", id, "_phasedGR.rds")
#	saveRDS(hap$geno, file = outFile)
								
#}