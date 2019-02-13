configfile: "config/config.yaml"
configfile: "config/samples.yaml"

include: "ichorCNA.snakefile"
include: "getAlleleCounts.snakefile"
import os.path

CLUST = {1:[1], 2:[1,2], 3:[1,2,3], 4:[1,2,3,4], 5:[1,2,3,4,5], 6:[1,2,3,4,5,6], 7:[1,2,3,4,5,6,7], 8:[1,2,3,4,5,6,7,8], 9:[1,2,3,4,5,6,7,8,9], 10:[1,2,3,4,5,6,7,8,9,10]}
PLOIDY = {2:[2], 3:[2,3], 4:[2,3,4]}


rule all:
	input: 
		expand("results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.titan.txt", tumor=config["pairings"], clustNum=CLUST[config["TitanCNA_maxNumClonalClusters"]], ploidy=PLOIDY[config["TitanCNA_maxPloidy"]]),
		#expand("results/titan/hmm/titanCNA_ploidy{ploidy}/", ploidy=PLOIDY[config["TitanCNA_maxPloidy"]]),
		"results/titan/hmm/optimalClusterSolution.txt"
		
#rule makeOutDir:
#	output:
#		"results/titan/hmm/titanCNA_ploidy{ploidy}/"
#	params:
#	shell:
#		"mkdir -p {output}"
		
rule runTitanCNA:
	input:
		alleleCounts="results/titan/tumCounts/{tumor}.tumCounts.txt",
		corrDepth="results/ichorCNA/{tumor}/{tumor}.correctedDepth.txt"		
	output:		
		titan="results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.titan.txt",
		param="results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.params.txt",
		segTxt="results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.segs.txt",
		seg="results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.seg"
	params:
		outRoot="results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}/",
		titanRscript=config["TitanCNA_rscript"],
		libdir=config["TitanCNA_libdir"],
		numCores=config["TitanCNA_numCores"],
		normal=config["TitanCNA_normalInit"],
		chrs=config["TitanCNA_chrs"],
		genomeStyle=config["genomeStyle"],
		genomeBuild=config["genomeBuild"],
		cytobandFile=config["cytobandFile"],
		estimatePloidy=config["TitanCNA_estimatePloidy"],
		estimateClonality=config["TitanCNA_estimateClonality"],
		estimateNormal=config["TitanCNA_estimateNormal"],
		centromere=config["centromere"],
		alphaK=config["TitanCNA_alphaK"],
		#alphaR=config["TitanCNA_alphaR"],
		#alleleModel=config["TitanCNA_alleleModel"],
		txnExpLen=config["TitanCNA_txnExpLen"],
		plotYlim=config["TitanCNA_plotYlim"],
		sex=config["sex"]
	log:
		"logs/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.log"
	shell:
		"Rscript {params.titanRscript} --hetFile {input.alleleCounts} --cnFile {input.corrDepth} --outFile {output.titan} --outSeg {output.segTxt} --outParam {output.param} --outIGV {output.seg} --outPlotDir {params.outRoot} --libdir {params.libdir} --id {wildcards.tumor} --numClusters {wildcards.clustNum} --numCores {params.numCores} --normal_0 {params.normal} --ploidy_0 {wildcards.ploidy} --genomeStyle {params.genomeStyle} --genomeBuild {params.genomeBuild} --cytobandFile {params.cytobandFile} --chrs \"{params.chrs}\" --estimateNormal {params.estimateNormal} --estimatePloidy {params.estimatePloidy} --estimateClonality {params.estimateClonality}  --centromere {params.centromere} --alphaK {params.alphaK} --txnExpLen {params.txnExpLen} --gender {params.sex}  --plotYlim \"{params.plotYlim}\" > {log} 2> {log}"
	
#--alleleModel {params.alleleModel} --alphaR {params.alphaR}
	
				
rule selectSolution:
	input:
		#ploidyDirs=expand("results/titan/hmm/titanCNA_ploidy{ploidy}/", ploidy=PLOIDY[config["TitanCNA_maxPloidy"]]),
		resultFiles=expand("results/titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.titan.txt", tumor=config["pairings"], clustNum=CLUST[config["TitanCNA_maxNumClonalClusters"]], ploidy=PLOIDY[config["TitanCNA_maxPloidy"]])
	output:
		"results/titan/hmm/optimalClusterSolution.txt"
	params:
		solutionRscript=config["TitanCNA_selectSolutionRscript"],
		threshold=config["TitanCNA_solutionThreshold"]
	log:
		"logs/titan/selectSolution.log"
	shell:
		"""
		ploidyRun2=results/titan/hmm/titanCNA_ploidy2/
		if [ -d results/titan/hmm/titanCNA_ploidy3/ ]; then
			ploidyRun3=results/titan/hmm/titanCNA_ploidy3/
		else
			ploidyRun3=NULL
		fi
		if [ -d results/titan/hmm/titanCNA_ploidy4/ ]; then
			ploidyRun4=results/titan/hmm/titanCNA_ploidy4/
		else
			ploidyRun4=NULL
		fi
		Rscript {params.solutionRscript} --ploidyRun2 $ploidyRun2 --ploidyRun3 $ploidyRun3 --ploidyRun4 $ploidyRun4 --threshold {params.threshold} --outFile {output} > {log} 2> {log}
		"""
	
