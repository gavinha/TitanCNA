[![Build Status](https://travis-ci.org/gavinha/TitanCNA.svg?branch=master)](https://travis-ci.org/gavinha/TitanCNA)

# *TitanCNA*

TitanCNA a R/Bioconductor package for analyzing subclonal copy number alterations (CNA) and loss of heterozygosity (LOH) in whole genome and exome sequencing of tumours.  

Ha, G., et al. (2014). [TITAN: Inference of copy number architectures in clonal cell populations from tumour whole genome sequence data. Genome Research, 24: 1881-1893.](http://genome.cshlp.org/content/24/11/1881) (PMID: 25060187)

## Contact
Gavin Ha  
Fred Hutchinson Cancer Research Center  
contact: <gavinha@gmail.com> or <gha@fredhutch.org>  
Date: May 30, 2019  
Website: [GavinHaLab.org](https://gavinhalab.org/)

## Table of Contents
* [Links](#links)
* [News](#news)
* [Installation](#installation)
* [Usage](#usage)
* [Vignette in TitanCNA R package](#vignette-in-titancna-r-package)
* [Acknowledgements](#acknowledgements)
* [License](#software-license)

## Links
**Snakemake Workflow:** https://github.com/gavinha/TitanCNA/tree/master/scripts/snakemake  
**10X Snakemake Workflow:** https://github.com/gavinha/TitanCNA_10X_snakemake  
Google Groups: https://groups.google.com/forum/#!forum/titancna  
Publication in Genome Research: http://genome.cshlp.org/content/24/11/1881

## News
(See [NEWS](NEWS) for previous version notes)

### May 30, 2019 - TitanCNA version 1.23.1
Addressed the issue of `RangedData` being deprecated by converting code to use `GRanges` from the `GenomicRanges` package. 
New or modified functions:  
	- `wigToGRanges`: to load WIG files and store in `GRanges` object.
	- `correctIntegerCN()`: performs allelic copy number (major/minor CN) adjustment.

### August 9, 2018
Improved parameter inference by handling errors and allowing EM to continue until convergence. This fixes runs that previously would fail because samples had very low tumor content. 

### July 26, 2018
Snakemake workflow for 10X Genomics whole genome sequencing data is now included in another Git repo.
https://github.com/gavinha/TitanCNA_10X_snakemake

### TitanCNA version 1.17.1 changes 
1)  New functions: 
	- `correctIntegerCN()`: recomputes high-level copy number that is capped by the maximum CN state. 
	Performs two tasks - (1) correct log ratio based on purity and ploidy, and then convert to decimal CN value; (2) Correct bins and segments in which the original predicted integer copy number was assigned the maximum CN state; bins and segments for all of chromosome X are also corrected, if provided in the input.

2) Modified functions:
	- `plotSegmentMedians()` and `plotCNlogRByChr()`: includes argument to show color-coding for corrected copy number; defaults to TRUE for this argument.
	
3) Removed functions/manual/dependencies:
	- `extractAlleleReadCounts()`	
	- `Rsamtools` dependency

### TitanCNA version 1.15.0 changes 
1) 10X Genomics analysis
  - Please see [scripts](scripts/) for instructions on running the 10X Genomics analysis.

2) New script to help **select optimal solutions**.  
	Please see [scripts/R_scripts](https://github.com/gavinha/TitanCNA/tree/master/scripts/R_scripts)

3) Added snakemake pipeline for entire TITAN workflow  
	Please see [scripts/snakemake](scripts/snakemake).

4) New function:
	- `plotSegmentMedians()`
	- `loadHaplotypeAlleleCounts()`: loads input allele counts with phasing information
	- `plotHaplotypeFraction()`: results from 10X Genomics WGS data with phasing of haplotype blocks
  
5) Modified features (no changes for user-accessible functions):
	- updateParameters: coordinate descent estimate of ploidy update uses previously estimated normal parameter from the same corodinate descent iteration ; leads to faster convergence
  

## Installation
### Install TitanCNA R package from github

From within R-3.3.2 or higher,  
```
install.packages("devtools")
library(devtools)
install_github("gavinha/TitanCNA")
```

### Install TitanCNA from Bioconductor
From within R-3.3.2 or higher,  
```
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("TitanCNA")
```

### Install other dependencies  
1. Install the HMMcopy suite
Please follow instructions on the HMMcopy GitHub <https://github.com/shahcompbio/hmmcopy_utils>.

2. Install ichorCNA
Please follow instructions on the ichorCNA GitHub Wiki <https://github.com/broadinstitute/ichorCNA>.

## Usage
R scripts are provided to run the R component of the TITAN analysis using the TitanCNA R/Bioconductor package.  
Please go to the [scripts](scripts/) directory and look at the README there for more details.

### Snakemake workflow
A [snakemake](scripts/snakemake) is also provided in this repo.  
This workflow will run the TITAN a set of tumour-normal pairs, starting from the BAM files and generating TitanCNA outputs. It will also perform model selection at the end of the workflow to choose the optimal ploidy and clonal cluster solutions. 

## Vignette in TitanCNA R package
The PDF of the vignette can be accessed from R
```
library(TitanCNA)
browseVignettes(package = "TitanCNA")
```
The path of the file can also be located using
```
pathToInstall <- system.file(package = "TitanCNA")
pathToPdf <- paste0(pathToInstall, "/int/doc/TitanCNA.pdf)
```
The example provided will reproduce Figure 1 in the manuscript. However, it will be slightly different because the example is only based on the analysis of chr2, not genome-wide.

## Acknowledgements
TitanCNA was developed by Gavin Ha while in the laboratories of Sohrab Shah (sshah@bccrc.ca) and Sam Aparicio (saparicio@bccrc.ca) at the Dept of Molecular Oncology, BC Cancer Agency, Vancouver, Canada.  
Yikan Wang and Daniel Lai have contributed code and discussions to this project.  
The KRONOS TITAN workflow was developed by Diljot Grewal (<dgrewal@bccrc.ca>) and Jafar Taghiyar (<jtaghiyar@bccrc.ca>).  
HMMcopy was co-developed by Daniel Lai and Gavin Ha.  
 
TitanCNA was inspired by existing methods including [OncoSNP](https://sites.google.com/site/oncosnp/)  and [PyClone](https://bitbucket.org/aroth85/pyclone/wiki/Home)  

## Software License
License: [GPLv3](LICENSE)

TitanCNA R code is open source and R/Bioconductor package is under [GPLv3](LICENSE).  This applies to the v1.9.0 and all subsequent versions within and obtained from Bioconductor.  

Users who are using TitanCNA earlier than v1.9.0 not for the purpose of academic research should contact gha@fredhutch.org, sshah@bccrc.ca, and prebstein@bccancer.bc.ca to inquire about previous licensing.
