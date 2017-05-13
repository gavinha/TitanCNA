[![Build Status](https://travis-ci.org/gavinha/TitanCNA.svg?branch=master)](https://travis-ci.org/gavinha/TitanCNA)

# *TitanCNA R/Bioconductor package*

TitanCNA is a probabilistic framework for analyzing subclonal copy number alterations (CNA) and loss of heterozygosity (LOH) in whole genome sequencing of tumours. The main software is implemented as a R Bioconductor package with preprocessing tools implemented in Python and BAMtools (C++). 

## Contact
author: Gavin Ha 
Dana-Farber Cancer Institute
Broad Institute
contact: <gavinha@gmail.com> or <gavinha@broadinstitute.org>
Date: May 13, 2017

## Table of Contents
* [Description](#description)
* [Installation](#installation)
* [Example Usage](#usage)
* [Acknowledgements](#acknowledgements)

## Links
TitanCNA GitHub: https://github.com/gavinha/TitanCNA
TitanCNA utils: https://github.com/gavinha/TitanCNA-utils
Google Groups: https://groups.google.com/forum/#!forum/titancna 
TitanCNA website: http://compbio.bccrc.ca/software/titan/
Publication in Genome Research: http://genome.cshlp.org/content/24/11/1881

## License
The license for the R Bioconductor package is now open source under GPLv3.  This applies to the v1.9.0 and all subsequent versions within and obtained from Bioconductor.  
Because the Bioconductor SVN is bridged to the GitHub TitanCNA repository "master" branch, the GitHub master branch for TitanCNA will also be under GPLv3.  
Please note that any other branch in the GitHub TitanCNA repository is still under the previous academic license. Users who are using TitanCNA for the purpose of academic research are free to use all branches and versions of the software. If this does not apply to you, please contact gavinha@broadinstitute.org, sshah@bccrc.ca, and prebstein@bccancer.bc.ca.

Thank you for your interest in the TitanCNA software.

## Latest version release notes 
(See NEWS for previous version notes)

### TitanCNA version 1.15.0 changes 
1) 10X Genomics analysis
  - Please see https://github.com/gavinha/TitanCNA-utils for instructions on running the 10X Genomics data

2) New function
  - plotSegmentMedians()
  - loadHaplotypeAlleleCounts(): loads input allele counts with phasing information
  - plotHaplotypeFraction(): results from 10X Genomics WGS data with phasing of haplotype blocks
  

3) Modified features (no changes for user-accessible functions)
  - updateParameters: coordinate descent estimate of ploidy update uses previously estimated normal parameter from the same corodinate descent iteration ; leads to faster convergence
  - fwd_back_clonalCN.c: returns 2-slice marginals 

## Example 

The example provided will reproduce Figure 1 in the manuscript. However, it will be slightly different because the example is only based the analysis of chr2, not genome-wide.

