## 1. Snakemake workflow in another Git repo:
https://github.com/gavinha/TitanCNA_10X_snakemake

## 2. 10X Genomics Whole Genome Sequencing Analysis (deprecated)
*The instructions in this section remains here for additional information of each analysis step. The full snakemake pipeline will contain all these steps, including additional post-processing steps.*

R script (`titanCNA_v1.15.0_TenX.R`) for running TitanCNA analysis on sequencing data with haplotype phasing such as 10X Genomics whole genome data

**Input files**  
1. GC-corrected, normalized read coverage using the HMMcopy suite  
  a. Get molecule coverage by counting unique barcodes in non-overlapping bins.  To do this, you will need the following installed:
    * [bxtools](https://github.com/walaj/bxtools#tile)
    * samtools
    To extract the molecule coverage at 10kb bins across chromosome 1 for both tumour and normal samples:  
    ```
    samtools view -h -F 0x4 -q 60 tumour.bam 1 | bxtools tile - -b chr1.10kb.bed > tumour_bxTile/chr1.10kb.bxTile.bed  
    ```
      The chromosome bed files for 10kb bins is provided in `TenX_scripts/data/10kb.bed.tar.gz`. Just `tar xvzf 10kb.bed.tar.gz` to extract the files.  
      
    b. Normalize the molecule coverage  
    This analysis uses ichorCNA to correct molecular coverage by GC and mappability biases.  
    See https://github.com/broadinstitute/ichorCNA/tree/master/scripts
    
    ```
    # from the command line
    >Rscript TenX_scripts/getMoleculeCoverage.R --help
    Usage: TenX_scripts/getMoleculeCoverage.R [options]


    Options:
            -t TUMORBXDIR, --tumorBXDir=TUMORBXDIR
                    Path to directory containing tumor bed files for each chromosome containing BX tags.

            -n NORMALBXDIR, --normalBXDir=NORMALBXDIR
                    Path to directory containing normal bed files for each chromosome containing BX tags.

            --minReadsPerBX=MINREADSPERBX
                    Minimum number of reads per barcode.
            
            --gcWig=GCWIG
                    Path to GC-content WIG file; Required
                    
            --libdir=LIBDIR
                    Script library path to include new or modified source R files (optional). Default: [NULL]
                
            (... plus other options)  

    ```
    Here is an example  
    
     ```
     Rscript getMoleculeCoverage.R --id test  \
        --tumorBXDir tumour_bxTile/ --normalBXDir normal_bxTile/ --minReadsPerBX 2 \
        --chrs "c(1:22, \"X\")" --outDir ../ \
        --centromere TenX_scripts/data/GRCh37.p13_centromere_UCSC-gapTable.txt \
        --gcWig ichorCNA/inst/extdata/gc_hg38_1000kb.wig \
        --libdir TitanCNA/R/
     ```
        
  
2. Tumour sample allelic read counts at heterozygous SNPs (identifed from the matched normal sample).  
  a. Identify the phased heterozygous SNPs (excluding indels) from the LongRanger 2.1 analysis on the matched normal sample.  
    Use this R script from the command line to process the phased variant file (`*phased_variants.vcf.gz`) and extracts heterozygous sites after filtering.  The output VCF file suffix `*_phasedHets.vcf`
      
    ```
    # from the command line
    > Rscript TenX_scripts/getPhasedHETSitesFromLLRVCF.R --help
    Usage: Rscript getPhasedHETSitesFromLLRVCF.R [options]


    Options:
            -i INVCF, --inVCF=INVCF
                    LongRanger 2.1 phased variant result VCF file; typically has filename suffix "*phased_variants.vcf.gz". [Required]

            -q MINQUALITY, --minQuality=MINQUALITY
                    Heterozygous variants with QUAL greater than or equal to this value are considered. [Default: 100]

            -d MINDEPTH, --minDepth=MINDEPTH
                    Heterozygous variants with read depth greater than or equal to this value are considered. [Default: 10]

            -v MINVAF, --minVAF=MINVAF
                    Heterozygous variants with variant allele fraction or reference allele fraction greater than this value are considered. [Default: 0.25]

            -o OUTVCF, --outVCF=OUTVCF
                    Output VCF file with suffix "_phasedHets.vcf" [Required]

            -h, --help
                    Show this help message and exit
    ```  
    b. Extract allelic read counts from the tumour sample  
      Using this python script, extract the allele read counts at the heterozygous SNP sites identified in the previous step.  This script is meant to be run per chromosome so that users can parallelize this step.  Users will need to combine the chromosome files into a single tab-delimited file for input into the R script.  
 
    ```
    # from the command line
    # run per chromosome
    # The `pysam` library will need to be installed in python.

    usage: python TenX_scripts/getTumourAlleleCountsAtHETSites.py $chr $phasedHetsVCF $bam $refFasta $baseQuality $mapQuality $vcfQuality > $tumCountsFile

    # $chr - Chromosome to analyze (e.g. 1)
    # $phasedHetsVCF - Output VCF file from getPhasedHETSitesFromLLRVCF.R
    # $bam - BAM file path of the tumour sample
    # $refFasta - File path to the reference FASTA file
    # $baseQuality - Minimum base quality to include in counts
    # $mapuality - Minimum mapping quality of reads to include in counts
    # $vcfQuality - Exclude HET sites with lower QUAL than this value
    # $tumCountsFile - Output tab-delimited file 
    ```  
      
**Running the R script**  
1. Look at the usage of the R script  

  ```
  # from the command line
  > Rscript titanCNA_v1.15.0_TenX.R --help
  Usage: Rscript titanCNA_v1.15.0_TenX.R [options]


  Options:
        --id=ID
                Sample ID

        --hetFile=HETFILE
                File containing allelic read counts at HET sites. (Required)
        
        (... other arguments similar to titanCNA.R)
        
        --alleleModel=ALLELEMODEL
                Emission density to use for allelic input data (binomial or Gaussian). [Default: binomial]
                
        --alphaR=ALPHAR
                Hyperparaemter on the Gaussian variance for allelic fraction data; used if --alleleModel="Gaussian". [Defaule: 5000]
        
        --haplotypeBinSize=HAPLOTYPEBINSIZE
                Bin size for summarizing phased haplotypes. [Default: 100000]

        --phaseSummarizeFun=PHASESUMMARIZEFUN
                Strategy to summarize specified haplotype bins: mean, sum, SNP. [Default: sum]
  ```
  The specific arguments that are different compared to `titanCNA.R` are listed above.  `HETFILE` should be the output of `getTumourAlleleCountsAtHETSites.py`. 
  
2. Example usage of R script  

  ```
  # normalized coverage file: test.cn.txt
  # allelic read count file from getTumourAlleleCountsAtHETSites.py: test.het.txt
  Rscript titanCNA.R --id test --hetFile test.het.txt --cnFile test.cn.txt \
    --numClusters 1 --numCores 1 --normal_0 0.5 --ploidy_0 2 \
    --chrs "c(1:22, \"X\")" --estimatePloidy TRUE --outDir ./ \
    --haplotypeBinSize 1e5 --phaseSummarizeFun sum --alleleModel Gaussian --alphaR 5000
  ```
3. Running TitanCNA for multiple restarts and model selection  
  Use the same approach as for Step 3 of [Standard Whole Genome/Exome Sequencing Analysis](#wgs).
