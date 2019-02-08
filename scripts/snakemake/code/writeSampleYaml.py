import pandas as pd
import yaml
import argparse
import os


def main():

    parser = argparse.ArgumentParser(
    description='Write a samples.yaml file in the input asked by the TitanCNA snakemake pipeline from a standard csv file containing columns of tumor sample name, tumor bam file location, normal bam name, and normal bam file location. csv should be in wide format with each row having a tumor-normal pair. Give the column numbers for the appropriate cols.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--sampleCSV', dest='sampleCSV', action='append', required=True,
        help='Path to csv file containing sample information, should have header')
    parser.add_argument('--tnc', dest='tumorNameCol', action='append', required=True,
        help='The column index with the names of the tumor samples')
    parser.add_argument('--tbc', dest='tumorBamCol', action='append', required=True,
        help='The column index with the filepaths of the tumor sample bams')
    parser.add_argument('--nnc', dest='normalNameCol', action='append', required=True,
        help='The column index with the names of the tumor sample bams')
    parser.add_argument('--nbc', dest='normalBamCol', action='append', required=True,
        help='The column index with the filepaths of the tumor sample bames')
    parser.add_argument('--r', dest='resultsFile', default = "samples.yaml", required=False,
        help='File path for output - default is current working directory samples.yaml')


    args = parser.parse_args()
    #read in the sample csv file with pandas
    data = pd.read_csv(args.sampleCSV[0])
    tumorBamCol = int(args.tumorBamCol[0])
    tumorNameCol = int(args.tumorNameCol[0])
    normalBamCol = int(args.normalBamCol[0])
    normalNameCol = int(args.normalNameCol[0])

    #tumor samples, using pd.Series to get a unique dictionary
    sample_bam = pd.Series(data.iloc[:,tumorBamCol].values,index=data.iloc[:,tumorNameCol]).to_dict()
    #update this dictionary with the normal samples dictionary
    sample_bam.update(pd.Series(data.iloc[:,normalBamCol].values,index=data.iloc[:,normalNameCol]).to_dict())
    #grab the normals and the samples and zip them into a pretty dict for writing
    tumor_pairings = pd.Series(data.iloc[:,normalNameCol].values,index=data.iloc[:,tumorNameCol]).to_dict() 
    with open(args.resultsFile, 'w') as outfile:
        yaml.dump({"samples": sample_bam,"pairings": tumor_pairings}, outfile, default_flow_style=False)


if __name__ == '__main__':
  main()
