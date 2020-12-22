# apples
Code pertinent to processing untargeted metabolomics data and integrating with SNP array data to see significant associations and visualize results. The code follows a pipeline for integrating genomic data with multiplatform metabolomics data.

## Processing Genomics Data
Dividing into separate chromosomes for parallel execution of mGWAS in the supercomputer.

## Processing Metabolomics Data
A file is designated for each metabolomics platform (LC-MS (+), LC-MS (-), NMR). It includes missing value imputation, log2 transformation, principal components analysis (PCA) with pooled quality control samples and without, and boxplots to check for outlier samples.

## Integrating Genomics and Metabolomics via Metabolite Genome-Wide Association Studies (mGWAS)
Code developed combining bash and R to execute multivariate mGWAS via batch script on a supercomputer. Additional code is included for collating and extracting data for each analysis.

## Prioritization of SNP-Feature Associations
To reduce the scope of the data to that of most interest, steps were developed to filter SNP-Feature associations to obtain a core set of features with putative metabolite quantitative trait loci (mQTL).