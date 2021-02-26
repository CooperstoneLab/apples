# Integrating genomics and multi-platform metabolomics enables metabolite QTL detection in breeding-relevant apple germplasm
Authors: Emma A Bilbrey, Kathryn Williamson, Emmanuel Hatzakis, Diane Doud Miller, Jonathan Fresnedo Ramiriez, Jessica L Cooperstone.
https://doi.org/10.1101/2021.02.18.431481

# Code Repository
Code pertinent to processing untargeted metabolomics data and integrating with SNP array data to see significant associations and visualize results. The code follows a pipeline for integrating genomic data with multiplatform metabolomics data. The code is separated into 3 parts. Part 1 details data processing and visualization for metabolomics and genomics data. Part 2 is a folder of 9 mGWAS batch scripts. Part 3 concerns result processing and visualization from mGWAS and PBA with a concentration on chlorogenic acid, our proof-of-concept metabolite.

Click the .md files for a browsable preview of the code that includes the figures produced.

## Part 1: Data Processing and Visualization of Genomics & Metabolomics Data
This code includes missing value imputation, log2 transformation, principal components analysis (PCA) with pooled quality control samples and without, and boxplots to check for outlier samples in the metabolomcis experiments. Also, code is included for prepping data for mGWAS, including A, G, and H matrix calculations via the AGHmatrix package, SNP PCA scree plots, and dividing SNP data into separate chromosomes for parallel execution of mGWAS in the supercomputer.

## Part 2: Integrating Genomics and Metabolomics via Metabolite Genome-Wide Association Studies (mGWAS)
Code developed combining bash and R to execute multivariate mGWAS via batch script that were processed on Owens at the Ohio Supercomputer Center (OSC).

## Part 3: mGWAS & PBA Results Processing & Visualizations
To reduce the scope of the data to that of most interest, steps were developed to filter SNP-Feature associations to obtain a core set of features with putative metabolite quantitative trait loci (mQTL). Broad and specific data visualization strategies were developed and are included here: number of mQTL per chromosome bar chart, composite mQTL chromosome map, number of metabolomic features associated with each SNP, mQTL across the NMR spectrum, manhattan plots, and mQTL genotype boxplot.
