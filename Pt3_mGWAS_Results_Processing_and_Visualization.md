Integrating genomics and multi-platform metabolomics enables metabolite
QTL detection in breeding-relevant apple germplasm
================
Emma Bilbrey
2/26/2021

**Publication Authors:** Emma A Bilbrey, Kathryn Wiliamson, Emmanuel
Hatzakis, Diane Doud Miller, Jonathan Fresnedo-Ramirez, Jessica L.
Cooperstone

To understand genotype-phytochemical associations in apple fruit, we
have developed a high-throughput integration strategy for genomic and
multi-platform metabolomics data.

**Context**: 124 apple genotypes, including members of three
pedigree-connected breeding families alongside diverse cultivars and
wild selections, were genotyped and phenotyped. Metabolite genome-wide
association studies (mGWAS) were conducted with ~10,000 single
nucleotide polymorphisms alongside phenotypic data acquired via liquid
chromatography mass spectrometry (LC-MS) and 1H nuclear magnetic
resonance (NMR) untargeted metabolomics. Putative metabolite
quantitative trait loci (mQTL) were then validated via pedigree-based
analyses (PBA).

Code in this markdown represents the last third of the work for our
publication. This third is dedicated to mGWAS and PBA results processing
and visualization. Data inputs can be found in the supplemental material
for the publication.

# Part 3: Results Processing and Visualization

# mGWAS Results Processing

## Filter for Significance

SNP-feature associations were filtered for significance with -log10(p)
threshold. For LC-MS (+) and (-) data sets, a threshold of greater than
or equal to 4 (P \< .0001) was used for the progeny results and greater
than or equal to 5 (P \< .00001) for the pedigree and diverse results. A
greater than or equal to 4 filter was used for each of three populations
of NMR analyses due to the fewer comparisons relative to the LC-MS data
sets. Values were not subjected to a multiple test correction, but
significance thresholds displayed on all Manhattan plots constructed
below correspond to p = .05 with a false discovery rate (FDR) correction
for that feature.

### LCMS(+)

#### Diverse

Resulting from the GWAS run in the OSC, we have a data matrix of
-log10(pvalues) for each SNP (rows) x Metabolite (columns) combination.
What we want to know is which are significant. The GWAS() code used
included P3D=TRUE, n.PC=10, and a pedigree based kinship matrix
corrected by genetic markers developed with AGHmatrix.

##### Read mGWAS Result

Columns are metabolomic features and rows are SNPs. The first 3 columns
are metadata concerning the SNPs: Index (the number assigned to each
SNP), Linkage\_Group (chromosome that the SNP is from), and
Genetic\_Distance (the distance in cM assigned from the iGLmap)

``` r
mGWASResultPosDiv<- read_csv("TableS13_LCMS_Pos_Div_mGWASresults.csv",
                             col_names=TRUE) # be patient, this takes a minute
dim(mGWASResultPosDiv)
```

    ## [1] 11165  4869

``` r
head_short(mGWASResultPosDiv)
```

    ## # A tibble: 5 x 5
    ##   Index Linkage_Group Genetic_Distance X235.16921223973_4.… X251.16328937515_4.…
    ##   <dbl>         <dbl>            <dbl>                <dbl>                <dbl>
    ## 1     1             1            0                   0.268                0.271 
    ## 2     3             1            0.002               0.0180               0.109 
    ## 3     4             1            0.003               0.0554               0.0231
    ## 4     5             1            0.004               1.13                 0.892 
    ## 5     6             1            0.005               0.143                0.335

Get rid of all rows with 0s (SNPs that did not meet the MAF minimum),
this takes a second. We are saying take only the rows from the data
frame that have no values equal to
0.

``` r
mGWASResultPosDiv0 <- mGWASResultPosDiv[apply(mGWASResultPosDiv[4:ncol(mGWASResultPosDiv)],
                                              1,
                                              function(z) any(z!=0)),]
dim(mGWASResultPosDiv0) # we see that 10294 SNPs were actually above the MAF minimum
```

    ## [1] 10294  4869

``` r
head_short(mGWASResultPosDiv0)
```

    ## # A tibble: 5 x 5
    ##   Index Linkage_Group Genetic_Distance X235.16921223973_4.… X251.16328937515_4.…
    ##   <dbl>         <dbl>            <dbl>                <dbl>                <dbl>
    ## 1     1             1            0                   0.268                0.271 
    ## 2     3             1            0.002               0.0180               0.109 
    ## 3     4             1            0.003               0.0554               0.0231
    ## 4     5             1            0.004               1.13                 0.892 
    ## 5     6             1            0.005               0.143                0.335

Some were not above the MAF minimum because we must have been filtering
for MAF minimum of .05 with some other individuals involved.

##### Filter -log(pvalue) \>= 5

We want to see which metabolomic features have strong significant
associations with at least one SNP. In this case, we are wanting an
association to have a -log10(pvalue) of greater than or equal to 5 for a
metabolomic feature to remain in the list. Column 4 is the first column
of a metabolomic
feature.

``` r
mGWASResultPosDivFilt <- Filter(function(x) max(x) >= 5, # the max of the column must be >=5
                                mGWASResultPosDiv0[,4:ncol(mGWASResultPosDiv0)])
dim(mGWASResultPosDivFilt) # we keep the same # of SNPs but reduce our # of metabolomic features
```

    ## [1] 10294  1186

``` r
head_short(mGWASResultPosDivFilt)
```

    ## # A tibble: 5 x 5
    ##   X130.1591723148… X227.1278530284… X297.1961726115… X284.3315326656…
    ##              <dbl>            <dbl>            <dbl>            <dbl>
    ## 1            0.225          0.422             0.0402           0.658 
    ## 2            0.116          0.00156           0.0983           0.310 
    ## 3            0.442          0.559             0.0752           0.465 
    ## 4            0.912          0.0611            1.86             0.0635
    ## 5            0.221          0.880             0.188            0.285 
    ## # … with 1 more variable: X138.055231846115_1.82607766106439 <dbl>

Here, the number of columns is the number of metabolomic features with
putative mQTL from the LCMS Pos, Diverse mGWAS –\> *1,186*.

#### Pedigree

Resulting from the GWAS run in the OSC, we have a data matrix of
-log10(pvalues) for each SNP (rows) x Metabolomic feature (columns)
combination. What we want to know is which are significant. The GWAS()
code used included P3D=TRUE, n.PC=6, and a pedigree based kinship matrix
corrected by genetic markers developed with AGHmatrix.

##### Read mGWAS Result

Columns are metabolomic features and rows are SNPs. The first 3 columns
are metadata concerning the SNPs: Index (the number assigned to each
SNP), Linkage\_Group (chromosome that the SNP is from), and
Genetic\_Distance (the distance in cM assigned from the iGLmap)

``` r
mGWASResultPosPed <- read_csv("TableS14_LCMS_Pos_Ped_mGWASresults.csv",
                              col_names=TRUE) # this takes a minute
dim(mGWASResultPosPed)
```

    ## [1] 11165  4858

``` r
head_short(mGWASResultPosPed)
```

    ## # A tibble: 5 x 5
    ##   Index Linkage_Group Genetic_Distance X235.16921223973_4.… X251.16328937515_4.…
    ##   <dbl>         <dbl>            <dbl>                <dbl>                <dbl>
    ## 1     1             1            0                   0.219                0.194 
    ## 2     3             1            0.002               0.127                0.0393
    ## 3     4             1            0.003               0.0731               0.0197
    ## 4     5             1            0.004               0.442                0.374 
    ## 5     6             1            0.005               0.0465               0.169

Get rid of all rows with 0s (SNPs that did not meet the MAF minimum).
Some were not above the MAF minimum because we must have been filtering
for MAF minimum of .05 with some other individuals involved. Also, there
are more removed in the pedigree set because more individuals have been
removed

``` r
mGWASResultPosPed0 <- mGWASResultPosPed[apply(mGWASResultPosPed[4],
                                              1,
                                              function(z) any(z!=0)),]
dim(mGWASResultPosPed0)
```

    ## [1] 9859 4858

``` r
head_short(mGWASResultPosPed0)
```

    ## # A tibble: 5 x 5
    ##   Index Linkage_Group Genetic_Distance X235.16921223973_4.… X251.16328937515_4.…
    ##   <dbl>         <dbl>            <dbl>                <dbl>                <dbl>
    ## 1     1             1            0                   0.219                0.194 
    ## 2     3             1            0.002               0.127                0.0393
    ## 3     4             1            0.003               0.0731               0.0197
    ## 4     5             1            0.004               0.442                0.374 
    ## 5     6             1            0.005               0.0465               0.169

``` r
mGWASResultPosPed0 <- mGWASResultPosPed[apply(mGWASResultPosPed[4:ncol(mGWASResultPosPed)],
                                              1,
                                              function(z) any(z!=0)),]
dim(mGWASResultPosPed0)
```

    ## [1] 9859 4858

``` r
head_short(mGWASResultPosPed0)
```

    ## # A tibble: 5 x 5
    ##   Index Linkage_Group Genetic_Distance X235.16921223973_4.… X251.16328937515_4.…
    ##   <dbl>         <dbl>            <dbl>                <dbl>                <dbl>
    ## 1     1             1            0                   0.219                0.194 
    ## 2     3             1            0.002               0.127                0.0393
    ## 3     4             1            0.003               0.0731               0.0197
    ## 4     5             1            0.004               0.442                0.374 
    ## 5     6             1            0.005               0.0465               0.169

##### Filter -log(pvalue) \>= 5

We want to see which metabolomic features have strong significant
associations with at least one SNP. In this case, we are wanting an
association to have a -log10(pvalue) of greater than or equal to 5 for a
metabolomic feature to remain in the list. Column 4 is the first column
of a metabolomic
feature.

``` r
mGWASResultPosPedFilt <- Filter(function(x) max(x) >= 5, # the max of the column must be >=5
                                mGWASResultPosPed0[,4:ncol(mGWASResultPosPed0)])
dim(mGWASResultPosPedFilt) # we keep the same # of SNPs but reduce our # of metabolomic features
```

    ## [1] 9859  953

``` r
head_short(mGWASResultPosPedFilt)
```

    ## # A tibble: 5 x 5
    ##   X423.2357579386… X130.1591723148… X356.2841832005… X379.2348339363…
    ##              <dbl>            <dbl>            <dbl>            <dbl>
    ## 1           0.274            0.0497           0.0633           0.267 
    ## 2           0.228            0.0586           0.357            0.870 
    ## 3           0.549            0.141            0.426            0.0510
    ## 4           0.617            0.438            0.379            1.05  
    ## 5           0.0274           0.229            0.301            0.495 
    ## # … with 1 more variable: X227.127853028416_2.90515959383752 <dbl>

Here, the number of columns is the number of metabolomic features with
putative mQTL from the LCMS Pos, Pedigree mGWAS –\> *953*

#### Progeny

Resulting from the GWAS run in the OSC, we have a data matrix of
-log10(pvalues) for each SNP (rows) x Metabolomic feature (columns)
combination. What we want to know is which are significant. The GWAS()
code used included P3D=TRUE, n.PC=3, and a pedigree based kinship matrix
corrected by genetic markers developed with AGHmatrix.

##### Read mGWAS Result

Columns are metabolomic features and rows are SNPs. The first 3 columns
are metadata concerning the SNPs: Index (the number assigned to each
SNP), Linkage\_Group (chromosome that the SNP is from), and
Genetic\_Distance (the distance in cM assigned from the
iGLmap)

``` r
mGWASResultPosProg <- read_csv("TableS15_LCMS_Pos_Prog_mGWASresults.csv",
                               col_names=TRUE) # this takes a minute
dim(mGWASResultPosProg)
```

    ## [1] 11165  4847

``` r
head_short(mGWASResultPosProg)
```

    ## # A tibble: 5 x 5
    ##   Index Linkage_Group Genetic_Distance X235.16921223973_4.… X251.16328937515_4.…
    ##   <dbl>         <dbl>            <dbl>                <dbl>                <dbl>
    ## 1     1             1            0                   0.493                0.290 
    ## 2     3             1            0.002               0.0696               0.106 
    ## 3     4             1            0.003               0.270                0.0944
    ## 4     5             1            0.004               0                    0     
    ## 5     6             1            0.005               0.347                0.485

Get rid of all rows with 0s (SNPs that did not meet the MAF
minimum)

``` r
mGWASResultPosProg0 <- mGWASResultPosProg[apply(mGWASResultPosProg[4:ncol(mGWASResultPosProg)],
                                                1,
                                                function(z) any(z!=0)),]
dim(mGWASResultPosProg0)
```

    ## [1] 9737 4847

``` r
head_short(mGWASResultPosProg0)
```

    ## # A tibble: 5 x 5
    ##   Index Linkage_Group Genetic_Distance X235.16921223973_4.… X251.16328937515_4.…
    ##   <dbl>         <dbl>            <dbl>                <dbl>                <dbl>
    ## 1     1             1            0                   0.493                0.290 
    ## 2     3             1            0.002               0.0696               0.106 
    ## 3     4             1            0.003               0.270                0.0944
    ## 4     6             1            0.005               0.347                0.485 
    ## 5     7             1            0.006               0.652                0.994

##### Filter -log(pvalue) \>= 4

We want to see which metabolomic features have strong significant
associations with at least one SNP. In this case, we are wanting an
association to have a -log10(pvalue) of greater than or equal to 4 for a
metabolomic feature to remain in the list. Column 4 is the first column
of a metabolomic
feature.

``` r
mGWASResultPosProgFilt <- Filter(function(x) max(x) >= 4, # the max of the column must be >=4
                                 mGWASResultPosProg0[,4:ncol(mGWASResultPosProg0)])
dim(mGWASResultPosProgFilt)
```

    ## [1] 9737 1787

``` r
head_short(mGWASResultPosProgFilt)
```

    ## # A tibble: 5 x 5
    ##   X235.1692122397… X251.1632893751… X377.2656580913… X305.1740957514…
    ##              <dbl>            <dbl>            <dbl>            <dbl>
    ## 1           0.493            0.290            0.129            0.0831
    ## 2           0.0696           0.106            0.261            0.324 
    ## 3           0.270            0.0944           0.0232           0.153 
    ## 4           0.347            0.485            0.123            0.106 
    ## 5           0.652            0.994            0.0389           0.982 
    ## # … with 1 more variable: X269.137908436327_3.08924992997199 <dbl>

Here, the number of columns is the number of metabolomic features with
putative mQTL from the LCMS Pos, Progeny mGWAS –\> *1,787*

### LCMS (-)

#### Diverse

Resulting from the GWAS run in the OSC, we have a data matrix of
-log10(pvalues) for each SNP (rows) x Metabolite (columns) combination.
What we want to know is which are significant. The GWAS() code used
included P3D=TRUE, n.PC=10, and a pedigree based kinship matrix
corrected by genetic markers developed with AGHmatrix.

##### Read mGWAS Results

Columns are metabolomic features and rows are SNPs. The first 3 columns
are metadata concerning the SNPs.

``` r
mGWASResultNegDiv<- read_csv("TableS16_LCMS_Neg_Div_mGWASresults.csv",
                             col_names=TRUE) # this takes a minute
dim(mGWASResultNegDiv)
```

    ## [1] 11165  4706

``` r
head_short(mGWASResultNegDiv)
```

    ## # A tibble: 5 x 5
    ##   Index Linkage_Group Genetic_Distance X885.2037_2.98177 X525.1583_3.24969
    ##   <dbl>         <dbl>            <dbl>             <dbl>             <dbl>
    ## 1     1             1            0                0.0440             0.382
    ## 2     3             1            0.002            0.0283             0.153
    ## 3     4             1            0.003            0.215              0.267
    ## 4     5             1            0.004            0.263              0.617
    ## 5     6             1            0.005            0.0467             0.161

Get rid of all rows with 0s (SNPs that did not meet the MAF
minimum)

``` r
mGWASResultNegDiv0 <- mGWASResultNegDiv[apply(mGWASResultNegDiv[4:ncol(mGWASResultNegDiv)],
                                              1,
                                              function(z) any(z!=0)),]
dim(mGWASResultNegDiv0)
```

    ## [1] 10294  4706

``` r
head_short(mGWASResultNegDiv0)
```

    ## # A tibble: 5 x 5
    ##   Index Linkage_Group Genetic_Distance X885.2037_2.98177 X525.1583_3.24969
    ##   <dbl>         <dbl>            <dbl>             <dbl>             <dbl>
    ## 1     1             1            0                0.0440             0.382
    ## 2     3             1            0.002            0.0283             0.153
    ## 3     4             1            0.003            0.215              0.267
    ## 4     5             1            0.004            0.263              0.617
    ## 5     6             1            0.005            0.0467             0.161

##### Filter -log(pvalue) \>= 5

We want to see which metabolomic features have strong significant
associations with at least one SNP. In this case, we are wanting an
association to have a -log10(pvalue) of greater than or equal to 5 for a
metabolomic feature to remain in the list. Column 4 is the first column
of a metabolomic
feature.

``` r
mGWASResultNegDivFilt <- Filter(function(x) max(x) >= 5, # the max of the column must be >=5
                                mGWASResultNegDiv0[,4:ncol(mGWASResultNegDiv0)])
dim(mGWASResultNegDivFilt)
```

    ## [1] 10294  1370

``` r
head_short(mGWASResultNegDivFilt)
```

    ## # A tibble: 5 x 5
    ##   X739.17477_2.71… X600.12641_2.10… X349.0664_1.932… X599.11324_2.33…
    ##              <dbl>            <dbl>            <dbl>            <dbl>
    ## 1           0.0614           0.0590          0.125              0.388
    ## 2           0.0808           0.242           0.235              0.403
    ## 3           0.265            0.109           0.00692            0.754
    ## 4           0.334            0.142           0.788              0.222
    ## 5           0.343            0.0199          0.306              0.236
    ## # … with 1 more variable: X1027.23518_2.83278 <dbl>

Here, the number of columns is the number of metabolomic features with
putative mQTL from the LCMS Neg, Diverse mGWAS –\> *1,370*

#### Pedigree

Resulting from the GWAS run in the OSC, we have a data matrix of
-log10(pvalues) for each SNP (rows) x Metabolomic feature (columns)
combination. What we want to know is which are significant. The GWAS()
code used included P3D=TRUE, n.PC=6, and a pedigree based kinship matrix
corrected by genetic markers developed with AGHmatrix.

##### Read mGWAS Result

Columns are metabolomic features and rows are SNPs. The first 3 columns
are metadata concerning the SNPs.

``` r
mGWASResultNegPed<- read_csv("TableS17_LCMS_Neg_Ped_mGWASresults.csv",
                             col_names=TRUE)
dim(mGWASResultNegPed)
```

    ## [1] 11165  4704

``` r
head_short(mGWASResultNegPed)
```

    ## # A tibble: 5 x 5
    ##   Index Linkage_Group Genetic_Distance X885.2037_2.98177 X525.1583_3.24969
    ##   <dbl>         <dbl>            <dbl>             <dbl>             <dbl>
    ## 1     1             1            0                0.176             0.117 
    ## 2     3             1            0.002            0.0978            0.0166
    ## 3     4             1            0.003            0.169             0.0519
    ## 4     5             1            0.004            0.217             0.309 
    ## 5     6             1            0.005            0.0548            0.110

Get rid of all rows with 0s (SNPs that did not meet the MAF
minimum)

``` r
mGWASResultNegPed0 <- mGWASResultNegPed[apply(mGWASResultNegPed[4:ncol(mGWASResultNegPed)],
                                              1,
                                              function(z) any(z!=0)),]
dim(mGWASResultNegPed0)
```

    ## [1] 9859 4704

``` r
head_short(mGWASResultNegPed0)
```

    ## # A tibble: 5 x 5
    ##   Index Linkage_Group Genetic_Distance X885.2037_2.98177 X525.1583_3.24969
    ##   <dbl>         <dbl>            <dbl>             <dbl>             <dbl>
    ## 1     1             1            0                0.176             0.117 
    ## 2     3             1            0.002            0.0978            0.0166
    ## 3     4             1            0.003            0.169             0.0519
    ## 4     5             1            0.004            0.217             0.309 
    ## 5     6             1            0.005            0.0548            0.110

##### Filter -log(pvalue) \>= 5

We want to see which metabolomic features have strong significant
associations with at least one SNP. In this case, we are wanting an
association to have a -log10(pvalue) of greater than or equal to 5 for a
metabolomic feature to remain in the list. Column 4 is the first column
of a metabolomic
feature.

``` r
mGWASResultNegPedFilt <- Filter(function(x) max(x) >= 5, # the max of the column must be >=5
                                mGWASResultNegPed0[,4:ncol(mGWASResultNegPed0)])
dim(mGWASResultNegPedFilt)
```

    ## [1] 9859 1187

``` r
head_short(mGWASResultNegPedFilt)
```

    ## # A tibble: 5 x 5
    ##   X600.12641_2.10… X599.11324_2.33… X666.02291_2.33… X599.12186_2.10…
    ##              <dbl>            <dbl>            <dbl>            <dbl>
    ## 1           0.323            0.143            0.0125          0.266  
    ## 2           0.173            0.0328           0.317           0.0346 
    ## 3           0.128            0.414            0.0250          0.120  
    ## 4           0.0503           0.198            0.325           0.497  
    ## 5           0.0784           0.349            0.240           0.00356
    ## # … with 1 more variable: X695.20221_2.25108 <dbl>

Here, the number of columns is the number of metabolomic features with
putative mQTL from the LCMS Neg, Pedigree mGWAS –\> *1,187*

#### Progeny

Resulting from the GWAS run in the OSC, we have a data matrix of
-log10(pvalues) for each SNP (rows) x Metabolomic feature (columns)
combination. What we want to know is which are significant. The GWAS()
code used included P3D=TRUE, n.PC=3, and a pedigree based kinship matrix
corrected by genetic markers developed with AGHmatrix.

##### Read mGWAS Result

Columns are metabolomic features and rows are SNPs. The first 3 columns
are metadata concerning the SNPs.

``` r
mGWASResultNegProg<- read_csv("TableS18_LCMS_Neg_Prog_mGWASresults.csv",
                              col_names=TRUE) # this takes a minute
dim(mGWASResultNegProg)
```

    ## [1] 11165  4704

``` r
head_short(mGWASResultNegProg)
```

    ## # A tibble: 5 x 5
    ##   Index Linkage_Group Genetic_Distance X885.2037_2.98177 X525.1583_3.24969
    ##   <dbl>         <dbl>            <dbl>             <dbl>             <dbl>
    ## 1     1             1            0                 0.329            0.178 
    ## 2     3             1            0.002             0.116            0.0703
    ## 3     4             1            0.003             0.322            0.0596
    ## 4     5             1            0.004             0                0     
    ## 5     6             1            0.005             0.161            0.0766

Get rid of all rows with 0s (SNPs that did not meet the MAF
minimum)

``` r
mGWASResultNegProg0 <- mGWASResultNegProg[apply(mGWASResultNegProg[4:ncol(mGWASResultNegProg)],
                                                1,
                                                function(z) any(z!=0)),]
dim(mGWASResultNegProg0)
```

    ## [1] 9737 4704

``` r
head_short(mGWASResultNegProg0)
```

    ## # A tibble: 5 x 5
    ##   Index Linkage_Group Genetic_Distance X885.2037_2.98177 X525.1583_3.24969
    ##   <dbl>         <dbl>            <dbl>             <dbl>             <dbl>
    ## 1     1             1            0                 0.329            0.178 
    ## 2     3             1            0.002             0.116            0.0703
    ## 3     4             1            0.003             0.322            0.0596
    ## 4     6             1            0.005             0.161            0.0766
    ## 5     7             1            0.006             0.297            0.190

##### Filter -log(pvalue) \>= 4

We want to see which metabolomic features have strong significant
associations with at least one SNP. In this case, we are wanting an
association to have a -log10(pvalue) of greater than or equal to 4 for a
metabolomic feature to remain in the list. Column 4 is the first column
of a metabolomic
feature.

``` r
mGWASResultNegProgFilt <- Filter(function(x) max(x) >= 4, # the max of the column must be >=4
                                 mGWASResultNegProg0[,4:ncol(mGWASResultNegProg0)])
dim(mGWASResultNegProgFilt)
```

    ## [1] 9737 1962

``` r
head_short(mGWASResultNegProgFilt)
```

    ## # A tibble: 5 x 5
    ##   X525.1583_3.249… X600.12641_2.10… X758.54656_8.644 X599.11324_2.33…
    ##              <dbl>            <dbl>            <dbl>            <dbl>
    ## 1           0.178            0.142            0.166            0.355 
    ## 2           0.0703           0.0680           0.0615           0.0688
    ## 3           0.0596           0.113            0.140            0.522 
    ## 4           0.0766           0.0881           0.0373           0.527 
    ## 5           0.190            0.400            0.422            0.986 
    ## # … with 1 more variable: X666.02291_2.33742 <dbl>

Here, the number of columns is the number of metabolomic features with
putative mQTL from the LCMS Neg, Progeny mGWAS –\> *1,962*

### NMR

#### Diverse

Resulting from the GWAS run in the OSC, we have a data matrix of
-log10(pvalues) for each SNP (rows) x Metabolite (columns) combination.
What we want to know is which are significant. The GWAS() code used
included P3D=TRUE, n.PC=10, and a pedigree based kinship matrix
corrected by genetic markers developed with AGHmatrix.

##### Read mGWAS Result

Columns are metabolomic features (bins) and rows are SNPs. The first 3
columns are metadata concerning the SNPs.

``` r
mGWASResultNMRDiv<- read_csv("TableS19_LCMS_NMR_Div_mGWASresults.csv",
                             col_names=TRUE) # this takes a minute
dim(mGWASResultNMRDiv)
```

    ## [1] 11165   759

``` r
head_short(mGWASResultNMRDiv)
```

    ## # A tibble: 5 x 5
    ##   Index Linkage_Group Genetic_Distance X9.45.9.44 X9.44.9.43
    ##   <dbl>         <dbl>            <dbl>      <dbl>      <dbl>
    ## 1     1             1            0        0.0972      0.320 
    ## 2     3             1            0.002    0.840       0.153 
    ## 3     4             1            0.003    0.705       0.0264
    ## 4     5             1            0.004    0.00280     0.0855
    ## 5     6             1            0.005    0.455       0.139

Get rid of all rows with 0s (SNPs that did not meet the MAF
minimum)

``` r
mGWASResultNMRDiv0 <- mGWASResultNMRDiv[apply(mGWASResultNMRDiv[4:ncol(mGWASResultNMRDiv)],
                                              1,
                                              function(z) any(z!=0)),]
dim(mGWASResultNMRDiv0)
```

    ## [1] 10294   759

``` r
head_short(mGWASResultNMRDiv0)
```

    ## # A tibble: 5 x 5
    ##   Index Linkage_Group Genetic_Distance X9.45.9.44 X9.44.9.43
    ##   <dbl>         <dbl>            <dbl>      <dbl>      <dbl>
    ## 1     1             1            0        0.0972      0.320 
    ## 2     3             1            0.002    0.840       0.153 
    ## 3     4             1            0.003    0.705       0.0264
    ## 4     5             1            0.004    0.00280     0.0855
    ## 5     6             1            0.005    0.455       0.139

##### Filter -log(pvalue) \>= 4

We want to see which metabolomic features have strong significant
associations with at least one SNP. In this case, we are wanting an
association to have a -log10(pvalue) of greater than or equal to 4 for a
metabolomic feature to remain in the list. Column 4 is the first column
of a metabolomic
feature.

``` r
mGWASResultNMRDivFilt <- Filter(function(x) max(x) >= 4, # the max of the column must be >=4
                                mGWASResultNMRDiv0[,4:ncol(mGWASResultNMRDiv0)])
dim(mGWASResultNMRDivFilt)
```

    ## [1] 10294   374

``` r
head_short(mGWASResultNMRDivFilt)
```

    ## # A tibble: 5 x 5
    ##   X9.45.9.44 X9.22.9.21 X9.18.9.17 X8.91.8.9 X8.9.8.89
    ##        <dbl>      <dbl>      <dbl>     <dbl>     <dbl>
    ## 1    0.0972       0.329     0.570      0.236     0.226
    ## 2    0.840        0.588     0.0388     0.875     0.103
    ## 3    0.705        0.146     0.221      0.341     1.04 
    ## 4    0.00280      0.965     0.400      0.501     0.112
    ## 5    0.455        0.255     0.827      0.291     0.140

Here, the number of columns is the number of metabolomic features with
putative mQTL from the NMR, Diverse mGWAS –\> *374*

#### Pedigree

Resulting from the GWAS run in the OSC, we have a data matrix of
-log10(pvalues) for each SNP (rows) x Metabolomic feature (columns)
combination. What we want to know is which are significant. The GWAS()
code used included P3D=TRUE, n.PC=6, and a pedigree based kinship matrix
corrected by genetic markers developed with AGHmatrix.

##### Read mGWAS Result

Columns are metabolomic features and rows are SNPs. The first 3 columns
are metadata concerning the SNPs.

``` r
mGWASResultNMRPed<- read_csv("TableS20_LCMS_NMR_Ped_mGWASresults.csv",
                             col_names=TRUE) # this takes a minute
dim(mGWASResultNMRPed)
```

    ## [1] 11165   759

``` r
head_short(mGWASResultNMRPed)
```

    ## # A tibble: 5 x 5
    ##   Index Linkage_Group Genetic_Distance X9.45.9.44 X9.44.9.43
    ##   <dbl>         <dbl>            <dbl>      <dbl>      <dbl>
    ## 1     1             1            0         0.270      0.272 
    ## 2     3             1            0.002     0.678      0.0916
    ## 3     4             1            0.003     0.899      0.0682
    ## 4     5             1            0.004     0.454      0.323 
    ## 5     6             1            0.005     0.0398     0.0445

Get rid of all rows with 0s (SNPs that did not meet the MAF
minimum)

``` r
mGWASResultNMRPed0 <- mGWASResultNMRPed[apply(mGWASResultNMRPed[4:ncol(mGWASResultNMRPed)],
                                              1,
                                              function(z) any(z!=0)),]
dim(mGWASResultNMRPed0)
```

    ## [1] 9859  759

``` r
head_short(mGWASResultNMRPed0)
```

    ## # A tibble: 5 x 5
    ##   Index Linkage_Group Genetic_Distance X9.45.9.44 X9.44.9.43
    ##   <dbl>         <dbl>            <dbl>      <dbl>      <dbl>
    ## 1     1             1            0         0.270      0.272 
    ## 2     3             1            0.002     0.678      0.0916
    ## 3     4             1            0.003     0.899      0.0682
    ## 4     5             1            0.004     0.454      0.323 
    ## 5     6             1            0.005     0.0398     0.0445

##### Filter -log(pvalue) \>= 4

We want to see which metabolomic features have strong significant
associations with at least one SNP. In this case, we are wanting an
association to have a -log10(pvalue) of greater than or equal to 4 for a
metabolomic feature to remain in the list. Column 4 is the first column
of a metabolomic
feature.

``` r
mGWASResultNMRPedFilt <- Filter(function(x) max(x) >= 4, # the max of the column must be >=4
                                mGWASResultNMRPed0[,4:ncol(mGWASResultNMRPed0)])
dim(mGWASResultNMRPedFilt)
```

    ## [1] 9859  385

``` r
head_short(mGWASResultNMRPedFilt)
```

    ## # A tibble: 5 x 5
    ##   X9.43.9.42 X9.18.9.17 X8.91.8.9 X8.9.8.89 X8.85.8.84
    ##        <dbl>      <dbl>     <dbl>     <dbl>      <dbl>
    ## 1     0.0921     0.508     0.280     0.210       0.495
    ## 2     0.409      0.305     0.626     0.0276      0.419
    ## 3     0.0487     0.0462    0.412     0.841       0.828
    ## 4     0.174      0.290     0.0958    0.673       0.286
    ## 5     0.236      0.115     0.0524    0.554       0.288

Here, the number of columns is the number of metabolomic features with
putative mQTL from the NMR, Pedigree mGWAS –\> *385*

#### Progeny

Resulting from the GWAS run in the OSC, we have a data matrix of
-log10(pvalues) for each SNP (rows) x Metabolomic feature (columns)
combination. What we want to know is which are significant. The GWAS()
code used included P3D=TRUE, n.PC=3, and a pedigree based kinship matrix
corrected by genetic markers developed with AGHmatrix.

##### Read mGWAS Result

Columns are metabolomic features and rows are SNPs. The first 3 columns
are metadata concerning the SNPs.

``` r
mGWASResultNMRProg<- read_csv("TableS21_LCMS_NMR_Prog_mGWASresults.csv",
                              col_names=TRUE) # this takes a minute
dim(mGWASResultNMRProg)
```

    ## [1] 11165   759

``` r
head_short(mGWASResultNMRProg)
```

    ## # A tibble: 5 x 5
    ##   Index Linkage_Group Genetic_Distance X9.45.9.44 X9.44.9.43
    ##   <dbl>         <dbl>            <dbl>      <dbl>      <dbl>
    ## 1     1             1            0         0.593      0.187 
    ## 2     3             1            0.002     0.862      0.0882
    ## 3     4             1            0.003     0.848      0.0284
    ## 4     5             1            0.004     0          0     
    ## 5     6             1            0.005     0.0966     0.0183

Get rid of all rows with 0s (SNPs that did not meet the MAF
minimum)

``` r
mGWASResultNMRProg0 <- mGWASResultNMRProg[apply(mGWASResultNMRProg[4:ncol(mGWASResultNMRProg)],
                                                1,
                                                function(z) any(z!=0)),]
dim(mGWASResultNMRProg0)
```

    ## [1] 9737  759

``` r
head_short(mGWASResultNMRProg0)
```

    ## # A tibble: 5 x 5
    ##   Index Linkage_Group Genetic_Distance X9.45.9.44 X9.44.9.43
    ##   <dbl>         <dbl>            <dbl>      <dbl>      <dbl>
    ## 1     1             1            0         0.593      0.187 
    ## 2     3             1            0.002     0.862      0.0882
    ## 3     4             1            0.003     0.848      0.0284
    ## 4     6             1            0.005     0.0966     0.0183
    ## 5     7             1            0.006     0.373      0.0783

##### Filter -log(pvalue) \>= 4

We want to see which metabolomic features have strong significant
associations with at least one SNP. In this case, we are wanting an
association to have a -log10(pvalue) of greater than or equal to 4 for a
metabolomic feature to remain in the list. Column 4 is the first column
of a metabolomic
feature.

``` r
mGWASResultNMRProgFilt <- Filter(function(x) max(x) >= 4, # the max of the column must be >=4
                                 mGWASResultNMRProg0[,4:ncol(mGWASResultNMRProg0)])
dim(mGWASResultNMRProgFilt)
```

    ## [1] 9737  281

``` r
head_short(mGWASResultNMRProgFilt)
```

    ## # A tibble: 5 x 5
    ##   X9.44.9.43 X9.43.9.42 X9.39.9.38 X9.37.9.36 X8.91.8.9
    ##        <dbl>      <dbl>      <dbl>      <dbl>     <dbl>
    ## 1     0.187      0.123      0.190      0.410      0.228
    ## 2     0.0882     0.513      0.0634     0.0697     0.205
    ## 3     0.0284     0.0541     0.384      0.361      0.134
    ## 4     0.0183     0.410      0.308      0.265      0.142
    ## 5     0.0783     0.200      1.12       0.151      0.360

Here, the number of columns is the number of metabolomic features with
putative mQTL from the NMR, Progeny mGWAS –\> *281*

## Venn Diagrams - Sig. Results

From the significance-filtered mGWAS results, sets of features were
compared via Venn diagrams (Progeny ∪ Pedigree ∪ Diverse). The
intersection (Progeny ∩ Pedigree ∩ Diverse) was then extracted. This
process resulted in three core data sets, one per metabolomics platform,
of -log10(p) values for features significantly associated with at least
one SNP in each population set. The features in these core data sets
were then considered as having putative mQTL.

### LCMS(+)

First we need to get column names of the metabolites

``` r
PosDivMetabs <- colnames(mGWASResultPosDivFilt)
PosPedMetabs <- colnames(mGWASResultPosPedFilt)
PosProgMetabs <- colnames(mGWASResultPosProgFilt)
```

#### Extract Intersection

Next, get intersection data for each element of the Venn
Diagram

``` r
PosProgCount <- length(PosProgMetabs) # number of sig metabs in the pos progeny
PosPedCount <- length(PosPedMetabs) # number of sig metabs in the pos pedigree
PosDivCount <- length (PosDivMetabs) # number of sig metabs in the pos diverse
PosProgPed <- intersect(PosProgMetabs,PosPedMetabs) # overlap of the pos progeny and pedigree features
PosProgPedCount <- length(PosProgPed) # number of sig metabs in the pos progeny/pedigree overlap
PosProgDiv <- intersect(PosProgMetabs,PosDivMetabs) # overlap of the pos progeny and diverse features
PosProgDivCount <- length(PosProgDiv) # number of sig metabs in the pos progeny/diverse overlap
PosPedDiv <- intersect(PosPedMetabs,PosDivMetabs) # overlap of the pos pedigree and diverse features
PosPedDivCount <- length(PosPedDiv) # number of sig metabs in the pos pedigree/diverse overlap
PosAllintersect <- intersect(PosProgPed,PosDivMetabs) # overlap of the pos prog/ped overlap and diverse
PosAllintersectCount <- length(PosAllintersect) # number of sig metabs in the pos prog/ped/div overlap
```

Then we can use these values to plot the venn diagrams. Colors were
determined by seeing what colors went well with the coral, cyan, and
gold that I wanted to use for the 3 metabolomics platforms. Then once I
chose those, I did a color mix calculator to determine the hex code for
the colors that should go in the overlapping sections.

<https://www.sessions.edu/color-calculator/>
<https://meyerweb.com/eric/tools/color-blend/#8BEC46:CC9933:1:hex>

``` r
plot(venn(c("Progenies" = PosProgCount, # This will give names to each section of the diagram
            "Pedigree" = PosPedCount, # This will give names to each section of the diagram
            "Diverse" = PosDivCount, # This will give names to each section of the diagram
            "Progenies&Pedigree" = PosProgPedCount,
            "Pedigree&Diverse" = PosPedDivCount,
            "Progenies&Diverse" = PosProgDivCount,
            "Progenies&Pedigree&Diverse" = PosAllintersectCount),
          input = "union"),
     labels = list(fontfamily = "Calibri", cex = 1.4), # change font and size of the labels
     quantities = list(fontfamily = "Calibri", cex = 1.4), # change font and size of the numbers
     lwd = 3, # gives the width of the outline of the circles
     col = "white", # color of the circle outlines
     fills = c("#E67CEF",
               "#FC9B23",
               "#8BEC46",
               "#F18C89",
               "#B9B49B",
               "#C4C435",
               "cyan3")) # this will color the center
```

![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-32-1.png)<!-- -->

519 features were found in the intersection

### LCMS(-)

First we need to get column names of the metabolites

``` r
NegDivMetabs <- colnames(mGWASResultNegDivFilt)
NegPedMetabs <- colnames(mGWASResultNegPedFilt)
NegProgMetabs <- colnames(mGWASResultNegProgFilt)
```

#### Extract Intersection

Next, get intersection data for each element of the Venn
Diagram

``` r
NegProgCount <- length(NegProgMetabs) # number of sig metabs in the Neg progeny
NegPedCount <- length(NegPedMetabs) # number of sig metabs in the Neg pedigree
NegDivCount <- length (NegDivMetabs) # number of sig metabs in the Neg diverse
NegProgPed <- intersect(NegProgMetabs,NegPedMetabs) # overlap of the Neg progeny and pedigree features
NegProgPedCount <- length(NegProgPed) # number of sig metabs in the Neg progeny/pedigree overlap
NegProgDiv <- intersect(NegProgMetabs,NegDivMetabs) # overlap of the Neg progeny and diverse features
NegProgDivCount <- length(NegProgDiv) # number of sig metabs in the Neg progeny/diverse overlap
NegPedDiv <- intersect(NegPedMetabs,NegDivMetabs) # overlap of the Neg pedigree and diverse features
NegPedDivCount <- length(NegPedDiv) # number of sig metabs in the Neg pedigree/diverse overlap
NegAllintersect <- intersect(NegProgPed,NegDivMetabs) # overlap of the Neg prog/ped overlap and diverse
NegAllintersectCount <- length(NegAllintersect) # number of sig metabs in the Neg prog/ped/div overlap
```

Then we can use these values to plot the venn diagrams. Colors were
determined by seeing what colors went well with the coral, cyan, and
gold that I wanted to use for the 3 metabolomics platforms. Then once I
chose those, I did a color mix calculator to determine the hex code for
the colors that should go in the overlapping sections.

<https://www.sessions.edu/color-calculator/>
<https://meyerweb.com/eric/tools/color-blend/#8BEC46:CC9933:1:hex>

``` r
plot(venn(c("Progenies" = NegProgCount, # This will give names to each section of the diagram
            "Pedigree" = NegPedCount, # This will give names to each section of the diagram
            "Diverse" = NegDivCount, # This will give names to each section of the diagram
            "Progenies&Pedigree" = NegProgPedCount,
            "Pedigree&Diverse" = NegPedDivCount,
            "Progenies&Diverse" = NegProgDivCount,
            "Progenies&Pedigree&Diverse" = NegAllintersectCount),
          input = "union"),
     labels = list(fontfamily = "Calibri", cex = 1.4), # change font and size of the labels
     quantities = list(fontfamily = "Calibri", cex = 1.4), # change font and size of the numbers
     lwd = 3, # gives the width of the outline of the circles
     col = "white", # color of the circle outlines
     fills = c("#E67CEF",
               "#FC9B23",
               "#8BEC46",
               "#F18C89",
               "#B9B49B",
               "#C4C435",
               "tomato")) # this will color the center
```

![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-35-1.png)<!-- -->
726 features were found in the interesection This is interesting because
it is ~200 more mQTL than detected via LCMS. It is possible that this
difference is due to the fact that LCMS (-) is the preferred ionization
mode for sensitive detection of phenolic compounds.

### NMR

First we need to get column names of the metabolites

``` r
NMRDivMetabs <- colnames(mGWASResultNMRDivFilt)
NMRPedMetabs <- colnames(mGWASResultNMRPedFilt)
NMRProgMetabs <- colnames(mGWASResultNMRProgFilt)
```

#### Extract Intersection

Next, get intersection data for each element of the Venn
Diagram

``` r
NMRProgCount <- length(NMRProgMetabs) # number of sig metabs in the NMR progeny
NMRPedCount <- length(NMRPedMetabs) # number of sig metabs in the NMR pedigree
NMRDivCount <- length (NMRDivMetabs) # number of sig metabs in the NMR diverse
NMRProgPed <- intersect(NMRProgMetabs,NMRPedMetabs) # overlap of the NMR progeny and pedigree features
NMRProgPedCount <- length(NMRProgPed) # number of sig metabs in the NMR progeny/pedigree overlap
NMRProgDiv <- intersect(NMRProgMetabs,NMRDivMetabs) # overlap of the NMR progeny and diverse features
NMRProgDivCount <- length(NMRProgDiv) # number of sig metabs in the NMR progeny/diverse overlap
NMRPedDiv <- intersect(NMRPedMetabs,NMRDivMetabs) # overlap of the NMR pedigree and diverse features
NMRPedDivCount <- length(NMRPedDiv) # number of sig metabs in the NMR pedigree/diverse overlap
NMRAllintersect <- intersect(NMRProgPed,NMRDivMetabs) # overlap of the NMR prog/ped overlap and diverse
NMRAllintersectCount <- length(NMRAllintersect) # number of sig metabs in the NMR prog/ped/div overlap
```

Then we can use these values to plot the venn diagrams. Colors were
determined by seeing what colors went well with the coral, cyan, and
gold that I wanted to use for the 3 metabolomics platforms. Then once I
chose those, I did a color mix calculator to determine the hex code for
the colors that should go in the overlapping sections.

<https://www.sessions.edu/color-calculator/>
<https://meyerweb.com/eric/tools/color-blend/#8BEC46:CC9933:1:hex>

``` r
plot(venn(c("Progenies" = NMRProgCount, # This will give names to each section of the diagram
            "Pedigree" = NMRPedCount, # This will give names to each section of the diagram
            "Diverse" = NMRDivCount, # This will give names to each section of the diagram
            "Progenies&Pedigree" = NMRProgPedCount,
            "Pedigree&Diverse" = NMRPedDivCount,
            "Progenies&Diverse" = NMRProgDivCount,
            "Progenies&Pedigree&Diverse" = NMRAllintersectCount),
          input = "union"),
     labels = list(fontfamily = "Calibri", cex = 1.4), # change font and size of the labels
     quantities = list(fontfamily = "Calibri", cex = 1.4), # change font and size of the numbers
     lwd = 3, # gives the width of the outline of the circles
     col = "white", # color of the circle outlines
     fills = c("#E67CEF",
               "#FC9B23",
               "#8BEC46",
               "#F18C89",
               "#B9B49B",
               "#C4C435",
               "gold1")) # this will color the center
```

![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-38-1.png)<!-- -->
177 features were found in the
intersection

### DF for Sig Metabs Only

``` r
PosDivSigMetabsSubset <- mGWASResultPosDivFilt[,colnames(mGWASResultPosDivFilt) %in% PosAllintersect]
dim(PosDivSigMetabsSubset) # scolumns hould be number in the middle of Pos intersection
```

    ## [1] 10294   519

``` r
head_short(PosDivSigMetabsSubset)
```

    ## # A tibble: 5 x 5
    ##   X130.1591723148… X227.1278530284… X475.3033473554… X229.1434426144…
    ##              <dbl>            <dbl>            <dbl>            <dbl>
    ## 1            0.225          0.422             0.0763            0.229
    ## 2            0.116          0.00156           0.0381            0.174
    ## 3            0.442          0.559             0.182             0.134
    ## 4            0.912          0.0611            0.780             0.182
    ## 5            0.221          0.880             0.0664            0.365
    ## # … with 1 more variable: X385.179963469991_3.61913942577029 <dbl>

``` r
NegDivSigMetabsSubset <- mGWASResultNegDivFilt[,colnames(mGWASResultNegDivFilt) %in% NegAllintersect]
dim(NegDivSigMetabsSubset) # scolumns hould be number in the middle of Neg intersection
```

    ## [1] 10294   726

``` r
head_short(NegDivSigMetabsSubset)
```

    ## # A tibble: 5 x 5
    ##   X600.12641_2.10… X599.11324_2.33… X666.02291_2.33… X599.12186_2.10…
    ##              <dbl>            <dbl>            <dbl>            <dbl>
    ## 1           0.0590            0.388            0.258          0.0135 
    ## 2           0.242             0.403            0.895          0.385  
    ## 3           0.109             0.754            0.372          0.123  
    ## 4           0.142             0.222            0.472          0.477  
    ## 5           0.0199            0.236            0.338          0.00492
    ## # … with 1 more variable: X616.11819_1.97414 <dbl>

``` r
NMRDivSigMetabsSubset <- mGWASResultNMRDivFilt[,colnames(mGWASResultNMRDivFilt) %in% NMRAllintersect]
dim(NMRDivSigMetabsSubset) # scolumns hould be number in the middle of NMR intersection
```

    ## [1] 10294   177

``` r
head_short(NMRDivSigMetabsSubset)
```

    ## # A tibble: 5 x 5
    ##   X8.91.8.9 X8.9.8.89 X8.67.8.66 X8.54.8.53 X8.42.8.41
    ##       <dbl>     <dbl>      <dbl>      <dbl>      <dbl>
    ## 1     0.236     0.226    0.463       1.07       0.0202
    ## 2     0.875     0.103    0.602       1.83       0.118 
    ## 3     0.341     1.04     1.00        2.49       0.489 
    ## 4     0.501     0.112    0.0455      0.0943     0.237 
    ## 5     0.291     0.140    0.00436     0.298      0.141

# mGWAS Results Visualizations

Despite filtering based on -log10(pvalue) thresholds and filtering by
only extracting the features from each metabolomics platform that were
found to be significantly associated with at least 1 SNP in each of the
population subsets, the large number of putative mQTL required further
visualization strategies to narrow focus to determine metabolomic
features of interest.

## Reattach SNP Metadata

When we were filtering for the significant metabolomic features, we lost
the SNP metadata (first 3 columns of the original data we read in). Now
we want it back. We will be using the diverse mGWAS results from now on
so we will use the Div data frames but the same could be done with the
pedigree and progeny data
frames.

``` r
PosDivSigMetabsSubsetSNP <- PosDivSigMetabsSubset %>% # the df we want to add columns to
  add_column(Index = mGWASResultPosDiv0$Index, # we add a column called Index
             Linkage_Group = mGWASResultPosDiv0$Linkage_Group, # add a column called Linkage_Group
             Genetic_Distance = mGWASResultPosDiv0$Genetic_Distance, # add column Genetic_Distance
             .before = 1) # put these three columns before the first column in the data frame
dim(PosDivSigMetabsSubsetSNP)
```

    ## [1] 10294   522

``` r
head_short(PosDivSigMetabsSubsetSNP)
```

    ## # A tibble: 5 x 5
    ##   Index Linkage_Group Genetic_Distance X130.159172314837_3… X227.127853028416_2…
    ##   <dbl>         <dbl>            <dbl>                <dbl>                <dbl>
    ## 1     1             1            0                    0.225              0.422  
    ## 2     3             1            0.002                0.116              0.00156
    ## 3     4             1            0.003                0.442              0.559  
    ## 4     5             1            0.004                0.912              0.0611 
    ## 5     6             1            0.005                0.221              0.880

``` r
NegDivSigMetabsSubsetSNP <- NegDivSigMetabsSubset %>% # the df we want to add columns to
  add_column(Index = mGWASResultNegDiv0$Index, # we add a column called Index
             Linkage_Group = mGWASResultNegDiv0$Linkage_Group, # add a column called Linkage_Group
             Genetic_Distance = mGWASResultNegDiv0$Genetic_Distance, # add column Genetic_Distance
             .before = 1) # put these three columns before the first column in the data frame
dim(NegDivSigMetabsSubsetSNP)
```

    ## [1] 10294   729

``` r
head_short(NegDivSigMetabsSubsetSNP)
```

    ## # A tibble: 5 x 5
    ##   Index Linkage_Group Genetic_Distance X600.12641_2.10625 X599.11324_2.3381
    ##   <dbl>         <dbl>            <dbl>              <dbl>             <dbl>
    ## 1     1             1            0                 0.0590             0.388
    ## 2     3             1            0.002             0.242              0.403
    ## 3     4             1            0.003             0.109              0.754
    ## 4     5             1            0.004             0.142              0.222
    ## 5     6             1            0.005             0.0199             0.236

``` r
NMRDivSigMetabsSubsetSNP <- NMRDivSigMetabsSubset %>% # the df we want to add columns to
  add_column(Index = mGWASResultNMRDiv0$Index, # we add a column called Index
             Linkage_Group = mGWASResultNMRDiv0$Linkage_Group, # add a column called Linkage_Group
             Genetic_Distance = mGWASResultNMRDiv0$Genetic_Distance, # add column Genetic_Distance
             .before = 1) # put these three columns before the first column in the data frame
dim(NMRDivSigMetabsSubsetSNP)
```

    ## [1] 10294   180

``` r
head_short(NMRDivSigMetabsSubsetSNP)
```

    ## # A tibble: 5 x 5
    ##   Index Linkage_Group Genetic_Distance X8.91.8.9 X8.9.8.89
    ##   <dbl>         <dbl>            <dbl>     <dbl>     <dbl>
    ## 1     1             1            0         0.236     0.226
    ## 2     3             1            0.002     0.875     0.103
    ## 3     4             1            0.003     0.341     1.04 
    ## 4     5             1            0.004     0.501     0.112
    ## 5     6             1            0.005     0.291     0.140

## mQTL per LG Bar Chart

First, we wanted to investigate how many putative mQTL were being
detected per chromosome in each of the metabolomic platforms. This is
displayed through a simple bar chart where each bar is a chromosome and
the y axis indicates the number of metabolomic features with a putative
mQTL on each chromosome.

### Prep Data

Split up each dataset of significant features by linkage
group

``` r
PosChromoSplit <- split(PosDivSigMetabsSubsetSNP, PosDivSigMetabsSubsetSNP$Linkage_Group)
NegChromoSplit <- split(NegDivSigMetabsSubsetSNP, NegDivSigMetabsSubsetSNP$Linkage_Group)
NMRChromoSplit <- split(NMRDivSigMetabsSubsetSNP, NMRDivSigMetabsSubsetSNP$Linkage_Group)
```

Loop through all chromosome groups

``` r
QTLcountsPos <- NULL
for(i in PosChromoSplit){
  QTLcountsPos <- cbind(QTLcountsPos, ncol(Filter(function(x) max(x) >=5, i[,-c(1:3)])))
}
QTLcountsPos
```

    ##      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14]
    ## [1,]   27   47   28   38   24   14   22   29   31    25    25    29    25    20
    ##      [,15] [,16] [,17]
    ## [1,]    42   284    97

``` r
QTLcountsNeg <- NULL
for(i in NegChromoSplit){
  QTLcountsNeg <- cbind(QTLcountsNeg, ncol(Filter(function(x) max(x) >=5, i[,-c(1:3)])))
}
QTLcountsNeg
```

    ##      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14]
    ## [1,]   34   48   14   32   30   17   29   13   37    29    24    32    13     8
    ##      [,15] [,16] [,17]
    ## [1,]    34   443   131

``` r
QTLcountsNMR <- NULL
for(i in NMRChromoSplit){
  QTLcountsNMR <- cbind(QTLcountsNMR, ncol(Filter(function(x) max(x) >=4, i[,-c(1:3)])))
}
QTLcountsNMR
```

    ##      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14]
    ## [1,]   14    6    3   11    7   11    4    2    3    14     3     9     7    17
    ##      [,15] [,16] [,17]
    ## [1,]    23    95     7

### LCMS(+)

``` r
par(mgp=c(1.5,.4,0), family="Calibri")
PosBP <- barplot(QTLcountsPos,
        xlab="Chromosome",
        ylab="Number of Putative mQTL",
        ylim=c(0,325),
        names.arg=c(1:17),
        cex.axis = .8,
        col="cyan3")
title("Number of Putative mQTL per Chromosome - LC-MS (+)",
      adj=0.05,
      line=-.5,
      font.main=1)
text(x = PosBP,
     y = QTLcountsPos,
     label = QTLcountsPos,
     pos = 3,
     cex = 0.8,
     col = "black")
```

![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-50-1.png)<!-- -->

Chromosomes 16 & 17 showed up as mQTL hotspots

### LCMS(-)

``` r
par(mgp=c(1.5,.5,0), family="Calibri")
NegBP <- barplot(QTLcountsNeg,
        xlab="Chromosome",
        ylab="Number of Putative mQTL",
        ylim=c(0,500),
        cex.axis = .8,
        names.arg=c(1:17),
        col="tomato")
title("Number of Putative mQTL per Chromosome - LC-MS (-)",
      adj=0.05,
      font.main=1)
text(x = NegBP,
     y = QTLcountsNeg,
     label = QTLcountsNeg,
     pos = 3,
     cex = 0.8,
     col = "black")
```

![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-51-1.png)<!-- -->

Chromosomes 16 & 17 showed up as mQTL hotspots

### NMR

``` r
par(mgp=c(1.5,.4,0), family="Calibri")
NMRBP <- barplot(QTLcountsNMR,
        xlab="Chromosome",
        ylab="Number of Putative mQTL",
        ylim=c(0,110),
        cex.axis = .8,
        names.arg=c(1:17),
        col="gold1"
        )
title("Number of Putative mQTL per Chromosome - NMR",
      adj=0.05,
      line=-.5,
      font.main=1)
text(x = NMRBP,
     y = QTLcountsNMR,
     label = QTLcountsNMR,
     pos = 3,
     cex = 0.8,
     col = "black")
```

![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-52-1.png)<!-- -->

Chromosome 16 showed up as a clear mQTL hotspot. Chromosome 17 did not
stand out as much as it did on the LCMS platforms.

## Composite mQTL Chr Map

To investigate the distribution of these mQTL within each chromosome, a
composite map of mQTL across the genome using all three metabolomic
platforms was created.

### Filter for SNPS with Significant Assoc with at least one Feature

We want to filter the SNPs now to keep only those that are significantly
associated with at least one feature. First we get the counts for each
row of the number of values greater than the threshold.

``` r
PosSNPSwithSigAssoc <- PosDivSigMetabsSubsetSNP %>% 
  transmute(Index,
            Linkage_Group,
            Genetic_Distance,
            NumSigMetabs = rowSums(.[-c(1:3)] >= 5)) %>% 
  filter(NumSigMetabs >= 1)
dim(PosSNPSwithSigAssoc)
```

    ## [1] 467   4

``` r
head(PosSNPSwithSigAssoc)
```

    ## # A tibble: 6 x 4
    ##   Index Linkage_Group Genetic_Distance NumSigMetabs
    ##   <dbl>         <dbl>            <dbl>        <dbl>
    ## 1     5             1            0.004           14
    ## 2   101             1            1.9              1
    ## 3   228             1           15.9              1
    ## 4   294             1           23.0             11
    ## 5   339             1           27.8              1
    ## 6   364             1           28.9              1

``` r
NegSNPSwithSigAssoc <- NegDivSigMetabsSubsetSNP %>% 
  transmute(Index,
            Linkage_Group,
            Genetic_Distance,
            NumSigMetabs = rowSums(.[-c(1:3)] >= 5)) %>% 
  filter(NumSigMetabs >= 1)
dim(NegSNPSwithSigAssoc)
```

    ## [1] 487   4

``` r
head(NegSNPSwithSigAssoc)
```

    ## # A tibble: 6 x 4
    ##   Index Linkage_Group Genetic_Distance NumSigMetabs
    ##   <dbl>         <dbl>            <dbl>        <dbl>
    ## 1     5             1            0.004           17
    ## 2   101             1            1.9              2
    ## 3   126             1            3.02             4
    ## 4   170             1            7.61             1
    ## 5   216             1           13.2              1
    ## 6   226             1           15.8              1

``` r
NMRSNPSwithSigAssoc <- NMRDivSigMetabsSubsetSNP %>% 
  transmute(Index,
            Linkage_Group,
            Genetic_Distance,
            NumSigMetabs = rowSums(.[-c(1:3)] >= 4)) %>% 
  filter(NumSigMetabs >= 1)
dim(NMRSNPSwithSigAssoc)
```

    ## [1] 198   4

``` r
head(NMRSNPSwithSigAssoc)
```

    ## # A tibble: 6 x 4
    ##   Index Linkage_Group Genetic_Distance NumSigMetabs
    ##   <dbl>         <dbl>            <dbl>        <dbl>
    ## 1   101             1              1.9            2
    ## 2   250             1             18.3            1
    ## 3   339             1             27.8            1
    ## 4   478             1             44.1            2
    ## 5   477             1             44.1            1
    ## 6   480             1             44.4            1

### Object of SNP Metadata

This object should be the same for each of the metabolomics datasets
because it is just about the SNPS, so we will just take it from the pos
version.

``` r
AllSNP_LG_cM <- PosDivSigMetabsSubsetSNP %>%
  transmute(Index,
            Linkage_Group,
            Genetic_Distance)
dim(AllSNP_LG_cM) # should be the full 10294 SNPs and just 3 columns
```

    ## [1] 10294     3

``` r
head(AllSNP_LG_cM)
```

    ## # A tibble: 6 x 3
    ##   Index Linkage_Group Genetic_Distance
    ##   <dbl>         <dbl>            <dbl>
    ## 1     1             1            0    
    ## 2     3             1            0.002
    ## 3     4             1            0.003
    ## 4     5             1            0.004
    ## 5     6             1            0.005
    ## 6     7             1            0.006

Now we need to get the length in genetic distance for each of the
chromosomes so that we can plot them in our
map.

``` r
allmaxGDs <- NULL # make an empty object where we will put the output in our 'for' loop
for(i in 1:17){ # for every i in 1 through 17 (aka 1, 2, 3, 4 ..)
  maxGD <- max(AllSNP_LG_cM[AllSNP_LG_cM$Linkage_Group==i,"Genetic_Distance"]) # get the maximum genetic distance per chromosome and write it to the object maxGD
  allmaxGDs <- c(allmaxGDs,maxGD) # fill in your empty object with output
}
allmaxGDsDF <- as.data.frame(cbind(1:17,allmaxGDs)) # bind a column for linkage group to the maxes
colnames(allmaxGDsDF) <- c("chr","MaxGD") # rename the columns
allmaxGDsDF$chr <- as.factor(allmaxGDsDF$chr) # set the chromosome column as a factor.
```

Write data set names to equal color names so that they will come up in
the legend appropriately instead of just saying the color name

``` r
colors <- c("LC-MS (+)" = "cyan3",
            "LC-MS (-)" = "tomato",
            "NMR" = "gold1")
```

Create pdf of the plot with the legend. The legend has to be manual
because the color is based on each geom\_point() because the color
doesnt correspond to one variable, but instead to 3 data sets.
pdf(file=“compositeCHRsigSNPmap061420.pdf”, width=13,
height=8)

``` r
compMap<- ggplot(data = allmaxGDsDF, # this allows me not to have to have a data call for geom_bars
                 aes(x = chr)) +
  geom_bar(aes(y = MaxGD), # plot bars to represent each chromosome
           stat = "identity",
           width = .85,
           fill = "gray19",
           colour = "black") + # outline of grey bars
  geom_bar(aes(y = MaxGD), # plot a grey line in the middle to help tell chromosomes apart
           stat = "identity",
           width = .1) +
  geom_point(data = PosSNPSwithSigAssoc, # plot the Pos SNPs as blue lines
             aes(x = Linkage_Group, # this divides the SNPs onto their correct chromosome
                 y = Genetic_Distance, # this plots them at the correct genetic distance (vertical)
                 colour = "LC-MS (+)"), # this is our color because of the object we set above
             shape = 95, # this will give a flat, horizontal line
             size = 6, # this is the size of the point
             position = position_nudge(x = -.3)) + # nudge to the left side of the bar graph
  geom_point(data = NegSNPSwithSigAssoc, # plot the Neg SNPs as coral lines
             aes(x = Linkage_Group,
                 y = Genetic_Distance,
                 colour = "LC-MS (-)"),
             shape = 95,
             size = 6) + # dont need a nudge because this one will stay in the middle
  geom_point(data = NMRSNPSwithSigAssoc, # plot the NMR SNps as yellow Lines
             aes(x = Linkage_Group,
                 y = Genetic_Distance,
                 colour = "NMR"),
             shape = 95,
             size = 6,
             position = position_nudge(x = .3)) + # nudge to the right side of the bar graph
  labs(x = "Linkage Group", # xaxis label
       y = "Genetic Distance (cM)", # Y axis label
       fill = "Legend", # 
       title = "SNPs Significantly Associated with at least one Metabolomic Feature") + # title
  scale_y_continuous(expand = c(0, 0)) +
  scale_y_reverse() + # put 0 at the top
  scale_x_discrete(name ="Linkage Group",
                   limits = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17")) +
  scale_color_manual(name = "Metabolomic Approach", # title for the legend
                     values = colors, # the values we assigned above
                     labels = c("LC-MS (+)","LC-MS (-)", "NMR"),
                     breaks = c("LC-MS (+)", "LC-MS (-)", "NMR")) +
  theme(axis.ticks.x = element_blank(), # remove ticks on xaxis
        axis.text = element_text(size = 15,
                                 vjust = 1),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 20),
        legend.key = element_rect(fill = 'gray19'), # fill behind the colored lines
        legend.position = c(.1,.1), # put the legend in the lower left corner of the plot
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        panel.grid = element_blank(), # remove the grid in the plot background
        text = element_text(family = "Calibri"))
compMap
```

![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-59-1.png)<!-- -->

Horizontal lines indicate the location of a SNP found to have a
significant association with at least one metabolomic feature in the
filtered and intersected mGWAS results. Lines are colored based on the
origin of the metabolomic feature. mQTL associated with features from
the LC-MS (+) are blue, LC-MS (-) are coral, and NMR are yellow.

### Plotting for Chlorogenic Acid

I want to select SNPS to plot only for the ones that are over the
threshold for chlorogenic acid. Based on the lines drawn for FDR on the
Manhattan plot, I am going to use 4 as the threshold.

``` r
PosSNPSwithSigAssocCGA <- PosDivSigMetabsSubsetSNP %>% 
  transmute(Index, # keep Index
            Linkage_Group, # keep linkage group
            Genetic_Distance, # keep genetic distance
            CGApos = X355.103099792586_2.2386339055794) %>% 
  filter(CGApos >= 4)
dim(PosSNPSwithSigAssocCGA)
```

    ## [1] 26  4

``` r
head(PosSNPSwithSigAssocCGA)
```

    ## # A tibble: 6 x 4
    ##   Index Linkage_Group Genetic_Distance CGApos
    ##   <dbl>         <dbl>            <dbl>  <dbl>
    ## 1 14854            17             35.8   4.53
    ## 2 14874            17             39.3   5.71
    ## 3 15090            17             58.2   5.85
    ## 4 15095            17             59.2   6.10
    ## 5 15096            17             59.2   4.60
    ## 6 15099            17             59.6   5.61

``` r
NegSNPSwithSigAssocCGA <- NegDivSigMetabsSubsetSNP %>% 
  transmute(Index, # keep Index
            Linkage_Group, # keep linkage group
            Genetic_Distance, # keep genetic distance
            CGANeg = X353.09194_2.23795) %>% 
  filter(CGANeg >= 4)
dim(NegSNPSwithSigAssocCGA)
```

    ## [1] 27  4

``` r
head(NegSNPSwithSigAssocCGA)
```

    ## # A tibble: 6 x 4
    ##   Index Linkage_Group Genetic_Distance CGANeg
    ##   <dbl>         <dbl>            <dbl>  <dbl>
    ## 1 14874            17             39.3   5.33
    ## 2 15090            17             58.2   6.24
    ## 3 15095            17             59.2   6.79
    ## 4 15096            17             59.2   5.96
    ## 5 15099            17             59.6   5.86
    ## 6 15098            17             59.6   4.76

``` r
NMRSNPSwithSigAssocCGA <- NMRDivSigMetabsSubsetSNP %>% 
  transmute(Index, # keep Index
            Linkage_Group, # keep linkage group
            Genetic_Distance, # keep genetic distance
            CGANMR = X2.15.2.14) %>% 
  filter(CGANMR >= 4)
dim(NMRSNPSwithSigAssocCGA)
```

    ## [1] 20  4

``` r
head(NMRSNPSwithSigAssocCGA)
```

    ## # A tibble: 6 x 4
    ##   Index Linkage_Group Genetic_Distance CGANMR
    ##   <dbl>         <dbl>            <dbl>  <dbl>
    ## 1 14737            17             27.5   4.99
    ## 2 14748            17             30.1   4.38
    ## 3 14756            17             30.1   4.04
    ## 4 14765            17             31.5   5.37
    ## 5 14854            17             35.8   4.45
    ## 6 14874            17             39.3   4.56

``` r
compMapCGA <- ggplot(data = allmaxGDsDF, # this allows me to skip a data call for geom_bars
                     aes(x = chr)) +
  geom_bar(aes(y = MaxGD), # plot bars to represent each chromosome
           stat = "identity",
           width = .85,
           fill = "gray19",
           colour = "black") + # outline of grey bars
  geom_bar(aes(y = MaxGD), # plot a grey line in the middle to help tell chromosomes apart
           stat = "identity",
           width = .1) +
  geom_point(data = PosSNPSwithSigAssocCGA, # plot the Pos SNPs as blue lines
             aes(x = Linkage_Group, # this divides the SNPs onto their correct chromosome
                 y = Genetic_Distance, # this plots them at the correct genetic distance (vertical)
                 colour = "LC-MS (+)"), # this is our color because of the object we set above
             shape = 95, # this will give a flat, horizontal line
             size = 6, # this is the size of the point
             position = position_nudge(x = -.3)) + # nudge to the left side of the bar graph
  geom_point(data = NegSNPSwithSigAssocCGA, # plot the Neg SNPs as coral lines
             aes(x = Linkage_Group,
                 y = Genetic_Distance,
                 colour = "LC-MS (-)"),
             shape = 95,
             size = 6) + # dont need a nudge because this one will stay in the middle
  geom_point(data = NMRSNPSwithSigAssocCGA, # plot the NMR SNps as yellow Lines
             aes(x = Linkage_Group,
                 y = Genetic_Distance,
                 colour = "NMR"),
             shape = 95,
             size = 6,
             position = position_nudge(x = .3)) + # nudge to the right side of the bar graph
  labs(x = "Linkage Group", # xaxis label
       y = "Genetic Distance (cM)", # Y axis label
       fill = "Legend", # 
       title = "SNPs Significantly Associated with Chlorogenic Acid") + # title
  scale_y_continuous(expand = c(0, 0)) +
  scale_y_reverse() + # put 0 at the top
  scale_x_discrete(name ="Linkage Group",
                   limits = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17")) +
  scale_color_manual(name = "mGWAS Significant SNPs", # title for the legend
                     values = colors, # the values we assigned above
                     labels = c("LC-MS (+)","LC-MS (-)", "NMR"),
                     breaks = c("LC-MS (+)", "LC-MS (-)", "NMR")) +
  theme(axis.ticks.x = element_blank(), # remove ticks on xaxis
        axis.text = element_text(size = 15,
                                 vjust = 1),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 20),
        legend.key = element_rect(fill = 'gray19'), # fill behind the colored lines
        legend.position = c(.1,.1), # put the legend in the lower left corner of the plot
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        panel.grid = element_blank(), # remove the grid in the plot background
        text = element_text(family = "Calibri"))
compMapCGA
```

![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-64-1.png)<!-- -->

## Number of Metabolomic Features Significantly Associated with each SNP

Additionally, to understand if many disparate SNPs were eliciting signal
or if a few SNPs were associated with many metabolomic features, plots
were constructed to visualize the number of significantly associated
features per SNP within each chromosome. This uses the core significant
features, not the full list of orginal metabolomic features.

### Prep Data

Subset For Chromosome and Genetic Distance subsetting the significant
rows for their linkage group. this will allow us to plot them on a
separate plot for each chromosome.

``` r
SigMetabsPerSNPPos <- PosDivSigMetabsSubsetSNP %>% 
  transmute(Index,
            Linkage_Group,
            Genetic_Distance,
            NumSigMetabs = rowSums(.[-c(1:3)] >= 5)) %>% 
  remove_rownames %>%
  column_to_rownames(var="Index") %>%
  as.data.frame()
dim(SigMetabsPerSNPPos)
```

    ## [1] 10294     3

``` r
head(SigMetabsPerSNPPos)
```

    ##   Linkage_Group Genetic_Distance NumSigMetabs
    ## 1             1            0.000            0
    ## 3             1            0.002            0
    ## 4             1            0.003            0
    ## 5             1            0.004           14
    ## 6             1            0.005            0
    ## 7             1            0.006            0

``` r
SigMetabsPerSNPNeg <- NegDivSigMetabsSubsetSNP %>% 
  transmute(Index,
            Linkage_Group,
            Genetic_Distance,
            NumSigMetabs = rowSums(.[-c(1:3)] >= 5)) %>% 
  remove_rownames %>%
  column_to_rownames(var="Index") %>%
  as.data.frame()
dim(SigMetabsPerSNPNeg)
```

    ## [1] 10294     3

``` r
head(SigMetabsPerSNPNeg)
```

    ##   Linkage_Group Genetic_Distance NumSigMetabs
    ## 1             1            0.000            0
    ## 3             1            0.002            0
    ## 4             1            0.003            0
    ## 5             1            0.004           17
    ## 6             1            0.005            0
    ## 7             1            0.006            0

``` r
SigMetabsPerSNPNMR <- NMRDivSigMetabsSubsetSNP %>% 
  transmute(Index,
            Linkage_Group,
            Genetic_Distance,
            NumSigMetabs = rowSums(.[-c(1:3)] >= 4)) %>% 
  remove_rownames %>%
  column_to_rownames(var="Index") %>%
  as.data.frame()
dim(SigMetabsPerSNPNMR)
```

    ## [1] 10294     3

``` r
head(SigMetabsPerSNPNMR)
```

    ##   Linkage_Group Genetic_Distance NumSigMetabs
    ## 1             1            0.000            0
    ## 3             1            0.002            0
    ## 4             1            0.003            0
    ## 5             1            0.004            0
    ## 6             1            0.005            0
    ## 7             1            0.006            0

### LCMS(+)

``` r
for(i in 1:17){
  x <- SigMetabsPerSNPPos[SigMetabsPerSNPPos$Linkage_Group==i,]
  y <- x[x$NumSigMetabs>0,]
plot <- ggplot() +
  geom_point(aes(x=Genetic_Distance,
                 y=NumSigMetabs),
             colour="cyan3",
             data=x) +
  geom_text_repel(aes(x=Genetic_Distance, y=NumSigMetabs,
                label=ifelse(NumSigMetabs >= 1,row.names(y),'')),
            size=2,
            segment.size=.2,
            point.padding=0,
            force=2,
            min.segment.length = .25,
            max.overlaps = 40,
            data=y) +
  scale_x_continuous(name="Genetic Distance (cM)",
                     breaks=scales::pretty_breaks(n = 7)) +
  scale_y_continuous(name="Number of Significantly Associated Features",
                     breaks= function(w) unique(floor(pretty(seq(0, (max(w) + 1) * 1.1))))) +
  ggtitle("Number of Significantly Associated Metabolomic Features per SNP",
          subtitle=paste("LC-MS (+) - Chromosome",i)) +
  theme_bw() +
  theme(text=element_text(family="Calibri"))+
  font("title", size = 14) +
  font("subtitle", size = 12) +
  font("xlab", size = 14) +
  font("ylab", size = 14) +
  font("xy.text", size = 12)
print(plot)
}
```

![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-69-1.png)<!-- -->![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-69-2.png)<!-- -->![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-69-3.png)<!-- -->![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-69-4.png)<!-- -->![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-69-5.png)<!-- -->![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-69-6.png)<!-- -->![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-69-7.png)<!-- -->![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-69-8.png)<!-- -->![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-69-9.png)<!-- -->![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-69-10.png)<!-- -->![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-69-11.png)<!-- -->![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-69-12.png)<!-- -->![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-69-13.png)<!-- -->![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-69-14.png)<!-- -->![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-69-15.png)<!-- -->![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-69-16.png)<!-- -->![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-69-17.png)<!-- -->

### LCMS(-)

``` r
for(i in 1:17){
  x <- SigMetabsPerSNPNeg[SigMetabsPerSNPNeg$Linkage_Group==i,]
  y <- x[x$NumSigMetabs>0,]
plot <- ggplot() +
  geom_point(aes(x=Genetic_Distance,
                 y=NumSigMetabs),
             colour="tomato",
             data=x) +
  geom_text_repel(aes(x=Genetic_Distance, y=NumSigMetabs,
                label=ifelse(NumSigMetabs >= 1,row.names(y),'')),
            size=2,
            segment.size=.2,
            point.padding=0,
            force=2,
            min.segment.length = .25,
            max.overlaps = 40,
            data=y) +
  scale_x_continuous(name="Genetic Distance (cM)",
                     breaks=scales::pretty_breaks(n = 7)) +
  scale_y_continuous(name="Number of Significantly Associated Features",
                     breaks= function(w) unique(floor(pretty(seq(0, (max(w) + 1) * 1.1))))) +
  ggtitle("Number of Significantly Associated Metabolomic Features per SNP",
          subtitle=paste("LC-MS (-) - Chromosome",i)) +
  theme_bw() +
  theme(text=element_text(family="Calibri"))+
  font("title", size = 14) +
  font("subtitle", size = 12) +
  font("xlab", size = 14) +
  font("ylab", size = 14) +
  font("xy.text", size = 12)
print(plot)
}
```

![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-70-1.png)<!-- -->![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-70-2.png)<!-- -->![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-70-3.png)<!-- -->![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-70-4.png)<!-- -->![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-70-5.png)<!-- -->![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-70-6.png)<!-- -->![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-70-7.png)<!-- -->![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-70-8.png)<!-- -->![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-70-9.png)<!-- -->![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-70-10.png)<!-- -->![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-70-11.png)<!-- -->![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-70-12.png)<!-- -->![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-70-13.png)<!-- -->![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-70-14.png)<!-- -->![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-70-15.png)<!-- -->![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-70-16.png)<!-- -->![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-70-17.png)<!-- -->

### NMR

``` r
for(i in 1:17){
  x <- SigMetabsPerSNPNMR[SigMetabsPerSNPNMR$Linkage_Group==i,]
  y <- x[x$NumSigMetabs>0,]
plot <- ggplot() +
  geom_point(aes(x=Genetic_Distance,
                 y=NumSigMetabs),
             colour="gold1",
             data=x) +
  geom_text_repel(aes(x=Genetic_Distance, y=NumSigMetabs,
                label=ifelse(NumSigMetabs >= 1,row.names(y),'')),
            size=2,
            segment.size=.2,
            point.padding=0,
            force=2,
            min.segment.length = .25,
            max.overlaps = 40,
            data=y) +
  scale_x_continuous(name="Genetic Distance (cM)",
                     breaks=scales::pretty_breaks(n = 7)) +
  scale_y_continuous(name="Number of Significantly Associated Bins",
                     breaks= function(w) unique(floor(pretty(seq(0, (max(w) + 1) * 1.1))))) +
  ggtitle("Number of Significantly Associated Metabolomic Features per SNP",
          subtitle=paste("NMR - Chromosome",i)) +
  theme_bw() +
  theme(text=element_text(family="Calibri"))+
  font("title", size = 14) +
  font("subtitle", size = 12) +
  font("xlab", size = 14) +
  font("ylab", size = 14) +
  font("xy.text", size = 12)
print(plot)
}
```

![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-71-1.png)<!-- -->![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-71-2.png)<!-- -->![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-71-3.png)<!-- -->![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-71-4.png)<!-- -->![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-71-5.png)<!-- -->![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-71-6.png)<!-- -->![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-71-7.png)<!-- -->![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-71-8.png)<!-- -->![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-71-9.png)<!-- -->![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-71-10.png)<!-- -->![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-71-11.png)<!-- -->![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-71-12.png)<!-- -->![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-71-13.png)<!-- -->![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-71-14.png)<!-- -->![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-71-15.png)<!-- -->![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-71-16.png)<!-- -->![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-71-17.png)<!-- -->

## Significant NMR Bins

1D 1H NMR spectrum of the apple extract pooled quality control. Yellow
lines indicate each bin that was significantly associated with at least
one SNP. Dashed lines approximately divide the spectrum according to the
type of compounds that elicit peaks at that chemical shift. The aromatic
region and aliphatic region are in much lower abundance than the sugar
region, so magnified inserts are also
presented.

### Prepare data

``` r
sigbinsdf <- data.frame(c(8.905,8.895,8.665,8.535,8.415,8.285,7.805,7.755,7.745,7.695,7.685,7.635,7.275,7.265,7.175,7.155,7.145,7.135,7.125,7.115,7.095,7.055,7.045,7.035,7.015,6.985,6.975,6.945,6.935,6.925,6.915,6.905,6.895,6.885,6.875,6.865,6.855,6.845,6.825,6.815,6.805,6.795,6.785,6.775,6.765,6.755,6.745,6.735,6.725,6.715,6.705,6.695,6.685,6.675,6.665,6.655,6.645,6.635,6.625,6.615,6.605,6.595,6.585,6.575,6.565,6.555,6.545,6.535,6.525,6.515,6.505,6.495,6.445,6.435,6.355,6.325,6.285,6.175,6.165,6.155,6.145,6.135,6.125,6.115,6.105,6.095,6.085,6.075,6.065,6.055,6.045,6.035,6.025,6.015,6.005,5.995,5.985,5.975,5.965,5.955,5.945,5.935,5.925,5.915,5.905,5.895,5.885,5.865,5.855,5.845,5.795,5.785,5.525,5.415,5.405,5.395,5.375,5.315,5.195,5.045,5.015,4.165,4.065,4.045,3.995,3.985,3.975,3.675,3.365,3.235,2.935,2.905,2.165,2.145,2.135,2.125,2.115,1.825,1.755,1.745,1.735,1.725,1.665,1.655,1.645,1.595,1.585,1.575,1.565,1.555,1.545,1.465,1.445,1.435,1.425,1.415,1.345,1.315,1.305,1.275,1.265,1.225,1.165,1.025,0.885,0.675,0.655,0.645,0.635,0.625,0.615,0.585,0.575,0.555,0.545,0.535,0.525))
colnames(sigbinsdf) <- "Bin"
head(sigbinsdf)
```

    ##     Bin
    ## 1 8.905
    ## 2 8.895
    ## 3 8.665
    ## 4 8.535
    ## 5 8.415
    ## 6 8.285

We need to do something similar to get the bins that will draw the full
spectrum of the QC sample analyzed by
NMR

``` r
allbins <- c(9.45,9.44,9.43,9.39,9.38,9.37,9.29,9.25,9.23,9.22,9.18,9.17,9.1,9.06,9.04,9.03,9.02,9.01,9.00,8.91,8.9,8.89,8.88,8.87,8.86,8.85,8.84,8.82,8.8,8.79,8.78,8.77,8.73,8.72,8.71,8.7,8.69,8.67,8.66,8.65,8.64,8.63,8.62,8.61,8.6,8.59,8.58,8.57,8.56,8.55,8.54,8.53,8.52,8.51,8.5,8.49,8.48,8.47,8.46,8.45,8.44,8.43,8.42,8.41,8.4,8.39,8.38,8.37,8.36,8.35,8.34,8.33,8.32,8.31,8.3,8.29,8.28,8.27,8.26,8.25,8.24,8.23,8.22,8.21,8.2,8.19,8.18,8.17,8.16,8.15,8.14,8.13,8.12,8.11,8.1,8.09,8.08,8.07,8.06,8.05,8.04,8.03,8.02,8.01,8.00,7.99,7.98,7.97,7.96,7.95,7.94,7.93,7.92,7.91,7.9,7.89,7.88,7.87,7.86,7.85,7.84,7.83,7.82,7.81,7.8,7.79,7.78,7.77,7.76,7.75,7.74,7.73,7.72,7.71,7.7,7.69,7.68,7.67,7.66,7.65,7.64,7.63,7.62,7.61,7.6,7.59,7.58,7.57,7.56,7.55,7.54,7.53,7.52,7.51,7.5,7.49,7.48,7.47,7.46,7.45,7.44,7.43,7.42,7.41,7.4,7.39,7.38,7.37,7.36,7.35,7.34,7.33,7.32,7.31,7.3,7.29,7.28,7.27,7.26,7.25,7.24,7.23,7.22,7.21,7.2,7.19,7.18,7.17,7.16,7.15,7.14,7.13,7.12,7.11,7.1,7.09,7.08,7.07,7.06,7.05,7.04,7.03,7.02,7.01,7.00,6.99,6.98,6.97,6.96,6.95,6.94,6.93,6.92,6.91,6.9,6.89,6.88,6.87,6.86,6.85,6.84,6.83,6.82,6.81,6.8,6.79,6.78,6.77,6.76,6.75,6.74,6.73,6.72,6.71,6.7,6.69,6.68,6.67,6.66,6.65,6.64,6.63,6.62,6.61,6.6,6.59,6.58,6.57,6.56,6.55,6.54,6.53,6.52,6.51,6.5,6.49,6.48,6.47,6.46,6.45,6.44,6.43,6.42,6.41,6.4,6.39,6.38,6.37,6.36,6.35,6.34,6.33,6.32,6.31,6.3,6.29,6.28,6.27,6.26,6.25,6.24,6.23,6.22,6.21,6.2,6.19,6.18,6.17,6.16,6.15,6.14,6.13,6.12,6.11,6.1,6.09,6.08,6.07,6.06,6.05,6.04,6.03,6.02,6.01,6.00,5.99,5.98,5.97,5.96,5.95,5.94,5.93,5.92,5.91,5.9,5.89,5.88,5.87,5.86,5.85,5.84,5.83,5.82,5.81,5.8,5.79,5.78,5.77,5.76,5.75,5.74,5.73,5.72,5.71,5.7,5.69,5.68,5.67,5.66,5.65,5.64,5.63,5.62,5.61,5.6,5.59,5.58,5.57,5.56,5.55,5.54,5.53,5.52,5.51,5.5,5.49,5.48,5.47,5.46,5.45,5.44,5.43,5.42,5.41,5.4,5.39,5.38,5.37,5.36,5.35,5.34,5.33,5.32,5.31,5.3,5.29,5.28,5.27,5.26,5.25,5.24,5.23,5.22,5.21,5.2,5.19,5.18,5.17,5.16,5.15,5.14,5.13,5.12,5.11,5.1,5.09,5.08,5.07,5.06,5.05,5.04,5.03,5.02,4.59,4.58,4.57,4.56,4.55,4.54,4.53,4.52,4.51,4.5,4.49,4.48,4.47,4.46,4.45,4.25,4.24,4.23,4.22,4.21,4.2,4.19,4.18,4.17,4.16,4.15,4.14,4.13,4.12,4.11,4.1,4.09,4.08,4.07,4.06,4.05,4.04,4.03,4.02,4.01,4.00,3.99,3.98,3.97,3.96,3.95,3.94,3.93,3.92,3.91,3.9,3.89,3.88,3.87,3.86,3.85,3.84,3.83,3.82,3.81,3.8,3.79,3.78,3.77,3.76,3.75,3.74,3.73,3.72,3.71,3.7,3.69,3.68,3.67,3.66,3.65,3.64,3.63,3.62,3.61,3.6,3.59,3.58,3.57,3.56,3.55,3.54,3.53,3.52,3.51,3.5,3.49,3.48,3.47,3.46,3.45,3.44,3.43,3.42,3.41,3.4,3.39,3.38,3.37,3.36,3.35,3.34,3.31,3.3,3.29,3.28,3.27,3.26,3.25,3.24,3.23,3.22,3.21,3.2,3.19,3.18,3.17,3.16,3.15,3.14,3.13,3.12,3.11,3.1,3.09,3.08,3.07,3.06,3.05,3.04,3.03,3.02,3.01,3.00,2.99,2.98,2.97,2.96,2.95,2.94,2.93,2.92,2.91,2.9,2.89,2.88,2.87,2.73,2.72,2.71,2.7,2.69,2.68,2.45,2.44,2.43,2.42,2.41,2.4,2.39,2.38,2.37,2.36,2.35,2.34,2.33,2.32,2.31,2.3,2.29,2.28,2.27,2.26,2.25,2.24,2.23,2.22,2.21,2.2,2.19,2.18,2.17,2.16,2.15,2.14,2.13,2.12,2.11,2.1,2.09,2.08,2.07,2.06,2.05,2.04,2.03,2.02,2.01,2.00,1.99,1.98,1.97,1.96,1.95,1.94,1.93,1.92,1.91,1.9,1.89,1.88,1.87,1.86,1.85,1.84,1.83,1.82,1.81,1.8,1.79,1.78,1.77,1.76,1.75,1.74,1.73,1.72,1.71,1.7,1.69,1.68,1.67,1.66,1.65,1.64,1.63,1.62,1.61,1.6,1.59,1.58,1.57,1.56,1.55,1.54,1.53,1.52,1.51,1.5,1.49,1.48,1.47,1.46,1.45,1.44,1.43,1.42,1.41,1.4,1.39,1.38,1.37,1.36,1.35,1.34,1.33,1.32,1.31,1.3,1.29,1.28,1.27,1.26,1.25,1.24,1.23,1.22,1.21,1.2,1.19,1.18,1.17,1.16,1.15,1.14,1.13,1.12,1.11,1.1,1.09,1.08,1.07,1.06,1.05,1.04,1.03,1.02,1.01,1.00,0.99,0.98,0.97,0.96,0.95,0.94,0.93,0.92,0.91,0.9,0.89,0.88,0.87,0.86,0.85,0.84,0.83,0.82,0.81,0.80,0.79,0.78,0.77,0.76,0.75,0.74,0.73,0.72,0.71,0.70,0.69,0.68,0.67,0.66,0.65,0.64,0.63,0.62,0.61,0.6,0.59,0.58,0.57,0.56,0.55,0.54,0.53,0.52,0.51)
length(allbins)
```

    ## [1] 756

Need to read in the raw NMR data for the pooled qc from the table
supplement.

``` r
NMRDataQCRaw <- read_excel("TableSupplement.xlsx", # Excel file with data
                            sheet = "Table S8 1D NMR Data", # sheet within the excel file with data
                            col_names = TRUE, # we want to keep the first row as column names
                            trim_ws = TRUE, # we want to trim any white space
                            na = "0") # we want all zero values to be changed to 'NA'
dim(NMRDataQCRaw)
```

    ## [1] 125 759

``` r
head_short(NMRDataQCRaw)
```

    ## # A tibble: 5 x 5
    ##   Genotype Number Label    `9.45,9.44` `9.44,9.43`
    ##   <chr>    <chr>  <chr>          <dbl>       <dbl>
    ## 1 HC       1      Pedigree       77.3       0.935 
    ## 2 PR       2      Diverse         1.53      0.0135
    ## 3 W03      3      Diverse       366.      323.    
    ## 4 GR       4      Pedigree       86.7      64.8   
    ## 5 FJ       5      Pedigree        2.06      1.08

``` r
NMRdataforSpectrum <- NMRDataQCRaw %>% 
  remove_rownames %>%
  column_to_rownames(var="Genotype") %>% 
  dplyr::select(-c(Number, Label)) %>% 
  filter(row.names(.)=="QC")
colnames(NMRdataforSpectrum) <- allbins
head_short(NMRdataforSpectrum)
```

    ##          9.45     9.44     9.43     9.39 9.38
    ## QC   2.066059 60.86719 43.18182 246.9062  326
    ## NA         NA       NA       NA       NA   NA
    ## NA.1       NA       NA       NA       NA   NA
    ## NA.2       NA       NA       NA       NA   NA
    ## NA.3       NA       NA       NA       NA   NA

### Full Spectrum

``` r
par(mar=c(3,3,4,.5),
    mgp=c(1.5,.5,0),
    family="Calibri",
    cex.axis=1.3,
    cex.lab=1.5)
plot(colnames(NMRdataforSpectrum),NMRdataforSpectrum[1,],
     col="black",
     type="l",
     xlim=c(9,0.5),
     xlab="Chemical Shift (ppm)",
     ylab="")
abline(v=sigbinsdf$Bin, 
       lwd=.75,
       col="gold1")
abline(v=c(5.49,3), 
       lwd=2,
       lty=2,
       col="black")
title(ylab="Intensity",
      line=2)
mytitle = "Pooled QC NMR Spectrum"
mysubtitle = "Bins with Putative mQTL detected by mGWAS"
mtext(side=3, line=2, adj=0, cex=1.8, mytitle)
mtext(side=3, line=.7, adj=0, cex=1.5, mysubtitle)
```

![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-76-1.png)<!-- -->

### Aromatic Region

``` r
par(mar=c(2.75,1.5,1.4,.5),
    mgp=c(1.5,.5,0))
plot(colnames(NMRdataforSpectrum),NMRdataforSpectrum[1,],
     col="black",
     type="l",
     ylim=c(0,75000),
     xlim=c(9,5.6),
     xlab="Chemical Shift (ppm)",
     ylab="")
abline(v=sigbinsdf$Bin, 
       lwd=.75,
       col="gold1")
title(main="Aromatic Region",
      adj=0,
      line=0.25,
      font.main=1)
```

![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-78-1.png)<!-- -->

### Amino Acid Region

``` r
par(mar=c(2.75,1.5,1.4,.5),
    mgp=c(1.5,.5,0))
plot(colnames(NMRdataforSpectrum),NMRdataforSpectrum[1,],
     col="black",
     type="l",
     ylim=c(0,550000),
     xlim=c(3,.5),
     xlab="Chemical Shift (ppm)",
     ylab="")
abline(v=sigbinsdf$Bin, 
       lwd=.75,
       col="gold1")
title(main="Amino Acid Region",
      adj=0,
      line=0.25,
      font.main=1)
```

![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-80-1.png)<!-- -->

## Chlorogenic Acid Manhattan Plots

Manhattan plots of chlorogenic acid phenotypic measurements from LC-MS
(+), (-), and NMR. Alternating colors were used to help delineate
neighboring chromosomes. The dashed line indicates an FDR-corrected
q-value equivalent to p =.05. The source code for these manhattan plots
comes from the package rrblup. I have modified it to get my desired
colors and other aesthetics so I have included it here inline.

Manhattan Plot Source Code

### LCMS(+)

``` r
par(mar = c(2.5,2.5,1.5,.5),
    mgp = c(1.5,.4,0),
    family = "Calibri",
    cex.axis = 1.3,
    cex.lab = 1.5,
    cex.main = 1.8)
manhattanPos(cbind(PosDivSigMetabsSubsetSNP[,c(1,2,3)], # SNP metadata required
                   PosDivSigMetabsSubsetSNP$X355.103099792586_2.2386339055794), # chlorogenic acid
             fdr.level = 0.05)
title("Chlorogenic Acid Manhattan Plot - LC-MS (+)",
      adj = .01,
      font.main = 1)
```

![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-82-1.png)<!-- -->

### LCMS (-)

``` r
par(mar = c(2.5,2.5,1.5,.5),
    mgp = c(1.5,.4,0),
    family = "Calibri",
    cex.axis = 1.3,
    cex.lab = 1.5,
    cex.main = 1.8)
manhattanNeg(cbind(NegDivSigMetabsSubsetSNP[,c(1,2,3)], # SNP metadata required
                   NegDivSigMetabsSubsetSNP$X353.09194_2.23795), # chlorogenic acid
             fdr.level = 0.05)
title("Chlorogenic Acid Manhattan Plot - LC-MS (-)",
      adj = .01,
      font.main = 1)
```

![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-84-1.png)<!-- -->

### NMR

``` r
par(mar = c(2.5,2.5,1.5,.5),
    mgp = c(1.5,.4,0),
    family = "Calibri",
    cex.axis = 1.3,
    cex.lab = 1.5,
    cex.main = 1.8)
manhattanNMR(cbind(NMRDivSigMetabsSubsetSNP[,c(1,2,3)], # SNP metadata required
                   NMRDivSigMetabsSubsetSNP$X2.15.2.14), # chlorogenic acid
             fdr.level = 0.05)
title("Chlorogenic Acid Manhattan Plot - NMR",
      adj = .01,
      font.main = 1)
```

![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-86-1.png)<!-- -->

## Find Top SNPs on Chr

If you want to see the top list of SNPs for a certain chromosome from
the mGWAS you can adapt this code for different chromosomes and
different metabolomic features.

### Chlorogenic Acid on Chr 17

``` r
CGApos <- PosDivSigMetabsSubsetSNP %>%
  transmute(Index,
            Linkage_Group,
            Genetic_Distance,
            X355.103099792586_2.2386339055794) %>%
  filter(Linkage_Group==17) %>%
  arrange(desc(X355.103099792586_2.2386339055794))
head(CGApos)
```

    ## # A tibble: 6 x 4
    ##   Index Linkage_Group Genetic_Distance X355.103099792586_2.2386339055794
    ##   <dbl>         <dbl>            <dbl>                             <dbl>
    ## 1 15109            17             59.9                              7.66
    ## 2 15123            17             60.2                              7.54
    ## 3 15110            17             59.9                              7.52
    ## 4 15111            17             59.9                              6.58
    ## 5 15102            17             59.6                              6.52
    ## 6 15113            17             60.1                              6.42

``` r
CGAneg <- NegDivSigMetabsSubsetSNP %>%
  transmute(Index,
            Linkage_Group,
            Genetic_Distance,
            X353.09194_2.23795) %>%
  filter(Linkage_Group==17) %>%
  arrange(desc(X353.09194_2.23795))
head(CGAneg)
```

    ## # A tibble: 6 x 4
    ##   Index Linkage_Group Genetic_Distance X353.09194_2.23795
    ##   <dbl>         <dbl>            <dbl>              <dbl>
    ## 1 15109            17             59.9               7.65
    ## 2 15123            17             60.2               7.55
    ## 3 15133            17             60.6               7.44
    ## 4 15110            17             59.9               7.29
    ## 5 15113            17             60.1               7.07
    ## 6 15124            17             60.2               7.07

``` r
CGAnmr <- NMRDivSigMetabsSubsetSNP %>%
  transmute(Index,
            Linkage_Group,
            Genetic_Distance,
            X2.15.2.14) %>%
  filter(Linkage_Group==17) %>%
  arrange(desc(X2.15.2.14))
head(CGAnmr)
```

    ## # A tibble: 6 x 4
    ##   Index Linkage_Group Genetic_Distance X2.15.2.14
    ##   <dbl>         <dbl>            <dbl>      <dbl>
    ## 1 15111            17             59.9       6.07
    ## 2 15123            17             60.2       6.02
    ## 3 15096            17             59.2       5.88
    ## 4 15133            17             60.6       5.72
    ## 5 15109            17             59.9       5.71
    ## 6 15114            17             60.1       5.70

# PBA Results Visualization

## Read in CGA Haplotype Info

``` r
CGAHaplotypeAll <- read_excel("TableSupplement.xlsx", # Excel file with data
                              sheet = "Table S24 CGA Haplotypes Ped",
                              col_names = FALSE,
                              trim_ws = TRUE,
                              na = "--")
head_short(CGAHaplotypeAll)
```

    ## # A tibble: 5 x 5
    ##   ...1      ...2             ...3              ...4            ...5             
    ##   <chr>     <chr>            <chr>             <chr>           <chr>            
    ## 1 <NA>      <NA>             Additive value    <NA>            Breeding value   
    ## 2 Individu… ChlAcid_Pos      MQTRa1            MQTRd1          SUMQTR           
    ## 3 Jonathan  18.937949540000… -0.1170000000000… 0.212999999999… 9.60000000000000…
    ## 4 D02       17.4242189       -0.6530000000000… 1E-3            -0.6520000000000…
    ## 5 D03       16.568941859999… -0.6530000000000… 1E-3            -0.6520000000000…

## Prep Data

``` r
# remove row that is helpful in excel viewing but not here (1) and row that will become column names (2)
CGAHaplotypeAll1 <- CGAHaplotypeAll[-c(1:2),]
head_short(CGAHaplotypeAll1)
```

    ## # A tibble: 5 x 5
    ##   ...1     ...2             ...3              ...4             ...5             
    ##   <chr>    <chr>            <chr>             <chr>            <chr>            
    ## 1 Jonathan 18.937949540000… -0.1170000000000… 0.2129999999999… 9.60000000000000…
    ## 2 D02      17.4242189       -0.6530000000000… 1E-3             -0.6520000000000…
    ## 3 D03      16.568941859999… -0.6530000000000… 1E-3             -0.6520000000000…
    ## 4 D07      17.490109530000… -0.6530000000000… 1E-3             -0.6520000000000…
    ## 5 D08      17.36657568      -0.6530000000000… 1E-3             -0.6520000000000…

``` r
colnames(CGAHaplotypeAll1) <- CGAHaplotypeAll[2,] # make column names row 2 from first df
head_short(CGAHaplotypeAll1)
```

    ## # A tibble: 5 x 5
    ##   Individual ChlAcid_Pos      MQTRa1           MQTRd1           SUMQTR          
    ##   <chr>      <chr>            <chr>            <chr>            <chr>           
    ## 1 Jonathan   18.937949540000… -0.117000000000… 0.2129999999999… 9.6000000000000…
    ## 2 D02        17.4242189       -0.653000000000… 1E-3             -0.652000000000…
    ## 3 D03        16.568941859999… -0.653000000000… 1E-3             -0.652000000000…
    ## 4 D07        17.490109530000… -0.653000000000… 1E-3             -0.652000000000…
    ## 5 D08        17.36657568      -0.653000000000… 1E-3             -0.652000000000…

``` r
colnames(CGAHaplotypeAll1)
```

    ##  [1] "Individual"     "ChlAcid_Pos"    "MQTRa1"         "MQTRd1"        
    ##  [5] "SUMQTR"         "std_sum"        "mce_sum"        "QQ"            
    ##  [9] "Qq"             "qq"             "QTL_genotype"   "SNP_FB_0193910"
    ## [13] "SNP_FB_0193921" "SNP_FB_0102606" "SNP_FB_0102615" "SNP_FB_0102619"
    ## [17] "SNP_FB_0438354"

``` r
# subset columns that we will need for the boxplot
CGAHaplotypeAll2 <- CGAHaplotypeAll1 %>%
  select(Individual,
         QTL_genotype,
         ChlAcid_Pos)
head(CGAHaplotypeAll2) # dont need to use head_short here because there are only 2 columns
```

    ## # A tibble: 6 x 3
    ##   Individual QTL_genotype ChlAcid_Pos       
    ##   <chr>      <chr>        <chr>             
    ## 1 Jonathan   _q           18.937949540000002
    ## 2 D02        qq           17.4242189        
    ## 3 D03        qq           16.568941859999999
    ## 4 D07        qq           17.490109530000002
    ## 5 D08        qq           17.36657568       
    ## 6 D10        qq           17.921560360000001

Notice that the values for Chlorogenic acid are seen as characters. this
is because when we read it in, there were some values that were
characters in those first two lines that we have since dealt with. Now
we need to convert that column to numeric before we can get summary
statistics on
it.

``` r
CGAHaplotypeAll2$ChlAcid_Pos <- as.numeric(as.character(CGAHaplotypeAll2$ChlAcid_Pos))
```

    ## Warning: NAs introduced by coercion

``` r
head(CGAHaplotypeAll2) # now the second column should be <dbl>
```

    ## # A tibble: 6 x 3
    ##   Individual QTL_genotype ChlAcid_Pos
    ##   <chr>      <chr>              <dbl>
    ## 1 Jonathan   _q                  18.9
    ## 2 D02        qq                  17.4
    ## 3 D03        qq                  16.6
    ## 4 D07        qq                  17.5
    ## 5 D08        qq                  17.4
    ## 6 D10        qq                  17.9

Some NAs were introduced because some of the individuals in our table do
not have any phenotypic measurements (ie measurements of chlorogenic
acid in LCMS positive mode). This is because they are included to track
haplotypes in the genetic data but were not available for metabolomic
analysis. So we need to remove the rows with NAs in the second column.

``` r
CGAHaplotypeAllNaOmit <- na.omit(CGAHaplotypeAll2)
dim(CGAHaplotypeAllNaOmit) # this should give us 98 rows because that is how many individuals we have in our pedigree group - the ones that have both metabolomics and genomics data
```

    ## [1] 98  3

We also want to make QTL\_genotype a factor because that is how we want
our boxplots
organized

``` r
CGAHaplotypeAllNaOmit$QTL_genotype <- as.factor(CGAHaplotypeAllNaOmit$QTL_genotype)
head(CGAHaplotypeAllNaOmit) # now we see column one labeled as a factor
```

    ## # A tibble: 6 x 3
    ##   Individual QTL_genotype ChlAcid_Pos
    ##   <chr>      <fct>              <dbl>
    ## 1 Jonathan   _q                  18.9
    ## 2 D02        qq                  17.4
    ## 3 D03        qq                  16.6
    ## 4 D07        qq                  17.5
    ## 5 D08        qq                  17.4
    ## 6 D10        qq                  17.9

Now we want to combine levels of this factor because the \_q is likely a
qq but FlexQTL was unable to distinguish clearly between qq and Qq

``` r
levels(CGAHaplotypeAllNaOmit$QTL_genotype) <- c("qq","qq","Qq","QQ")
head(CGAHaplotypeAllNaOmit)
```

    ## # A tibble: 6 x 3
    ##   Individual QTL_genotype ChlAcid_Pos
    ##   <chr>      <fct>              <dbl>
    ## 1 Jonathan   qq                  18.9
    ## 2 D02        qq                  17.4
    ## 3 D03        qq                  16.6
    ## 4 D07        qq                  17.5
    ## 5 D08        qq                  17.4
    ## 6 D10        qq                  17.9

## Make Plot of mQTL Genotypes

We set “QTL\_genotype” as a factor, so it goes in the order we specify
in “levels”. This will be important for our plot order and assigning
colors

``` r
CGAHaplotypeAllNaOmit$QTL_genotype <- factor(CGAHaplotypeAllNaOmit$QTL_genotype,
                                             levels = c("QQ", "Qq", "qq"))
```

Then it is helpful to reorder the samples by their “Label” level as well
as setting “Genotype” as a factor with levels in the order we set for
“Label”

``` r
CGAHaplotypeAllNaOmitOrdered <- CGAHaplotypeAllNaOmit %>%
  arrange(QTL_genotype) %>% # sort by label
  mutate(Individual = fct_inorder(Individual)) # set Genotype as a factor in the order we have arranged label
```

``` r
CGAmQTLgenotypeBoxplot <- CGAHaplotypeAllNaOmit %>%
  ggplot(aes(x = QTL_genotype, 
             y = ChlAcid_Pos,
             fill = QTL_genotype)) +
  geom_boxplot(outlier.shape = NA,
               width = .35) +
  geom_point(position = position_jitterdodge()) + 
  scale_fill_manual(values=c("#006262","#00CECC","#B5F1EF")) +
  scale_y_continuous(trans = log2_trans(),
                     breaks=seq(16,22,2)) + # this makes the y-axis log2 transformed as the data is
  theme_classic() + # gives the black x and y axis lines
  theme(text=element_text(family = "Calibri"), # makes Calibri the font
        panel.background = element_blank(), # removes color from plot background
        legend.position = "none", #removes the legend
        axis.ticks.x = element_blank(), # removes the x axis ticks
        axis.text.x = element_text(face = "italic")) + # make the mQTL genotypes italic
  labs(title = "Chlorogenic Acid Abundance",
       subtitle = "Log2-transformed abundance data from LC-MS (+)",
       x = "mQTL Genotype", # the \n here allows an extra blank line to push the label lower
       y = "Chlorogenic Acid Abundance (log2)") +
  font("title", size = 14) +
  font("subtitle", size = 12) +
  font("xlab", size = 14) +
  font("ylab", size = 14) +
  font("xy.text", size = 12)
CGAmQTLgenotypeBoxplot
```

![](Pt3_mGWAS_Results_Processing_and_Visualization_files/figure-gfm/unnamed-chunk-101-1.png)<!-- -->
