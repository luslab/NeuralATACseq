# Metzis et.al. (2018) analysis scripts

## Abstract 
Neural induction in vertebrates generates a central nervous system that extends the rostral-caudal length of the body. The prevailing view is that neural cells are initially induced with anterior (forebrain) identity, with caudalising signals then converting a proportion to posterior fates (spinal cord). To test this model, we used chromatin accessibility assays to define how cells adopt region-specific neural fates. Together with genetic and biochemical perturbations this identified a developmental time window in which genome-wide chromatin remodeling events preconfigure epiblast cells for neural induction. Contrary to the established model, this revealed that cells commit to a regional identity before acquiring neural identity. This 'primary regionalization' allocates cells to anterior or posterior regions of the nervous system, explaining how cranial and spinal neurons are generated at appropriate axial positions. These findings prompt a revision to models of neural induction and support the proposed dual evolutionary origin of the vertebrate central nervous system.

## Requirements 

* Install following R packages from CRAN: 
som, dplyr, parallel, readr, tidyr, purrr, stringr, broom, knitr,
ggplot2, RColorBrewer, ggrepel, ggpubr, ggsci, cowplot

* Install following R packages from bioconductor:
GenomicRanges, BSgenome.Mmusculus.UCSC.mm10, DESeq2, GenomicFeatures, GenomicInteractions, org.Mm.eg.db, tximport, ChIPseeker, LOLA, rGREAT, motifmatchr, TFBSTools, JASPAR2016, limma, oligo, mouse4302.db, wiggleplotr, ComplexHeatmap

* Install following R packages from R-forge: 
CENTIPEDE

* Install atacR from github:

    atacR contains useful functions to analyse (bulk) ATAC-seq data with R.

    ``` r
    devtools::install_github("luslab/atacR")
    ``` 


## Figure to code map 

* Figure 1 [C](../analysis/regulatoryRegionPlots/01_Pou5f1Example.R), [D](../master/analysis/regulatoryRegionPlots/02_0lig2Example.R),
[E](../master/analysis/02_atacseq_DESeq2-WT-analysis.R), [F](../master/analysis/02_atacseq_DESeq2-WT-analysis.R), [G](../master/analysis/04_atacseq_WT-downstreamAnalysis.R)

* Figure 2

* Figure 3

* Figure 4

* Figure 5

* Figure S1

* Figure S2

* Figure S3

* Figure S4

* Figure S5


## Reference 
**Regionalization of the nervous system requires axial allocation prior to neural lineage commitment**

[biorxiv preprint](https://www.biorxiv.org/content/early/2017/12/04/229203,
