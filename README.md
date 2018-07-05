# Regionalization of the nervous system requires axial allocation prior to neural lineage commitment 

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

* Figure 1 [C](../master/analysis/regulatoryRegionPlots/01_Pou5f1Example.R), [D](../master/analysis/regulatoryRegionPlots/02_0lig2Example.R),
[E](../master/analysis/02_atacseq_DESeq2-WT-analysis.R), [F](../master/analysis/02_atacseq_DESeq2-WT-analysis.R), [G](../master/analysis/04_atacseq_WT-downstreamAnalysis.R)

* Figure 2
A
B
C
D
F
G
H
I
J

* Figure 3
A
B
C
D
P

* Figure 4
A
B
C
D
I
J
N

* Figure 5
D
H

* Figure S1
A
B
C
D
E
F

* Figure S2
A
B
C
D
E
F
G
H
I
J
K

* Figure S3
A
B
C

* Figure S4
A
B
C
F

* Figure S5 [A](../master/analysis/07_atacseq_Cdx2-analysis.R), [A'](../master/analysis/regulatoryRegionPlots/08_Phox2bCdx2Example.R), [A''](../master/analysis/regulatoryRegionPlots/09_MafbCdx2Example.R), [B](../master/sh/plotCdx2Heatmap.sh), [C](../master/analysis/07_atacseq_Cdx2-analysis.R), [D](../master/analysis/07_atacseq_Cdx2-analysis.R), [E](../master/analysis/07_atacseq_Cdx2-analysis.R), [F](../master/analysis/07_atacseq_Cdx2-analysis.R), [G](../master/analysis/xx_microarray_analysis.R), [H](../master/analysis/07_atacseq_Cdx2-analysis.R)


## Reference 
**Regionalization of the nervous system requires axial allocation prior to neural lineage commitment**

[biorxiv preprint](https://www.biorxiv.org/content/early/2017/12/04/229203)
