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

* Figure 1
 [C](../master/analysis/regulatoryRegionPlots/01_Pou5f1Example.R),
 [D](../master/analysis/regulatoryRegionPlots/02_0lig2Example.R),
 [E](../master/analysis/02_atacseq_DESeq2-WT-analysis.R),
 [F](../master/analysis/02_atacseq_DESeq2-WT-analysis.R),
 [G](../master/analysis/04_atacseq_WT-downstreamAnalysis.R)

* Figure 2
 [A](../master/analysis/03_atacseq_SOM-WT-analysis.R),
 [B](../master/analysis/regulatoryRegionPlots/03_ShhAnteriorExample.R),
 [C](../master/analysis/regulatoryRegionPlots/04_Phox2bHindbrainExample.R),
 [E](../master/analysis/regulatoryRegionPlots/05_Hoxc8SpinalCordExample.R),
 [F](../master/analysis/04_atacseq_WT-downstreamAnalysis.R),
 [G](../master/analysis/04_atacseq_WT-downstreamAnalysis.R),
 [H](../master/analysis/04_atacseq_WT-downstreamAnalysis.R),
 [I](../master/analysis/04_atacseq_WT-downstreamAnalysis.R),
 [J](../master/analysis/04_atacseq_WT-downstreamAnalysis.R),

* Figure 3
 [A](../master/analysis/04_atacseq_WT-downstreamAnalysis.R),
 [B](../master/analysis/04_atacseq_WT-downstreamAnalysis.R),
 [C](../master/analysis/04_atacseq_WT-downstreamAnalysis.R),
 [D](../master/analysis/04_atacseq_WT-downstreamAnalysis.R),
 [P](../master/analysis/xx_rnaseq_analysis.R)

* Figure 4
 [A](../master/analysis/05_atacseq_NMP-analysis.R),
 [B](../master/analysis/05_atacseq_NMP-analysis.R),
 [C](../master/analysis/05_atacseq_NMP-analysis.R),
 [D](../master/sh/plotCdx2Heatmap.sh),
 [I](../master/analysis/04_atacseq_WT-downstreamAnalysis.R),
 [J](../master/analysis/07_atacseq_Cdx2-analysis.R),
 [N](../master/analysis/05_atacseq_NMP-analysis.R)

* Figure 5
 [D](../master/analysis/xx_D5SCind-D5Hrep_analysis.R),
 [H](../master/analysis/xx_D5SCind-D5Hrep_analysis.R)

* Figure S1
 [A](../master/analysis/01_atacseq_qc.R),
 [B](../master/analysis/01_atacseq_qc.R),
 [C](../master/analysis/01_atacseq_qc.R),
 [D](../master/analysis/01_atacseq_qc.R),
 [E](../master/analysis/01_atacseq_qc.R),
 [F](../master/analysis/01_atacseq_qc.R)

* Figure S2
 [A](../master/analysis/04_atacseq_WT-downstreamAnalysis.R),
 [B](../master/analysis/04_atacseq_WT-downstreamAnalysis.R),
 [C](../master/analysis/04_atacseq_WT-downstreamAnalysis.R),
 [D](../master/analysis/04_atacseq_WT-downstreamAnalysis.R),
 [E](../master/analysis/04_atacseq_WT-downstreamAnalysis.R),
 [F](../master/analysis/04_atacseq_WT-downstreamAnalysis.R),
 [G](../master/analysis/04_atacseq_WT-downstreamAnalysis.R),
 [H](../master/analysis/04_atacseq_WT-downstreamAnalysis.R),
 [I](../master/analysis/04_atacseq_WT-downstreamAnalysis.R),
 [J](../master/analysis/regulatoryRegionPlots/06_ShhInvivoExample.R),
 [K](../master/analysis/regulatoryRegionPlots/07_Olig2InvivoExample.R)

* Figure S3
 [A](../master/analysis/04_atacseq_WT-downstreamAnalysis.R),
 [B](../master/analysis/04_atacseq_WT-downstreamAnalysis.R),
 [C](../master/analysis/04_atacseq_WT-downstreamAnalysis.R)

* Figure S4
 [A](../master/analysis/xx_atacseq_hbPlus-analysis.R),
 [B](../master/analysis/xx_atacseq_hbPlus-analysis.R),
 [C](../master/analysis/xx_atacseq_hbPlus-analysis.R),
 [F](../master/analysis/xx_rnaseq_analysis.R)

* Figure S5
 [A](../master/analysis/07_atacseq_Cdx2-analysis.R),
 [A'](../master/analysis/regulatoryRegionPlots/08_Phox2bCdx2Example.R),
 [A''](../master/analysis/regulatoryRegionPlots/09_MafbCdx2Example.R),
 [B](../master/sh/plotCdx2Heatmap.sh),
 [C](../master/analysis/07_atacseq_Cdx2-analysis.R),
 [D](../master/analysis/07_atacseq_Cdx2-analysis.R),
 [E](../master/analysis/07_atacseq_Cdx2-analysis.R),
 [F](../master/analysis/07_atacseq_Cdx2-analysis.R),
 [G](../master/analysis/xx_microarray_analysis.R),
 [H](../master/analysis/07_atacseq_Cdx2-analysis.R)

## Additional information 

* Pre-processing pipelines are available upon request.

* Processed data (*.bam, *.bw, etc.) are available upon request.

## Reference 
**Regionalization of the nervous system requires axial allocation prior to neural lineage commitment**

[biorxiv preprint](https://www.biorxiv.org/content/early/2017/12/04/229203)
