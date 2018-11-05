#' ---
#' title: 'Analysis scripts: Regionalization of the nervous system requires axial allocation prior to neural lineage commitment, Metzis et. al. (2018)'
#' author: 'Sebastian Steinhauser'
#' date: '14/07/2017'
#' output:
#'   html_document:
#'     number_sections: yes
#'     toc: true
#'     toc_float: true
#'     fig_caption: yes
#'     code_folding: hide
#'   pdf_document:
#'     number_sections: yes
#'     toc: true
#'     fig_caption: yes
#' ---

#/*==========================================================================#*/
#' # Libraries, paths and parameters
#+ chunk_preparation, cache=F, error=F, warning=F, message=F, echo=T
#/*==========================================================================#*/
# Libraries
library(GenomicRanges)
library(BSgenome.Mmusculus.UCSC.mm10)
library(blutwuRst)
library(DESeq2)
library(GenomicFeatures)
library(GenomicInteractions)
library(org.Mm.eg.db)
library(tximport)
library(ChIPseeker)
library(LOLA)
library(rGREAT)
library(motifmatchr)
library(TFBSTools)
library(JASPAR2016)
library(CENTIPEDE)
library(limma)
library(oligo)
library(mouse4302.db)

library(som)

library(dplyr)
library(readr)
library(tidyr)
library(purrr)
library(stringr)
library(broom)
library(knitr)

library(ggplot2)
library(RColorBrewer)
library(wiggleplotr)
library(ggrepel)
library(ggpubr)
library(ggsci)
library(cowplot)
library(ComplexHeatmap)

# Define paths
project.path <- '/Users/steinhs/projects/atac_vicki/'
data.path <- file.path(project.path, 'data')
atac.path <- file.path(data.path, 'atacseq/vicki_final')
rnaseq.path <- file.path(data.path, 'rnaseq')
results.path <- file.path(project.path, 'results')
figures.path <- file.path(project.path, 'figs')
src.path <- file.path(project.path, 'src')
notebook.path <- file.path(project.path, 'src/publication_scripts/metzis2018')

# Additional parameters
set.seed(1234)
figure.suffix <- c('pdf', 'eps')

## Define color code for sample names
color.pal <- c('D0' = '#B4B4B4', 'D1' = '#B4B4B4', 'D2' = '#B4B4B4',
               'D2.5' = '#B4B4B4',
               'D3A' = '#A7DDFF', 'D4A' = '#74ADD1', 'D5A' = '#74ADD1',
               'D4H' = '#FDCD61', 'D5H' = '#FDCD61', 'D2.5NMP' = '#F79992',
               'D3NMP' = '#F79992', 'D4SC' = '#FB564A', 'D5SC' = '#FB564A',
               'D5A-invivo' = '#484F9E', 'D5SC-invivo' = '#AE2C22',
               'D3NMP-BRA' = '#9957BE', 'D4SC-BRA' = '#9957BE', 'D5SC-BRA' = '#9957BE',
               'D3NMP-CDX' = '#39C051', 'D4SC-CDX' = '#39C051', 'D5SC-CDX' = '#39C051',
               'D5H+' = 'black',
               'D5SCind' = '#b784a7', 'D5Hrep' =  '#E29ECD')

shape.pal <- c('D0' = 21, 'D1' = 23, 'D2' = 24,
               'D3NMP' = 21, 'D4SC' = 21, 'D5SC' = 24,
               'D3A' = 21, 'D4H' = 21, 'D5H' = 24,
               'D4A' = 21, 'D5A' = 24)

bioCluster.color <- c('Anterior' = '#74ADD1', 'Hindbrain' = '#FDCD61',
                      'Spinal cord' = '#FB564A', 'Neural' = 'black',
                      'ESC'= '#B4B4B4')

### Load GENCODE annotation as TxDB file
gencode.gtf <- '/Users/steinhs/genomes/mm10/gencode.vM14.annotation.gtf.gz'
gencode.rdata <- '/Users/steinhs/genomes/mm10/gencode.vM14.annotation.RData'

if(!file.exists(gencode.rdata)) {
  gencode.txdb <- makeTxDbFromGFF(file = gencode.gtf,
                                  format = 'gtf',
                                  dataSource = 'gencode.vM14.annotation',
                                  organism = 'Mus musculus',
                                  chrominfo = seqinfo(BSgenome.Mmusculus.UCSC.mm10))
  saveDb(gencode.txdb, file = gencode.rdata)
} else {
  gencode.txdb <- loadDb(gencode.rdata)
  gencode.metadata <- import.gff(gencode.gtf)
  gencode.metadata <- mcols(gencode.metadata) %>%
    tbl_df %>%
    #filter(type == 'gene') %>%
    select(transcript_id, gene_id, gene_type, gene_name, exon_number) %>%
    mutate(strand = as.character(strand(gencode.metadata)))
}

### Load atacR interaction data
if(!exists('interactions')) {
  interactions <- loadInteractionDb('mm10')
}

#/*==========================================================================#*/
#' # Load RNA-seq data from SALMON files
#+ chunk_loadRNA, cache=F, echo=F, warning=F, message=F
#/*==========================================================================#*/
# Define *.Rds file to save SALMON RNAseq quantification
rnaseq.rds <- file.path(results.path, 'RNAseq_txiSalmon.Rds')

# Load SALMON RNA-seq quantification
if(!file.exists(rnaseq.rds)) {
  # Get transcript id to GeneSymbol table
  keys <- keys(gencode.txdb, keytype = 'GENEID')
  df <- AnnotationDbi::select(gencode.txdb, keys = keys, keytype = 'GENEID',
                              columns = c('TXNAME'))
  tx2gene <- df[, 2:1]  # tx ID, then gene ID

  # Load GENCODE metadata from *.gtf file
  gencode.gtf <- '/Users/steinhs/genomes/mm10/gencode.vM14.annotation.gtf.gz'
  gencode.gr <- rtracklayer::import.gff(gencode.gtf)
  gencode.metadata <- mcols(gencode.gr) %>%
    tbl_df() %>%
    select(gene_id, transcript_id, gene_name, exon_number) %>%
    mutate(strand = as.character(strand(gencode.gr)))

  # Find salmon files
  salmon.files <- atacR::listFiles(rnaseq.path, pattern = 'quant.sf') %>% sort(.)
  salmon.files <- salmon.files[-grep(salmon.files, pattern = 'Amin')]
  sample.names <- lapply(str_split(salmon.files, pattern = '/'), function(x) rev(x)[3]) %>% unlist()
  sample.names[grep(sample.names, pattern = 'ESC')] <- c('D0_1', 'D0_2')

  # Read salmon files
  txi.salmon <- tximport(salmon.files, type = 'salmon', tx2gene = tx2gene)

  # Assign sample names as colnames
  colnames(txi.salmon$abundance) <- sample.names
  colnames(txi.salmon$counts) <- sample.names

  # Save *.Rds
  saveRDS(txi.salmon, file = rnaseq.rds)
}else{
  # Load RNA-seq quantification (by SALMON) from *.Rds file
  txi.salmon <- readRDS(rnaseq.rds)
}

# Get TPMS from import
tpms <- txi.salmon$abundance
rownames(tpms) <- gsub('\\..*', '', rownames(tpms))

#/*==========================================================================#*/
#' # ATAC-seq quality control
#+ chunk_qc, cache=F, echo=F, warning=F, message=F
#/*==========================================================================#*/
spin_child(file.path(notebook.path, 'analysis/01_atacseq_qc.R'))

#/*==========================================================================#*/
#' # Differential analysis of WT ATAC-seq using DESeq2 and SOM clustering
#+ chunk_deseq2WT, cache=F, echo=F, warning=F, message=F
#/*==========================================================================#*/
spin_child(file.path(notebook.path, 'analysis/02_atacseq_DESeq2-WT-analysis.R'))
spin_child(file.path(notebook.path, 'analysis/03_atacseq_SOM-WT-analysis.R'))

#/*==========================================================================#*/
#' # Analysis of Neuro-mesodermal progentior (NMP) specific regions
#+ chunk_NMP, cache=F, echo=F, warning=F, message=F
#/*==========================================================================#*/

#/*==========================================================================#*/
#' # Analysis of day5 hindbrain repressed and spinal cord induced
#+ chunk_d5RepInd, cache=F, echo=F, warning=F, message=F
#/*==========================================================================#*/
#spin_child(file.path(notebook.path, 'xx_D5SCind-D5Hrep_analysis.R'))

#/*==========================================================================#*/
#' # Session Info
#+ chunk_session_info, cache=F, echo=F, warning=F, message=F, tidy=T
#/*==========================================================================#*/
sessionInfo()
