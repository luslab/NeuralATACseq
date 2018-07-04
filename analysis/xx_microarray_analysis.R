#/*==========================================================================#*/
#' # Microarray analysis: MNP WT vs Hindbrain Cdx2 induced
#+ chunk_microarray, cache=T, echo=F, warning=F, message=F
#/*==========================================================================#*/
# Define path to microarrys
microarray.path <- file.path(data.path, 'microarray')
cel.files <- list.files(microarray.path, pattern = 'Cdx2', full.names = T)

# Real CEL files
mnp.exp <- read.celfiles(cel.files)

# Normalise using RMA
norm.mnpExp <- rma(mnp.exp)

# Create design matrix
f <- colnames(norm.mnpExp) %>% basename %>%
  gsub('.*_iC', 'C', .) %>% gsub('5_.*', '5', .) %>%
  factor
design <- model.matrix(~0+f)
colnames(design) <- levels(f)

# Fit linear model with LIMMA
fit <- lmFit(exprs(norm.mnpExp), design = design)

# Create constrast matrix
cont.matrix <- makeContrasts(WTvsCdxInd = 'Cdx2_RAHh_Day5-Cdx2_RAHhDoxFGF_Day5', levels = design)

# Extract the linear model fit for the contrasts
fit2  <- contrasts.fit(fit, cont.matrix)
fit2  <- eBayes(fit2)

# Extract LIMMA results
mnpCdx.diffGenes <- topTable(fit2, coef = 1, adjust = 'fdr', number = nrow(norm.mnpExp))

# Map affy ids to gene names
affy.symbols <- as.list(mouse4302SYMBOL) %>% unlist

affy.symbols <- affy.symbols %>%
  tbl_df %>%
  mutate(affy_id = names(affy.symbols)) %>%
  rename(gene_name = value)

# Assign gene_names to affyIds
mnpCdx.diffGenes <- mnpCdx.diffGenes %>%
  tbl_df() %>%
  mutate(affy_id = rownames(mnpCdx.diffGenes)) %>%
  left_join(., affy.symbols, by = 'affy_id') %>%
  filter(!is.na(gene_name))

# Write diff genes to file to run beta
#write_tsv(mnpCdx.diffGenes[,c(8, 1, 5)], path = file.path(results.path, 'MNP_CdxDiffGenes.txt'), col_names = F)
#/*==========================================================================#*/
#' ## Differential gene expression analysis (DESeq2)
#+ chunk_rnaseq_deseq2, cache=T, echo=F, warning=F, message=F
#/*==========================================================================#*/
# Prepare DESeqData object
sample.names <- colnames(tpms)
txi.colData <- data_frame(sample = sample.names) %>%
  separate(sample, c('cond', 'rep'), sep = '_', remove = F)
deseq.rna <- DESeqDataSetFromTximport(txi.salmon, colData = txi.colData, design = ~ cond)
deseq.rna <- deseq.rna[,grep(txi.colData$cond, pattern = '(D5SC|D5H)')]
colData(deseq.rna)  <- colData(deseq.rna) %>% droplevels()

# Perform diff analysis with DESeq2
deseq.rna <- DESeq(deseq.rna)

# Get and annotate results
deseqRna.results <- results(deseq.rna) %>%
  tbl_df() %>%
  mutate(gene_id = rownames(deseq.rna)) %>%
  left_join(., gencode.metadata, by = 'gene_id') %>%
  select(-transcript_id, -exon_number, -strand) %>%
  unique() %>%
  rename(logFC_rnaseq = log2FoldChange)

#/*==========================================================================#*/
#' # Compare log2(WT/Cdx2 induced) with other published gene expression data (Figure S5 G)
#+ chunk_microarray_comp, cache=T, echo=F, warning=F, message=F
#/*==========================================================================#*/
mnpCdx.cor <- deseqRna.results %>%
  filter(padj <= 0.1) %>%
  left_join(., mnpCdx.diffGenes, by = 'gene_name') %>%
  filter(!is.na(logFC))

mnpCdx.label <- mnpCdx.cor %>%
  filter(abs(logFC) > 3.5 | abs(logFC_rnaseq) > 4)

gg.mnpCdx <- mnpCdx.cor %>%
  mutate(padj = -log10(padj)) %>%
  ggplot(aes(x = logFC, y =logFC_rnaseq, col= padj)) +  geom_point() +
  xlab('Microarray log2(WT/Cdx2 induced)') + ylab('RNA-seq log2(D5SC/D5H)') +
  geom_text_repel(data = mnpCdx.label, aes(label = gene_name), col = 'black') +
  scale_color_gradientn(colours = c('#fff5f0', '#fee0d2', '#fcbba1', '#fc9272',
                                    '#fb6a4a', '#ef3b2c', '#cb181d', '#a50f15',
                                    '#67000d'), name = '-log10(adj pvalue)') +
  geom_hline(yintercept = 0, col = 'black', linetype = 'dashed') +
  geom_vline(xintercept = 0, col = 'black', linetype = 'dashed') +
  theme(legend.position = 'bottom')
gg.mnpCdx

#savePlots(gg.mnpCdx, file.name = file.path(figures.path, 'HBvsSC_WTvsCdxInduced'),
#          base_height = 8, base_width = 9)
