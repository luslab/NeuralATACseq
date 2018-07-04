#/*==========================================================================#*/
#' ## Load and preprocess WT ATAC-seq data
#+ chunk_deseq2WT_load, cache=T, echo=F, warning=F, message=F
#/*==========================================================================#*/
regions.path <- file.path(results.path, 'atac_conregions')
sample.txt <- list.files(regions.path, pattern = 'samples', full.names = T)
count.txts <- list.files(regions.path, pattern = 'counts.txt', full.names = T)

# Read selected high quality samples
samples <- read_tsv(sample.txt, col_names = F) %>%
  select(X1, X2) %>%
  unlist() %>%
  sort()
names(samples) <- NULL

# Convert sample IDs from different runs to same format (rename replicates to 1/2)
new.samples <- gsub('-', '', samples) %>%
  tbl_df() %>%
  mutate(exp = gsub('[1-3]$', '', value),
         rep = str_extract(value, pattern = '[1-3]$')) %>%
  group_by(exp) %>%
  mutate(r = 1:length(exp),
         sample = sprintf('%s-%s', exp, r)) %>%
  ungroup() %>% select(sample) %>% unlist(.)

# Read con region counts by featureCounts
i <- grep(count.txts, pattern = '.*conPeaks_counts.txt$')
region.counts <- read_tsv(count.txts[i], skip = 1, col_names = T)
colnames(region.counts) <- gsub('_.*', '', basename(colnames(region.counts)))

# Parse regions from featureCount table
regions <- GRanges(seqnames = region.counts$Chr,
                   ranges = IRanges(region.counts$Start, region.counts$End),
                   id = region.counts$Geneid)

# Annotate regions with TSS dist/genomic features
summits <- getSummit(regions)
summits <- annotatePeak(summits, TxDb = gencode.txdb, annoDb = 'org.Mm.eg.db')
mcols(regions) <- cbind(mcols(regions), mcols(summits@anno))

# Filter low quality samples from count table
region.counts <- region.counts[,colnames(region.counts) %in% samples]
colnames(region.counts) <- new.samples

### Pre-filtering --> filter outlier regions
# Plot max read count per region distribution
max.counts <- apply(region.counts, 1, max)
data_frame(x = log2(max.counts)) %>%
  #ggplot(aes(x = x)) + geom_histogram(bins = 5*10^2) + xlab('log2(max counts)')
  ggplot(aes(x = x)) + geom_density() + xlab('log2(max counts)') +
  geom_vline(xintercept = 4, col = 'red', linetype = 'dashed')

# Filter extrem small and large valued regions
i.rmv <- which(max.counts < 16 | max.counts > 2500)
region.counts <- region.counts[-i.rmv,]
regions <- regions[-i.rmv,]

### Combine metadata, peaks and count matrix in DESeq2 object
# Create columnData
col.data <- data_frame(sample = new.samples) %>%
  mutate(sample = gsub('-', '_', sample),
         cond = gsub('_[1-3]', '', sample),
         rep = gsub('.*_', '', sample))
colnames(region.counts) <- col.data$sample %>% unlist()

# Create 'DESeqData' object
deseq.data <- DESeqDataSetFromMatrix(countData = region.counts,
                                     rowRanges = regions,
                                     colData = col.data,
                                     design = ~ cond)

#/*==========================================================================#*/
#' ## Diff. accessibility analysis with DESeq2
#+ chunk_deseq2WT_diffPeak, cache=T, echo=F, warning=F, message=F
#/*==========================================================================#*/
# DESeq2 output path
deseq2.path <- file.path(results.path, 'deseq2')
if(!dir.exists(deseq2.path)) dir.create(deseq2.path)

# RUN DESeq2
deseq.data <- DESeq(deseq.data, parallel = F)

# Get sig. upregulated regions during differentation process
result.names <- resultsNames(deseq.data)[-1]
i.diff <- sapply(result.names, function(n) {
  r <- results(deseq.data, name = n)
  which(r$padj < 0.01 & r$log2FoldChange > 1)
})

#/*==========================================================================#*/
#' ## DESeq2 analysis visualization
#+ chunk_deseq2WT_diffPeakPlots, cache=T, echo=F, warning=F, message=F
#/*==========================================================================#*/
### Plot pvalue distribution
pvalue.dists <- lapply(result.names, function(rn) {
  data_frame(comparison = rn, pvalue = results(deseq.data, name = rn)$pvalue)
})

gg.pvalue <- bind_rows(pvalue.dists) %>%
  tbl_df() %>%
  ggplot(aes(x = pvalue)) + geom_histogram(bins = 10^2) +
  facet_wrap(~comparison, scales = 'free', ncol = 4) +
  xlab('p-value')
gg.pvalue

#/*==========================================================================#*/
#' ### Number of differential peaks (Figure 1 F)
#+ chunk_deseq2WT_diffPeakPlots_nPeaks, cache=T, echo=F, warning=F, message=F
#/*==========================================================================#*/
# Preapre df with number of differential regions over conditons vs D0
n.diffRegions <- lapply(i.diff, length) %>%
  tbl_df() %>%
  gather(comp, n, 1:10) %>%
  mutate(sample = gsub('(cond_|_vs.*)', '', comp),
         sample = gsub('_', '-', sample),
         day = gsub('[A-Z]', '', sample),
         cond = gsub('D[0-9]', '', sample),
         cond = if_else(cond == '', 'epiblast', cond),
         group = cond,
         group = if_else(cond == 'NMP', 'SC', group))

# Plot number of differential regions in each comparison
gg.nDiffRegions <- n.diffRegions %>%
  mutate(sample = factor(sample,  levels = names(color.pal))) %>%
  ggplot(aes(x = sample, y = n, fill = sample)) + geom_bar(stat = 'identity') +
  xlab('') + ylab('Number of diff. regions (vs. D0)') +
  scale_fill_manual(values = color.pal) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = 'none')
gg.nDiffRegions

# Save number of diff regions barplots
#savePlots(gg.nDiffRegions, file.name = file.path(figures.path, 'nDiffRegionsAll_bar'))

#/*==========================================================================#*/
#' ### Differential peaks - MA-plot (Figure 1 E)
#+ chunk_deseq2WT_diffPeakPlots_MAPlot, cache=T, echo=F, warning=F, message=F
#/*==========================================================================#*/
# Plot MA-plots for all comparison
gg.maplots <- lapply(result.names, function(rn){
  # Prepare labels for MAplot
  cond <- str_split(rn, pattern = '_') %>% unlist(.)
  cond <- cond[c(2,4)]
  nDiff.df <- results(deseq.data, name = rn) %>%
    tbl_df() %>%
    filter(padj <= 0.01 & abs(log2FoldChange) > 1) %>%
    mutate(sample = if_else(log2FoldChange > 0, cond[1], cond[2])) %>%
    group_by(sample) %>%
    summarise(n = n(),
              x = log2(max(baseMean)) - 0.5,
              y = max(abs(log2FoldChange))) %>%
    mutate(x = round(max(x), 1),
           y = if_else(sample == 'D0', -y, y),
           label = sprintf('%s regions (n = %s)', sample, n))

  # Plot MA-plot
  gg.maplot <- results(deseq.data, name = rn) %>%
    tbl_df() %>%
    mutate(padj_cat = '> 0.1',
           padj_cat = if_else(padj < 0.1, '< 0.1', padj_cat),
           padj_cat = if_else(padj < 0.01, '< 0.01', padj_cat),
           padj_cat = if_else(padj < 0.001, '< 0.001', padj_cat)) %>%
    ggplot(aes(x = log2(baseMean), y = log2FoldChange, col = padj_cat)) +
    geom_point() +
    geom_hline(yintercept = 0, col = 'red', linetype = 'dashed') +
    scale_color_manual(values = rev(c('grey', '#fee5d9', '#fb6a4a', '#a50f15')), name = 'FDR') +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 4))) +
    geom_text(data = nDiff.df, aes(x = x, y = y, label = label), inherit.aes = F) +
    xlab('log2(mean)') + ylab('log2(FC)') + #  theme_cowplot(font_size = 10)
    theme(legend.position = 'bottom')
    #theme(legend.position = c(0.75, 0.15))
  # Add histogram to MAplot
  gg.maplot <- ggExtra::ggMarginal(gg.maplot, type = 'histogram')
  return(gg.maplot)
})
names(gg.maplots) <- result.names

# Plot MA-plots in report
gg.maplots[['cond_D5SC_vs_D0']]

# Save MA-plots
#lapply(result.names, function(rn) {
#  maplot.file <-  file.path(figures.path, sprintf('%s_MAplot', rn))
#  savePlots(gg.maplots[[rn]], file.name = maplot.file, base_height = 5, base_width = 6,
#            figure.suffix = c('pdf', 'eps', 'png'))
#})

#/*==========================================================================#*/
#' ## Post-process DESeq2 results
#+ chunk_deseq2WT_postprocess, cache=T, echo=F, warning=F, message=F
#/*==========================================================================#*/
### Get ALL variable regions during diff process
# Combine all sig differential regions over the differentation process
i.diff <- sort(unique(unlist(i.diff)))
#length(i.diff)

# Get log2(normalised counts + 1) of the variable regions
#norm.counts <- log2(counts(deseq.data, normalized = T)[i.diff,] + 1)
norm.counts <- counts(deseq.data, normalized = T)[i.diff,]

# Get index of ESC specific regions
i.diffEsc <- lapply(result.names[grep(result.names, pattern = 'D5')], function(rn) {
  r <- results(deseq.data, name = rn)
  which(r$padj < 0.01 & r$log2FoldChange < -1)
})

# Write ESC specific sites to *.bed file
i.diffEsc <- unique(unlist(i.diffEsc))
diffRegions.esc <- regions[i.diffEsc,] %>% tbl_df(.)
#write_tsv(x = diffRegions.esc[,1:3],
#          path = file.path(deseq2.path, 'ESC_diffRegions_D0.bed'))
