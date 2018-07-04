#/*==========================================================================#*/
#' ## Load and preprocess D5H-repressed and D5SC-induced ATAC-seq data
#+ chunk_d5RepInd_load, cache=T, echo=F, warning=F, message=F
#/*==========================================================================#*/
# Define paths and find files
regions.path <- file.path(results.path, 'atac_D5Hrep-D5SCind/conPeaks')
d5HrepSCind.path <- file.path(figures.path, 'D5Hrep-D5SCind')
sample.txt <- list.files(regions.path, pattern = 'samples', full.names = T)
count.txts <- list.files(regions.path, pattern = 'counts.txt', full.names = T)

# Read selected high quality samples
samples <- read_tsv(sample.txt, col_names = F) %>%
  select(X1, X2) %>%
  unlist() %>%
  sort()
names(samples) <- NULL
samples <- gsub('_rep', '-', samples)

# Convert sample IDs from different runs to same format (rename replicates to 1/2)
new.samples <- gsub('-', '', samples) %>%
  tbl_df() %>%
  mutate(exp = gsub('[1-3]$', '', value),
         rep = str_extract(value, pattern = '[1-3]$')) %>%
  group_by(exp) %>%
  mutate(r = 1:length(exp),
         sample = sprintf('%s-%s', exp, r),
         sample = gsub('_rep', '', sample)) %>%
  ungroup() %>% select(sample) %>% unlist(.)

# Read con region counts by featureCounts
i <- grep(count.txts, pattern = '.*conPeaks_counts.txt$')
region.counts <- read_tsv(count.txts[i], skip = 1, col_names = T)
colnames(region.counts) <- gsub('_rep', '-', gsub('_rmbqr.*', '', basename(colnames(region.counts))))

# Parse regions from featureCount table
regions <- GRanges(seqnames = region.counts$Chr,
                   ranges = IRanges(region.counts$Start, region.counts$End),
                   id = region.counts$Geneid)

# Filter low quality samples from count table
region.counts <- region.counts[,colnames(region.counts) %in% samples]

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
col.data <- data_frame(sample = colnames(region.counts)) %>%
  mutate(sample = gsub('-', '_', sample),
         cond = gsub('_[1-3]', '', sample),
         rep = gsub('.*_', '', sample),
         cond = gsub('_', '', cond))
colnames(region.counts) <- col.data$sample %>% unlist()

# Create 'DESeqData' object
deseq.data <- DESeqDataSetFromMatrix(countData = region.counts,
                                     rowRanges = regions,
                                     colData = col.data,
                                     design = ~ cond)

#/*==========================================================================#*/
#' ## Comparison of counts on old consensus peak set D5H/SC vs D5Hrep/D5SCind (Figure 5D,H)
#+ chunk_d5RepInd_comparison, cache=T, echo=F, warning=F, message=F
#/*==========================================================================#*/
# Normalise count data
deseq.data <- estimateSizeFactors(deseq.data)
mcols(regions) <- counts(deseq.data, normalized = T)

# Overlap diffRegions with D5Hrep and D5SCind regions
diff.regions <- readRDS('/Users/steinhs/projects/atac_vicki/results/diffRegions.Rds')
ov.diff <- findOverlaps(regions, diff.regions, type = 'equal')
regions <- regions[queryHits(ov.diff),]
regions$bio_cluster <- diff.regions$bio_cluster[subjectHits(ov.diff)]

# Summarise replicates and format sample names
region.counts <- mcols(regions) %>%
  tbl_df %>%
  select(contains('D5'), 'bio_cluster') %>%
  mutate(peak_id = 1:length(regions)) %>%
  gather(sample, counts, 1:10) %>%
  mutate(sample = if_else(str_count(sample, pattern = '_') == 2, sub('_', '', sample), sample)) %>%
  separate(sample, c('cond', 'rep'), sep = '_') %>%
  group_by(peak_id, bio_cluster, cond) %>%
  summarise(counts = mean(counts)) %>%
  mutate(counts = log2(counts)) %>%
  spread(cond, counts)

# List of sample comparisons
rc.comps <- list(c('D5Hrep', 'D5H'),
                 c('D5SCind', 'D5H'),
                 c('D5SC', 'D5H'),
                 c('D5Hrep', 'D5SC'),
                 c('D5SCind', 'D5SC'),
                 c('D5SCind', 'D5Hrep'))

# Plot spinal cord/hindbrain scatterplots
lapply(rc.comps, function(rc.comp){
  print(rc.comp)
  gg.rcComp <- region.counts %>%
    filter(bio_cluster %in% c('Hindbrain', 'Spinal cord')) %>%
    ggplot(aes_string(x = rc.comp[2], y = rc.comp[1], col = 'bio_cluster')) +
    geom_point() + geom_abline(col = 'red', linetype = 'dashed') +
    stat_cor(show.legend = F) +
    scale_color_manual(values = bioCluster.color, name = 'Cluster') +
    xlim(0, 11) + ylim(0, 11) +
    xlab(sprintf('%s log2(norm. counts)', rc.comp[2])) +
    ylab(sprintf('%s log2(norm. counts)', rc.comp[1])) +
    theme(legend.position = 'bottom')
  gg.rcComp
})
