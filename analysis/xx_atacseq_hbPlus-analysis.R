#/*==========================================================================#*/
#' ## Load and preprocess D5/D5-plus ATAC-seq data
#+ chunk_hbPlus_load, cache=T, echo=F, warning=F, message=F
#/*==========================================================================#*/
# Define path to count table
regions.path <- file.path(results.path, 'atac_conregions_hbPlus')
sample.txt <- list.files(regions.path, pattern = 'samples', full.names = T)
count.txts <- list.files(regions.path, pattern = 'counts.txt', full.names = T)

# Read con region counts by featureCounts
i <- grep(count.txts, pattern = '.*conPeaks_counts.txt$')
region.counts <- read_tsv(count.txts[i], skip = 1, col_names = T)
colnames(region.counts) <- gsub('_.*', '', gsub('_rep', '-', basename(colnames(region.counts))))

# Parse regions from featureCount table
regions <- GRanges(seqnames = region.counts$Chr,
                   ranges = IRanges(region.counts$Start, region.counts$End),
                   id = region.counts$Geneid)

## Filter for D5H and D5H-plus
region.counts <- region.counts[,grep(colnames(region.counts), pattern = 'D')]
colnames(region.counts) <- gsub('-', '', colnames(region.counts))

### Pre-filtering --> filter outlier regions
# Plot max read count per region distribution
max.counts <- apply(region.counts, 1, max)
data_frame(x = log2(max.counts)) %>%
  ggplot(aes(x = x)) + geom_density() + xlab('log2(max counts)') +
  geom_vline(xintercept = 4, col = 'red', linetype = 'dashed')

# Filter extrem small and large valued regions
i.rmv <- which(max.counts < 16 | max.counts > 2500)
region.counts <- region.counts[-i.rmv,]
regions <- regions[-i.rmv,]

### Combine metadata, peaks and count matrix in DESeq2 object
# Create columnData
col.data <- data_frame(sample = colnames(region.counts)) %>%
  mutate(cond = gsub('[1-3]$', '', sample),
         rep = str_extract(sample, pattern = '[1-3]$'))

# Create 'DESeqData' object
deseq.data <- DESeqDataSetFromMatrix(countData = region.counts,
                                     rowRanges = regions,
                                     colData = col.data,
                                     design = ~ cond)

# Estimate size factors and normalise counts
deseq.data <- estimateSizeFactors(deseq.data)
normCounts.all<- counts(deseq.data, normalized = T)

#/*==========================================================================#*/
#' ## Diff. accessibility analysis with DESeq2 between D5H vs D5Hplus
#+ chunk_hbPlus_deseq2, cache=T, echo=F, warning=F, message=F
#/*==========================================================================#*/
# DESeq2 output path
deseq2.path <- file.path(regions.path, 'deseq2')
if(!dir.exists(deseq2.path)) dir.create(deseq2.path)

# Filter for samples of interest
deseq.data <- deseq.data[,grep(colData(deseq.data)$cond, pattern = 'D5H')]
colData(deseq.data)$cond <- droplevels(colData(deseq.data)$cond)

# RUN DESeq2
deseq.data <- DESeq(deseq.data, parallel = F)

# Get DESeq2 results
result.names <- resultsNames(deseq.data)[-1]
deseq.results <- results(deseq.data)

# Get index for diff regions
i.diff <- deseq.results %>%
  tbl_df %>%
  mutate(region_id = 1:length(padj)) %>%
  filter(!is.na(padj) & padj < 0.1 & abs(log2FoldChange) > 1) %>%
  select(region_id) %>% unlist

# Define differential regions
diff.regions <- regions[i.diff,]

# Get normalised counts of the variable regions
norm.counts <- counts(deseq.data, normalized = T)[i.diff,] %>%
  tbl_df %>%
  mutate(region_id = i.diff) %>%
  gather(sample, count, 1:4) %>%
  mutate(sample = gsub('[1-3]$', '', sample)) %>%
  group_by(region_id, sample) %>%
  summarise(count = mean(count)) %>%
  spread(sample, count) %>%
  mutate(logFc = log2(D5Hplus/D5H))

# Add normalised counts to diffRegions
mcols(diff.regions) <- cbind(mcols(diff.regions), norm.counts[,-1])

# Add cluster annotation
diff.regions$cluster <- 'D5H'
diff.regions$cluster[diff.regions$logFc > 0] <- 'D5Hplus'

#/*==========================================================================#*/
#' ## Plot D5Hplus scatterplot for Hindbrain/Spinalcord regions (Figure S4A)
#+ chunk_hbPlus_scatter, cache=T, echo=F, warning=F, message=F
#/*==========================================================================#*/
# Load SOM regions from *.Rds
som.regions <- readRDS(file.path(results.path, 'diffRegions.Rds'))

# Overlap all regions  with SOM regions
ov.all <- findOverlaps(regions, som.regions)

# Log2-foldchange over SC
i.counts <- grep(colnames(normCounts.all), pattern = '(D5H(1|3)$|D5SC(1|3)$|D5Hplus)')

# Plot log2(D5H/SC) vs log2(D5Hplus/D5SC)
gg.fcScatter <- normCounts.all[queryHits(ov.all), i.counts] %>%
  tbl_df %>%
  mutate(id = 1:length(ov.all),
         bio_cluster = som.regions$bio_cluster[subjectHits(ov.all)],
         is_cdx2 = som.regions$is_cdx2[subjectHits(ov.all)]) %>%
  gather(sample, count, 1:4) %>%
  mutate(cond = gsub('(1|2|3)$', '', sample)) %>%
  group_by(id, bio_cluster, cond) %>%
  summarise(count = mean(count)) %>%
  spread(cond, count) %>%
  mutate(logFC = log2(D5H/D5SC),
         logFC_plus = log2(D5Hplus/D5SC)) %>%
  filter(bio_cluster %in% c('Hindbrain', 'Spinal cord')) %>%
  ggplot(aes(x = logFC, y = logFC_plus, col = bio_cluster)) + geom_point() +
  #geom_point(shape = 16, alpha = 0.5, size = 0.6) +
  geom_hline(yintercept = 0, col = 'black', linetype = 'dashed') +
  geom_vline(xintercept = 0, col = 'black', linetype = 'dashed') +
  scale_color_manual(values = c(bioCluster.color, 'H-SC' = 'black'), name = 'Cluster') +
  xlim(-5.5, 5.5) + ylim(-5.5, 5.5) + #stat_cor(show.legend = F) +
  xlab('log2(D5H/D5SC)') + ylab('log2(D5Hplus/D5SC)')

gg.fcScatter

#/*==========================================================================#*/
#' ## Overlap D5/D5plus diff regions with SOM regions (Figure S4B)
#+ chunk_hbPlus_somOV, cache=T, echo=F, warning=F, message=F
#/*==========================================================================#*/
# Load SOM regions from *.Rds
som.regions <- readRDS(file.path(results.path, 'diffRegions.Rds'))
#diff.regions <- readRDS(file.path(regions.path, 'diffRegions.Rds'))

# Overlap SOM regions and differential regions
ov.som <- findOverlaps(som.regions, diff.regions)

# Plot gained/lost D5H-plus regions relativ to SOM clusters
gg.somOvBar <- data_frame(region_id = 1:length(diff.regions),
           cond = diff.regions$cluster) %>%
  left_join(., tbl_df(ov.som), by = c('region_id' = 'subjectHits')) %>%
  rowwise() %>%
  mutate(som_cluster = som.regions$bio_cluster[queryHits],
         som_cluster = if_else(is.na(som_cluster), 'D5-plus specific', som_cluster)) %>%
  group_by(cond, som_cluster) %>%
  summarise(n = n()) %>%
  mutate(freq = n/sum(n)) %>%
  ungroup() %>%
  mutate(cond = if_else(cond == 'D5H', 'D5H-plus lost', 'D5H-plus gained')) %>%
  ggplot(aes(x = reorder(som_cluster, freq), y = freq)) + # fill = som_cluster)) +
  geom_bar(stat = 'identity') +  #scale_fill_manual(values = color.pal) +
  facet_wrap(~cond) + xlab('') + ylab('Proportion of diff. regions (%)') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = 'none')

gg.somOvBar

#/*==========================================================================#*/
#' ## Motif enrichment with iCisTarget - day5 hindbrain plus (Figure S4C)
#+ chunk_hbPlus_iCisTarget, cache=T, echo=F, warning=F, message=F
#/*==========================================================================#*/
# Define path and find files
icistarget.path <- file.path(regions.path, 'icistarget')
icistarget.files <- listFiles(icistarget.path, pattern = 'tbl')

# Read statistics.tbl
icistarget.results <- data_frame(files = icistarget.files[grep(icistarget.files, pattern = 'stat.*tbl')]) %>%
  mutate(data = map(files, ~ read_tsv(.)),
         sample = str_extract(files, pattern = 'D5.*(H|plus|Som)')) %>%
  select(sample, data)

# Get featureIds for top 10 motifs
feature.ids <- icistarget.results %>%
  unnest()  %>%
  filter(FeatureDatabase == 'PWMs') %>%
  group_by(sample) %>%
  mutate(r = rank(-NES)) %>%
  filter(r <= 10) %>% ungroup %>%
  select(FeatureID) %>% unlist

# Plot NES from icistarget analysis
gg.icistarget <- icistarget.results %>%
  unnest() %>%
  filter(FeatureID %in% feature.ids) %>%
  select(sample, FeatureID, FeatureAnnotations, NES) %>%
  filter(sample %in% c('D5H', 'D5Hplus')) %>%
  filter(!is.na(FeatureAnnotations)) %>%
  ggplot(aes(x = sample, y = reorder(FeatureAnnotations, NES), fill = NES)) + geom_tile() +
  scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0) +
  xlab('') + ylab('TFs') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
gg.icistarget

