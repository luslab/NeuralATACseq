#/*==========================================================================#*/
#' ## Cluster variable regions using self-organising maps (SOMs) (Figure 2 A)
#+ chunk_deseq2WT_som, cache=T, echo=F, warning=F, message=F
#/*==========================================================================#*/
### SOM Pre-processing
#  Average replicates and build count matrix
count.matrix <- log2(norm.counts + 1)%>%
  tbl_df() %>%
  mutate(id = 1:nrow(norm.counts)) %>%
  gather(sample, counts, 1:ncol(norm.counts)) %>%
  mutate(sample = gsub('_[1-3]', '', sample)) %>%
  group_by(id, sample) %>%
  summarise(counts = mean(counts)) %>%
  spread(sample, counts) %>%
  ungroup() %>%
  select(-id)

# Re-scale count matrix that row mean ~ 0 and var = 1
zscore.matrix <- t(apply(count.matrix, 1, function(r) (r - mean(r))/sd(r)))

# Sort columns (according to lineages not timepoints)
i.sort <- c(1:4, 6, 9, 7, 10, 5, 8, 11)
zscore.matrix <- zscore.matrix[,i.sort]
sorted.samples <- colnames(zscore.matrix)

### SOM cluserting
# Self-organising map clustering with 'som' package
som.x <- 5
som.y <- 5
som.results <- som::som(zscore.matrix, som.x, som.y, neigh = 'gaussian', topol = 'hex')

# Get cluster definitions
som.cluster <- paste(som.results$visual$x, som.results$visual$y, sep = "_")

# Compute number of clusters
n.cluster <- data_frame(x = som.results$visual$x, y = som.results$visual$y) %>%
  group_by(x, y) %>%
  summarise(n = n()) %>%
  mutate(label = sprintf('n = %s', n))

# Add cluster to matrix and plot
clustered.matrix <- zscore.matrix %>%
  tbl_df() %>%
  mutate(x = som.results$visual$x, y = som.results$visual$y)

### Cluster SOM into bigger meta-clusters
# Compute distance measure: pearson corr based on z-score medians
#som.codes <- som.results$code %>% tbl_df()
som.codes <- zscore.matrix %>%
  tbl_df() %>%
  mutate(cluster = som.cluster) %>%
  gather('sample', 'zscore', 1:ncol(zscore.matrix)) %>%
  group_by(cluster, sample) %>%
  summarise(zscore = median(zscore)) %>%
  spread(sample, zscore)

# Plot SOM grid as violin plots
gg.som <- clustered.matrix %>%
  gather(sample, norm_count, 1:11) %>%
  mutate(sample = factor(sample, levels = sorted.samples)) %>%
  ggplot(aes(x = sample, y = norm_count)) + geom_violin() +
  facet_grid(y~x) + xlab('') + ylab('Z-Score') +
  geom_text_repel(data = n.cluster, aes(x = 0.75, y = 2.5, label = label), inherit.aes = F) +
  geom_hline(yintercept = 0, col = 'red', linetype = 'dashed') +
  #scale_fill_manual(values = color.pal) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = 'none')
gg.som

#save_plot(filename = som.pdf, plot = gg.som, base_height = 7, base_width = 9)
# Save genomic loci plot
#som.file <- file.path(figures.path, sprintf('conRegions_deseq2ESCdiff_som-x%s-y%s_uncolored', som.x, som.y))
#savePlots(gg.som, file.name = som.file, base_width = 9, base_height = 8)

# Get differential regions and assign SOM location
diff.regions <- regions[i.diff,]
diff.regions$cluster <- som.cluster

# Write regions to file
#som.path <- file.path(deseq2.path)
#lapply(unique(som.cluster), function(c) {
#  print(c)
#  cluster.diffRegions <- diff.regions[diff.regions$cluster == c]
#  cluster.diffRegions <- as.data.frame(cluster.diffRegions)[,1:3]
#  diffRegions.bed <- file.path(som.path, sprintf('diffRegions_%s.bed', c))
#  write_tsv(cluster.diffRegions, diffRegions.bed, col_names = F)
#})

### Manually define SOM clusters
som.bioCluster <- list('Hindbrain' = '0_0',
                       'Neural' = c('1_0', '2_0', '3_0', '1_1', '2_1', '3_1'),
                       'Anterior' = c('4_0', '4_1'),
                       'H-SC' = c('0_1', '0_2', '0_3', '1_2', '1_3'),
                       'Spinal cord' = c('0_4', '1_4'),
                       'NMP-SC' = c('2_4'),
                       'NMP' = c('3_4'),
                       'Epi' = c('4_4'),
                       'A-H' = c('4_2', '4_3'),
                       'Undefined' = c('2_2', '2_3', '3_2', '3_3'))

i.som <- match(diff.regions$cluster, unlist(som.bioCluster))
diff.regions$bio_cluster <- names(unlist(som.bioCluster))[i.som] %>% gsub('[0-9]', '', .)

# Write biological clusters to file
#lapply(names(som.bioCluster), function(c) {
#  print(c)
#  cluster.diffRegions <- diff.regions[diff.regions$bio_cluster == c]
#  cluster.diffRegions <- as.data.frame(cluster.diffRegions)[,1:3]
#  diffRegions.bed <- file.path(som.path, sprintf('diffRegions_%s.bed', c))
#  write_tsv(cluster.diffRegions, diffRegions.bed, col_names = F)
#})

# Add normalised count data to diff. regions and save GRange as RData file
mcols(diff.regions) <- cbind(mcols(diff.regions), count.matrix)

#saveRDS(diff.regions, file = file.path(results.path, 'diffRegions.Rds'))
