#/*==========================================================================#*/
#' ## Perform Cdx2 footprint analysis (CENTIPEDE) (Figure 4J)
#+ chunk_cdx2_footprinting, cache=T, echo=F, warning=F, message=F
#/*==========================================================================#*/
# Define path to insertion RData files
rdata.path <- file.path(data.path, 'atacseq/*/*/rdata')
ins.rdata <- listFiles(rdata.path, pattern = 'RData')

# Subset for samples of interest
ins.rdata <- ins.rdata[grep(ins.rdata, pattern = '(NMP|D4SC|D5SC)')]
ins.rdata <- ins.rdata[grep(ins.rdata, pattern = '(CDX-2|NMP-3)')]

# Load D3 NMP WT best replicate
load(ins.rdata[1])

# Find Cdx2 motif matches in and around differential regions
opts <- list('taxon' = 'vertebrates', name = 'CDX2')
cdx2.pwm <- getMatrixSet(JASPAR2016, opts = opts)
footprints <- matchMotifs(pwms = cdx2.pwm,
                          subject = getSummit(diff.regions) + 5000,
                          genome = BSgenome.Mmusculus.UCSC.mm10,
                          out = 'position', p.cutoff = 10^-4)
footprints <- footprints[[1]]

# Compute Tn5 Insertion matrix
centipede.matrix <- computeCentipedeMatrix(insertions.gr[insertions.gr$fraglen <= 100], footprints, motif.ext = 150)
score.matrix <- cbind(rep(1, length(footprints)), footprints$score)

# Run CENTIPEDE
centipede.fit <- CENTIPEDE::fitCentipede(Xlist = list('atac' = centipede.matrix),
                                         Y = score.matrix)

# Plot CENTIPE profile
CENTIPEDE::plotProfile(centipede.fit$LambdaParList[[1]])

# Prepare bound/unbound Cdx2 motifs from footprinting
footprints$prob <- centipede.fit$PostPr[,1]
footprints <- subsetByOverlaps(footprints, diff.regions)
cdx2.footprints <- list('bound' = footprints[footprints$prob >= 0.9],
                        'unbound' = footprints[footprints$prob < 0.9])

# Load insertion GRange from RData file and compute insertion profile around motif
cdx2.insProfiles <- lapply(ins.rdata, function(ins) {
  print(ins)
  load(ins)
  cdx2.insProfile <- lapply(names(cdx2.footprints), function(fn) {
    print(fn)
    cdx2.insProfile <- computeInsertionProfile(insertions.gr = insertions.gr[insertions.gr$fraglen <= 100],
                                               regions = cdx2.footprints[[fn]] + 175,
                                               strand.spec =  T) %>%
      mutate(type = fn, sample = gsub('_.*', '', basename(ins)))
    return(cdx2.insProfile)
  })
  cdx2.insProfile <- bind_rows(cdx2.insProfile)
  return(cdx2.insProfile)
})

gg.cdx2Footprint <- bind_rows(cdx2.insProfiles) %>%
  select(sample, rel_pos, freq, type) %>%
  #spread(type, freq) %>%
  #group_by(sample) %>%
  #mutate(freq = log2(bound/unbound)) %>%
  ggplot(aes(x = rel_pos, y = freq, col = type)) + geom_line() + facet_wrap(~sample) +
  geom_vline(xintercept = c(-6,6), col = 'red', linetype = 'dashed') +
  #geom_hline(yintercept = 0, linetype = 'dashed') +
  scale_color_manual(values = c('red', 'blue'), name = '') +
  xlab('Distance realtive to motif (bp)') +
  ylab('Insertion frequency') + xlim(-150, 150) + theme(legend.position = 'bottom')

gg.cdx2Footprint

#savePlots(gg.cdx2Footprint, file.name = file.path(figures.path, 'Cdx2_footprint'),
#          base_height = 4, base_width = 5.5)

#/*==========================================================================#*/
#' ## Cdx2 ChIP-seq peak overlap with SOM clustering (Figure S5A)
#+ chunk_cdx2_peakOv, cache=T, echo=F, warning=F, message=F
#/*==========================================================================#*/
# Define path to Cdx2 ChIP-seq peaks
chipseq.path <- file.path(data.path, 'chipseq')
cdx2.beds <- atacR::listFiles(file.path(chipseq.path, '*/*/*/macs2'), pattern = '.*narrow')
cdx2.beds <- cdx2.beds[grep(cdx2.beds, pattern = '(Cdx)')] #(Cdx|Sox2)')]
cdx2.beds <- cdx2.beds[-grep(cdx2.beds, pattern = 'noFGF')]

# Read peak from files
cdx2.peaks <- lapply(cdx2.beds, importNarrowPeak)
cdx2.peaks <- GRangesList(cdx2.peaks)
names(cdx2.peaks) <- gsub('_peaks.*', '', basename(cdx2.beds))

# Combine replicates by keeping only overlaping peaks
cdx2.samples <- gsub('_rep.*', '', names(cdx2.peaks)) %>% unique
cdx2.peaks <- lapply(cdx2.samples, function(cdx2.sample) {
  i.cdx2 <- grep(names(cdx2.peaks), pattern = cdx2.sample)
  if(length(i.cdx2) > 1) {
    subsetByOverlaps(cdx2.peaks[[i.cdx2[1]]], cdx2.peaks[[i.cdx2[2]]])
  } else {
    cdx2.peaks[[i.cdx2]]
  }
})
names(cdx2.peaks) <- cdx2.samples
cdx2.peaks <- GRangesList(cdx2.peaks)

# Compute SOM cluster size
n.som <- diff.regions$cluster %>%
  table() %>% tbl_df() %>%
  rename('som_cluster' = '.', 'n_regions' = 'n')

# Overlap SOM regions with Cdx2 ChIP-seq peaks
cdx2.ov <- findOverlaps(diff.regions, cdx2.peaks)

# Plot Cdx2 peak overlap over SOM
gg.somCdx <- data_frame(peak_id = queryHits(cdx2.ov),
                        som_cluster = diff.regions$cluster[queryHits(cdx2.ov)],
                        bio_cluster = diff.regions$bio_cluster[queryHits(cdx2.ov)],
                        chip_ov = names(cdx2.peaks)[subjectHits(cdx2.ov)]) %>%
  group_by(som_cluster, chip_ov) %>%
  summarise(n = n()) %>%
  left_join(., n.som, by = 'som_cluster') %>%
  mutate(freq = (n/n_regions)*100) %>%
  separate(som_cluster, c('x', 'y'), sep = '_', remove = F) %>%
  mutate(chip_ov = gsub('_', ' ', chip_ov)) %>%
  ggplot(aes(x = chip_ov, y = freq, fill = chip_ov)) +
  geom_bar(stat = 'identity') + facet_grid(y~x) +
  scale_fill_brewer(palette = 'Dark2', name = '') +
  guides(fill = guide_legend(nrow = 2)) +
  ylab('Proportion of regions overlapping peak (%)') +
  xlab('') +  theme(legend.position = 'bottom',
                    axis.text.x = element_blank(),
                    axis.ticks.x = element_blank())
gg.somCdx

#savePlots(gg.somCdx, file.name = file.path(figures.path, 'CdxOverlapSom'),
#          base_width = 5, base_height = 5)

# Write Cdx2 peaks to files --> plot region heatmap
#ov.cdxMnp <- findOverlaps(cdx2.peaks[['MNP_Cdx2-FGF_rep1']], diff.regions)
#cdxMnp.peaks <- cdx2.peaks[['MNP_Cdx2-FGF_rep1']][queryHits(ov.cdxMnp)]
#cdxMnp.peaks$cluster <- diff.regions$bio_cluster[subjectHits(ov.cdxMnp)]
#lapply(c('NMP', 'NMP-SC', 'Spinal cord', 'H-SC', 'Hindbrain'), function(bc) {
#  print(bc)
#  cmp <- cdxMnp.peaks[cdxMnp.peaks$cluster == bc] %>%
#    as.data.frame
#  cmp.bed <- file.path('/Users/steinhs/projects/atac_vicki/results/cdx2/regions', sprintf('%s_MnpCdx2.bed', bc))
#  write_tsv(cmp[,1:3], path = cmp.bed, col_names = F)
#})

#/*==========================================================================#*/
#' ## GREAT analysis of Cdx2 peaks overlapping with SOM (Figure S5C,D,E,F)
#+ chunk_cdx2_great, cache=T, echo=F, warning=F, message=F
#/*==========================================================================#*/
# Set GREAT path and find GO enrichment files
great.path <- file.path(results.path, 'cdx2/great')
go.files <- list.files(great.path, full.names = T, pattern = 'GO')

# Parse GREAT GO results
great.results <- data_frame(files = go.files) %>%
  mutate(bio_cluster = gsub('-Cdx2_.*', '', basename(files)),
         bio_cluster = gsub('Spinalcord', 'Spinal cord', bio_cluster),
         data = map(files, ~ read_tsv(., skip = 1))) %>%
  select(-files) %>% unnest() %>%
  rename(padj = `Binom FDR Q-Val`, term = `# Term Name`)

# Plot GREAT GO enrichment as barplot and save as pdf
lapply(unique(great.results$bio_cluster), function(bc) {
  print(bc)
  gg.great <- great.results %>%
    filter(bio_cluster == bc) %>%
    mutate(padj = -log10(padj)) %>%
    ggplot(aes(x  = reorder(term, padj), y = padj, fill = bio_cluster)) +
    geom_bar(stat = 'identity') +
    scale_fill_manual(values = c(bioCluster.color, 'NMP-SC' = '#F79992',
                                 'NMP' = '#F79992', 'H-SC' = 'black')) +
    coord_flip() + ylab('-log10(adj. pvalue)') + xlab('') +
    #geom_hline(yintercept = -log10(0.01), col = 'red', linetype = 'dashed') +
    theme_cowplot(font_size = 10) + theme(legend.position = 'none')
  gg.great
  #savePlots(gg.plot = gg.great,
  #          file.name = file.path(figures.path, sprintf('%s-Cdx2_greatEnrichment', bc)),
  #          base_width = 6, base_height = 4)
})

#/*==========================================================================#*/
#' ## Cdx2 peaks overlapping with SOM - iCisTarget motif enrichment (Figure S5H)
#+ chunk_cdx2_icistarget, cache=T, echo=F, warning=F, message=F
#/*==========================================================================#*/
# Define path to icistarget results
icistarget.path <- '/Users/steinhs/projects/atac_vicki/results/cdx2/icistarget/'
icis.files <- listFiles(icistarget.path, pattern = 'tbl$')

# Read icis target results
icis.results <- lapply(icis.files[grep(icis.files, pattern = 'stat')], function(cf) {
  read_tsv(cf)
})

# Get top ten motifs from each conditions
icis.features <- lapply(icis.results, function(icr) {
  icr %>%
  filter(FeatureDatabase == 'PWMs') %>%
  arrange(desc(NES)) %>%
  mutate(r = 1:length(NES)) %>%
  filter(r <= 20) %>%
  select(FeatureID, FeatureAnnotations)
})
icis.features <- bind_rows(icis.features) %>% unique() %>%
  filter(!is.na(FeatureAnnotations))

# Tidy iCisTarget output for plotting
icis.matrix <- bind_rows(icis.results) %>%
  filter(FeatureID %in% icis.features$FeatureID) %>%
  select( `#GeneSignatureID`, FeatureID, FeatureAnnotations, NES) %>%
  group_by(FeatureID) %>%
  mutate(max_NES = max(NES),
         cv = sd(NES)/mean(NES),
         entropy  = shannon.entropy(NES)) %>%
  spread(`#GeneSignatureID`, NES) %>%
  arrange(desc(cv)) %>%
  filter(!is.na(FeatureAnnotations))

# Plot heatmap with conditons over top 20 motifs
gg.icisHeatmap <- icis.matrix %>%
  filter(str_count(FeatureAnnotations, ',') < 2) %>%
  group_by(FeatureAnnotations) %>%
  filter(max_NES == max(max_NES)) %>%
  arrange(desc(Hindbrain)) %>%
  gather(sample, NES, 6:10) %>%
  mutate(sample = factor(sample, levels = c('Hindbrain', 'H-SC', 'SC', 'NMP-SC', 'NMP'))) %>%
  ggplot(aes(x = sample, y = FeatureAnnotations, fill = NES)) + geom_tile() +
  scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0) +
  xlab('') + ylab('') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
gg.icisHeatmap

#savePlots(gg.plot = gg.icisHeatmap,
#          file.name = file.path(figures.path, 'Cdx2Regions_iCisTargetHeatmap'),
#          base_height = 6, base_width = 5)
