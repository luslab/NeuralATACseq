#/*==========================================================================#*/
#' ## Tn5 chrM bias (Figure S1 A)
#+ chunk_qc_chrM, cache=T, echo=T, warning=F, message=F
#/*==========================================================================#*/
# Create qc dir if not exists
if(!dir.exists(file.path(figures.path, 'qc'))) dir.create(file.path(figures.path, 'qc'))

# Get paths to selected samples
samplePath.txt <- list.files(atac.path, full.names = T, pattern = '2.*txt$')
sample.paths <- read_tsv(samplePath.txt, col_names = F) %>%
  mutate(X2 = sub('atacseq', 'atac', gsub('/home/camp', '/Users', X2)))

# Reformat sample names
sample.paths <- sample.paths %>%
  mutate(cond = if_else(str_count(X1, pattern = '_rep') > 0,
                        sub('-', '', gsub('_rep.*', '', X1)),
                        sub('-', '', gsub('-[1-3]', '', X1)))) %>%
  group_by(cond) %>%
  mutate(rep = 1:length(cond)) %>%
  mutate(sample = sprintf('%s_%s', cond, rep))

# Find chr distribution files for relevant samples
chrDist.files <- lapply(unlist(sample.paths$X2), function(p) {
  list.files(file.path(p, 'QC'), full.names = T, pattern = 'chrDist.txt')
}) %>% unlist(.)

# Read chr distribution from files
chr.dists <- data_frame(files = chrDist.files) %>%
  filter(!grepl(files, pattern = '(1|2)r')) %>%
  mutate(sample = gsub('_chrDist.*', '', basename(files)),
         data = map(files, ~ read_lines(.x))) %>%
  unnest() %>%
  mutate(chr = gsub('.*chr', 'chr', data),
         n = as.numeric(gsub(' c.*', '', data))) %>%
  select(-files, -data)

# Plot read chr distribution
gg.chrDist <- chr.dists %>%
  left_join(., sample.paths, by = c('sample' = 'X1')) %>%
  group_by(sample) %>%
  mutate(freq = (n/sum(n))*100) %>%
  filter(chr == 'chrM') %>%
  #filter(!grepl(cond, pattern = '(CDX|BRA)')) %>%
  mutate(cond = gsub('D5Hplus', 'D5H+', cond)) %>%
  mutate(cond = factor(cond, levels = names(color.pal))) %>%
  ggplot(aes(x = cond, y = freq, fill = cond)) + geom_bar(stat = 'identity') +
  scale_fill_manual(values = color.pal) +
  facet_grid(rep~.) + xlab('') + ylab('Proportion of chrM fragments (%)') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = 'none')

gg.chrDist

#savePlots(gg.plot = gg.chrDist, file.name = file.path(figures.path, 'qc/chrM_bias'),
#          base_width = 5, base_height = 5)

#/*==========================================================================#*/
#' ## Fragment length distribution (Figure S1 B)
#+ chunk_qc_fragLengthDist, cache=T, echo=T, warning=F, message=F
#/*==========================================================================#*/
# Find and load insertion *.RData
ins.rdata <- listFiles('/Users/steinhs/projects/atac_vicki/data/atacseq/*/*/rdata/', pattern = 'D5H-3')
load(ins.rdata)

# Plot fragment length distribution
gg.flength <- plotFragmentLengthDist(insertions.gr)
gg.flength

#/*==========================================================================#*/
#' ## TSS enrichment scores (Figure S1 C)
#+ chunk_qc_tssScores, cache=T, echo=T, warning=F, message=F
#/*==========================================================================#*/
# Plot ATAC-seq enrichment scores summary
enrichment.files <- lapply(unlist(sample.paths$X2), function(p) {
  list.files(file.path(p, 'QC'), full.names = T, pattern = 'Enrich')
}) %>% unlist(.)
#enrichment.files <- atacR:::listFiles(qc.path, pattern = 'Enrich')

# Read enrichment scores from files
enrichment.scores <- data_frame(files = enrichment.files) %>%
  filter(!grepl(files, pattern = '(1|2)r')) %>%
  mutate(sample = gsub('_atacEnr.*', '', basename(files)),
         data = map(files, read_tsv)) %>%
  select(-files) %>%
  unnest()

# Plot TSS enrichment scores
gg.enrichmentScores <- enrichment.scores %>%
  left_join(., sample.paths, by = c('sample' = 'X1')) %>%
  mutate(cond = gsub('D5Hplus', 'D5H+', cond)) %>%
  mutate(cond = factor(cond, levels = names(color.pal))) %>%
  ggplot(aes(x = cond, y = score, fill = cond)) + geom_bar(stat = 'identity') +
  facet_grid(rep~.) + xlab('') + ylab('TSS enrichment score') +
  scale_fill_manual(values = color.pal) +
  #geom_hline(yintercept = 5, col = 'Red', linetype = 'dashed') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = 'none')

gg.enrichmentScores
#savePlots(gg.plot = gg.enrichmentScores,
#          file.name = file.path(figures.path, 'qc/tssEnrichmentScores'),
#          base_width = 5, base_height = 5)

#/*==========================================================================#*/
#' ## TSS enrichment metaprofile (Figure S1 D)
#+ chunk_qc_tssEnrichMeta, cache=T, echo=T, warning=F, message=F
#/*==========================================================================#*/
# Load genes for given genome
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# Compute TSS enrichment score & plot TSS enrichment
gg.tssEnrichMeta <- plotTssEnrichment(insertions.gr, txdb)
gg.tssEnrichMeta

#/*==========================================================================#*/
#' ## Fraction of reads/fragments in peaks (Figure S1 E)
#+ chunk_qc_frip, cache=T, echo=T, warning=F, message=F
#/*==========================================================================#*/
# Find featureCount table file
summary.files <- list.files(file.path(results.path, 'atac_conRegions_all'),
                            full.names = T, pattern = 'summary')

# Compute fraction of read in peaks from featureCounts
frip.df <- read_tsv(summary.files[2], col_names = T) %>%
  filter(Status %in% c('Assigned', 'Unassigned_NoFeatures')) %>%
  gather(sample, n, 2:46) %>%
  mutate(sample = gsub('_rmbqr.*', '', basename(sample))) %>%
  spread(Status, n) %>%
  mutate(frip = Assigned/(Assigned + Unassigned_NoFeatures)) %>%
  filter(!grepl(sample, pattern = 'Sox2neg'))

# Plot FRiP
gg.frip <- frip.df %>%
  left_join(., sample.paths, by = c('sample' = 'X1')) %>%
  mutate(cond = gsub('D5Hplus', 'D5H+', cond)) %>%
  mutate(cond = factor(cond, levels = names(color.pal))) %>%
  ggplot(aes(x = cond, y = frip*100, fill = cond)) + geom_bar(stat = 'identity') +
  facet_grid(rep~.) + xlab('') + ylab('Fraction of fragments in consensus peaks (%)') +
  scale_fill_manual(values = color.pal) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = 'none')

gg.frip
#savePlots(gg.plot = gg.frip,
#          file.name = file.path(figures.path, 'qc/frip_barplot'),
#          base_width = 5, base_height = 5)

#/*==========================================================================#*/
#' ## CTCF footprinting benchmark (Figure S1 F)
#+ chunk_qc_footprintCTCF, cache=T, echo=T, warning=F, message=F
#/*==========================================================================#*/
# Define path to centipede dirs and find RData files with footprint motifs
centipede.path <- file.path(data.path, 'atacseq/vicki*/*/centipede')
pwm.id <- 'CTCF-JASPAR2014.MA0139.1'
footprint.files <- atacR::listFiles(file.path(centipede.path, pwm.id), pattern = 'RData$')
footprint.files <- footprint.files[grep(footprint.files, pattern = 'D0')]

# Load motifs inc. centipede probability
footprints <- lapply(footprint.files, function(footprint.file) {
  load(footprint.file)
  return(motifs)
})
names(footprints) <- gsub('_.*', '', basename(footprint.files))

# Keep only motifs which are in both replicates
footprints$`D0-3` <- subsetByOverlaps(footprints$`D0-3`, footprints$`D0-1`, type = 'equal')
footprints$`D0-1` <- subsetByOverlaps(footprints$`D0-1`, footprints$`D0-3`, type = 'equal')

### Plot footprints
# Classify footprints given the binding probability into bound/unbound
footprints <- lapply(footprints, function(footprint) {
  list('bound' = footprint[footprint$prob >= 0.99],
       'unbound' = footprint[footprint$prob < 0.99])
})

# Load insertions GRanges
ins.files <- atacR::listFiles(file.path(data.path, 'atacseq/vicki*/*/rdata'), pattern = 'RData$')
ins.files <- ins.files[grep(ins.files, pattern = 'D0')]
ins.files <- ins.files[c(3, 4)]

insertions.list <- lapply(ins.files, function(ins.file) {
  load(ins.file)
  return(insertions.gr)
})
names(insertions.list) <- gsub('_rmbqr.*', '', basename(ins.files))

# Compute insertion profile around CTCF motifs
ctcf.insProfiles <- lapply(names(footprints), function(fn) {
  print(fn)
  ctcf.insProfile <- lapply(names(footprints[[fn]]), function(type) {
    print(type)
    ctcf.insProfile <- computeInsertionProfile(insertions.gr = insertions.list[[fn]],
                                               regions = footprints[[fn]][[type]] + 150,
                                               strand.spec = T) %>%
      mutate(type = type, exp = fn)
    return(ctcf.insProfile)
  })
  return(bind_rows(ctcf.insProfile))
})

# Normalise insertion profile by subtracting Tn5 insertions from unbound motif (background)
ctcf.normProfile <- bind_rows(ctcf.insProfiles) %>%
  select(rel_pos, freq, type, exp) %>%
  spread(type, freq) %>%
  mutate(freq = bound - unbound,
         #freq = log2(bound/unbound),
         type = 'Normalised cuts') %>%
  select(rel_pos, freq, type, exp)

# Get motif width
motif.width <- width(footprints[[1]]$bound[1])

# Plot CTCF footprint Tn5 insertion profile
gg.ctcfProfile <- bind_rows(ctcf.insProfiles) %>%
  select(rel_pos, freq, type, exp) %>%
  bind_rows(., ctcf.normProfile) %>%
  #filter(grepl(type, pattern = 'Norm') &  exp == 'D0-3') %>%
  filter(!grepl(type, pattern = 'Norm') &  exp == 'D0-3') %>%
  ggplot(aes(x = rel_pos, y = freq, col = type)) + geom_line(size = 1) +
  xlab('Distance to motif mid (bp)') +
  ylab('Insertion frequency') +
  #ylab('Norm. insertion freq (bound - unbound)') +
  #ylab('log2(bound/unbound)') +
  #geom_hline(yintercept = 0, col = 'black', size = 0.5, linetype = 'dashed') +
  geom_vline(xintercept = c(-1,1)*floor(motif.width/2), col = 'red',
             linetype = 'dashed', size = 0.75) +
  theme(legend.position = 'none')

# Compute CTCF motif and higlight it in footprint
ctcf.seq <- getSeq(BSgenome.Mmusculus.UCSC.mm10, footprints[[1]]$bound)
ctcf.pwm <- consensusMatrix(ctcf.seq)[1:4,]

# Plot CTCF motif
ctcf.seqlogo <- ggseqlogo::ggseqlogo(ctcf.pwm) + ylim(0, 2) + theme_cowplot()

# Plot v-plot to add seqLogo on top of Tn5 insertion plot
gg.motifZoom <- data_frame(rel_pos = -125:125,
                           y = if_else(rel_pos < -9, rel_pos/min(rel_pos),
                                       if_else(rel_pos > 9, rel_pos/max(rel_pos), 0))) %>%
  ggplot(aes(x = rel_pos, y = y)) + geom_line(col = 'red', linetype = 'dashed', size = 0.75) +
  ylim(0.1, 1) + theme_transparent() +
  theme(plot.margin = margin(t = -0.75, b = -0.5, unit = 'cm'))

gg.footprint <- plot_grid(ctcf.seqlogo, gg.motifZoom, gg.ctcfProfile, nrow = 3,
                          align = 'v', rel_heights = c(2, 0.1, 4.5))
gg.footprint
#savePlots(gg.plot = gg.footprint,
#          file.name = file.path(figures.path, 'D0-3_CTCF-footprint'),
#          base_width = 4.75)
