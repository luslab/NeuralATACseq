#/*==========================================================================#*/
#' ## Plot Mafb loci - Cdx2 binding example (Figure S5A'')
#+ chunk_regionPlots_MafbCdx2, cache=T, echo=F, warning=F, message=F
#/*==========================================================================#*/
# Define gene of interest and region size
gene.ofi <- 'Mafb'

# Get ENSEMBL gene id for symbol
gene.id <- gencode.metadata%>%
  filter(gene_name %in% gene.ofi) %>%
  arrange(gene_name) %>% unique(.) %>%
  mutate(gene_id = gsub('\\.[0-9]*', '', gene_id)) %>%
  select(gene_id) %>% unlist() %>% unique()

### Plot ATAC-signal and gene expression for given genomic loci
bw.path <- '/Users/steinhs/projects/atac_vicki/results/atac_conregions_KO/bw/'
bw.files <- atacR::listFiles(bw.path, pattern = 'fpm.bw')
bw.files <- bw.files[grep(bw.files, pattern = '(D3|D4|D5)')]
bw.files <- bw.files[-grep(bw.files, pattern = '(D4-A|D5-A|CDX|BRA)')]

# Get promoter region for gene
promoter.region <- GRanges(seqnames = 'chr2', IRanges(160356793, 160538238)) # MAFB

promoter.chr <- seqnames(promoter.region) %>% as.character()
region.limits <- c(start(promoter.region), end(promoter.region))

# Parse FPM for bigWig files for given promoter region
pregion.signal <- lapply(bw.files, function(bw.file) {
  pregion.signal <- rtracklayer::import.bw(bw.file, which = promoter.region)
  pregion.signal$sample <- gsub('_.*', '', basename(bw.file))
  return(as.data.frame(pregion.signal))
})
pregion.signal <- bind_rows(pregion.signal)

# Plot gene models for promoter region
promoter.genes <- subsetByOverlaps(genes(gencode.txdb), promoter.region)
promoter.genes <- promoter.genes[grep(promoter.genes$gene_id, pattern = gene.id)]
tx.id <- gencode.metadata %>%
  filter(gene_id %in% unlist(promoter.genes$gene_id)) %>%
  group_by(gene_name) %>%
  filter(exon_number == max(exon_number, na.rm = T)) %>%
  ungroup() %>%
  unique(.) %>% select(transcript_id) %>% unlist(.)

exons <- exonsBy(gencode.txdb, by = "tx", use.names = TRUE)[tx.id]
cds <- cdsBy(gencode.txdb, by = "tx", use.names = TRUE)
i.cds <- lapply(tx.id, function(tid) which(names(cds) == tid)) %>% unlist(.)
cds <- cds[i.cds]

# Plot gene model
gg.gene <- plotTranscripts(exons, cds, gencode.metadata, rescale_introns = F,
                           region_coords = region.limits) +
  xlab('') + ylab('GENE') + theme_bw() + theme_cowplot() +
  theme(plot.margin = margin(r = -0.5, l = 1, t = 0, b = 0, unit = 'cm'),
        legend.position = 'none',
        strip.background = element_rect(fill = 'white'),
        axis.title.y = element_text(color = 'white'),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())

# Plot ATAC-seq signal tracks
signal.max <- ceiling(pregion.signal %>% summarise(max(score))) %>% unlist()
n.by <- signal.max/2
gg.tracks <- pregion.signal %>%
  tbl_df() %>%
  mutate(sample = sub('-', '',  sample)) %>%
  ggplot(aes(xmin = start, xmax = end, ymin = 0, ymax = score, fill = sample)) +
  geom_rect() + scale_fill_manual(values = color.pal) +
  scale_y_continuous(breaks = seq(0, signal.max, by = n.by), limits = c(0, signal.max)) +
  scale_x_continuous(limits = region.limits, expand = c(0,0)) +
  facet_grid(sample~., scales = 'free_y') + xlab('Genomic position (bp)') + ylab('FPM') +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(colour = 'white'),
        axis.title.x = element_text(colour = 'white'),
        plot.margin = margin(r = -1, l = 1, t = -0.5, b = 0, unit = 'cm'))

### Plot Cdx2 logFC signal
cdx2.bws <- listFiles(file.path(chipseq.path, '*/*/*/bw'), pattern = 'Cdx.*sub')
cdx2.bws <- cdx2.bws[c(1,3)]

# Parse FPM for bigWig files for given promoter region
cdx2.signal <- lapply(cdx2.bws, function(bw) {
  pregion.signal <- rtracklayer::import.bw(bw, which = promoter.region)
  pregion.signal$sample <- gsub('_rep.*', '', basename(bw))
  return(as.data.frame(pregion.signal))
})
cdx2.signal <- bind_rows(cdx2.signal)

#
cdx2.max <- ceiling(cdx2.signal %>% summarise(max(score))) %>% unlist()
cdx2.by <- cdx2.max/2

# Cdx2 ChIP-seq tracks plot
gg.cdx2Tracks <- cdx2.signal %>%
  tbl_df() %>%
  mutate(sample = sub('-', '',  sample)) %>%
  filter(score > 0) %>%
  ggplot(aes(xmin = start, xmax = end, ymin = 0, ymax = score, fill = sample)) +
  geom_rect() +  scale_fill_brewer(palette = 'Paired', name = '') +
  scale_y_continuous(breaks = seq(0, cdx2.max, by = cdx2.by), limits = c(0, cdx2.max)) +
  scale_x_continuous(limits = region.limits, expand = c(0,0)) +
  facet_grid(sample~., scales = 'free_x') + xlab('Genomic position (bp)') + ylab('log2(FC)') +
  theme(legend.position =  c(0.25, 0.6),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(colour = 'white'),
        axis.title.x = element_text(colour = 'white'),
        plot.margin = margin(r = -0.5, l = 1, t = -1, b = 0, unit = 'cm'))

### Plot expression data as barplot
# Get TPM for gene of interest
tpm.ofi <- log2(tpms[gene.id,] + 1)

# Gene expression as barplot
gg.tpmBar <- tpm.ofi %>%
  tbl_df() %>%
  mutate(sample = names(tpm.ofi)) %>%
  separate(sample, c('cond', 'rep'), sep = '_') %>%
  mutate(cond = factor(cond, levels = sort(unique(cond), decreasing = F))) %>%
  filter(cond %in% names(color.pal)) %>%
  group_by(cond) %>%
  summarise(mean_tpm = mean(value),
            sd = sd(value),
            std_error = atacR:::stdError(value)) %>%
  filter(grepl(cond, pattern = '(D3A|NMP|SC|H)')) %>%
  filter(!grepl(cond, pattern = '2.5')) %>%
  ggplot(aes(x = cond, y = mean_tpm, fill = cond)) +
  geom_bar(stat = 'identity', width = 0.5) +
  geom_errorbar(aes(ymin = mean_tpm - sd, ymax = mean_tpm + sd), col = 'darkgrey', size = 1, width = 0.15) +
  xlab('') + ylab('log2(TPM + 1)') + scale_fill_manual(values = color.pal) +
  coord_flip() + facet_grid(cond~., scales = 'free_y') +
  theme(legend.position = 'none',
        axis.text.y = element_blank(),
        axis.line.y = element_line(linetype = 'solid'),
        axis.ticks.y = element_blank(),
        plot.margin = margin(l = -0.5, unit = 'cm'))

# Get Cdx2 track legend
legend.cdx2 <- get_legend(gg.cdx2Tracks)
gg.cdx2Tracks <- gg.cdx2Tracks + theme(legend.position = 'none')

# Combine genemodel, tracks, TPM barplot and interactions
gg.loci <- plot_grid(gg.gene, gg.tracks, gg.cdx2Tracks, nrow = 3, rel_heights = c(1.5, 6, 2), align = 'v', axis = 'lr')
gg.exp <- plot_grid(NULL, gg.tpmBar, legend.cdx2, nrow = 3, rel_heights = c(1.5, 6, 2))
gg.genomicLoci <- plot_grid(gg.loci, gg.exp, ncol = 2, rel_widths = c(3, 1))

gg.genomicLoci
