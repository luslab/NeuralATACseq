#/*==========================================================================#*/
#' ## Plot Phox2b loci - Hindbrain example (Figure 2C)
#+ chunk_regionPlots_Phox2b, cache=T, echo=F, warning=F, message=F
#/*==========================================================================#*/
# Define gene of interest and region size
gene.ofi <- 'Phox2b'

# Get ENSEMBL gene id for symbol
gene.id <- gencode.metadata%>%
  filter(gene_name %in% gene.ofi) %>%
  arrange(gene_name) %>% unique(.) %>%
  mutate(gene_id = gsub('\\.[0-9]*', '', gene_id)) %>%
  select(gene_id) %>% unlist() %>% unique()

### Plot ATAC-signal and gene expression for given genomic loci
bw.path <- '/Users/steinhs/projects/atac_vicki/results/atac_conregions/bw/'
bw.files <- atacR::listFiles(bw.path, pattern = 'fpm.bw')
bw.files <- bw.files[grep(bw.files, pattern = '(D0|D5)')]

# Get promoter region for gene
promoter.region <- GRanges(seqnames = 'chr5',
                           ranges = IRanges(67079562, 67221323))
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
if(n.by < 1) n.by <- 1
gg.tracks <- pregion.signal %>%
  tbl_df() %>%
  mutate(sample = gsub('-', '',  sample)) %>%
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
  filter(grepl(cond, pattern = '(D0|D5)')) %>%
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

### Interaction data
# Define set of regions of interest
regions.promoter <- subsetByOverlaps(diff.regions, promoter.region)

# Find interactions overlapping with the promoter region
promoter.inter <- subsetByOverlaps(interactions, promoter.region, type = 'within')
promoter.inter <- subsetInteractionsByRegion(promoter.inter, regions.promoter)
promoter.inter <- subsetByOverlaps(promoter.inter, regions.promoter)

# Remove to large interacting regions
i.toLarge <- which(c(width(anchorOne(promoter.inter)) > 10^4 |
  width(anchorTwo(promoter.inter)) > 10^4))
if(length(i.toLarge) > 0) promoter.inter <- promoter.inter[-i.toLarge,]

# Filter for trans-interactions
promoterInter.df <- promoter.inter %>% tbl_df() %>%
  filter(seqnames1 == seqnames2) %>%
  filter((Cell.Tissue %in% c('Cortex', 'Cerebellum', 'ESC', 'NPC', 'NSC', 'Whole brain'))) %>%
  #filter(!(Cell.Tissue %in% c('Th1'))) %>%
  mutate(ymin = 0:(length(seqnames1) - 1) + 0.5, ymax = ymin)

inter.max <- max(promoterInter.df$ymax)

regionInter.df <- promoterInter.df[,c(1:3, 14)] %>%
  mutate(type = 'A') %>%
  rename(seqnames = seqnames1, start = start1, end = end1, sample = Cell.Tissue)

regionInter.df <- promoterInter.df[,c(6:8, 14)] %>%
  mutate(type = 'B') %>%
  rename(seqnames = seqnames2, start = start2, end = end2,  sample = Cell.Tissue) %>%
  bind_rows(., regionInter.df) %>%
  group_by(type) %>%
  mutate(ymin = 0:(length(seqnames) - 1),
         ymax = ymin + 1)

gg.inter <- regionInter.df %>%
  ggplot(aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax, fill = sample)) +
  geom_segment(data = promoterInter.df, aes(x = start1, xend = end2, y = ymax, yend = ymax), inherit.aes = F) +
  geom_rect(col = 'black') + ylab('Interactions') +
  scale_fill_brewer(palette = 'Dark2', name = '') +
  geom_rect(data = tbl_df(regions.promoter),
            aes(xmin = start, xmax = end, ymin = inter.max + 2, ymax = inter.max + 4),
            fill = 'black', inherit.aes = F) +
  xlab(sprintf('%s', promoter.chr)) +
  scale_x_continuous(limits = region.limits, expand = c(0,0)) +
  theme(legend.position = c(0.25, 0.6),
        axis.ticks.y = element_blank(),
        plot.margin = margin(r = -0.5, l = 1, t = -1, b = 0, unit = 'cm'),
        axis.text.y = element_blank())

# Get interaction plot legend
legend.inter <- get_legend(gg.inter)
gg.inter <- gg.inter + theme(legend.position = 'none')

# Combine genemodel, tracks, TPM barplot and interactions
gg.loci <- plot_grid(gg.gene, gg.tracks, gg.inter, nrow = 3, rel_heights = c(1.5, 8, 2), align = 'v', axis = 'lr')
gg.exp <- plot_grid(NULL, gg.tpmBar, legend.inter, nrow = 3, rel_heights = c(1.5, 8, 2))
gg.genomicLoci <- plot_grid(gg.loci, gg.exp, ncol = 2, rel_widths = c(3, 1))

gg.genomicLoci
