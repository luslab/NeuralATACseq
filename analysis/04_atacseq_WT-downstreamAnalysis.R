#/*==========================================================================#*/
#' ## Genomic annotation using ChIPseeker (Figure S2 F,H,I,G)
#+ chunk_deseq2WT_genomeAnno, cache=T, echo=F, warning=F, message=F
#/*==========================================================================#*/
if(exists('diff.regions')) {
  diff.regions <- readRDS(file.path(results.path, 'diffRegions.Rds'))
}
# Define colors and bioclusters of interest
bio.cluster  <- c('Anterior', 'Hindbrain', 'Spinal cord', 'Neural')

# Annotate SOM clusters -ChIPseeker and Gencode
som.anno <- lapply(bio.cluster, function(bc){
  som.anno <- getSummit(diff.regions[diff.regions$bio_cluster == bc])
  som.anno <- annotatePeak(som.anno, TxDb = gencode.txdb, annoDb = 'org.Mm.eg.db')
  som.anno <- som.anno@annoStat %>%
    tbl_df() %>%
    mutate(cluster = bc)
  return(som.anno)
})
som.anno <- bind_rows(som.anno)

# Plot genomic annotation for each bio cluster
lapply(bio.cluster, function(bc) {
  print(bc)
  # Plot genomic annotation for a given cluster
  gg.anno <- som.anno %>%
    dplyr::filter(cluster == bc) %>%
    ggplot(aes(x = Feature, y = Frequency, fill = cluster)) + geom_bar(stat = 'identity') +
    scale_fill_manual(values = bioCluster.color, name = 'Cluster') +
    xlab('') + ylab('Proportion of regions (%)') +
    theme(legend.position = 'none',
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

  # Save genomic annotation on SOM plot
  #anno.file <- file.path(figures.path, sprintf('genomicAnno_%s', bc))
  #savePlots(gg.anno, file.name = anno.file, base_height = 4.5)
})

#/*==========================================================================#*/
#' ##  Metaprofile of NPC H3K27ac over biocluster regions (Figure 1 G)
#+ chunk_deseq2WT_k27acProfile, cache=T, echo=F, warning=F, message=F
#/*==========================================================================#*/
# Define paths
k27ac.path <- file.path(data.path, 'chipseq/nishi2012_GSE42132')
k27ac.bws <- list.files(k27ac.path, recursive = T, full.names = T, pattern = 'K27ac.*log2')

# Comput consensus regions around summit
profile.summits <- lapply(bio.cluster, function(bc) {
  getSummit(diff.regions[diff.regions$bio_cluster == bc]) + 2700
})
names(profile.summits) <- bio.cluster
profile.summits <- GRangesList(profile.summits)
profile.summits$ESC <- getSummit(regions[i.diffEsc,]) + 2700

# Compute H3K27ac metaprofiles
k27ac.profiles <- lapply(k27ac.bws[2], function(k27ac.bw) {
  print(k27ac.bw)
  k27ac.profile <- lapply(names(profile.summits), function(rn) {
    print(rn)
    k27ac.profile <- computeProfile(profile.summits[[rn]], k27ac.bw, region.names = rn)
    return(k27ac.profile)
  })
  k27ac.profile <- bind_rows(k27ac.profile)
  k27ac.profile$sample <- gsub('_rep.*', '', basename(k27ac.bw))
  return(k27ac.profile)
})
k27ac.profiles <- bind_rows(k27ac.profiles)

# Plot H3K27ac metaprofiles
gg.k27acProfiles <- k27ac.profiles %>%
  mutate(sample = gsub('_', ' ', sample)) %>%
  ggplot(aes(x = rel_pos/10^3, y = mean_score, col = name, fill = name)) +
  geom_line() + scale_color_manual(values = bioCluster.color, name = 'Cluster') +
  scale_fill_manual(values = bioCluster.color, name = 'Cluster') +
  geom_ribbon(aes(ymin = mean_score - std_error, ymax = mean_score + std_error), alpha = 0.25) +
  xlab('Distance to summit (kbp)') + geom_vline(xintercept = 0, linetype = 'dashed') +
  xlim(-2.5, 2.5) + ylab('Mean log2(FC)') + facet_wrap(~sample, scales = 'free_x') +
  theme(legend.position = 'bottom')

gg.k27acProfiles
#savePlots(gg.k27acProfiles,
#          file.name = file.path(figures.path, 'NEB-H3K27ac_metaprofiles'),
#          base_height = 4.5, base_width = 5)

#/*==========================================================================#*/
#' ## ChIP-seq enrichment using LOLA (Figure 2 H,I,J & Figure 4 I)
#+ chunk_deseq2WT_lola, cache=T, echo=F, warning=F, message=F
#/*==========================================================================#*/
# Load LOLA database
lolaDb.path <- '/Users/steinhs/public_data/LOLACoreCaches/scratch/ns5bc/resources/regions/LOLACore/mm10/'
lola.db <- loadRegionDB(lolaDb.path)

# LOLA figure path
lola.path <- file.path(figures.path, 'lola')
if(!dir.exists(lola.path)) dir.create(lola.path)

### Run LOLA enrichment
# For all biological clusters
lola.results <- lapply(unique(diff.regions$bio_cluster), function(c.name) {
  lola.result <- runLOLA(userSets = diff.regions[diff.regions$bio_cluster == c.name],
                         userUniverse = diff.regions,
                         regionDB = lola.db)
  lola.result$cluster_id <- c.name
  return(lola.result)
})
names(lola.results) <- unique(diff.regions$bio_cluster)

# Define cellType dict and corresponding colours
cell.type <- c('neural progenitor cells', 'Embryonic Stem Cell', 'Embryonic Stem cell',
               'Embryonic stem cells', 'Embryonic stem cell', 'Neural Tube Cells',
               'Embryonic Stem Cells', 'Motor neuron progenitors', 'ES-E14',
               'Neuron', 'WholeBrain', 'D5-H', 'D5-SC', 'EpiSCs', 'Motorneuron',
               'NEB', 'NPC', 'MNP', 'NMP', 'neuraltube', 'NSC')
cell.short <- c('NPC', 'ESC', 'ESC', 'ESC', 'ESC', 'NT', 'ESC', 'MNP', 'ESC',
                'Neuron', 'WholeBrain', 'D5-H', 'D5-SC', 'EpiSCs', 'MN',
                'NEB', 'NPC', 'MNP', 'NMP', 'NT', 'NSC')

# Define colours for cellType
lola.cells <- c('EpiSCs', 'NEB', 'NPC', 'MNP', 'NT', 'NSC', 'Neuron', 'WholeBrain')
lola.colPal <- brewer.pal(length(lola.cells), 'Paired')
names(lola.colPal) <- lola.cells
lola.colPal <- c('D5-H' = '#FDCD61', 'D5-SC' = '#FB564A', 'ESC' = '#B4B4B4',
                 'NMP'= '#F79992', 'other cell type' = 'black', lola.colPal)
cell.dict <- data_frame(cell_type = cell.type, cell_short = cell.short)

# Plot for each som cluster the enrichment seperatly
gg.lola <- lapply(lola.results, function(lola.result) {

  # Parse name of significantly enriched TFs
  antibody.ids <- lola.result %>%
    filter(-log10(qValue) >= -log10(0.01)) %>%
    mutate(antibody = gsub('_.*', '', antibody),
           antibody = gsub('Ctcf', 'CTCF', antibody),
           antibody = gsub('SMC3', 'Smc3', antibody)) %>%
    filter(!is.na(antibody)) %>%
    filter(!is.infinite(oddsRatio)) %>%
    select(antibody) %>% unlist() %>% unique()

  # If more then 20 different TFs are enriched, take top 20
  if(length(antibody.ids) > 20) {
    print(unique(lola.result$cluster_id))
    antibody.ids <- antibody.ids[1:15]
  }

  # Plot LOLA enrichment
  # Remove GSM968709_V5.bed redundant Cdx2 CODEX vs our own database
  lola.result %>%
    tbl_df() %>%
    filter(-log10(qValue) >= -log10(0.01)) %>%
    filter(filename != 'GSM968709_V5.bed') %>%
    filter(antibody %in% antibody.ids) %>%
    mutate(x = gsub('_.*', '', cluster_id),
           y = gsub('.*_', '', cluster_id),
           antibody = gsub('_.*', '', antibody),
           antibody = gsub('Ctcf', 'CTCF', antibody),
           antibody = gsub('SMC3', 'Smc3', antibody)) %>%
    filter(antibody %in% antibody.ids) %>%
    left_join(., cell.dict, by = c('cellType' = 'cell_type')) %>%
    mutate(cell_short = if_else(is.na(cell_short), 'other cell type', cell_short)) %>%
    ggplot(aes(x = reorder(antibody, -log10(qValue)), y = -log10(qValue),
               size = oddsRatio, col = cell_short)) +
    geom_hline(yintercept = -log10(0.01), col = 'red', linetype = 'dashed') +
    geom_jitter(width = 0.1) + xlab('') +
    #scale_color_brewer(palette = 'Paired', name = 'Cell types') +
    scale_color_manual(values = lola.colPal, name = 'Cell types')+
    guides(colour = guide_legend(override.aes = list(size = 4))) +
    scale_size(name = 'Odds ratio') + ylab('-log10(adj. pvalue)') +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
})
names(gg.lola) <- names(som.bioCluster) #unique(som.cluster)

# Save LOLA enrichment plot for each SOM cluster
#lapply(names(gg.lola), function(c.name){
#  print(c.name)
#  lola.file <- file.path(figures.path, sprintf('lolaEnrichment_%s', c.name))
#  savePlots(gg.lola[[c.name]], file.name = lola.file)
#})

#/*==========================================================================#*/
#' ##  Invitro vs invivo comparison of SC/A ATAC-seq (Figure 2 F)
#+ chunk_deseq2WT_invivo, cache=T, echo=F, warning=F, message=F
#/*==========================================================================#*/
# Define path to invitro peaks count table (with invivo samples)
invivoResults.path <- '/Users/steinhs/projects/atac_vicki/results/atac_conregions/'
invivoCount.txt <- list.files(path = invivoResults.path, full.names = T,
                              pattern = 'conPeaks_counts.txt$')

# Read con region counts by featureCounts
invivo.regionCounts <- read_tsv(invivoCount.txt, skip = 1, col_names = T)
colnames(invivo.regionCounts) <- gsub('_.*', '', basename(colnames(invivo.regionCounts)))

# Parse regions from count table file
invivo.regions <- GRanges(seqnames = invivo.regionCounts$Chr,
                          IRanges(start = invivo.regionCounts$Start,
                                  end = invivo.regionCounts$End))

# Subset count table for only invivo/invitro D5A/SC samples
i.invivo <- grep(colnames(invivo.regionCounts), pattern = '(D5-?A|D5-?SC|invivo)')
invivo.regionCounts <- invivo.regionCounts[,i.invivo]
i.sox2neg <- grep(colnames(invivo.regionCounts), pattern = '(CDX|BRA|Sox2neg)')
if(length(i.sox2neg) > 0) invivo.regionCounts <- invivo.regionCounts[,-i.sox2neg]
invivo.regionCounts <- invivo.regionCounts[-i.rmv,]

# Combine metadata, peaks and count matrix in DESeq2 object
invivo.colData <- data_frame(sample = colnames(invivo.regionCounts)) %>%
  mutate(sample = gsub('-', '_', sample),
         cond = gsub('_[1-3]', '', sample),
         rep = gsub('.*_', '', sample))

# Create 'DESeqData' object
invivo.deseq <- DESeqDataSetFromMatrix(countData = invivo.regionCounts,
                                       colData = invivo.colData,
                                       design = ~ cond)

# Normalise counts using DESeq2
invivo.deseq <- estimateSizeFactors(invivo.deseq)

### Plot LogFC plot comparing invivo/invitro
# Prepare logFC invivo/invitro from norm count table
invivoInvitro.fc <- counts(invivo.deseq, normalized = T)[i.diff,] %>%
  tbl_df() %>%
  mutate(peak_id = i.diff, bio_cluster = diff.regions$bio_cluster) %>%
  gather(sample, norm_counts, 1:8) %>%
  mutate(exp = gsub('[0-9]$', '', gsub('-', '', sample))) %>%
  group_by(peak_id, bio_cluster, exp) %>%
  summarise(norm_counts = mean(norm_counts)) %>%
  spread(exp, norm_counts) %>%
  mutate(logFC = log2(D5SC/D5A),
         logFC_invivo = log2(D5SCinvivo/D5Ainvivo))

# Add labels baed on logFC threshold
invivoInvitro.fc <- invivoInvitro.fc %>%
  filter(bio_cluster %in% c('Anterior', 'Spinal cord', 'Neural')) %>%
  mutate('invivo_status' = if_else(logFC_invivo > 0.5, 'SC',
                                   if_else(logFC_invivo < -0.5, 'A', '-')),
         'invitro_status' = if_else(logFC> 0.5, 'SC',
                                    if_else(logFC< -0.5, 'A', '-')))

# Plot invivo vs invitro scatterplot
gg.invivoVsInvitro <- invivoInvitro.fc %>%
  filter(bio_cluster %in% c('Anterior', 'Spinal cord')) %>%
  ggplot(aes(x = logFC, y = logFC_invivo, col = bio_cluster)) + geom_point() +
  geom_point(data = filter(invivoInvitro.fc, bio_cluster == 'Neural'),
             aes(x = logFC, y = logFC_invivo, col = bio_cluster), inherit.aes = F) +
  scale_color_manual(values = bioCluster.color, name = 'Cluster') +
  xlim(-4, 4.5) + ylim(-4, 4.5) + xlab('invitro log2(D5SC/D5A)') +
  ylab('invivo log2(D5SC/D5A)') +
  geom_hline(yintercept = 0, col = 'black', linetype = 'dashed') +
  geom_vline(xintercept = 0, col = 'black', linetype = 'dashed') +
  theme(legend.position = c(0.15, 0.75),
        plot.margin = margin(t = -0.75, r = -0.75, unit = 'cm'))

# Define wilcox test schematic
box.comp <- list(c('Anterior', 'Neural'),
                 c('Spinal cord', 'Neural'),
                 c('Anterior', 'Spinal cord'))

# Plot invitro SC/A foldchange as boxplot and perform wilcox.test comparison
gg.invitroBox <- invivoInvitro.fc %>%
  ggplot(aes(x = bio_cluster, y = logFC, col = bio_cluster)) + geom_boxplot(notch = T) +
  scale_color_manual(values = bioCluster.color, name = 'Cluster') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  stat_compare_means(comparisons = box.comp, method = 'wilcox.test',
                     label.y = c(2.5, 4, 4.5)) +
  xlab('') + ylab('') + ylim(-4, 4.5) + coord_flip() +
  theme(legend.position = 'none',
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(r = -0.75, unit = 'cm'))

# Plot invivo SC/A foldchange as boxplot and perform wilcox.test comparison
gg.invivoBox <- invivoInvitro.fc %>%
  ggplot(aes(x = bio_cluster, y = logFC_invivo, col = bio_cluster)) + geom_boxplot(notch = T) +
  scale_color_manual(values = bioCluster.color, name = 'Cluster') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  stat_compare_means(comparisons = box.comp, method = 'wilcox.test',
                     label.y = c(2.35, 2.85, 3.35)) +
  ylab('') + ylim(-4, 4.5) +
  theme(legend.position = 'none',
        axis.text.x = element_text(color = 'white'),
        axis.title.x = element_text(color = 'white'),
        axis.ticks.x = element_line(color = 'white'),
        plot.margin = margin(t = -0.75, unit = 'cm'))

# Combine scatter and boxplots
gg.scatter <- plot_grid(gg.invitroBox, gg.invivoVsInvitro, nrow = 2,
                        rel_heights = c(1, 3.25), align = 'v')
gg.scatter2 <- plot_grid(NULL, gg.invivoBox, nrow = 2,
                         rel_heights = c(1, 3.25), align = 'v')
gg.invivoVsInvitro <- plot_grid(gg.scatter, gg.scatter2, ncol = 2,
                                rel_widths = c(3, 1), align = 'h', axis = 'b')
gg.invivoVsInvitro

# Save invivo/invitro scatterplot
#savePlots(gg.invivoVsInvitro,
#          file.name = file.path(figures.path, 'InvivoVsInvitro_scatterPlot'),
#          base_height = 7.5, base_width = 7.5)

#/*==========================================================================#*/
#' ##  Invitro vs invivo ATAC-seq metaprofiles (Figure 2 G)
#+ chunk_deseq2WT_invivoMetaProfile, cache=T, echo=F, warning=F, message=F
#/*==========================================================================#*/
### Plot meta-profile for selected regions
i.a <- invivoInvitro.fc %>%
  ungroup() %>%
  filter(bio_cluster == 'Anterior' &
           invivo_status == 'A' &
           invitro_status == 'A') %>%
  select(peak_id) %>% unlist %>% as.numeric

i.sc <- invivoInvitro.fc %>%
  ungroup() %>%
  filter(bio_cluster == 'Spinal cord' &
           invivo_status == 'SC' &
           invitro_status == 'SC') %>%
  select(peak_id) %>% unlist %>% as.numeric

### Plot in-vivo data for clusters
# Merged invitro/invivo *.bigWigs
bw.files <- atacR::listFiles(file.path(results.path, '*/bw'), pattern = 'D5.*bw$')
bw.files <- bw.files[!grepl(bw.files, pattern = 'KO')]
bw.files <- bw.files[!grepl(bw.files, pattern = '(rep|ind|plus|-H)')]

# Comput consensus regions around summit
profile.summits <- lapply(list('A' = i.a, 'SC' = i.sc), function(i) {
  getSummit(regions[i,]) + 2700
})
names(profile.summits) <- c('Anterior', 'Spinal cord')

# Compute metaprofile for *.bw over con. regions
invivo.atacProfiles <- lapply(bw.files, function(bw.file) {
  print(bw.file)
  invivo.atacProfile <- lapply(names(profile.summits), function(rn) {
    print(rn)
    invivo.profile <- computeProfile(profile.summits[[rn]], bw.file, region.names = rn)
    return(invivo.profile)
  })
  invivo.atacProfile <- bind_rows(invivo.atacProfile)
  invivo.atacProfile$sample <- gsub('_rep.*', '', basename(bw.file))
  return(invivo.atacProfile)
})
invivo.atacProfiles <- bind_rows(invivo.atacProfiles)

# Plot invivo metaprofiles over invitro/invivo shared regions
gg.invivoProfiles <- invivo.atacProfiles %>%
  mutate(sample = gsub('_merged.*', ' ', sample)) %>%
  ggplot(aes(x = rel_pos/10^3, y = mean_score, col = name, fill = name)) +
  geom_line() +  scale_color_manual(values = bioCluster.color, name = 'Cluster') +
  scale_fill_manual(values = bioCluster.color, name = 'Cluster') +
  geom_ribbon(aes(ymin = mean_score - std_error, ymax = mean_score + std_error), alpha = 0.25) +
  xlab('Distance to summit (kbp)') + geom_vline(xintercept = 0, linetype = 'dashed') +
  xlim(-2.5, 2.5) + ylab('Mean FPM') + facet_wrap(~sample, scales = 'free') +
  theme(legend.position = 'bottom')

gg.invivoProfiles

# Save metaprofile of ATAC-seq data over consistent hindbrain/spinalcord regions
#savePlots(gg.invivoProfiles,
#          file.name = file.path(figures.path, 'invivoMetaProfile_logFCThreshold'),
#          base_height = 8, base_width = 9)

#/*==========================================================================#*/
#' ##  Motif enrichment analysis with HOMER (Figure S3)
#+ chunk_deseq2WT_homer, cache=T, echo=F, warning=F, message=F
#/*==========================================================================#*/
# Define and create motif figures path
motifFigure.path <- file.path(figures.path, 'motif_enrichment')
if(!dir.exists(motifFigure.path)) dir.create(motifFigure.path)

# Define path to homer enrichment files
homer.path <- file.path(results.path, 'denovo_motifs/homer')
homer.txts <- atacR::listFiles(homer.path, pattern = 'knownResults.txt')

# Read HOMER known motif enrichment
motif.enrichment <- data_frame(files = homer.txts) %>%
  mutate(sample = gsub('.*Regions_', '', gsub('\\/known.*', '', files)),
         data = map(files, ~ read_tsv(.))) %>%
  select(-files) %>% unnest() %>%
  rename(motif_name = `Motif Name`,
         pvalue =  `P-value`,
         log_pvalue =  `Log P-value`,
         qvalue = `q-value (Benjamini)`)

# Filter enrichment of known motifs for top 30
bc <- 'Hindbrain'
motifEnrich.df <- motif.enrichment %>%
  filter(sample == bc) %>%
  mutate(r = 1:length(pvalue),
         motif_id = motif_name,
         motif_name = gsub('\\/.*', '', motif_name)) %>%
  filter(r <= 30) %>%
  filter(!grepl(motif_name, pattern = '(Unknown|X-|\\:|half|COUP)')) %>%
  group_by(motif_name) %>%
  filter(max(-log_pvalue) == -log_pvalue) %>%
  ungroup() %>%
  mutate(tf_type = gsub('(.*\\(|\\).*)', '', motif_name),
         motif_name = str_to_title(gsub('\\(.*', '', motif_name)))

# Define TF synonyme dict and replace it motif enrichment table
tf.dict <- c('Scl' = 'Tal1', 'E2a' = 'Tcf3', 'Erra' = 'Esrra', 'Myod' = 'Myod1',
             'Rxr' = 'Rxra', 'Ap4' = 'Tfap4', 'Heb' = 'Tcf12', 'Ppare' = 'Pparg',
             'Ronin' = 'Thap11', 'Znf467' = 'Zfp467')

# Filter human only Znf and rename TFs that they fit gencode
motifEnrich.df <- motifEnrich.df %>%
  mutate(motif_name = if_else(motif_name %in% names(tf.dict),
                              tf.dict[motif_name], motif_name)) %>%
  group_by(motif_name) %>%
  filter(max(-log_pvalue) == -log_pvalue) %>%
  ungroup() %>%
  filter(!(motif_name %in% c('Znf263', 'Znf416')))

# Plot known motif enrichment for a given cluster
gg.homer <- motifEnrich.df %>%
  ggplot(aes(x = reorder(motif_name, -log_pvalue), y = -log_pvalue, fill = tf_type)) +
  geom_bar(stat = 'identity') + coord_flip() + ylab('-log(pvalue)') +
  xlab('') + scale_fill_brewer(palette = 'Dark2', name = 'TF class') +
  theme(legend.position = 'bottom')

# Parse gene_ids by overlapping gene_name and TF names
motif.tfs <- motifEnrich.df %>%
  mutate(motif_name = gsub('\\(.*', '', motif_name),
         motif_name = str_to_title(motif_name)) %>%
  select(motif_name) %>% unlist()

motif.tfs <- gencode.metadata  %>%
  filter(gene_name %in% motif.tfs) %>%
  select(gene_id, gene_name) %>% unique() %>%
  mutate(gene_id = gsub('\\..*', '', gene_id),
         gene_name = factor(gene_name, levels = rev(motif.tfs)))

# Get TPM for D5A|H|SC for TF of interest
i.tf <- if_else(bc == 'NMP', 'D3',
        if_else(bc == 'Epi', '(D2_|D3)', 'D5(A|H|SC)'))
i.tf <- grep(colnames(tpms), pattern = i.tf)
motif.tpms <- log2(tpms[motif.tfs$gene_id, i.tf] + 1)

# Plot TPM for corresponding TFs
gg.tfExp <- motif.tpms %>%
  tbl_df() %>%
  bind_cols(., motif.tfs) %>%
  gather(sample, tpm, 1:length(i.tf)) %>%
  separate(sample, c('exp', 'rep'), sep = '_') %>%
  group_by(gene_name, exp) %>%
  summarise(mean_tpm = mean(tpm),
            stderr_tpm = atacR:::stdError(tpm)) %>%
  ggplot(aes(x = gene_name, y = mean_tpm, fill = exp, group = exp)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_tpm - stderr_tpm,
                    ymax = mean_tpm + stderr_tpm, group = exp),
                width = 0.3, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = color.pal, name = '') +
  coord_flip() + xlab('') + ylab('log2(TPM + 1)') +
  theme(legend.position = 'bottom')

# Plot top 6 HOMER motifs
homerMotif.path <- '/Users/steinhs/software/homer/motifs/'
homerMotif.meta <- read_tsv('/Users/steinhs/software/homer/motifs/table.txt')

# Parse Homer metadata file and filter with top 6 motif ids
homerMotif.files <- homerMotif.meta %>%
  filter(Name %in% motifEnrich.df$motif_id[1:6]) %>%
  left_join(., motifEnrich.df[,c('motif_name','motif_id')], by = c('Name' = 'motif_id')) %>%
  select(Filename, Name, motif_name) %>%
  mutate(Filename = file.path(homerMotif.path, Filename))

# Read homer motifs
homer.motifs <- lapply(homerMotif.files$Filename, function(homerMotif.file) {
  homer.motif <- read_tsv(homerMotif.file, skip = 1, col_names = F) %>%
    rename(A = X1, C = X2, G = X3, T = X4) %>% t(.)
  return(homer.motif)
})
names(homer.motifs) <- homerMotif.files$motif_name

# Plot top 6 homer motifs as seqLogos
gg.seqlogos <- ggseqlogo::ggseqlogo(homer.motifs, ncol = 2) +
  ylim(0, 2) + theme_cowplot()

gg.homerSummary <- plot_grid(gg.homer, gg.tfExp, align = 'hv', ncol = 2)
gg.homerSummary <- plot_grid(gg.homerSummary, gg.seqlogos, align = 'hv',
                             ncol = 2, rel_widths = c(1.75, 1))
gg.homerSummary

#savePlots(gg.plot = gg.homerSummary,
#          file.name = file.path(figures.path, sprintf('%s_homerMotifPlot', bc)),
#          base_height = 6, base_width = 9)

#/*==========================================================================#*/
#' ## Overlap ENCODE DNAse with SOM/Bio.Cluster (Figure S2 A,B,C,D)
#+ chunk_deseq2WT_encodeDHS, cache=T, echo=F, warning=F, message=F
#/*==========================================================================#*/
# Define path to ENCODE DNAse data and list bed files
dnase.path <- '/Users/steinhs/public_data/encode/mm10/DNAse/'
dnase.beds <- list.files(dnase.path, full.names = T, pattern = 'gz')

# Read meta-data to annotate IDs with sample names
dnase.meta <- read_tsv(file.path(dnase.path, 'metadata.tsv'))
dnase.meta <- dnase.meta[,c(1,4,7,11,30)] %>%
  rename(encode_id = `File accession`,
         exp_id = `Experiment accession`,
         sample = `Biosample term name`,
         time = `Biosample Age`,
         rep = `Biological replicate(s)`)

# Import all DNAse peaks
dnase.regions <- lapply(dnase.beds, function(dnase.bed) {
  print(dnase.bed)
  importBroadPeak(dnase.bed)
})
dnase.regions <- GRangesList(dnase.regions)
names(dnase.regions) <- gsub('.bed.gz', '', basename(dnase.beds))

# Combine DHS regions by keeping only regions in all replciates (per exp_id)
dnase.regions <- lapply(unique(dnase.meta$exp_id), function(exp.id) {
  print(exp.id)
  encode.ids <- dnase.meta %>% filter(exp_id == exp.id) %>% select(encode_id) %>% unlist()
  if(length(encode.ids) > 1) {
    robust.peaks <- subsetByOverlaps(dnase.regions[[encode.ids[1]]], dnase.regions[encode.ids[-1]])
    return(robust.peaks)
  } else {
    return(dnase.regions[[encode.ids]])
  }
})
names(dnase.regions) <- unique(dnase.meta$exp_id)
dnase.regions <- GRangesList(dnase.regions)

# Overlap DNAse with SOM regions
ov.dnase <- findOverlaps(diff.regions, dnase.regions)

# Prepare number of biocluster regions
n.biocluster <- table(diff.regions$bio_cluster) %>%
  tbl_df %>%
  rename(bio_cluster = Var1, n_cluster = n)

# Compute DNAse overlap dataFrame - tidy overlap
dnaseOv.df <- data_frame(peak_id = queryHits(ov.dnase),
                         bio_cluster = diff.regions$bio_cluster[queryHits(ov.dnase)],
                         exp_id = names(dnase.regions)[subjectHits(ov.dnase)]) %>%
  unique %>%
  group_by(bio_cluster, exp_id) %>%
  summarise(n = n()) %>%
  left_join(., n.biocluster, by = 'bio_cluster') %>%
  mutate(freq = n/n_cluster) %>% arrange(desc(freq)) %>%
  left_join(., dnase.meta, by = 'exp_id') %>%
  select(-encode_id, -rep) %>% unique %>%
  filter(!grepl(sample, pattern = 'CD71')) %>%
  mutate(sample = str_to_title(sample),
         time = if_else(grepl(time, pattern = 'unknown'), 'unknown', time),
         time = if_else(grepl(time, pattern = '50'), gsub('50', '5', time), time),
         time = if_else(time == '0 day', 'postnatal', time))

# Get ENCODE time categories and sort in logical order
dnase.sortedTime <- dnaseOv.df$time %>% unique %>% sort
dnase.sortedTime <- dnase.sortedTime[c(8, 10, 1:7, 9, 11:12)]

dnase.colPal <- c('#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99',
                  '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a',
                  '#ffff99', '#b15928')

n.dnase <- dnaseOv.df$sample %>% unique %>% length
n.dnase <-  seq(1.5, n.dnase + 0.5, by = 1)

# Plot DHS overlap as jitter plot for each biocluster
lapply(names(som.bioCluster), function(bc) {
  print(bc)

  gg.dnaseOv <- dnaseOv.df %>%
    filter(bio_cluster == bc) %>%
    group_by(sample) %>%
    mutate(time = factor(time, levels = dnase.sortedTime),
           max_freq = max(freq)) %>%
    ggplot(aes(x = reorder(sample, max_freq), y = freq*100, fill = time)) +
    geom_jitter(size = 4, shape = 21, width = 0.3) + xlab('') + #ylim(0, 100) +
    #geom_point(size = 4, shape = 21) + xlab('') +
    scale_fill_manual(values = dnase.colPal, name = '') +
    ylab('Proportion of regions overlap DHS (%)') +
    geom_vline(xintercept = n.dnase, linetype = 'dashed') +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

  # Save ENCODE DHS SOM overlap plot
  #savePlots(gg.plot = gg.dnaseOv,
  #          file.name = file.path(figures.path, sprintf('%s_encodeDHSOverlap', bc)),
  #          base_height = 6, base_width = 9)
})

# DNAseq overlap proportion heatmap
#dnaseOv.df %>%
#  filter(bio_cluster == 'Neural') %>%
#  mutate(time = factor(time, levels = dnase.sortedTime)) %>%
#  group_by(bio_cluster, sample, time) %>%
#  summarise(freq = max(freq)) %>%
#  ggplot(aes(x = reorder(sample, freq), y = time, fill = freq)) +
#  geom_tile() + facet_wrap(~bio_cluster) +
#  scale_fill_gradient(low = '#feedde', high = '#a63603',
#                      name = 'Proportion (%)', limits = c(0, 1)) +
#  ylab('Proportion of regions overlap DHS (%)') +
#  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

#/*==========================================================================#*/
#' ##  VISTA enhancer enrichment of SOM regions (Figure S2 E)
#+ chunk_deseq2WT_vista, cache=T, echo=F, warning=F, message=F
#/*==========================================================================#*/
# Path to Vista enhancer *.bed
vista.bed <- '/Users/steinhs/public_data/vista/vista_enhancer_mm10_ENCFF239XXL.bed'

# Preprocess Vista Enhancers
vista.enhancers <- read_tsv(vista.bed, col_names = F)
vista.enhancers <- vista.enhancers %>%
  mutate(anno = gsub('\\[.*?\\]', '', vista.enhancers$X10),
         anno = gsub('\\(.*?\\)', '', anno),
         anno = str_split(anno, pattern = ',')) %>%
  unnest(anno) %>%
  mutate(anno = gsub('(^ | $)', '', anno)) %>%
  filter(!grepl(anno, pattern = 'negative'))

# Convert to GenomicRanges
vista.gr <- GRanges(seqnames = vista.enhancers$X1,
                    ranges = IRanges(vista.enhancers$X2, vista.enhancers$X3),
                    anno = vista.enhancers$anno)

# DEF. function to compute the proportion the 'annotated genome' given an
# annotation GenomicRanges and a BSgenome object
coveredGenome <- function(gr, bs.genome) {
  genome.size <- seqinfo(bs.genome) %>% seqlengths()
  genome.size <- genome.size[-grep(names(genome.size), pattern = '(_|M)')]
  genome.size <- sum(as.numeric(genome.size))
  return(sum(width(gr))/genome.size)
}

# Test VISTA enhancer enrichment using a binomial Test
vista.enrichment <- lapply(n.biocluster$bio_cluster, function(bc) {
  print(bc)
  # Perform enrichment testing for given bio_cluster
  vista.enrich <- lapply(unique(vista.gr$anno), function(vista.anno) {
    # Compute proportion of annotated genome
    p <- coveredGenome(gr = vista.gr[vista.gr$anno == vista.anno],
                       bs.genome = BSgenome.Mmusculus.UCSC.mm10)
    # Compute number of regions within SOM cluster
    n.bc <- length(diff.regions[diff.regions$bio_cluster == bc])
    # Overlap SOM regions with VISTA enhancers
    ov.vista <- findOverlaps(diff.regions[diff.regions$bio_cluster == bc],
                             vista.gr[vista.gr$anno == vista.anno])
    # Test enrichment using a binomial test
    binom.test(x = length(ov.vista), n = n.bc, p = p, alternative = 'greater') %>%
    tidy() %>%
    mutate(vista_anno = vista.anno,
           prop = length(ov.vista)/n.bc,
           bio_cluster = bc)
  })
  return(bind_rows(vista.enrich))
})

# PLot VISTA enhancer enrichment
gg.vista <- bind_rows(vista.enrichment) %>%
  mutate(vista_anno = if_else(grepl(vista_anno, pattern = 'crest'), 'mesenchyme (neural crest)', vista_anno)) %>%
  mutate(padj = p.adjust(p.value, method = 'fdr')) %>%
  filter(bio_cluster %in% c(names(bioCluster.color), 'Epi')) %>%
  ggplot(aes(x = reorder(vista_anno, prop), y = 100*prop,
             fill = bio_cluster, size = -log10(padj))) +
  geom_point(shape = 21) + xlab('') +
  scale_fill_manual(values = bioCluster.color, name = 'Cluster') +
  guides(fill = guide_legend(override.aes = list(size = 4))) +
  scale_size(name = '-log10(adj pvalue)') +
  ylab('Proportion of regions \n overlapping VISTA enhancer (%)') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
gg.vista

#savePlots(gg.vista, file.name = file.path(figures.path, 'VistaEnhancer_enrichment'),
#          base_width = 8, base_height = 6)

#/*==========================================================================#*/
#' ##  Regional/Axis identity regions over time plot (Figure 3 A,B,C,D)
#+ chunk_deseq2WT_regionIdentity, cache=T, echo=F, warning=F, message=F
#/*==========================================================================#*/
# Plot ATAC-seq signal vs time for A/H/SC bioclusters
axisIdentiy.cluster <- c('Anterior', 'Hindbrain', 'Spinal cord', 'Neural')

# Compute relative accesibilty as z-Scores
axisIdentity.matrix <- zscore.matrix[diff.regions$bio_cluster %in% axisIdentiy.cluster,] %>%
  tbl_df() %>%
  mutate(bio_cluster = diff.regions$bio_cluster[diff.regions$bio_cluster %in% axisIdentiy.cluster],
         peak_id = diff.regions$id[diff.regions$bio_cluster %in% axisIdentiy.cluster])

axisIdentity.matrix <- axisIdentity.matrix %>%
  gather(exp, value, 1:11) %>%
  mutate(day = str_extract(exp, pattern = 'D[0-9]'))

# Reformat matrix to dataframe
axisIdentity.df <- lapply(c('A', 'H', 'SC'), function(ef) {
  exp.filter <- sprintf('(D0|D1|D2|D3A|D4%s|D5%s)', ef, ef)
  if(ef == 'SC') exp.filter <- sprintf('(D0|D1|D2|D3NMP|D4%s|D5%s)', ef, ef)
  axisIdentity.df <- lapply(axisIdentiy.cluster, function(axis.cluster) {
    axisIdentity.matrix %>%
      filter(bio_cluster == axis.cluster) %>%
      filter(grepl(exp, pattern = exp.filter)) %>%
      mutate(exp_cond = ef)
  })
  axisIdentity.df <- bind_rows(axisIdentity.df)
  return(axisIdentity.df)
})
axisIdentity.df <- bind_rows(axisIdentity.df)

# Plot zscores over days for each cluster seperatly
lapply(c('Anterior', 'Hindbrain', 'Spinal cord', 'Neural'), function(bc) {
  gg.axisIdentity <- axisIdentity.df %>%
    group_by(bio_cluster, day, exp, exp_cond) %>%
    summarise(std_error = sd(value), value = mean(value)) %>%
    filter(day != 'D0') %>%
    ungroup() %>%
    mutate(day = gsub('D', '', day)) %>%
    filter(bio_cluster == bc) %>%
    ggplot(aes(x = day, y = value, group = 1, shape = exp_cond)) + #geom_line(size = 1) +
    geom_errorbar(aes(ymin = value - std_error, ymax = value + std_error), width = 0.5) +
    geom_point(size = 2, fill = 'black') +
    scale_shape_manual(values = c(21, 22, 24)) +
    geom_hline(yintercept = 0, col = 'red', linetype = 'dashed') +
    facet_grid(. ~ exp_cond) + xlab('Day') +
    ylab('Relative chromatin accessibility') + theme(legend.position = 'none')

  #savePlots(gg.plot = gg.axisIdentity,
  #          file.name = file.path(figures.path, sprintf('axisIdentity_atacSignal_withoutD0_%s', bc)))
})
