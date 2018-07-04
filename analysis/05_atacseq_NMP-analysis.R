#/*==========================================================================#*/
#' ## Overlap NMP sides with ATAC peaks from published studies (Figure 4 A)
#+ chunk_NMP_peakOV, cache=T, echo=F, warning=F, message=F
#/*==========================================================================#*/
# Find peak files for inhouse invivo data
pubAtac.files <- atacR::listFiles(file.path(data.path, 'atacseq/vicki/*invivo*/macs2'),
                                  pattern = 'broadPeak')

# Define path published atacseq/chipseq data and get path to relevant files
pubAtac.files <- lapply(c('amin2017_GSE84899', 'neijts2016_GSE81203'), function(study) {
  study.path <- file.path(data.path, '*', study, '*/macs2')
  atacR::listFiles(study.path, pattern = 'broadPeak')
}) %>% unlist %>% c(., pubAtac.files)

# Import ATAC-seq/ChIP-seq regions (broad|narrowPeak)
pubAtac.regions <- lapply(pubAtac.files, function(pubAtac.file) {
  print(pubAtac.file)
  if(grepl(pubAtac.file, pattern = 'broadPeak')) {
    importBroadPeak(pubAtac.file)
  } else {
    importNarrowPeak(pubAtac.file)
  }
})
names(pubAtac.regions) <- basename(gsub('_pea.*', '', basename(pubAtac.files)))
names(pubAtac.regions) <- gsub('-1', '_rep1', names(pubAtac.regions))
names(pubAtac.regions) <- gsub('-2', '_rep2', names(pubAtac.regions))
pubAtac.regions <- GRangesList(pubAtac.regions)

# Merge broadPeaks from replicates
pubAtac.names <- gsub('_rep.*', '', names(pubAtac.regions))
pubAtac.regions <- lapply(pubAtac.names, function(pn) {
  i.rep <- grep(names(pubAtac.regions), pattern = pn)
  if(length(i.rep) > 1) {
    #GenomicRanges::reduce(unlist(pubAtac.regions[i.rep]))
    subsetByOverlaps(pubAtac.regions[[i.rep[1]]],
                     pubAtac.regions[[i.rep[2]]])
  } else {
    pubAtac.regions[[i.rep]]
  }
})
names(pubAtac.regions) <- pubAtac.names

# Compute number of overlaps between peaks from public data and NMP regions
ov.pubAtac <- lapply(names(pubAtac.regions), function(pr.n){
  pubAtac.region <- pubAtac.regions[[pr.n]]
  ov.pubAtac <- lapply(c('2_4', '3_4', '4_4'), function(sc) {
    n <- sum(countOverlaps(diff.regions[diff.regions$cluster == sc], pubAtac.region) != 0)
    data_frame(cluster = sc, exp = pr.n, n_ov = n, n = sum(diff.regions$cluster == sc))
  })
  return(bind_rows(ov.pubAtac))
})
ov.pubAtac <- bind_rows(ov.pubAtac)

# Rename sample names with more scientific names
pubAtac.dict <- c('EpiSCs' = 'EpiSCs',
                  'EpiSCs-Chiron48h' = 'EpiSCs-Chiron48h',
                  'embryo-day6.0' = 'E6.0',
                  'embryo-day7.2' = 'E7.2',
                  'embryo-day7.5-post' = 'E7.5-post',
                  'D5-SC-invivo' = 'E9.5-SC',
                  'D5-A-invivo' = 'E9.5-A')

# Plot region public dataset overlap
gg.pubAtacOv <-  ov.pubAtac %>%
  mutate(freq = (n_ov/n)*100) %>%
  filter(!grepl(exp, pattern = '(Sox2|Cdx)')) %>%
  mutate(exp = factor(pubAtac.dict[exp], levels = pubAtac.dict)) %>%
  ggplot(aes(x = exp, y = freq, fill = cluster)) +
  geom_bar(stat = 'identity', position = 'dodge') + xlab('') +
  scale_fill_brewer(palette = 'Dark2', name = 'Cluster') +
  ylab('Proportion of NMP regions (%)') + ylim(0, 100) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = c(0.8, 0.5))

gg.pubAtacOv

#savePlots(gg.pubAtacOv,
#          file.path(figures.path, 'NMPCluster_pubAtacOV.pdf'),
#          base_height = 5, base_width = 5)

#/*==========================================================================#*/
#' ## D3 samples meta-profile over NMP sites (Figure 4 B)
#+ chunk_nmp_d3MetaProfile, cache=T, echo=F, warning=F, message=F
#/*==========================================================================#*/
# Define path to bigWig files
d3.bws <- atacR::listFiles(file.path(results.path, '*KO/bw'), pattern = 'D3')

# Compute +/-2700bp around summits for NMP regions
#nmp.summits <- getSummit(diff.regions[diff.regions$bio_cluster == 'NMP']) + 2700
nmp.summits <-  lapply(c('2_4', '3_4'), function(sc) {
  getSummit(diff.regions[diff.regions$cluster == sc]) + 2700
})
names(nmp.summits) <- c('2_4', '3_4')

# Compute metaprofile for all D3 *.bw over NMP regions
d3.atacProfiles <- lapply(d3.bws, function(d3.bw) {
  print(d3.bw)
  d3.profiles <- lapply(names(nmp.summits), function(sc) {
    print(sc)
    d3.profile <- computeProfile(nmp.summits[[sc]], d3.bw, region.names = sc)
    d3.profile$sample <- gsub('_merged.*', '', basename(d3.bw))
    return(d3.profile)
  })
  d3.profiles <- bind_rows(d3.profiles)
  return(d3.profiles)
})
d3.atacProfiles <- bind_rows(d3.atacProfiles)

# Plot meta-profile for 2_4/3_4 nmp sites for d3 samples including Cdx ko
ko.design <- c('CDX' = 'BRA', 'BRA' = 'CDX')
lapply(ko.design, function(ko) {
  print(ko)
  gg.d3KoAtacProfiles <- d3.atacProfiles %>%
    mutate(sample = sub('-', '', sample)) %>%
    filter(!grepl(sample, pattern = ko)) %>%
    ggplot(aes(x = rel_pos/10^3, y = mean_score, col = sample, fill = sample)) +
    geom_line() + scale_color_manual(values = color.pal, name = 'Sample') +
    scale_fill_manual(values = color.pal, name = 'Sample') +
    geom_ribbon(aes(ymin = mean_score - std_error, ymax = mean_score + std_error), alpha = 0.25) +
    xlab('Distance to summit (kbp)') + geom_vline(xintercept = 0, linetype = 'dashed') +
    xlim(-2.5, 2.5) + ylab('Mean FPM') + facet_wrap(~name, scales = 'free_x') +
    theme(legend.position = 'bottom')

  gg.d3KoAtacProfiles
  # Save meta-profile to file
  #out.file <- file.path(figures.path, sprintf('d3%sKO_metaprofile', ko.design[ko]))
  #savePlots(gg.plot = gg.d3KoAtacProfiles, file.name = out.file)
})

#/*==========================================================================#*/
#' ## D5 samples meta-profile over HB/SC sites (Figure 4C, 4N)
#+ chunk_nmp_d5MetaProfile, cache=T, echo=F, warning=F, message=F
#/*==========================================================================#*/
# Define path to bigWig files
d5.bws <- atacR::listFiles(file.path(results.path, '*KO/bw'), pattern = 'D5')

# Compute +/-2700bp around summits for Hindbrain/Spinal cord regions
d5.summits <-  lapply(c('Hindbrain', 'Spinal cord'), function(bc) {
  d5.summits <- getSummit(diff.regions[diff.regions$bio_cluster == bc]) + 2700
})
names(d5.summits) <- c('Hindbrain', 'Spinal cord')

# Compute metaprofile for all D5 *.bw over Hindbrain/Spinal cord regions
d5.atacProfiles <- lapply(d5.bws, function(d5.bw) {
  print(d5.bw)
  d5.profiles <- lapply(names(d5.summits), function(bc) {
    print(bc)
    d5.profile <- computeProfile(d5.summits[[bc]], d5.bw, region.names = bc)
    d5.profile$sample <- gsub('_merged.*', '', basename(d5.bw))
    return(d5.profile)
  })
  d5.profiles <- bind_rows(d5.profiles)
  return(d5.profiles)
})
d5.atacProfiles <- bind_rows(d5.atacProfiles)

### SPINAL CORD
# Plot Spinal cord metaprofiles for all D5 inc KO samples
gg.d5SCAtacProfiles <- d5.atacProfiles %>%
  filter(name == 'Spinal cord' & !grepl(sample, pattern = 'BRA')) %>%
  mutate(sample = sub('-', '', sample)) %>%
  ggplot(aes(x = rel_pos/10^3, y = mean_score, col = sample, fill = sample)) +
  geom_line() + scale_color_manual(values = color.pal, name = 'Sample') +
  scale_fill_manual(values = color.pal, name = 'Sample') +
  geom_ribbon(aes(ymin = mean_score - std_error, ymax = mean_score + std_error), alpha = 0.25) +
  xlab('Distance to summit (kbp)') + geom_vline(xintercept = 0, linetype = 'dashed') +
  xlim(-2.5, 2.5) + ylab('Mean FPM')
gg.d5SCAtacProfiles

# Save metaprofile to file
#savePlots(gg.plot = gg.d5SCAtacProfiles,
#          file.name = file.path(figures.path, 'D5KO_SC_metaprofile'))

# Plot Spinal cord metaprofiles for all D5 inc Bra KO samples
gg.d5BraSCAtacProfiles <- d5.atacProfiles %>%
  filter(name == 'Spinal cord' & !grepl(sample, pattern = 'CDX')) %>%
  mutate(sample = sub('-', '', sample)) %>%
  ggplot(aes(x = rel_pos/10^3, y = mean_score, col = sample, fill = sample)) +
  geom_line() + scale_color_manual(values = color.pal, name = 'Sample') +
  scale_fill_manual(values = color.pal, name = 'Sample') +
  geom_ribbon(aes(ymin = mean_score - std_error, ymax = mean_score + std_error), alpha = 0.25) +
  xlab('Distance to summit (kbp)') + geom_vline(xintercept = 0, linetype = 'dashed') +
  xlim(-2.5, 2.5) + ylab('Mean FPM')
gg.d5BraSCAtacProfiles

# Save metaprofile to file
#savePlots(gg.plot = gg.d5BraSCAtacProfiles,
#          file.name = file.path(figures.path, 'D5BraKO_SC_metaprofile'))

### HINDBRAIN
# Plot Hindbrain metaprofiles for all D5 inc KO samples
gg.d5HbAtacProfiles <- d5.atacProfiles %>%
  filter(name == 'Hindbrain' & !grepl(sample, pattern = 'BRA')) %>%
  mutate(sample = sub('-', '', sample)) %>%
  ggplot(aes(x = rel_pos/10^3, y = mean_score, col = sample, fill = sample)) +
  geom_line() + scale_color_manual(values = color.pal, name = 'Sample') +
  scale_fill_manual(values = color.pal, name = 'Sample') +
  geom_ribbon(aes(ymin = mean_score - std_error, ymax = mean_score + std_error), alpha = 0.25) +
  xlab('Distance to summit (kbp)') + geom_vline(xintercept = 0, linetype = 'dashed') +
  xlim(-2.5, 2.5) + ylab('Mean FPM')
gg.d5HbAtacProfiles

#savePlots(gg.plot = gg.d5HbAtacProfiles,
#          file.name = file.path(figures.path, 'D5KO_Hb_metaprofile'))
