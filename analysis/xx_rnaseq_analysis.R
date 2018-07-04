#/*==========================================================================#*/
#' ## Plot gene expression of SOX1 (Figure 3 P)
#+ chunk_rnaseq_expBar, cache=T, echo=F, warning=F, message=F
#/*==========================================================================#*/
# Check if necessary count data is existing
if(!exists('txi.salmon')) {
  txi.salmon <- readRDS(file.path(results.path, 'RNAseq_txiSalmon.Rds'))
}

# Define gene of interest and region size
gene.ofi <- 'Sox1'

# Get ENSEMBL gene id for symbol
gene.id <- gencode.metadata%>%
  filter(gene_name %in% gene.ofi) %>%
  arrange(gene_name) %>% unique(.) %>%
  mutate(gene_id = gsub('\\.[0-9]*', '', gene_id)) %>%
  select(gene_id) %>% unlist() %>% unique()

### Plot expression data as barplot
# Get TPM for gene of interest
tpm.ofi <- log2(tpms[gene.id,] + 1)

# Gene expression as barplot
gg.expBar <- tpm.ofi %>%
  tbl_df() %>%
  mutate(sample = names(tpm.ofi)) %>%
  separate(sample, c('cond', 'rep'), sep = '_') %>%
  mutate(cond = factor(cond, levels = sort(unique(cond), decreasing = F))) %>%
  filter(cond %in% names(color.pal)) %>%
  group_by(cond) %>%
  summarise(mean_tpm = mean(value),
            sd = sd(value),
            std_error = atacR:::stdError(value)) %>%
  mutate(cond = factor(cond, levels = names(color.pal))) %>%
  filter(cond != 'D0') %>%
  ggplot(aes(x = cond, y = mean_tpm, fill = cond)) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymin = mean_tpm - sd, ymax = mean_tpm + sd), col = 'darkgrey', size = 1, width = 0.15) +
  xlab('') + ylab('log2(TPM + 1)') + scale_fill_manual(values = color.pal) +
  coord_flip() +
  theme(legend.position = 'none')
gg.expBar

#savePlots(gg.plot = gg.expBar,
#          file.name = file.path(figures.path, sprintf('ExpBar_%s', gene.ofi)),
#          base_width = 3.5, base_height = 4.5)
#          #base_width = 5.5)

#/*==========================================================================#*/
#' ## Plot Cdx expression profile in spinal cord condition (Figure S4 F)
#+ chunk_rnaseq_expBar, cache=T, echo=F, warning=F, message=F
#/*==========================================================================#*/
# Get Cdx1/2/4 gene_ids
cdx.geneId <- gencode.metadata %>%
  filter(grepl(gene_name, pattern = 'Cdx')) %>%
  arrange(gene_name) %>% unique(.) %>%
  mutate(gene_id = gsub('\\.[0-9]*', '', gene_id)) %>%
  select(gene_id, gene_name) %>% unique()

# Plot expression of Cdx1/2/4 from D0 to D5SC
gg.cdxBar <- log2(tpms[cdx.geneId$gene_id,] + 1) %>%
  reshape2::melt() %>%
  rename(gene_id = Var1, sample = Var2, tpm = value) %>%
  left_join(., cdx.geneId, by = 'gene_id') %>%
  separate(sample, c('exp', 'rep'), sep = '_') %>%
  group_by(gene_name, exp) %>%
  summarise(mean_tpm = mean(tpm), std_error = atacR:::stdError(tpm)) %>%
  filter(!grepl(exp, pattern = '(H|A|D2.5$)')) %>%
  mutate(exp = factor(exp, levels = names(color.pal))) %>%
  ggplot(aes(x = exp, y = mean_tpm, fill = exp)) +
  geom_bar(stat = 'identity') + facet_grid(.~gene_name) + # , scales = 'fre') +
  scale_fill_manual(values = color.pal) +
  geom_errorbar(aes(ymin = mean_tpm - std_error, ymax = mean_tpm + std_error), width = 0.3) +
  xlab('') + ylab('log2(TPM + 1)') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = 'none')

gg.cdxBar
