# Following 
# https://bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html
# https://github.com/hbctraining/DGE_workshop_salmon_online
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# 
# Import files using tximeta
library(DESeq2)
library(tximeta)
library(tidyverse)

dir <- "../../../../linux/mycmax/salmon_alignments"

files <- file.path(list.files(dir, "quant_orig_", full.names = T), "quant.sf")
file.exists(files)

newNames <- data.frame("RNASEQ" = gsub("_quant_orig_index", "", list.files(dir, "quant_orig")),
                       "genotype" = c("lmax_A","lmax_B","lmax_C","lmax_hmyc_A","lmax_hmyc_B","lmax_hmyc_C","neg_A","neg_B","neg_C"))
# RNASEQ              genotype
# 1  RNASEQ1A_S3      lmax_A
# 2  RNASEQ1B_S8      lmax_B
# 3  RNASEQ1C_S6      lmax_C
# 4  RNASEQ2A_S2 lmax_hmyc_A
# 5 RNASEQ2B_S10 lmax_hmyc_B
# 6  RNASEQ2C_S5 lmax_hmyc_C
# 7  RNASEQ3A_S4       neg_A
# 8  RNASEQ3B_S1       neg_B
# 9 RNASEQ3C_S11       neg_C

coldata <- data.frame(
  "names" = c("lmax_A","lmax_B","lmax_C","lmax_hmyc_A","lmax_hmyc_B","lmax_hmyc_C","neg_A","neg_B","neg_C"),
  "files" = files)

fname <- "data/002_gse.rds"
if(!file.exists(fname)){
  
  se <- tximeta(coldata)
  gse <- summarizeToGene(se)
  saveRDS(gse, file = fname)
  
}else{gse <- readRDS(fname)}

colData(gse)
rowRanges(gse)  # this object contains column "symbol" with all genes' symbols
assayNames(gse) # names of the assays (counts = rawcounts, abundance = TPMs, length = gene length)
assays(gse)[["length"]]   
assay(gse)      # genes raw counts
seqinfo(rowRanges(gse))
round(colSums(assay(gse))/1e6, 1) # Ms of fragments mapped

# Data preparation and QC
library(ComplexHeatmap)
library(RColorBrewer)
library(Rtsne)
set.seed(3)

fname <- "data/002_dds.rds"
if(!file.exists(fname)){
  
  gse$gt <- factor(c(rep("lmax", 3), rep("lmax_hmyc", 3), rep ("neg", 3))) %>% relevel("neg")
  dds <- DESeqDataSet(gse, design = ~gt)
  saveRDS(dds, fname)
  nrow(dds)   # 14094
  
}else{dds <- readRDS(fname)}

# keep <- rowSums(counts(dds)) > 1
# dds <- dds[keep,]
# nrow(dds)   # 11191
fname <- ("data/002_dds_rlog_normalized.rds")
if(!file.exists(fname)){
  rld <- rlog(dds, blind = FALSE)   # n<30 so we use rlog instead of vst
  saveRDS(rld, fname)
  all.equal(rowRanges(rld)$gene_id, rownames(rld))    # TRUE
  rownames(rld) <- rowRanges(rld)$symbol
  saveRDS(assay(rld), "data/002_dds_rlog_normalized_assayonly.rds")
  }else{rld <- readRDS(fname)}

sampleDists <- dist(t(assay(rld)))    # Euclidean distance between samples
sampleDistMatrix <- as.matrix(sampleDists)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

png("plots/002_1_Euclidean_distance.png", h = 3000, w = 4200, res = 600)
Heatmap(sampleDistMatrix,
        col = colors,
        column_title = "Overall similarity between samples",
        name = "Euclidean \ndistance\n",
        rect_gp = gpar(col = "grey60", lwd = 1)
)
dev.off()

png("plots/002_2_PCA_plot.png", h = 2000, w = 3500, res = 600)
plotPCA(rld, intgroup = c("names", "gt"), ntop = nrow(rld)) +
  labs(title = "Principal Component Analysis") +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(xlim = c(-20, 25), ylim = c(-20, 25)) + 
  theme_bw()
dev.off()

rldMat <- as.matrix(t(assay(rld)))
tsne <- Rtsne(rldMat, perplexity = floor((nrow(rldMat)-1)/3))
tsnePlot <- data.frame(tSNE_1 = tsne$Y[,1], 
                       tSNE_2 = tsne$Y[,2],
                       Genotype = gsub(".{2}$", "", colnames(rld)))

png("plots/002_3_tSNE.png", h = 3000, w = 4000, res = 600)
ggplot(tsnePlot, aes(x = tSNE_1, y = tSNE_2, color = Genotype)) +
  geom_point(size = 5) +
  labs(title = "tSNE") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# Differential expression analysis
fname <- "data/002_des.rds"
if(!file.exists(fname)){
  des <- DESeq(dds)
  saveRDS(des, fname)
}else{des <- readRDS(fname)}

resultsNames(des)
all.equal(rowRanges(dds)$gene_id, rownames(des))    # TRUE means that genes are all in the same order
write.csv(res[[1]], "results/DESeq_results_lmaxVsNeg.csv")
write.csv(res[[2]], "results/DESeq_results_lmaxhmycVsNeg.csv")
write.csv(res[[3]], "results/DESeq_results_lmaxhmycVslmax.csv")
rownames(des) <- rowRanges(des)$symbol           # Change gene ids to symbols

# Unshrunk results
fname <- "data/002_DESeq2_results.rds"
if(!file.exists(fname)){
  
  # We set alpha (FDR threshold) to the same threshold we usually keep for padj (0.05). 
  # From the argument description of ?results:
  # alpha	
  # the significance cutoff used for optimizing the independent filtering (by default 0.1). 
  # If the adjusted p-value cutoff (FDR) will be a value other than 0.1, alpha should be set to that value.
  # We can save as GRanges by specifying in the results argument "format = "GRanges""
  # We can also save a column with gene symbols by specifying tidy = T, but this will give an error with plotMA
  
    res <- list(
    lmax_vs_neg = results(des, contrast = c("gt", "lmax", "neg"), alpha = 0.05),
    lmax_hmyc_vs_neg = results(des, contrast = c("gt", "lmax_hmyc", "neg"), alpha = 0.05),
    lmax_hmyc_vs_lmax = results(des, contrast = c("gt", "lmax_hmyc", "lmax"), alpha = 0.05)
  )
  saveRDS(res, fname)

}else(res <- readRDS(fname))



lapply(res, summary)

lapply(seq(1:length(res)), function(i){
  
  # Need to create this df if we have a GRanges object or tidy = T
  # df <- data.frame("baseMean" = res[[i]]$baseMean,
  #                  "log2FoldChange" = res[[i]]$log2FoldChange,
  #                  "padj" = res[[i]]$padj < 0.1)
  
  png(paste0("plots/002_4_", names(res)[i], "_MAplot.png"), w = 4000, h = 2500, res = 600)
  plotMA(res[[i]], main = paste("MA plot", names(res)[i], "unshrunk"),
         colSig = "#1E20FF80",
         colNonSig = "#88888870",
         colLine = "#ff000080"
  )
  dev.off()
  
})

# Shrunken results
fname <- "data/002_DESeq2_shrunk_results.rds"
if(!file.exists(fname)){
  
  # For the third contrast (lmaxhmyc vs lmax) we must create a coef by releveling "gt"factors and using lmax as the reference level
  desLmax <- des 
  desLmax$gt <- relevel(desLmax$gt, "lmax")
  desLmax <- nbinomWaldTest(desLmax)
  resLmax <- results(desLmax, contrast = c("gt", "lmax_hmyc", "lmax"), alpha = 0.05)
  coef <- resultsNames(des)           # "Intercept"  "gt_lmax_vs_neg" "gt_lmax_hmyc_vs_neg"
  coefLmax <- resultsNames(desLmax)   # "Intercept" "gt_neg_vs_lmax" "gt_lmax_hmyc_vs_lmax"
  
  shrRes <- list(lmaxVsNeg = lfcShrink(dds = des, res = res[[1]], coef = coef[2], type = "apeglm"),
                 lmaxhmycVsNeg = lfcShrink(dds = des, res = res[[2]],coef = coef[3], type = "apeglm"),
                 lmaxhmycVslmax = lfcShrink(dds = desLmax, res = resLmax, coef = coefLmax[3], type = "apeglm"))
  
  saveRDS(shrRes, fname)
}else{shrRes <- readRDS(fname)}

lapply(seq(1:length(shrRes)), function(i){
  
  png(paste0("plots/002_5_", names(shrRes)[i], "_MAplot.png"), w = 4000, h = 2500, res = 600)
  plotMA(shrRes[[i]], main = paste("MA plot", names(shrRes)[i], "shrunk"),
         colSig = "#1E20FF80",
         colNonSig = "#88888870",
         colLine = "#ff000080"
  )
  dev.off()
  
})

# Check if the dispersion estimates for the 2 DESeq objects are identical (they should, otherwise something is wrong)
png("plots/002_6_Dispersion_estimates_unshrunk.png", w = 4000, h = 3000, res = 600)
plotDispEsts(des, main = "Dispersion estimates")
dev.off()

png("plots/002_6_Dispersion_estimates_lfcShrink.png", w = 4000, h = 3000, res = 600)
plotDispEsts(des2, main = "Dispersion estimates releveled dataset")
dev.off()

# Check summaries for the results. If we saved the results in GenomicRanges, in order to access the DESeqResults 
# object we must use mcols(shrRes[[i]]).  

sapply(shrRes, summary)

# Plot gene counts
best_list <- map(seq(1:length(shrRes)), function(i, ngenes = 4){
  
  res <- shrRes[[i]] %>% 
    as.data.frame() %>% 
    rownames_to_column("gene_name") %>% 
    mutate(gene_id = rowRanges(des)$gene_id) # We can do this because because all.equal(rowRanges(des)$gene_name, rownames(shrRes[[i]])) #[1] TRUE
  
  ord <- res %>% 
    filter(padj < 0.05) %>% 
    arrange(desc(abs(log2FoldChange))) %>% 
    slice_head(n = ngenes) %>% 
    bind_rows(res %>% filter(gene_name %in% c("nub", "pros", "sd", "vg")))
  
  best <- ord %>%
    pull(gene_id) %>% 
    map_dfr(function(x){
    
    geneCounts <- plotCounts(dds = des, 
                             gene = x, 
                             intgroup = c("gt", "names"), 
                             returnData = T) %>%
      mutate(gene = x)
    
    }) %>% 
    inner_join(ord %>% select(gene_name, gene_id), by = c("gene" = "gene_id")) %>%
    mutate(gene_name = factor(gene_name, levels = unique(ord$gene_name)))
  
      png(paste0("plots/002_7.", names(shrRes)[i], "_counts.png"), h = 2000, w = 4000, res = 600)
      
       p <-  ggplot(best, aes(x = gt, y = count, color = names)) +
          ggbeeswarm::geom_beeswarm(size = 3, cex = 3, alpha = 0.7) +
          facet_wrap(~gene_name, nrow = 2) +
          labs(title = paste0("Best DE genes and differentiation markers"),
               x = "Genotype",
               y = "Normalized counts") +
          scale_x_discrete(labels = function(x) str_replace(x, "_", "\n")) +
          scale_y_log10(expand = expansion(mult = 0.1)) +
          theme(plot.title = element_text(hjust = 0.5))
       print(p)
      
      dev.off()
      best
})



# Plot best scoring genes
source("code/p_star.R")
map(seq(1:length(shrRes)), function(i, ngenes = 10){

  ends <- shrRes[[i]] %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    filter(padj < 0.05) %>%
    filter(str_detect(rowname, "^CG|^CR", negate = T)) %>%    # Keep only known genes
    arrange(log2FoldChange) %>%
    dplyr::slice(-((ngenes+1):(n()-ngenes))) %>%              # Keep head and tail simultaneously (by removing what's in between)
    mutate(Sign = if_else(log2FoldChange > 0, "Up","Down") %>%
             as.factor()) %>%
    left_join(rowRanges(des) %>% 
                as.data.frame() %>% 
                select(description, symbol), by = c("rowname" = "symbol"))

  png(paste0("plots/002_8.", i,"_best_log2FC_genes_", names(shrRes)[i],".png"), h = 3500, w = 4500, res = 600)
  
  p <- ggplot(ends, aes(x = c(1:20), y = log2FoldChange, fill = Sign)) +
    geom_bar(stat = "identity") +
    theme_light() +
    theme(legend.position = "none") +
    labs(title = paste0("Most significant differentially expressed genes in ", names(shrRes)[i]), 
         subtitle = paste0("padj < ", signif(max(ends$padj), 4), " - only known genes are shown"),
         y = "log2FoldChange", 
         x = "") +
    geom_text(aes(label = p_star(padj)), 
              hjust = ifelse(ends$Sign == "Up", -0.3, 1.3),
              vjust = 0.75) +
    theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          axis.title.y = element_text(size = 12),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 10),
          axis.ticks.y = element_blank()) +
    coord_flip() +
    scale_fill_manual(values = c("#00BFC4","#F8766D")) +
    geom_text(aes(label = paste(rowname, "-", description) %>%
                    str_wrap(width = 60),
                  y = ifelse(sign(log2FoldChange) > 0, -0.5, 0.5),
                  hjust = ifelse(sign(log2FoldChange) > 0, 1, 0)),
              position = position_dodge(width = 0.2),
              size = 3.2,
              lineheight = 0.85)
  print(p)
  
  dev.off()
})

# Volcano plots with inset plot
library(EnhancedVolcano)

map(seq(1:length(shrRes)), function(i){
  
  name <- names(shrRes)[i]
  plotres <- shrRes[[i]] %>% 
    as.data.frame() %>% 
    mutate(padj = ifelse(padj < 1e-50, 1e-50, padj),
           pvalue = padj,
           log2FoldChange = ifelse(log2FoldChange < -8, -8, ifelse(
             log2FoldChange > 8, 8, log2FoldChange
           ))) %>% 
    rownames_to_column()
  goi <- c("Antp", "ftz", "nub", "sd", "vg", "pros", "Myc", "Max", "ana", "elav", "retn",
           best_list[[i]] %>% pull(gene_name) %>% unique() %>% as.character())
  
  png(paste0("plots/002_9.", i, "Volcano_plot_", name, "_shr.png"), h = 3500, w = 6500, res = 700)
  
  # labs <- plotres %>%
  #   filter(abs(log2FoldChange) > 7 | padj < 1e-150) %>%
  #   slice_max(order_by = -log10(padj)*abs(log2FoldChange), n = 15) %>% 
  #   pull(rowname)
    gglabs <- ifelse(plotres$rowname %in% goi, plotres$rowname, "")
  
  p <- EnhancedVolcano(plotres, subtitle = "padj cutoff = 0.05 \nlog2FC cutoff = 1",
                       lab = "",
                       x = 'log2FoldChange',
                       y = 'padj',
                       title = paste0('Differential expression in ', name),
                       pCutoff = 0.05, #0.05 cutoff
                       FCcutoff = 1,
                       ylim = c(0, 50),
                       xlim = c(-8, 8),
                       axisLabSize = 14,
                       titleLabSize = 14,
                       subtitleLabSize = 10,
                       captionLabSize = 10,
                       labSize = 4,
                       legendLabSize = 10,
                       legendPosition = "right",
                       drawConnectors = T,
                       maxoverlapsConnectors = Inf,
                       caption = paste0("Downregulated genes = ", 
                                        plotres %>% 
                                          filter(log2FoldChange < -1 & padj < 0.05) %>%
                                          nrow(),
                                        "     ",
                                        "Upregulated genes = ",
                                        plotres %>% 
                                          filter(log2FoldChange > 1 & padj < 0.05) %>% 
                                          nrow()),
                       #selectLab = labs
  ) +
    
    theme(plot.caption = element_text(hjust = 0.5)) +
    geom_text_repel(label = gglabs, 
                    size = 3.5,
                    min.segment.length = 0,
                    max.overlaps = Inf,
                    arrow = arrow(length = unit(0.015, "npc"), ends = "first", type = "closed"),
                    ylim = c(0, 50),
                    #max.time = 5,
                    box.padding = 0.3,
                    point.padding = 0.2,
                    segment.size = 0.3,
                    segment.color = "#00000070",
                    nudge_x = .15,
                    segment.curvature = -1e-20
                    
    )
 
  
  print(p)
  dev.off()

})
  