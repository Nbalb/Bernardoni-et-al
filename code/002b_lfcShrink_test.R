dresults <- readRDS("../../mycmax/results/001_results_deseqobj.rds")   #check script 001_DE.R in mycmax/results folder
class(dresults[[1]])    #"DESeqResults"
des <- readRDS("../../mycmax/results/001_Deseq_obj.rds")
resultsNames(des)   #kept same names in the newer script
#[1] "Intercept"          "gt_lmax_vs_neg"     "gt_lmaxhmyc_vs_neg"
coef1 <- resultsNames(des)[2]
freshr <- lfcShrink(des, coef = coef1, res = dresults[[1]])

png("plots/002b_1_MAplot_unshrunk_fede.png", w = 5000, h = 3000, res = 600)
plotMA(dresults[[1]],
       main = paste0("MAplot ", names(shrRes)[1]),
       colSig = "#1E20FF80",
       colNonSig = "#88888870",
       colLine = "#ff000080")
dev.off()

png("plots/002b_2_MAplot_shrunk_fede.png", w = 5000, h = 3000, res = 600)
plotMA(freshr,
       main = paste0("MAplot ", names(shrRes)[1], " after shrinkage"),
       colSig = "#1E20FF80",
       colNonSig = "#88888870",
       colLine = "#ff000080")
dev.off()

# Volcano Plot
library(EnhancedVolcano) 
library(airway)
library(magrittr)


  ###### Enhanced Volcano
  
  res = dresults[[1]] #nrow = 55291
  res = res[!is.na(res$log2FoldChange) & !is.na(res$padj),] #nrow = 18441
  res$padj = res$padj+1e-300
  upreg = res[res$log2FoldChange > 1 & res$padj < 0.01,]
  downreg = res[res$log2FoldChange < -1 & res$padj < 0.01,]
  
  png(paste0("plots/Volcano_", names(dresults)[[1]], ".png"),width=6000,height=3500,res=600,pointsize=40)
  EnhancedVolcano(res, subtitle = "pvalue cutoff = 0.01 \nlog2FC cutoff = 1",
                       lab = rownames(res),
                       selectLab = "FBgn0001077",
                       x = 'log2FoldChange',
                       y = 'padj',
                       # xlim = c(-6, 6),
                       # ylim = c(0,5),
                       title = paste0('Differential expression in ', names(dresults[[1]])),
                       pCutoff = 0.01, #0.05 cutoff
                       FCcutoff = 1,
                       axisLabSize = 14,
                       titleLabSize = 14,
                       subtitleLabSize = 10,
                       captionLabSize = 10,
                       pointSize = 2.5,
                       labSize = 3.5,
                       legendLabSize = 10,
                       drawConnectors = T,
                       caption = paste0("Downregulated genes = ",nrow(downreg), "     ","Upregulated genes = ",nrow(upreg))
  )+
    theme(plot.title = element_text(hjust = 0.5), 
          plot.subtitle = element_text(hjust = 0.5),
          plot.caption = element_text(hjust = 0.5))
  
  dev.off()
  
# Plot shrunk
  res = freshr #nrow = 55291
  res = res[!is.na(res$log2FoldChange) & !is.na(res$padj),] #nrow = 18441
  res$padj = res$padj+1e-300
  upreg = res[res$log2FoldChange > 1 & res$padj < 0.01,]
  downreg = res[res$log2FoldChange < -1 & res$padj < 0.01,]
  
  png(paste0("plots/Volcano_", names(dresults)[[1]], "_shrunk.png"),width=6000,height=3500,res=600,pointsize=40)
  EnhancedVolcano(res, subtitle = "pvalue cutoff = 0.01 \nlog2FC cutoff = 1",
                  lab = rownames(res),
                  selectLab = "FBgn0001077",
                  x = 'log2FoldChange',
                  y = 'padj',
                  # xlim = c(-6, 6),
                  # ylim = c(0,5),
                  title = paste0('Differential expression in ', names(dresults[[1]])),
                  pCutoff = 0.01, #0.05 cutoff
                  FCcutoff = 1,
                  axisLabSize = 14,
                  titleLabSize = 14,
                  subtitleLabSize = 10,
                  captionLabSize = 10,
                  pointSize = 2.5,
                  labSize = 3.5,
                  legendLabSize = 10,
                  drawConnectors = T,
                  caption = paste0("Downregulated genes = ",nrow(downreg), "     ","Upregulated genes = ",nrow(upreg))
  )+
    theme(plot.title = element_text(hjust = 0.5), 
          plot.subtitle = element_text(hjust = 0.5),
          plot.caption = element_text(hjust = 0.5))
  
  dev.off()
  