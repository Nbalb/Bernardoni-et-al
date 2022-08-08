library(DESeq2)
library(tidyverse)
getwd()   # D:/OneDrive - Alma Mater Studiorum Università di Bologna/PhD/projects/mycmax2
dir <- "../../../../linux/mycmax/salmon_alignments"

files <- file.path(list.files(dir, "quant_orig_", full.names = T), "quant.sf")
file.exists(files)

quants <- lapply(files, read_delim, delim = "\t") %>%
  setNames(str_extract(files, "RNASEQ[1-3][A-C]_S[0-9]{1,2}"))    # Salmon transcritps

gse <- readRDS("data/002_gse.rds")    # Salmon genes after importing with tximeta and summarizeToGene

newNames <- data.frame("RNASEQ" = gsub("_quant_orig_index", "", list.files(dir, "quant_orig")),
                       "genotype" = c("lmax_A","lmax_B","lmax_C","lmax_hmyc_A","lmax_hmyc_B","lmax_hmyc_C","neg_A","neg_B","neg_C"))

fcc <- read_delim("../../../../linux/mycmax/bams/fc_counts/featureCounts_gene.txt",
                  delim = "\t",
                  skip = 1) %>% 
  select(starts_with("RNASEQ"), Geneid) %>% 
  rename_with(~(newNames %>% pull(genotype)), starts_with("RNASEQ"))

# Plot common genes ----
common <- intersect(rownames(assay(gse)), fcc$Geneid)
salmon <- assay(gse)[common,]
fc <- as.data.frame(fcc)
rownames(fc) <- fcc$Geneid
fc <- fc[common,]

png("plots/003_1_Comparison_common_genes.png", h = 4000, w = 6000, res = 600)
par(mfrow = c(3,3))
for (i in 1:ncol(assay(gse))) {
  
  x <- salmon[,i]
  y <- setNames(fc[,i], common)
  
  plot(x = salmon[,i],
       y = fc[,i],
       main = colnames(fc)[i],
       pch = 16,
       cex = 1.5,
       col = "#66666660",
       xlab = "Salmon counts",
       ylab = "FeatureCounts counts")
  
  points(x = x["FBgn0001077"],
         y = y["FBgn0001077"],
         col = "#f94144",
         pch = 16)

}

dev.off()

# Single plot common ----

png("plots/003_2_Comparison_common_genes_lmax_A.png", h = 3000, w = 3500 , res = 600)

x <- salmon[,"lmax_A"]
y <- setNames(fc[,"lmax_A"], common)

plot(x = salmon[,"lmax_A"],
     y = fc[,"lmax_A"],
     main = "Salmon vs featureCounts common genes' counts",
     xlim = c(0, 3e4),
     ylim = c(0, 3e4),
     pch = 16,
     cex = 1.5,
     col = "#66666660",
     xlab = "Salmon counts",
     ylab = "FeatureCounts counts")

points(x = x["FBgn0001077"],
       y = y["FBgn0001077"],
       col = "#f94144",
       pch = 16)

text(x = x["FBgn0001077"],
     y = y["FBgn0001077"]+1200,
     "ftz",
     font = 2)

dev.off()


# Barplot total counts with different minOverlap vs ftz ----
dir <- "../../../../linux/mycmax/bams/fc_counts"
countsFiles <- list.files(dir, "minOv[0-9]{1,2}.txt.summary", full.names = T)

totCounts <- map(countsFiles, function(x) {
  
  read_delim(x, delim = "\t") %>% 
    mutate(minOv = str_extract(x, "minOv[0-9]{1,2}") %>% 
             str_extract("[0-9]{1,2}") %>% 
             as.numeric(), .after = Status)
  
}) %>% 
  map_dfr(dplyr::filter, Status == "Assigned") %>% 
  arrange(minOv) %>% 
  rename_with(~str_replace(., ".sorted.bam", ""), starts_with("RNASEQ")) %>% 
  gather(starts_with("RNASEQ"), key = "Sample", value = "totCounts") %>% 
  select(-Status)


ftzFiles <- str_replace(countsFiles, ".summary", "")

ftzCounts <- map(ftzFiles, function(x){
  
  read_delim(x, delim = "\t", skip = 1) %>% 
    select(Geneid, starts_with("RNASEQ")) %>% 
    mutate(minOv = str_extract(x, "minOv[0-9]{1,2}") %>% 
             str_extract("[0-9]{1,2}") %>% 
             as.numeric(), .after = Geneid)
  
}) %>% 
  map_dfr(dplyr::filter, Geneid == "FBgn0001077") %>% 
  arrange(as.numeric(minOv)) %>% 
  rename_with(~str_replace(., ".sorted.bam", ""), starts_with("RNASEQ")) %>% 
  gather(starts_with("RNASEQ"), key = "Sample", value = "ftzCounts") %>% 
  select(-Geneid) 

finalCounts <- totCounts %>% 
  left_join(ftzCounts, by = c("minOv", "Sample")) %>% 
  mutate(minOv = factor(minOv))

coeff <- mean(finalCounts %>% dplyr::filter(ftzCounts > 1) %>% pull())/mean(finalCounts$totCounts)

png("plots/003_3_ftz_vs_total_ggplotmore.png", h = 5000, w = 6500, res = 600)
ggplot(finalCounts, aes(x = minOv, fill = minOv)) +
  geom_bar(aes(y = totCounts), stat = "identity", colour = "black") +
  geom_line(aes(y = ftzCounts/coeff, x = minOv, group = 1), size = 1.2, colour = "#F71735") +
  geom_text(aes(y = totCounts, label = scales::scientific(totCounts, 4) %>% str_replace("e","\ne")), 
            # position = position_dodge(width=0.9), 
            vjust = -0.25,
            size = 3) +
  scale_y_continuous(name = "Total Counts",
                     #expand = expansion(mult = c(0, 0)),
                     labels = scales::scientific,
                     breaks = scales::pretty_breaks(3),
                     sec.axis = sec_axis(trans = ~.*coeff,
                                         name = "ftz Counts")) +
  facet_wrap(~Sample, ncol = 3) +
  scale_fill_brewer(palette = "Blues") +
  labs(title = "Effect of different minOverlap values in featureCounts",
       subtitle = "ftz counts vs total gene counts") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.title.y.left = element_text(color = "#0247A1"),
        axis.text.y.left = element_text(color = "#0247A1"),
        axis.title.y.right = element_text(color = "#F71735"),
        axis.text.y.right = element_text(color = "#F71735")
        )
dev.off()

# 2 y axes plot in R base (wrong scales) ----
png("plots/003_3_ftz_vs_total_base.png", h = 4500, w = 6000, res = 600)
par(mar = c(5.1, 4.1, 4.1, 5.1), mfrow = c(3,3))
for (samp in unique(finalCounts$Sample)) {
  
  plotData <- finalCounts %>% 
    dplyr::filter(Sample == samp) %>% 
    mutate(minOv = as.numeric(levels(minOv)))
  
  barplot(plotData$totCounts,
          names = plotData$minOv,
          col = RColorBrewer::brewer.pal(5, "Greens"),
          xlab = "Minimum Overlap",
          ylab = "Total Counts")
  par(new = T)
  plot(x = plotData$minOv,
       y = plotData$ftzCounts,
       ylim = c(0, 50),
       main = samp,
       type = "l",
       lwd = 2,
       axes = F,
       xlab = "",
       ylab = ""
       )
  axis(side = 4, at = seq(0,50,10))
  mtext("ftzCounts", side = 4, line = 3, cex = 0.7)
  
}
dev.off()

# log2fc comparison ----
shrRes <- readRDS("data/002_DESeq2_shrunk_results.rds")
load("../../mycmax/results/001_results.rda")   #dresults
names(dresults)   # "lmax_vs_neg"      "lmaxhmyc_vs_neg"  "lmaxhmyc_vs_lmax"
names(shrRes)     #"lmaxVsNeg"      "lmaxhmycVsNeg"  "lmaxhmycVslmax"

names(dresults) <- names(shrRes)

for(contr in names(shrRes)){
  
  splot <- shrRes[[contr]] %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    left_join(dresults[[contr]] %>%
                as.data.frame() %>%
                select(Gene_Symbol, log2FoldChange, padj), 
              by = c("rowname" = "Gene_Symbol"), suffix = c(".new", ".old")) %>%
    filter(!is.na(log2FoldChange.new) & !is.na(log2FoldChange.old))  %>%
    mutate(log2FCdiff = log2FoldChange.new - log2FoldChange.old)
  
  png(paste("plots/003_4_log2FC_comparison_old_new_", contr, "_vir.png"), h = 2000, w = 3500, res = 600)
  p <- ggplot(splot, aes(x = log2FoldChange.new, y = log2FoldChange.old)) +
    ggpointdensity::geom_pointdensity(size = 2, alpha = 0.8) +
    labs(title = paste("Log2FC comparison new vs old pipeline in", contr),
           subtitle = paste0("Number of common significant genes: ", nrow(splot), "\n",
             "Number of common significant genes with |log2FC difference| > 1: ", 
             sum(abs(splot$log2FCdiff) > 1), 
             " (", signif((sum(abs(splot$log2FCdiff) > 1))/(nrow(splot))*100, 3), " %)")) +
    scale_color_viridis_c() +
    theme(plot.subtitle = element_text(size = 9))
  print(p)  
  dev.off()
  
  
  
}

# Significant genes comparison
tot <- map(seq(1:length(shrRes)), function(i){
  
  as.data.frame(shrRes[[i]]) %>% 
    filter(!is.na(padj)) %>% 
    summarise(total = n(),
              signif = sum(padj < 0.05)) %>% 
    mutate(contrast = names(shrRes)[i])
  
}) %>% 
  bind_rows() %>% 
  mutate(pipeline = "New") %>% 
  bind_rows(map(seq(1:length(dresults)), function(i){
    
    as.data.frame(dresults[[i]]) %>% 
      filter(!is.na(padj)) %>% 
      summarise(total = n(),
                signif = sum(padj < 0.05)) %>% 
      mutate(contrast = names(dresults)[i])
    
  }) %>% 
    bind_rows() %>% 
    mutate(pipeline = "Old"))
  
png("plots/003_5_Signifcant_genes_comparison.png", h = 2500, w = 5000, res = 600)

p1 <- ggplot(tot, aes(x = contrast, y = total, fill = pipeline)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black") +
  geom_text(aes(y = total, label = total), 
            position = position_dodge(0.9),
            vjust = -0.4,
            size = 3) +
  labs(title = "Total measured genes New vs Old pipeline",
       fill = "Pipeline") +
  scale_y_continuous(name = "Number of total genes") +
  scale_x_discrete(name = "Contrast") +
  theme(plot.title =  element_text(hjust = 0.5, size = 10),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 7))


p2 <- ggplot(tot, aes(x = contrast, y = signif, fill = pipeline)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black") +
  geom_text(aes(y = signif, label = signif), 
            position = position_dodge(0.9),
            vjust = -0.4,
            size = 3) +
  labs(title = "Signifcant genes New vs Old pipeline",
       subtitle = "padj < 0.05",
       fill = "Pipeline") +
  scale_y_continuous(name = "Number of significant genes") +
  scale_x_discrete(name = "Contrast") +
  theme(plot.title =  element_text(hjust = 0.5, size = 10),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 7))

cowplot::plot_grid(p1, p2, ncol = 2)

dev.off()


