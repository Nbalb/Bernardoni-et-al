# New enrichment analyses ----
# (https://hbctraining.github.io/Training-modules/DGE-functional-analysis/lessons/functional_analysis_2019.html)
# http://yulab-smu.top/clusterProfiler-book/chapter3.html#wikipathways-analysis
library(tidyverse)
library(fgsea)
library(pathview)
library(clusterProfiler)
library(org.Dm.eg.db)
library(enrichplot)
source("../../../functions/pretty_path_label.R")

shrRes <- readRDS("data/002_DESeq2_shrunk_results.rds")
mdf <- readRDS("data/005_mdf_dm.rds")
mlist <- split(x = mdf$gene_symbol, f = mdf$gs_name)

## Complete analysis working on the list
set.seed(3)
map(seq(1:length(shrRes)), function(i){
  
  df <- shrRes[[i]] %>% as.data.frame() %>% 
    rownames_to_column("Gene_Symbols") %>% 
    mutate(padj = replace(padj, padj == 0, 2.225074e-308))
  
  entrezid <- bitr(df$Gene_Symbols, fromType = "SYMBOL",toType = "ENTREZID", OrgDb = "org.Dm.eg.db")
  
  df_id <- left_join(df, entrezid, by = c("Gene_Symbols" = "SYMBOL"))
  
  ### ClusterProfiler
  
  ### Msigdb
  #### GSEA 
  sig <- setNames(as.numeric(df_id$log2FoldChange), as.character(df_id$ENTREZID))
  sig <- sort(sig, decreasing = T)
  
  nmdf <- mdf %>%
    dplyr::select(gs_name, entrez_gene) %>% 
    dplyr::rename(term = gs_name, gene = entrez_gene) %>% 
    as.data.frame()
  
  gsemdf <- GSEA(sig, TERM2GENE = nmdf)
  
  gsemdf@result <- gsemdf@result %>% 
    mutate(Description = pretty_path_label(Description, 50))
  
  #### ORA (Over-Representation Analysis)
  sig_genes <- df_id %>% 
    filter(padj < 0.05) %>% 
    pull(ENTREZID)
  
  all_genes <- pull(df_id, ENTREZID)
  
  emdf <- enricher(gene = sig_genes,
                   universe = all_genes,
                   TERM2GENE = nmdf)
  
  emdf@result <- emdf@result %>% 
    mutate(Description = pretty_path_label(Description, 50))
  
  ### Plots
  p1 <- dotplot(gsemdf, showCategory = 15) + ggtitle("dotplot for GSEA")
  p2 <- dotplot(emdf, showCategory = 15) + ggtitle("dotplot for ORA")
  
  png(paste0("plots/005b_1_", names(shrRes)[i],"Dotplot_ORA_vs_GSEA_msigdb.png"), h = 5000, w = 8500, res = 600)
  p3 <- cowplot::plot_grid(p1, p2, ncol = 2)
  print(p3)
  dev.off()
  
  # heatplot(gsemdf, foldChange = sig, showCategory = 10)
  # 
  # emdf_emap <- pairwise_termsim(emdf)
  # 
  # emapplot(emdf_emap)
  
  ### GSEA WP (only 2 results)
  # gsewp <- gseWP(sig, organism = "Drosophila melanogaster")
  # gsewp@result <- gsewp@result %>% 
  #   mutate(Description = pretty_path_label(Description))
  # 
  # dotplot(gsewp, showCategory = 30)
  
  ### EnrichGO
  res <- shrRes[[i]] %>% 
    as.data.frame()
  
  all_genes <- as.character(rownames(res))
  sig_res <- res %>% 
    dplyr::filter(!is.na(padj) & padj < 0.05)
  sig_genes <- as.character(rownames(sig_res))
  
  ego <- enrichGO(gene = sig_genes, 
                  universe = all_genes,
                  keyType = "SYMBOL",
                  OrgDb = org.Dm.eg.db, 
                  ont = "BP", 
                  pAdjustMethod = "BH", 
                  qvalueCutoff = 0.05)
  
  egos <- clusterProfiler::simplify(ego, cutoff = 0.65)
  
  ###Plots
  # dotplot(ego, showCategory = 20)
  # dotplot(egos, showCategory = 20)
  # 
  # fcs <- sig_res$log2FoldChange
  # names(fcs) <- rownames(sig_res)
  # 
  # plotegos <- egos
  # plotegos@result <- egos@result[c(3,6,7,10),]
  
  # cnetplot(plotegos, 
  #          categorySize = "pvalue", 
  #          showCategory = 5, 
  #          foldChange = fcs, 
  #          vertex.label.font = 3,
  #          colorEdge = T,
  #          node_label = "category",
  #          color_category = "black") +
  #   scale_color_gradientn(name = "fold change",
  #                         colours = c("blue","white","red"), 
  #                         values = scales::rescale(c(-1, 0, 1)),
  #                         limits = c(-1, 1),
  #                         breaks = c(-1 , 0, 1))
  # 
  # goplot(egos)
  
  
  png(paste0("plots/005b_2_", names(shrRes)[i],"_Emapplot_ORA_GO_terms.png"), h = 4000, w = 5000, res = 600)
  p <- emapplot(pairwise_termsim(egos), 
           showCategory = 15,
           layout = "graphopt",
           cex_label_category = 0.8,
           repel = T,
           label_format = 10,
           cex_pie2axis = 0.8)
  
  print(p)
  dev.off()
  
  gsemdf_long <- gsemdf
  gsemdf_long@result <- gsemdf_long@result %>% 
    mutate(Description = str_replace_all(Description, "\n", " ")) %>% 
    mutate(Description = str_replace_all(Description, "Endoplasmic Reticulum", "ER"))
  
  png(paste0("plots/005b_3_", names(shrRes)[i],"_Upsetplot_GSEA_GO_terms.png"), h = 4000, w = 8000, res = 600)
  p <- upsetplot(gsemdf_long) +
    ylab("log2FoldChange")
  print(p)
  dev.off()
  
  emdf_long <- emdf
  emdf_long@result <- emdf_long@result %>% 
    mutate(Description = str_replace_all(Description, "\n", " ")) %>% 
    mutate(Description = str_replace_all(Description, "Endoplasmic Reticulum", "ER"))
  
  png(paste0("plots/005b_4_", names(shrRes)[i],"_Upsetplot_ORA_GO_terms.png"), h = 4000, w = 8000, res = 600)
  p <- upsetplot(emdf_long) +
    ylab("Nunmber of genes")
  print(p)
  dev.off()
  
  ### Ridgeplot
  png(paste0("plots/005b_5_", names(shrRes)[i],"_Ridgeplot_msigdb_gsea.png"), h = 6000, w = 9000, res = 600)
  p <- ridgeplot(gsemdf, showCategory = 20) +
    ggtitle(paste0("GSEA results distribution for ", names(shrRes[i])))
  print(p)
  dev.off()
  
  ### GSEA KEGG
  # sig_cg <- sig
  # cgid <- bitr(names(sig), fromType = "ENTREZID", toType = "FLYBASECG", OrgDb = "org.Dm.eg.db")
  # cgid <- map_chr(cgid$FLYBASECG, ~paste0("Dmel_", .))
  # names(sig_cg) <- cgid
  # 
  # gseaKEGG <- gseKEGG(geneList = sig_cg,
  #                     organism = "dme",
  #                     keyType = "kegg",
  #                     minGSSize = 20
  # )
  # 
  # gseaKEGG@result$Description
  # gseaKEGG@result$ID
  # 
  # pathview(gene.data = sig,
  #          pathway.id = "dme03030",
  #          species = "dme",
  #          limit = list(gene = 2,
  #                       cpd = 1),
  #          kegg.dir = "./plots")
  
})
