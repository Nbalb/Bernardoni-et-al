library(DESeq2)
library(tidyverse)
source("../../functions/plotEnrichment2.R")
source("../../functions/pretty_path_label.R")


# Following msigdbr and fgsea vignettes ----
fname <- "data/005_mdf_dm.rds"
if(!file.exists(fname)){
  library(msigdbr)  
  msigdbr_species()
  msigdbr_collections() %>% print(n = Inf)
  mdf2 <- msigdbr(species = "Drosophila melanogaster") %>%     # Retrieve all Dm gene sets
    filter(gs_cat %in% c("C1", "C2", "C5", "C7", "H")) %>% 
    filter(gs_subcat %in% c("", "GO:MF", "GO:BP", "GO:CC", "CP:KEGG", "CP:WIKIPATHWAYS", "CP:BIOCARTA"))
  saveRDS(mdf, file = fname)
}else(mdf <- readRDS(fname))

mlist <- split(x = mdf$gene_symbol, f = mdf$gs_name)

# Get signatures
res <- readRDS("data/002_DESeq2_shrunk_results.rds")

fname <- "data/005_fgsea_results.rds"
set.seed(3)

if(!file.exists(fname)){
  
  library(fgsea)
  library(data.table)
  
  fgseaRes <- list()
  sig <- list()
  
  for(i in 1:length(res)){
    
    df <- res[[i]] %>% as.data.frame() %>% 
      mutate(padj = replace(padj, padj == 0,  2.225074e-308)) %>% 
      filter(!is.na(padj) & !is.na(log2FoldChange))
      
    signature <- setNames(-log10(df$padj)*sign(df$log2FoldChange), rownames(df))
    sig[[i]] <- signature
    
    cat("Running fgsea for", names(res)[i], "\n")
    #Use fgseaMultilevel for better accuracy than fgseaSimple (https://bioconductor.riken.jp/packages/3.9/bioc/vignettes/fgsea/inst/doc/fgseaMultilevel-tutorial.html)
    fgseaRes[[i]] <- fgseaMultilevel(pathways = mlist,
                                stats = sig[[i]],
                                eps = 0,
                                nproc = parallel::detectCores()-2,
                                )  
    names(fgseaRes)[i] <- names(res)[i]
    
  }
   
  saveRDS(fgseaRes, fname)
  
}else{fgseaRes <- readRDS(fname)}

# Plots ----
des <- readRDS("data/002_des.rds")
map(seq(1:length(fgseaRes)), function(i){
   
  contr <- names(fgseaRes)[i]
  collapsedPathways <- collapsePathways(fgseaRes[[i]][order(pval)][padj < 0.01], 
                                        mlist, sig[[i]])
  
  # mainPathways <- fgseaRes[[i]][pathway %in% collapsedPathways$mainPathways][
  #    order(-NES), pathway]
  
  mainPathways <- fgseaRes[[i]] %>% 
    filter(pathway %in% collapsedPathways$mainPathways) %>% 
    arrange(-NES) %>% 
    pull(pathway)
  
  png(paste0("plots/005_01.",i,"_", contr, "_top20_enriched_pathways.png"), h = 5000, w = 6500, res = 600)
  plotGseaTable(mlist[mainPathways[1:20]], sig[[i]], fgseaRes[[i]], 
                gseaParam = 0.5, colwidths = c(12,4,1.8,2,2))
  dev.off()
  
  ## Plot single GSEAs
  pdf(paste0("plots/005_02.",i,"_", contr, "_main_pathways_enrichment.pdf"), w = 5, h = 3)
  for(o in 1:length(mainPathways)){
    
    NES <- signif(fgseaRes[[i]][pathway == mainPathways[o]]$NES, 3)
    padj <- signif(fgseaRes[[i]][pathway == mainPathways[o]]$padj, 3)
    title <- paste0("GSEA ", 
                    contr, 
                    " in \n", 
                    pretty_path_label(mainPathways[o]) %>% str_wrap(50))
    
    p <- plotEnrichment2(mlist[[mainPathways[o]]], stats = sig[[i]]) +
      labs(title = title, 
              subtitle = paste0("NES = ", NES, "  p.adj = ", padj)) +
      theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
            plot.subtitle = element_text(size = 8, hjust = 0.5)) +
      xlim(0, 10500)
    print(p)
  }
  dev.off()
})

## Barplot gsea
map(seq(1:length(fgseaRes)), function(i){
  
  contr <- names(fgseaRes)[i]
  collapsedPathways <- collapsePathways(fgseaRes[[i]][order(pval)][padj < 0.01], 
                                        mlist, sig[[i]])
  
  # mainPathways <- fgseaRes[[i]][pathway %in% collapsedPathways$mainPathways][
  #    order(-NES), pathway]
  
  mainPathways <- fgseaRes[[i]] %>% 
    filter(pathway %in% collapsedPathways$mainPathways) %>% 
    arrange(-NES) %>% 
    pull(pathway)
  
  topPathwaysUp <- fgseaRes[[i]][ES > 0][pathway %in% mainPathways][head(order(pval), n=10), pathway]
  topPathwaysDown <- fgseaRes[[i]][ES < 0][pathway %in% mainPathways][head(order(pval), n=10), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  
  ends <- fgseaRes[[i]] %>%
    filter(pathway %in% topPathways) %>%
    slice_max(abs(NES), n = floor(nrow(.)/2)*2) %>%       # New way to slice ends, order by absolute value of NES, then keep an equal number of rows from both ends
    arrange(NES) %>%
    mutate(pathlab = pretty_path_label(pathway) %>% str_wrap(width = 50))
  
  ends$p_star <- sapply(ends$padj, 
                        function(x) {
                          if(x <= 10^(-4)){"****"
                          }else if(x <= 10^(-3)){"***"
                          }else if(x <= 10^(-2)){"**"
                          }else if(x <= 5*10^(-2)){"*"
                          }else{""}
                        })  
  
  png(paste0("plots/005_03.",i,"_", contr,"GSEA_barplot.png"), h = 3500, w = 4500, res = 600)
  
  plt <- ggplot(ends, aes(x = c(1:nrow(ends)), y = NES, fill = as.factor(sign(-NES)))) +
    geom_bar(stat='identity') +
    theme_light() +
    theme(legend.position = "none") +
    labs(title = paste0("Best scoring pathways for ", contr), y = "Combined Score", x="") +
    theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
          axis.title.y = element_text(size = 12),
          axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 10)) +
    coord_flip() +
    geom_text(aes(label = pathlab, 
                  y = ifelse(NES < 0, 0.05, -0.05),
                  hjust = ifelse(NES < 0, 0, 1)),
              position = position_dodge(width = 0),
              size = 3,
              lineheight = 0.85) +
    geom_text(aes(label = p_star), 
            hjust = ifelse(sign(ends$NES) < 0, 1.3, -0.3),
            vjust = 0.75)
  
  
  print(plt)
  dev.off()

})
# Save in an excel format all the significant pathways
library(openxlsx)
fres <- map(fgseaRes, function(x){
  
  x %>% 
    filter(padj < 0.05 & size > 10) %>% 
    arrange(padj)
  
})
write.xlsx(fres, file = "results/Dm_mycmax_pathway_analysis.xlsx")


# Plot LeadingEdge counts and related pathway
## Take only pathways large and significant, extract all the leading edges from all the fgsea contrasts, count which are the most recurrent genes (Spoiler: Rp, Ribosomal proteins - removed them since not informative enough) and select the top 20 
top_le <- fres %>% 
  map(select, leadingEdge) %>% 
  unlist(use.names = F) %>% 
  as_tibble() %>% 
  dplyr::count(value, sort = T) %>% 
  filter(str_detect(value, "^Rp", negate = T)) %>% 
  slice_head(n = 30)    # Count how many times a gene is found in best significant pathways


top_le_count <- map_dfr(top_le$value, function(x){
  
  geneCounts <- plotCounts(dds = des,
                           gene = rowRanges(des) %>%
                             as.data.frame() %>%
                             rownames_to_column("rowname") %>%
                             filter(gene_name == x) %>%
                             pull(rowname),
                           intgroup = c("gt", "names"),
                           returnData = T) %>%
    mutate(gene = x) 
  
}) %>%
  mutate(gene = factor(gene, levels = unique(top_le$value)))

# map(fres, pull, leadingEdge, pathway) %>% flatten()

# Leading edge analysis
# Create a list of all the significant pathways with their LE as object (for all contrasts)
# Flatten removes a level from a list, discard keeps all the FALSE elements from a logical test, this way we remove duplicated pathways (which should not have been there)
s_res <- map(fres, pull, leadingEdge, pathway) %>% 
  flatten() %>% 
  discard(map(fres, pull, leadingEdge, pathway) %>% 
            flatten() %>% 
            duplicated())


pathway_x <- vector(mode = "list", length = nrow(top_le))   # Preallocate list space, names of the list are the genes, elements of the list are vectors containing all the significant pathways the gene is found in
names(pathway_x) <- top_le$value

for(i in 1:length(top_le$value)){
  
  pathway_x[[i]] <- names(s_res)[map_lgl(s_res, ~ any(str_detect(., top_le$value[i])))]  # search in all the elements of the list with LEs if there is one of the genes of interest. str_detect converts all the elements in the list as logicals, any transforms it to a single value (gene present = T), map_lgl transforms the list to a named (with pathway name) logical vector that is then used to extract all the pathways that have T as logical from the vector with all the pathways names. These pathways are saved to a list
  
}

plist <- vector(mode = "list", length = length(unique(top_le_count$gene)))

# top_le_count = all the counts for the selected genes
# pathway_x = list of all the pathways in which each gene is involved, names of the list are the genes
# Every loop I create a dataset for a gene where NES from each fgsea contrast are extracted, then group by contrast and remove duplicated pathways (which I don't understand why they duplicate), then create tot_NES in order to find pathways with largest abs(NES)

for (i in unique(top_le_count$gene)) {
  
  genedf <- top_le_count %>% 
    filter(gene == i)
  
  assoc_path <- map_dfr(pathway_x[[i]], function(x){
    
    df <- data.frame(
      pathway = x,
      lmaxVsNeg = fgseaRes[[1]] %>% 
        filter(pathway == x) %>% 
        pull(NES),
      lmaxhmycVsNeg = fgseaRes[[2]] %>% 
        filter(pathway == x) %>% 
        pull(NES),
      lmaxhmycVslmax = fgseaRes[[3]] %>% 
        filter(pathway == x) %>% 
        pull(NES)
    )
      
      df
  }) %>% gather(`lmaxVsNeg`, `lmaxhmycVsNeg`, `lmaxhmycVslmax`, key = "contrast", value = "NES") %>% 
    mutate(contrast = factor(contrast, levels = names(fres)),
           pathway = pretty_path_label(pathway)) %>%
    group_by(contrast) %>% 
    distinct(pathway, .keep_all = T) %>% 
    ungroup() %>% 
    group_by(pathway) %>% 
    mutate(tot_NES = sum(abs(NES))) %>% 
    ungroup() %>% 
    slice_max(order_by = tot_NES, n = 72)
    
    
  p1 <- ggplot(genedf, aes(x = gt, y = count, color = names)) +
    scale_y_log10() +
    ggbeeswarm::geom_beeswarm(size = 5, cex = 3, alpha = 0.7) +
    labs(title = paste0(i, " expression"),
         x = "Genotype",
         y = parse(text=paste("log[10]","~normalized ~counts"))) +
    scale_x_discrete(labels = function(x) str_replace(x, "_", "\n")) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  
  p2 <- ggplot(assoc_path, aes(x = reorder(pathway, NES), y = NES, fill = as.factor(sign(NES)))) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = 0, color = "lightgrey", size = 1.1) +
    theme_light() +
    theme(legend.position = "none") +
    labs(title = paste0("Best scoring pathways for ", i), y = "NES", x="") +
    facet_wrap(~ contrast, ncol = 3) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          axis.title.x = element_text(size = 12),
          axis.text.y = element_text(size = 9, lineheight = 0.8),
          axis.text.x = element_text(size = 10)) +
    coord_flip() +
    scale_x_discrete(position = "top",
                     labels = function(x) str_wrap(x, width = 55)) +
    theme(strip.background = element_rect(fill = "white")) +
    theme(strip.text = element_text(colour = "black", face = "bold", size = 10)) + 
    scale_fill_manual(values = c("#00BFC4", "#F8766D"))
  
  p3 <- cowplot::plot_grid(p1, p2, rel_widths = c(1,2))
  
  plist[[which(i == unique(top_le_count$gene))]] <- p3
}

ggsave("005_4_LE_analysis.pdf",
       plot = gridExtra::marrangeGrob(plist, nrow=1, ncol=1, top = NULL),
       device = "pdf",
       path = "plots",
       width = 8500,
       height = 3500,
       units = "px",
       dpi = 600)
