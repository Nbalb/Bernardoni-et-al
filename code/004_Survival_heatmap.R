library(DESeq2)
library(ComplexHeatmap)
library(RColorBrewer)
library(tidyverse)
library(viridis)

rld <- readRDS("data/002_dds_rlog_normalized.rds")
res <- readRDS("data/002_DESeq2_results.rds")
ort <- readRDS("data/ort_flyrnai_TARGET.rds")        
surv_target <- readRDS("data/Survival_analysis_target.rds") 

# Filter genes that are all increasing or all decreasing in ge in all 3 lists
## Function that takes in a list, coerces to a tibble, 
## creates a column with those names then filters only those genes that are 
## either increasing or decreasing in all 3 contrasts
stouffer <- function(p1, p2){
  st <- (qnorm(p1) + qnorm(p2))/(sqrt(2))
  return(st)
}

filter_lfc_sign <- function(res, sign) {
  map_dfr(names(res), ~ res[[.x]] %>% 
            as_tibble(rownames = "gene_name") %>%
            mutate(contrast = .x)
  ) %>% 
    filter(sign(log2FoldChange) == sign) %>%
    group_by(gene_name) %>% 
    mutate(count = n()) %>% 
    ungroup() %>% 
    filter(count == 3)
}

gn_up <- filter_lfc_sign(res, 1)  # 6831 genes
gn_dn <- filter_lfc_sign(res, -1) # 7206 genes

gn_up_int <- gn_up %>%
  mutate(padj_contr = rep(gn_up %>% 
                            filter(contrast == "lmax_hmyc_vs_neg") %>% 
                            pull(padj), 3),
         stouff_z = stouffer(padj, padj_contr)) %>% 
  filter(!(contrast == "lmax_hmyc_vs_neg"))

gn_dn_int <- gn_dn %>%
  mutate(padj_contr = rep(gn_dn %>% 
                            filter(contrast == "lmax_hmyc_vs_neg") %>% 
                            pull(padj), 3),
         stouff_z = stouffer(padj, padj_contr)) %>% 
  filter(!(contrast == "lmax_hmyc_vs_neg"))

# Take the z-scores of the 2 log2fc, intersect the significant ones 
expratio <- assay(rld) %>% 
  as.data.frame() %>% 
  rownames_to_column("Dm_symbol") %>%
  inner_join(ort, by = c("Dm_symbol" = "FlyBaseID")) %>% 
  relocate(where(is.character))


expsurv <- expratio %>% 
  left_join(surv_target, by = c("Human.Symbol" = "Symbol")) %>% 
  filter(!is.na(Cox_p) & Cox_p < 0.05)

# Plot Heatmap ----
fname <- "data/msigdb_Dm.rds"
if(!file.exists(fname)){
  
  library(msigdbr)
  mdf <- msigdbr(species = "Drosophila melanogaster")
  saveRDS(mdf, fname)

  }else{mdf <- readRDS(fname)}

## Subset for genes present only for GO:BP
mdf %>%
  filter(gs_subcat == "GO:BP") %>% 
  filter(str_detect(gs_name, "NERVOUS|BRAIN")) %>% 
  group_by(gs_name) %>% 
  summarise(n = n()) %>% 
  arrange(-n) %>% 
  print(n = Inf)

ns_genes <- mdf %>% 
  filter(gs_subcat == "GO:BP", str_detect(gs_name, "NERVOUS|BRAIN")) %>% 
  select(gene_symbol, human_gene_symbol) %>% 
  distinct(gene_symbol, human_gene_symbol)

## Plot heatmap
plotdf <- expsurv %>% 
  inner_join(gn_up_int %>% 
               select(gene_name, stouff_z) %>% 
               group_by(gene_name) %>% 
               summarise(mean_z = mean(stouff_z)), by = c("Fly.Symbol" = "gene_name")) %>% 
  slice_min(mean_z*-log10(Cox_padj), n = 50)

plotmat <- plotdf %>% 
  select(matches("lmax|neg", ignore.case = F), Fly.Symbol) %>%
  relocate(matches("neg"), .before = matches("lmax")) %>% 
  column_to_rownames("Fly.Symbol") %>%
  as.matrix()

coxdf <- plotdf %>% 
  filter(Fly.Symbol %in% rownames(plotmat)) %>% 
  select(Cox_coef, Fly.Symbol, Human.Symbol, Cox_padj, Cox_p) %>% 
  mutate(NS = ifelse(Fly.Symbol %in% ns_genes$gene_symbol, "yes", "no"),
         Cox_star = case_when(
           Cox_padj <= 0.001 ~ "<= 0.001",
           Cox_padj > 0.001 & Cox_padj <= 0.01 ~ "<= 0.01",
           Cox_padj > 0.01 & Cox_padj <= 0.05 ~ "<= 0.05",
           Cox_padj > 0.05 ~ "> 0.05",
         )) %>% 
  as.data.frame()

rownames(coxdf) <- paste0(coxdf$Fly.Symbol, " - ", coxdf$Human.Symbol)

htmp_cox <- function(plotmat, coxdf, pathway = NULL) {
  coxdf_plot <- coxdf %>%
    select(Cox_coef, NS, Cox_padj, Cox_star)
  
  Heatmap(
    t(scale(t(plotmat))),
    name = "Row-scaled \ngene expression",
    column_title = "Genotype",
    row_title = "Genes",
    show_row_names = T,
    cluster_columns = F,
    col = viridis::viridis_pal(option = "C")(12),
    
    column_labels = colnames(plotmat) %>%
      str_replace_all("_", " ") %>%
      str_to_title() %>%
      str_wrap(w = 8),
    
    column_names_centered = T,
    column_names_rot = 315
  ) +
    
    Heatmap(
      coxdf_plot[, "Cox_coef"] %>% as.matrix(),
      name = "Prognosis",
      show_row_names = T,
      row_names_gp = gpar(fontsize = 7),
      col = viridis::viridis_pal(option = "D", dir = -1)(12),
      color_space = "LAB",
      width = unit(5, "mm"),
      column_labels = "",
      column_names_centered = T,
      column_names_rot = -45,
      
      # All these parameters can be found in ?Legend()
      heatmap_legend_param = list(
        title = "Prognosis for \nTARGET NBL",
        at = c(-1, 1),
        labels = c("Good", "Bad"),
        legend_height = unit(1.5, "cm")
      )
    ) +
    
    Heatmap(
      coxdf_plot[, "Cox_star"] %>% as.matrix(),
      name = "Cox FDR",
      col = viridis_pal(option = "G")(n_distinct(coxdf_plot[, "Cox_star"])),
      row_names_gp = gpar(fontsize = 7),
      column_labels = "",
      column_names_centered = T,
      column_names_rot = -45,
      width = unit(5, "mm")
      
    ) +
    
    Heatmap(
      coxdf_plot[, "NS", drop = F] %>% as.matrix(),
      name = "Included in a \nNervous\\Brain \nPathway",
      col = c("#D5C2FF", "#580AFF"),
      row_names_gp = gpar(fontsize = 7),
      column_labels = "",
      column_names_centered = T,
      column_names_rot = -45,
      width = unit(5, "mm")
    )
  
  
}

png("plots/004_1.1_expratio_complex_positive_best_fdr.png", w = 5500, h = 4000, res = 600)
draw(htmp_cox(plotmat, coxdf), ht_gap = unit(c(3, 0, 0), "mm"))
dev.off()

# Heatmap opposite expression
plotdf <- expsurv %>% 
  inner_join(gn_dn_int %>% 
               select(gene_name, stouff_z) %>% 
               group_by(gene_name) %>% 
               summarise(mean_z = mean(stouff_z)), by = c("Fly.Symbol" = "gene_name")) %>% 
  slice_min(mean_z*-log10(Cox_padj), n = 50)

plotmat <- plotdf %>% 
  select(matches("lmax|neg", ignore.case = F), Fly.Symbol) %>% 
  relocate(matches("neg"), .before = matches("lmax")) %>% 
  column_to_rownames("Fly.Symbol") %>% 
  as.matrix()

coxdf <- plotdf %>% 
  filter(Fly.Symbol %in% rownames(plotmat)) %>% 
  select(Cox_coef, Fly.Symbol, Human.Symbol, Cox_padj) %>% 
  mutate(NS = ifelse(Fly.Symbol %in% ns_genes$gene_symbol, "yes", "no"),
         Cox_star = case_when(
           Cox_padj <= 0.001 ~ "<= 0.001",
           Cox_padj > 0.001 & Cox_padj <= 0.01 ~ "<= 0.01",
           Cox_padj > 0.01 & Cox_padj <= 0.05 ~ "<= 0.05",
           Cox_padj > 0.05 ~ "> 0.05",
         )) %>% 
  as.data.frame()

rownames(coxdf) <- paste0(coxdf$Fly.Symbol, " - ", coxdf$Human.Symbol)

png("plots/004_1.2_expratio_complex_negative_best.png", w = 5500, h = 4000, res = 600)

draw(htmp_cox(plotmat, coxdf), ht_gap = unit(c(3, 0, 0), "mm"))

dev.off()


# Step survival KM curve
source("code/stepsurvival_cox.R")
expmat <- readRDS("data/target_NBL-expmat.rds")
survival <- readRDS("data/target_NBL-survival.rds")

colnames(expmat) <- gsub(pattern = ".{8}$", replacement =  "", colnames(expmat))

expvar <- apply(expmat, 1, var)
expvar <- names(expvar[expvar > 0.1])

expmat <- expmat[expvar,]

common <- intersect(colnames(expmat),names(survival))
expmat <- expmat[,common]
survival <- survival[common]

mygene <- "UTP3"
oritrack <- expmat[mygene,]

png(paste0("plots/004_2.1_Stepsurvival_", mygene, ".png"), h = 3000, w = 4500, res = 600)
plotstepsurv(oritrack, survival, ngroups = 6, title = paste0("KM survival for ", mygene), cox_survival = TRUE)
dev.off()

mygene <- "DST"
oritrack <- expmat[mygene,]

png(paste0("plots/004_2.2_Stepsurvival_", mygene, ".png"), h = 3000, w = 4500, res = 600)
plotstepsurv(oritrack, survival, ngroups = 6, title = paste0("KM survival for ", mygene), cox_survival = TRUE)
dev.off()
