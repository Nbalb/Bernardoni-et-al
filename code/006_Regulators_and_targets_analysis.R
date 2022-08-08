library(DESeq2)
library(tidyverse)

# Find target of Antp, Nub, Vg, Sd
res <- readRDS("data/002_DESeq2_results.rds")
res_df <- map_dfr(names(res), function(x) {
  
  res[[x]] %>% 
    as_tibble(rownames = "gene") %>% 
    mutate(contrast = x) %>% 
    filter(!(contrast == "lmax_hmyc_vs_lmax")) %>% 
    rename(fc = log2FoldChange) %>% 
    select(gene, contrast, fc, padj)

  }) 

# Get interaction networks from FlyBase (http://ftp.flybase.org/releases/FB2022_04/precomputed_files/genes/gene_genetic_interactions_fb_2022_04.tsv.gz)
gint <- read_delim("data/gene_genetic_interactions_fb_2022_03.tsv.gz",
                        comment = "##", 
                   col_names = c("target_gn", 
                                 "target_id", 
                                 "regulator_gn", 
                                 "regulator_id",
                                 "Interaction_type", 
                                 "Publication_FBrf"))

gint_int <- map_dfr(c("Antp", "nub", "vg", "sd"), function(x){
  
  gint %>% 
    filter(regulator_gn == x) %>% 
    mutate(primary = x,
           secondary = target_gn,
           primary_role = "regulator") %>% 
    bind_rows(gint %>% 
                filter(str_detect(target_gn, x)) %>% 
                mutate(primary = x,
                       secondary = regulator_gn,
                       primary_role = "target")) %>% 
    select(primary, secondary, primary_role , Interaction_type)
  
})

# Plot distributions of the interactions of selected genes
gint_int %>% 
  group_by(primary) %>% 
  count(primary_role) %>% 
  ungroup() %>% 
  mutate(primary_role = as_factor(primary_role)) %>% 
  ggplot(aes(primary, n, fill = primary_role)) +
  geom_bar(stat = "identity", 
           position = "dodge", 
           alpha = 0.8, 
           size = 0.8, 
           color =  "black") + 
  geom_text(aes(label = n), 
            position = position_dodge(0.9), 
            vjust = -0.5,
            size = 5) +
  labs(title = "Genes of interest interactions distribution",
       x = "Genes of interest",
       y = "count") +
  scale_fill_viridis_d(option = "C") +
  coord_cartesian(clip = "off")
ggsave("plots/006_Interactions_distributions.png", h = 1600, w = 2400, units = "px")

# Plot FC of the targets and regulators of the selected genes
# First we have to take into account that there might be discordant evidences
# for a seconadry gene (being either a regulator or a target) effect on one of
# the primary genes. We choose the correct interaction by selecting the one
# that has most reported Interaction_type. If the score is even we select the 
# "enhanceable" interaction
gint_no_dup <- gint_int %>% 
  group_by(primary, secondary, primary_role, Interaction_type) %>% 
  mutate(n = n()) %>% 
  ungroup(Interaction_type) %>% 
  arrange(desc(n)) %>% 
  slice_head() %>% 
  ungroup() %>%
  select(-n)

# Add FC and remove genes that are not found in the DE analysis
gint_fc <- gint_no_dup %>% 
  inner_join(res_df, by = c("secondary" = "gene")) %>% 
  mutate(fc = ifelse(padj > 0.05, 0, fc)) %>% 
  group_by(primary, secondary) %>% 
  mutate(mean_fc = mean(fc)) %>% 
  ungroup()

for(i in unique(gint_fc$primary)){
  
  gint_fc_reg <- gint_fc %>% 
    filter(primary == i, primary_role == "regulator") %>% 
    mutate(secondary = fct_reorder(secondary, mean_fc))
  
   gint_fc_reg %>% 
     filter(!(mean_fc == 0)) %>% 
     ggplot(aes(secondary, fc, fill = contrast)) +
     geom_bar(stat = "identity", position = "dodge", color = "black", alpha = 0.8) +
     labs(title = paste0(i, " targets log2FC"),
         y = "Target gene",
         x = "log2FoldChange") +
     coord_flip() +
     scale_fill_manual(values = c("lmax_hmyc_vs_neg" = "#153AE0",
                                 "lmax_vs_neg" = "#A2B1F6"),
                      labels = c("LowMaxHighMyc vs Neg", "LowMax vs Neg")) +
     theme_light()
  ggsave(paste0("plots/011_", i, "_targets.png"), h = 2000, w = 2500, units = "px")
  
  # Same plot for the regulators of the genes of interest
  gint_fc_tar <- gint_fc %>% 
    filter(primary == i, primary_role == "target", !(mean_fc == 0)) %>% 
    mutate(secondary = fct_reorder(secondary, mean_fc))
  
  gint_fc_tar %>% 
    filter(!(mean_fc) == 0) %>% 
    ggplot(aes(secondary, fc, fill = contrast)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", alpha = 0.8) +
    labs(title = paste0(i, " regulators log2FC"),
         y = "Regulator gene",
         x = "log2FoldChange") +
    coord_flip() +
    scale_fill_manual(values = c("lmax_hmyc_vs_neg" = "#FB690E",
                                 "lmax_vs_neg" = "#FDC09B"),
                      labels = c("LowMaxHighMyc vs Neg", "LowMax vs Neg")) +
    theme_light()
  ggsave(paste0("plots/011_", i, "_regulators.png"), h = 2000, w = 2500, units = "px")
  
  # Mosaic plots and Fisher exact test
  # The fcs look to be higher when the regulation is positive when myc is 
  # overexpressed and lower when the regulation is negative, compared to Low Max
  reg_mos <- gint_fc_reg %>% 
    filter(!(mean_fc == 0), !(padj > 0.05)) %>% 
    select(primary, secondary, mean_fc, contrast, fc) %>% 
    pivot_wider(names_from = contrast, values_from = fc) %>% 
    filter(!if_any(starts_with("lmax"), is.na)) %>% 
    mutate(obs_regulation = ifelse(mean_fc > 0, "positive", "negative")) %>% 
    mutate(HMyc_greater = lmax_hmyc_vs_neg > lmax_vs_neg) 
  
  if(nrow(reg_mos) > 0){
    
  pos_hmyc_greater <- reg_mos %>% filter(obs_regulation == "positive") %>% pull(HMyc_greater) %>% sum()
  pos_hmyc_smaller <- nrow(filter(reg_mos, obs_regulation == "positive")) - pos_hmyc_greater
  neg_hmyc_greater <- sum(reg_mos[reg_mos$obs_regulation == "negative",]$HMyc_greater == TRUE)
  neg_hmyc_smaller <- nrow(filter(reg_mos, obs_regulation == "negative")) - neg_hmyc_greater
  
  mos <- data.frame(row.names = c("Positive", "Negative"),
                "HMyc_greater" =  c(pos_hmyc_greater, neg_hmyc_greater),
                "HMyc_smaller" = c(pos_hmyc_smaller, neg_hmyc_smaller)) 
  fet <- fisher.test(mos)
  png(paste0("plots/011_", i, "_mosaic_targets.png"), h = 1800, w = 1800, res = 300)
  mosaicplot(mos, 
             main = paste0(i, " targets regulation"),
             color = TRUE,
             cex.axis = 1)
  mtext(side=3, paste0("FET p = ", round(fet$p.value, 3)))
  dev.off()
  }
  
  # Mosaic plots for regulators
  tar_mos <- gint_fc_tar %>% 
    filter(!(mean_fc == 0), !(padj > 0.05)) %>% 
    select(primary, secondary, mean_fc, contrast, fc) %>% 
    pivot_wider(names_from = contrast, values_from = fc) %>% 
    filter(!if_any(starts_with("lmax"), is.na)) %>% 
    mutate(obs_regulation = ifelse(mean_fc > 0, "positive", "negative")) %>% 
    mutate(HMyc_greater = lmax_hmyc_vs_neg > lmax_vs_neg) 
  
  if(nrow(tar_mos) > 0){
    
    pos_hmyc_greater <- tar_mos %>% filter(obs_regulation == "positive") %>% pull(HMyc_greater) %>% sum()
    pos_hmyc_smaller <- nrow(filter(tar_mos, obs_regulation == "positive")) - pos_hmyc_greater
    neg_hmyc_greater <- sum(tar_mos[tar_mos$obs_regulation == "negative",]$HMyc_greater == TRUE)
    neg_hmyc_smaller <- nrow(filter(tar_mos, obs_regulation == "negative")) - neg_hmyc_greater
    
    
    mos <- data.frame(row.names = c("Positive", "Negative"),
                      "HMyc_greater" =  c(pos_hmyc_greater, neg_hmyc_greater),
                      "HMyc_smaller" = c(pos_hmyc_smaller, neg_hmyc_smaller)) 
    fet <- fisher.test(mos)
    png(paste0("plots/011_", i, "_mosaic_regulators.png"), h = 1800, w = 1800, res = 300)
    mosaicplot(mos, 
               main = paste0(i, " regulators expression profiles"),
               color = TRUE,
               cex.axis = 1)
    mtext(side=3, paste0("FET p = ", round(fet$p.value, 3)))
    dev.off()
  }
  
}

  

