library(tidyverse)

# TARGET survival ----
expmat_orig <- readRDS("data/target_NBL-expmat.rds")
survival <- readRDS("data/target_NBL-survival.rds")

# Make sure that the patients name match between expmat and clinical
colnames(expmat_orig) <- gsub(pattern = ".{8}$", replacement =  "", colnames(expmat_orig))

## Remove genes with low expression
expvar <- apply(expmat_orig, 1, var)
expvar <- names(expvar[expvar > 0.1])
expmat <- expmat_orig[expvar,]

## Intersect: analysis must be done for the same patients (they must have both survival info and the biomarkers)
common <- intersect(colnames(expmat_orig), names(survival))
expmat <- expmat[,common]
survival <- survival[common]

surv_target <- data.frame(Symbol = character(),
                          Cox_coef = numeric(),
                          Cox_p = numeric(),
                          Survival_p = numeric(),
                          Pred_sign = character(),
                          stringsAsFactors = F)

pb <- txtProgressBar(min = 0, max = nrow(expmat), initial = 0, style = 3)

if(!file.exists("data/Survival_analysis_target.rds")){
  
  for (i in 1:nrow(expmat)) {
    
    library(survival)
    
    mygene <- rownames(expmat)[i]
    
    #### Fit Cox model
    res.cox <- summary(coxph(survival ~ expmat[mygene,]))
    
    #### Survfit
    oritrack <- expmat[as.character(mygene),]
    
    #### Define types for survfit
    track <- oritrack
    track[] <- "Other"
    track[oritrack>median(oritrack)] <- paste0(mygene,"up")
    track[oritrack<=median(oritrack)] <- paste0(mygene,"dn")
    track <- as.factor(track)
    track <- relevel(track,ref=paste0(mygene,"dn"))
    
    #### Survdiff for all vs. all p-value
    sdiff <- survdiff(survival~track)
    p <- 1-pchisq(sdiff$chisq, df=length(sdiff$n) - 1)
    
    
    #### Function to determine predictor sign
    events_difference <- sdiff$obs-sdiff$exp
    if(events_difference[1] > events_difference[length(events_difference)]){
      sign<-"neg"
    } else {
      sign<-"pos"
    }
    
    #### Group results
    surv_target[i,"Symbol"] <- mygene
    surv_target[i,"Cox_coef"] <- res.cox$coefficients[1,1]
    surv_target[i,"Cox_p"] <- res.cox$coefficients[1,5] 
    surv_target[i,"Survival_p"] <- p
    surv_target[i,"Pred_sign"] <- sign
    
    setTxtProgressBar(pb, value = i)
  }
  
  surv_target <- surv_target %>% 
    mutate(Survival_padj = p.adjust(Survival_p, method = "BH"),
           Cox_padj = p.adjust(Cox_p, method = "BH"))
  saveRDS(surv_target, file = "data/Survival_analysis_target.rds")
  
}else{surv_target <- readRDS("data/Survival_analysis_target.rds")}


# Kocak survival ----
# load('data/kocak_NBL-survival.rda')
# load('data/kocak_NBL-expmat.rda')
# 
# common <- intersect(colnames(expmat),names(survival))
# expmat <- expmat[,common]
# survival <- survival[common]
# 
# surv_kocak <- data.frame(Symbol = character(),
#                          Cox_coef = numeric(),
#                          Cox_p = numeric(),
#                          Survival_p = numeric(),
#                          Pred_sign = character(),
#                          stringsAsFactors = F)
# 
# pb <- txtProgressBar(min = 0, max = nrow(expmat), initial = 0, style = 3)
# 
# if(!file.exists("data/Survival_analysis_kocak.rda")){
#   for (i in 1:nrow(expmat)) {
#     
#     library(survival)
#     
#     mygene <- rownames(expmat)[i]
#     
#     #### Fit Cox model
#     res.cox <- summary(coxph(survival ~ expmat[mygene,]))
#     
#     #### Survfit
#     oritrack <- expmat[as.character(mygene),]
#     
#     #### Define types for survfit
#     track <- oritrack
#     track[] <- "Other"
#     track[oritrack>median(oritrack)] <- paste0(mygene,"up")
#     track[oritrack<=median(oritrack)] <- paste0(mygene,"dn")
#     track <- as.factor(track)
#     track <- relevel(track,ref=paste0(mygene,"dn"))
#     
#     #### Survdiff for all vs. all p-value
#     sdiff <- survdiff(survival~track)
#     p <- 1-pchisq(sdiff$chisq, df=length(sdiff$n) - 1)
#     
#     
#     #### Function to determine predictor sign
#     events_difference < -sdiff$obs-sdiff$exp
#     if(events_difference[1]>events_difference[length(events_difference)]){
#       sign<-"neg"
#     } else {
#       sign<-"pos"
#     }
#     
#     #### Group results
#     surv_kocak[i,"Symbol"] <- mygene
#     surv_kocak[i,"Cox_coef"] <- res.cox$coefficients[1,1]
#     surv_kocak[i,"Cox_p"] <- res.cox$coefficients[1,5] 
#     surv_kocak[i,"Survival_p"] <- p
#     surv_kocak[i,"Pred_sign"] <- sign
#     
#     setTxtProgressBar(pb, value = i)
#     
#   }
#   save(surv_kocak, file = "data/Survival_analysis_kocak.rda")
#   
# }else{load("data/Survival_analysis_kocak.rda")}


# Tried the code below, but there too few genes in common (~400 after filtering for orthologs)
# surv <- surv_kocak %>% 
#   left_join(surv_target, by = "Symbol", suffix = c("_kocak", "_target")) %>% 
#   filter(Survival_p_target < 0.05 & Survival_p_kocak < 0.05) %>% 
#   filter(!is.na(Survival_p_target))


# Create conversion map based on TARGET expression ----
ids <- read.delim("data/diopt_drosophila_human_03feb2021.tsv", as.is = TRUE) # Downloaded from https://www.flyrnai.org/diopt

orth_orig <- ids %>% 
  distinct(Fly.Gene.ID, .keep_all = T)

tarmeans <- data.frame(nbl_symbol = rownames(expmat),
                       expression = rowMeans(expmat))

if(!file.exists("data/ort_flyrnai_TARGET.rds")){
  ort <- orth_orig %>% 
    left_join(tarmeans, by = c("Human.Symbol" = "nbl_symbol")) %>% 
    group_by(FlyBaseID) %>% 
    slice_max(Weighted.Score) %>% 
    filter(!is.na(expression)) %>%
    slice_max(expression) %>% 
    
    # If the orthologs have the same expression, select only the first
    filter(row_number(FlyBaseID) == 1) %>% 
    select(FlyBaseID, Fly.Symbol, Human.Symbol) %>% 
    ungroup()
  
  saveRDS(ort, "data/ort_flyrnai_TARGET.rds")
}else{ort <- readRDS("ort_flyrnai_TARGET.rds")}


