library(tidyverse)
library(openxlsx)
library(conflicted)
library(oligo)
conflict_prefer_all("dplyr")

Metadata.GSE2603 <- read.xlsx("TNBC Metadata GSE2603.xlsx")
Metadata.GSE5327 <- read.xlsx("Metadata GSE5327.xlsx")
Metadata.GSE12276 <- read.xlsx("TNBC Metadata GSE12276.xlsx")
Metadata.GSE2034 <- read.xlsx("TNBC Metadata GSE2034.xlsx")


#Metadata combined
Metadata.GSE2603.filt <-  Metadata.GSE2603 %>%
  select(GEO_id, Cohort, LMEvent, LMFS_years)

Metadata.GSE5327.filt <- Metadata.GSE5327 %>%
  select(GEO_id, Cohort, LMEvent, LMFS_years)

# LMS in months
Metadata.GSE2034.filt <- Metadata.GSE2034 %>%
  select(GEO.array, Cohort,  Lung.relapse, MFS) %>%
  rowwise() %>%
  mutate(MFS = MFS/12)
colnames(Metadata.GSE2034.filt) <- colnames(Metadata.GSE2603.filt)

Metadata.GSE12276.filt <- Metadata.GSE12276 %>%
  select(GEO.array, Cohort, Lung.relapse, MFS) %>%
  rowwise() %>%
  mutate(MFS = MFS/12)
colnames(Metadata.GSE12276.filt) <- colnames(Metadata.GSE2603.filt)


Metadata.all <- rbind(Metadata.GSE2603.filt, Metadata.GSE5327.filt, Metadata.GSE2034.filt, Metadata.GSE12276.filt)
table(Metadata.all$LMEvent)


# GSE2603
setwd("GSE2603 Data")
data.raw.2603 <- read.celfiles(list.celfiles(listGzipped = TRUE))
sampleNames(data.raw.2603) <- gsub("\\.CEL\\.gz$", "", sampleNames(data.raw.2603)) 
data.raw.2603 <- data.raw.2603[, sampleNames(data.raw.2603) %in% Metadata.GSE2603$GEO_id]

# GSE5327
setwd("GSE5327 Data")
data.raw.5327 <- read.celfiles(list.celfiles(listGzipped = TRUE))
sampleNames(data.raw.5327) <- gsub("\\.CEL\\.gz$", "", sampleNames(data.raw.5327)) 

# GSE 2034
setwd("GSE2034 Data")
data.raw.2034 <- read.celfiles(list.celfiles(listGzipped = TRUE))
sampleNames(data.raw.2034) <- gsub("\\.CEL\\.gz$", "", sampleNames(data.raw.2034)) 
data.raw.2034 <- data.raw.2034[, sampleNames(data.raw.2034) %in% Metadata.GSE2034$GEO.array]

# GSE 12276
setwd("GSE12276 Data")
data.raw.12276 <- read.celfiles(list.celfiles(listGzipped = TRUE))
sampleNames(data.raw.12276) <- gsub("\\.CEL\\.gz$", "", sampleNames(data.raw.12276))
data.raw.12276 <- data.raw.12276[, sampleNames(data.raw.12276) %in% Metadata.GSE12276$GEO.array]



setwd("C:/Users/idisb/Desktop/Andrés/02 - Projects/05 - DNAm organ specific/MethMice LMFS")

# Data normalization
data.norm.2603 <- rma(data.raw.2603, background = T, normalize = T)
data.RNA.2603 <- exprs(data.norm.2603)
data.RNA.2603 <- as.data.frame(data.RNA.2603)
boxplot(data.raw.2603, target = "core", 
        main = "Boxplot of log2-intensitites for the raw data")
boxplot(data.norm.2603, target = "core", 
        main = "Boxplot of log2-intensitites for the normalized data", col = "gray")



data.norm.5327 <- rma(data.raw.5327, background = T, normalize = T)
data.RNA.5327 <- exprs(data.norm.5327)
data.RNA.5327 <- as.data.frame(data.RNA.5327)
boxplot(data.raw.5327, target = "core", 
        main = "Boxplot of log2-intensitites for the raw data")
boxplot(data.norm.5327, target = "core", 
        main = "Boxplot of log2-intensitites for the normalized data", col = "gray")



data.norm.2034 <- rma(data.raw.2034, background = T, normalize = T)
data.RNA.2034 <- exprs(data.norm.2034)
data.RNA.2034 <- as.data.frame(data.RNA.2034)
boxplot(data.raw.2034, target = "core", 
        main = "Boxplot of log2-intensitites for the raw data")
boxplot(data.norm.2034, target = "core", 
        main = "Boxplot of log2-intensitites for the normalized data", col = "gray")


data.norm.12276 <- rma(data.raw.12276, background = T, normalize = T)
data.RNA.12276 <- exprs(data.norm.12276)
data.RNA.12276 <- as.data.frame(data.RNA.12276)
boxplot(data.raw.12276, target = "core", 
        main = "Boxplot of log2-intensitites for the raw data")
boxplot(data.norm.12276, target = "core", 
        main = "Boxplot of log2-intensitites for the normalized data", col = "gray")



library(Biobase)
library(sva)

eset1 <- data.norm.2603
eset2 <- data.norm.5327
eset3 <- data.norm.2034
eset4 <- data.norm.12276


common_features <- Reduce(
  intersect,
  list(
    featureNames(eset1),
    featureNames(eset2),
    featureNames(eset3),
    featureNames(eset4)))

#
eset1 <- eset1[common_features, ]
eset2 <- eset2[common_features, ]
eset3 <- eset3[common_features, ]
eset4 <- eset4[common_features, ]


exprs_combined <- cbind(
  exprs(eset1),
  exprs(eset2),
  exprs(eset3),
  exprs(eset4))

batch <- c(
  rep(1, ncol(eset1)),
  rep(2, ncol(eset2)),
  rep(3, ncol(eset3)),
  rep(4, ncol(eset4)))

# ComBat
mod <- model.matrix(
  ~ LMEvent,
  data = Metadata.all)

exprs_corrected <- ComBat(dat=exprs_combined, batch=batch, mod = mod)

head(exprs_corrected[1:5, 1:5])


#BiocManager::install("hgu133a.db")
library(hgu133a.db)
library(biomaRt)

RNA.anot <- hgu133aALIAS2PROBE %>%
  as.data.frame() %>%
  inner_join(hgu133aCHR %>%
               as.data.frame(), by = "probe_id") %>%
  mutate(chromosome = paste0("chr", chromosome))

Genes.interest <- read.xlsx("Final Intersect MM, AURORA DNAm, RNAseq.xlsx")


# Differential expression between LM2 and WT
library(limma)
library(ggpubr)
data.all <- exprs_corrected[,Metadata.all$GEO_id]

genes <- unique(KM.data$probe_id)
gene_symbol <- unique(KM.data$alias_symbol)

for (i in 1:length(genes)) {
  gene <- genes[i]
  
  KM.genes <- KM.data %>%
    drop_na(Exprs) %>%
    select(GEO_id, LMFS_years, LMEvent, Exprs, probe_id, alias_symbol, Cohort) %>%
    filter(probe_id == gene) %>%
    filter(!is.na(LMFS_years)) %>%
    filter(!is.na(LMEvent)) %>%
    mutate(LMFS_years = as.numeric(LMFS_years)) %>%
    mutate(LMFS_years = LMFS_years * 12)
  
  
  gene_symbol <- unique(KM.genes$alias_symbol)

  KM.genes <- KM.genes %>%
    select(-probe_id, - alias_symbol)
  colnames(KM.genes) <- c("Sample ID", "Survival time", "Survival event", gene_symbol, "Cohort")
  
file_name <- paste0("KM data expanded/KM_", gene_symbol,"-", gene, ".csv")
  
write_csv(KM.genes, file = file_name)
}

# FORESTPLOT de los KM
library(grid)
library(forestploter)

Gene <- c("AK1", "AUTS2", "CAPS", "CD5", "CHI3L2", "CYB5R3",
          "EYA1", "FKBP4", "HOXA5", "MEIS2", "MICAL3", "PRKCZ", "PRPH", 
          "SLC2A5", "TCIRG1", "TPI1", "ZBTB17")

HR <- c(
  2.01, 0.52, 0.70, 0.52, 1.25, 1.41, 0.59,
  0.49, 1.33, 0.72, 2.12, 0.57, 1.91,
  2.37, 1.62, 2.92, 1.73)

lower_CI <- c(
  1.04, 0.25, 0.36, 0.22, 0.63, 0.71, 0.31,
  0.26, 0.65, 0.37, 1.09, 0.28, 1.00,
  1.22, 0.85, 1.51, 0.87
)

upper_CI <- c(
  3.88, 1.08, 1.37, 1.26, 2.48, 2.82, 1.13,
  0.94, 2.76, 1.38, 4.12, 1.19, 3.67,
  4.61, 3.11, 5.64, 3.45
)

logrank_P <- c(
  0.033, 0.076, 0.30, 0.14, 0.53, 0.32, 0.11,
  0.028, 0.44, 0.32, 0.024, 0.13, 0.047,
  0.0087, 0.14, 0.00082, 0.11
)


low <- c(
  95, 101, 94, 130, 63, 129, 68,
  48, 129, 53, 130, 104, 122,
  100, 114, 106, 130
)

high <- c(
  78, 72, 79, 43, 110, 44, 105,
  125, 44, 120, 43, 69, 51,
  73, 59, 67, 43
)

test.forest <- data.frame(
  Genes = Gene,
  low_n = low,
  high_n = high,
  logrank_P = logrank_P)

test.forest$` ` <- paste(rep(" ", 30), collapse = " ") 


test.forest <- test.forest %>%
  select(Genes, low_n, high_n, logrank_P, ` `) %>%
  mutate("   HR (95% CI)" = sprintf("   %.2f (%.2f to %.2f)",
                                    HR, lower_CI, upper_CI))

test.forest$logrank_P <- round(test.forest$logrank_P, 3)
colnames(test.forest) <- c(
  "Gene",
  "Low (n)",
  "High (n)",
  "Logrank P",
  " ",
  "HR (95% CI)"
)

test.forest$`Logrank P` <- signif(test.forest$`Logrank P`, 3)

tm <- forest_theme(base_size = 10,
                   ci_pch = 15,
                   ci_col = "#31A331",
                   ci_fill = "black",
                   ci_alpha = 0.8,
                   ci_lty = 1,
                   ci_lwd = 2,
                   ci_Theight = 0.3, # Set a T end at the end of CI 
                   refline_gp  = gpar(lwd = 1, lty = "dashed", col = "grey20"),
                   vertline_lwd = 1,
                   vertline_lty = "dashed",
                   vertline_col = "grey20",
                   summary_fill = "#4575b4",
                   summary_col = "#4575b4",
                   footnote_gp = gpar(cex = 1, fontface = "italic", col = "blue"),
                   core = list(padding = unit(c(10, 6), "mm")))


p <- forest(test.forest,
            est = HR,
            lower = lower_CI, 
            upper = upper_CI,
            sizes = 0.75,
            ci_column = 5,
            ref_line = 1,
            x_trans = "log2",
            theme = tm, xlim = c(0 , 8),
            ticks_at = c(0.125, 0.25, 0.5, 1, 2, 4, 8))

g <- edit_plot(p, row = c(1,8, 11, 13, 14, 16), which = "background",
               gp = gpar(fill = "#E9FFE9"))
f <- edit_plot(g, row = c(2,9), which = "background",
               gp = gpar(fill = "tomato"))
plot(g) 





# GENES INTEREST ENHANCERS
setwd("C:/Users/idisb/Desktop/Andrés/02 - Projects/05 - DNAm organ specific")

Genes.interest.enh <- read.xlsx("Genes dysregulated by DNAm changes.xlsx")
Genes.interest

Genes.interest.both <- c(Genes.interest.enh$symbol, Genes.interest$geneNames)


library(limma)
library(ggpubr)
data.all <- exprs_corrected[,Metadata.all$GEO_id]



# Boxplot probes/genes of interest
KM.data.enh <- as.data.frame(data.all)  %>%
  rownames_to_column(var = "probe_id") %>%
  pivot_longer(-probe_id, names_to = "GEO_id", values_to = "Exprs") %>%
  inner_join(Metadata.all, by = "GEO_id") %>%
  inner_join(RNA.anot, by = c("probe_id")) %>%
  filter(alias_symbol %in% Genes.interest.both) %>%
  mutate(Group = case_when(LMEvent == 0 ~ "No",
                           LMEvent == 1 ~ "Yes"))

KM.gene <- KM.data.enh %>%
  group_by(GEO_id, alias_symbol) %>%
  summarise(Exprs = mean(Exprs, na.rm = TRUE), .groups = "drop") %>%
  left_join(Metadata.all, by = "GEO_id") %>%
  mutate(Group = case_when(LMEvent == 0 ~ "No",
                           LMEvent == 1 ~ "Yes"))


library(dplyr)
library(survival)
library(survminer)
library(purrr)


get_best_cutoff <- function(df) {
  
  q1 <- quantile(df$Exprs, 0.25, na.rm = TRUE)
  q3 <- quantile(df$Exprs, 0.75, na.rm = TRUE)
  
  cutoffs <- unique(df$Exprs[df$Exprs >= q1 & df$Exprs <= q3])
  
  if (length(cutoffs) < 5) return(NULL)
  
  res <- lapply(cutoffs, function(cut) {
    
    df$group <- ifelse(df$Exprs > cut, "high", "low")
    df$group <- factor(df$group, levels = c("low", "high"))
    
    if (min(table(df$group)) < 5) return(NULL)
    
    # Log-rank
    sd <- survdiff(Surv(LMFS_years, LMEvent) ~ group, data = df)
    pval <- 1 - pchisq(sd$chisq, df = 1)
    
    # Cox
    cox <- coxph(Surv(LMFS_years, LMEvent) ~ group, data = df)
    HR <- summary(cox)$coefficients[1, "exp(coef)"]
    
    data.frame(
      cutoff = cut,
      pval = pval,
      HR = HR
    )
  })
  
  res <- bind_rows(res)
  if (nrow(res) == 0) return(NULL)
  
  # FDR
  res$FDR <- p.adjust(res$pval, method = "BH")
  
  best <- res %>%
    filter(FDR == min(FDR, na.rm = TRUE)) %>%
    arrange(desc(HR)) %>%
    slice(1)
  
  return(best)
}


#
results <- map_dfr(unique(KM.gene$alias_symbol), function(gene) {
  
  df <- KM.gene %>%
    filter(alias_symbol == gene) %>%
    mutate(LMFS_years = as.numeric(LMFS_years)) %>%
    filter(!is.na(LMFS_years), !is.na(LMEvent), !is.na(Exprs))
  
  if (nrow(df) < 20) return(NULL)
  
  best <- get_best_cutoff(df)
  if (is.null(best)) return(NULL)
  
  # Aplicar cutoff final
  df$group <- ifelse(df$Exprs > best$cutoff, "high", "low")
  df$group <- factor(df$group, levels = c("low", "high"))
  
  n_high <- sum(df$group == "high")
  n_low  <- sum(df$group == "low")
  
  # Cox final
  cox <- coxph(Surv(LMFS_years, LMEvent) ~ group, data = df)
  cox_sum <- summary(cox)
  
  data.frame(
    gene = gene,
    cutoff = best$cutoff,
    n_high = n_high,
    n_low = n_low,
    HR = cox_sum$coefficients[1, "exp(coef)"],
    CI_low = cox_sum$conf.int[1, "lower .95"],
    CI_high = cox_sum$conf.int[1, "upper .95"],
    p_logrank = best$pval,
    FDR = best$FDR
  )
})



analyze_probe <- function(df) {
  q1 <- quantile(df$Exprs, 0.25, na.rm = TRUE)
  q3 <- quantile(df$Exprs, 0.75, na.rm = TRUE)
  
  cutoffs <- unique(df$Exprs[df$Exprs > q1 & df$Exprs < q3])
  
  if (length(cutoffs) < 5) return(NULL)
  
  res <- lapply(cutoffs, function(cut) {
    
    df$group <- ifelse(df$Exprs > cut, "high", "low")
    df$group <- factor(df$group, levels = c("low", "high"))
    
    if (min(table(df$group)) < 5) return(NULL)
    
    # Log-rank
    sd <- survdiff(Surv(LMFS_months, LMEvent) ~ group, data = df)
    pval <- 1 - pchisq(sd$chisq, df = 1)
    
    # Cox
    cox <- coxph(Surv(LMFS_months, LMEvent) ~ group, data = df)
    cox_sum <- summary(cox)
    
    HR <- cox_sum$coefficients[1, "exp(coef)"]
    CI_low  <- cox_sum$conf.int[1, "lower .95"]
    CI_high <- cox_sum$conf.int[1, "upper .95"]
    
    data.frame(
      cutoff = cut,
      pval = pval,
      HR = HR,
      CI_low = CI_low,
      CI_high = CI_high,
      n_high = sum(df$group == "high"),
      n_low  = sum(df$group == "low")
    )
  })
  
  res <- bind_rows(res)
  if (nrow(res) == 0) return(NULL)
  
  res$FDR <- p.adjust(res$pval, method = "BH")
  
  best <- res %>%
    filter(FDR == min(FDR, na.rm = TRUE)) %>%
    arrange(desc(abs(log(HR)))) %>%
    slice(1)
  
  return(best)
}



library(purrr)

results_probe <- map_dfr(unique(KM.data.enh$probe_id), function(probe) {
  
  df <- KM.data.enh %>%
    filter(probe_id == probe) %>%
    mutate(LMFS_years = as.numeric(LMFS_years)) %>%
    filter(!is.na(LMFS_years), !is.na(LMEvent), !is.na(Exprs))
  
  if (nrow(df) < 20) return(NULL)
  
  best <- analyze_probe(df)
  if (is.null(best)) return(NULL)
  
  data.frame(
    probe_id = probe,
    gene = unique(df$alias_symbol),
    cutoff = best$cutoff,
    n_high = best$n_high,
    n_low = best$n_low,
    HR = best$HR,
    CI_low = best$CI_low,
    CI_high = best$CI_high,
    p_logrank = best$pval,
    FDR_probe = best$FDR)})


best_probes <- results_probe %>%
  group_by(gene) %>%
  mutate(FDR_within_gene = p.adjust(p_logrank, method = "BH")) %>%
  filter(FDR_within_gene == min(FDR_within_gene, na.rm = TRUE)) %>%
  arrange(desc(abs(log(HR)))) %>%
  slice(1) %>%
  ungroup()


best_probes.final <- best_probes %>%
  mutate(group = case_when(gene %in% Genes.interest$geneNames & gene %in% Genes.interest.enh$symbol ~ "Promoter and enhancer",
                           gene %in% Genes.interest$geneNames ~ "Promoter",
                           gene %in% Genes.interest.enh$symbol ~ "Enhancer"
                           ))

table(best_probes.final$group)


write.xlsx(best_probes.final, file = "LMFS for both enhancer and promoter DNAm regulated.xlsx")


best_probes.final.sig <- best_probes.final %>%
  filter(FDR_within_gene < 0.05) %>%
  arrange(FDR_within_gene)


# AURORA 
AURORA.res.naive <- read_tsv("C:/Users/idisb/Desktop/Andrés/02 - Projects/05 - DNAm organ specific/AURORA RNA-seq LU vs PT.txt")

best_probes.final.sig <- best_probes.final.sig %>%
  inner_join(AURORA.res.naive, by = c("gene" = "gene_name")) 

best_probes.final.sig.filt <- best_probes.final.sig %>%
  mutate(HR.direction = case_when(HR > 1 ~ "Up", T ~ "Down")) %>%
  filter(Expression %in% "Upregulated" & HR.direction %in% "Up"|
           Expression %in% "Downregulated" & HR.direction %in% "Down")

best_probes.final.sig.filt.no <- best_probes.final.sig %>%
  mutate(HR.direction = case_when(HR > 1 ~ "Up", T ~ "Down")) %>%
  filter(!Expression %in% "Upregulated" & HR.direction %in% "Up"|
           !Expression %in% "Downregulated" & HR.direction %in% "Down")




##

library(dplyr)
library(survival)
library(survminer)

plot_KM_probe <- function(probe, data, show_CI = TRUE) {
  
  df <- data %>%
    filter(probe_id %in% probe) %>%
    mutate(
      LMFS_years = as.numeric(LMFS_years),
      LMFS_months = LMFS_years * 12) %>%
    filter(!is.na(LMFS_months), !is.na(LMEvent), !is.na(Exprs)) %>%
    mutate(LMFS_months = case_when(LMFS_months > 60 ~ 61, T ~ LMFS_months))
  
  gene_name <- unique(df$alias_symbol)
  
  cutoff_value <- best_probes %>%
    filter(probe_id == probe) %>%
    pull(cutoff)
  
  df$Exprs_group <- ifelse(df$Exprs > cutoff_value, "high", "low")
  df$Exprs_group <- factor(df$Exprs_group, levels = c("low", "high"))
  
  fit <- survfit(Surv(LMFS_months, LMEvent) ~ Exprs_group, data = df)
  
  survdiff_res <- survdiff(Surv(LMFS_months, LMEvent) ~ Exprs_group, data = df)
  pval <- 1 - pchisq(survdiff_res$chisq, df = 1)
  
  cox <- coxph(Surv(LMFS_months, LMEvent) ~ Exprs_group, data = df)
  cox_sum <- summary(cox)
  
  HR  <- round(cox_sum$coefficients[1, "exp(coef)"], 2)
  CI_low  <- round(cox_sum$conf.int[1, "lower .95"], 2)
  CI_high <- round(cox_sum$conf.int[1, "upper .95"], 2)
  
  p <- ggsurvplot(
    fit,
    data = df,
    xlim = c(0, 60),
    break.time.by = 12, 
    risk.table = TRUE,
    conf.int = show_CI, 
    
    legend.labs = c("Low", "High"),
    palette = c("black", "red"),
    
    title = gene_name,
    xlab = "Time (months)",
    ylab = "LMFS probability"
  )
  
  annotation_text <- paste0(
    "HR = ", HR, " (", CI_low, " - ", CI_high, ")\n",
    "Log-rank p = ", signif(pval, 3)
  )
  
  p$plot <- p$plot +
    annotate(
      "text",
      x = 60, y = 0.1,     
      label = annotation_text,
      hjust = 1,
      size = 4.5, color = "black"
    ) + theme_bw() + 
    theme(
      text = element_text(size = 12, color = "black"),   
      axis.text = element_text(size = 11, color = "black"),
      axis.title = element_text(size = 11,color = "black" ),
      legend.position = "bottom", 
    )
  
  p$table <- p$table +
    theme(
      text = element_text(size = 10, color = "black"),  
      axis.text = element_text(size = 10, color = "black"),
      axis.title = element_text(size = 10,color = "black" ))
  
  return(p)}

plot_KM_probe("203574_at", KM.data.enh, show_CI = F) 


## PROMOTERS
AK1 <- plot_KM_probe(best_probes.final.sig %>% filter(gene %in% "AK1") %>% pull(probe_id), 
                     KM.data.enh, show_CI = F) 

TPI1 <- plot_KM_probe(best_probes.final.sig %>% filter(gene %in% "TPI1") %>% pull(probe_id), 
                     KM.data.enh, show_CI = F) 

SLC2A5 <- plot_KM_probe(best_probes.final.sig %>% filter(gene %in% "SLC2A5") %>% pull(probe_id), 
                     KM.data.enh, show_CI = F) 
MICAL3 <- plot_KM_probe(best_probes.final.sig %>% filter(gene %in% "MICAL3") %>% pull(probe_id), 
                        KM.data.enh, show_CI = F) 

library(patchwork)
combined <- arrange_ggsurvplots(
  list(AK1, TPI1, SLC2A5, MICAL3),
  ncol = 2,
  nrow = 2)

combined


## ENHANCERS
ENO1 <- plot_KM_probe(best_probes.final.sig %>% filter(gene %in% "ENO1") %>% pull(probe_id), 
                     KM.data.enh, show_CI = F) 

SLC16A3 <- plot_KM_probe(best_probes.final.sig %>% filter(gene %in% "SLC16A3") %>% pull(probe_id), 
                      KM.data.enh, show_CI = F) 

CDK5RAP3 <- plot_KM_probe(best_probes.final.sig %>% filter(gene %in% "CDK5RAP3") %>% pull(probe_id), 
                        KM.data.enh, show_CI = F) 

PAX5 <- plot_KM_probe(best_probes.final.sig %>% filter(gene %in% "PAX5") %>% pull(probe_id), 
                          KM.data.enh, show_CI = F) 

combined2 <- arrange_ggsurvplots(
  list(ENO1, SLC16A3, CDK5RAP3, PAX5),
  ncol = 2,
  nrow = 2)

combined2


## FORESTPLOT 
library(grid)
library(forestploter)

Figure6 <- best_probes.final.sig.filt %>%
  filter(gene %in% c("AK1", "MICAL3", "SLC2A5", "TPI1",  "ENO1", "NFIL3", "PAX5", "SLC16A3"))
Figure6$gene <- factor(Figure6$gene, levels = c("AK1", "MICAL3", "SLC2A5", "TPI1",  "ENO1", "NFIL3", "PAX5", "SLC16A3"))

Figure6 <- Figure6 %>%
  select(1:11) %>% select(-p_logrank, -FDR_probe, -cutoff) %>%
  arrange(gene)

HR. <- Figure6$HR
CI.low <- Figure6$CI_low
CI.high <- Figure6$CI_high

Figure6$` ` <- paste(rep(" ", 30), collapse = " ") 

Figure6 <- Figure6 %>%
  select(probe_id, gene, n_low, n_high, FDR_within_gene, ` `, everything()) %>%
  mutate("   HR (95% CI)" = sprintf("   %.2f (%.2f to %.2f)",
                                    HR, CI_low, CI_high)) %>%
  select(-HR, -CI_low, -CI_high)

Figure6$FDR_within_gene <- formatC(Figure6$FDR_within_gene, format = "e", digits = 2)

tm <- forest_theme(base_size = 10,
                   ci_pch = 15,
                   ci_col = "#31A331",
                   ci_fill = "black",
                   ci_alpha = 0.8,
                   ci_lty = 1,
                   ci_lwd = 2,
                   ci_Theight = 0.3, # Set a T end at the end of CI 
                   refline_gp  = gpar(lwd = 1, lty = "dashed", col = "grey20"),
                   vertline_lwd = 1,
                   vertline_lty = "dashed",
                   vertline_col = "grey20",
                   summary_fill = "#4575b4",
                   summary_col = "#4575b4",
                   footnote_gp = gpar(cex = 1, fontface = "italic", col = "blue"),
                   core = list(padding = unit(c(12, 12), "mm")))


p <- forest(Figure6,
            est = HR.,
            lower = CI.low, 
            upper = CI.high,
            sizes = 0.75,
            ci_column = 6,
            ref_line = 1,
            x_trans = "log2",
            theme = tm, xlim = c(0 , 16),
            ticks_at = c(0.125, 0.25, 0.5, 1, 2, 4, 8))


plot(p) 


plot_KM_probe <- function(probe, data, show_CI = TRUE) {
  
  df <- data %>%
    filter(probe_id %in% probe) %>%
    mutate(
      LMFS_years = as.numeric(LMFS_years),
      LMFS_months = LMFS_years * 12) %>%
    filter(!is.na(LMFS_months), !is.na(LMEvent), !is.na(Exprs)) %>%
    mutate(LMFS_months = case_when(LMFS_months > 60 ~ 61, T ~ LMFS_months))
  
  gene_name <- unique(df$alias_symbol)
  
  cutoff_value <- best_probes %>%
    filter(probe_id == probe) %>%
    pull(cutoff)
  
  df$Exprs_group <- ifelse(df$Exprs > cutoff_value, "high", "low")
  df$Exprs_group <- factor(df$Exprs_group, levels = c("low", "high"))
  
  # KM
  fit <- survfit(Surv(LMFS_months, LMEvent) ~ Exprs_group, data = df)
  
  # Log-rank
  survdiff_res <- survdiff(Surv(LMFS_months, LMEvent) ~ Exprs_group, data = df)
  pval <- 1 - pchisq(survdiff_res$chisq, df = 1)
  
  # Cox
  cox <- coxph(Surv(LMFS_months, LMEvent) ~ Exprs_group, data = df)
  cox_sum <- summary(cox)
  
  HR  <- round(cox_sum$coefficients[1, "exp(coef)"], 2)
  CI_low  <- round(cox_sum$conf.int[1, "lower .95"], 2)
  CI_high <- round(cox_sum$conf.int[1, "upper .95"], 2)
  
  # Base plot
  p <- ggsurvplot(
    fit,
    data = df,
    xlim = c(0, 60),
    break.time.by = 12,
    risk.table = F,
    conf.int = show_CI, 
    
    legend.labs = c("Low", "High"),
    palette = c("black", "red"),
    
    title = gene_name,
    xlab = "Time (months)",
    ylab = "LMFS probability"
  )
    annotation_text <- paste0(
    "HR = ", HR, " (", CI_low, " - ", CI_high, ")\n",
    "Log-rank p = ", signif(pval, 3)
  )
  
  p$plot <- p$plot +
    annotate(
      "text",
      x = 60, y = 0.1,
      label = annotation_text,
      hjust = 1,
      size = 4.5, color = "black"
    ) + theme_bw() + 
    theme(
      text = element_text(size = 12, color = "black"),   
      axis.text = element_text(size = 11, color = "black"),
      axis.title = element_text(size = 11,color = "black" ),
      legend.position = "bottom")
  
  return(p)}


FKBP4 <- plot_KM_probe(best_probes.final.sig %>% filter(gene %in% "FKBP4") %>% pull(probe_id), 
                     KM.data.enh, show_CI = F)

PPFIBP2 <- plot_KM_probe(best_probes.final.sig %>% filter(gene %in% "PPFIBP2") %>% pull(probe_id), 
                      KM.data.enh, show_CI = F)

IRX5 <- plot_KM_probe(best_probes.final.sig %>% filter(gene %in% "IRX5") %>% pull(probe_id), 
                        KM.data.enh, show_CI = F)

OBSL1 <- plot_KM_probe(best_probes.final.sig %>% filter(gene %in% "OBSL1") %>% pull(probe_id), 
                        KM.data.enh, show_CI = F)


combined.x <- arrange_ggsurvplots(
  list(FKBP4, PPFIBP2, IRX5, OBSL1),
  ncol = 4,
  nrow = 1)

combined.x
