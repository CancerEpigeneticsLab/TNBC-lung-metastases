
library(tidyverse)
library(openxlsx)
library(DESeq2)
library(conflicted)
conflict_prefer_all("dplyr")

metadata.AURORA <- read.xlsx("Metadata AURORA USA cohort.xlsx") %>%
  filter(Sample.Type %in% c("Primary", "Metastasis")) %>%
  filter(RNA.Seq.FreezeSet.123 %in% "Sequenced") %>%
  filter((Inferred.ER.from.RNAseq %in% "Negative" & 
       Inferred.PR.from.RNAseq %in% "Negative" & 
       Inferred.HER2.from.RNAseq %in% "Negative") |
      (Profiled.ER %in% "Negative" & 
         Profiled.PR %in% "Negative" & 
         Profiled.HER2 %in% "Negative")) %>%
  filter(Tissue.primary.type.treatment %in% c("Pre-treatment", "Metastasis")) %>%
  filter(DNA.Methylation.FreezeSet.131 %in% "Sequenced")

write.xlsx(metadata.AURORA, file = "Supplementary Table 1. AURORA US cohort.xlsx")

table(metadata.AURORA$Sample.Type)


table <- table(metadata.AURORA$Anatomic.Site.Simplified)[table(metadata.AURORA$Anatomic.Site.Simplified) > 2]
names(table)


metadata.AURORA <- metadata.AURORA %>%
  filter(Anatomic.Site.Simplified %in% names(table))

table(metadata.AURORA$Anatomic.Site.Simplified)


BRCA.raw <- read_tsv("AURORA GSE209998_AUR_129_raw_counts.txt")
colnames(BRCA.raw)[1] <- "gene_name"

BRCA.raw <- BRCA.raw %>%
  filter(!duplicated(gene_name)) %>%
  column_to_rownames(var = "gene_name")


# Ajustamos los barcodes
colnames(BRCA.raw) <- case_when(
  grepl("TTM10", colnames(BRCA.raw)) ~ substr(colnames(BRCA.raw), 1, 16),
  grepl("NT", colnames(BRCA.raw)) ~ substr(colnames(BRCA.raw), 1, 14),
  TRUE ~ substr(colnames(BRCA.raw), 1, 15))

setdiff(colnames(BRCA.raw), metadata.AURORA$BCR.Sample.barcode)

BRCA.raw <- round(BRCA.raw)
BRCA.raw <- as.matrix(BRCA.raw)

colData <- metadata.AURORA %>%
  column_to_rownames(var = "BCR.Sample.barcode")

BRCA.raw <- BRCA.raw[, rownames(colData)]


##
dds <- DESeqDataSetFromMatrix(countData = BRCA.raw,
                              colData = colData,
                              design = ~ Anatomic.Site.Simplified)

dds <- estimateSizeFactors(dds)

AURORA.Matrix <- counts(dds, normalized = TRUE) 

AURORA.Matrix <- AURORA.Matrix %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene_name")



### Load genes of interest
Genes.interest <- read.xlsx("Final Intersect MM, AURORA DNAm, RNAseq.xlsx")

AURORA.Matrix <- AURORA.Matrix %>%
  filter(gene_name %in% c(Genes.interest$geneNames)) %>%
  column_to_rownames(var = "gene_name")


AURORA.Matrix.log <- log2(AURORA.Matrix + 1)

Aurora.long <- AURORA.Matrix.log %>%
  rownames_to_column(var = "gene_name") %>%
  pivot_longer(-gene_name, values_to = "log2Expr", names_to = "BCR.Sample.barcode") %>%
  inner_join(metadata.AURORA, by = "BCR.Sample.barcode") %>%
  filter(Anatomic.Site.Simplified %in% c("Breast", "Brain","Lymph node", "Liver", "Lung")) %>% 
  mutate(Anatomic.Site.Simplified = factor(Anatomic.Site.Simplified, levels = c("Breast", "Lymph node", "Lung",  "Brain", "Liver")))


library(ggpubr)
set.seed(46)
unique(Aurora.long$Anatomic.Site.Simplified)
Aurora.long %>%
  filter(log2Expr > 1) %>%
  ggplot(aes(x = Anatomic.Site.Simplified, y = log2Expr, fill = Anatomic.Site.Simplified)) +
  geom_boxplot(width = 0.5, outlier.shape = T) +
  geom_jitter(size = 1.5, position = position_jitterdodge(dodge.width = 0.1, jitter.width = 0.8)) +
  theme_bw() + 
  facet_wrap(~ gene_name, scales = "free_y", nrow = 7) + 
  labs(y = "log2(Normalized counts + 1)", x = "Tumor type") +
  theme(legend.position = "none") +
  stat_compare_means(label = "p.signif",
                     comparisons = list(c("Breast", "Lymph node"),
                                        c("Breast", "Lung"),
                                        c("Breast", "Brain"),
                                        c("Breast", "Liver")))




library(tidyverse)
library(DESeq2)

# Comparaciones a realizar
sites <- c("Lung", "Brain", "Liver", "Lymph node")

# Lista para guardar resultados
deseq_results <- list()

# Loop por cada sitio metastásico
for (site in sites) {
  
  message("Comparing Breast vs ", site)
  
  # Subset de colData
  subset_colData <- colData %>%
    filter(Anatomic.Site.Simplified %in% c("Breast", site))
  
  # Subset de matriz de expresión cruda
  subset_counts <- BRCA.raw[, rownames(subset_colData)]
  
  # Crear objeto DESeq2
  dds <- DESeqDataSetFromMatrix(countData = subset_counts,
                                colData = subset_colData,
                                design = ~ Anatomic.Site.Simplified)
  
  # Normalización y filtrado
  dds <- estimateSizeFactors(dds)
  keep <- rowMeans(counts(dds)) >= 20
  dds <- dds[keep,]
  
  # Relevel para que "Breast" sea la referencia
  dds$Anatomic.Site.Simplified <- relevel(dds$Anatomic.Site.Simplified, ref = "Breast")
  
  # DESeq
  dds <- DESeq(dds)
  
  # Resultados
  res <- results(dds, contrast = c("Anatomic.Site.Simplified", site, "Breast"), pAdjustMethod = "BH") %>%
    as.data.frame() %>%
    rownames_to_column(var = "gene_name") %>%
    drop_na(pvalue, padj) %>%
    mutate(Expression = case_when(
      log2FoldChange > 1 & pvalue < 0.05 ~ "Upregulated",
      log2FoldChange < -1 & pvalue < 0.05 ~ "Downregulated",
      TRUE ~ "NS"
    )) %>%
    filter(baseMean > 20) %>%
    mutate(Expression = factor(Expression, levels = c("Downregulated", "NS", "Upregulated")))
  
  # Guardar
  deseq_results[[site]] <- res
}

#
Lung <- deseq_results$Lung %>%
  filter(gene_name %in% Genes.interest$geneNames)

Lymph <- deseq_results$`Lymph node` %>%
  filter(gene_name %in% Genes.interest$geneNames)
  
LBrain <- deseq_results$Brain %>%
  filter(gene_name %in% Genes.interest$geneNames)

Liver <- deseq_results$Liver %>%
  filter(gene_name %in% Genes.interest$geneNames)
