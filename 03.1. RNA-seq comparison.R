library(tidyverse)
library(openxlsx)
library(DESeq2)
library(conflicted)
conflict_prefer_all("dplyr")

metadata.AURORA <- read.xlsx("Curated Metadata AURORA USA cohort.xlsx") %>%
  filter(!Anatomic.Site.Simplified %in% "Lymph node")
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

keep <- rowMeans(counts(dds)) >=20
dds <- dds[keep,]

# as reference:
dds$Anatomic.Site.Simplified <- relevel(dds$Anatomic.Site.Simplified, ref = "Breast")
dds <- DESeq(dds)


# results
res.LU.PT <- results(dds, contrast=c("Anatomic.Site.Simplified","Lung","Breast"),  pAdjustMethod="BH")


# results
res.naive.LU.PT <- as.data.frame(res.LU.PT) %>%
  rownames_to_column(var = "gene_name") %>%
  mutate(baseMean = round(baseMean, 3)) %>%
  drop_na(pvalue, padj) %>%
  mutate(Expression = case_when(log2FoldChange > 1 & pvalue < 0.05 ~ "Upregulated",
                                log2FoldChange < -1 & pvalue < 0.05 ~ "Downregulated",
                                T ~ "NS")) %>%
  filter(!grepl("RP11-", gene_name)) %>%
  filter(!grepl("RP5-", gene_name)) %>%
  filter(!grepl("RP4-", gene_name)) %>%
  filter(!grepl("RP1-", gene_name)) %>%
  filter(!grepl("RP6-", gene_name)) %>%
  filter(!grepl("RP3-", gene_name)) %>%
  filter(!grepl("RP13-", gene_name)) %>%
  mutate(Expression = factor(Expression, levels = c("Downregulated", "NS", "Upregulated"))) %>%
  filter(baseMean > 20)

table(res.naive.LU.PT$Expression)


write_tsv(res.naive.LU.PT, file = "AURORA RNA-seq LU vs PT.txt")


#
library(ggpubr)
ggplot(res.naive.LU.PT, aes(x = log2FoldChange, y = -log10(pvalue), color = Expression)) + 
  geom_point(size = 1) + 
  theme_pubr() +
  scale_color_manual(values=c("blue", "grey", "red")) +
  geom_vline(xintercept = 1, color = "black", linetype = "dashed", linewidth = 0.75) + 
  geom_vline(xintercept = -1, color = "black", linetype = "dashed", linewidth = 0.75) + 
  geom_hline(yintercept = -log10(0.05), color = "black", linetype = "dashed", linewidth = 0.75)


Genes.interest <- read.xlsx("Final Intersect MM, AURORA DNAm, RNAseq.xlsx")

res.interest <- res.naive.LU.PT %>%
  filter(gene_name %in% Genes.interest$geneNames)
