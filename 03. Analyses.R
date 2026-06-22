library(tidyverse)
library(conflicted)
library(openxlsx)
library(data.table)
library(ggpubr)
conflict_prefer_all("dplyr")

setwd("C:/Users/idisb/Desktop/Andrés/02 - Projects/05 - DNAm organ specific")


## FIGTREE
myCombat <- read.table("AURORA Filtered Beta-Values corrected.txt", sep = "\t")
SD <- apply(myCombat, 1, sd)
top.SD <- sort(SD, decreasing = T)[1:5000]
myCombat.SD <- myCombat[names(top.SD), ]

Fig.Tree <- all_probes.LU.PT %>%
  inner_join(all_probes.LN.PT, by = "probeID") %>%
  filter(probeID %in% names(top.SD)) %>%
  select(LU.mean, LN.mean, PT.mean.x)

Fig.Tree <- t(Fig.Tree)
rownames(Fig.Tree) <- c("Lung Metastases", "Lymph Node Metastases", "Primary Tumor")

dist_matrix <- dist(Fig.Tree, method = "euclidean")
hc <- hclust(dist_matrix, method = "average")  # UPGMA
library(ape)
tree <- as.phylo(hc)
write.tree(tree, file = "AURORA tree.nwk")




#
all_probes.LU.PT <- read_tsv("AURORA All_probes.LU.PT.txt") 

all_probes.LU.PT <- all_probes.LU.PT %>%
  mutate(Meth_status = case_when(pvalue < 0.05 & fold.LU.PT > 0.15 ~ "Hypermethylated",
                                 pvalue < 0.05  & fold.LU.PT < -0.15 ~ "Hypomethylated",
                                 T ~ "NS"))

LU.PT.hyper <- all_probes.LU.PT %>%
  filter(Meth_status %in% "Hypermethylated") %>%
  filter(pvalue < 0.05)
LU.PT.hypo <- all_probes.LU.PT %>%
  filter(Meth_status %in% "Hypomethylated") %>%
  filter(pvalue < 0.05)

#
all_probes.LN.PT <- read_tsv("AURORA All_probes.LN.PT.txt")

all_probes.LN.PT <- all_probes.LN.PT %>%
  mutate(Meth_status = case_when(pvalue < 0.05 & fold.LN.PT > 0.15 ~ "Hypermethylated",
                                 pvalue < 0.05  & fold.LN.PT < -0.15 ~ "Hypomethylated",
                                 T ~ "NS"))

LN.PT.hyper <- all_probes.LN.PT %>%
  filter(Meth_status %in% "Hypermethylated")

LN.PT.hypo <- all_probes.LN.PT %>%
  filter(Meth_status %in% "Hypomethylated")

#
all_probes.LU.LN <- read_tsv("AURORA All_probes.LU.LN.txt")

all_probes.LU.LN <- all_probes.LU.LN %>%
  mutate(Meth_status = case_when(pvalue < 0.05 & fold.LU.LN > 0.15 ~ "Hypermethylated",
                                 pvalue < 0.05  & fold.LU.LN < -0.15 ~ "Hypomethylated",
                                 T ~ "NS"))


LU.LN.hyper <- all_probes.LU.LN %>%
  filter(Meth_status %in% "Hypermethylated") 

LU.LN.hypo <- all_probes.LU.LN %>%
  filter(Meth_status %in% "Hypomethylated")


Barplot <- data.frame(Condition = c("LU.PT", "LU.PT", "LN.PT", "LN.PT", "LU.LN", "LU.LN"),
                      Methylation = c("Hypermethylated", "Hypomethylated", "Hypermethylated", 
                                      "Hypomethylated", "Hypermethylated", "Hypomethylated"),
                      DMS = c(nrow(LU.PT.hyper), nrow(LU.PT.hypo), nrow(LN.PT.hyper), nrow(LN.PT.hypo),
                              nrow(LU.LN.hyper), nrow(LU.LN.hypo))) %>%
  mutate(Condition = factor(Condition, levels = c("LU.PT", "LN.PT", "LU.LN")))


ggplot(Barplot, aes(x = Condition, y = DMS, fill = Methylation))  +
  geom_col(position = "dodge", width = 0.8) +
  theme_bw()


###############
RNA.anot <- read_tsv("Gene_metadata_annotations.txt")
RNA.anot <- RNA.anot %>%
  filter(gene_type %in% "protein_coding")

AURORA.all.probes.LU.PT <- all_probes.LU.PT %>%
  filter(pvalue < 0.05,
         abs(fold.LU.PT) > 0.15) %>%
  separate_rows(geneNames, distToTSS, transcriptTypes, sep = ";") %>%
  filter(geneNames %in% RNA.anot$gene_name)  %>%
  mutate(distToTSS = as.numeric(distToTSS)) %>%
  filter(transcriptTypes %in% "protein_coding") %>%
  filter(between(distToTSS, -2000, 0)) %>%
  filter(!duplicated(paste(geneNames, probeID)))

AURORA.count.per.gene.LU <- AURORA.all.probes.LU.PT %>%
  group_by(geneNames) %>%
  summarize(total = n(),
            status = ifelse(all(fold.LU.PT > 0), "Hypermethylated",
                            ifelse(all(fold.LU.PT < 0), "Hypomethylated", "Mixed"))) %>%
  arrange(desc(total)) %>%
  filter(total > 1)


table(AURORA.count.per.gene.LU$status)



### CANCER HALLMARKS INTEGRATED ###
Cancer_hallmarks <- read_tsv("Menyhart_JPA_CancerHallmarks_integrated.txt")

Cancer_hallmarks.curated <- Cancer_hallmarks %>%
  mutate(ID = rownames(Cancer_hallmarks)) %>%
  pivot_longer(-ID, names_to = "gs_name", values_to = "gene_symbol") %>%
  select(-ID) %>%
  arrange(gs_name) %>%
  drop_na(gene_symbol)


table(Cancer_hallmarks.curated$gs_name)

library(clusterProfiler)
Cancer.hallmark.LU <- enricher(gene = AURORA.count.per.gene.LU$geneNames,
                                 TERM2GENE =  Cancer_hallmarks.curated,
                                 pAdjustMethod = "BH",
                                 pvalueCutoff = 1,
                                 minGSSize = 1,
                                 maxGSSize = 1000,
                                 qvalueCutoff = 1)
Cancer.hallmark.LU.df <- as.data.frame(Cancer.hallmark.LU@result)





library(dplyr)
library(ggplot2)
library(stringr)


plot.df.our <- Cancer.hallmark.LU.df %>%
  select(Description, p.adjust) %>%
  mutate(p.adjust2 = case_when(p.adjust < 0.0001 ~ 0.00005,
                               T ~ p.adjust))


# Plot
plt <- ggplot(plot.df.our) +
  geom_hline(
    aes(yintercept = y), 
    data.frame(y = -log10(c(0.1, 0.01, 0.001,0.0001, 0.00005))),
    color = "lightgrey"
  ) + 
  geom_hline(
    aes(yintercept = z), 
    data.frame(z = -log10(c(0.05))),
    color = "red",
    linetype = "dashed"
  ) +
  geom_col(
    aes(
      x = Description,
      y = -log10(p.adjust2)),
    position = "dodge2",
    fill = "steelblue",
    color = "black",
    show.legend = FALSE) +
  geom_segment(
    aes(
      x = Description,
      y = 0,
      xend = Description,
      yend = -log10(0.0001)
    ),
    linetype = "dashed",
    color = "gray12"
  ) +
  # Make it circular with theme adjustments
  coord_polar() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      size = 8, 
      angle = 0, 
      hjust = 1, 
      vjust = 1),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    plot.margin = margin(10, 10, 10, 10)) +
  scale_y_continuous(
    limits = c(-1, 5),
    expand = c(0, 0),
    breaks = c(0, -log10(0.1), -log10(0.01), -log10(0.001), -log10(0.0001))) +
  # Wrap long descriptions
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))
plt



#### KEGG ENRICHMENT RESULTS ####
library(msigdbr)
KEGG <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>%
  select(gs_name, gene_symbol) %>%
  filter(!duplicated(paste(gs_name, gene_symbol)))


#
KEGG.AUR <- enricher(gene = AURORA.count.per.gene.LU$geneNames,
                    TERM2GENE = KEGG,
                    pAdjustMethod = "BH",
                    pvalueCutoff = 1,
                    minGSSize = 1,
                    maxGSSize = 1000,
                    qvalueCutoff = 1)
KEGG.AUR.DF <- as.data.frame(KEGG.AUR@result) %>%
  filter(pvalue < 0.05)


KEGG.AUR.DF %>%
  arrange(pvalue) %>%
  slice(1:10) %>%
  mutate(Description = fct_reorder(Description, -log10(pvalue))) %>%
  ggplot(aes(x = Description, y = -log10(pvalue))) +
  geom_bar(stat = "identity", width = 0.8, fill = "tomato") +
  coord_flip() +
  theme_bw() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")







#### NOW INTERSECT GENES OUR MM AND AURORA

# OUR cohort
Meth.Mice.LU.LNandPT <- read_tsv("All_probes.LUvsLN+PT.txt") %>%
  filter(!duplicated(probeID))


Meth.Mice.LU.LNandPT.filt <- Meth.Mice.LU.LNandPT %>%
  filter(pvalue < 0.05,
         abs(fold.LU.LN_PT) > 0.15) %>%
  separate_rows(geneNames, distToTSS, transcriptTypes, sep = ";") %>%
  filter(geneNames %in% RNA.anot$gene_name)  %>%
  mutate(distToTSS = as.numeric(distToTSS)) %>%
  filter(transcriptTypes %in% "protein_coding") %>%
  filter(between(distToTSS, -2000, 0)) %>%
  filter(!duplicated(paste(geneNames, probeID)))

Meth.Mice.LU.LNandPT.count <- Meth.Mice.LU.LNandPT.filt %>%
  group_by(geneNames) %>%
  summarize(total_MM = n(),
            status_MM = ifelse(all(fold.LU.LN_PT > 0), "Hypermethylated",
                               ifelse(all(fold.LU.LN_PT < 0), "Hypomethylated", "Mixed"))) %>%
  arrange(desc(total_MM)) %>%
  filter(total_MM > 1) 


table(Meth.Mice.LU.LNandPT.count$status_MM)


## Intersect both
Intersect.DNAm.all <-  Meth.Mice.LU.LNandPT.count%>%
  inner_join(AURORA.count.per.gene.LU, by = c("geneNames" = "geneNames")) %>%
  filter(total_MM > 1 & total > 1)


Intersect.DNAm.all.filt <- Intersect.DNAm.all %>%
  filter(!(status_MM == "Hypermethylated" & status == "Hypomethylated"),
         !(status_MM == "Hypomethylated" & status == "Hypermethylated"))




library(ggvenn)

genes_mice <- unique(Meth.Mice.LU.LNandPT.count$geneNames)
genes_aurora <- unique(AURORA.count.per.gene.LU$geneNames)
genes_filtered <- unique(Intersect.DNAm.all.filt$geneNames)

gene_lists <- list(
  Mice = genes_mice,
  AURORA = genes_aurora)

ggvenn(gene_lists,
       fill_color = c("#E69F00", "#56B4E9"),
       stroke_size = 0.5,
       set_name_size = 5)


gene_lists <- list(
  Mice = genes_mice,
  AURORA = genes_aurora,
  Filtered_Intersect = genes_filtered)

ggvenn(gene_lists,
       fill_color = c("#E69F00", "#56B4E9", "#009E73"),
       stroke_size = 0.5,
       set_name_size = 5)







# After RNA-seq data analysis...
AURORA.RNA.LU <- read_tsv("AURORA RNA-seq LU vs PT.txt") %>%
  mutate(Expression = case_when(log2FoldChange > 1 & pvalue < 0.05 ~ "Upregulated",
                                log2FoldChange < -1 & pvalue < 0.05 ~ "Downregulated",
                                T ~ "NS")) %>%
  filter(baseMean > 20)



library(ggpubr)
ggplot(AURORA.RNA.LU, aes(x = log2FoldChange, y = -log10(pvalue), color = Expression)) + 
  geom_point(size = 1.5) + 
  theme_pubr() +
  xlim(-10, 10) +
  scale_color_manual(values=c("blue", "grey", "red")) +
  geom_vline(xintercept = 1, color = "black", linetype = "dashed", linewidth = 0.5) + 
  geom_vline(xintercept = -1, color = "black", linetype = "dashed", linewidth = 0.5) + 
  geom_hline(yintercept = -log10(0.05), color = "black", linetype = "dashed", linewidth = 0.5)



AURORA.RNA.LU.intersect <- Intersect.DNAm.all.filt %>%
  inner_join(AURORA.RNA.LU, by = c("geneNames" = "gene_name")) 


AURORA.RNA.LU.intersect.filt <- AURORA.RNA.LU.intersect %>%
  filter(pvalue < 0.05) %>%
  filter(status_MM %in% c("Hypomethylated", "Mixed") & status %in% c("Hypomethylated", "Mixed") & Expression %in% "Upregulated"|
         status_MM %in% c("Hypermethylated", "Mixed") & status %in% c("Hypermethylated", "Mixed") & Expression %in% "Downregulated")

table(AURORA.RNA.LU.intersect.filt$Expression)




library(msigdbr)
KEGG <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>%
  select(gs_name, gene_symbol) %>%
  filter(!duplicated(paste(gs_name, gene_symbol)))

library(UpSetR)


gene_sets <- list(
  DNAm.MM = Meth.Mice.LU.LNandPT.count$geneNames,
  DNAm.AUR = Intersect.DNAm.all.filt$geneNames,
  AUR.RNA.LU = AURORA.RNA.LU.intersect.filt$geneNames)

all_genes <- unique(unlist(gene_sets))
upset_data <- data.frame(Gene = all_genes)

for (set_name in names(gene_sets)) {
  upset_data[[set_name]] <- ifelse(upset_data$Gene %in% gene_sets[[set_name]], 1, 0)
}

upset(upset_data, 
      sets = names(gene_sets), 
      sets.bar.color = "steelblue", 
      keep.order = TRUE)


Winners <- upset_data %>%
  filter(if_all(-Gene, ~ . == 1))








######## Figure 4

Interest.genes.4B.MM <- Meth.Mice.LU.LNandPT %>%
  filter(ProbeV1 %in% all_probes.LU.PT$probeID) %>%
  separate_rows(geneNames, distToTSS, transcriptTypes, sep = ";") %>%
  mutate(distToTSS = as.numeric(distToTSS)) %>%
  filter(transcriptTypes %in% "protein_coding") %>%
  filter(between(distToTSS, -2000, 0)) %>%
  filter(!duplicated(paste(probeID, geneNames))) %>%
  filter(geneNames %in% Winners$Gene) %>%
  arrange(CpG_beg) 


Interest.genes.4B.MM <- Interest.genes.4B.MM %>%
  group_by(geneNames) %>%
  mutate(pos = row_number()) %>%
  mutate(pos = as.numeric(pos)) %>%
  filter(geneNames %in% c("SLC2A5", "CD5", "HOXA5"))


Fig4b.1 <- ggplot(Interest.genes.4B.MM, aes(x = pos)) +
  geom_area(aes(y = LU.mean, fill = "LU"), alpha = 0.3) +
  geom_area(aes(y = LN.PT.mean, fill = "LN and PT"), alpha = 0.3) +
  theme_bw()  +
  scale_x_continuous(
    name = "",
    breaks = Interest.genes.4B.MM$pos,
    labels = Interest.genes.4B.MM$ProbeV1
  )  +
  labs(y = "Methylation B-value", x = "") +
  facet_wrap(~geneNames, scales = "free_x", nrow = 3) +
  scale_fill_manual(values = c("LU" = "blue", "LN and PT" = "black")) +
  theme(legend.title = element_blank(),
        legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0, size = 6)) +
  ylim(0, 1)




#AUR
Interest.genes.4B.AUR <- all_probes.LU.PT %>%

  filter(probeID %in% Interest.genes.4B.MM$ProbeV1) %>%
  separate_rows(geneNames, distToTSS, transcriptTypes, sep = ";") %>%
  filter(geneNames %in% RNA.anot$gene_name)  %>%
  mutate(distToTSS = as.numeric(distToTSS)) %>%
  filter(transcriptTypes %in% "protein_coding") %>%
  filter(between(distToTSS, -2000, 0)) %>% 
  filter(!duplicated(paste(geneNames, probeID)))  %>%
  filter(geneNames %in% Winners$Gene) %>%
  arrange(CpGbeg)

table(Interest.genes.4B.AUR$geneNames, Interest.genes.4B.AUR$MyMeth)


Interest.genes.4B.AUR <- Interest.genes.4B.AUR %>%
  group_by(geneNames) %>%
  mutate(pos = row_number()) %>%
  mutate(pos = as.numeric(pos)) %>%
  filter(geneNames %in% c("SLC2A5", "CD5", "HOXA5"))


library(ggpubr)
Fig4b.2 <- ggplot(Interest.genes.4B.AUR, aes(x = pos)) +
  geom_area(aes(y = LU.mean, fill = "LU"), alpha = 0.3) +
  geom_area(aes(y = PT.mean, fill = "PT"), alpha = 0.3) +
  theme_bw() +
  scale_x_continuous(
    name = "",
    breaks = Interest.genes.4B.AUR$pos,
    labels = Interest.genes.4B.AUR$probeID
  )  +
  labs(y = "", x = "") +
  facet_wrap(~geneNames, scales = "free_x", nrow = 3) +
  scale_fill_manual(values = c("LU" = "blue", "PT" = "black")) +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0, size = 6)) +
  ylim(0, 1)



## ## Gene expression 4C
AURORA.Matrix <- read_tsv("Normalized counts AURORA Lung and PT.tsv")
metadata.AURORA <- read_tsv("Metadata AURORA for plot RNA.tsv")

AURORA.Matrix <- AURORA.Matrix %>%
  filter(!duplicated(gene_name)) %>%
  column_to_rownames(var = "gene_name")


colnames(AURORA.Matrix) <- case_when(
  grepl("TTM10", colnames(AURORA.Matrix)) ~ substr(colnames(AURORA.Matrix), 1, 16),
  grepl("NT", colnames(AURORA.Matrix)) ~ substr(colnames(AURORA.Matrix), 1, 14),
  TRUE ~ substr(colnames(AURORA.Matrix), 1, 15))

setdiff(colnames(AURORA.Matrix), metadata.AURORA$BCR.Sample.barcode)

AURORA.Matrix <- AURORA.Matrix[, metadata.AURORA$BCR.Sample.barcode]


##
metadata.metastasis <- metadata.AURORA %>%
  filter(Sample.Type %in% "Metastasis")

metadata.primary <- metadata.AURORA %>%
  filter(Sample.Type %in% "Primary")

table(metadata.metastasis$Anatomic.Site.Simplified)


AURORA.Matrix.log <- log2(AURORA.Matrix + 1)


Genes.interest <- AURORA.Matrix.log %>%
  rownames_to_column(var = "gene_name") %>%
  filter(gene_name %in% c("SLC2A5", "CD5", "HOXA5")) %>%
  pivot_longer(-gene_name, values_to = "log2Expr", names_to = "BCR.Sample.barcode") %>%
  inner_join(metadata.AURORA, by = "BCR.Sample.barcode") %>%
  filter(Anatomic.Site.Simplified %in% c("Breast", "Lung"))



library(ggpubr)
set.seed(46)
Fig4b.3 <- ggplot(Genes.interest, aes(x = Anatomic.Site.Simplified, y = log2Expr, fill = Anatomic.Site.Simplified)) +
  geom_boxplot(width = 0.4, outlier.shape = T) +
  geom_jitter(size = 1.5, position = position_jitterdodge(dodge.width = 0, jitter.width = 0.1)) +
  theme_bw() + 
  facet_wrap(~ gene_name, nrow = 3, scales = "free_y") +
  stat_compare_means(comparisons = list(c("Breast", "Lung")),
                     method = "t.test", label = "p.signif") + 
  labs(y = "log2(Normalized counts + 1)", x = "Tumor type") +
  theme(legend.position = "none")


library(patchwork)

Fig4b.1 | Fig4b.2 | Fig4b.3 +
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")


