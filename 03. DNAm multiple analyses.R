library(tidyverse)
library(openxlsx)
library(conflicted)
library(ggpubr)
library(stringr)
conflict_prefer_all("dplyr")


## 01. Looking how samples are distributed ##
Metadata <- read_tsv("XEN metadata.txt") %>%
  mutate(Color = case_when(Organ %in% "PT" ~ "black",
                           Organ %in% "LU" ~ "red",
                           Organ %in% "LN" ~ "pink"))

myCombat <- read.table("XEN Beta-values.txt", header = T, sep = "\t")
colnames(myCombat) <- gsub("\\.", "-", colnames(myCombat))
Manifest.EPIC <- read_tsv("EPICv2 hg38 Final Manifest.tsv")

myCombat <- myCombat[, Metadata$Sample_Name]

# Apply SD to extract top 5000, 10000 and 25000 variable sites
SD <- apply(myCombat, 1, sd)
top.SD <- sort(SD, decreasing = T)[1:25000]
myCombat.SD <- myCombat[names(top.SD), ]

myCombat.M2 <- myCombat.SD %>%
  select(Metadata$Sample_Name[Metadata$Animal %in% "M2"])
mean.M2 <- apply(myCombat.M2, 1 , mean)
myCombat.M2 <- myCombat.M2/mean.M2


myCombat.M5 <- myCombat.SD %>%
  select(Metadata$Sample_Name[Metadata$Animal %in% "M5"])
mean.M5 <- apply(myCombat.M5, 1 , mean)
myCombat.M5 <- myCombat.M5/mean.M5


myCombat.M6 <- myCombat.SD %>%
  select(Metadata$Sample_Name[Metadata$Animal %in% "M6"])
mean.M6 <- apply(myCombat.M6, 1 , mean)
myCombat.M6 <- myCombat.M6/mean.M6


myCombat.M7 <- myCombat.SD %>%
  select(Metadata$Sample_Name[Metadata$Animal %in% "M7"])
mean.M7 <- apply(myCombat.M7, 1 , mean)
myCombat.M7 <- myCombat.M7/mean.M7



myCombat.normalized <- cbind(myCombat.M2, myCombat.M5, myCombat.M6, myCombat.M7)


# Pretty heatmap
library(pheatmap)
my_sample_col <- Metadata %>%
  column_to_rownames(var = "Sample_Name") %>%
  select(Organ, Animal, Array_number)
pheatmap(myCombat.normalized, 
         annotation_col = my_sample_col,
         scale = "row",
         show_rownames = F,
         cluster_rows = T,
         treeheight_row = 0)

# FIGTREE
Fig.Tree <- all_probes.LUvsLNPT %>%
  filter(probeID %in% names(top.SD)) %>%
  select(LU.mean, LN.mean, PT.mean)

Fig.Tree <- t(Fig.Tree)
rownames(Fig.Tree) <- c("Lung Metastases", "Lymph Node Metastases", "Primary Tumor")

dist_matrix <- dist(Fig.Tree, method = "euclidean")
hc <- hclust(dist_matrix, method = "average")  # UPGMA

library(ape)
tree <- as.phylo(hc)
write.tree(tree, file = "tree.nwk")



# Cargamos los datos de All probes
all_probes.LU.PT <- read_tsv("All_probes.LU.PT.txt") %>%
  mutate(ProbeV1 = gsub("_.*", "", probeID)) %>%
  inner_join(Manifest.EPIC, by = "probeID") 

all_probes.LN.PT <- read_tsv("All_probes.LN.PT.txt") %>%
  mutate(ProbeV1 = gsub("_.*", "", probeID)) %>%
  inner_join(Manifest.EPIC, by = "probeID") 



# DMSs for each

# 1 - LN vs PT #
Metadata.LN.PT <- Metadata %>%
  filter(Organ %in% c("LN", "PT"))

all_probes.LN.PT <- all_probes.LN.PT %>%
  mutate(Meth_status = case_when(fold.LN.PT > 0.15 & pvalue < 0.1 ~ "Hypermethylated",
                                 fold.LN.PT < -0.15 & pvalue < 0.1 ~ "Hypomethylated",
                                 T ~ "NS")) %>%
  mutate(Meth_status = factor(Meth_status, levels = c("Hypomethylated", "NS", "Hypermethylated"))) %>%
  filter(!duplicated(probeID))

table(all_probes.LN.PT$Meth_status)

probes_with_changes.LN.PT <- all_probes.LN.PT %>%
  filter(Meth_status %in% c("Hypermethylated", "Hypomethylated"))

LN.PT.hyper <- probes_with_changes.LN.PT %>%
  filter(Meth_status %in% "Hypermethylated")
LN.PT.hypo <- probes_with_changes.LN.PT %>%
  filter(Meth_status %in% "Hypomethylated")



# 2 - LU vs PT #
Metadata.LU.PT <- Metadata %>%
  filter(Organ %in% c("LU", "PT"))
all_probes.LU.PT <- all_probes.LU.PT %>%
  mutate(Meth_status = case_when(fold.LU.PT > 0.15 & pvalue < 0.1 ~ "Hypermethylated",
                                 fold.LU.PT < -0.15 & pvalue < 0.1 ~ "Hypomethylated",
                                 T ~ "NS")) %>%
  mutate(Meth_status = factor(Meth_status, levels = c("Hypomethylated", "NS", "Hypermethylated")))  %>%
  filter(!duplicated(probeID))
table(all_probes.LU.PT$Meth_status)

probes_with_changes.LU.PT <- all_probes.LU.PT %>%
  filter(Meth_status %in% c("Hypermethylated", "Hypomethylated"))


LU.PT.hyper <- probes_with_changes.LU.PT %>%
  filter(Meth_status %in% "Hypermethylated")
LU.PT.hypo <- probes_with_changes.LU.PT %>%
  filter(Meth_status %in% "Hypomethylated")


# 3 - LU vs LN
all_probes.LU.LN <- read_tsv("All_probes.LN.LU.txt") %>%
  rowwise() %>%
  mutate(fold.LU.LN = LU.mean - LN.mean) %>%
  select(-fold.LN.LU)  %>%
  mutate(ProbeV1 = gsub("_.*", "", probeID)) %>%
  inner_join(Manifest.EPIC, by = "probeID") %>%
  filter(!duplicated(probeID))


all_probes.LU.LN <- all_probes.LU.LN %>%
  mutate(Meth_status = case_when(fold.LU.LN > 0.15 & pvalue < 0.1 ~ "Hypermethylated",
                                 fold.LU.LN < -0.15 & pvalue < 0.1 ~ "Hypomethylated",
                                 T ~ "NS")) %>%
  mutate(Meth_status = factor(Meth_status, levels = c("Hypomethylated", "NS", "Hypermethylated"))) 

table(all_probes.LU.LN$Meth_status)

probes_with_changes.LU.LN <- all_probes.LU.LN %>%
  filter(Meth_status %in% c("Hypermethylated", "Hypomethylated"))

LU.LN.hyper <- probes_with_changes.LU.LN %>%
  filter(Meth_status %in% "Hypermethylated")
LU.LN.hypo <- probes_with_changes.LU.LN %>%
  filter(Meth_status %in% "Hypomethylated")


#
Barplot <- data.frame(Condition = c("LU.PT", "LU.PT", "LN.PT", "LN.PT", "LU.LN", "LU.LN"),
                      Methylation = c("Hypermethylated", "Hypomethylated", "Hypermethylated", 
                                      "Hypomethylated", "Hypermethylated", "Hypomethylated"),
                      DMS = c(nrow(LU.PT.hyper), nrow(LU.PT.hypo), nrow(LN.PT.hyper), nrow(LN.PT.hypo),
                              nrow(LU.LN.hyper), nrow(LU.LN.hypo))) %>%
  mutate(Condition = factor(Condition, levels = c("LU.PT", "LU.LN", "LN.PT")))


ggplot(Barplot, aes(x = Condition, y = DMS, fill = Methylation))  +
  geom_col(position = "dodge", width = 0.8) +
  theme_bw()




## 02. LU vs LN and PT together analyses #

all_probes.LUvsLNPT <- read_tsv("All_probes.LUvsLN+PT.txt")

all_probes.LUvsLNPT <- all_probes.LUvsLNPT %>%
  mutate(Meth_status = case_when(fold.LU.LN_PT > 0.15 & pvalue < 0.05 ~ "Hypermethylated",
                                 fold.LU.LN_PT < -0.15 & pvalue < 0.05 ~ "Hypomethylated",
                                 T ~ "NS")) %>%
  mutate(Meth_status = factor(Meth_status, levels = c("Hypomethylated", "NS", "Hypermethylated"))) %>%
  #inner_join(Manifest.EPIC, by = "probeID") %>%
  filter(!duplicated(probeID)) %>%
  mutate(ProbeV1 = gsub("_.*", "", probeID))


table(all_probes.LUvsLNPT$Meth_status)


probes_with_changes.LUvsLNPT <- all_probes.LUvsLNPT %>%
  filter(Meth_status %in% c("Hypomethylated", "Hypermethylated"))
myCombat.DMS <- myCombat[probes_with_changes.LUvsLNPT$probeID, ]


# 
hyper.LUvsLNPT <- probes_with_changes.LUvsLNPT %>%
  filter(Meth_status %in% c("Hypermethylated"))

hypo.LUvsLNPT <- probes_with_changes.LUvsLNPT %>%
  filter(Meth_status %in% c("Hypomethylated"))


# Circus
library(tidyverse)
library(circlize)

Hyper.circ <- hyper.LUvsLNPT %>%
  select(CpG_chrm, CpG_beg, CpG_end) %>%
  filter(!duplicated(paste(CpG_chrm,CpG_beg, CpG_end)))

Hypo.circ <- hypo.LUvsLNPT %>%
  select(CpG_chrm, CpG_beg, CpG_end) %>%
  filter(!duplicated(paste(CpG_chrm,CpG_beg, CpG_end)))


circos.clear()
circos.par("start.degree" = 90)
circos.initializeWithIdeogram(species = "hg38", chromosome.index = paste0("chr", c(1:22, "X")))
circos.genomicDensity(Hypo.circ, col = "limegreen", track.height = 0.15)
circos.genomicDensity(Hyper.circ, col = "red",track.height = 0.15)





# 03. Characterization of DMS #


# Calculating Odds Ratio CpG location reference to genome
# Hypermethylated
library(questionr)

table(all_probes.LUvsLNPT$Meth_status, all_probes.LUvsLNPT$CGIposition)

all_probes.LUvsLNPT$CGIposition <- factor(all_probes.LUvsLNPT$CGIposition, levels = c("Island", "Shore", "Shelf", "OpenSea"))

nrow(all_probes.LUvsLNPT %>% filter(CGIposition %in% "Shelf", 
                                    Meth_status %in% "Hypermethylated"))

hyper.island <- matrix(c(415, 4209, 148307, 715896),
                       nrow = 2, byrow = TRUE)

hyper.shore <- matrix(c(597,4027, 151269,712934), 
                      nrow = 2, byrow = TRUE)

hyper.shelf <- matrix(c(298, 4326, 57673, 806530), 
                      nrow = 2, byrow = TRUE)

hyper.opensea <- matrix(c(3314, 1310, 506954, 357249), 
                        nrow = 2, byrow = TRUE)

hyper.island.od <- odds.ratio(hyper.island) %>%
  mutate(Type = "Island")
hyper.shore.od <- odds.ratio(hyper.shore) %>%
  mutate(Type = "Shore")
hyper.shelf.od <- odds.ratio(hyper.shelf)%>%
  mutate(Type = "Shelf")
hyper.opensea.od <- odds.ratio(hyper.opensea)%>%
  mutate(Type = "OpenSea")

hyper.all.od <- rbind(hyper.island.od,hyper.shore.od,hyper.shelf.od,hyper.opensea.od)

library(forestplot)
library(tidyverse)
hyper.all.od$Type <- factor(hyper.all.od$Type, levels = c("OpenSea", "Shelf", "Shore", "Island"))

hyper.all.od %>%
  arrange(Type) %>%
  ggplot(aes(x = OR, y = as.factor(Type))) +
  geom_point(size = 3.5, color = "red") +
  geom_vline(aes(xintercept = 1), linewidth = .25, linetype = "dashed") + 
  geom_errorbarh(aes(xmax =`97.5 %` , xmin = `2.5 %`), size = 0.5, height = .2, color = "black") + 
  theme_bw() +
  geom_text(aes(label = paste(round(OR, 2), 
                              "[" , round(`2.5 %`, 2), 
                              " - ", round(`97.5 %`,2), 
                              "]")),
            vjust = -3, size = 3, color = "black") +
  geom_text(aes(label = paste("P =", round(p, 5))),
            vjust = -1.5, size = 3, color = "black") + 
  ylab("Relation To Genome Position") + 
  xlim(0,2)



# Hypomethylated
nrow(all_probes.LUvsLNPT %>% filter(CGIposition %in% "Island", 
                                    Meth_status %in% "Hypomethylated"))



hypo.island <- matrix(c(4563, 31017, 144159, 689088), 
                      nrow = 2, byrow = TRUE)

hypo.shore <- matrix(c(5001, 30579, 146865, 686382), 
                     nrow = 2, byrow = TRUE)

hypo.shelf <- matrix(c(2764, 32816, 55207, 778040), 
                     nrow = 2, byrow = TRUE)

hypo.opensea <- matrix(c(23252, 12328, 487016, 346231), 
                       nrow = 2, byrow = TRUE)


hypo.island.od <- odds.ratio(hypo.island) %>%
  mutate(Type = "Island")
hypo.shore.od <- odds.ratio(hypo.shore) %>%
  mutate(Type = "Shore")
hypo.shelf.od <- odds.ratio(hypo.shelf)%>%
  mutate(Type = "Shelf")
hypo.opensea.od <- odds.ratio(hypo.opensea)%>%
  mutate(Type = "OpenSea")

hypo.all.od <- rbind(hypo.island.od,hypo.shore.od,hypo.shelf.od,hypo.opensea.od)

library(forestplot)
library(tidyverse)
hypo.all.od$Type <- factor(hypo.all.od$Type, levels = c("OpenSea", "Shelf", "Shore", "Island"))

hypo.all.od %>%
  arrange(Type) %>%
  ggplot(aes(x = OR, y = as.factor(Type))) +
  geom_point(size = 3.5, color = "green") +
  geom_vline(aes(xintercept = 1), linewidth = .25, linetype = "dashed") + 
  geom_errorbarh(aes(xmax =`97.5 %` , xmin = `2.5 %`), size = 0.5, height = .2, color = "black") + 
  theme_bw() +
  geom_text(aes(label = paste(round(OR, 2), 
                              "[" , round(`2.5 %`, 2), 
                              " - ", round(`97.5 %`,2), 
                              "]")),
            vjust = -3, size = 3, color = "black") +
  geom_text(aes(label = paste("P =", round(p, 2))),
            vjust = -1.5, size = 3, color = "black") + 
  ylab("Relation To Genome Position") + 
  xlim(0,2)










# DMS in Promoter Regions
RNA.anot <- read_tsv("Gene_metadata_annotations.txt")
RNA.anot <- RNA.anot %>%
  filter(gene_type %in% "protein_coding")

Meth.Mice.LU.LNandPT.filt <- all_probes.LUvsLNPT %>%
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




# CANCER HALLMARKS INTEGRATED #
Cancer_hallmarks <- read_tsv("Menyhart_JPA_CancerHallmarks_integrated.txt")

Cancer_hallmarks.curated <- Cancer_hallmarks %>%
  mutate(ID = rownames(Cancer_hallmarks)) %>%
  pivot_longer(-ID, names_to = "gs_name", values_to = "gene_symbol") %>%
  select(-ID) %>%
  arrange(gs_name) %>%
  drop_na(gene_symbol)


table(Cancer_hallmarks.curated$gs_name)


library(clusterProfiler)
Cancer.hallmark.MM <- enricher(gene = Meth.Mice.LU.LNandPT.count$geneNames,
                               TERM2GENE =  Cancer_hallmarks.curated,
                               pAdjustMethod = "BH",
                               pvalueCutoff = 1,
                               minGSSize = 1,
                               maxGSSize = 1000,
                               qvalueCutoff = 1)
Cancer.hallmark.MM.df <- as.data.frame(Cancer.hallmark.MM@result)



plot.df.our <- Cancer.hallmark.MM.df %>%
  select(Description, p.adjust) %>%
  mutate(p.adjust2 = case_when(p.adjust < 0.0001 ~ 0.00005,
                               T ~ p.adjust))


# Plot
plt <- ggplot(plot.df.our) +
  geom_hline(
    aes(yintercept = y), 
    data.frame(y = -log10(c(0.1, 0.01, 0.001,0.0001, 0.00005))),
    color = "lightgrey") + 
  geom_hline(
    aes(yintercept = z), 
    data.frame(z = -log10(c(0.05))),
    color = "red",
    linetype = "dashed") +
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
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))
plt




# KEGG ENRICHMENT RESULTS #
library(msigdbr)
KEGG <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>%
  select(gs_name, gene_symbol) %>%
  filter(!duplicated(paste(gs_name, gene_symbol)))


#
KEGG.MM <- enricher(gene = Meth.Mice.LU.LNandPT.count$geneNames,
                        TERM2GENE = KEGG,
                        pAdjustMethod = "BH",
                        pvalueCutoff = 1,
                        minGSSize = 1,
                        maxGSSize = 1000,
                        qvalueCutoff = 1)
KEGG.MM.df <- as.data.frame(KEGG.MM@result) %>%
  filter(pvalue < 0.05)


KEGG.MM.df %>%
  arrange(pvalue) %>%
  slice(1:10) %>%
  mutate(Description = fct_reorder(Description, -log10(pvalue))) %>%
  ggplot(aes(x = Description, y = -log10(pvalue))) +
  geom_bar(stat = "identity", width = 0.8, fill = "tomato") +
  coord_flip() +
  theme_bw() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")

