library(tidyverse)
library(conflicted)
library(openxlsx)
library(data.table)
conflict_prefer_all("dplyr")

Metadata <- read.xlsx("Curated Metadata AURORA USA cohort.xlsx")
Manifest.EPIC <- read.xlsx("Manifest EPIC h38 more annotation.xlsx")

table(Metadata$Anatomic.Site.Simplified)
table(Metadata$Sample.Type)


myCombat <- fread("AURORA Filtered Beta-Values corrected.txt")
myCombat <- myCombat %>%
  column_to_rownames(var = "V1")

myCombat <- myCombat[, Metadata$BCR.Sample.barcode]

### I WILL CHECK COMPARISON
Metadata.PT <- Metadata %>%
  filter(Anatomic.Site.Simplified %in% "Breast")
Metadata.LN <- Metadata %>%
  filter(Anatomic.Site.Simplified %in% "Lymph node")
Metadata.LU <- Metadata %>%
  filter(Anatomic.Site.Simplified %in% "Lung")

## Datas and means
data.PT <- myCombat %>%
  select(all_of(Metadata.PT$BCR.Sample.barcode))
data.LN <- myCombat %>%
  select(all_of(Metadata.LN$BCR.Sample.barcode))
data.LU <- myCombat %>%
  select(all_of(Metadata.LU$BCR.Sample.barcode))

PT.mean <- apply(data.PT, 1, mean)
LU.mean <- apply(data.LU, 1, mean)
LN.mean <- apply(data.LN, 1, mean)


## COMPARISON 1. LU vs PT ####
fold.LU.PT <- LU.mean - PT.mean

# Compute statistical significance #
pvalue = NULL
tstat = NULL
for(i in 1 : nrow(data.PT)) {
  x = data.PT[i,]
  y = data.LU[i,]
  
  t = wilcox.test(as.numeric(x), as.numeric(y)) 
  pvalue[i] = t$p.value
  tstat[i] = t$statistic
  print(i)
}

all_probes.LU.PT <- cbind(pvalue, LU.mean, PT.mean, fold.LU.PT) %>%
  as.data.frame() %>%
  rownames_to_column(var = "probeID") %>%
  mutate(Meth_status = case_when(fold.LU.PT > 0.20 & pvalue < 0.05 ~ "Hypermethylated",
                                 fold.LU.PT < -0.20 & pvalue < 0.05 ~ "Hypomethylated",
                                 T ~ "NS")) %>%
  mutate(Meth_status = factor(Meth_status, levels = c("Hypomethylated", "NS", "Hypermethylated")))
table(all_probes.LU.PT$Meth_status)

write_tsv(all_probes.LU.PT, file = "New All_probes.LU.PT.txt")


## COMPARISON 2. LN vs PT ####
fold.LN.PT <- LN.mean - PT.mean

# Compute statistical significance #
pvalue = NULL
tstat = NULL
for(i in 1 : nrow(data.PT)) {
  x = data.PT[i,]
  y = data.LN[i,]
  
  t = wilcox.test(as.numeric(x), as.numeric(y)) 
  pvalue[i] = t$p.value
  tstat[i] = t$statistic
  print(i)
}

all_probes.LN.PT <- cbind(pvalue, LN.mean, PT.mean, fold.LN.PT) %>%
  as.data.frame() %>%
  rownames_to_column(var = "probeID") %>%
  mutate(Meth_status = case_when(fold.LN.PT > 0.20 & pvalue < 0.05 ~ "Hypermethylated",
                                 fold.LN.PT < -0.20 & pvalue < 0.05 ~ "Hypomethylated",
                                 T ~ "NS")) %>%
  mutate(Meth_status = factor(Meth_status, levels = c("Hypomethylated", "NS", "Hypermethylated")))
table(all_probes.LN.PT$Meth_status)


write_tsv(all_probes.LN.PT, file = "New All_probes.LN.PT.txt")


## COMPARISON 3. LN vs LU ####
fold.LU.LN <- LU.mean - LN.mean

# Compute statistical significance #
pvalue = NULL
tstat = NULL
for(i in 1 : nrow(data.LU)) {
  x = data.LU[i,]
  y = data.LN[i,]
  
  t = wilcox.test(as.numeric(x), as.numeric(y)) 
  pvalue[i] = t$p.value
  tstat[i] = t$statistic
  print(i)
}

all_probes.LU.LN <- cbind(pvalue, LU.mean, LN.mean, fold.LU.LN) %>%
  as.data.frame() %>%
  rownames_to_column(var = "probeID") %>%
  mutate(Meth_status = case_when(fold.LU.LN > 0.20 & pvalue < 0.05 ~ "Hypermethylated",
                                 fold.LU.LN < -0.20 & pvalue < 0.05 ~ "Hypomethylated",
                                 T ~ "NS")) %>%
  mutate(Meth_status = factor(Meth_status, levels = c("Hypomethylated", "NS", "Hypermethylated")))
table(all_probes.LU.LN$Meth_status)

write_tsv(all_probes.LU.LN, file = "New All_probes.LU.LN.txt")



