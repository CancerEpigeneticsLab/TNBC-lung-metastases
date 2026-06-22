library(tidyverse)
library(openxlsx)
library(conflicted)
conflict_prefer_all("dplyr")


# ANALYSIS PROBES
#5. DNA METHYLATION ANALYSIS POST NORMALIZATION ####
Annotation <- read_tsv("XEN metadata.txt")
myCombat <- read_tsv("XEN Beta-values.txt")
Manifest.EPIC <- read_tsv("EPICv2.hg38.manifest.gencode.v41.tsv")
Manifest.EPIC <- Manifest.EPIC %>%
  filter(probeID %in% myCombat$Probe_ID)

myCombat <- myCombat %>%
  column_to_rownames(var = "probeID")
# Identify DMS (3 organs referencing to PT; LU vs PT; LN vs PT)
Annotation.PT <- Annotation %>%
  filter(Organ %in% "PT")
Annotation.LU <- Annotation %>%
  filter(Organ %in% "LU")
Annotation.LN <- Annotation %>%
  filter(Organ %in% "LN")
  
  
## Datas and means
data.PT <- myCombat %>%
  select(all_of(Annotation.PT$Sample_Name))
data.LU <- myCombat %>%
  select(all_of(Annotation.LU$Sample_Name))
data.LN <- myCombat %>%
  select(all_of(Annotation.LN$Sample_Name))

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
  rownames_to_column(var = "probeID")

write_tsv(all_probes.LU.PT, file = "All_probes.LU.PT.txt")


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
  rownames_to_column(var = "probeID")

write_tsv(all_probes.LN.PT, file = "All_probes.LN.PT.txt")



## COMPARISON 3. LN vs LU ####
fold.LN.LU <- LN.mean - LU.mean

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

all_probes.LN.LU <- cbind(pvalue, LN.mean, LU.mean, fold.LN.LU) %>%
  as.data.frame() %>%
  rownames_to_column(var = "probeID")

write_tsv(all_probes.LN.LU, file = "All_probes.LN.LU.txt")

table(all_probes.LN.LU$pvalue)






# LU vs LN and PT

data.LN.PT <- cbind(data.LN, data.PT)
LN.PT.mean <- apply(data.LN.PT, 1, mean)

fold.LU.LN_PT <- LU.mean - LN.PT.mean

# Compute statistical significance #
pvalue = NULL
tstat = NULL
for(i in 1 : nrow(data.LU)) {
  x = data.LU[i,]
  y = data.LN.PT[i,]
  
  t = wilcox.test(as.numeric(x), as.numeric(y)) 
  pvalue[i] = t$p.value
  tstat[i] = t$statistic
  print(i)
}

all_probes.LU.LN_PT <- cbind(pvalue, fold.LU.LN_PT, LU.mean, LN.mean, PT.mean, LN.PT.mean) %>%
  as.data.frame() %>%
  rownames_to_column(var = "probeID")

table(all_probes.LU.LN_PT$pvalue)

write_tsv(all_probes.LU.LN_PT, file = "All_probes.LUvsLN+PT.txt")