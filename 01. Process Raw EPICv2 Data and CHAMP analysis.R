library(tidyverse)
library(ChAMP)
library(openxlsx)
library(conflicted)
conflict_prefer_all("dplyr")


Annotation <- read_tsv(".../1_SampleAnnotation.txt")


#1. IMPORT RAW DATA #
set.seed(46)
testDir <- "..."
myImport <- champ.import(directory = testDir, arraytype = "EPICv2")

#2. FILTER RAW DATA
myLoad <- champ.filter(arraytype = "EPICv2", intensity = myImport$intensity,
                       Meth = myImport$Meth, UnMeth = myImport$UnMeth, 
                       detP = myImport$detP, filterDetP = T, filterXY = FALSE)

myLoad$pd$Sample_Name <- as.character(myLoad$pd$Sample_Name)
myLoad$pd$Slide <- as.factor(myLoad$pd$Slide)
myLoad$pd$Array <- as.factor(myLoad$pd$Array)
myLoad$pd$Array_number <- as.factor(myLoad$pd$Array_number)
myLoad$pd$Plate_well <- factor(myLoad$pd$Plate_well)
myLoad$pd$Animal <- as.factor(myLoad$pd$Animal)
myLoad$pd$Organ <- as.factor(myLoad$pd$Organ)

myLoad$pd <- myLoad$pd %>%
  select(-X) %>%
  mutate(Array_number = paste0("Array", Array_number))

str(myLoad$pd)

#3. NORMALIZE DATA
myNorm <- champ.norm(beta = myLoad$beta, arraytype="EPICv2", plotBMIQ = TRUE)
sapply(as.data.frame(myNorm), function(x) sum(is.na(x)))

#4. REMOVE BATCH EFFECT AND COVARIATES. COMBAT.
Covar <- champ.SVD(beta = myNorm %>% as.matrix(), pd = myLoad$pd, RGEffect = TRUE)

myCombat <- myNorm %>%
  as.data.frame() %>%
  rownames_to_column(var = "Probe_ID")

write_tsv(myCombat, file = "XEN Beta-values.txt")