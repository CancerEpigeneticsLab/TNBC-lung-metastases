library(tidyverse)
library(openxlsx)
library(conflicted)
conflict_prefer_all("dplyr")


metadata.AURORA <- read.xlsx("Metadata AURORA USA cohort.xlsx")
metadata.AURORA <- metadata.AURORA %>%
  filter(DNA.Methylation.FreezeSet.131 %in% "Sequenced") %>%
  filter(Sample.Type %in% c("Primary", "Metastasis")) %>%
  filter(
    (Inferred.ER.from.RNAseq %in% "Negative" & 
       Inferred.PR.from.RNAseq %in% "Negative" & 
       Inferred.HER2.from.RNAseq %in% "Negative") |
      (Profiled.ER %in% "Negative" & 
         Profiled.PR %in% "Negative" & 
         Profiled.HER2 %in% "Negative")) %>%
  filter(Anatomic.Site.Simplified %in% c("Breast", "Lung", "Lymph node"))


table(metadata.AURORA$Anatomic.Site.Simplified)
table(metadata.AURORA$Sample.Type)


## CHAMP process ##
library(ChAMP)

set.seed(46)
testDir <- "..."
myImport <- champ.import(directory = testDir, arraytype = "EPIC")


#2. FILTER RAW DATA
myLoad <- champ.filter(arraytype = "EPIC", intensity = myImport$intensity,
                       Meth = myImport$Meth, UnMeth = myImport$UnMeth, 
                       detP = myImport$detP, filterDetP = T, filterXY = FALSE)


myLoad$pd$Slide <- as.factor(myLoad$pd$Slide)
myLoad$pd$Array <- as.factor(myLoad$pd$Array)

#3. NORMALIZE DATA
myNorm <- champ.norm(beta=myLoad$beta,arraytype="EPIC", plotBMIQ = TRUE)


#4. REMOVE BATCH EFFECT AND COVARIATES. COMBAT.
Covar <- champ.SVD(beta=myNorm %>% as.data.frame(),pd=myLoad$pd,RGEffect = TRUE)

myCombat <- as.data.frame(myNorm)
write.table(myCombat, file = "AURORA Filtered Beta-Values corrected.txt", sep = "\t")
