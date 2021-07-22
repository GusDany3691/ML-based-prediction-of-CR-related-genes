#LOADING AGEING AND CR GENES
AgeGenes = read.csv("E:/gusdany/Documents/AgeingOMA.csv")$AgeGenes
CrGenes = read.csv("E:/gusdany/Documents/CrOMA.csv")$CrGenes

# LOAD GTEXT DATA
Gtex = read.csv("E:/gusdany/Desktop/El Proyecto/Databases/DataSetsHagr/Gtex.csv")

AgeingGtex = Gtex[Gtex$Gene %in% AgeGenes,]
AgeingGtex$Class = ifelse(AgeingGtex$Gene %in% CrGenes, "CR", "NotCR")

# SAVE
write.csv(AgeingGtex,"D:/gusdany/Desktop/El Proyecto/Databases/DataSetsHagr/Final/OMA/GtexDataset.csv")
