#LOADING AGEING AND CR GENES
AgeGenes = read.csv("/ML-based-prediction-of-CR-related-genes/Data/Ageing and CR genes/AgeingGenes.csv")$AgeGenes
CrGenes = read.csv("/ML-based-prediction-of-CR-related-genes/Data/Ageing and CR genes/CrGenes.csv.csv")$CrGenes

# LOAD GTEXT DATA
Gtex = read.csv("/ML-based-prediction-of-CR-related-genes/Data/Datasets/External/Gtex.csv")

AgeingGtex = Gtex[Gtex$Gene %in% AgeGenes,]
AgeingGtex$Class = ifelse(AgeingGtex$Gene %in% CrGenes, "CR", "NotCR")

# SAVE
write.csv(AgeingGtex,"/ML-based-prediction-of-CR-related-genes/Data/Datasets/Gtex.csv")
