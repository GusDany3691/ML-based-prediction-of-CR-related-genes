library(scales)

DIP = read.csv('/ML-based-prediction-of-CR-related-genes/Data/CR-probabilities/Strongest models/CAT_DIP.csv')
GO = read.csv('/ML-based-prediction-of-CR-related-genes/Data/CR-probabilities/Strongest models/BRF_GO.csv')

GO = GO[order(GO$Prob, decreasing = TRUE),]
DIP = DIP[order(DIP$Prob, decreasing = TRUE),]

GO$Ranking = seq(1,length(GO$Label))
DIP$Ranking = seq(1,length(DIP$Label))

CommonGenes = intersect(GO$Label,DIP$Label)

commonGO = GO[GO$Label %in% CommonGenes,]
commonDIP = DIP[DIP$Label %in% CommonGenes,]

row.names(commonGO) = commonGO$Label
row.names(commonDIP) = commonDIP$Label

commonDIP = commonDIP[CommonGenes,]
commonGO = commonGO[CommonGenes,]

Label = CommonGenes
Test = commonDIP$Test
DIP_Ranking = commonDIP$Ranking
GO_Ranking = commonGO$Ranking
Average_Ranking = (DIP_Ranking + GO_Ranking) / 2
DIP_Prob_Original = commonDIP$Prob
GO_Prob_Original = commonGO$Prob
Average_Prob_Original = ( DIP_Prob_Original + GO_Prob_Original ) / 2
DIP_Prob_Normalised = rescale(DIP_Prob_Original)
GO_Prob_Normalised = rescale(GO_Prob_Original)
Average_Prob_Normalised = ( DIP_Prob_Normalised + GO_Prob_Normalised ) / 2

JointAnalysis = data.frame(Label, Test, DIP_Ranking, GO_Ranking, Average_Ranking, DIP_Prob_Original, GO_Prob_Original, Average_Prob_Original, DIP_Prob_Normalised, GO_Prob_Normalised, Average_Prob_Normalised)

write.csv(JointAnalysis, '/ML-based-prediction-of-CR-related-genes/Data/CR-probabilities/Strongest models/JointAnalysis.csv')
