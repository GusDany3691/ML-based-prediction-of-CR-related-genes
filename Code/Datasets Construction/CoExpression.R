RawFrame =  read.csv('/ML-based-prediction-of-CR-related-genes/Data/Ageing and CR genes/ensemblid_conversion.csv')   #GET COEXPRESSION FRAME

AgeingGenes = readRDS("/ML-based-prediction-of-CR-related-genes/Data/Ageing and CR genes/AgeingGenes.rds")
CrGenes = readRDS("/ML-based-prediction-of-CR-related-genes/Data/Ageing and CR genes/CrGenes.rds")

Frame = data.frame()
for ( i in 1:dim(RawFrame)[1] ){
  print(i)
  ConversionElements = unlist(strsplit(RawFrame[i,1],';'))
  sFrame = t(data.frame(ConversionElements))
  Frame = rbind(Frame,sFrame)
}
colnames(Frame) = c('Symbol','Ensemble')
row.names(Frame) = Frame$Ensemble
Frame
EnsInstances = Frame$Ensemble


for(InstanceIndex in 1:length(EnsInstances)){
  print(InstanceIndex)
  ReadingText = paste('/ML-based-prediction-of-CR-related-genes/Data/Datasets/External/Coexpression', EnsInstances[InstanceIndex], sep = "", collapse = NULL) # These files are not available in GitHub (there are more than 1000 coexpression profiles files from Ageing-related genes in Genefriends)
  GeneCoexpression = read.csv(ReadingText)
  GeneCoexpression[1,1]
  
  colnames(GeneCoexpression)
  
  Feature = c()
  FeatureValue = c()
  for ( i in 1:dim(GeneCoexpression)[1] ){
    #print(i)
    GeneValuePair = unlist(strsplit(GeneCoexpression[i,1],'\t'))
    Feature[i] = GeneValuePair[1]
    FeatureValue[i] = GeneValuePair[2]
  }
  CoexpressionFrameColumn = data.frame(FeatureValue, row.names = Feature)
  colnames(CoexpressionFrameColumn) = colnames(GeneCoexpression)
  CoexpressionFrameRow = t(CoexpressionFrameColumn)
  row.names(CoexpressionFrameRow)
  
  CoexpressionFrame = rbind(CoexpressionFrame, CoexpressionFrameRow)
}
dim(CoexpressionFrame)
row.names(CoexpressionFrame)
row.names(CoexpressionFrame) = Frame[row.names(CoexpressionFrame),]$Symbol


Class = ifelse(row.names(CoexpressionFrame) %in% CrGenes,'CR','NotCR')
CoexpressionFrame$Class = Class

write.csv(CoexpressionFrame,'/ML-based-prediction-of-CR-related-genes/Data/Datasets/Coexpression.csv')
