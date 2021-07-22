RawFrame =  read.csv('GitHub/Data/Ageing and CR genes/ensemblid_conversion.csv')   #GET COEXPRESSION FRAME
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
  ReadingText = paste('C:\\Users\\Usuario\\OneDrive\\Escritorio\\Coexpression\\Coexpression\\', EnsInstances[InstanceIndex], sep = "", collapse = NULL)
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


AgeingGenes = readRDS("G:/gusdany/Documents/AgeingFinal.rds")
CrGenes = readRDS("G:/gusdany/Documents/CrFinal.rds")


Class = ifelse(row.names(CoexpressionFrame) %in% CrGenes,'CR','NotCR')

CoexpressionFrame$Class = Class

write.csv(CoexpressionFrame,'G:\\The Project\\OMA\\CoexpressionDataset.csv')

tail(colnames(CoexpressionFrame))
ReadingText = paste('C:\\Users\\Usuario\\OneDrive\\Escritorio\\Coexpression\\Coexpression\\', EnsInstances[929], sep = "", collapse = NULL)
GeneCoexpression = read.csv(ReadingText)
