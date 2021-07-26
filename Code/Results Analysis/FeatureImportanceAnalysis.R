suppressMessages(library('rattle'))

library('tibble')
library('dplyr') # for data manipulation
library('caret') # for model-building
library('DMwR') # for smote implementation
library('purrr') # for functional programming
library('pROC')
library('data.table') # For dataframe managing
library('DMwR')
library('unbalanced')
library('textshape')
library("KEGGlincs")
library("KEGGgraph")
library("igraph")
library("StarBioTrek")
library("limma")
library("stringr")
library("plyr")
library('org.Hs.eg.db')
library("GenomicRanges")
library("GenomicFiles")
library("gwascat")
library("ensembldb")
library("VariantAnnotation")
library("GenomicAlignments")
library("SummarizedExperiment")
library("AnnotationHub")
library("erma")
library("biomaRt")
library("igraph")


# RETRIEVE DATASETS FOR FEATURES DESCRIPTION
GoDataSet = read.csv("/ML-based-prediction-of-CR-related-genes/Data/Datasets/GoDataset.csv")
PathDipDataset = read.csv("/ML-based-prediction-of-CR-related-genes/Data/Datasets/PathDip.csv")

##### GET GO TERMS DEFINITIONS ###################################################################################################

## Bimap interface: Convert the object to a list
GoTerms <- as.list(GOTERM)
if(length(GoTerms) > 0){
  GOID(GoTerms[[1]])
  Term(GoTerms[[1]])
  Synonym(GoTerms[[1]])
  Secondary(GoTerms[[1]])
  Definition(GoTerms[[1]])
  Ontology(GoTerms[[1]])
}

##### GET FEATURES NAMES FOR GO, KEGG AND DIP (USED FOR WHOLEDATASET) ##############################################################

GoFeatures = colnames(GoDataSet)[2:(length(GoDataSet)-1)]
head(GoFeatures)
tail(GoFeatures)

KEGGpFeatures = colnames(sKeggPathDataset)[2:(length(sKeggPathDataset)-1)]
head(KEGGpFeatures)
tail(KEGGpFeatures)

DipFeatures = colnames(PathDipDataset)[2:(length(PathDipDataset)-1)]
head(DipFeatures)
tail(DipFeatures)

##### SELECT DATA OF INTEREST ######################################################################################################

Class1 = "CR"
Class2 = "NotCR"

# SELECT DATASET FOR QUERYING FEATURE IMPORTANCE AND CR-RELATED INSTANCES
Dataset =  GoDataSet # GoDataSet or PathDipDataset 
head(colnames(Dataset))
Dataset = column_to_rownames(Dataset, "X")
row.names(Dataset)
dim(Dataset)

Dataset$Class = ifelse(Dataset$Class == Class1, Class1, Class2)
sum(Dataset$Class %in% Class1)
sum(Dataset$Class %in% Class2)

# SELECT CR INSTANCES
Instances = row.names(Dataset)
CrInstances = Instances[Dataset$Class %in% Class1]

# RETRIEVE FEATURES IMPORTANCE DATAFRAME
VarImp = read.csv('/ML-based-prediction-of-CR-related-genes/Data/Features_Importance/Original/fi_DIP_CAT.csv')   # fi_DIP_CAT.csv or fi_GO_BRF.csv
VarImp <- VarImp %>% column_to_rownames("Feature")

####################################################################################################################
###### CRRATE STATISTICAL FEATURES IMPORTANCE DATAFRAME ############################################################
####################################################################################################################

method = "DIP"   # DIP or GO 
Description = FALSE
FeatureType = 'Discrete' # Discrete or Continuous

VarImpStatistic = data.frame()
i = 1
for (i in 1:length(VarImp$Feature)){
  print(i)
  #i = 5
  Feature = VarImp$Feature[i]
  Feature = gsub("`", "", Feature)                           # For GO Terms whose Feature is written as GO.XXXXX instead of GO:XXXXX
  Feature = gsub(":", ".", Feature)
  Score =  VarImp$Score[i]
  CrKeyFeatureValues = Dataset[CrInstances,Feature]
  NonCrInstances = setdiff(Instances, CrInstances)
  NonCrKeyFeatureValues = Dataset[NonCrInstances,Feature]
  
  x = CrKeyFeatureValues
  y = NonCrKeyFeatureValues
  
  if (FeatureType == 'Continuous'){
    ttest = t.test(x,y)
    CrMean = mean(x) 
    NonCrMean = mean(y)  
    FoldChange = ((CrMean - NonCrMean) / NonCrMean)
    pvalue = ifelse(CrMean == NonCrMean, 0, ttest$p.value) 
  }
  else{
  res = prop.test(x = c(sum(x), sum(y)), n = c(length(x), length(y)), alternative = "greater")
  pvalue = res$p.value
  }
  Significative = (pvalue <= 0.05)
  
  
  # FANCY NUMBERS
  Score = round(Score,3)
  if (CrMean != 0) CrMean = as.character(if (CrMean != 0) ifelse(abs(CrMean) < 0.001 || abs(CrMean) > 99999, formatC(CrMean, format = "e", digits = 2),  round(CrMean,2) )) else CrMean = "0"
  if (NonCrMean != 0) NonCrMean = as.character(if (NonCrMean != 0) ifelse(abs(NonCrMean) < 0.001 || abs(NonCrMean) > 99999, formatC(NonCrMean, format = "e", digits = 2),  round(NonCrMean,2) )) else NonCrMean = "0"
  if (pvalue != 0 ) pvalue = as.character(if (pvalue != 0) ifelse(pvalue < 0.001, formatC(pvalue, format = "e", digits = 2), round(pvalue,3))) else pvalue = "0"
  NonCrMean
  
  # ---  RETRIEVE DESCRIPTION -------------------------------------------------------------------------------
  
  # GO TERMS
  if (method == "GO"){
    Description = TRUE
    Feature = gsub("[.]", ":", Feature)
    Definition = Term(GoTerms[[Feature]])
    SelectedCrAmount = as.character(sum(x))
    CrAmount = as.character(length(x))
    SelectedNonCrAmount = as.character(sum(y))
    NonCrAmount = as.character(length(y))
    CrMean = paste(CrMean, "%", sep = "")
    NonCrMean = paste(NonCrMean, "%", sep = "")
  } 
  
  # PATHDIP
  if (method == "Dip"){
    Description = TRUE
    Definition = strsplit(Feature, "[.][.][.][.]")[[1]][2]
    Feature = strsplit(Feature, "[.][.][.][.]")[[1]][1]
  }
  
  
  if (method == "None"){
    Description = TRUE
    Definition = "-"
  }
  
  # ---------------------------------------------------------------------------------------------------------
  
  if (Description)
    cVarImpStatistic = data.frame(Feature, Definition, Score, 'CR Mean' = CrMean, "NonCR mean" = NonCrMean, pvalue)
  else
    cVarImpStatistic = data.frame(Feature, Score, 'CR %' = CrMean, "NonCR %" = NonCrMean, pvalue)
  VarImpStatistic = rbind(VarImpStatistic, cVarImpStatistic)
  
}

VarImpStatistic = VarImpStatistic[order(VarImpStatistic$Score, decreasing = TRUE),]

# MULTIPLE TESTS CORRECTION
VarImpStatistic$padjust = p.adjust(VarImpStatistic$pvalue, method = 'bonferroni')

# SAVING
write.csv('/ML-based-prediction-of-CR-related-genes/Data/Features_Importance/Scaled_and_Discussed/FeaturesImportance(PathDIP_CAT).csv', pathfile)
