##########################################################################################################################
###### LIBRARIES #########################################################################################################
##########################################################################################################################

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

##########################################################################################################################
###### GET GENES #########################################################################################################
##########################################################################################################################
hCR == CrGenes

AgeingGenes = read.csv("/GitHub/Data/Ageing and CR genes/AgeingFinal.csv")
CrGenes = read.csv("/GitHub/Data/Ageing and CR genes/CrFinal.csv")


Class = ifelse((AgeingGenes %in% CrGenes),"CR","NotCR")

ClassFrame = data.frame('Class' = Class)
row.names(ClassFrame) = AgeingGenes 

Class = Data.frame()
#########################################################################################################################
###### GET GO TERMS FROM GENES ##########################################################################################
#########################################################################################################################

AllInstances = AgeingGenes
sAllGoPerInstanceList = list()
sBpGoPerIntanceList = list()
BpGoPerIntanceList = list()
sBpGoTerms = c()
BpGoTerms = c()
AllGoAncestorsPerInstance = c()

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

###### QUERING GO TERMS AND ANCESTORS ######################################################################################
for (i in 1:length(AllInstances)){
  print(i)
  # ----- RETRIEVING GOs ----------------------------------------------------------------------------------------------------
  Instance = AllInstances[i]
  sAllGoPerInstance = getBM(attributes = c('hgnc_symbol','go_id', 'namespace_1003','go_linkage_type'),
                         filters = c('hgnc_symbol'),
                         values = Instance,
                         mart = human)
  sAllGoPerInstanceList[[Instance]] = sAllGoPerInstance
  
  # ----- FILTERING GOs -----------------------------------------------------------------------------------------------------
  BpIndices = sAllGoPerInstance$namespace_1003=="biological_process" 
  FilteredBpIndices = sAllGoPerInstance$namespace_1003=="biological_process" & sAllGoPerInstance$go_linkage_type=="IEA"
  InstanceBpTerms = sAllGoPerInstance[BpIndices,]$go_id
  
  if(length(InstanceBpTerms) > 0){
    if(!is.na(InstanceBpTerms)){
      
      # ----- SAVING GOs ----------------------------------------------------------------------------------------------------
      sBpGoPerIntanceList[[Instance]] = InstanceBpTerms
      sBpGoTerms = unique(c(sBpGoTerms, InstanceBpTerms))
      
      # ----- ANCESTORS -----------------------------------------------------------------------------------------------------
      InstanceGoAncestors = c()
      CurrentGoes = InstanceBpTerms
      for (j in 1:length(CurrentGoes)){
        CurrentGo = CurrentGoes[j]
        CurrentAncestors = GOBPANCESTOR[[CurrentGo]]
        if (!is.null(CurrentAncestors))
        InstanceGoAncestors = unique(c(InstanceGoAncestors, CurrentAncestors, CurrentGo))   # Get Current Go Ancestors and Current Go
      }
      InstanceGoAncestors = setdiff(InstanceGoAncestors, "all")
      BpGoPerIntanceList[[Instance]] = InstanceGoAncestors                                        # List of all go terms with ancestros per instance
      BpGoTerms = unique(c(BpGoTerms, InstanceGoAncestors))                                       # Character of all go terms with ancestros
    }
  }
}

CompleteFeaturesPerInstanceList = BpGoPerIntanceList
IncompleteFeaturesPerInstanceList = sBpGoPerIntanceList
  
###### DATASETS #########################################################################################################

Instances = names(sBpGoPerIntanceList)
GoDataSet = data.frame()
sGoDataSet = data.frame()

for (i in 1:length(Instances)){
  print(i)
  
  # GO TERMS WITHOUT ANCESTORS FRAME
  Instance = Instances[i]
  ColNames = c(sBpGoTerms, "Class")
  RowName = Instance
  FeatureValues = sBpGoTerms %in% sBpGoPerIntanceList[[Instance]]
  Class = Instance %in% CrGenes
  Row = t(c(FeatureValues, Class))
  sGoFrame = data.frame(Row)                                                                          # InCompleteFrame
  colnames(sGoFrame) = ColNames
  row.names(sGoFrame) = RowName
  
  # GO TERMS WITH ANCESTORS FRAME
  ColNames = c(BpGoTerms, "Class")
  RowName = Instance
  FeatureValues = BpGoTerms %in% BpGoPerIntanceList[[Instance]]
  Class = Instance %in% CrGenes
  Row = t(c(FeatureValues, Class))
  GoFrame = data.frame(Row)                                                                          # CompleteFrame
  colnames(GoFrame) = ColNames
  row.names(GoFrame) = RowName
  
  GoDataSet = rbind(GoDataSet, GoFrame)                                                              # CompleteDataSet                                             
  sGoDataSet = rbind(sGoDataSet, sGoFrame)                                                           # InCompleteDataSet   
}

setdiff(BpGoTerms,unique(ParentsChildFrame$Go))


save.csv(GoDataSet,"/GitHub/Data/Datasets/GoDataSet.rds")
save.csv(sGoDataSet,"/GitHub/Data/Datasets/sGoDataSet.rds")
