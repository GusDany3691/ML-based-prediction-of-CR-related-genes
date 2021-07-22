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


AgeingGenes = readRDS("/GitHub/Data/Ageing and CR genes/AgeingFinal.rds")
CrGenes = readRDS("/GitHub/Data/Ageing and CR genes/CrFinal.rds")


# - FOR EXPRESSION ANALYSIS ----------------------------------------------------------------------------------
# MANUAL
AgeingGenes = c('MAP3K6', 'ADH1A', 'ADH1B', 'ADH1C', 'SC5D')
ClassFrame = c("CR","CR","CR","CR","CR")

# NEW CR-RELATED GENES
UpCrHuman = read.csv("C:/Users/gusdany/Desktop/El Proyecto/Databases/Datasets/UpCR_HM.csv")
DownCrHuman = read.csv("C:/Users/gusdany/Desktop/El Proyecto/Databases/Datasets/DownCR_HM.csv")

# RETRIEVE CR-RELATED HUMAN GENES
AgeingGenes  = c(UpCrHuman$Human, DownCrHuman$Human)



# ------------------------------------------------------------------------------------------------------------

ClassFrame = data.frame(ifelse(AgeingGenes %in% CrGenes,"CR","NotCR"),row.names = AgeingGenes)

AgeCrGenes = readRDS("AgeCrGenes.rds")
KnownAgeinGenes = readRDS("KnownAgeinGenes.rds")
UnKnownAgeingGenes = readRDS("UnKnownAgeingGenes.rds")

write.csv(as.data.frame(AgeingGenes),"D:/gusdany/Desktop/El Proyecto/Databases/DataSetsHagr/Final/OMA/AgeingFinal.csv")
write.csv(as.data.frame(CrGenes),"D:/gusdany/Desktop/El Proyecto/Databases/DataSetsHagr/Final/OMA/CrFinal.csv")

write.csv(ClassFrame,"D:/gusdany/Desktop/El Proyecto/Databases/DataSetsHagr/Final/OMA/ClassFrame.csv")

Class = ifelse((AgeingGenes %in% CrGenes),"CR","NotCR")

ClassFrame = data.frame('Class' = Class)
row.names(ClassFrame) = AgeingGenes 
write.csv(ClassFrame,"C:/Users/gusdany/Desktop/El Proyecto/Databases/DataSetsHagr/Final/OMA/ClassFrame.csv")


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

GoDataset2 = GoDataSet
for (i in 1:(length(GoDataset2)-1)){
  GoDataset2[[i]] = ifelse(GoDataSet[[i]] == TRUE, 1, 0)
  print(i)
}
GoDataset2[['Class']] = ifelse(GoDataSet2[['Class']] == 1, TRUE, FALSE)
sGoDataset2 = sGoDataSet

for (i in 1:(length(sGoDataset2)-1)){
  sGoDataset2[[i]] = ifelse(sGoDataSet[[i]] == TRUE, 1, 0)
  print(i)
}
sGoDataset2[['Class']] = ifelse(sGoDataset2[['Class']] == 1, TRUE, FALSE)

dim(GoDataset2)
dim(sGoDataset2)

write.csv(GoDataset2,"C:/Users/gusdany/Desktop/El Proyecto/Databases/DataSetsHagr/Final/OMA/noGoDataset.csv")
write.csv(sGoDataset2,"C:/Users/gusdany/Desktop/El Proyecto/Databases/DataSetsHagr/Final/OMA/nosGoDataset.csv")

#gDB = read.csv("D:/gusdany/Desktop/El Proyecto/Databases/DataSetsHagr/Final/OMA/GoDataset.csv")
#GoCR = gDB[gDB["Class"] == 'CR',]$X


saveRDS(GoDataSet,"C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/GO/Final/OMA/GoDataSet.rds")
saveRDS(sGoDataSet,"C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/GO/Final/OMA/sGoDataSet.rds")
##### GRAPH ############################################################################################################

GoDataSet = readRDS("C:/Users/gusdany/Desktop/Classification/RDS/GO/GoDataSet.rds")
Instances = row.names(GoDataSet)
Features = colnames(GoDataSet)
BpGoTerms = Features
# GENERATE A PARENT-CHILD DATAFRAME FOR ALL ANCESTOR SET
ParentsChildFrame = data.frame()
for (i in 1:length(BpGoTerms)){  # All Go term with ancestors
  print(i)
  CurrentGo = BpGoTerms[i]
  CurrentParents = GOBPPARENTS[[CurrentGo]]
  RepeatedCurrentGo = rep(CurrentGo, length(CurrentParents))
  cParentsChildFrame = data.frame( 'Parent' = CurrentParents, 'Go' = RepeatedCurrentGo )
  ParentsChildFrame = rbind(ParentsChildFrame, cParentsChildFrame)
}
ParentsChildFrame = ParentsChildFrame[ParentsChildFrame$Parent != "all",]

# GENERATE A PARENT-CHILD DATAFRAME FOR ALL NON-ANCESTOR SET
sParentsChildFrame = data.frame()
for (i in 1:length(sBpGoTerms)){  # All Go term without ancestors
  print(i)
  CurrentGo = BpGoTerms[i]
  CurrentParents = GOBPPARENTS[[CurrentGo]]
  RepeatedCurrentGo = rep(CurrentGo, length(CurrentParents))
  cParentsChildFrame = data.frame( 'Parent' = CurrentParents, 'Go' = RepeatedCurrentGo )
  sParentsChildFrame = rbind(sParentsChildFrame, cParentsChildFrame)
}
sParentsChildFrame = sParentsChildFrame[sParentsChildFrame$Parent != "all",]

# GENERATE GRPHS
GoGraph = graph_from_data_frame(ParentsChildFrame, directed = TRUE, vertices = NULL)
sGoGraph = graph_from_data_frame(sParentsChildFrame, directed = TRUE, vertices = NULL)

GoGraph = simplify(GoGraph, remove.multiple = TRUE, remove.loops = TRUE)
sGoGraph = simplify(sGoGraph, remove.multiple = TRUE, remove.loops = TRUE)

CompleteGraph = GoGraph
IncompleteGraph = sGoGraph

saveRDS(CompleteGraph,"C:/Users/gusdany/Desktop/Classification/RDS/GO/GoCompleteGraph.rds")
GoDataSet = readRDS("C:/Users/gusdany/Desktop/Classification/RDS/GO/GoDataSet.rds")

#########################################################################################################################
##### GRAPH-BASED FEATURES SELECTION PREPROCESSING ######################################################################
#########################################################################################################################

###### RELEVANCE CALCULATION ############################################################################################

CompleteFeatures = BpGoTerms
InompleteFeatures = sBpGoTerms
CompleteFeaturesFrame = GoDataSet
InompleteFeaturesFrame = sGoDataSet

Rlist = list()
NumInstances = length(Instances)
NumFeatures = length(CompleteFeatures)
iNumFeatures = length(IncompleteFeatures)
for (i in 1:NumFeatures){
  print(i)
  CyGy = 0
  CyGn = 0
  CnGy = 0
  CnGn = 0
  for (j in 1:NumInstances){
    Feature = CompleteFeatures[[i]]
    if (CompleteFeaturesFrame[["Class"]][j] == TRUE){
      if(CompleteFeaturesFrame[[Feature]][j] == TRUE)
        CyGy = CyGy + 1
      else
        CyGn = CyGn + 1
    } else {
      if(CompleteFeaturesFrame[[Feature]][j])
        CnGy = CnGy + 1
      else
        CnGn = CnGn + 1
    }
  }
  Pyy = CyGy / NumInstances
  Pyn = CyGn / NumInstances
  Pny = CnGy / NumInstances
  Pnn = CnGn / NumInstances
  Rlist[[Feature]] = (Pyy - Pyn)^2 + (Pyn - Pnn)^2
}
Rchar = unlist(Rlist)
R = data.frame('Feature' = names(Rchar), 'Score' = unname(Rchar))

write.csv(CompleteFeaturesFrame,"CompleteFeaturesFrame.csv")

##### CREATE LISTS OF ANCESTORS, DESCENDANTS AND STATUS #################################################################

# CREATEM EMPTY LISTS
Anc = list()
Dec = list()
Status = list()
iAnc = list()
iDec = list()
iStatus = list()

# FOR COMPLETE GRAPH: FILL ANCESTORS AND DESCENDANTS LISTS. INITIALIZE STATUS LIST 
for (i in 1 :length(CompleteFeatures)){
  key <- CompleteFeatures[i]                                                        
  Anc[[key]] = subcomponent(CompleteGraph, CompleteFeatures[i], mode = "in")[-1]  
  Dec[[key]] = subcomponent(CompleteGraph, CompleteFeatures[i], mode = "out")[-1]
  Status[[key]] = "Selected"
  print(i)
}

# FOR INCOMPLETE GRAPH: FILL ANCESTORS AND DESCENDANTS LISTS. INITIALIZE STATUS LIST 
for (i in 1 :length(IncompleteFeatures)){
  key <- IncompleteFeatures[i]                                                        
  iAnc[[key]] = subcomponent(IncompleteGraph, IncompleteFeatures[i], mode = "in")[-1]  
  iDec[[key]] = subcomponent(IncompleteGraph, IncompleteFeatures[i], mode = "out")[-1]
  iStatus[[key]] = "Selected"
  print(i)
}

#########################################################################################################################
##### SAVE AND LOAD VARIABLES ###########################################################################################
#########################################################################################################################
saveRDS(Instances, file = "Instances.rds")
saveRDS(CompleteGraph, file = "CompleteGraph.rds")
saveRDS(IncompleteGraph, file = "IncompleteGraph.rds")
saveRDS(CompleteFeatures, file = "CompleteFeatures.rds")
saveRDS(IncompleteFeatures, file = "IncompleteFeatures.rds")
saveRDS(CompleteFeaturesFrame, file = "CompleteFeaturesFrame.rds")
saveRDS(IncompleteFeaturesFrame, file = "IncompleteFeaturesFrame.rds")
saveRDS(CompleteFeaturesPerInstanceList, file = "CompleteFeaturesPerInstanceList.rds")
saveRDS(IncompleteFeaturesPerInstanceList, file = "IncompleteFeaturesPerInstanceList.rds")

Instances = readRDS(file = "Instances.rds")
CompleteGraph = readRDS(file = "CompleteGraph.rds")
IncompleteGraph = readRDS(file = "IncompleteGraph.rds")
CompleteFeatures = readRDS(file = "CompleteFeatures.rds")
IncompleteFeatures = readRDS(file = "IncompleteFeatures.rds")
CompleteFeaturesFrame = readRDS(file = "CompleteFeaturesFrame.rds")
IncompleteFeaturesFrame = readRDS(file = "IncompleteFeaturesFrame.rds")
CompleteFeaturesPerInstanceList = readRDS(file = "CompleteFeaturesPerInstanceList.rds")
IncompleteFeaturesPerInstanceList = readRDS(file = "IncompleteFeaturesPerInstanceList.rds")

#########################################################################################################################
##### GRAPH-BASED FEATURES SELECTION ALGORITHMS #########################################################################
#########################################################################################################################

NumInstances = length(Instances)
NumFeatures = length(CompleteFeatures)

##### # CREATE FEATURES FREQUENCY LIST ##################################################################################

FrequencyList1 = list()
FrequencyList2 = list()
FrequencyList3 = list()
for (i in 1:length(CompleteFeatures)){
  print(i)
  FrequencyList1[[CompleteFeatures[i]]] = 0
  FrequencyList2[[CompleteFeatures[i]]] = 0
  FrequencyList3[[CompleteFeatures[i]]] = 0
}

##### ALGORITHM 1 #######################################################################################################

# INITIALIZE TRAN AND TEST DATASETS *

# INITIALIZE LISTS
FrequencyFrame1 = as.data.frame(unlist(FrequencyList1))
colnames(FrequencyFrame1) = "Frequency"
FrequencyFrame1 = tibble::rownames_to_column(FrequencyFrame1, "Feature")
row.names(FrequencyFrame1) = FrequencyFrame1$Feature
i = 1

SelectedFeaturesPerIntance1 = list()
# BEGIN FEATURES FILTERING
for (i in 1:NumInstances){
  print(i)
  print(dim(FrequencyFrame1))
  Instance = Instances[i]
  for (j in 1:NumFeatures){
    Feature = CompleteFeatures[j]
    if (Feature %in% CompleteFeaturesPerInstanceList[[Instance]]){
      FeatureAncestors = names(Anc[[Feature]])
      NumFeatureAnc = length(FeatureAncestors)
      if (NumFeatureAnc > 0)
        Status[FeatureAncestors] = "Removed"
    } 
    else {
      FeatureDescendants = Dec[[Feature]]
      NumFeatureDec = length(FeatureDescendants)
      if (NumFeatureDec > 0)
        Status[FeatureDescendants] = "Removed"
    }
  }
  
  # FILTER ONLY SELECTED FEATURES FOR THE INSTANCE
  SelectedFeatures = names(Status[Status == "Selected"])
  FrequencyFrame1[SelectedFeatures,]$Frequency = FrequencyFrame1[SelectedFeatures,]$Frequency + 1
  
  # ASSIGN SELECTED ITEMS TO INSTANCE
  SelectedFeaturesPerIntance1[[Instance]] = SelectedFeatures
  
  # SAVE ACCUMULATED VARIABLES 
  saveRDS(FrequencyFrame1, file = "C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/GO/GoFrequencyFrame1.rds")
  saveRDS(SelectedFeaturesPerIntance1, file = "C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/GO/GoSelectedFeaturesPerIntance1.rds")
  
  # CLASSIFICATION
  
  # REASSIGN SELECTED
  Status[] = "Selected"
  
}

###### ALGORITHM 2 #####################################################################################################

# CREATE LISTS OF ANCESTORS, DESCENDANTS AND STATUS
Anc2 = list()
Dec2 = list()

# FOR COMPLETE GRAPH: FILL ANCESTORS AND DESCENDANTS LISTS. INITIALIZE STATUS LIST
for (i in 1 :length(CompleteFeatures)){
  Feature <- CompleteFeatures[i]                                                        
  Anc2[[key]] = all_simple_paths(CompleteGraph, from = Root, to = Feature)  
  Dec2[[key]] = all_simple_paths(CompleteGraph, from = Feature, to = Leaves)
  print(i)
}

SelectedFeaturesPerIntance2 = list()

for (i in 1:NumInstances){
  Instance = Instances[i]
  for (j in 1:NumFeatures){
    Feature = CompleteFeatures[j]
    print(j)
    if (Feature %in% CompleteFeaturesPerInstanceList[[Instance]]){
      # Get paths from root to Current Feature
      PathsToRoot = all_simple_paths(CompleteGraph, from = Root, to = Feature) 
      # Keep the Most Relevant Feature for every path from Current Feature to Root
      if (length(PathsToRoot) > 0 ) {
        for (k in 1:length(PathsToRoot)){
          FeaturesInPath = names(PathsToRoot[[k]])
          FeaturesInPathScoresFrame = R[R$Feature %in% FeaturesInPath,]
          MRFindex = FeaturesInPathScoresFrame$Score %in%  max(FeaturesInPathScoresFrame$Score)
          MRF = as.character(FeaturesInPathScoresFrame[MRFindex,]$Feature)
          ToRemoveFeatures = setdiff(FeaturesInPath, MRF)
          Status[ToRemoveFeatures] = "Removed"
        }
      }
    } 
    else {
      
      # Get paths from Current Feature to leaves
      PathsToLeaves = all_simple_paths(CompleteGraph, from = Feature, to = Leaves) 
      # Keep the Most Relevant Feature for every path from Current Feature to Root
      if (length(PathsToLeaves) > 0 ) {
        for (k in 1:length(PathsToLeaves)){
          FeaturesInPath = names(PathsToLeaves[[k]])
          FeaturesInPathScoresFrame = R[R$Feature %in% FeaturesInPath,]
          MRFindex = FeaturesInPathScoresFrame$Score %in%  max(FeaturesInPathScoresFrame$Score)
          MRF = as.character(FeaturesInPathScoresFrame[MRFindex,]$Feature)
          ToRemoveFeatures = setdiff(FeaturesInPath, MRF)
          Status[ToRemoveFeatures] = "Removed"
        }
      }
    }
  }
  
  # FILTER ONLY SELECTED FEATURES FOR THE INSTANCE
  SelectedFeatures = names(Status[Status == "Selected"])
  FrequencyFrame2[SelectedFeatures,]$Frequency = FrequencyFrame2[SelectedFeatures,]$Frequency + 1
  
  # ASSIGN SELECTED ITEMS TO INSTANCE
  SelectedFeaturesPerIntance2[[Instance]] = SelectedFeatures
  
  # SAVE ACCUMULATED VARIABLES 
  saveRDS(FrequencyFrame2, file = "FrequencyFrame2.rds")
  saveRDS(SelectedFeaturesPerIntance2, file = "SelectedFeaturesPerIntance2.rds")
  
  # CLASSIFICATION
  
  # REASSIGN SELECTED
  Status[] = "Selected"
  #print("ESTO ES")
  print(i)
  
}


###### ALGORITHM 3 ######################################################################################################

FrequencyFrame3 = as.data.frame(unlist(FrequencyList3))
colnames(FrequencyFrame3) = "Frequency"
FrequencyFrame3 = tibble::rownames_to_column(FrequencyFrame3, "Feature")
row.names(FrequencyFrame3) = FrequencyFrame3$Feature

SelectedFeaturesPerIntance3 = list()
# BEGIN FEATURES FILTERING
for (i in 1:NumInstances){
  Instance = Instances[i]
  print(i)
  print(dim(FrequencyFrame3))
  for (j in 1:NumFeatures){
    Feature = CompleteFeatures[j]
    Rfeature = R[R$Feature %in% Feature,]
    if (Feature %in% CompleteFeaturesPerInstanceList[[Instance]]){
      FeatureAncestors = names(Anc[[Feature]])
      NumFeatureAnc = length(FeatureAncestors)
      if (NumFeatureAnc > 0) {
        for (k in 1:NumFeatureAnc){
          Ancestor = FeatureAncestors[k]
          Rancestor = R[R$Feature %in% Ancestor,]
          if (Rancestor$Score <= Rfeature$Score)
            Status[[Ancestor]] = "Removed"
        }
      }
    } 
    else {
      FeatureDescendants = names(Dec[[Feature]])
      NumFeatureDec = length(FeatureDescendants)
      if (NumFeatureDec > 0){
        for (k in 1:NumFeatureDec){
          Descendant = FeatureDescendants[k]
          Rdescendant = R[R$Feature %in% Descendant,]
          if (Rdescendant$Score <= Rfeature$Score)
            Status[[Descendant]] = "Removed"
        }
      }
    }
  }
  
  # FILTER ONLY SELECTED FEATURES FOR THE INSTANCE
  SelectedFeatures = names(Status[Status == "Selected"])
  FrequencyFrame3[SelectedFeatures,]$Frequency = FrequencyFrame3[SelectedFeatures,]$Frequency + 1
  
  # ASSIGN SELECTED ITEMS TO INSTANCE
  SelectedFeaturesPerIntance3[[Instance]] = SelectedFeatures
  
  # SAVE ACCUMULATED VARIABLES 
  saveRDS(FrequencyFrame3, file = "C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/GO/GoFrequencyFrame3.rds")
  saveRDS(SelectedFeaturesPerIntance3, file = "C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/GO/GoSelectedFeaturesPerIntance3.rds")
  
  # CLASSIFICATION
  
  # REASSIGN SELECTED
  Status[] = "Selected"
  
  
}

#########################################################################################################################
##### GENERATE GRAPH-BASED FILTERED FEATURE DATAFRAMES ##################################################################
#########################################################################################################################

# FRAME FROM ALGORITHM 1
FrequencyFrame1 = as.data.frame(unlist(FrequencyList1))
colnames(FrequencyFrame1) = "Frequency"
FrequencyFrame1 = tibble::rownames_to_column(FrequencyFrame1, "Feature")
row.names(FrequencyFrame1) = FrequencyFrame1$Feature

# FRAME FROM ALGORITHM 2
FrequencyFrame2 = as.data.frame(unlist(FrequencyList2))
colnames(FrequencyFrame2) = "Frequency"
FrequencyFrame2 = tibble::rownames_to_column(FrequencyFrame2, "Feature")
row.names(FrequencyFrame2) = FrequencyFrame2$Feature

# FRAME FROM ALGORITHM 3
FrequencyFrame3 = as.data.frame(unlist(FrequencyList3))
colnames(FrequencyFrame3) = "Frequency"
FrequencyFrame3 = tibble::rownames_to_column(FrequencyFrame3, "Feature")
row.names(FrequencyFrame3) = FrequencyFrame3$Feature


GoDataSet[1,]
