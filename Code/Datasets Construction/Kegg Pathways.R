# LOADING PACKAGES AND LIBRARIES
BiocManager::install("KEGGlincs")
BiocManager::install("KEGGgraph")
BiocManager::install("StarBioTrek")
BiocManager::install("limma")
BiocManager::install("stringr")
BiocManager::install("DO.db")

library("KEGGlincs")
library("KEGGgraph")
library("igraph")
library("StarBioTrek")
library("limma")
library("stringr")
library("plyr")
library('org.Hs.eg.db')

# DOWNLOAD LIST OF GENES-RELATED KEGG PATHWAYS                                                                                                                                               
DfAllKegg <- getGeneKEGGLinks(species="hsa")                                                                                 
DfAllKegg$Symbol <-  as.character(mapIds(org.Hs.eg.db, DfAllKegg$GeneID, 'SYMBOL', 'ENTREZID'))#mapIds(org.Hs.eg.db, tab$GeneID, column="SYMBOL", keytype="ENTREZID")
DfAllKegg$PathwayID = str_remove(DfAllKegg$PathwayID, "path:")
DfAllKegg$GeneID =  paste("hsa:", DfAllKegg$GeneID ,sep="")
head(DfAllKegg)



# GET AGEING AND CR GENES
AgeGenes = readRDS("/ML-based-prediction-of-CR-related-genes/Data/Ageing and CR genes/AgeingGenes.csv")
CrGenes = readRDS("/ML-based-prediction-of-CR-related-genes/Data/Ageing and CR genes/CrGenes.csv")

# GET AGE AND CR GENES PRESENT IN KEGG DATABASE
AllKeggGenes = DfAllKegg$Symbol
AgeKeggGenes = intersect(AllKeggGenes, AgeGenes)

InstancesLabel = AgeKeggGenes 
Instances = unique(DfAllKegg[DfAllKegg$Symbol %in% AgeKeggGenes,]$GeneID)

# #########################################################################################################
# ################################### GENERATE DATAFRAME LISTS CODE ####################################### 
# #########################################################################################################

InstancesList = list()
ListFrames = list()
GraphInfoList = list()
InstanceGraphInfoList = list()
#AcceptedFrameIndices = c()

# CREATE LIST
for ( i in 1:length(Instances)){   
  #print(i)   
  AcceptedFrameIndices = c()
  ListFrames = list()
  GraphInfoList = list()
  
  Instance = Instances[i] 
  InstanceLabel = InstancesLabel[i]                                                                  # Define Instance
  GeneKeggFrame = DfAllKegg[DfAllKegg$GeneID == Instance,]                                           # Retrieve All Pathways and Genes related to this Gene
  GenePaths = GeneKeggFrame$PathwayID                                                                # Select only pathways of the preceding query
  
  # CLEAR PATHWAYS 
  for ( j in 1:length(GenePaths)){                                                                   # Clear Past Pathway Dataframes for further analysis
    if (exists(GenePaths[j])) {                                                                      # Check if a nprevious pathway dataframe with this name exists. If so, erase it
      rm(list = c(GenePaths[j]))
    }
  }
  
  # CREATE DATAFRAMES CONTAINING GENES FOR EACH PATH
  for ( j in 1:length(GenePaths)){                                                                  # For each path
    Pathway = GenePaths[j]                                                                          # Geth the jth path name
    graphAndMore = GenerateGraph(Instance, Pathway)                                                 # HERE I GET GRAPH AND EDGE INFO!!!
    print(paste(i,":" ,j))
    if (!is.character(graphAndMore)){                                                               # If I could get a query for the pathway
      FeaturesAndGraphInfo = InfluenceByDeletion(graphAndMore, Instance)                            # Get the features' scores for the current Instance and Pathway
      GeneScoreFrame = FeaturesAndGraphInfo[[1]]
      EdgeRels = FeaturesAndGraphInfo[[2]]
      InstanceNames = FeaturesAndGraphInfo[[3]]
      graphAndMore[[2]] = cbind(EdgeRels, graphAndMore[[2]])
      GraphInfoList[[Pathway]] = graphAndMore
      
      if (length(GeneScoreFrame) != 0) {                                                             # If the Dataframe contains features
        AcceptedFrameIndices = c(AcceptedFrameIndices, j)
        gene = GeneScoreFrame[,1]                                                                    # Gene features names
        score <- GeneScoreFrame[,2]                                                                  # Gene Features values 
        path = rep(Pathway, times = length(gene))                                                   # Introduce the pathname as one of the future dataframe columns
        assign(GenePaths[j], data.frame(path, gene, score))                                          # Create a DataFrame indireclty called 'hsa30000': GenePaths[1] = "hsa30000" containign 'gene', 'score' and 'path' for the current pathway 
        #print(j)
      }
    } 
  }

  # GENERATE LIST OF DATAFRAMES FOR THE QUERY GENE
  GeneFrame = data.frame()                                                                           # Data fraim containing all the path frames fora single gene
  ListFrames = list()                                                                                # List of feature Frames corresponding to a specific gene
  MixedFinalFrame = data.frame()
  
  # FRAMES UNIFICATION
  for (j in 1:length(GenePaths)){                                                                    # For each of the pathways correponding to instance gene
    if (j %in% AcceptedFrameIndices ){
      #print(j)
      ListFrames[[GenePaths[j]]] = get(GenePaths[j])                                                   # Set Path key to list.. Eg. get return the variable hsaXX000 = Dataframe.. from the first path: GenePaths[i] = "hsaXX000"
      MixedFinalFrame =  rbind(MixedFinalFrame, get(GenePaths[j]) )                                    # Bind all the gene-paths's dataframes
    }
  }
  
  # TAKE DATAFRAME LIST AS AN ELEMENT FOR THE INSTANCES LISTS
  InstancesList[[Instance]] = ListFrames                                                             # Attach the previous list of dataframes as the gene-instance element of ListGenes
  repeatedGenes = unique(as.character(MixedFinalFrame$gene[duplicated(MixedFinalFrame$gene)]))       # Look for repeated gene features in the gene instance
  NonRepeatedGenes = setdiff(as.character(unique(MixedFinalFrame$gene)), repeatedGenes)                             # Extract non repeated gene features in the gene instance
  FinalFrame = data.frame()                                                                          # Create FinalDataFrame variable
  
  # COMBINE THE FEATURES REPEATED THORUG MULTIPLE PATHWAYS AND ADD TO FINAL DATAFRAME
  if (length(repeatedGenes > 0))
  {
    for (i in 1:length(repeatedGenes)){
      GeneScores = MixedFinalFrame[MixedFinalFrame$gene == repeatedGenes[i],]$score                    # Get all scores for repeated feature i
      score = sum(GeneScores)                                                                          # Sum all scores for repeated feature i
      RepeatedDF = data.frame('gene' = repeatedGenes[i], 'score' = score)                              # Generate data frame for gene feature i containing its name and its operated score
      FinalFrame = rbind(FinalFrame, RepeatedDF)                                                       # Bind operated feature i dataframe to Final Datframe.
    }
  }
  
  # ADD THE NON REPEATED FEATURES TO FINALFRAME AND ATTACH THIS AND THE ORIGINAL MIXED DATAFRAME TO LISTS
  if (length(MixedFinalFrame) > 0){
    NonRepeatedDF = subset(MixedFinalFrame[as.character(MixedFinalFrame$gene) %in% NonRepeatedGenes,], select = -c(path) )    # Extract a gene-score dataframe containing only the original non-repeated features
  }
  FinalFrame = rbind(FinalFrame, NonRepeatedDF)                                                      # Bind the non-repeated features dataframe to the repeated-features dataframe 
  #InstancesList[[Instance]][["MixedFinalFrame"]] = MixedFinalFrame                                   # Append the current instance's "MixedFinalFrame" to its index in ListGenes
  InstancesList[[Instance]][["FinalFrame"]] = FinalFrame                                             # Append the current instance's "FinalFrame" to its index in ListGenes
  InstancesList[[Instance]][["GraphInfo"]] = GraphInfoList
  InstancesList[[Instance]][["InstanceNames"]] = InstanceNames 
  print(FinalFrame)
}



InstancesList[[10]]$FinalFrame$gene
Features = c()
# VISUALIZING LIST PATHWAYS LENGTHS
for (i in 1:length(Instances)){
  Features = unique(c(Features, as.character(InstancesList[[i]]$FinalFrame$gene)))
}

# ------------------------------- PATHS' GENES FEATURES DATAFRAME ----------------------------------------
FeatureValue = matrix(0, nrow=length(Instances), ncol = length(Features))
for (i in 1:length(Instances)){
   LogicFeatures = Features %in% InstancesList[[i]]$FinalFrame$gene
  for (j in 1:length(Features)){
    if(LogicFeatures[j] == TRUE)
      FeatureValue[i,j] = InstancesList[[i]]$FinalFrame$score[InstancesList[[i]]$FinalFrame$gene == Features[j]]
    else
      FeatureValue[i,j] = 0
  }
}

# GENERATE GENES FEATURE DATAFRAME AND WRITE IT
FeatureValue = as.data.frame(FeatureValue)
GeneFeatures = RenameFeatures(Features)
row.names(FeatureValue) = InstancesLabel    # Instances
colnames(FeatureValue) = GeneFeatures       # Features
Class = row.names(FeatureValue) %in% CRgenes
FeatureValue$Class = Class
KeggGenes = FeatureValue                   # Define the KeggGenes dataset


# -------------------------------------- PATHS' FEATURES DATAFRAME --------------------------------------
# GENERATE PATHS FEATURES DATAFRAME AND WRITE IT
PathFeatures = unique(DfAllKegg[DfAllKegg$GeneID %in% Instances,]$PathwayID)

# GENERATE LIST OF PATHS FOR EVERY INSTANCE
GenePathsList = list()
for (i in 1:length(Instances)) {
  Instance = Instances[i]                                                               # Define Instance
  GeneKeggFrame = DfAllKegg[DfAllKegg$GeneID == Instance,]                                           # Retrieve All Pathways and Genes related to this Gene
  GenePathsList[[Instance]] = GeneKeggFrame$PathwayID   
}

# CREATE MATRIX CONTAINING THE PATHS DATA FOR INSTANCE
PathFeatureValue = matrix(0, nrow=length(Instances), ncol = length(PathFeatures))
for (i in 1:length(Instances)){
  LogicFeatures = PathFeatures %in% GenePathsList[[i]]    ############### GenePath
    PathFeatureValue[i,] = PathFeatures %in% GenePathsList[[i]]
}

# GENERATE PATHS FEATURE DATAFRAME AND WRITE IT
PathFeatureValue = as.data.frame(PathFeatureValue)
#GeneFeatures = RenameFeatures(Features)
row.names(PathFeatureValue) = InstancesLabel    # Instances
colnames(PathFeatureValue) = PathFeatures       # Features
Class = row.names(PathFeatureValue)
PathFeatureValue$Class = Class                 
KeggPaths = PathFeatureValue                        # Define the KeggPaths dataset


# SAVE DATASETS
write.csv(KeggPaths,"/ML-based-prediction-of-CR-related-genes/Data/Datasets/KeggPaths.csv")
write.csv(KeggGenes,"/ML-based-prediction-of-CR-related-genes/Data/Datasets/KeggGenes.csv")

# ###################################################################################################################
# ################################# INFLUENCE METHOD 1 (INTERMEDIATES DELETION) ##################################### 
# ################################# RECEIVE INSTANCE AND PATHWAY FOR INFLUENCE  #####################################
# ###################################################################################################################

InfluenceByDeletion <- function(graphAndMore, Instance) {
  
  g = graphAndMore[[1]]                                                                               # Get the graph from list
  GraphInfo = graphAndMore[[2]]
  gGeneMap = graphAndMore[[3]]
  
  for (i in 1:length(V(g))){
    if (Instance %in% strsplit(names(V(g))[i], '~')[[1]]) {                                                 # If instance is within hsa Entries IDs
      Instance = names(V(g))[i]
    }
  }

  # DESCENDANTS
  DescInstanceInclusive = names(subcomponent(g, Instance, mode = "out"))                              # Get all Instance's descendants, including itself
  DescInstance = DescInstanceInclusive[-1]                                                            # Exlude the same Instance in the list of descendants
  NumDescInstance = length(DescInstance)                                                              # Get number of Intsance's descendants
  DescInstanceInclusiveMap = names(subcomponent(gGeneMap, Instance, mode = "out"))                              # Get all Instance's descendants, including itself
  DescInstanceMap = DescInstanceInclusiveMap[-1]                                                            # Exlude the same Instance in the list of descendants
  DescFrame = data.frame('Names' = DescInstanceMap, 'Relation' = rep("Des",length(DescInstanceMap)))
  
  # ANCESTORS
  AscInstanceInclusive = names(subcomponent(gGeneMap, Instance, mode = "in"))                                # Get all Instance's ascendants, including itself 
  AscInstance = AscInstanceInclusive[-1]                                                              # Get all Instance's ascendants, including itself
  AscInstanceInclusiveMap = names(subcomponent(gGeneMap, Instance, mode = "in"))                                # Get all Instance's ascendants, including itself 
  AscInstanceMap = AscInstanceInclusiveMap[-1]                                                              # Get all Instance's ascendants, including itself
  AscFrame = data.frame('Names' = AscInstanceMap, 'Relation' = rep("Asc",length(AscInstanceMap)))
  
  # ALL
  AscDescInstance = unique(c(DescInstanceInclusive,AscInstanceInclusive))                             # Get all ascendants and descendants of Instance
  AscDescInstanceMap = unique(c(DescInstanceInclusiveMap,AscInstanceInclusiveMap))                             # Get all ascendants and descendants of Instance
  
  # NEIGHBOORS
  InstanceTree = c()
  if (length(AscInstanceMap) > 0){
    for ( j in 1:length(AscInstanceMap)){
      AscDown = names(subcomponent(gGeneMap, AscInstanceMap[j], mode = "out"))
      InstanceTree = unique(c(InstanceTree, AscDown))
    }
  } else {
    InstanceTree = c()
  }
  FarNeighboors = setdiff(InstanceTree, AscDescInstanceMap)
  FarNeighboorsFrame = data.frame('Names' = FarNeighboors, 'Relation' = rep("Ngb", length(FarNeighboors)))
  
  # INSTANCE 
  InstanceFrame = data.frame('Names' = Instance, 'Relation' = rep("INS",length(Instance)))
  
  # OTHER (Not including Neighboors)
  OtherNoNgb = setdiff(names(V(gGeneMap)), unique(c(AscDescInstanceMap, FarNeighboors)))
  OtherNoNgbFrame = data.frame('Names' = OtherNoNgb, 'Relation' = rep("Oth",length(OtherNoNgb)))
  
  # OTHER (Including NeighboorS)
  Others = setdiff(names(V(g)), AscDescInstance)
  
  # BIND ITERACTION DATAFRAMES ()
  InstanceRelationsFrame = do.call("rbind", list(DescFrame, AscFrame, FarNeighboorsFrame, InstanceFrame, OtherNoNgbFrame))
  
  #print("Hola mundo 2")
  FromRel = c()
  ToRel = c()
  EdgeRel = c()
  
  i = 1
  FromRel[i] = unique(as.character(InstanceRelationsFrame[InstanceRelationsFrame$Names %in% as.character(GraphInfo$From)[i],]$Relation))
  
  
  for (i in 1: length(GraphInfo$From)){                                                   # For edge type
    #print(i)
    FR = unique(as.character(InstanceRelationsFrame[InstanceRelationsFrame$Names %in% as.character(GraphInfo$From)[i],]$Relation))  # From type ########################################
    if (length(FR) == 1)
      FromRel[i] = FR
    else
      FromRel[i] = "AD"
    TR =   unique(as.character(InstanceRelationsFrame[InstanceRelationsFrame$Names %in% as.character(GraphInfo$To)[i],]$Relation))    # To Type ##########################################
    if (length(TR) == 1)
      ToRel[i] = TR
    else
      ToRel[i] = "AD"
    EdgeRel[i] = paste( c(FromRel[i], ToRel[i]), collapse="->")                                    # Bind source and destination types
  }
  EdgeRels = data.frame(EdgeRel)                                                                   # Generate Multiplicity and Loop Edge info Frame
  
  
  # ----------------------------------- DOWN INFLUENCE CALCULATION -----------------------------------
  Influence = c()
  #print(DescInstance)                                                                                   #####################################################################################################################
  if (length(DescInstance) > 0) {
    # DOWN INFLUENCE (TRASCENDENT)
    for (i in 1:length(DescInstance)){                                                                 # For each Feature (Instance descendant)
      #i = 8
      DescFeature = DescInstance[i]                                                                    # Get current downstream feature
      DescFeature
      AscFeatureInclusive = names(subcomponent(g, DescFeature, mode = "in"))                           # Get ancsestors for the current downstream feature
      AscFeature = AscFeatureInclusive[-1]                                                             # Get feature ancestors without the feature includig itself
      
      IntermediateNodes = intersect(DescInstanceInclusive, AscFeature)                                 # Get intermediate nodes between Feature and Instance 
      IntermediateNodes
      gd = delete_vertices(g, IntermediateNodes)                                                       # Delete intermediate nodes
      AscFeatureInclusiveAfter = names(subcomponent(gd, DescFeature, mode = "in"))                     # Get ancsestors of the current downstream feature that are not intermediate to the instance
      AscFeatureAfter = AscFeatureInclusiveAfter[-1]                                                     # Exclude the downstream gene itself from the non-intermediate ancestors list of the feature
      
      # INFLUENCE CALCULATION
      NumAscFeatureAfter = length(AscFeatureAfter)                                                     # Get num of non downstream asncestors
      Influence[i] = 1/(1+NumAscFeatureAfter)                                                          # Calculate influence
      Influence[i]
    }
    DownFeatures = data.frame(DescInstance, Influence)                                                 # Generate DataFrame of the features ans its score
  } else{
    DownFeatures = data.frame()                                                                        # Void dataframe if no descendants
  }
  
  # ------------------------------------ UP INFLUENCE CALCULATION ------------------------------------
  Influence = c()
  if (length(AscInstance) > 0) {
    # UP INFLUENCE (TRASCENDENT)
    for (i in 1:length(AscInstance)){                                                                  # For each Feature (Instance ascendat)                                                    
      Influence[i] = 0                                                                                 # Calculate up influence
    }
    UpFeatures = data.frame(AscInstance, Influence)                                                    # Generate DataFrame of the features ans its score
  } else{
    UpFeatures = data.frame()                                                                          # Void dataframe if no others
  }
  
  # ----------------------------------- OTHER INFLUENCE CALCULATION -----------------------------------
  Influence = c()
  if (length(Others) > 0) {
    # UP INFLUENCE (TRASCENDENT)
    for (i in 1:length(Others)){                                                                      # For each Feature (Instance ascendat)                                                    
      Influence[i] = 0                                                                                # Calculate other influence
    }
    OtherFeatures = data.frame(Others, Influence)                                                     # Generate DataFrame of the features ans its score
  } else{
    OtherFeatures = data.frame()                                                                      # Void dataframe if no others
  }
  
  # ------------------------------------ SELF INFLUENCE CALCULATION ------------------------------------
  SelfFeature = data.frame(Instance, 0) 
  
  ########################################## DECIDE WHAT TO RETURN #########################################
  
  # ONLY DOWN NODES
  features = DownFeatures
  
  FeaturesAndGraphInfo = list(features, EdgeRels, Instance)
  
  return(FeaturesAndGraphInfo)
}


# #########################################################################################################
# ############################################ FUNCTION ###################################################
# ############################### RECEIVE PATHWAY AND INSTANCE TO GRAPH ###################################
# #########################################################################################################

`%notin%` <- Negate(`%in%`)                                                                           # Declare function notin: to select elemnts that are not in a set
NNodes = c()

GenerateGraph <- function(Instance, Pathway) {
  
  KGML <- get_KGML(Pathway)
  
  if (!is.logical(KGML)){
    # OBTAIN GRAPH NODES AND EDGES
    Ngraph = graph::nodes(KGML)
    Egraph = graph::edges(KGML)
    
    # CREATE NODES' LISTS AND DATAFRAMES
    Entry = c()
    NodeAttList = list()
    NodeAttListKey = list()
    Nlist = list()
    Label = c()
    Type = c()
    DisplayNames = list()
    Names = c()
    Llist = list()
    
    for (i in 1:length(names(Ngraph))) {                                                                # For each of the Ngraph nodes
      Entry[i] = names(Ngraph)[i]                                                                       # Get ith node entry
      DisplayNames[i] = strsplit( str_remove_all(getDisplayName(Ngraph[[i]]), "[.]") , ", ")            # Get all the display names associated to ith node
      Label[i] = DisplayNames[[i]][1]                                                                   # Get ith node label
      Names[i] = paste(getName(Ngraph[[i]]), collapse='~')                                              # Get Name: eg.. "hsa:217~hsa:219"
      Type[i] = getType(Ngraph[[i]])                                                                    # Get type: eg.. "gene", "compound"
      NodeAttList[[i]] = list(Entry = Entry[i], DisplayNames = DisplayNames[[i]], Label = Label[i], Names = Names[i], Type = Type[i])             # Generate list with all attributes
      NodeAttListKey[[Entry[i]]] = list(Entry = Entry[i], DisplayNames = DisplayNames[[i]], Label = Label[i], Names = Names[i], Type = Type[i])   # Generate list with all attributes accesable from entry key
      Nlist[[Entry[i]]] = Names[i]    #list(Label[i])[[1]]   
      Llist[[Entry[i]]] = Label[i]                                                                      # Generate list of label accesable form Entry Key
    }
    
    
    DfNodes = data.frame(Entry, Label, Type, Names)                                                    # Generate Dataframe containing the nodes' entries, label and type
    DfNodes

    # CREATE EDGES' DATAFRAMES
    SourceNode = c()
    GoalNode = c()
    SourceNodeLabel = c()
    GoalNodeLabel = c()
    SourceNodeNames = c()
    GoalNodeNames = c()
    FromType = c()
    ToType = c()
    EdgeType = c()
    
    for (i in 1:length(names(Egraph))) {                                                                # For each Edge
      EdgeNodes = strsplit(getName(Egraph[[i]]), "~")[[1]]                                              # Get the entries related to the two nodes of the edge
      SourceNode[i] = EdgeNodes[1]                                                                      # Get the source node
      GoalNode[i] = EdgeNodes[2]                                                                        # Get the destination node
      SourceNodeLabel[i] = Llist[[EdgeNodes[1]]]                                                        # Get the label of the source node
      GoalNodeLabel[i] = Llist[[EdgeNodes[2]]]                                                          # Get the label of the destination node
      SourceNodeNames[i] = Nlist[[EdgeNodes[1]]]                                                        # Get the label of the source node
      GoalNodeNames[i] = Nlist[[EdgeNodes[2]]]                                                          # Get the label of the destination node
      }
    
    DfEdges = data.frame(SourceNode, GoalNode)                                                          # Create edges dataframe with entreis as nodes
    DfEdgesLabels = data.frame(SourceNodeLabel, GoalNodeLabel)                                          # create edges dataframe with labels as nodes
    DfEdgesNames = data.frame(SourceNodeNames, GoalNodeNames)                                           # create edges dataframe with labels as nodes
    
    # CREATE GRAPH FROM NODES' AND EDGES' DATA
    nodes <- data.frame(name=unique(Names))                                                             # Generate all the label vertices (nodes)
    from = SourceNodeNames                                                                              # Generate a sources nodes vector
    to = GoalNodeNames                                                                                  # Generate a destination nodes vector
    relations <- data.frame(from, to)                                                                   # Append the sorces and destination nodes in a dataframe
    G <- graph.data.frame(relations, directed=TRUE, vertices=nodes)                                     # Generate non-simplified graph of the pathway
    IntermediateRelations = ddply(relations,.(from,to),nrow)                                            # Get Intermediater dataframe for edge analysis
    multiplicity = IntermediateRelations$V1                                                             # Extract edges multiplicity from the intermediatr dataframe
    loop <- as.character(IntermediateRelations$from) == as.character(IntermediateRelations$to)          # Get edges loops vector
    for (i in 1: length(IntermediateRelations$from)){                                                   # For edge type
      FromType[i] = unique(as.character(DfNodes[DfNodes$Names == as.character(IntermediateRelations$from)[i],]$Type))  # From type
      ToType[i] = unique(as.character(DfNodes[DfNodes$Label == as.character(IntermediateRelations$to)[i],]$Type))    # To Type
      EdgeType[i] = paste( c(FromType[i], ToType[i]), collapse="->")                                    # Bind source and destination types
    }
    EdgeInfo = data.frame(EdgeType, 'From' = as.character(IntermediateRelations$from), 'To' = as.character(IntermediateRelations$to), 'Multiplicity' = as.character(multiplicity), 'Loop' = loop) # Generate Multiplicity and Loop Edge info Frame
    
    g = igraph::simplify(G, remove.multiple=TRUE, remove.loops=TRUE)     
    
    # SELECT ONLY GENE-NODES IN THE GRAPH
    ToKeepNodes = as.character(DfNodes[DfNodes[,3] %in% "gene",]$Names)                                 # Look for gene-type nodes 
    ToDeleteNodes = as.character(DfNodes[DfNodes[,3] != "gene",]$Names)                                 # Look for non-gene-type nodes
    gGenes = delete_vertices(g, ToDeleteNodes)                                                          # Delete non-gene nodes
    
    for (i in 1:length(V(gGenes))){
      if (Instance %in% strsplit(names(V(gGenes))[i], '~')[[1]]) {                                      # If instance is within hsa Entries IDs
        #print(i)
        Instance = names(V(gGenes))[i]
      }
    }
    
    # SELECT CONNECTED NODES IN THE GRAPH
    ToKeepNodes = unique((c(Instance, as.character(EdgeInfo$From), as.character(EdgeInfo$To))))
    
    #unique(as.character(DfNodes[DfNodes[,3] %in% c("gene", "map"),]$Label))                           # Look for gene-type nodes 
    ToDeleteNodes = setdiff(names(V(g)), ToKeepNodes)                                                  # as.character(DfNodes[DfNodes[,3] %notin% c("gene", "map"),]$Label)                                 # Look for non-gene-type nodes
    gMapsGenes = delete_vertices(g, ToDeleteNodes)                                                     # Delete non-gene nodes
    # plot(gd, layout = layout_as_star(gd))                                                     
    
    graphAndMore = list('Graph' = gGenes, 'Info' = EdgeInfo, 'MapGeneGraph' = gMapsGenes)
    return(graphAndMore)                                                                               # Return the simplified graph
  } else {
    return("Esto valio madre seniores")
  }
}


# #########################################################################################################
# ############################################ FUNCTION ###################################################
# ############################### RECEIVE DATAFRAME AND RETURN BETTER DF ##################################
# #########################################################################################################


CompleteFrame <- function(Instance, FinalFrame, DfAllKegg, InstanceLabel) {
  # SPLIT THE FEATURES NAME
  GeneList = list()
  BindedGeneKeggId = FinalFrame$gene
  GeneKeggIdList = strsplit(as.character(FinalFrame$gene), '~')
  BindedGeneList = list()
  MixedGeneScoresFrame = data.frame()
  Label = c()
  LabelId = c()
  GeneScoresFrame = data.frame()
  
  # CREATE COMPLETE DATAFRAME
  for( i in 1:length(GeneKeggIdList)){
    GeneList[[i]] = unique(DfAllKegg[ DfAllKegg$GeneID %in% GeneKeggIdList[[i]], ]$Symbol)
    Label[i] = GeneList[[i]][1]
    LabelId[i] = GeneKeggIdList[[i]][1]
    BindedGeneList[[i]] = paste(GeneList[[i]], collapse =  '~')
    GeneScores = data.frame('gene' = GeneList[[i]], 'score' = rep(FinalFrame$score[i], length(GeneList[[i]])))
    MixedGeneScoresFrame = rbind(MixedGeneScoresFrame, GeneScores)
  }
  
  # COMBINE REPETAED GENES IN MixedGeneScoresFrame
  
  repeatedGenes = unique(as.character(MixedGeneScoresFrame$gene[duplicated(MixedGeneScoresFrame$gene)]))       # Look for repeated gene features in the gene instance
  NonRepeatedGenes = setdiff(as.character(unique(MixedGeneScoresFrame$gene)), repeatedGenes)                             # Extract non repeated gene features in the gene instance
  #FinalFrame = data.frame()                                                                          # Create FinalDataFrame variable
  
  # COMBINE THE FEATURES REPEATED THORUG MULTIPLE PATHWAYS AND ADD TO FINAL DATAFRAME
  if (length(repeatedGenes > 0))
  {
    for (i in 1:length(repeatedGenes)){
      GeneScores = MixedGeneScoresFrame[as.character(MixedGeneScoresFrame$gene) %in% repeatedGenes[i],]$score                    # Get all scores for repeated feature i
       score = sum(GeneScores)                                                                          # Sum all scores for repeated feature i
       RepeatedDF = data.frame('gene' = repeatedGenes[i], 'score' = score)                              # Generate data frame for gene feature i containing its name and its operated score
       GeneScoresFrame = rbind(GeneScoresFrame, RepeatedDF)                                                       # Bind operated feature i dataframe to Final Datframe.
    }
  }
  
    # ADD THE NON REPEATED FEATURES TO FINALFRAME AND ATTACH THIS AND THE ORIGINAL MIXED DATAFRAME TO LISTS
  NonRepeatedDF = MixedGeneScoresFrame[as.character(MixedGeneScoresFrame$gene) %in% NonRepeatedGenes,]   # Extract a gene-score dataframe containing only the original non-repeated features
  GeneScoresFrame = rbind(GeneScoresFrame, NonRepeatedDF) 
  row.names(GeneScoresFrame) = 1:length(GeneScoresFrame[,1])
  
  #DfAllKegg[ DfAllKegg$GeneID == GeneKeggIdList[[1]], ]$Symbol
  #AllFinalFrame = 1
  
}


##################################################################################################################################


RenameFeatures <- function(Features) {
  GeneKeggIdList = strsplit(as.character(Features), '~')
  BindedGeneList = list()
  Label = c()
  LabelId = c()
  GeneScoresFrame = data.frame()
  GeneList = list()
  
  # CREATE COMPLETE DATAFRAME
  for( i in 1:length(GeneKeggIdList)){
    GeneList[[i]] = unique(DfAllKegg[ DfAllKegg$GeneID %in% GeneKeggIdList[[i]], ]$Symbol)
    Label[i] = GeneList[[i]][1]
    LabelId[i] = GeneKeggIdList[[i]][1]
    BindedGeneList[[i]] = paste(GeneList[[i]], collapse =  '~')
  }
  return(BindedGeneList)
}

InstancesList$ADH1A$GraphInfo[[3]]$Info
InstancesList$ADH1A$hsa00350



##### GET PATHWAY GENES ############################################################################################################

KGML <- get_KGML("hsa04730")
Ngraph = graph::nodes(KGML)
Egraph = graph::edges(KGML)

Entry = c()
NodeAttList = list()
NodeAttListKey = list()
Nlist = list()
Label = c()
Type = c()
DisplayNames = list()
Names = c()
Llist = list()

for (i in 1:length(names(Ngraph))) {                                                                # For each of the Ngraph nodes
  Entry[i] = names(Ngraph)[i]                                                                       # Get ith node entry
  DisplayNames[i] = strsplit( str_remove_all(getDisplayName(Ngraph[[i]]), "[.]") , ", ")            # Get all the display names associated to ith node
  Label[i] = DisplayNames[[i]][1]                                                                   # Get ith node label
  Names[i] = paste(getName(Ngraph[[i]]), collapse='~')                                              # Get Name: eg.. "hsa:217~hsa:219"
  Type[i] = getType(Ngraph[[i]])                                                                    # Get type: eg.. "gene", "compound"
  NodeAttList[[i]] = list(Entry = Entry[i], DisplayNames = DisplayNames[[i]], Label = Label[i], Names = Names[i], Type = Type[i])             # Generate list with all attributes
  NodeAttListKey[[Entry[i]]] = list(Entry = Entry[i], DisplayNames = DisplayNames[[i]], Label = Label[i], Names = Names[i], Type = Type[i])   # Generate list with all attributes accesable from entry key
  Nlist[[Entry[i]]] = Names[i]    #list(Label[i])[[1]]   
  Llist[[Entry[i]]] = Label[i]                                                                      # Generate list of label accesable form Entry Key
}

Label[Label %in% AgeGenes]
Label[Label %in% CrGenes]