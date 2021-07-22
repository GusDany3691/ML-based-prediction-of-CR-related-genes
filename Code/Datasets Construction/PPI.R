library(igraph)
library(dplyr) 
library(expm)
library(textshape)

install.packages('textshape')


#######################################################################################################################################################################
##### PREPARE GRAPH ###################################################################################################################################################
#######################################################################################################################################################################

# READ PPI FILE
Biogrid=read.delim('C:/Users/gusdany/Desktop/El Proyecto/Files/BIOGRID-MV-Physical-3.5.181.tab2.txt'
                   , sep='\t',header = TRUE,stringsAsFactors = FALSE)%>%dplyr::select(Official.Symbol.Interactor.A
                   , Official.Symbol.Interactor.B, Experimental.System, Experimental.System.Type, Organism.Interactor.A
                   , Organism.Interactor.B,Pubmed.ID,Throughput,Source.Database)

# READ PPI FILE
Biogrid=read.delim('M:/gusdany/Documents/Desmadre/BIOGRID-MV-Physical-4.2.191.tab2.txt'
                   , sep='\t',header = TRUE,stringsAsFactors = FALSE)%>%dplyr::select(Official.Symbol.Interactor.A
                   , Official.Symbol.Interactor.B, Experimental.System, Experimental.System.Type, Organism.Interactor.A
                   , Organism.Interactor.B,Pubmed.ID,Throughput,Source.Database)


#:/Users/gusdany/Documents/Desmadre/BIOGRID-MV-Physical-4.2.191.mitab.txt

Biogrid=read.delim('C:/Users/gusdany/Documents/Desmadre/BIOGRID-MV-Physical-4.2.191.mitab.txt'
                   , sep='\t',header = TRUE,stringsAsFactors = FALSE)

#BIOGRID-MV-Physical-4.2.191.mitab.txt

'%!in%' <- function(x,y)!('%in%'(x,y))

# FILTER FOR HUMAN GENES AND CREATE GRAPH
HumanBiogridFrame = Biogrid%>%dplyr::filter(Organism.Interactor.A==9606,Organism.Interactor.B==9606)
HumanBiogridGraph = graph.data.frame(data.frame(Biogrid$Official.Symbol.Interactor.A, Biogrid$Official.Symbol.Interactor.B),directed = FALSE)

# READ AGEING AND CR GENES
# GET AGEING AND CR GENES
AgeGenes = readRDS("D:/gusdany/Documents/AgeingFinal.rds")
CrGenes = readRDS("D:/gusdany/Documents/CrFinal.rds")#("C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/CrGenes.rds")

# MAKE A VECTOR OF ALL GENES IN BIOGRID
HumanBiogridGenes = unique(c(HumanBiogridFrame$Official.Symbol.Interactor.A, HumanBiogridFrame$Official.Symbol.Interactor.B))

# KEEP AGEINS AND CR GENES THAT ARE IN BIOGRID
BiogridAgeGenes = intersect(AgeGenes, HumanBiogridGenes)
BiogridCrGenes  = intersect(CrGenes, HumanBiogridGenes)

# SIMPLIFY GRAPH
sHumanBiogridGraph = igraph::simplify(HumanBiogridGraph, remove.multiple=TRUE, remove.loops=TRUE)

########################################################################################################################################################################
##### GRAPH MEASURES ###################################################################################################################################################
########################################################################################################################################################################

g = sHumanBiogridGraph
Nodes = BiogridAgeGenes

##### CENTRALITIES #####################################################################################################################################################

# --- LOAD PROVIUS MEASURES (IF SAVED) ---------------------------------------------------------------------------------------------------------------------------------

Degree = readRDS("C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/PPI/Degree.rds")
Betweenness = readRDS("C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/PPI/Betweenness.rds")
Closeness = readRDS("C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/PPI/Closeness.rds")
EigenCentrality = readRDS("C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/PPI/EigenCentrality.rds")
Eccentricity = readRDS("C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/PPI/Eccentricity.rds")

SubgraphCentrality = readRDS("C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/PPI/SubgraphCentrality.rds")

#Flowbet = readRDS("C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/PPI/Flowbet.rds")
#Gilschmidt = readRDS("C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/PPI/Gilschmidt.rds")
#LoadCent = readRDS("C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/PPI/LoadCent.rds")
#Infocent = readRDS("C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/PPI/Infocent.rds")
#Stresscent = readRDS("C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/PPI/Stresscent.rds")

#Averagedis = readRDS("C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/PPI/Averagedis.rds")
#Barycenter = readRDS("C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/PPI/Barycenter.rds")
#ClosenessCurrentFlow = readRDS("C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/PPI/ClosenessCurrentFlow.rds")
ClosenessLatora = readRDS("C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/PPI/ClosenessLatora.rds")
ClosenessResidual = readRDS("C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/PPI/ClosenessResidual.rds")
#Communibet = readRDS("C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/PPI/Communibet.rds")
#Crossclique = readRDS("C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/PPI/Crossclique.rds")
Decay = readRDS("C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/PPI/Decay.rds")
Diffusion = readRDS("C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/PPI/Diffusion.rds")
#*Entropy = readRDS("C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/PPI/Entropy.rds")
Geokpath  = readRDS("C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/PPI/Geokpath.rds")
Laplacian  = readRDS("C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/PPI/Laplacian.rds")
Leverage  = readRDS("C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/PPI/Leverage.rds")
Licnent  = readRDS("C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/PPI/Licnent.rds")
Lobby  = readRDS("C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/PPI/Lobby.rds")
Markovcent  = readRDS("C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/PPI/Markovcent.rds")
Mnc  = readRDS("C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/PPI/Mnc.rds")
Radiallity  = readRDS("C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/PPI/Radiallity.rds")
Semilocal  = readRDS("C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/PPI/Semilocal.rds")
Topocoefficient = readRDS("C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/PPI/Topocoefficient.rds")

#DangalchevClosenes = readRDS("C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/PPI/DangalchevClosenes.rds")
#* Harmonic = readRDS("C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/PPI/Harmonic.rds")
#* LocalBriding = readRDS("C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/PPI/LocalBriding.rds")

# --- COMPUTE CENTRALITIES -------------------------------------------------------------------------------------------------------------------------------------------------

Degree <- igraph::degree(g, v = Nodes)                         # Independient and name
saveRDS(Degree,"Degree.rds")
print(1)
Betweenness <- igraph::betweenness(g, v = Nodes)               # Independient and name
saveRDS(Betweenness,"Betweenness.rds")
print(2)
Closeness <- igraph::closeness(g, v = Nodes)                   # Independient and name
saveRDS(Closeness,"Closeness.rds")
print(3)
EigenCentrality <- igraph::eigen_centrality(g)$vector          # All and Name 
saveRDS(EigenCentrality,"EigenCentrality.rds")
print(4)
Eccentricity <- 1/igraph::eccentricity(g, v = Nodes)           # Independient and name
saveRDS(Eccentricity,"Eccentricity.rds")
print(5)
SubgraphCentrality <- igraph::subgraph_centrality(g)           # All and Name  
saveRDS(SubgraphCentrality,"SubgraphCentrality.rds")
print(6)

#A <- get.adjacency(g,sparse=F)
#Flowbet <- sna::flowbet(A)                                     # All and No Name
#saveRDS(Flowbet,"Flowbet.rds")
#print(7)
#LoadCent <- sna::loadcent(A)                                   # All and No Name
#saveRDS(LoadCent,"Loadcent.rds")
#print(8)
#Gilschmidt <- sna::gilschmidt(A)                               # All and No Name
#saveRDS(Gilschmidt,"Gilschmidt.rds")
#print(9)
#Infocent <- sna::infocent(A)                                   # All and No Name
#saveRDS(Infocent,"Infocent.rds")
#print(10)
#Stresscent <- sna::stresscent(A)                               # All and no Name
#saveRDS(Stresscent,"Stresscent.rds")
#print(11)

Averagedis <- 1/centiserve::averagedis(g, vids = Nodes)        # Independient and name                             
saveRDS(Averagedis,"Averagedis.rds")
print(12)  # Must be loopfree
Barycenter <- centiserve::barycenter(g, vids = Nodes)          # Independient and name
saveRDS(Barycenter,"Barycenter.rds")
print(13) # Must be loop free
ClosenessCurrentFlow <- centiserve::closeness.currentflow(g, vids = Nodes)     # Independient and name
saveRDS(ClosenessCurrentFlow,"ClosenessCurrentflow.rds")
print(14)
ClosenessLatora <- centiserve::closeness.latora(g, vids = Nodes)               # Independient and name
saveRDS(ClosenessLatora,"ClosenessLatora.rds")
print(15)
ClosenessResidual <- centiserve::closeness.residual(g, vids = Nodes)           # Independient and name
saveRDS(ClosenessResidual,"ClosenessResidual.rds")
print(16)
Communibet <- centiserve::communibet(g, vids = Nodes)          # Independient and name
saveRDS(Communibet,"Communibet.rds")
print(17)
Crossclique <- centiserve::crossclique(g, vids = Nodes)        # Independient and name
saveRDS(Crossclique,"Crossclique.rds")
print(18)
Decay <- centiserve::decay(g, vids = Nodes)                    # Independient and name
saveRDS(Decay,"Degree.rds")#
print(19)
Diffusion <- centiserve::diffusion.degree(g, vids = Nodes)     # Independient and name
saveRDS(Diffusion,"Decay.rds")#
print(20)
Entropy <- 1/centiserve::entropy(g, vids = Nodes)              # Independient and name
saveRDS(Entropy,"Entropy.rds")
print(21)
Geokpath <- centiserve::geokpath(g, vids = Nodes)              # Independient and name
saveRDS(Geokpath,"Geokpath.rds")
print(22)
Katzcent <- centiserve::katzcent(g, vids = Nodes)              # Independient and name
saveRDS(Katzcent,"Katzcent.rds")
print(23)
Laplacian <- centiserve::laplacian(g, vids = Nodes)            # Independient and name
saveRDS(Laplacian,"Laplacian.rds")
print(24)
Leverage <- centiserve::leverage(g, vids = Nodes)              # Independient and name
saveRDS(Leverage,"Leverage.rds")
print(25)
Licnent <- centiserve::lincent(g, vids = Nodes)                # Independient and name
saveRDS(Licnent,"Lincent.rds")
print(26)
Lobby <- centiserve::lobby(g, vids = Nodes)                    # Independient and name
saveRDS(Lobby,"Lobby.rds")
print(27)
Markovcent <- centiserve::markovcent(g, vids = Nodes)          # Independient and name
saveRDS(Markovcent,"Markovcent.rds")
print(28)
Mnc <- centiserve::mnc(g, vids = Nodes)                        # Independient and name
saveRDS(Mnc,"Mnc.rds")
print(29)
Radiality <- centiserve::radiality(g, vids = Nodes)            # Independient and name
saveRDS(Radiality,"Radiality.rds")
print(30) 
Semilocal <- centiserve::semilocal(g, vids = Nodes)            # Independient and name
saveRDS(Semilocal,"Semilocal.rds")
print(31)  
Topocoefficient <- 1/centiserve::topocoefficient(g, vids = Nodes)            # Independient and name
saveRDS(Topocoefficient,"Topocoefficient.rds")
print(32)

#DangalchevCloseness <- CINNA::dangalchev_closeness_centrality(g)             # All and name
#saveRDS(DangalchevCloseness,"DangalchevCloseness.rds")
#print(33)
#Harmonic <- CINNA::harmonic_centrality(g)                      # All and name
#saveRDS(Harmonic,"Harmonic.rds")
#print(34)
#LocalBriding <- 1/CINNA::local_bridging_centrality(g)          # All and name
#saveRDS(LocalBriding,"LocalBriding.rds")
#print(35)

# --- PRE PROCESSING --------------------------------------------------------------------------------------------------------------------------------------------------

SubgraphCentrality = SubgraphCentrality[Nodes]
EigenCentrality = EigenCentrality[Nodes]

##### OTHER MEASURES ##################################################################################################################################################

# --- UPLOAD OTHER MEASURES -------------------------------------------------------------------------------------------------------------------------------------------

Kcore = readRDS("C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/PPI/Kcore.rds")
ClusteringCoefficient = readRDS("C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/PPI/ClusteringCoefficient.rds")
CrRatio = readRDS("C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/PPI/CrRatio.rds")
AdjacencyMatrix = readRDS("C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/PPI/AgeAdjacencyMatrix.rds")

# --- COMPUTE OTHER MEASURES ------------------------------------------------------------------------------------------------------------------------------------------

# K CORE
HumanKcore = 'Kcore' = coreness(g, mode = c("all"))
Kcore = HumanKcore[Nodes]
saveRDS(Kcore, "Kcore.rds")#(),"C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/PPI/Kcore.rds")

# CLUSTERING COEFFICIENT
ClusteringCoefficientFrame = transitivity(as.undirected(g),"local", isolates = c("zero"), vids = Nodes)
ClusteringCoefficientFrame = transitivity(as.undirected(g),"local", isolates = c("zero"))
#ClusteringCoefficient = as.numeric(ClusteringCoefficientFrame$ClusteringCoefficient)
ClusteringCoefficient = as.numeric(ClusteringCoefficientFrame)
names(ClusteringCoefficientFrame) = Nodes
saveRDS(ClusteringCoefficient,"ClusteringCoefficient.rds")#,"C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/PPI/ClusteringCoefficient.rds")

# CR RATIO
CrRatio = data.frame()
BiogriCrGenes = CrGenes[CrGenes %in% names(V(g))]
for (i in 1:length(Nodes)){
  # Cr Ratio
  print(i)
  Node = Nodes[i]
  NeighborhoodGenes = neighborhood(graph=g, order=1, nodes=Node, mode = c("all"), mindist = 1)
  NeighborhoodGenes = names(NeighborhoodGenes[[1]])
  TotalNeighbors = length(NeighborhoodGenes)
  TotalCrNeighbors = sum(NeighborhoodGenes %in% CrGenes)
  Ratio = TotalCrNeighbors / TotalNeighbors
  cCrRatio = data.frame('Gene' = Node, 'CrRatio' = Ratio)
  CrRatioFrame = rbind(CrRatio, cCrRatio)
  
  # Other indices 
  #toLeaves = all_simple_paths(g, from = Nodes[1], to = CrGenes[1])
}
saveRDS(CrRatio,"C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/PPI/CrRatio.rds")

# --- ADJACENCIES TO CR AND AGE ---------------------------------------------------------------------------------------------------------------------------------------

ReachableCr = c() 
Adjacent1Cr = c() 
Adjacent2Cr = c() 
Adjacent3Cr = c() 
Adjacent4Cr = c() 
Adjacent5Cr = c()
Adjacent6Cr = c() 
Adjacent7Cr = c() 
Adjacent8Cr = c() 
Adjacent9Cr = c() 
Adjacent10Cr = c() 
MeanDistanceCr = c()

ReachableAge = c() 
Adjacent1Age = c() 
Adjacent2Age = c() 
Adjacent3Age = c() 
Adjacent4Age = c() 
Adjacent5Age = c()
Adjacent6Age = c() 
Adjacent7Age = c() 
Adjacent8Age = c() 
Adjacent9Age = c()
Adjacent10Age = c()
MeanDistanceAge = c()

for (i in 1:length(Nodes)){
  print(i)
  Node = Nodes[i]
  # CR PATHS
  ShortestPaths = get.shortest.paths(g,Node,BiogriCrGenes)$vpath
  PathsLength = c()
  for (j in 1:length(ShortestPaths))
    PathsLength = c(PathsLength, length(ShortestPaths[[j]])) 
  ReachableCr[i] = sum(PathsLength != 0)
  #Adjacent1Cr[i] = sum(PathsLength == 1)
  Adjacent1Cr[i] = sum(PathsLength == 2)
  Adjacent2Cr[i] = sum(PathsLength == 3)
  Adjacent3Cr[i] = sum(PathsLength == 4)
  Adjacent4Cr[i] = sum(PathsLength == 5)
  Adjacent5Cr[i] = sum(PathsLength == 6)
  Adjacent6Cr[i] = sum(PathsLength == 7)
  Adjacent7Cr[i] = sum(PathsLength == 8)
  Adjacent8Cr[i] = sum(PathsLength == 9)
  Adjacent9Cr[i] = sum(PathsLength == 10)
  MeanDistanceCr[i] = mean(PathsLength)
  
  # AGE PATHS
  ShortestPaths = get.shortest.paths(g,Node,Nodes)$vpath
  PathsLength = c()
  for (j in 1:length(ShortestPaths))
    PathsLength = c(PathsLength, length(ShortestPaths[[j]])) 
  ReachableAge[i] = sum(PathsLength != 0)
  #Adjacent1Age[i] = sum(PathsLength == 1)
  Adjacent1Age[i] = sum(PathsLength == 2)
  Adjacent2Age[i] = sum(PathsLength == 3)
  Adjacent3Age[i] = sum(PathsLength == 4)
  Adjacent4Age[i] = sum(PathsLength == 5)
  Adjacent5Age[i] = sum(PathsLength == 6)
  Adjacent6Age[i] = sum(PathsLength == 7)
  Adjacent7Age[i] = sum(PathsLength == 8)
  Adjacent8Age[i] = sum(PathsLength == 9)
  Adjacent9Age[i] = sum(PathsLength == 10)
  MeanDistanceAge[i] = mean(PathsLength)
  
}

CrNeighboursFrame = data.frame(ReachableCr, Adjacent1Cr, Adjacent2Cr, Adjacent3Cr, Adjacent4Cr, Adjacent5Cr, 
                               Adjacent6Cr, Adjacent7Cr, Adjacent8Cr, Adjacent9Cr, MeanDistanceCr)

AgeNeighboursFrame = data.frame(ReachableAge, Adjacent1Age, Adjacent2Age, Adjacent3Age, Adjacent4Age, Adjacent5Age, 
                                Adjacent6Age, Adjacent7Age, Adjacent8Age, Adjacent9Age, MeanDistanceAge)

CrRatio = as.numeric(CrRatioFrame$CrRatio)
names(CrRatio) = CrRatioFrame$Gene




# --- ADJACENCY MATRIX ---------------------------------------------------------------------------------------------------------------------
HumanAdjacencyMatrix = as_adjacency_matrix(g)
AdjencyMatrix = HumanAdjacencyMatrix[Nodes,]
AdjencyDataSet = as.data.frame(as.matrix(AdjencyMatrix))
saveRDS(AdjencyDataSet,"C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/PPI/AdjencyDataSet.rds")
#######################################################################################################################################################################
#####  CREATE DATABASES FRAMES ########################################################################################################################################
#######################################################################################################################################################################

Class = Nodes %in% CrGenes
Class = ifelse(Class == "TRUE", "CR", "NOtCR") 
Class = data.frame('Class' = Class)

###### NETWORK ATTRIBUTES #############################################################################################################################################

PPImeasuresFrame = data.frame(
                              'Degree' = Degree,
                              'Betweenness' = Betweenness,
                              'Closeness' = Closeness,
                              'EigenCentrality' = EigenCentrality,
                              'Eccentricity' = Eccentricity,
                              'SubgraphCentrality' = SubgraphCentrality,
                              
                              #'Flowbet' = Flowbet,
                              #'Gilschmidt' = Gilschmidt,
                              #'LoadCent' = LoadCent,
                              #'Infocent' = Infocent,
                              #'Stresscent' = Stresscent,
                              
                              #'Averagedis' = Averagedis,
                              #'Barycenter' = Barycenter,
                              #'ClosenessCurrentFlow' = ClosenessCurrentFlow,
                              #'ClosenessLatora' = ClosenessLatora,
                              #'ClosenessResidual' = ClosenessResidual,
                              #'Communibet' = Communibet,
                              #'Crossclique' = Crossclique,
                              #'Decay' = Decay,
                              'Diffusion' = Diffusion,
                              #'Entropy' = Entropy,
                              'Geokpath' = Geokpath,
                              'Laplacian' = Laplacian,
                              'Leverage' = Leverage,
                              #'Licnent' = Licnent,
                              'Lobby' = Lobby,
                              'Markovcent' = Markovcent,
                              'Mnc' = Mnc,
                              #'Radiallity' = Radiallity,
                              'Semilocal' = Semilocal,
                              'Topocoefficient' = Topocoefficient,
                              
                              #'DangalchevClosenes' = DangalchevClosenes,
                              #'Harmonic' = Harmonic,
                              #'LocalBriding' = LocalBriding,
                              
                              Kcore,
                              ClusteringCoefficient,
                              CrRatio
          
)

PPImeasuresDataset = cbind(PPImeasuresFrame, CrNeighboursFrame, AgeNeighboursFrame, Class)
saveRDS(PPImeasuresDataset,"C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/PPI/PPImeasuresDataset.rds")

#PPImeasures = readRDS("C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/PPI/PPImeasuresDataset.rds")
#Class = PPImeasures$Class
#Measures = PPImeasures[,-length(PPImeasures)]
#PPImeasuresDataset = cbind(Measures,CrNeighboursFrame, AgeNeighboursFrame, Class)
#saveRDS(PPImeasuresDataset,"C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/PPI/PPImeasuresDataset2.rds")

##### ADJENCY ATTRIBUTES ##################################################################################################################################################



PPIa = as_adjacency_matrix(g)
nPPIa = PPIa[row.names(PPIa) %in% AgeGenes,]
nPPIa = as.data.frame(as.matrix(nPPIa))

for (i in 1:length(nPPIa)){
  nPPIa[[i]] = ifelse(nPPIa[[i]] == 1, 1, 0)
}

UsefulCols = as.logical(colSums(nPPIa, na.rm = FALSE, dims = 1) > 0)

nPPIa = nPPIa[,UsefulCols]

Class = ifelse(row.names(nPPIa) %in% CrGenes,"CR", "NotCR")

nPPIa$Class = Class

PPIadjencyDataset = cbind(AdjencyDataSet, Class)
saveRDS(nPPIa,"D:/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/PPI/PPIadjencyDataset.rds")
write.csv(nPPIa,"D:/gusdany/Desktop/El Proyecto/Databases/DataSetsHagr/Final/OMA/PPIadjencyDataset.csv")
PPIadjencyDataset = read.csv("D:/gusdany/Desktop/El Proyecto/Databases/DataSetsHagr/Final/OMA/PPIadjencyDataset.csv")

PPIadjencyDataset = column_to_rownames(PPIadjencyDataset, loc = 1)

GoodFeatureIndex = c()
for (i in 1:(length(PPIadjencyDataset)-1)){
  NumFeaturedInstances = sum(PPIadjencyDataset[,i])
  if (NumFeaturedInstances >= 2)
    GoodFeatureIndex = c(GoodFeatureIndex,i)
}
length(GoodFeatureIndex)
GoodFeatureIndex = c(GoodFeatureIndex, length(PPIadjencyDataset))


PPIfilteredAdjencyDataSet = PPIadjencyDataset[,GoodFeatureIndex]

saveRDS(PPIfilteredAdjencyDataSet,"C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/PPI/PPIfilteredAdjencyDataSet.rds")
write.csv(PPIfilteredAdjencyDataSet, "C:/Users/gusdany/Desktop/El Proyecto/Databases/DataSetsHagr/PPIfilteredAdjencyDataSet.csv")

write.csv(nPPIa, "C:/Users/gusdany/Desktop/El Proyecto/Databases/DataSetsHagr/Final/PPIAdjencyDataSet.csv")
