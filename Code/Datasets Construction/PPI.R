library(igraph)
library(dplyr) 
library(expm)
library(textshape)

install.packages('textshape')


#######################################################################################################################################################################
##### PREPARE GRAPH ###################################################################################################################################################
#######################################################################################################################################################################

# READ PPI FILE
Biogrid=read.delim('/ML-based-prediction-of-CR-related-genes/Data/Datasets/External/BIOGRID-MV-Physical-4.2.191.mitab.txt'
                   , sep='\t',header = TRUE,stringsAsFactors = FALSE)


'%!in%' <- function(x,y)!('%in%'(x,y))

# FILTER FOR HUMAN GENES AND CREATE GRAPH
HumanBiogridFrame = Biogrid%>%dplyr::filter(Organism.Interactor.A==9606,Organism.Interactor.B==9606)
HumanBiogridGraph = graph.data.frame(data.frame(Biogrid$Official.Symbol.Interactor.A, Biogrid$Official.Symbol.Interactor.B),directed = FALSE)

# READ AGEING AND CR GENES
AgeGenes = readRDS("/ML-based-prediction-of-CR-related-genes/Data/Ageing and CR genes/AgeingGenes.rds")
CrGenes = readRDS("/ML-based-prediction-of-CR-related-genes/Data/Ageing and CR genes/CrGenes.rds")

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


# --- COMPUTE CENTRALITIES -------------------------------------------------------------------------------------------------------------------------------------------------

Degree <- igraph::degree(g, v = Nodes)                      
Betweenness <- igraph::betweenness(g, v = Nodes)               
Closeness <- igraph::closeness(g, v = Nodes)                 
EigenCentrality <- igraph::eigen_centrality(g)$vector         
Eccentricity <- 1/igraph::eccentricity(g, v = Nodes)          
SubgraphCentrality <- igraph::subgraph_centrality(g)          
Averagedis <- 1/centiserve::averagedis(g, vids = Nodes)                                 
Barycenter <- centiserve::barycenter(g, vids = Nodes)        
ClosenessCurrentFlow <- centiserve::closeness.currentflow(g, vids = Nodes)     
ClosenessLatora <- centiserve::closeness.latora(g, vids = Nodes)               
ClosenessResidual <- centiserve::closeness.residual(g, vids = Nodes)     
Communibet <- centiserve::communibet(g, vids = Nodes)       
Crossclique <- centiserve::crossclique(g, vids = Nodes)    
Decay <- centiserve::decay(g, vids = Nodes)                  
Diffusion <- centiserve::diffusion.degree(g, vids = Nodes)    
Entropy <- 1/centiserve::entropy(g, vids = Nodes)             
Geokpath <- centiserve::geokpath(g, vids = Nodes)            
Katzcent <- centiserve::katzcent(g, vids = Nodes) 
Laplacian <- centiserve::laplacian(g, vids = Nodes)   
Leverage <- centiserve::leverage(g, vids = Nodes)           
Licnent <- centiserve::lincent(g, vids = Nodes)            
Lobby <- centiserve::lobby(g, vids = Nodes)                   
Markovcent <- centiserve::markovcent(g, vids = Nodes)         
Mnc <- centiserve::mnc(g, vids = Nodes)                      
Radiality <- centiserve::radiality(g, vids = Nodes)         
Semilocal <- centiserve::semilocal(g, vids = Nodes)        
Topocoefficient <- 1/centiserve::topocoefficient(g, vids = Nodes)         

SubgraphCentrality = SubgraphCentrality[Nodes]
EigenCentrality = EigenCentrality[Nodes]

##### OTHER MEASURES ##################################################################################################################################################


# K CORE
HumanKcore = 'Kcore' = coreness(g, mode = c("all"))
Kcore = HumanKcore[Nodes]


# CLUSTERING COEFFICIENT
ClusteringCoefficientFrame = transitivity(as.undirected(g),"local", isolates = c("zero"), vids = Nodes)
ClusteringCoefficient = as.numeric(ClusteringCoefficientFrame)
names(ClusteringCoefficientFrame) = Nodes
saveRDS(ClusteringCoefficient,"ClusteringCoefficient.rds")

# CR RATIO
CrRatio = data.frame()
BiogriCrGenes = CrGenes[CrGenes %in% names(V(g))]
for (i in 1:length(Nodes)){
  print(i)
  Node = Nodes[i]
  NeighborhoodGenes = neighborhood(graph=g, order=1, nodes=Node, mode = c("all"), mindist = 1)
  NeighborhoodGenes = names(NeighborhoodGenes[[1]])
  TotalNeighbors = length(NeighborhoodGenes)
  TotalCrNeighbors = sum(NeighborhoodGenes %in% CrGenes)
  Ratio = TotalCrNeighbors / TotalNeighbors
  cCrRatio = data.frame('Gene' = Node, 'CrRatio' = Ratio)
  CrRatioFrame = rbind(CrRatio, cCrRatio)
}


###########################################################################################################################################################################
##### ADJENCY ATTRIBUTES ##################################################################################################################################################
###########################################################################################################################################################################

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


#######################################################################################################################################################################
#####  CREATE DATABASES FRAMES ########################################################################################################################################
#######################################################################################################################################################################

Class = Nodes %in% CrGenes
Class = ifelse(Class == "TRUE", "CR", "NOtCR") 
Class = data.frame('Class' = Class)

###### NETWORK ATTRIBUTES #############################################################################################################################################

PPImeasuresFrame = data.frame('Degree' = Degree,
                              'Betweenness' = Betweenness,
                              'Closeness' = Closeness,
                              'EigenCentrality' = EigenCentrality,
                              'Eccentricity' = Eccentricity,
                              'SubgraphCentrality' = SubgraphCentrality,
                              'Diffusion' = Diffusion,
                              'Geokpath' = Geokpath,
                              'Laplacian' = Laplacian,
                              'Leverage' = Leverage,
                              'Lobby' = Lobby,
                              'Markovcent' = Markovcent,
                              'Mnc' = Mnc,
                              'Semilocal' = Semilocal,
                              'Topocoefficient' = Topocoefficient,
                              Kcore,
                              ClusteringCoefficient,
                              CrRatio)

PPImeasuresDataset = cbind(PPImeasuresFrame, CrNeighboursFrame, AgeNeighboursFrame, Class)
PPIadjencyDataset = cbind(AdjencyDataSet, Class)

write.csv(PPImeasuresDataset,"/ML-based-prediction-of-CR-related-genes/Data/Datasets/PPImeasuresDataset.csv")
write.csv(PPIadjencyDataset,"/ML-based-prediction-of-CR-related-genes/Data/Datasets/PPIadjencyDataset.csv")
