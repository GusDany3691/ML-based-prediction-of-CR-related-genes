#LOADING AGEING AND CR GENES
AgeGenes = read.csv("E:/gusdany/Documents/AgeingOMA.csv")$AgeGenes
CrGenes = read.csv("E:/gusdany/Documents/CrOMA.csv")$CrGenes

# LOAD PATHDIP DATA
PathDipDataset = read.csv("C:/Users/gusdany/Desktop/El Proyecto/Databases/DataSetsHagr/PathDipOriginal.csv") 

# GET PATHNAMES
PathNames = c()
for (i in 1:length(PathDip)){
  Path = as.character(PathDip[1,i])
  PathNames = c(PathNames, Path)
}

# GET PATH TYPES
Colnames = PathDip[2,]

# COUPLE PATH TYPE AND PATH NAME
NewColnames = paste(Colnames,PathNames, sep = " -> ")

# GENERATE NEW FRAME WITH COLUMN NAMES CONSIDERING BOTH PATH TYPE AND NAME
NiceColsPathDip = PathDip
colnames(NiceColsPathDip) = NewColnames
NiceColsPathDip = NiceColsPathDip[-1,]
row.names(NiceColsPathDip) = NiceColsPathDip[[1]]
NiceColsPathDip[[1]] = NULL

Class = ifelse(NiceColsPathDip$X %in% CrGenes, "CR", "NotCR")

# FILL VOID WITH ZEROS AND Y's WITH ONES
NumRows = length(row.names(NiceColsPathDip))
NumCols = length(colnames(NiceColsPathDip))

NicePathDip =  data.frame()
  for (i in 1:NumRows){
    print(i)
    States = c()
    for (j in 1:NumCols){
      States[j] = ifelse(NiceColsPathDip[i,j] == "y", "1", "0")
    }
    cNicePathDip = as.data.frame(t(States))
    colnames(cNicePathDip) = colnames(NiceColsPathDip)
    row.names(cNicePathDip) = row.names(NiceColsPathDip)[i]
    NicePathDip = rbind(NicePathDip, cNicePathDip)
  }

PathDipDataset = cbind(NicePathDip, Class)
write.csv(PathDipDataset, "C:/Users/gusdany/Desktop/El Proyecto/Databases/DataSetsHagr/Final/PathDip.csv")
