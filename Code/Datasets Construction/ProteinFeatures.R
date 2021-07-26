library("ensembldb")
library("EnsDb.Hsapiens.v86")
library("Peptides")
library("stringr")
library("protr")

# ASSIGN THE HUMAN DATABASE TO A VARIABLE
edb <- EnsDb.Hsapiens.v86        

# EVALUATE WETHER IT EXISTS PROTEINT ANNOTATION AVAIALABLE
hasProteinData(edb)

#LOADING AGEING AND CR GENES
AgeingGenes = readRDS("/ML-based-prediction-of-CR-related-genes/Data/Ageing and CR genes/AgeingGenes.rds")
CrGenes = readRDS("/ML-based-prediction-of-CR-related-genes/Data/Ageing and CR genes/CrGenes.rds")

# LOAD AGEING GENES
InstanceGenes = AgeGenes

IntanceTranscriptsList = list()
IntanceTranscriptsPerGeneList = list()
Protenis = list()
for(i in 1:length(InstanceGenes)){                  # For all Age Genes
#for(i in 1:2){  
  TranscriptsList = list()
  # GET THE CURRENT GENE INSTANCE
  InstanceGene = InstanceGenes[i]                   # Get Gen name
  Class = InstanceGene %in% CrGenes                 # Get class
  
  # GET TRANSCRIPTS INFO (TRANSCRIPT ID, LENGHT AND SEQUENCE) FROM GENE NAME
  Transcripts <- proteins(edb, filter = GeneNameFilter(InstanceGene), return.type = "AAStringSet")  # Get transcripts from current gene
  
  # CREATE LIST OF TRANSCRIPT INSTANCES
  if (length(Transcripts) > 0)
  {
    for (j in 1:length(Transcripts)){               # For all Transcripts of current gene
      print(paste(i,j, sep = "-"))
      Transcript = Transcripts[[j]]                 # Select transcript
      Seq = as.character(Transcript)                # Get transcript sequence
      Seq <- Seq[(sapply(Seq, protcheck))]          # Quit non aminoacid terms
      if (length(Seq) > 0) {
        if (width(Seq) >= 31) {
          TransciptName = names(Transcripts[j])         # Get transcript name
          Instance = paste(InstanceGene, TransciptName, sep = "-")  # Define instance name considering GeneName and TranscriptName
          Length = lengthpep(Transcript)                 # Get sequcne length (1)
          Weigth = mw(Seq)                              # Get sequence weight (1)
          Zscales = zScales(Seq)[[1]]                   # Get Zscales (5)
          AAcomp = aaComp(Seq)                          # Get aminoacid composition (First librarie)
          
          # AMINOACID COMPOSITION
          AminoAcidComposition = extractAAC(Seq)        # Aminoacid compositon (20)
          DipeptideComposition = extractDC(Seq)         # Dipeptide composition (400)
          TripeptideComposition = extractTC(Seq)        # Tripeptide Composition (8000)

          
          # CTD (Aminoacids classified by attributes)
          Composition = extractCTDC(Seq)                # Composition (21)
          Transition = extractCTDT(Seq)                 # Transition (21)
          Distribution = extractCTDD(Seq)               # Distribution (105)
          
          # -------------- DATA COLLECTION ----------------
          
          TranscriptsList[[TransciptName]] = list('Name' = Instance,
                                                  'Seq' = Seq,                                                                       # Seq  Info
                                                  'Length' = Length,                                                                 # 1    Feature
                                                  'Weigth' = Weigth,                                                                 # 1    Feature
                                                  'Zscales' = Zscales,                                                               # 5    Features
                                                  'AminoAcidComposition' = AminoAcidComposition,                                     # 20   Features
                                                  'DipeptideComposition' = DipeptideComposition,                                     # 400  Features
                                                  'TripeptideComposition' = TripeptideComposition,                                   # 8000 Features
                                                  'Composition' = Composition,                                                       # 21   Features
                                                  'Transition' = Transition,                                                         # 21   Features                                 
                                                  'Distribution' = Distribution,                                                     # 105  Features
                                                  'Class' = Class)                                                                   # 1    Class (T/F)
        }
      }
    }
  }
  Proteins[[InstanceGene]] = TranscriptsList
}


######################################################################################################################################################
## CREATE DATAFRAME (BETTER RUN DATASETSPREPARATION) #################################################################################################
######################################################################################################################################################

ProteinsFrame = data.frame()
for (i in 1 : length(Proteins)){
  print(i)
  if (length(Proteins[[i]]) > 0 ){
    for (j in 1:length(Proteins[[i]])){
      CurrentProteinsFrame =  data.frame(
                              'Name' = Proteins[[i]][[j]]$Name,
                              'Length' = Proteins[[i]][[j]]$Length,
                              'Weigth' = Proteins[[i]][[j]]$Weigth,
                              t(Proteins[[i]][[j]]$Zscales),
                              t(Proteins[[i]][[j]]$AminoAcidComposition),
                              t(Proteins[[i]][[j]]$DipeptideComposition),
                              t(Proteins[[i]][[j]]$TripeptideComposition),
                              t(Proteins[[i]][[j]]$Composition),
                              t(Proteins[[i]][[j]]$Transition),
                              t(Proteins[[i]][[j]]$Distribution),
                              'Class' = Proteins[[i]][[j]]$Class)   
      ProteinsFrame = rbind(ProteinsFrame, CurrentProteinsFrame)
    }
  }
}


write.csv(ProteinsFrame, "/ML-based-prediction-of-CR-related-genes/Data/Datasets/ProteinsDataset.csv")

