library("ensembldb")
library("EnsDb.Hsapiens.v86")
library("Peptides")
library("stringr")
library("protr")
# FOR MORE INFO CHECK HERE
# https://rdrr.io/cran/Peptides/man/
# https://www.genome.jp/dbget-bin/www_bget?pf:RNA_pol_Rpb2_6
# file:///C:/Users/gusdany/Downloads/Selecting_different_protein_representations_and_cl.pdf

GeneInstancesList30 = readRDS("C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/Proteins/GeneInstancesList30.rds")

# ASSIGN THE HUMAN DATABASE TO A VARIABLE
edb <- EnsDb.Hsapiens.v86        

# EVALUATE WETHER IT EXISTS PROTEINT ANNOTATION AVAIALABLE
hasProteinData(edb)

# LOAD AGE AND CR GENES

AgeGenes = readRDS("AgeingFinal.rds")
CrGenes = readRDS("CrFinal.rds")

# LOAD AGEING GENES
InstanceGenes = AgeGenes

IntanceTranscriptsList = list()
IntanceTranscriptsPerGeneList = list()
GeneInstancesList30 = list()
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
          
          # AUTOCORRELATION
          NormalizedMoreauBrotoAutocorrelation = extractMoreauBroto(Seq)  # NormalizedMoreauBrotoAutocorrelation (240)
          MoranAutocorrelation = extractMoran(Seq)      # Moran Autocorrelation (240)
          GearyAutocorrelation = extractGeary(Seq)      # Geary Autocorrelation (240)
          
          # CTD (Aminoacids classified by attributes)
          Composition = extractCTDC(Seq)                # Composition (21)
          #extractCTDCClass()                           # Composition (21)
          Transition = extractCTDT(Seq)                 # Transition (21)
          # Transition = extractCTDTClass()             # Transition (21)
          Distribution = extractCTDD(Seq)               # Distribution (105)
          # Distribution = extractCTDDClass()           # Distribution (105)
          
          # CONJOINT TRIAD
          ConjointTriad = extractCTriad(Seq)            # ConjointTriad (343)
          # ConjointTriad = extractCTriadClass(Seq)     # ConjointTriad (343)
          
          # QUASI-SEQUENCE-ORDER
          SequenceOrderCouplingNumber = extractSOCN(Seq)   # Quasi-Sequence-Order number (60)
          QuasiSequenceOrderDescriptors = extractQSO(Seq)  # Quasi-Sequence-Order Descriptors (100)
          
          # PSEUDO AMINOACID COMPOSITION 
          PseudoAminoAcidComposition = extractPAAC(Seq)    # PseudoAminoAcidComposition (50)
          AmphiphilicPseudoAminoAcidComposition = extractAPAAC(Seq)  # Amphiphilic Pseudo-Amino Acid Composition (80)
          
          # -------------- DATA COLLECTION ----------------
          
          TranscriptsList[[TransciptName]] = list('Name' = Instance,
                                                  'Seq' = Seq,                                                                       # Seq  Info
                                                  'Length' = Length,                                                                 # 1    Feature
                                                  'Weigth' = Weigth,                                                                 # 1    Feature
                                                  'Zscales' = Zscales,                                                               # 5    Features
                                                  'AminoAcidComposition' = AminoAcidComposition,                                     # 20   Features
                                                  'DipeptideComposition' = DipeptideComposition,                                     # 400  Features
                                                  'TripeptideComposition' = TripeptideComposition,                                   # 8000 Features
                                                  'NormalizedMoreauBrotoAutocorrelation' = NormalizedMoreauBrotoAutocorrelation,     # 240  Features
                                                  'MoranAutocorrelation' = MoranAutocorrelation,                                     # 240  Features
                                                  'GearyAutocorrelation' = GearyAutocorrelation,                                     # 240  Features
                                                  'Composition' = Composition,                                                       # 21   Features
                                                  'Transition' = Transition,                                                         # 21   Features                                 
                                                  'Distribution' = Distribution,                                                     # 105  Features
                                                  'ConjointTriad' = ConjointTriad,                                                   # 343  Features
                                                  'SequenceOrderCouplingNumber' = SequenceOrderCouplingNumber,                       # 60   Features
                                                  'QuasiSequenceOrderDescriptors' = QuasiSequenceOrderDescriptors,                   # 100  Features
                                                  'PseudoAminoAcidComposition' = PseudoAminoAcidComposition,                         # 50   Features
                                                  'AmphiphilicPseudoAminoAcidComposition' = AmphiphilicPseudoAminoAcidComposition,   # 80   Features
                                                  'Class' = Class)                                                                   # 1    Class (T/F)
        }
      }
    }
  }
  GeneInstancesList30[[InstanceGene]] = TranscriptsList
  saveRDS(GeneInstancesList30, "C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/Proteins/GeneInstancesList30.rds")
}

saveRDS(GeneInstancesList30, "C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/Proteins/GeneInstancesList30.rds")

######################################################################################################################################################
## CREATE DATAFRAME (BETTER RUN DATASETSPREPARATION) #################################################################################################
######################################################################################################################################################

Proteins = GeneInstancesList30#readRDS("GeneInstancesList30.rds")

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
                              #t(Proteins[[i]][[j]]$NormalizedMoreauBrotoAutocorrelation),
                              #t(Proteins[[i]][[j]]$MoranAutocorrelation),
                              #t(Proteins[[i]][[j]]$GearyAutocorrelation),
                              t(Proteins[[i]][[j]]$Composition),
                              t(Proteins[[i]][[j]]$Transition),
                              t(Proteins[[i]][[j]]$Distribution),
                              #t(Proteins[[i]][[j]]$ConjointTriad),
                              #t(Proteins[[i]][[j]]$SequenceOrderCouplingNumber),
                              #t(Proteins[[i]][[j]]$QuasiSequenceOrderDescriptors),
                              #t(Proteins[[i]][[j]]$PseudoAminoAcidComposition),
                              #t(Proteins[[i]][[j]]$AmphiphilicPseudoAminoAcidComposition),
                              'Class' = Proteins[[i]][[j]]$Class)   
      ProteinsFrame = rbind(ProteinsFrame, CurrentProteinsFrame)
    }
  }
}

#saveRDS(ProteinsFrame, "ProteinsFrame.rds")
#ProteinsFrame = readRDS("ProteinsFrame.rds")

saveRDS(ProteinsFrame, "C:/Users/gusdany/Desktop/El Proyecto/Code/FinalCode/RDS/Proteins/ProteinsDataset.rds")
write.csv(ProteinsFrame, "C:/Users/gusdany/Desktop/El Proyecto/Databases/DataSetsHagr/Final/OMA/ProteinsDataset.csv")


column_to_rownames(ProteinsFrame, loc = 1)
