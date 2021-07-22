#### PATH DIP MODEL #########################################################################################################

Labels = read.csv("E:\\The Project\\NewPred\\y_pred_test\\labels_pathdip_cat.csv")
pProb = read.csv("E:\\The Project\\NewPred\\y_pred_test\\y_pred_pathdip_cat.csv")
Test = read.csv("E:\\The Project\\NewPred\\y_pred_test\\y_test_pathdip_cat.csv")

dim(Test)

Labels$X = NULL
Test$X = NULL
pProb$X = NULL

Prob = pProb[,c(2,4,6,8,10,12,14,16,18,20)]
Pred = as.data.frame(ifelse(Prob >= 0.5, 1, 0))

uPred = as.numeric(unlist(Pred))
uProb = as.numeric(unlist(Prob))
uTest = as.numeric(unlist(Test))
uLabels = as.character(unlist(Labels))

length(uPred)

AvailableIndices = !is.na(uPred)

aPred = uPred[AvailableIndices]
aProb = uProb[AvailableIndices]
aTest = uTest[AvailableIndices]
aLabels = uLabels[AvailableIndices]

length(aPred)

oDIP = data.frame('Label' = aLabels, 'Test' = aTest, 'Pred' = aPred, 'Prob' = aProb)
row.names(oDIP) = oDIP$Label


nDIP = oDIP[oDIP$Test == 0 & oDIP$Pred == 1,]

DIP = nDIP[order(nDIP$Prob,decreasing = TRUE),]

dim(nDIP)

row.names(DIP) = DIP$Label

DIP$Rank = seq(1, length(row.names(DIP)), by=1)
DIP


######################################################################################
# GO MODEL ###########################################################################
######################################################################################

oGO = read.csv('E:/The Project/BRF_GO.csv')
oGO = oGO[order(oGO$Prob),] # FOR BRF
oGO$Test = ifelse(oGO$Test == 1, 0 ,1 )
oGO$Pred = ifelse(oGO$Pred == 1, 0 ,1 )
oGO$Prob = 1 - oGO$Prob

oGO$X = NULL
oGO$Rank = 1:dim(oGO)[1]
row.names(oGO) = oGO$Label

bGO = oGO

oGO$Class = ifelse(oGO$Test == 0, 'NotCR', 'CR')


#########################################################################################################################
#### JOINT ANALYSIS #####################################################################################################
#########################################################################################################################
bGO = bGO[bGO$Test == 0,  ]


#fDIP = oDIP
fDIP = oDIP[oDIP$Test == 0,]
sDIP = fDIP[order(fDIP$Prob,decreasing = TRUE),]
sDIP$Rank = 1:dim(sDIP)[1]


row.names(oGO) = oGO$Label
oGO$X = NULL
#fGO = oGO
fGO = oGO[oGO$Test == 0,]
sGO = fGO[order(fGO$Prob,decreasing = TRUE),]
sGO$Rank = 1:dim(sGO)[1]

CommonLabels = intersect(sGO$Label, sDIP$Label)


cGO = sGO[CommonLabels,]
cDIP = sDIP[CommonLabels,]

colnames(cGO) = c("gLabel", "gTest",  "gPred",  "gProb",  "gRank", "gClass")
colnames(cDIP) = c("dLabel", "dTest",  "dPred",  "dProb",  "dRank")


cGO$gRank = 1:dim(cGO)[1]
cDIP$gRank = 1:dim(cDIP)[1] 


cDIP = cDIP[order(cDIP$Prob, decreasing = TRUE),]
cDIP$dRank = 1:length(cDIP$Rank)


cSET = cbind(cGO,cDIP[cGO$gLabel,])
cSET$MeanRank = (cSET$gRank + cSET$dRank) / 2
cSET$GeometricRank = (cSET$gRank * cSET$dRank)^0.5
cSET$GeometricScore = (cSET$gProb * cSET$dProb)^0.5
cSET$MeanScore = (cSET$gProb + cSET$dProb) / 2
gRank = cSET$gRank
dRank = cSET$dRank
cSET = cSET[order(cSET$MeanScore, decreasing = TRUE),]

ClassLabel =  ifelse(cSET$gTest >= 0.50, 'CR', 'NotCR')
MeanClass = ifelse(cSET$MeanScore >= 0.50, 1, 0)
MeanSet = data.frame('Label' = cSET$gLabel, 'Test' = cSET$gTest, 'Pred' = MeanClass, 'Prob' = cSET$MeanScore, 'Class' = ClassLabel)
GeomClass = ifelse(cSET$GeometricScore >= 0.50, 1, 0)
GeomSet = data.frame('Label' = cSET$gLabel, 'Test' = cSET$gTest, 'Pred' =GeomClass, 'Prob' = cSET$GeometricScore, 'Class' = ClassLabel)


data.frame('GS' = cSET$GeometricScore, 'MS' = cSET$MeanScore)
write.csv(MeanSet,"E:\\The Project\\MeanSetForGraphCommon.csv")
write.csv(GeomSet,"E:\\The Project\\GeomSetForGraph2.csv")

head(MeanSet)
head(GeomSet)

gANC = paste( paste( paste(as.character(round(100*gRank/872,1)), "% {", sep = ''), gRank , sep = ''), "th/872}",  sep = '')
dANC = paste( paste( paste(as.character(round(100*dRank/872,1)), "% {", sep = ''), dRank , sep = ''), "th/872}",  sep = '')

gFP = paste( paste( paste(as.character(round(100*gRank/273,1)), "% {", sep = ''), gRank , sep = ''), "th/273}",  sep = '')
dFP = paste( paste( paste(as.character(round(100*dRank/201,1)), "% {", sep = ''), dRank , sep = ''), "th/201}",  sep = '')

DF = data.frame('Label' = cSET$gLabel, 'mScore' = cSET$MeanScore, 'geScore' = cSET$GeometricScore, 'gRank' = gRank, 'dRank' = dRank,  'gANC' = gANC, 'dANC' = dANC, 'gFP' = gFP, 'dFP' = dFP, 'gCR score' = round(cSET$gProb,2), 'dCR score' = round(cSET$dProb,2))

write.csv(DF,"E:\\The Project\\FrameFinal_NotCR_Common.csv")


###################################################################################################
# FURTHER ANALYSIS ################################################################################
###################################################################################################

library(scales)

# CALLING DATASETS
Go_original = oGO
Dip_original = oDIP

dim(oDIP)

Go_original = oGO[oGO$Test == 0,]
Dip_original = oDIP[oDIP$Test == 0,]

CommonLabels = intersect(Go_original$Label, Dip_original$Label)

# s = single, j =joint / o = original, r = rescaled

# GO TERMS PARAMETERS
soGO = Go_original
soGO$sRank = as.integer(length(soGO$Prob) - rank(soGO$Prob) + 1)

srGO = soGO
srGO$Prob = rescale(soGO$Prob)
srGO$Pred = 1*(srGO$Prob >= 0.5)

jsrGO = srGO[CommonLabels,] # Rescale original distribution, compute classification and then select only the common terms

joGO = soGO[CommonLabels,]
joGO$jRank = as.integer(length(joGO$Prob) - rank(joGO$Prob) + 1) # Take common terms of original diustribution and rank

jrGO = joGO
jrGOProb = rescale(joGO$Prob)
jrGO$Pred = 1*(jrGO$Prob >= 0.5)           # Take common terms of original distribution, rescale, recompute classification


# DIP TERMS PARAMETERS

dim(joGO)


soDIP = Dip_original
soDIP$sRank = as.integer(length(soDIP$Prob) - rank(soDIP$Prob) + 1)

srDIP = soDIP
srDIP$Prob = rescale(soDIP$Prob)
srDIP$Pred = 1*(srDIP$Prob >= 0.5)

jsrDIP = srDIP[CommonLabels,] #

joDIP = soDIP[CommonLabels,]
joDIP$jRank = as.integer(length(joDIP$Prob) - rank(joDIP$Prob) + 1)

jrDIP = joDIP
jrDIP$Prob = rescale(joDIP$Prob)
jrDIP$Pred = 1*(jrDIP$Prob >= 0.5)


Test = jsrDIP$Test

# DIP RANKING
spRankDIP = round( 100 * (jrDIP$sRank / max(jrDIP$sRank)) , 2)
jpRankDIP =  round( 100 * (jrDIP$jRank / max(jrDIP$jRank)), 2)
# GO RANKING
spRankGO = round( 100 * (jrGO$sRank / max(jrGO$sRank)) , 2)
jpRankGO =  round( 100 * (jrGO$jRank / max(jrGO$jRank)), 2)

# ALL JOIN GENES FRAME CREATION

#####################################################################################
# ALL JOIN GENES FRAME CREATION
#####################################################################################

# ORIGINAL
AJGF = data.frame('Label' = CommonLabels, 'Test' = Test,
                  # DIP RANKS
                  'sRankDIP_986' = jrDIP$sRank, 'jRankDIP_981' = jrDIP$jRank,
                  'spRankDIP_986' = spRankDIP, 'jpRankDIP_981' = jpRankDIP,
                  # DIP PROB AND PRED
                  'joProbDIP' = joDIP$Prob, 'jrProbDIP' = jrDIP$Prob, 'jsrProbDIP' = jsrDIP$Prob,
                  'joPredDIP' = joDIP$Pred, 'jrPredDIP' = jrDIP$Pred, 'jsrPredDIP' = jsrDIP$Pred,
                  # GO RANKS
                  'sRankGO_1124' = jrGO$sRank, 'jRankGO_981' = jrGO$jRank,
                  'spRankGO_1124' = spRankGO, 'jpRankGO_981' = jpRankGO,
                  # GO PROB AND PRED
                  'joProbGO' = joGO$Prob, 'jrProbGO' = jrGO$Prob, 'jsrProbGO' = jsrGO$Prob,
                  'joPredGO' = joGO$Pred, 'jrPredGO' = jrGO$Pred, 'jsrPredGO' = jsrGO$Pred)

# ALL
AJGF = data.frame('Label' = CommonLabels, 'Test' = Test,
                  # DIP RANKS
                  'AsRankDIP_986' = jrDIP$sRank, 'AjRankDIP_981' = jrDIP$jRank,
                  'AspRankDIP_986' = spRankDIP, 'AjpRankDIP_981' = jpRankDIP,
                  # DIP PROB AND PRED
                  'AjoProbDIP' = joDIP$Prob, 'AjrProbDIP' = jrDIP$Prob, 'AjsrProbDIP' = jsrDIP$Prob,
                  'AjoPredDIP' = joDIP$Pred, 'AjrPredDIP' = jrDIP$Pred, 'AjsrPredDIP' = jsrDIP$Pred,
                  # GO RANKS
                  'AsRankGO_1124' = jrGO$sRank, 'AjRankGO_981' = jrGO$jRank,
                  'AspRankGO_1124' = spRankGO, 'AjpRankGO_981' = jpRankGO,
                  # GO PROB AND PRED
                  'AjoProbGO' = joGO$Prob, 'AjrProbGO' = jrGO$Prob, 'AjsrProbGO' = jsrGO$Prob,
                  'AjoPredGO' = joGO$Pred, 'AjrPredGO' = jrGO$Pred, 'AjsrPredGO' = jsrGO$Pred)


# NOT CR
NJGF = data.frame('Label' = CommonLabels, 'Test' = Test,
                  # DIP RANKS
                  'NsRankDIP_876' = jrDIP$sRank, 'NjRankDIP_872' = jrDIP$jRank,
                  'NspRankDIP_876' = spRankDIP, 'NjpRankDIP_872' = jpRankDIP,
                  # DIP PROB AND PRED
                  'NjoProbDIP' = joDIP$Prob, 'NjrProbDIP' = jrDIP$Prob, 'NjsrProbDIP' = jsrDIP$Prob,
                  'NjoPredDIP' = joDIP$Pred, 'NjrPredDIP' = jrDIP$Pred, 'NjsrPredDIP' = jsrDIP$Pred,
                  # GO RANKS
                  'NsRankGO_1010' = jrGO$sRank, 'NjRankGO_872' = jrGO$jRank,
                  'NspRankGO_1010' = spRankGO, 'NjpRankGO_872' = jpRankGO,
                  # GO PROB AND PRED
                  'NjoProbGO' = joGO$Prob, 'NjrProbGO' = jrGO$Prob, 'NjsrProbGO' = jsrGO$Prob,
                  'NjoPredGO' = joGO$Pred, 'NjrPredGO' = jrGO$Pred, 'NjsrPredGO' = jsrGO$Pred)


##################################################################################################
#### FURTHER ON FRAME
##################################################################################################

#### ORIGINAL ####################################################################################
# --- AVERAGEING RANKS ----------------------------------------
# Numeric Rank - Joint
AJGF$Am_jRank_981 = (AJGF$jRankDIP_981 + AJGF$jRankGO_981) / 2
AJGF$Gm_jRank_981 = sqrt(AJGF$jRankDIP_981 * AJGF$jRankGO_981)
# Numeric Rank - Single
AJGF$Am_sRank_986_1124 = (AJGF$sRankDIP_986 + AJGF$sRankGO_1124) / 2
AJGF$Gm_sRank_986_1124 = sqrt(AJGF$sRankDIP_986 * AJGF$sRankGO_1124)

# Probabilistic Rank
AJGF$Am_jpRank_981 = (AJGF$jpRankDIP_981 + AJGF$jpRankGO_981) / 2
AJGF$Gm_jpRank_981 = sqrt(AJGF$jpRankDIP_981 * AJGF$jpRankGO_981)
# Probabilistic Single
AJGF$Am_spRank_986_1124 = (AJGF$spRankDIP_986 + AJGF$spRankGO_1124) / 2
AJGF$Gm_spRank_986_1124 = sqrt(AJGF$spRankDIP_986 * AJGF$spRankGO_1124)

# --- AVERAGEING PROBS ----------------------------------------
# Joint Original
AJGF$Am_joProb = (AJGF$joProbDIP + AJGF$joProbGO) / 2
AJGF$Gm_joProb = sqrt(AJGF$joProbDIP * AJGF$joProbGO)          
# Joint Rescaled
AJGF$Am_jrProb = (AJGF$jrProbDIP + AJGF$jrProbGO) / 2
AJGF$Gm_jrProb = sqrt(AJGF$jrProbDIP * AJGF$jrProbGO)  
# Single-Rescaled - Joint
AJGF$Am_jsrProb = (AJGF$jsrProbDIP + AJGF$jsrProbGO) / 2
AJGF$Gm_jsrProb = sqrt(AJGF$jsrProbDIP * AJGF$jsrProbGO)  

AJGF[order(AJGF$Am_jsrProb, decreasing = TRUE),]

AJGF[order(AJGF$Gm_jsrProb, decreasing = TRUE),]

# ALL #########################################################################################


# Numeric Rank - Joint
AJGF$Aam_jRank_981 = (AJGF$AjRankDIP_981 + AJGF$AjRankGO_981) / 2
AJGF$Agm_jRank_981 = sqrt(AJGF$AjRankDIP_981 * AJGF$AjRankGO_981)
# Numeric Rank - Single
AJGF$Aam_sRank_986_1124 = (AJGF$AsRankDIP_986 + AJGF$AsRankGO_1124) / 2
AJGF$Agm_sRank_986_1124 = sqrt(AJGF$AsRankDIP_986 * AJGF$AsRankGO_1124)

# Probabilistic Rank
AJGF$Aam_jpRank_981 = (AJGF$AjpRankDIP_981 + AJGF$AjpRankGO_981) / 2
AJGF$Agm_jpRank_981 = sqrt(AJGF$AjpRankDIP_981 * AJGF$AjpRankGO_981)
# Probabilistic Single
AJGF$Aam_spRank_986_1124 = (AJGF$AspRankDIP_986 + AJGF$AspRankGO_1124) / 2
AJGF$Agm_spRank_986_1124 = sqrt(AJGF$AspRankDIP_986 * AJGF$AspRankGO_1124)

# --- AVERAGEING PROBS ----------------------------------------
# Joint Original
AJGF$Aam_joProb = (AJGF$AjoProbDIP + AJGF$AjoProbGO) / 2
AJGF$Agm_joProb = sqrt(AJGF$AjoProbDIP * AJGF$AjoProbGO)          
# Joint Rescaled
AJGF$Aam_jrProb = (AJGF$AjrProbDIP + AJGF$AjrProbGO) / 2
AJGF$Agm_jrProb = sqrt(AJGF$AjrProbDIP * AJGF$AjrProbGO)  
# Single-Rescaled - Joint
AJGF$Aam_jsrProb = (AJGF$AjsrProbDIP + AJGF$AjsrProbGO) / 2
AJGF$Agm_jsrProb = sqrt(AJGF$AjsrProbDIP * AJGF$AjsrProbGO)  



# NOT #########################################################################################


# Numeric Rank - Joint
NJGF$Nam_jRank_872 = (NJGF$NjRankDIP_872 + NJGF$NjRankGO_872) / 2
NJGF$Ngm_jRank_872 = sqrt(NJGF$NjRankDIP_872 * NJGF$NjRankGO_872)
# Numeric Rank - Single
NJGF$Nam_sRank_876_1010 = (NJGF$NsRankDIP_876 + NJGF$NsRankGO_1010) / 2
NJGF$Ngm_sRank_876_1010 = sqrt(NJGF$NsRankDIP_876 * NJGF$NsRankGO_1010)

# Probabilistic Rank - Joint
NJGF$Nam_jpRank_872 = (NJGF$NjpRankDIP_872 + NJGF$NjpRankGO_872) / 2
NJGF$Ngm_jpRank_872 = sqrt(NJGF$NjpRankDIP_872 * NJGF$NjpRankGO_872)
# Probabilistic Single - Single
NJGF$Nam_spRank_876_1010 = (NJGF$NspRankDIP_876 + NJGF$NspRankGO_1010) / 2
NJGF$Ngm_spRank_876_1010 = sqrt(NJGF$NspRankDIP_876 * NJGF$NspRankGO_1010)

# --- AVERAGEING PROBS ----------------------------------------
# Joint Original
NJGF$Nam_joProb = (NJGF$NjoProbDIP + NJGF$NjoProbGO) / 2
NJGF$Ngm_joProb = sqrt(NJGF$NjoProbDIP * NJGF$NjoProbGO)          
# Joint Rescaled
NJGF$Nam_jrProb = (NJGF$NjrProbDIP + NJGF$NjrProbGO) / 2
NJGF$Ngm_jrProb = sqrt(NJGF$NjrProbDIP * NJGF$NjrProbGO)  
# Single-Rescaled - Joint
NJGF$Nam_jsrProb = (NJGF$NjsrProbDIP + NJGF$NjsrProbGO) / 2
NJGF$Ngm_jsrProb = sqrt(NJGF$NjsrProbDIP * NJGF$NjsrProbGO)  




### SAVE ######################################################################################
write.csv(AJGF,'E:\\The Project\\AJGF_BRF.csv')
write.csv(NJGF,'E:\\The Project\\NJGF_BRF.csv')


dim(AJGF)
dim(NJGF)


row.names(AJGF) = AJGF$Label
row.names(NJGF) = NJGF$Label

FJGF = cbind(AJGF[NJGF$Label,], NJGF)

# FALSE POSITIVES
# All Dip
FJGF$AjoRankDIP_981_FP200 = FJGF$AjRankDIP_981 / sum(FJGF$AjoPredDIP)    # # Take common terms of original distribution and rank
FJGF$AjrRankDIP_981_FP223 = FJGF$AjRankDIP_981 / sum(FJGF$AjrPredDIP)    # # Take common terms of original distribution, rescale, recompute classification
FJGF$AjsrRankDIP_986_FP223 = FJGF$AsRankDIP_986 / sum(FJGF$AjsrPredDIP)   # # Rescale original distribution, compute classification and then select only the common terms
# N Dip
FJGF$NjoRankDIP_872_FP200 = FJGF$NjRankDIP_872 / sum(FJGF$NjoPredDIP)
FJGF$NjrRankDIP_872_FP270 = FJGF$NjRankDIP_872 / sum(FJGF$NjrPredDIP)
FJGF$NjsrRankDIP_876_FP270 = FJGF$NsRankDIP_876 / sum(FJGF$NjsrPredDIP)
# All GO
FJGF$AjoRankGO_981_FP309 = FJGF$AjRankGO_981 / sum(FJGF$AjoPredGO)
FJGF$AjrRankGO_981_FP309 = FJGF$AjRankGO_981 / sum(FJGF$AjrPredGO)
FJGF$AjsrRankGO_1124_FP195 = FJGF$AsRankGO_1124 / sum(FJGF$AjsrPredGO)
# N GO
FJGF$NjoRankGO_872_FP309 = FJGF$NjRankGO_872 / sum(FJGF$NjoPredGO)
FJGF$NjrRankGO_872_FP309 = FJGF$NjRankGO_872 / sum(FJGF$NjrPredGO)
FJGF$NjsrRankGO_1010_FP195 = FJGF$NsRankGO_1010 / sum(FJGF$NjsrPredGO)


write.csv(FJGF,'E:\\The Project\\FJGF_BRF.csv')

FJGF = AJGF[AJGF$Test==0,]
