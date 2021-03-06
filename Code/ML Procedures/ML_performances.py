################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################

# IMPORT LIBRARIES
import pickle
import math
import numpy as np 
import pandas as pd
import seaborn as sns
import scikitplot as skplt
from itertools import chain
import statsmodels.api as sm
import inspect, math, itertools
from xgboost import XGBClassifier
from sklearn import preprocessing
from varname import Wrapper, nameof
from matplotlib import pyplot as plt
from catboost import CatBoostClassifier
from sklearn.impute import SimpleImputer, KNNImputer
#import statsmodels.formula.api as sm
from imblearn.combine import SMOTEENN, SMOTETomek
from sklearn.ensemble import RandomForestClassifier
from imblearn.under_sampling import RandomUnderSampler
from plot_metric.functions import BinaryClassification
from sklearn.feature_selection import VarianceThreshold
from sklearn.preprocessing import OneHotEncoder, label_binarize, StandardScaler
from imblearn.metrics import geometric_mean_score, sensitivity_score, specificity_score
from sklearn.model_selection import GridSearchCV, train_test_split, StratifiedKFold, cross_validate
from imblearn.over_sampling import BorderlineSMOTE, SMOTENC, SMOTE, SVMSMOTE, ADASYN, RandomOverSampler
from imblearn.ensemble import BalancedRandomForestClassifier, EasyEnsembleClassifier, BalancedBaggingClassifier
from sklearn.metrics import roc_auc_score, roc_curve, auc, confusion_matrix, make_scorer, classification_report, plot_confusion_matrix

from MLfunctions import BRFparamGridFunc, EECparamGridFunc, CATparamGridFunc, XGBparamGridFunc, GeometricMean, FiPreProcessing, FstatisticFilter
from MLfunctions import PreProcessing, backwardElimination,  Sampling, unique, NanInfs, simple_imputer, Knn_imputer, delete_nan_values, get_binary_columns, unique, intersection

################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################

# --- DEFINE CLASSIFICATION DATASET AND PARAMETERS ----------------------------------------------------------------------------------------------------------

# DATASET
oDataset = pd.read_csv('/ML-based-prediction-of-CR-related-genes/Data/Datasets/PathDip.csv') # Or any of the datasets in /Data/Datasets 
nDataset = 'nPathDip'   # Or any of the datasets in /Data/Datasets 
nnDataset = nDataset
SplittedDataset = False
NanInfFilter = False
Imputation = 'None'
#SimpleImputerStrategy = 'mean'
ClassLabel = 'Class'
Target = ClassLabel

# CLASSIFIER 
Classifier = 'CAT'
if Classifier == 'EEC':
    MLalgorithm = EasyEnsembleClassifier(random_state=42)
    ParamGridFunc = EECparamGridFunc
elif Classifier == 'BRF':
    MLalgorithm = BalancedRandomForestClassifier(random_state=42)
    ParamGridFunc = BRFparamGridFunc
    nJobs = None
elif Classifier == 'CAT':
    MLalgorithm = CatBoostClassifier(random_state=42)
    ParamGridFunc = CATparamGridFunc
    nJobs = None
elif XGboost == 'XGB':
    MLalgorithm = XGBClassifier(random_state=42)
    ParamGridFunc = XGBparamGridFunc
    nJobs = None

# SCORE
MyScore = make_scorer(GeometricMean, greater_is_better=True)
Score = 'Gmean'

# SAMPLING
SamplingMethod = 'normal'

# PREPROCESSING
VUCPS = [True, True, True, False, True]   # For preprocessing, V: Variance-based Filtering, U: Univariative, C: Correlation-based Filtering, P: Pvalue-based Filtering, S: SCALING
bMask = [True, False, False, False, False]
cMask = [False, True, True, False, True]

MinocFlag = True
if MinocFlag:
    #Minoc = 4                      # Minimun allowed ocurrence in minority feature
    #Minoc = Minoc-1
    MinocArray = [3,4,5]
else:
    NZV = 0.995

CorrArray = [0.99]
#CorrCoef = 0.99
SL = 0.5

# CROSS VALIDATIONS
Method = 'cv'
OuterSplits = 10
CvOuterMethod = 'o'
InnerSplits = 5
CvInnerMethod = 'i'
ImportanceSplits = 5


# SEPARATOR
Cntr = '_'


# ------------------------------------------------------------------------------------------------------------------------------------------------------------

# DEFINE AND PREPARE DATASET
if oDataset.columns[0] == "Unnamed: 0":
    oDataset = oDataset.rename(columns={"Unnamed: 0": "Genes"})
else:
    oDataset = oDataset.rename(columns={"X": "Genes"})

oDataset.index = oDataset.iloc[:,0]            # Set the first column as row names 
oDataset = oDataset.drop(columns=['Genes'])    # Drop first column
oDataset
oDataset.shape

# FIND INF AND NAN INDICES
if NanInfFilter:
    NanIndexName, NanColumnName, InfIndexName, InfColumnName, NanInfColumnName, NanInfIndexName, NanInfLocNames = NanInfs(oDataset)
    Dataset = oDataset.drop(index=InfIndexName.values)    # Drop first column
else:
    Dataset = oDataset


# DEFINE FEATURES AND TARGET

Features = Dataset.columns.tolist()
del Features[-1]                             # Remove class from features
oLen = len(Features)

# FEATURES AND TARGET ARRAYS
X = Dataset[Features]


# Y LABEL
Ylabel = Dataset[Target].values

# LABEL ENCODING
le = preprocessing.LabelEncoder()
le.fit(Ylabel)
le.classes_
y = le.transform(Ylabel)


####################################################################################################################################################################
##### OUTER CROSS VALIDATION #######################################################################################################################################
####################################################################################################################################################################

##### OUTER SPLITTING ##############################################################################################################################################

# PARAMETERS SPLITTING
kXtrain = []
kYtrain = []
kXtest = [] 
kYtest = []
kYlabel = []
kYtrainCountCr = []
kYtestCountCr = []
kYtrainCountNotCr = []
kYtestCountNotCr = []

cont = 0

# OUTER CV SPLITTING FOR PROTEINS DATASET
if nDataset == 'ProteinsDataset':
    OuterSkf, OuterIntersect = TranscriptsSpliting(X, y,  Splits = OuterSplits, random_state = 42)
else:
   OuterSkf = StratifiedKFold(n_splits=OuterSplits, shuffle=True, random_state = 42)


for TrainIndex, TestIndex in OuterSkf.split(X, y):
    kXtrain.append(X.iloc[TrainIndex,:])
    kXtest.append(X.iloc[TestIndex]) 
    kYtrain.append(y[TrainIndex])
    kYtest.append(y[TestIndex])
    kYtrainCountCr.append(sum(kYtrain[cont] == 0))
    kYtestCountCr.append(sum(kYtest[cont] == 0))
    kYtrainCountNotCr.append(sum(kYtrain[cont] == 1))
    kYtestCountNotCr.append(sum(kYtest[cont] == 1))
    kYlabel.append(kXtrain[0].index)
    cont = cont+1

    print("DATASET SET:", len(y), "\n")
    print("MAJORITY CLASS IN DATASET:", sum(y == 1), "\n")
    print("MINORITY CLASS IN DATASET:", sum(y == 0), "\n\n")

    print("MAJORITY CLASS IN TRAINIG SET SPLITS:", kYtrainCountNotCr, "\n")
    print("MINORITY CLASS IN TRAINING SET SPLITS:", kYtrainCountCr, "\n")
    print("MAJORITY CLASS IN TESTING SET SPLITS:", kYtestCountNotCr, "\n")
    print("MINORITY CLASS IN TESTING SET SPLITS:", kYtestCountCr, "\n")


ModelsList = {}
k = 0
for k in range(0, OuterSplits):
    oXtrain = kXtrain[k]
    oYtrain = kYtrain[k]
    oXtest = kXtest[k]
    Ytest = kYtest[k]

    print("-------------------------------------------------- OUTER ", k, " ----------------------------------------------------------\n")

    ##### PREPROCESSING ################################################################################################################################################

    # RESAMPLING  
    sXtrain, Ytrain = Sampling(oXtrain, oYtrain, opc = SamplingMethod) # opc = {normal, under, over, smote, blsmote, csmote, svmsmote, smoteenn, smotetomek, adasyn}

    # --- IMPUTATION --------------------------------------------------------------------------------------------------------------------------------------------------------
    if Imputation == 'Simple':
        iXtrain = simple_imputer(sXtrain)
        iXtest = simple_imputer(oXtest)
    elif Imputation == 'Interative':
        ImpMean = IterativeImputer(random_state=0)
        X = ImpMean.fit_transform(oX)
        X = pd.DataFrame(X, index = oX.index.values, columns = oX.columns.values)
    elif Imputation == 'Knn':
        iXtrain = Knn_imputer(sXtrain)
        iXtest = Knn_imputer(oXtest)
    else:
        iXtrain = sXtrain
        iXtest = oXtest


    # --- PREPROCESSING  -----------------------------------------------------------------------------------------------------------------------------------------
    # DEFINE MINIMUN NUMBER OF FEATURE OCURRENCES FOR NEAR ZERO VARIANCE

    GmeanTrainMax = 0
    for CorrCoefK in CorrArray:
        for MinocK in MinocArray:
            MinocK = MinocK - 1
            print("-------------------------------------------------- OUTER ", k, ": ", MinocK, " ----------------------------------------------------------\n")

            Nrows = iXtrain.shape[0]
            Nrows
            if MinocFlag:
                MCP = 1 - (MinocK / Nrows)
                MCP # Put in MaxClassPercent to get a determined repetition of features across instances
                print('MinocFlag')
            if MinocFlag == False:
                MCP = NZV
                MinocK = (1-MCP) * Nrows
                print('NoMinocFlag')

            print('NZV = ' +  str(MCP))

            # SPLITTING BINARY AND CONTINUOUS
            biXtrain, ciXtrain = get_binary_columns(iXtrain)
            biXtest = iXtest[biXtrain.columns]
            ciXtest = iXtest[ciXtrain.columns]

            #biXtest, ciXtest = get_binary_columns(iXtest)

            b_oLenK = len(biXtrain.columns)
            c_oLenK = len(ciXtrain.columns)

            oLenK = b_oLenK + c_oLenK

            # FEATURE SELECTION (NEAR ZERO VARIANCE, CORRELATION, PVALUE-BASED ELIMINATION, AND SCALING-CENTERING)
            if SplittedDataset:
                print('SPLITTED')
                bXtrainK, bXtestK, b_vLenK, b_uvLenK, b_cuvLenK, b_pcuvLenK, b_VUCPtextK = PreProcessing(biXtrain, biXtest, Ytrain, MaxClassPercent = MCP, Ftop = Ftop, UnivariateTechnique = UnivariateTechnique, CorrCoef = CorrCoefK, SL = SL, VUCPS = np.multiply(VUCPS,bMask), Minoc = MinocK+1, Mds = 'B')
                print("Original Size     :", b_oLenK , "\n")
                print("Near Zero Var Size:", b_vLenK , "\n")
                print("Nzv + Univariate:", b_uvLenK , "\n")
                print("Nzv + Univariate + Corr Size   :", b_cuvLenK , "\n")

                cXtrainK, cXtestK, c_vLenK, c_uvLenK, c_cuvLenK, c_pcuvLenK, c_VUCPtextK = PreProcessing(ciXtrain, ciXtest, Ytrain, MaxClassPercent = MCP, Ftop = Ftop, UnivariateTechnique = UnivariateTechnique, CorrCoef = CorrCoefK, SL = SL, VUCPS = np.multiply(VUCPS,cMask), Minoc = MinocK+1, Mds = 'C')
                print("Original Size     :", c_oLenK , "\n")
                print("Near Zero Var Size + Univariate:", b_uvLenK , "\n")
                print("Nzv + Univariate + Corr Size   :", b_cuvLenK , "\n")

                XtrainK = pd.concat((bXtrainK, cXtrainK), axis=1)
                XtestK = pd.concat((bXtestK, cXtestK), axis=1)
                VUCPtextK =  b_VUCPtextK + c_VUCPtextK

                vLenK = b_vLenK + c_vLenK
                uvLenK = b_uvLenK + c_uvLenK
                cuvLenK = b_cuvLenK + c_cuvLenK
                pcuvLenK = b_pcuvLenK + c_pcuvLenK

                print("Original Size     :", oLenK , "\n")
                print("Near Zero Var Size:", vLenK , "\n")
                print("Nzv + Univariate   :", uvLenK , "\n")
                print("Nzv + Univariate + Corr Size   :", cuvLenK , "\n")

                #ysXtrain = Xtrain
                #ysXtest = Xtest

            
            if SplittedDataset ==  False:
                print('NON-SPLITTED')
                XtrainK, XtestK, vLenK, uvLenK, cuvLenK, pcuvLenK, VUCPtextK = PreProcessing(iXtrain, iXtest, Ytrain, MaxClassPercent = MCP, Ftop = Ftop, CorrCoef = CorrCoefK, SL = SL, VUCPS = VUCPS, Minoc = MinocK+1, Mds = 'M')
                b_vLenK = vLenK
                b_uvLenK = uvLenK
                b_cuvLenK = cuvLenK
                b_pcuvLenK = pcuvLenK
                c_vLenK = vLenK
                c_uvLenK = uvLenK
                c_cuvLenK = cuvLenK
                c_pcuvLenK = pcuvLenK

                #nsXtrain = Xtrain
                #nsXtest = Xtest

                print("Original Size     :", oLenK , "\n")
                print("Near Zero Var Size:", vLenK , "\n")
                print("Nzv + Univariate   :", uvLenK , "\n")
                print("Nzv + Univariate + Corr Size   :", cuvLenK , "\n")
            # ------------------------------------------------------------------------------------------------------------------------------------------------------------

            ##### INNER CV ###############################################################################################################################################

            # PARAMETERS SPLITTING
            kXlearn = []
            kXlearn = []
            kYlearn = []
            kXval = [] 
            kYval = []
            kYlearnCountCr = []
            kYvalCountCr = []
            kYlearnCountNotCr = []
            kYvalCountNotCr = []

            cont = 0
            for LearnIndex, ValIndex in InnerSkf.split(XtrainK, Ytrain):
                #print("TRAIN:", TrainIndex, "TEST:", TestIndex)
                #print("TEST:", TestIndex)
                kXlearn.append(XtrainK.iloc[LearnIndex,:])
                kXval.append(XtrainK.iloc[ValIndex]) 
                kYlearn.append(Ytrain[LearnIndex])
                kYval.append(Ytrain[ValIndex])
                kYlearnCountCr.append(sum(kYlearn[cont] == 0))
                kYvalCountCr.append(sum(kYval[cont] == 0))
                kYlearnCountNotCr.append(sum(kYlearn[cont] == 1))
                kYvalCountNotCr.append(sum(kYval[cont] == 1))
                cont = cont+1

            print("---------------------------------------- INNER ------------------------------------------------\n")
            print("TRAINING SET:", len(Ytrain), "\n")
            print("TESTING SET:", len(Ytest), "\n\n")

            print("MAJORITY CLASS IN TRAINING SET:", sum(Ytrain == 1), "\n")
            print("MINORITY CLASS IN TRAINING SET:", sum(Ytrain == 0), "\n\n")

            print("MAJORITY CLASS IN TESTING SET:", sum(Ytest == 1), "\n")
            print("MINORITY CLASS IN TESTING SET:", sum(Ytest == 0), "\n\n")

            print("MAJORITY CLASS IN LEARNING SET SPLITS:", kYlearnCountNotCr, "\n")
            print("MINORITY CLASS IN LEARNING SET SPLITS:", kYlearnCountCr, "\n")
            print("MAJORITY CLASS IN VALIDATION SET SPLITS:", kYvalCountNotCr, "\n")
            print("MINORITY CLASS IN VALIDATION SET SPLITS:", kYvalCountCr, "\n")

            # --- CREATE CLASSIFIER --------------------------------------------------------------------------------------------------------------------------------------
            # RANDOM FORESTS
    

            # ------------------------------------------------------------------------------------------------------------------------------------------------------------

            # INNER CROSS VALIDATION AND FIT
            ParamGridK = ParamGridFunc(XtrainK)
            Ylabel = XtestK.index

            # OUTER CV SPLITTING FOR PROTEINS DATASET
            if nDataset == 'ProteinsDataset':
                InnerSkf, InnerIntersect = TranscriptsSpliting(XtrainK, Ytrain,  Splits = InnerSplits, random_state = 42)
            else:
                InnerSkf = StratifiedKFold(n_splits=InnerSplits, shuffle=True, random_state = 42)

            print('INNER INTERSECT: ', InnerIntersect)

            GridSearchK = GridSearchCV(estimator = MLalgorithm, param_grid = ParamGridK, scoring = MyScore, cv = InnerSkf, verbose = 1, return_train_score = True)
            GridSearchK.fit(XtrainK, Ytrain)
            XtrainK.shape
            len(Ytrain)

            # BEST PARAMETERS
            BestGridK = GridSearchK.best_estimator_
            YpredTestK = BestGridK.predict(XtestK)
            YprobTestK = BestGridK.predict_proba(XtestK)[:, 1]
            YpredTrainK = BestGridK.predict(XtrainK)
            YprobTrainK = BestGridK.predict_proba(XtrainK)[:, 1]

            ##### PERFORMANCE ############################################################################################################################################

            # TRAINING METRICS
            ConfusionTrainK = confusion_matrix(YpredTrainK, Ytrain)     # It is reversed for easy of interpretation (it is just the transpose)
            print(ConfusionTrainK)

            SensitivityTrainK = sensitivity_score(Ytrain, YpredTrainK) 
            SpecificityTrainK = specificity_score(Ytrain, YpredTrainK) 
            GmeanTrainK = geometric_mean_score(Ytrain, YpredTrainK) 
            AucTrainK = roc_auc_score(Ytrain, YprobTrainK)

            if GmeanTrainK > GmeanTrainMax:
                GmeanTrainMax = GmeanTrainK
                BestGrid = BestGridK
                GridSearch = GridSearchK
                YpredTest = YpredTestK
                YprobTest = YprobTestK
                YpredTrain = YpredTrainK
                YprobTrain = YprobTrainK
                ConfusionTrain = ConfusionTrainK
                SensitivityTrain = SensitivityTrainK
                SpecificityTrain = SpecificityTrainK
                GmeanTrain = GmeanTrainK
                AucTrain = AucTrainK
                ParamGrid = ParamGridK

                b_oLen = b_oLenK
                c_oLen = c_oLenK
                oLen = oLenK
                b_vLen = b_vLenK
                b_cvLen = b_cvLenK
                b_pcvLen = b_pcvLenK
                c_vLen = c_vLenK
                c_cvLen = c_cvLenK
                c_pcvLen = c_pcvLenK
                VCPtext = VCPtextK
                Xtrain = XtrainK
                Xtest = XtestK
                vLen = vLenK
                cvLen = cvLenK
                pcvLen = pcvLenK
                Minoc = MinocK
                CorrCoef = CorrCoefK
        
            print('SensitivityTrain: ', SensitivityTrainK)
            print('SpecificityTrain: ', SpecificityTrainK)
            print('GmeanTrain: ', GmeanTrainK)
            print('AucTrain: ', AucTrainK)
        MinocArrayEnd= 1


    print("-------------------------------------------------- BEST MINOC: ",  Minoc, " ----------------------------------------------------------\n")
    
    print('SensitivityTrain: ', SensitivityTrain)
    print('SpecificityTrain: ', SpecificityTrain)
    print('GmeanTrain: ', GmeanTrain)
    print('AucTrain: ', AucTrain)

    # TESTING METRICS 
    ConfusionTest = confusion_matrix(YpredTest, Ytest)        # It is reversed for easy of interpretation (it is just the transpose)
    print(ConfusionTest)

    SensitivityTest = sensitivity_score(Ytest, YpredTest) 
    SpecificityTest = specificity_score(Ytest, YpredTest) 
    GmeanTest = geometric_mean_score(Ytest, YpredTest) 
    AucTest = roc_auc_score(Ytest, YprobTest)

    print('SensitivityTest: ', SensitivityTest )
    print('SpecificityTest: ', SpecificityTest)
    print('GmeanTest: ', GmeanTest)
    print('AucTest: ', AucTest)

    # LISTING RESULTS
    Metrics = {'ConfusionTrain': ConfusionTrain, 'SensitivityTrain': SensitivityTrain, 'SpecificityTrain': SpecificityTrain, 'GmeanTrain': GmeanTrain, 'AucTrain': AucTrain,\
                'ConfusionTest': ConfusionTest, 'SensitivityTest': SensitivityTest, 'SpecificityTest': SpecificityTest, 'GmeanTest': GmeanTest, 'AucTest': AucTest}
    Models = {'BestModel': BestGrid, 'Grid': GridSearch, 'MLalgorithm': MLalgorithm, 'ParamGrid': ParamGrid, 'Score':MyScore}
    Targets = {'Ytest': Ytest, 'YtestPred': YpredTest, 'YtestProb': YprobTest, 'Ylabel': Ylabel}
    Preprocessing = {'VUCPS': VUCPS, 'Mcp': MCP, 'Minoc': Minoc, 'MinocArray': MinocArray, 'CorrArray': CorrArray, 'CorrCoef': CorrCoef, 'SL': SL, 'oLen': oLen, 'vLen': vLen, 'cvLen': cvLen, 'pcvLen': pcvLen,\
                        'b_oLen': b_oLen, 'b_vLen': b_vLen, 'b_cvLen': b_cvLen, 'b_pcvLen': b_pcvLen,\
                        'c_oLen': c_oLen, 'c_vLen': c_vLen, 'c_cvLen': c_cvLen, 'c_pcvLen': c_pcvLen}
    
    CrossValidation = {'InnerSkf': InnerSkf, 'OuterSkf': OuterSkf}
    ModelsList[k] = {'Metrics':Metrics, 'Models': Models, 'Targets':Targets, 'Preprocessing':Preprocessing, 'CrossValidation':CrossValidation}

# AVERAGING RESULTS
CumulatedConfusionTrain = ModelsList[0]['Metrics']['ConfusionTrain'] - ModelsList[0]['Metrics']['ConfusionTrain']  #Just to initialize a zeros confusion matrix
CumulatedSensitivityTrain = 0
CumulatedSpecificityTrain = 0
CumulatedGmeanTrain = 0
CumulatedAucTrain = 0
CumulatedConfusionTest = ModelsList[0]['Metrics']['ConfusionTest'] - ModelsList[0]['Metrics']['ConfusionTest']  #Just to initialize a zeros confusion matrix
CumulatedSensitivityTest = 0
CumulatedSpecificityTest = 0
CumulatedGmeanTest = 0
CumulatedAucTest = 0

for k in range(0, OuterSplits):
    CumulatedConfusionTrain = CumulatedConfusionTrain + ModelsList[k]['Metrics']['ConfusionTrain']
    CumulatedSensitivityTrain = CumulatedSensitivityTrain + ModelsList[k]['Metrics']['SensitivityTrain']
    CumulatedSpecificityTrain = CumulatedSpecificityTrain + ModelsList[k]['Metrics']['SpecificityTrain']
    CumulatedGmeanTrain = CumulatedGmeanTrain + ModelsList[k]['Metrics']['GmeanTrain']
    CumulatedAucTrain = CumulatedAucTrain + ModelsList[k]['Metrics']['AucTrain']

    CumulatedConfusionTest = CumulatedConfusionTest + ModelsList[k]['Metrics']['ConfusionTest']
    CumulatedSensitivityTest = CumulatedSensitivityTest + ModelsList[k]['Metrics']['SensitivityTest']
    CumulatedSpecificityTest = CumulatedSpecificityTest + ModelsList[k]['Metrics']['SpecificityTest']
    CumulatedGmeanTest = CumulatedGmeanTest + ModelsList[k]['Metrics']['GmeanTest']
    CumulatedAucTest = CumulatedAucTest + ModelsList[k]['Metrics']['AucTest']

AverageConfusionTrain = CumulatedConfusionTrain / OuterSplits
AverageSensitivityTrain = CumulatedSensitivityTrain / OuterSplits
AverageSpecificityTrain = CumulatedSpecificityTrain / OuterSplits
AverageGmeanTrain = CumulatedGmeanTrain / OuterSplits
AverageAucTrain = CumulatedAucTrain / OuterSplits
AverageConfusionTest = CumulatedConfusionTest / OuterSplits
AverageSensitivityTest = CumulatedSensitivityTest / OuterSplits
AverageSpecificityTest = CumulatedSpecificityTest / OuterSplits
AverageGmeanTest = CumulatedGmeanTest / OuterSplits
AverageAucTest = CumulatedAucTest / OuterSplits

AveragePerformances = {'AverageConfusionTrain': AverageConfusionTrain, \
                        'AverageSensitivityTrain': AverageSensitivityTrain, \
                        'AverageSpecificityTrain': AverageSpecificityTrain, \
                        'AverageGmeanTrain': AverageGmeanTrain, \
                        'AverageAucTrain': AverageAucTrain, \
                        'AverageConfusionTest': AverageConfusionTest, \
                        'AverageSensitivityTest': AverageSensitivityTest, \
                        'AverageSpecificityTest': AverageSpecificityTest, \
                        'AverageGmeanTest': AverageGmeanTest, \
                        'AverageAucTest': AverageAucTest \
                        }



print(DatasetModel['AveragePerformances']['AverageConfusionTest'], '\n\n')
print(DatasetModel['AveragePerformances']['AverageSensitivityTest'], '\n\n')
print(DatasetModel['AveragePerformances']['AverageSpecificityTest'], '\n\n')
print(DatasetModel['AveragePerformances']['AverageGmeanTest'], '\n\n')
print(DatasetModel['AveragePerformances']['AverageAucTest'], '\n\n')



