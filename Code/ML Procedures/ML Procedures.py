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

################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################

# IMPORT LIBRARIES
                                                          

##### LOAD PRE-COMPUTED DATA (IF NECCESARY) ########################################################################################################################

file = open('C:/Users/gusdany/Desktop/El Proyecto/brfGoModelsList.py', 'r')
LoadedModelsList = pickle.load(file)

pickle_in = open('C:/Users/gusdany/Desktop/El Proyecto/brfGoModelsList.py',"rb")
example_dict = pickle.load(pickle_in)

##### READING THE DATASETS #########################################################################################################################################
GoDataset = pd.read_csv('C:/Users/gusdany/Desktop/El Proyecto/Databases/DataSetsHagr/GoDataSet.csv') 
sKeggPathDataset = pd.read_csv('C:/Users/gusdany/Desktop/El Proyecto/Databases/DataSetsHagr/sKeggPathDataset.csv') 
KeggGeneDataset = pd.read_csv('C:/Users/gusdany/Desktop/El Proyecto/Databases/DataSetsHagr/KeggGeneDataset.csv') 
ProteinsDataset = pd.read_csv('C:/Users/gusdany/Desktop/El Proyecto/Databases/DataSetsHagr/ProteinsDataset.csv') 
PPImeasuresDataset = pd.read_csv('C:/Users/gusdany/Desktop/El Proyecto/Databases/DataSetsHagr/PPImeasuresDataset.csv') 
PPIadjencyDataset = pd.read_csv('C:/Users/gusdany/Desktop/El Proyecto/Databases/DataSetsHagr/PPIadjencyDataset.csv') 
PPIfilteredAdjencyDataSet = pd.read_csv('C:/Users/gusdany/Desktop/El Proyecto/Databases/DataSetsHagr/PPIfilteredAdjencyDataSet.csv') 
PathDipDataset = pd.read_csv('C:/Users/gusdany/Desktop/El Proyecto/Databases/DataSetsHagr/PathDipDataSet.csv') 
GtexDataset = pd.read_csv('C:/Users/gusdany/Desktop/El Proyecto/Databases/DataSetsHagr/GtexDataset.csv') 
WholeDatasetImputation847 = pd.read_csv('C:/Users/gusdany/Desktop/El Proyecto/Databases/DataSetsHagr/WholeDatasetImputation847.csv') 
WholeDatasetImputation888 = pd.read_csv('C:/Users/gusdany/Desktop/El Proyecto/Databases/DataSetsHagr/WholeDatasetImputation888.csv') 
WholeDatasetImputation887 = pd.read_csv('C:/Users/gusdany/Desktop/El Proyecto/Databases/DataSetsHagr/WholeDatasetsImputation887.csv') 
WholeDatasetImputation884 = pd.read_csv('C:/Users/gusdany/Desktop/El Proyecto/Databases/DataSetsHagr/WholeDatasetsImputation884.csv') #
WholeDatasetOverlap460 = pd.read_csv('C:/Users/gusdany/Desktop/El Proyecto/Databases/DataSetsHagr/WholeDatasetOverlap460.csv') #

# FLORA DATASETS
MetabolonCocomo = pd.read_csv('C:/Users/gusdany/Downloads/Metabolon_terms_COCOMO.csv')
ExpressionCocomo = pd.read_csv('C:/Users/gusdany/Downloads/expression_data_cocomo.csv')

ExpressionCameroon = pd.read_csv('C:/Users/gusdany/Desktop/Metabolites/ML/NewCameroonExpressionDataset.csv')
KeggCameroon = pd.read_csv('C:/Users/gusdany/Desktop/Metabolites/ML/CameroonKeggFdrDataset.csv')
MetabolonCameroon = pd.read_csv('C:/Users/gusdany/Desktop/Metabolites/ML/CameroonMetabolonFdrDataset.csv')
ExpressionIndia = pd.read_csv('C:/Users/gusdany/Desktop/Metabolites/ML/NewIndiaExpressionDataset.csv')
KeggIndia = pd.read_csv('C:/Users/gusdany/Desktop/Metabolites/ML/IndiaKeggFdrDataset.csv')
MetabolonIndia = pd.read_csv('C:/Users/gusdany/Desktop/Metabolites/ML/IndiaMetabolonDataset.csv')

# --- DEFINE CLASSIFICATION DATASET AND PARAMETERS ----------------------------------------------------------------------------------------------------------
# PATHS
PathModel = 'C:/Users/gusdany/Desktop/El Proyecto/PythonModels/'
PathVarimp = 'C:/Users/gusdany/Desktop/El Proyecto/PythonVarimp/'

# DATASET
oDataset = pd.read_csv('E:/The Project/OMA/PathDip.csv') 
#del WholeDatasetOverlap460
nDataset = 'nPathDip'
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

# SCORE
MyScore = make_scorer(GeometricMean, greater_is_better=True)
Score = 'Gmean'

# SAMPLING
SamplingMethod = 'normal'

# PREPROCESSING
VCPS = [True, False, False, False]   # For preprocessing, V: Variance-based Filtering, C: Correlation-based Filtering, P: Pvalue-based Filtering, S: SCALING
bMask = [True, False, False, False]
cMask = [False, False, False, False]
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
OuterSkf = StratifiedKFold(n_splits=OuterSplits, shuffle=True, random_state = 42)
InnerSplits = 5
CvInnerMethod = 'i'
InnerSkf = StratifiedKFold(n_splits=InnerSplits, shuffle=True, random_state = 42)

ImportanceSplits = 5
ImportanceSkf = StratifiedKFold(n_splits=ImportanceSplits, shuffle=True, random_state = 42)

# SEPARATOR
Cntr = '_'


# ------------------------------------------------------------------------------------------------------------------------------------------------------------

#oDataset = oDataset.rename(columns={1: "Genes"})
#oDataset
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
# REMOVE NANS AND INFS

Dataset.shape
oDataset.shape

#Dataset.to_csv('C:/Users/gusdany/Desktop/El Proyecto/Databases/DataSetsHagr/WholeDatasetsOverlap454.csv')

#Dataset['CrRatio']

#.drop(columns=[InfColumn])    # Drop first column

# DEFINE FEATURES AND TARGET
#Features = list(Dataset.columns)
Features = Dataset.columns.tolist()
del Features[-1]                             # Remove class from features
oLen = len(Features)

# FEATURES AND TARGET ARRAYS
#XX = GoDataSet1[Features].values
X = Dataset[Features]

###########################
# HERE THAT CODE USED TO BE
###########################

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
                bXtrainK, bXtestK, b_vLenK, b_cvLenK, b_pcvLenK, b_VCPtextK = PreProcessing(biXtrain, biXtest, MaxClassPercent = MCP, CorrCoef = CorrCoefK, SL = SL, VCPS = np.multiply(VCPS,bMask), Minoc = MinocK+1, Mds = 'B')
                print("Original Size     :", b_oLenK , "\n")
                print("Near Zero Var Size:", b_vLenK , "\n")
                print("Nzv + Corr Size   :", b_cvLenK , "\n")

                cXtrainK, cXtestK, c_vLenK, c_cvLenK, c_pcvLenK, c_VCPtextK = PreProcessing(ciXtrain, ciXtest, MaxClassPercent = MCP, CorrCoef = CorrCoefK, SL = SL, VCPS = np.multiply(VCPS,cMask), Minoc = MinocK+1, Mds = 'C')
                print("Original Size     :", c_oLenK , "\n")
                print("Near Zero Var Size:", c_vLenK , "\n")
                print("Nzv + Corr Size   :", c_cvLenK , "\n")

                Xtrain = pd.concat((bXtrainK, cXtrainK), axis=1)
                Xtest = pd.concat((bXtestK, cXtestK), axis=1)
                VCPtextK =  b_VCPtextK + c_VCPtextK

                vLenK = b_vLenK + c_vLenK
                cvLenK = b_cvLenK + c_cvLenK
                pcvLenK = b_pcvLenK + c_pcvLenK

                print("Original Size     :", oLenK , "\n")
                print("Near Zero Var Size:", vLenK , "\n")
                print("Nzv + Corr Size   :", cvLenK , "\n")

                #ysXtrain = Xtrain
                #ysXtest = Xtest

            
            if SplittedDataset ==  False:
                print('NON-SPLITTED')
                XtrainK, XtestK, vLenK, cvLenK, pcvLenK, VCPtextK = PreProcessing(iXtrain, iXtest, MaxClassPercent = MCP, CorrCoef = CorrCoefK, SL = SL, VCPS = VCPS, Minoc = MinocK+1, Mds = 'M')
                b_vLenK = vLenK
                b_cvLenK = cvLenK
                b_pcvLenK = pcvLenK
                c_vLenK = vLenK
                c_cvLenK = cvLenK
                c_pcvLenK = pcvLenK

                #nsXtrain = Xtrain
                #nsXtest = Xtest

                print("Original Size     :", oLenK , "\n")
                print("Near Zero Var Size:", vLenK , "\n")
                print("Nzv + Corr Size   :", cvLenK , "\n")
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
    Preprocessing = {'VCPS': VCPS, 'Mcp': MCP, 'Minoc': Minoc, 'MinocArray': MinocArray, 'CorrArray': CorrArray, 'CorrCoef': CorrCoef, 'SL': SL, 'oLen': oLen, 'vLen': vLen, 'cvLen': cvLen, 'pcvLen': pcvLen,\
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

Preprocessing = {'VCPS': VCPS, 'Minoc': Minoc, 'CorrCoef': CorrCoef, 'SL': SL, 'MinocArray': MinocArray}
CrossValidation = {'OuterSkf':OuterSkf}

DatasetModel = {'AveragePerformances': AveragePerformances, 'ModelsList': ModelsList, 'Preprocessing':Preprocessing, 'CrossValidation':CrossValidation, 'Score':MyScore}

AvSensString = str(int(round(AverageSensitivityTest*100)))#str(round(AverageSensitivityTrain, 2))
AvSpecString = str(int(round(AverageSpecificityTest*100)))#str(round(AverageSpecificityTrain, 2))

# --- SAVING --------------------------------------------------------------------------------------------------------------------------------------------------------------------

SavingText = PathModel + \
             nDataset + Cntr + \
            Classifier + Cntr + \
             Score + Cntr + \
            SamplingMethod + Cntr + \
             VCPtext + Cntr + \
             Method + OuterSplits + CvOuterMethod + InnerSplits + Cntr + CvInnerMethod + Cntr + \
             AvSensString + AvSpecString + \
             + '.py'


file = open("E:\The Project\sGO.py", 'wb')
pickle.dump(DatasetModel, file)
file.close()                                                                    

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# TAKING A LOOK AT EACH OF THE TEST PERFORMANCES
show = 'Confusion'
for k in range(0,OuterSplits):
    if(show == 'Confusion'):
        print('Split ', k, ':\n', DatasetModel['ModelsList'][k]['Metrics']['ConfusionTest'], '\n\n')
    elif(show == 'Sensitivity'):
        print('Split ', k, ':\n', DatasetModel['ModelsList'][k]['Metrics']['SensitivityTest'], '\n\n')
    elif(show == 'Specificity'):
        print('Split ', k, ':\n', DatasetModel['ModelsList'][k]['Metrics']['SpecificityTest'], '\n\n')
    elif(show == 'Gmean'):
        print('Split ', k, ':\n', DatasetModel['ModelsList'][k]['Metrics']['GmeanTest'], '\n\n')
    elif(show == 'Auc'):
        print('Split ', k, ':\n', DatasetModel['ModelsList'][k]['Metrics']['AucTest'], '\n\n')


print(DatasetModel['AveragePerformances']['AverageConfusionTest'], '\n\n')
print(DatasetModel['AveragePerformances']['AverageSensitivityTest'], '\n\n')
print(DatasetModel['AveragePerformances']['AverageSpecificityTest'], '\n\n')
print(DatasetModel['AveragePerformances']['AverageGmeanTest'], '\n\n')
print(DatasetModel['AveragePerformances']['AverageAucTest'], '\n\n')


print('Split ', k, ':\n', DatasetModel['ModelsList'][0]['Metrics']['ConfusionTest'], '\n\n')
print('Split ', k, ':\n', DatasetModel['ModelsList'][0]['Metrics']['SensitivityTest'], '\n\n')
print('Split ', k, ':\n', DatasetModel['ModelsList'][0]['Metrics']['SpecificityTest'], '\n\n')



##############################################################################################################################################################
##### FEATURE IMPORTANCE #####################################################################################################################################
##############################################################################################################################################################


sXtrain, Ytrain = Sampling(X, y, opc = SamplingMethod) # opc = {normal, under, over, smote, blsmote, csmote, svmsmote, smoteenn, smotetomek, adasyn}


# --- IMPUTATION --------------------------------------------------------------------------------------------------------------------------------------------

# --- IMPUTATION --------------------------------------------------------------------------------------------------------------------------------------------
if Imputation == 'Simple':
    iXtrain = simple_imputer(sXtrain)
elif Imputation == 'Interative':
    ImpMean = IterativeImputer(random_state=0)
    X = ImpMean.fit_transform(oX)
    X = pd.DataFrame(X, index = oX.index.values, columns = oX.columns.values)
elif Imputation == 'Knn':
    iXtrain = Knn_imputer(sXtrain)
else:
    iXtrain = sXtrain

BestGrids = []

GmeanTrainMax = 0
for CorrCoefK in CorrArray:
    for MinocK in MinocArray:
        MinocK = MinocK - 1
        print("-------------------------------------------------- IMP: ", MinocK, " ----------------------------------------------------------\n")
        oLenImpK = oLen
        # FEATURE SELECTION (NEAR ZERO VARIANCE, CORRELATION, PVALUE-BASED ELIMINATION, AND SCALING-CENTERING)


        if MinocFlag:
            Nrows = sXtrain.shape[0]
            Nrows
            MCP = 1 - (MinocK / Nrows)
            MCP # Put in MaxClassPercent to get a determined repetition of features across instances
        else:
            MCP = NZV


        biXtrain, ciXtrain = get_binary_columns(iXtrain)
        b_oLenImpK = len(biXtrain.columns)
        c_oLenImpK = len(ciXtrain.columns)

        if SplittedDataset:
            # SPLITTING BINARY AND CONTINUOUS

            bXtrainK, b_vLenImpK, b_cvLenImpK, b_pcvLenImpK = FiPreProcessing(biXtrain, MaxClassPercent = MCP, CorrCoef = CorrCoefK, SL = SL, VCPS = np.multiply(VCPS,bMask), Minoc = MinocK, Mds = 'B')
            print("Original Size     :", b_oLenImpK , "\n")
            print("Near Zero Var Size:", b_vLenImpK , "\n")
            print("Nzv + Corr Size   :", b_cvLenImpK , "\n")

            cXtrainK, c_vLenImpK, c_cvLenImpK, c_pcvLenImpK = FiPreProcessing(ciXtrain, MaxClassPercent = MCP, CorrCoef = CorrCoefK, SL = SL, VCPS = np.multiply(VCPS,cMask), Minoc = MinocK, Mds = 'C')
            print("Original Size     :", c_oLenImpK , "\n")
            print("Near Zero Var Size:", c_vLenImpK , "\n")
            print("Nzv + Corr Size   :", c_cvLenImpK , "\n")

            XtrainK = pd.concat((bXtrainK, cXtrainK), axis=1)

            oLenImpK = b_oLenImpK + c_oLenImpK
            vLenImpK = b_vLenImpK + c_vLenImpK
            cvLenImpK = b_cvLenImpK + c_cvLenImpK
            pcvLenImpK = b_pcvLenImpK + c_pcvLenImpK

        else:
            Mds = 'Nmds'
            XtrainK, vLenImpK, cvLenImpK, pcvLenImpK, VCPtextK = FiPreProcessing(iXtrain, MaxClassPercent = MCP, CorrCoef = CorrCoefK, SL = SL, VCPS = VCPS, Minoc = MinocK, Mds = 'M')
            b_vLenImpK = vLenImpK
            b_cvLenImpK = cvLenImpK
            b_pcvLenImpK = pcvLenImpK
            c_vLenImpK = vLenImpK
            c_cvLenImpK = cvLenImpK
            c_pcvLenImpK = pcvLenImpK

            print("Original Size     :", oLenImpK , "\n")
            print("Near Zero Var Size:", vLenImpK , "\n")
            print("Nzv + Corr Size   :", cvLenImpK , "\n")



        ##### CLASSIFIER ##################################################################################################################################################

        # --- CREATE CLASSIFIER ---------------------------------------------------------------------------------------------------------------------------------------
        # RANDOM FORESTS
        #BRF = BalancedRandomForestClassifier(random_state=0)

        # -------------------------------------------------------------------------------------------------------------------------------------------------------------
        #custom_cv
        # INNER CROSS VALIDATION AND FIT
        ParamGridK = ParamGridFunc(XtrainK)
        GridSearchK = GridSearchCV(estimator = MLalgorithm, param_grid = ParamGridK, scoring = MyScore, cv = ImportanceSkf, verbose = 1, return_train_score = True)
        GridSearchK.fit(XtrainK, Ytrain)
        XtrainK.shape
        len(Ytrain)

        # BEST PARAMETERS
        BestGridK = GridSearchK.best_estimator_
        BestGridK.fit(XtrainK, Ytrain)

        BestGrids.append(BestGridK)

        # PREDICT TRAING
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
                    XtrainBest = XtrainK
                    GmeanTrainMax = GmeanTrainK
                    BestGrid = BestGridK
                    GridSearch = GridSearchK
                    YpredTrain = YpredTrainK
                    YprobTrain = YprobTrainK
                    ConfusionTrain = ConfusionTrainK
                    SensitivityTrain = SensitivityTrainK
                    SpecificityTrain = SpecificityTrainK
                    GmeanTrain = GmeanTrainK
                    AucTrain = AucTrainK
                    VCPtext = VCPtextK

                    b_oLenImp = b_oLenImpK
                    c_oLenImp = c_oLenImpK
                    oLenImp = oLenImpK
                    b_vLenImp = b_vLenImpK
                    b_cvLenImp = b_cvLenImpK
                    b_pcvLenImp = b_pcvLenImpK
                    c_vLenImp = c_vLenImpK
                    c_cvLenImp = c_cvLenImpK
                    c_pcvLenImp = c_pcvLenImpK
                    Xtrain = XtrainK
                    vLenImp = vLenImpK
                    cvLenImp = cvLenImpK
                    pcvLenImp = pcvLenImpK
                    Minoc = MinocK
                    CorrCoef = CorrCoefK

        print('SensitivityTrain: ', SensitivityTrainK )
        print('SpecificityTrain: ', SpecificityTrainK)
        print('GmeanTrain: ', GmeanTrainK)
        print('AucTrain: ', AucTrainK)
    MinocArrayEnd = 1

print("-------------------------------------------------- BEST MINOC: ",  Minoc, " ----------------------------------------------------------\n")

print('SensitivityTrain: ', SensitivityTrain )
print('SpecificityTrain: ', SpecificityTrain)
print('GmeanTrain: ', GmeanTrain)
print('AucTrain: ', AucTrain)


# --------------------------------------------------------

BestGridK.predict(Xtr)

Xtr = XtrainBest
Genes = Xtr['Genes']
Xtr = Xtr.drop(['Genes'], axis=1)

Xtr.shape
XtrainBest.shape
XtrainBest.columns
X.columns

XtrainBest.to_csv('C:/Users/gusdany/Desktop/El Proyecto/Databases/DataSetsHagr/Final/OMA/FilteredDipOldDataset.csv')

XtrainBest = pd.read_csv('C:/Users/gusdany/Desktop/El Proyecto/Databases/DataSetsHagr/Final/OMA/FilteredGoOldDataset.csv')


X = pd.read_csv('C:/Users/gusdany/Desktop/El Proyecto/Databases/DataSetsHagr/Final/OMA/fDipOldDatasetTemplate3.csv')
Genes = X['Unnamed: 0']
#X = X.drop(['Unnamed: 0', 'Class'], axis=1)
X = X.drop(['Unnamed: 0'], axis=1)

X.columns = Xtrain.columns
X.columns


BG = BestGrid
BG = BestGrids
BG.predict(X)
probs = BG.predict_proba(X)[:, 0]
Frame = pd.DataFrame({'Gene': Genes, 'Proba':probs})
Frame.loc[Frame['Proba'] > 0.5,]
Frame.to_csv('C:/Users/gusdany/Desktop/El Proyecto/Databases/DataSetsHagr/Final/OMA/CrDipFrame.csv')

GeneScoreFrame = pd.DataFrame({'Gene': Genes, 'Score':probs})
GeneScoreFrame.to_csv('C:/Users/gusdany/Desktop/El Proyecto/Databases/DataSetsHagr/Final/OMA/eGeneScoreDip3.csv')

BestModel.predict(XtrainBest)
BestModel.predict_proba(XtrainBest)




##### GET FEATURES IMPORTANCE #################################################################################################################################

# EXTRACT FEATURE IMPORTANCES
fi = pd.DataFrame({'Feature': list(Xtrain.columns),
                    'Score': BestGrid.feature_importances_}).\
                    sort_values('Score', ascending = False)
fi.index = fi.Feature

# DISPLAY
OriginalImportance = fi

fi.Score = fi.Score * 100 / (max(fi.Score))
fi.head()

ScaledImportance = fi

# LISTING RESULTS
Metrics = {'ConfusionTrain': ConfusionTrain, 'SensitivityTrain': SensitivityTrain, 'SpecificityTrain': SpecificityTrain, 'GmeanTrain': GmeanTrain, 'AucTrain': AucTrain}
Models = {'BestModel': BestGrid, 'Grid': GridSearch, 'MLalgorithm': MLalgorithm, 'ParamGrid': ParamGrid, 'Score':MyScore}
Targets = {'Ytrain': Ytest, 'YtrainPred': YpredTrain, 'YtrainProb': YprobTrain}
Preprocessing = {'VCPS': VCPS, 'Minoc': Minoc, 'MinocArray': MinocArray, 'CorrArray': CorrArray, 'CorrCoef': CorrCoef, 'SL': SL, 'oLenImp': oLenImp, 'vLenImp': vLenImp, 'cvLenImp': cvLenImp, 'pcvLenImp': pcvLenImp,\
                    'b_oLenImp': b_oLenImp, 'b_vLenImp': b_vLenImp, 'b_cvLenImp': b_cvLenImp, 'b_pcvLenImp': b_pcvLenImp,\
                    'c_oLenImp': c_oLenImp, 'c_vLenImp': c_vLenImp, 'c_cvLenImp': c_cvLenImp, 'c_pcvLenImp': c_pcvLenImp}
CrossValidation = {'ImportanceSkf': ImportanceSkf}
Importance= {'OriginalImportance':OriginalImportance, 'ScaledImportance':ScaledImportance}
ImportanceModel = {'Metrics':Metrics, 'Models': Models, 'Targets':Targets, 'Preprocessing':Preprocessing, 'CrossValidation':CrossValidation, 'Importance':Importance}

Arrays = 'ma'
for Minoc in MinocArray:
    Arrays = Arrays + '-' + str(Minoc)
Arrays = Arrays + '-oa'
for CorrCoef in CorrArray:
    Arrays = Arrays + '-' + str(CorrCoef)

# CRATING AND SAVING FINAL MODEL
FinalModel = {'DatasetModel':DatasetModel, 'ImportanceModel':ImportanceModel}
ModelSavingText = PathModel + \
                    nDataset + Cntr + \
                    Classifier + Cntr + \
                    Score + Cntr + \
                    SamplingMethod + Cntr + \
                    Arrays + Cntr + \
                    VCPtext + Cntr + \
                    Method + str(OuterSplits) + CvOuterMethod + str(InnerSplits) + CvInnerMethod + Cntr + \
                    'sn' + AvSensString + 'sp' + AvSpecString + \
                    '.py'
ModelSavingText
file = open(ModelSavingText, 'wb')
pickle.dump(FinalModel, file)
file.close()
LoadedModelsList = pickle.load(file)

#FinalModel

# SAVING IMPORTANCES FRAME
ImportanceSavingText = PathVarimp + \
                        nDataset + Cntr + \
                        Classifier + Cntr + \
                        Score + Cntr + \
                        SamplingMethod + Cntr + \
                        VCPtext + Cntr + \
                        Method + str(OuterSplits) + CvOuterMethod + str(InnerSplits) + CvInnerMethod + Cntr + \
                        'sn' + AvSensString + 'sp' + AvSpecString + \
                        '_vi.csv'
ImportanceSavingText 
ScaledImportance.to_csv(ImportanceSavingText) 

print(DatasetModel['AveragePerformances']['AverageConfusionTest'], '\n\n')
print(DatasetModel['AveragePerformances']['AverageSensitivityTest'], '\n\n')
print(DatasetModel['AveragePerformances']['AverageSpecificityTest'], '\n\n')
print(DatasetModel['AveragePerformances']['AverageGmeanTest'], '\n\n')
print(DatasetModel['AveragePerformances']['AverageAucTest'], '\n\n')
print(ScaledImportance)



###############################################################################################################################################################
###############################################################################################################################################################
##### FUNCTIONS ###############################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################


##### GRID ####################################################################################################################################################

def BRFparamGridFunc(Xtrain):
    MaxFeaturesGridNum = 5
    PowDivisions = 1/MaxFeaturesGridNum
    MaxFeatures = []
    NumFeatures = Xtrain.shape[1] 
    for i in range(1, MaxFeaturesGridNum + 1):
        MaxFeatures.append(round(pow(NumFeatures,PowDivisions*i)))
    MaxFeatures
    BRFparamGrid = {
        'bootstrap': [True],
        'replacement': [True, False],
        'max_features': ['sqrt','log2'], 
        'n_estimators': [500],
        'class_weight': ['balanced', 'balanced_subsample', None],
        'sampling_strategy': [1],
        }
    return(BRFparamGrid)

ParamGridVlad = {
    'min_samples_split': [3, 5, 10, 20],
    'n_estimators': [100, 500],
    'max_depth': [1, 5, 20],
    'max_features': [1, 5, 10, 30]
}

def EECparamGridFunc(Xtrain):
     EECparamGrid = {
        'replacement': [True, False],
        'n_estimators': [500],
        'sampling_strategy': [1]
     }
     return EECparamGrid


def CATparamGridFunc(Xtrain):
     CATparamGrid = {
        'n_estimators': [500],
     }
     return CATparamGrid



def BBCparamGridFunc(Xtrain):
     BBCparamGrid = {
        'sampling_strategy':1,
        'n_estimators':5
     }
     return BBCparamGrid


def GridCreator(Min, Max, Num, Func):
    Num = Num - 1
    Grid = []
    if func == 'lin':
        for i in range(0,Num+1):
            Slope = (Max - Min)/Num
            iValue = Min + Slope*i
            Grid.append(round(iValue))
    elif func == 'exp':
        for i in range(0,Num+1):
            ExpTerm = pow(Max-Min,i)
            iValue = rpund(Min + ExpTerm - 1 + i)
            Grid.append(round(iValue))
    elif func == 'iexp':
        for i in range(0,Num+1):
            ExpTerm = pow(Max-Min,i)
            iValue = Max - Min - rpund(Min + ExpTerm - 1 + i)
            Grid.append(round(iValue))
    elif func == 'gauss':  # FAKE GAUSS PROJECTION MADE WITH POLYNOMIAL
            Xmidl = 0.5
            Xup = 0.8
            Yup = 0.6
            A = np.array([[0, 0, 0, 1], [1, 1, 1, 1], [pow(Xmidl,3), pow(Xmidl,2), Xmidl, 1], [pow(Xup,3), pow(Xup,2), Xup, 1]])
            B = np.array([Min, Max, (Min+Max)/2, (Yup*Min)+((1-Yup)*max)])
            X = np.linalg.inv(A).dot(B)
            iValue = X[0]*pow(i,3) + X[1]*pow(i,2) + X[2]*i + X[3]
            Grid.append(round(iValue))
    return Grid
##### SCORE ###################################################################################################################################################

def GeometricMean(Ytest,Ypred):
   GM = geometric_mean_score(Ytest, Ypred)
   return GM

##### CONFUSION MATRIX ########################################################################################################################################

def plot_confusion_matrix(cm, classes, normalize=False, title='Confusion matrix', cmap=plt.cm.Oranges): 

    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
        print("Normalized confusion matrix")
    else:
        print('Confusion matrix, without normalization')

    print(cm)

    # Plot the confusion matrix
    plt.figure(figsize = (10, 10))
    plt.imshow(cm, interpolation='nearest', cmap=cmap)
    plt.title(title, size = 24)
    plt.colorbar(aspect=4)
    tick_marks = np.arange(len(classes))
    plt.xticks(tick_marks, classes, rotation=45, size = 14)
    plt.yticks(tick_marks, classes, size = 14)

    fmt = '.2f' if normalize else 'd'
    thresh = cm.max() / 2.
    
    # Labeling the plot
    for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
        plt.text(j, i, format(cm[i, j], fmt), fontsize = 20,
                 horizontalalignment="center",
                 color="white" if cm[i, j] > thresh else "black")
        
    plt.grid(None)
    plt.tight_layout()
    plt.ylabel('True label', size = 18)
    plt.xlabel('Predicted label', size = 18)

##### PRE PROCESSING ################################################################################################################################################

def FiPreProcessing(sXtrain, MaxClassPercent = 1, CorrCoef = 1, SL = 0.5, VCPS = [True, False, False, True], Minoc = 1, Mds = 'M'):

    Mtext = Mds

    # --- INDICES AND COLUMNS ---------------------------------------------------------------------------------------------------------------------------------------
    Features = sXtrain.columns
    TrainIndices = sXtrain.index

    Vtext = ''
    Ctext = ''
    Ptext = ''
    Stext = ''

    # --- FEATURES SELECTION (NEAR ZERO VARIANCE) -------------------------------------------------------------------------------------------------------------------
    if (VCPS[0]):
        Vtext = 'm' + str(math.floor(Minoc))# + 'v' + str(round(MaxClassPercent*100))
        Variance = MaxClassPercent * (1 - MaxClassPercent)
        selector = VarianceThreshold(Variance)
        selector.fit_transform(sXtrain)
        vSelectedIndices = selector.get_support(indices=True)
        vSelectedColumns = Features[vSelectedIndices]
    else:
        vSelectedColumns= Features
    vXtrain = sXtrain[vSelectedColumns]
    vLen = len(vSelectedColumns)

    # --- FEATURES SELECTION (CORRELATION) ---------------------------------------------------------------------------------------------------------------------------
    # CORRELATION SELECTION
    if(VCPS[1]):
        Ctext = 'c' + str(round(CorrCoef*100))
        corr = vXtrain.corr()
        #sns.heatmap(corr)
        columns = np.full((corr.shape[0],), True, dtype=bool)
        for i in range(corr.shape[0]):
            for j in range(i+1, corr.shape[0]):
                if corr.iloc[i,j] >= CorrCoef:
                    if columns[j]:
                        columns[j] = False
        cvSelectedColumns = vXtrain.columns[columns]
    else:
        cvSelectedColumns = vSelectedColumns
    cvXtrain = vXtrain[cvSelectedColumns]
    cvLen = len(cvSelectedColumns)
    
    # PVALUES SELECTION (BACKWARD ELIMINATION)
    if (VCPS[2]):
        Ptext = 'p' + str(round(SL))
        pcvSelectedColumns = backwardElimination(cvXtrain.values, Ytrain, SL, cvSelectedColumns)
    else:
        pcvSelectedColumns = cvSelectedColumns
    pcvXtrain = cvXtrain[pcvSelectedColumns]
    pcvLen = len(pcvSelectedColumns)

    # --- SCALER -----------------------------------------------------------------------------------------------------------------------------------------------------
    if (VCPS[3]):
        Stext = 's' 
        sc = StandardScaler()
        spcvXtrain = sc.fit_transform(pcvXtrain)
        #spcvXtest = sc.transform(pcvXtest)
    else:
        spcvXtrain = pcvXtrain
        #spcvXtest = pcvXtest

    # --- DATA FRAMES ------------------------------------------------------------------------------------------------------------------------------------------------
    Xtrain = pd.DataFrame(data=spcvXtrain, columns = pcvSelectedColumns, index = TrainIndices)
    VCPtext = Mtext + Vtext + Ctext + Ptext + Stext

    return Xtrain, vLen, cvLen, pcvLen, VCPtext




def PreProcessing(sXtrain, oXtest, MaxClassPercent = 1, CorrCoef = 1, SL = 0.5, VCPS = [True, False, False, True], Minoc = 1,  Mds = 'M'):

    Mtext = Mds

    # --- INDICES AND COLUMNS ---------------------------------------------------------------------------------------------------------------------------------------
    Features = sXtrain.columns
    TrainIndices = sXtrain.index
    TestIndices = oXtest.index
    Vtext = ''
    Ctext = ''
    Ptext = ''
    Stext = ''

    # --- FEATURES SELECTION (NEAR ZERO VARIANCE) -------------------------------------------------------------------------------------------------------------------
    if (VCPS[0]):
        Vtext = 'm' + str(math.floor(Minoc))# + 'v' + str(round(MaxClassPercent*100))
        Variance = MaxClassPercent * (1 - MaxClassPercent)
        selector = VarianceThreshold(Variance)
        selector.fit_transform(sXtrain) # Originally oXtrain
        vSelectedIndices = selector.get_support(indices=True)
        vSelectedColumns = Features[vSelectedIndices]
    else:
        vSelectedColumns= Features
    vXtrain = sXtrain[vSelectedColumns]
    vLen = len(vSelectedColumns)

    # --- FEATURES SELECTION (CORRELATION) ---------------------------------------------------------------------------------------------------------------------------
    # CORRELATION SELECTION
    if(VCPS[1]):
        Ctext = 'c' + str(round(CorrCoef*100))
        corr = vXtrain.corr()
        #sns.heatmap(corr)
        columns = np.full((corr.shape[0],), True, dtype=bool)
        for i in range(corr.shape[0]):
            for j in range(i+1, corr.shape[0]):
                if corr.iloc[i,j] >= CorrCoef:
                    if columns[j]:
                        columns[j] = False
        cvSelectedColumns = vXtrain.columns[columns]
    else:
        cvSelectedColumns = vSelectedColumns
    cvXtrain = vXtrain[cvSelectedColumns]
    cvLen = len(cvSelectedColumns)
    
    # --- PVALUES SELECTION (BACKWARD ELIMINATION) -------------------------------------------------------------------------------------------------------------------
    if (VCPS[2]):   
        Ptext = 'p' + str(round(SL))
        pcvSelectedColumns = backwardElimination(cvXtrain.values, Ytrain, SL, cvSelectedColumns)
    else:
        pcvSelectedColumns = cvSelectedColumns
    pcvXtrain = cvXtrain[pcvSelectedColumns]
    pcvLen = len(pcvSelectedColumns)

    # --- XTEST ------------------------------------------------------------------------------------------------------------------------------------------------------
    pcvXtest = oXtest[pcvSelectedColumns]

    # --- SCALER -----------------------------------------------------------------------------------------------------------------------------------------------------

    if (VCPS[3]):
        Stext = 's' 
        sc = StandardScaler()
        spcvXtrain = sc.fit_transform(pcvXtrain)
        spcvXtest = sc.transform(pcvXtest)
    else:
        spcvXtrain = pcvXtrain
        spcvXtest = pcvXtest

    # --- DATA FRAMES ------------------------------------------------------------------------------------------------------------------------------------------------
    Xtrain = pd.DataFrame(data=spcvXtrain, columns = pcvSelectedColumns, index = TrainIndices)
    Xtest = pd.DataFrame(data=spcvXtest, columns = pcvSelectedColumns, index = TestIndices)
    VCPtext = Mtext + Vtext + Ctext + Ptext + Stext

    return Xtrain, Xtest, vLen, cvLen, pcvLen, VCPtext



# BACKWARD ELIMINATION FUCTION
def backwardElimination(x, Y, sl, columns):
    numVars = len(x[0])
    for i in range(0, numVars):
        regressor_OLS = sm.OLS(Y, x).fit()
        maxVar = max(regressor_OLS.pvalues).astype(float)
        if maxVar > sl:
            for j in range(0, numVars - i):
                if (regressor_OLS.pvalues[j].astype(float) == maxVar):
                    x = np.delete(x, j, 1)
                    columns = np.delete(columns, j)
                    
    regressor_OLS.summary()
    return columns


##### SAMPLING #####################################################################################################################################################

def Sampling(oXtrain, oYtrain, opc = "normal"):

    # NORMAL
    if (opc == "normal"):
        sXtrain = oXtrain
        Ytrain = oYtrain

    # RANDOM UNDER SAMPLING
    elif(opc == "under"):
        rus = RandomUnderSampler(random_state=42)
        sXtrain, Ytrain = rus.fit_resample(oXtrain, oYtrain)  

    # RANDOM OVER SAMPLER
    elif(opc == "over"):
        ros = RandomOverSampler(random_state=42)
        sXtrain, Ytrain = ros.fit_resample(oXtrain, oYtrain)               

    # SMOTE
    elif(opc == "smote"):
        sm = SMOTE(random_state=42)
        sXtrain, Ytrain = sm.fit_resample(oXtrain, oYtrain)

    # BORDER LINE SMOTE
    elif(opc == "blsmote"):
         blm = BorderlineSMOTE(random_state=42)                            
         sXtrain, Ytrain = blm.fit_resample(oXtrain, oYtrain)

    #CATEGORIAL SMOTE
    elif(opc == "csmote"):
        smote_nc = SMOTENC(categorical_features=[range(0,len(Features))], random_state=42)
        sXtrain, Ytrain = smote_nc.fit_resample(oXtrain, oYtrain)

    # SVM SMOTE
    elif(opc == "svmsmote"):
        sm = SVMSMOTE(random_state=42)
        sXtrain, Ytrain = sm.fit_resample(oXtrain, oYtrain)

    # SMOTE ENN
    elif(opc == "smoteenn"):
        sme = SMOTEENN(random_state=42)
        uXtrain, Ytrain = sme.fit_resample(oXtrain, oYtrain)

    # SMOTE TOMEK
    elif(opc == "smotetomek"):
        smt = SMOTETomek(random_state=42)
        sXtrain, Ytrain = smt.fit_resample(oXtrain, oYtrain)

    # ADASYN
    elif(opc == "adasyn"):
        ada = ADASYN(random_state=42)
        sXtrain, Ytrain = ada.fit_resample(oXtrain, oYtrain)

    else:
        sXtrain = oXtrain
        Ytrain = oYtrain

    return sXtrain, Ytrain

##### UNIQUE ######################################################################################################################################################


def unique(list1): 
    list_set = set(list1) 
    unique_list = (list(list_set)) 
    return(unique_list)

##### NANINFS #####################################################################################################################################################

def NanInfs(Dataset):
    Instances = Dataset.index
    Features = Dataset.columns
    Nrows = Dataset.shape[0]
    Ncols = Dataset.shape[1]
    NanLoc = []
    InfLoc = []
    NanLocNames = []
    InfLocNames = []
    for row in range(0,Nrows):
        for col in range(0,Ncols-1):
            if math.isnan(Dataset.iloc[row,col]):
                NanLoc.append([row,col])
                NanLocNames.append([Instances[row],Features[col]])
            if np.isinf(Dataset.iloc[row,col]):
                InfLoc.append([row,col])
                InfLocNames.append([Instances[row],Features[col]])
    Nnan = len(NanLoc)
    NanIndex = []
    NanColumn = []
    for i in range(0,Nnan):
        NanIndex.append(NanLoc[i][0])
        NanColumn.append(NanLoc[i][1])
    NanIndex = unique(NanIndex)
    NanColumn = unique(NanColumn)
    NanIndexName = Instances[NanIndex]
    NanColumnName = Features[NanColumn]

    Ninf = len(InfLoc)
    InfIndex = []
    InfColumn = []
    for i in range(0,Ninf):
        InfIndex.append(InfLoc[i][0])
        InfColumn.append(InfLoc[i][1])
    InfIndex = unique(InfIndex)
    InfColumn = unique(InfColumn)
    InfIndexName = Instances[InfIndex]
    InfColumnName = Features[InfColumn]

    NanInfColumn = unique(NanColumn + InfColumn)
    NanInfIndex = unique(NanIndex + InfIndex)

    NanInfColumnName = Features[NanInfColumn]
    NanInfIndexName = Instances[NanInfIndex]
    
    NanInfLoc = NanLoc + InfLoc
    NanInfLocNames = NanLocNames + InfLocNames

    return NanIndexName, NanColumnName, InfIndexName, InfColumnName, NanInfColumnName, NanInfIndexName, NanInfLocNames

##### IMPUTATION #################################################################################################################################################

def simple_imputer(features_array, strategy='mean'):
    #features_array = oX
    binary_columns_gaps, non_binary_columns_gaps = get_binary_columns(features_array)
    imputer_binary = SimpleImputer(strategy='most_frequent')#SimpleImputer(strategy='constant', fill_value=0)#SimpleImputer(strategy=strategy)
    imputer_binary.fit(binary_columns_gaps)
    binary_columns_imputed = imputer_binary.transform(binary_columns_gaps)
    binary_columns_imputed = pd.DataFrame(binary_columns_imputed, index = binary_columns_gaps.index, columns = binary_columns_gaps.columns)

    imputer = SimpleImputer(strategy=strategy)#SimpleImputer(strategy='constant', fill_value=0)
    imputer.fit(non_binary_columns_gaps)
    non_binary_columns_imputed = imputer.transform(non_binary_columns_gaps)
    non_binary_columns_imputed = pd.DataFrame(non_binary_columns_imputed, index = non_binary_columns_gaps.index, columns = non_binary_columns_gaps.columns)
    return pd.concat((binary_columns_imputed, non_binary_columns_imputed), axis=1)


def Knn_imputer(features_array, strategy='mean'):
    binary_columns_gaps, non_binary_columns_gaps = get_binary_columns(features_array)
    imputer_binary = KNNImputer(n_neighbors=5)#SimpleImputer(strategy='most_frequent')#SimpleImputer(strategy='constant', fill_value=0)#SimpleImputer(strategy=strategy)
    imputer_binary.fit(binary_columns_gaps)
    binary_columns_imputed = imputer_binary.transform(binary_columns_gaps)
    binary_columns_imputed = pd.DataFrame(binary_columns_imputed, index = binary_columns_gaps.index, columns = binary_columns_gaps.columns)
    binary_columns_imputed = round(binary_columns_imputed)

    imputer = KNNImputer(n_neighbors=5)#SimpleImputer(strategy='constant', fill_value=0)
    imputer.fit(non_binary_columns_gaps)
    non_binary_columns_imputed = imputer.transform(non_binary_columns_gaps)
    non_binary_columns_imputed = pd.DataFrame(non_binary_columns_imputed, index = non_binary_columns_gaps.index, columns = non_binary_columns_gaps.columns)
    return pd.concat((binary_columns_imputed, non_binary_columns_imputed), axis=1)



def delete_nan_values(features_array, y_labels):
    ar_nan = np.where(pd.isnull(features_array))
    unique_rows_with_inf_values = np.unique(ar_nan[0])
    features_array = np.delete(features_array, unique_rows_with_inf_values, axis=0)
    y_labels = np.delete(y_labels, unique_rows_with_inf_values, axis=0)
    return features_array, y_labels
    # print(ar_nan)


def get_binary_columns(feature_array):
    #feature_array = oX
    binary_columns_indices = []
    for i in range(feature_array.shape[1]):
        #print(i)
        column = feature_array.iloc[:, i]
        is_binary = np.all((column == 0) | (column == 1) | (np.isnan(column)))
        if is_binary:
            binary_columns_indices.append(i)
    binary_columns = feature_array.iloc[:, binary_columns_indices]
    mask = np.ones(feature_array.shape[1], np.bool)
    mask[binary_columns_indices] = 0
    non_binary_columns = feature_array.loc[:, mask]
    return binary_columns, non_binary_columns



# UNIQUE
def unique(list1): 
    list_set = set(list1) 
    unique_list = (list(list_set)) 
    return(unique_list)

# INTERSECTION
def intersection(List1, List2):
    Set1 = set(List1)
    Set2 = set(List2)
    return Set1.intersection(Set2)