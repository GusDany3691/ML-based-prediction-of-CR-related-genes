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

from MLfunctions import BRFparamGridFunc, EECparamGridFunc, CATparamGridFunc, BBCparamGridFunc, GridCreator, GeometricMean, FiPreProcessing
from MLfunctions import PreProcessing, backwardElimination,  Sampling, unique, NanInfs, simple_imputer, Knn_imputer, delete_nan_values, get_binary_columns, unique, intersection

################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################

# --- DEFINE CLASSIFICATION DATASET AND PARAMETERS ----------------------------------------------------------------------------------------------------------
# PATHS

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
ImportanceMethod = 'Gini' # Gini or Permutation

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


##############################################################################################################################################################
##### FEATURE IMPORTANCE #####################################################################################################################################
##############################################################################################################################################################


sXtrain, Ytrain = Sampling(X, y, opc = SamplingMethod) # opc = {normal, under, over, smote, blsmote, csmote, svmsmote, smoteenn, smotetomek, adasyn}

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

            Xtrain, vLen, uvLen, cuvLen, pcuvLen, VUCPtext

            bXtrainK, b_vLenImpK, b_cvLenImpK, b_pcvLenImpK, b_pcuvLenImpK, VUCPtext = FiPreProcessing(biXtrain, MaxClassPercent = MCP, CorrCoef = CorrCoefK, SL = SL, VUCPS = np.multiply(VCPS,bMask), Minoc = MinocK, Mds = 'B')
            print("Original Size     :", b_oLenImpK , "\n")
            print("Near Zero Var Size:", b_vLenImpK , "\n")
            print("Nzv + Corr Size   :", b_cvLenImpK , "\n")

            cXtrainK, c_vLenImpK, c_cvLenImpK, c_pcvLenImpK, c_pcuvLenImpK, VUCPtext = FiPreProcessing(ciXtrain, MaxClassPercent = MCP, CorrCoef = CorrCoefK, SL = SL, VUCPS = np.multiply(VCPS,cMask), Minoc = MinocK, Mds = 'C')
            print("Original Size     :", c_oLenImpK , "\n")
            print("Near Zero Var Size:", c_vLenImpK , "\n")
            print("Nzv + Corr Size   :", c_cvLenImpK , "\n")

            XtrainK = pd.concat((bXtrainK, cXtrainK), axis=1)

            oLenImpK = b_oLenImpK + c_oLenImpK
            vLenImpK = b_vLenImpK + c_vLenImpK
            cvLenImpK = b_cvLenImpK + c_cvLenImpK
            pcvLenImpK = b_pcvLenImpK + c_pcvLenImpK
            pcuvLenImpK =  b_pcuvLenImpK +  c_pcuvLenImpK

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

##### GET FEATURES IMPORTANCE #################################################################################################################################

# GINI
if ImportanceMethod == 'Gini':
    fi = pd.DataFrame({'Feature': list(Xtrain.columns), 'Score': BestModel.feature_importances_})
    fi = fi.sort_values('Score', ascending = False)
    fi.index = fi.Feature
    #fi

# PERMUTATION
if ImportanceMethod == 'Permutation':
    oi = permutation_importance(BestModel, Xtrain, Ytrain, n_repeats=30, random_state=0)
    Ftrs = []
    Impts = []
    for i in oi.importances_mean.argsort()[::-1]:
    #if oi.importances_mean[i] - 2 * oi.importances_std[i] > 0:
        print(f"{Xtrain.columns[i]:<8}"
                f"{oi.importances_mean[i]:.3f}"
                f" +/- {oi.importances_std[i]:.3f}")
        Ftrs.append(Xtrain.columns[i])
        Impts.append(oi.importances_mean[i])
    fi = pd.DataFrame({'Feature': Ftrs, 'Score': Impts}, index = Ftrs)


fi['ScaledScore'] = fi.Score * 100 / (max(fi.Score))

fi.to_csv('/ML-based-prediction-of-CR-related-genes/Data/Features_Importance/Original/fi_'+nDataset+'_'+Classifier+'.csv')

