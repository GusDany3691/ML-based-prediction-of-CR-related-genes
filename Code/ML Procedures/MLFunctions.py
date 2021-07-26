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


##### GRID ####################################################################################################################################################

def BRFparamGridFunc(Xtrain):
     BRFparamGrid = {
        'bootstrap': [True],
        'replacement': [True, False],
        'max_features': ['sqrt','log2'], 
        'n_estimators': [500],
        'class_weight': ['balanced', 'balanced_subsample', None],
        'sampling_strategy': [1],
     }
     return(BRFparamGrid)

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


def XGBparamGridFunc(Xtrain):
     XGBparamGrid = {
        'num_feature': ['sqrt','log2'],
     }
     return BBCparamGrid


##### SCORE ###################################################################################################################################################

def GeometricMean(Ytest,Ypred):
   GM = geometric_mean_score(Ytest, Ypred)
   return GM

##### PRE PROCESSING ################################################################################################################################################

def FiPreProcessing(sXtrain, Ytrain, MaxClassPercent = 1, Ftop = 2000, UnivariateTechnique = mutual_info_classif, CorrCoef = 1, SL = 0.5, VUCPS = [True, False, False, False, True], Minoc = 1, Mds = 'M'):

    Mtext = Mds

    # --- INDICES AND COLUMNS ---------------------------------------------------------------------------------------------------------------------------------------
    Features = sXtrain.columns
    TrainIndices = sXtrain.index

    Vtext = ''
    Ctext = ''
    Utext = ''
    Ptext = ''
    Stext = ''

    # --- FEATURES SELECTION (NEAR ZERO VARIANCE) -------------------------------------------------------------------------------------------------------------------
    if (VUCPS[0]):
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

    # --- FEATURES SELECTION (UNIVARIATE FILTER) ---------------------------------------------------------------------------------------------------------------------

    if (VUCPS[1]):
        Utext = 'u' + str(Ftop)# + 'v' + str(round(MaxClassPercent*100))
        fs = SelectKBest(score_func=UnivariateTechnique, k=Ftop)
        fs.fit(vXtrain, Ytrain)
        fs.transform(vXtrain)
        UnivariateFeatures = pd.DataFrame(fs.scores_, index =  vXtrain.columns, columns = ['Values'])
        SortedUnivariateFeatures = UnivariateFeatures.sort_values(by=['Values'], ascending=False)
        uvSelectedColumns = SortedUnivariateFeatures.index[0:Ftop]
    else:
        uvSelectedColumns = vSelectedColumns
    uvXtrain = vXtrain[uvSelectedColumns]
    uvLen = len(uvSelectedColumns)

    # --- FEATURES SELECTION (CORRELATION) ---------------------------------------------------------------------------------------------------------------------------
    # CORRELATION SELECTION
    if(VUCPS[2]):
        Ctext = 'c' + str(round(CorrCoef*100))
        corr = uvXtrain.corr().abs()
        #sns.heatmap(corr)
        columns = np.full((corr.shape[0],), True, dtype=bool)
        for i in range(corr.shape[0]):
            for j in range(i+1, corr.shape[0]):
                if corr.iloc[i,j] >= CorrCoef:
                    if columns[j]:
                        columns[j] = False
        cuvSelectedColumns = uvXtrain.columns[columns]
    else:
        cuvSelectedColumns = uvSelectedColumns
    cuvXtrain = vXtrain[cuvSelectedColumns]
    cuvLen = len(cuvSelectedColumns)
    
    # PVALUES SELECTION (BACKWARD ELIMINATION)
    if (VUCPS[3]):
        Ptext = 'p' + str(round(SL))
        pcuvSelectedColumns = backwardElimination(cuvXtrain.values, Ytrain, SL, cuvSelectedColumns)
    else:
        pcuvSelectedColumns = cuvSelectedColumns
    pcuvXtrain = cvXtrain[pcuvSelectedColumns]
    pcuvLen = len(pcuvSelectedColumns)

    # --- SCALER -----------------------------------------------------------------------------------------------------------------------------------------------------
    if (VUCPS[4]):
        Stext = 's' 
        sc = StandardScaler()
        spcuvXtrain = sc.fit_transform(pcuvXtrain)
        #spcvXtest = sc.transform(pcvXtest)
    else:
        spcuvXtrain = pcuvXtrain
        #spcvXtest = pcvXtest

    # --- DATA FRAMES ------------------------------------------------------------------------------------------------------------------------------------------------
    Xtrain = pd.DataFrame(data=spcuvXtrain, columns = pcuvSelectedColumns, index = TrainIndices)
    VUCPtext = Mtext + Vtext +  Utext + Ctext + Ptext + Stext

    return Xtrain, vLen, uvLen, cuvLen, pcuvLen, VUCPtext

sXtrain = iXtrain
oXtest = iXtest



def PreProcessing(sXtrain, oXtest, Ytrain, MaxClassPercent = 1, Ftop = 2000, UnivariateTechnique = mutual_info_classif, CorrCoef = 1, SL = 0.5, VUCPS = [True, False, False, False, True], Minoc = 1,  Mds = 'M'):

    Mtext = Mds

    # --- INDICES AND COLUMNS ---------------------------------------------------------------------------------------------------------------------------------------
    Features = sXtrain.columns
    TrainIndices = sXtrain.index
    TestIndices = oXtest.index
    Vtext = ''
    Utext = ''
    Ctext = ''
    Ptext = ''
    Stext = ''

    # --- FEATURES SELECTION (NEAR ZERO VARIANCE) -------------------------------------------------------------------------------------------------------------------
    if (VUCPS[0]):
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


    # --- FEATURES SELECTION (UNIVARIATE FILTER) ---------------------------------------------------------------------------------------------------------------------

    if (VUCPS[1]):
        Utext = 'u' + str(Ftop)# + 'v' + str(round(MaxClassPercent*100))
        fs = SelectKBest(score_func=UnivariateTechnique, k=Ftop)
        fs.fit(vXtrain, Ytrain)
        fs.transform(vXtrain)
        UnivariateFeatures = pd.DataFrame(fs.scores_, index =  vXtrain.columns, columns = ['Values'])
        SortedUnivariateFeatures = UnivariateFeatures.sort_values(by=['Values'], ascending=False)
        uvSelectedColumns = SortedUnivariateFeatures.index[0:Ftop]
    else:
        uvSelectedColumns = vSelectedColumns
    uvXtrain = vXtrain[uvSelectedColumns]
    uvLen = len(uvSelectedColumns)


    # --- FEATURES SELECTION (CORRELATION) ---------------------------------------------------------------------------------------------------------------------------
    # CORRELATION SELECTION
    if(VUCPS[2]):
        Ctext = 'c' + str(round(CorrCoef*100))
        corr = uvXtrain.corr().abs()
        #sns.heatmap(corr)
        columns = np.full((corr.shape[0],), True, dtype=bool)
        for i in range(corr.shape[0]):
            for j in range(i+1, corr.shape[0]):
                if corr.iloc[i,j] >= CorrCoef:
                    if columns[j]:
                        columns[j] = False
        cuvSelectedColumns = uvXtrain.columns[columns]
    else:
        cuvSelectedColumns = uvSelectedColumns
    cuvXtrain = vXtrain[cuvSelectedColumns]
    cuvLen = len(cuvSelectedColumns)
    
    # --- PVALUES SELECTION (BACKWARD ELIMINATION) -------------------------------------------------------------------------------------------------------------------
    if (VUCPS[3]):   
        Ptext = 'p' + str(round(SL))
        pcuvSelectedColumns = backwardElimination(cuvXtrain.values, Ytrain, SL, cuvSelectedColumns)
    else:
        pcuvSelectedColumns = cuvSelectedColumns
    pcuvXtrain = cuvXtrain[pcuvSelectedColumns]
    pcuvLen = len(pcuvSelectedColumns)

    # --- XTEST ------------------------------------------------------------------------------------------------------------------------------------------------------
    pcuvXtest = oXtest[pcuvSelectedColumns]

    # --- SCALER -----------------------------------------------------------------------------------------------------------------------------------------------------

    if (VUCPS[4]):
        Stext = 's' 
        sc = StandardScaler()
        spcuvXtrain = sc.fit_transform(pcuvXtrain)
        spcuvXtest = sc.transform(pcuvXtest)
    else:
        spcuvXtrain = pcuvXtrain
        spcuvXtest = pcuvXtest

    # --- DATA FRAMES ------------------------------------------------------------------------------------------------------------------------------------------------
    Xtrain = pd.DataFrame(data=spcuvXtrain, columns = pcuvSelectedColumns, index = TrainIndices)
    Xtest = pd.DataFrame(data=spcuvXtest, columns = pcuvSelectedColumns, index = TestIndices)
    VUCPtext = Mtext + Vtext + Utext + Ctext + Ptext + Stext

    return Xtrain, Xtest, vLen, uvLen, cuvLen, pcuvLen, VUCPtext



# feature selection
def FstatisticFilter(X_train, y_train, X_test):
	# configure to select all features
	fs = SelectKBest(score_func=f_classif, k=10)
	# learn relationship from training data
	fs.fit(X_train, y_train)
	# transform train input data
	X_train_fs = fs.transform(X_train)
	# transform test input data
	X_test_fs = fs.transform(X_test)
	return X_train_fs, X_test_fs, fs



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


##### PROTEINS SPLITTING ##################################################################################################################

def TranscriptsSpliting(X, y, Splits = 10, random_state = 50):
    #Dataset = TrainingSets[0]

    WholeInstances = X.index

    ##### GETTING GENES ####################################################################################################################

    # GENE INSTANCES
    GeneInstancesList = []
    for GeneProteinInstance in WholeInstances:
      GeneInstancesList.append(GeneProteinInstance.split('-')[0])
 
    # GENE INDICES
    GeneIndicesDict = {}
    for GeneInstance in GeneInstancesList:
        GeneIndicesDict[GeneInstance] = []
    for i, GeneInstance in enumerate(GeneInstancesList):
        GeneIndicesDict[GeneInstance].append(i)

    # GENE NAMES
    YgeneList = []
    YgeneDict = {}
    Ygene = []
    uGeneInstancesList = [*GeneIndicesDict]                # Get key names from dictionary
    ClassIndicesDict = {}
    for GeneName in uGeneInstancesList:
        ClassIndicesDict[GeneName] = []
    for GeneName in uGeneInstancesList:
        GeneIndex = GeneIndicesDict[GeneName][0]
        ClassIndicesDict[GeneName].append( GeneIndex )
        YgeneDict[GeneName] = y[GeneIndex]
        Ygene.append(y[GeneIndex])

    # --- SPLITTING ---------------------------------------------------------------------------------------------------------------------------------------

    gTrainingIndicesList = []
    gTestingIndicesList = []
    cont = 0

    Skf = StratifiedKFold(n_splits=Splits, shuffle=True, random_state = random_state)

    for TrainingIndices, TestingIndices in Skf.split(Ygene, Ygene):
        gTrainingIndicesList.append(TrainingIndices)
        gTestingIndicesList.append(TestingIndices)
        cont = cont+1

    #gene2pindex(Indices = gTestingIndicesList[0], GeneIndicesDict = GeneIndicesDict)
    TestingSplits = [gene2pindex(Indices = gTestingIndices, GeneIndicesDict = GeneIndicesDict) for gTestingIndices in gTestingIndicesList]
    TrainingSplits = [gene2pindex(Indices = gTrainingIndices, GeneIndicesDict = GeneIndicesDict) for gTrainingIndices in gTrainingIndicesList]

    #TrainingSets = [Dataset.iloc[pTrainingIndices] for pTrainingIndices in OuterTrainingSplits]
    #TestingSets = [Dataset.iloc[pTestingIndices] for pTestingIndices in OuterTestingSplits]

    FoldIterator = list(CustomCv(TrainingSplits, TestingSplits))
    
    Intersect = TestTrainIntersection(TrainingSplits, TestingSplits, GeneInstancesList)
    
    return FoldIterator, Intersect#OuterTrainingSplits, OuterTestingSplits, GeneInstancesList, TrainingSets, TestingSets



##### CUSTOM CV #######################################################################################################################
def CustomCv(TrainingSplits, TestingSplits):
    folds = len(TrainingSplits)
    i = 0
    while i < folds:
        yield np.array(TrainingSplits[i]), np.array(TestingSplits[i])
        i += 1
