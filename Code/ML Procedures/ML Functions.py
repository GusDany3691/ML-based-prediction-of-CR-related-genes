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
