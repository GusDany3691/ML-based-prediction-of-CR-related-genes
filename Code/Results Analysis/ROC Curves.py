from imblearn.metrics import geometric_mean_score, sensitivity_score, specificity_score
from sklearn.metrics import roc_auc_score, roc_curve, auc, confusion_matrix, make_scorer, classification_report, plot_confusion_matrix
import numpy as np
import pandas as pd

#### RETRIEVE DATA ######################################################################################################################

Frame = {}
Frame['DIP'] = {}
Frame['GO'] = {}
Frame['KEGG'] = {}
Frame['WHOLE'] = {}

Datasets= ['DIP','GO']

# - BRF MODELS --------------------------------------------------------------------------------------------------------------------------------

for dataset in Datasets:
    if dataset == 'DIP':
        oFrame = pd.read_csv('E:/The Project/BRF_DIP.csv')
    if dataset == 'GO':
        oFrame = pd.read_csv('E:/The Project/BRF_GO.csv')
    Frame[dataset]['BRF'] = FlatternGus(oFrame)


# - OTHER MODELS --------------------------------------------------------------------------------------------------------------------------------

Folder = "E:\\The Project\\NewPred\\y_pred_test\\"

Algorithms = ['EEC', 'XGB', 'CAT']

for i in range(4):
    for j in range(3):

        dataset = Datasets[i]        #'GO'
        algorithm = Algorithms[j]    #'EEC'

        # PATHDIP EEC
        if dataset == 'DIP' and algorithm == 'EEC':
            LabelsFile = "labels_pathdip_ee.csv"
            ProbFile = "y_pred_pathdip_ee.csv"
            TestFile = "y_test_pathdip_ee.csv"

        # PATHDIP XBG
        if dataset == 'DIP' and algorithm == 'XGB':
            LabelsFile = "labels_pathdip_xgboost.csv"
            ProbFile = "y_pred_pathdip_xgboost.csv"
            TestFile = "y_test_pathdip_xgboost.csv"

        # PATHDIP CAT
        if dataset == 'DIP' and algorithm == 'CAT':
            LabelsFile = "labels_pathdip_cat.csv"
            ProbFile = "y_pred_pathdip_cat.csv"
            TestFile = "y_test_pathdip_cat.csv"


        # GO EEC
        if dataset == 'GO' and algorithm == 'EEC':
            LabelsFile = "labels_go_ee.csv"
            ProbFile = "y_pred_go_ee.csv"
            TestFile = "y_test_go_ee.csv"

        # GO XBG
        if dataset == 'GO' and algorithm == 'XGB':
            LabelsFile = "labels_go_xgb.csv"
            ProbFile = "y_pred_go_xgb.csv"
            TestFile = "y_test_go_xgb.csv"

        # GO CAT
        if dataset == 'GO' and algorithm == 'CAT':
            LabelsFile = "labels_go_cat.csv"
            ProbFile = "y_pred_go_cat.csv"
            TestFile = "y_test_go_cat.csv"


         # KEGG EEC
        if dataset == 'KEGG' and algorithm == 'EEC':
            LabelsFile = "labels_kegg_path_ee.csv"
            ProbFile = "y_pred_kegg_path_ee.csv"
            TestFile = "y_test_kegg_path_ee.csv"

        # KEGG XBG
        if dataset == 'KEGG' and algorithm == 'XGB':
            LabelsFile = "labels_kegg_path_xgb.csv"
            ProbFile = "y_pred_kegg_path_xgb.csv"
            TestFile = "y_test_kegg_path_xgb.csv"

        # KEGG CAT
        if dataset == 'KEGG' and algorithm == 'CAT':
            LabelsFile = "labels_kegg_path_cat.csv"
            ProbFile = "y_pred_kegg_path_cat.csv"
            TestFile = "y_test_kegg_path_cat.csv"


        # WHOLE EEC
        if dataset == 'WHOLE' and algorithm == 'EEC':
            LabelsFile = "y_test_wholeimp_ee.csv"
            ProbFile = "y_pred_wholeimp_ee.csv"
            TestFile = "y_test_wholeimp_ee.csv"

        # KEGG XBG
        if dataset == 'WHOLE' and algorithm == 'XGB':
            LabelsFile = "y_test_wholeimp_xg.csv"
            ProbFile = "y_pred_wholeimp_xg.csv"
            TestFile = "y_test_wholeimp_xg.csv"


        # KEGG CAT
        if dataset == 'WHOLE' and algorithm == 'CAT':
            LabelsFile = "y_test_wholeimp_cat.csv"
            ProbFile = "y_pred_wholeimp_cat.csv"
            TestFile = "y_test_wholeimp_cat.csv"


        # ---  DATA RETRIEVAL -------------------------------------------------------------------------------------------------------------------

        oLabels = pd.read_csv(Folder + LabelsFile)
        oProb = pd.read_csv(Folder + ProbFile)
        oTest = pd.read_csv(Folder + TestFile)

        Labels, Prob, Pred, Test = FixData(oLabels, oProb, oTest)
        Frame[dataset][algorithm] = Flattern(Labels, Prob, Test)
        Frame[dataset][algorithm]['Frame']



dataset = 'DIP'
algorithm = 'CAT'


Frame[dataset][algorithm]['Frame']['Pred'] == 1.0 and (Frame[dataset][algorithm]['Frame']['Test'] == 0.0).values()
poss = (Frame[dataset][algorithm]['Frame']['Pred'] == 1.0) & (Frame[dataset][algorithm]['Frame']['Test'] == 0.0)

Frame[dataset][algorithm]['Frame'].to_csv('E:/The Project/EEC_DIP.csv')

range(311,332)

dataset = 'GO'
algorithm = 'EEC'

Frame[dataset][algorithm]['wMeasures']['Auc']
Frame[dataset][algorithm]['wMeasures']['Gmean']
Frame[dataset][algorithm]['wMeasures']['Conf']

Frame[dataset][algorithm]['aMeasures']['Gmean']
Frame[dataset][algorithm]['aMeasures']['Auc']

Frame[dataset][algorithm]['aMeasures']['Sens']
Frame[dataset][algorithm]['aMeasures']['Spec']

Frame[dataset][algorithm]['aMeasures']['Sens']
Frame[dataset][algorithm]['aMeasures']['Spec']

Frame[dataset][algorithm]['wMeasures']['Sens']
Frame[dataset][algorithm]['wMeasures']['Spec']

Frame[dataset][algorithm]['aMeasures']['Sens']
Frame[dataset][algorithm]['aMeasures']['Spec']

Folder + LabelsFile
Folder + ProbFile
Folder + TestFile

##### SAVE ############################################################################################################################

file = open("E:\The Project\Frame.py", 'wb')
pickle.dump(Frame, file)
file.close()   

##### TESTING #########################################################################################################################

dataset = 'GO'
algorithm = 'BRF'

Frame[dataset][algorithm]['wMeasures']['Auc']
Frame[dataset][algorithm]['wMeasures']['Gmean']


Frame['DIP']['BRF']['wMeasures']['Auc']
Frame['DIP']['BRF']['wMeasures']

Frame['DIP']['EEC']['aMeasures']['Auc']
Frame['DIP']['EEC']['wMeasures']

Frame['DIP']['XGB']['aMeasures']['Auc']
Frame['DIP']['XGB']['aMeasures']['Gmean']

Frame['DIP']['CAT']['aMeasures']['Auc']
Frame['DIP']['CAT']['aMeasures']['Gmean']


#######################################################################################################################################
##### PLOT DATA #######################################################################################################################
#######################################################################################################################################

import numpy as np
from scipy import interp
import pylab as pl
import matplotlib.pyplot as plt
from sklearn import svm, datasets
from sklearn.metrics import roc_curve, auc
import seaborn as sns

sns.set(context='notebook', style='darkgrid', palette='deep', font='sans-serif', font_scale=1, color_codes=False, rc=None)

fig, axs = plt.subplots(1,2)
axs[0]
axs[1]

Datasets= ['GO','DIP']
Algorithms = ['BRF','EEC', 'XGB', 'CAT']
Colors = ['y-','b-','g-','r-']


Pal = sns.color_palette()
Colors = Pal.as_hex()[0:4]

LegendSize = 10

for cont, dataset in enumerate(Datasets):
    for i, algorithm in enumerate(Algorithms):
        print(dataset + " " + str(algorithm))
        fpr, tpr, _ = roc_curve(Frame[dataset][algorithm]['Frame']['Test'], Frame[dataset][algorithm]['Frame']['Prob'])

        print(i)

        if dataset == 'GO' and algorithm == 'EEC':
            axs[cont].plot(fpr, tpr, Colors[i], label =  algorithm + '_AUC = %0.2f' % 0.83, lw=2)
        else:
            axs[cont].plot(fpr, tpr, Colors[i], label =  algorithm + '_AUC = %0.2f' % Frame[dataset][algorithm]['wMeasures']['Auc'], lw=2)

    axs[cont].plot([0, 1], [0, 1], '--', color=(0.6, 0.6, 0.6), label='Random predictions')       
    if dataset == 'GO':
        axs[cont].set_title('(A) ROC Curves - GO terms', fontweight = 'bold')
        
    if dataset == 'DIP':
        axs[cont].set_title('(B) ROC Curves - PathDIP', fontweight = 'bold')
    axs[cont].set_xlabel('1 - Specificity', size = 12)

    axs[cont].set_ylabel('Sensitivity', size = 12)
    axs[cont].legend(loc="lower right")


    axs[cont].legend(ncol=1,                    #loc = 'center',#loc='upper center', # upper center
        columnspacing=1.3, labelspacing=0.5,
        handletextpad=0.0, handlelength=1.5,
        fontsize = LegendSize)


plt.show()


#######################################################################################################################################
##### FUNCTIONS #######################################################################################################################
#######################################################################################################################################

##### FLATTERN ########################################################################################################################

# FIX DATA
def FixData(Labels, pProb, Test):
    Labels = Labels.drop(columns = ['Unnamed: 0'])
    pProb = pProb.drop(columns = ['Unnamed: 0'])
    Test = Test.drop(columns = ['Unnamed: 0'])
    Prob = pProb.iloc[:,[1,3,5,7,9,11,13,15,17,19]]
    Pred = Prob.iloc[:,] >= 0.5
    Pred = Pred.replace(True,1)
    return Labels, Prob, Pred, Test

# FLATTERN
def Flattern(Labels, Prob, Test):
    Sens = np.array([]) 
    Spec = np.array([]) 
    Gmean = np.array([]) 
    Auc = np.array([]) 
    fTest = np.array([]) 
    fPred = np.array([]) 
    fProb = np.array([]) 
    fLabel = np.array([]) 
    for i in range(10):
        #print(i)
        Te = Test.iloc[:,i]
        Pr = Pred.iloc[:,i]
        Pb = Prob.iloc[:,i]
        Lb = Labels.iloc[:,i]
        Te = Te.dropna()
        Pr = Pr.iloc[range(len(Te))]
        Pb = Pb.iloc[range(len(Te))]
        Lb = Lb.iloc[range(len(Te))]
        fTest = np.append (fTest, Te.values )
        fPred = np.append (fPred, Pr.values )
        fProb = np.append (fProb, Pb.values )
        fLabel = np.append (fLabel, Lb.values )
        Sens = np.append (Sens, sensitivity_score(Te, Pr) )
        Spec = np.append (Spec, specificity_score(Te, Pr) )
        Gmean = np.append (Gmean, geometric_mean_score(Te, Pr) )
        Auc = np.append( Auc, roc_auc_score(Te, Pb) )

    aSens = np.average(Sens)
    aSpec = np.average(Spec)
    aGmean = np.average(Gmean)
    aAuc = np.average(Auc)

    wSens = sensitivity_score(fTest, fPred)
    wSpec = specificity_score(fTest, fPred)
    wGmean = geometric_mean_score(fTest, fPred)
    wAuc = roc_auc_score(fTest, fProb)
    wConf = confusion_matrix(fTest, fPred)

    Frame = pd.DataFrame(data = {'Label': fLabel, 'Test': fTest, 'Pred': fPred, 'Prob': fProb}, index= fLabel)
    Frame.sort_values(by=['Prob'], inplace=True)

    aMeasures = {'Sens': aSens, 'Spec': aSpec, 'Gmean': aGmean, 'Auc': aAuc}
    wMeasures = {'Sens': wSens, 'Spec': wSpec, 'Gmean': wGmean, 'Auc': wAuc, 'Conf': wConf}

    Model = {'Frame': Frame, 'aMeasures': aMeasures, 'wMeasures': wMeasures}

    return Model


# FLATTERN GO
def FlatternGus(oFrame):
    wSens = sensitivity_score(oFrame['Test'], oFrame['Pred'])
    wSpec = specificity_score(oFrame['Test'], oFrame['Pred'])
    wGmean = geometric_mean_score(oFrame['Test'], oFrame['Pred'])
    wAuc = roc_auc_score(oFrame['Test'], oFrame['Prob'])
    wConf = confusion_matrix(oFrame['Test'], oFrame['Pred'])
    wMeasures = {'Sens': wSens, 'Spec': wSpec, 'Gmean': wGmean, 'Auc': wAuc, 'Conf': wConf}

    Model = {'Frame': oFrame, 'wMeasures': wMeasures}

    return Model

