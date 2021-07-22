import numpy as np
from scipy import interp
import pylab as pl
import matplotlib.pyplot as plt
from sklearn import svm, datasets
from sklearn.metrics import roc_curve, auc

from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm
from matplotlib.ticker import MultipleLocator
import itertools
import pandas as pd
import seaborn as sns
import matplotlib.ticker as ticker

import matplotlib as mpl
from scipy.interpolate import interp1d
import math

from sklearn.metrics import roc_curve, auc
from imblearn.metrics import geometric_mean_score, sensitivity_score, specificity_score
from sklearn.metrics import roc_auc_score, roc_curve, auc, confusion_matrix, make_scorer, classification_report, plot_confusion_matrix
import sklearn

from scipy.interpolate import interp1d


# GO AND DIP
GoDensity  = pd.read_csv('E:\\The Project\\NewPred\\eGO.csv')
DipDensity = pd.read_csv('E:\\The Project\\NewPred\\oDIP.csv')

# FIGURE CREATION
fig, axs = plt.subplots(1,2)

# PARAMETERS
sp = 2 # size parameter
TitleSize = 12
LabelSize = 12
TickSize = 10
ConfusionSize = 12
LegendSize = 10
normalize    = False
target_names = ['CR', 'NotCR']
cmap='Blues'
LegendSize = 10
BigTitle = 12
TitleWeight = 'bold'

ElementsInList = 10

TopGraph = 10
TopGraph = TopGraph - 1


#########################################################################################################################################################
### CONFUSION COMPUTATION ###############################################################################################################################
#########################################################################################################################################################

from imblearn.metrics import geometric_mean_score, sensitivity_score, specificity_score
from sklearn.metrics import roc_auc_score, roc_curve, auc, confusion_matrix, make_scorer, classification_report, plot_confusion_matrix
import numpy as np
import pandas as pd


##### PATHDIP - CAT #####################################################################################################################################

Labels = pd.read_csv("E:\\The Project\\NewPred\\y_pred_test\\labels_pathdip_cat.csv")
pProb = pd.read_csv("E:\\The Project\\NewPred\\y_pred_test\\y_pred_pathdip_cat.csv")
Test = pd.read_csv("E:\\The Project\\NewPred\\y_pred_test\\y_test_pathdip_cat.csv")

Labels = Labels.drop(columns = ['Unnamed: 0'])
pProb = pProb.drop(columns = ['Unnamed: 0'])
Test = Test.drop(columns = ['Unnamed: 0'])

Prob = pProb.iloc[:,[1,3,5,7,9,11,13,15,17,19]]

Pred = Prob.iloc[:,] >= 0.5
Pred = Pred.replace(True,1)

Test.iloc[:,0]
Pred.iloc[:,0]


Sens = np.array([]) 
Spec = np.array([]) 
Gmean = np.array([]) 
Auc = np.array([]) 
fTest = np.array([]) 
fPred = np.array([]) 
fProb = np.array([]) 
fLabel = np.array([]) 
for i in range(10):
    print(i)
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

cm_dip = confusion_matrix(fTest, fPred)



##### GO TERMS - BRF #####################################################################################################################################

oGO = pd.read_csv("E:/The Project/BRF_GO.csv")

cm_go = confusion_matrix(oGO[['Test']], oGO[['Pred']])


#########################################################################################################################################################
###### CONFUSION PLOTTING ##############################################################################################################################
#########################################################################################################################################################

row = 0

# ORIGINAL GO and DIP

cm_list = [cm_go, cm_dip]
title_list = ['{GO terms, BRF}', '{PathDip, CAT}']
intro_list = ['(A) Confusion Matrix:\n', '(B) Confusion Matrix:\n']

cont = 0
for cont in range(2):

    cm = cm_list[cont]
    thresh = cm.max() / 2

    title1 = intro_list[cont]
    title2 = title_list[cont]

    title = title1 + title2

    # SET PARAMETERS
    axs[cont].set_title(title, size = TitleSize, fontweight = TitleWeight)#, size = TitleSize, fontweight = 'bold')
    im = axs[cont].imshow(cm, interpolation='nearest', cmap=cmap)
    divider = make_axes_locatable(axs[cont])
    caxs = divider.append_axes("right", size="20%", pad=0.05)
    cbar = plt.colorbar(im, cax=caxs)
    axs[cont].set_yticks([0.0, 1.0])
    axs[cont].set_yticks([0.0, 1.0])
    axs[cont].set_yticklabels(target_names, size = TickSize)#, size = TickSize)    
    axs[cont].set_xticks([0.0, 1.0])
    axs[cont].set_xticklabels(target_names, size = TickSize)#, size = TickSize)   
    axs[cont].set_ylabel('Predicted', size = LabelSize)#, size = LabelSize)
    axs[cont].set_xlabel('Actual', size = LabelSize)#, size = LabelSize)
    axs[cont].grid(False)

    # SHOW NUMBERS
    thresh = cm.max() / 1.5 if normalize else cm.max() / 2
    for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
            axs[cont].text(j, i, "{:,}".format(cm[i, j]), horizontalalignment="center", size = ConfusionSize, color="white" if cm[i, j] > thresh else "black")

plt.tight_layout(pad=0.1, w_pad=0.5, h_pad=-1)

plt.show()


