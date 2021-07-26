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


# FIGURE CREATION
fig, axs = plt.subplots(1,2)

# FIGURE PARAMETERS
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

# GO AND DIP
GoDensity  = pd.read_csv('/ML-based-prediction-of-CR-related-genes/Data/CR-probabilities/Strongest models/BRF_GO.csv')
DipDensity = pd.read_csv('/ML-based-prediction-of-CR-related-genes/CR-probabilities/Strongest models/CAT_DIP.csv')

cm_dip = confusion_matrix(GoDensity['Test'], GoDensity['Pred'])
cm_go = confusion_matrix(DipDensity['Test'], DipDensity['Pred'])


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


