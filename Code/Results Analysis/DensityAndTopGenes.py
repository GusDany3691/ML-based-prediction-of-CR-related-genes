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
GoDensity  = pd.read_csv('/ML-based-prediction-of-CR-related-genes/Data/CR-probabilities/Strongest models/BRF_GO.csv')
DipDensity = pd.read_csv('/ML-based-prediction-of-CR-related-genes/Data/CR-probabilities/Strongest models/CAT_DIP.csv')

oDIP = DipDensity.sort_values(by='Prob', ascending=False)
oGO = GoDensity.sort_values(by='Prob', ascending=False)



###########################################################################
Top = [16.9, 2.75]            # CR-probability / NotCR probability text positions for {GO Temrs, BRF} and {PathDIP, CAT}

# FIGURE CREATION
fig, axs = plt.subplots(2,2)
axs[0,0]
axs[0,1]

# PARAMETERS
sp = 2 # size parameter
TitleSize = 12
LabelSize = 11
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

plt.yticks(style='italic')



title_list = ['{GO terms, BRF}', '{PathDip, CAT}']

#########################################################################################################################################################
### DENSITY #############################################################################################################################################
#########################################################################################################################################################

sns.set(context='notebook', style='darkgrid', palette='deep', font='sans-serif', font_scale=1, color_codes=False, rc=None)

xleyend = 0.77
yleyend = 0.76

row = 0
DensityList = [GoDensity, DipDensity]

NotCrDIP = oDIP[oDIP['Class'] == 'NotCR']
NotCrDIP = NotCrDIP.sort_values(by='Prob', ascending=False)
NotCrDIP.index = range(len(NotCrDIP))
NotCrDIP = NotCrDIP.loc[range(ElementsInList),]

oGO.index = range(len(oGO))
NotCrGO = oGO3[oGO['Class'] == 'NotCR']
NotCrGO = NotCrGOiloc[range(ElementsInList),]

top_list = [NotCrGO, NotCrDIP]

intro_list = ['(A) CR-probabilities density:\n', '(B) CR-probabilities density:\n']

cont = 0 
for cont in range(2):
 
    title1 = intro_list[cont]
    title2 = title_list[cont]

    title = title1 + title2

    top = top_list[cont]

    FPmin = top['Prob'].iloc[TopGraph]

    classes = ['CR', 'NotCR']
    Density = DensityList[cont]

    # Iterate through the five airlines

    colores = ['#D62728', '#1F77B4']
    cont1 = 0
    for clase in classes:

        # Subset to CR-relationship
        subset = Density[Density['Class'] == clase] 
        xmin = min(subset['Prob'])
        xmax = max(subset['Prob'])
        
        # Draw the density plot
        sns.distplot(subset['Prob'], hist = False, kde = True,
                     kde_kws = {'shade': True, 'linewidth': 2},
                     label = clase, ax=axs[row,cont], color=colores[cont1])
        cont1 = cont1+1
    axs[row,cont].axvline(0.5, color=(0, 0, 0), label='Threshold', linestyle='--')
    axs[row,cont].axvline(FPmin, color='#55a868', label='Top ten FP', linestyle='--')
    axs[row,cont].set_ylabel('Genetic density', size = LabelSize)#, size = LabelSize)
    axs[row,cont].set_xlabel('CR-probability', size = LabelSize)#, size = LabelSize)
    axs[row,cont].set_title(title, size = TitleSize, fontweight = TitleWeight)#, size = TitleSize, fontweight = 'bold')
    axs[row,cont].tick_params(axis='both', which='major', labelsize=TickSize)
    handles, _ = axs[row,cont].get_legend_handles_labels()


    axs[row,cont].legend(ncol=1, loc = 'center',#loc='upper center', # upper center
          bbox_to_anchor=[xleyend, yleyend], # [0.5, 1.25] # Here, I can play with the legend box location
          handles = handles[0:],# labels = labels,
          columnspacing=1.3, labelspacing=0.0,
          handletextpad=0.0, handlelength=1.5, fontsize = LegendSize)

    
    # Lims
    axs[row,cont].set_xlim([min(Density['Prob']), max(Density['Prob'])])

    # CR-NotCR Text: GO
    fs = 9
    Ltext_pos = ( min(DensityList[cont]['Prob']) + 0.5 ) / 2
    Rtext_pos = ( max(DensityList[cont]['Prob']) + 0.5 ) / 2

# both no mapped
    axs[row,cont].text(Ltext_pos, Top[cont], 'NotCR-predicted', fontsize = fs, horizontalalignment = 'center', weight='bold')#, ha='left', va='left')#,ha='center', va='center')
    axs[row,cont].text(Rtext_pos , Top[cont], 'CR-predicted', fontsize = fs, weight='bold', horizontalalignment = 'center')#, ha='left', va='left')#,ha='center', va='center')


#########################################################################################################################################################
### DISCOVERY ###########################################################################################################################################
#########################################################################################################################################################

xleyend = 0.25
yleyend = 0.9

row = 1
f = [0.02, 0.00]

intro_list = ['(C) Top 10 FP genes:\n', '(D) Top 10 FP genes:\n']

for cont in range(2):

    plt.yticks(style='italic')

    axs[row,cont].tick_params(axis='both', which='major', labelsize=TickSize)

    title1 = intro_list[cont]
    title2 = title_list[cont]

    title = title1 + title2

    top = top_list[cont]
    sns.barplot(data=top, x="Prob", y="Label", saturation=0.6, ax=axs[row,cont], color = ('#AC6668'))   # Skin: E6D0C9  #Brown: #C18C6E  #Red: #AC6668
    axs[row,cont].set_title(title, size = TitleSize, fontweight = TitleWeight)                          #, fontweight = 'bold')#, size = TitleSize, fontweight = 'bold')
    axs[row,cont].axvline(0.5, color=(0, 0, 0), label='Threshold', linestyle='--')
    for p in axs[row,cont].patches:
        width = p.get_width()
        axs[row,cont].text(p.get_width() - 0.08 + f[cont], p.get_y() + 0.55*p.get_height(), 
           '{:1.2f}'.format(width),
           ha='center', va='center', fontsize = 10, weight = 'bold')

        params = {'mathtext.default': 'regular' }          
        plt.rcParams.update(params)
        axs[row,cont].set_xlabel('CR-probability', size = LabelSize)
        axs[row,cont].set_ylabel('$Ageing_{NotCR}$-related\n genes', size = LabelSize)

    handles, _ = axs[row,cont].get_legend_handles_labels()

    axs[row,cont].legend(ncol=1, loc = 'center',  #loc='upper center', # upper center
        bbox_to_anchor=[xleyend, yleyend],        # [0.5, 1.25] # Here, I can play with the legend box location
        handles = handles[0:],                    # labels = labels,
        columnspacing=1.3, labelspacing=0.0,
        handletextpad=0.0, handlelength=1.5, 
        fontsize = LegendSize)


plt.tight_layout(pad=0.1, w_pad=0.5, h_pad=-1)


for tick in axs[1][0].yaxis.get_major_ticks():
  tick.label1.set_fontstyle('italic')


plt.show()


