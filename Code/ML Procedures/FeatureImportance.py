from sklearn.inspection import permutation_importance

BestModel = GmeanModel['Model']
ImportanceMethod = 'Permutation'

##### FEATURE IMPORTANCE #############################################################################################################
# GINI
if ImportanceMethod == 'Gini':
    fi = pd.DataFrame({'Feature': list(Xtrain.columns), 'OriginalScore': BestModel.feature_importances_})
    fi = fi.sort_values('OriginalScore', ascending = False)
    fi.index = fi.Feature
    #fi

# PERMUTATION
if ImportanceMethod == 'Permutation':
    oi = permutation_importance(BestModel, Xtest, Ytest, n_repeats=30, random_state=0)
    Ftrs = []
    Impts = []
    for i in oi.importances_mean.argsort()[::-1]:
        #print(oi.importances_mean[i] - 2 * oi.importances_std[i])
        #if oi.importances_mean[i] - 2 * oi.importances_std[i] > 0:
        print(f"{Xtrain.columns[i]:<8}"
                f"{oi.importances_mean[i]:.3f}"
                f" +/- {oi.importances_std[i]:.3f}")
        Ftrs.append(Xtrain.columns[i])
        Impts.append(oi.importances_mean[i])
    fi = pd.DataFrame({'Feature': Ftrs, 'OriginalScore': Impts}, index = Ftrs)
    #fi


# DISPLAY
OriginalImportance = fi
fi['Score'] = fi.OriginalScore * 100 / (max(fi.OriginalScore))
fi



##### SAVING
Folder = 'C:/Users/gusdany/Desktop/FloraFinal/Results/VariableImportance/'
Path = Folder + Document + '.csv'
fi.to_csv(Path)
