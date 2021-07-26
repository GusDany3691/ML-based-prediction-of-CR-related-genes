import pandas as pd
import numpy  as np

CrGenes = pd.read_csv('/ML-based-prediction-of-CR-related-genes/Data/Ageing and CR genes/CRgenes.csv')

GoDataset = DatasetAdapation(pd.read_csv('/ML-based-prediction-of-CR-related-genesb/Data/Datasets/GoDataset.csv'))
sKeggPathDataset = DatasetAdapation(pd.read_csv('/ML-based-prediction-of-CR-related-genes/Data/Datasets/KeggPaths.csv')) 
KeggGeneDataset = DatasetAdapation(pd.read_csv('/ML-based-prediction-of-CR-related-genes/Data/Datasets/KeggGenes.csv'), FirstCol = "X") 
ProteinsDataset = DatasetAdapation(pd.read_csv('/ML-based-prediction-of-CR-related-genes/Data/Datasets/ProteinsDataset.csv')) 
PPImeasuresDataset = DatasetAdapation(pd.read_csv('/ML-based-prediction-of-CR-related-genes/Data/Datasets/PPImeasuresDataset.csv')) 
PPIadjencyDataset = DatasetAdapation(pd.read_csv('/ML-based-prediction-of-CR-related-genes/Data/Datasets/PPIadjencyDataset.csv')) 
PPIfilteredAdjencyDataSet = DatasetAdapation(pd.read_csv('/ML-based-prediction-of-CR-related-genes/Data/Datasets/PPIAdjencyDataset.csv')) 
PathDipDataset = DatasetAdapation(pd.read_csv('/ML-based-prediction-of-CR-related-genes/Data/Datasets/PathdipDataset.csv')) 
GtexDataset = DatasetAdapation(pd.read_csv('/ML-based-prediction-of-CR-related-genes/Data/Datasets/GtexDataset.csv'))                                   # COMPLET
CoexpressionDataset = DatasetAdapation(pd.read_csv('/ML-based-prediction-of-CR-related-genes/Data/Datasets/CoexpressionDataset.csv'))  

FinalClassFrame

WholeDataset  = pd.concat([GoDataset, sKeggPathDataset, KeggGeneDataset, PPImeasuresDataset, PPIfilteredAdjencyDataSet, PathDipDataset, GtexDataset, CoexpressionDataset], axis=1, sort=False)
AgeingGenes = WholeDataset.index
LogicCR = Dataframe(AgeingGenes.isin(CrGenes['CrGenes']))
CrArray = np.where(LogicCR, 'CR','NotCR')
ClassFrame = pd.DataFrame({'Class':CrArray}, index = AgeingGenes)

WholeDataset  = pd.concat([WholeDataset, ClassFrame], axis=1, sort=False)
                
WholeDatasetIntersection = WholeDataset.dropna(thresh=1)


WholeDataset.to_csv('/ML-based-prediction-of-CR-related-genes/Data/Datasets/WholeDataset.csv')
WholeDatasetIntersection.to_csv('/ML-based-prediction-of-CR-related-genes/Data/Datasets/WholeDatasetIntersection.csv')

#############################################################################################################################################################################################
##### FUNCTIONS #############################################################################################################################################################################
#############################################################################################################################################################################################

def DatasetAdapation(Dataset, FirstCol = "Unnamed: 0"):
    Dataset = Dataset.rename(columns={FirstCol: "Genes"})
    Dataset.index = Dataset.iloc[:,0]            # Set the first column as row names 
    Dataset = Dataset.drop(columns=['Genes', 'Class'])    # Drop first column
    return(Dataset)

# INTERSECTION
def intersection(List1, List2):
    Set1 = set(List1)
    Set2 = set(List2)
    return Set1.intersection(Set2)

# Ifelse
def if_this_else_that(x, list_of_checks, yes_label, no_label):
    if x in list_of_checks:
        res = yes_label
    else: 
        res = no_label
    return(res)