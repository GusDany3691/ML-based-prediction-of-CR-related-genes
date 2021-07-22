# ML-based-prediction-of-CR-related-genes
Here you will find two folders: 'Data' and 'Code'.

The folder 'Data' contains other 4 folders:

 * 'Ageing and CR genes' : Contains ageing-related genes, CR-related genes (both, orthologs from 'GenDR - Gene Manipulations' and 'GenDR - Gene Expression').
 * 'Datasets'            : Contains most of the datasets created during the work. Not all of them are contained (Wholedatasets, Proteins dataset and Coexpression dataset) due to their big size. They can be, however, shared by request.
 * 'Features Importance' : Features importance from GO terms and PathDIP pathways.
 * 'CR-probabilites'     : Contains two folders: 'Strongest models' has features importances for {PathDIP, CAT} and {GO terms, BRF}. 'Nested CV outer folds' contains the outer fold data for each of the ML algorithms in GO terms and PathDIP pathways, the two most predictive datasets.

The 'Datasets' folder contains csv files of each one of the datasets used in this work.

The 'Code' folder is split  in three internal folders: 
  * 'Datasets construction' which contains R files each of wich is focused on retrieving different types of features and constructing the datasets.
  * 'Machine learning' which contains python files for conducting the ML predictions and features importance computations.
  * 'Results analysis' which contains R and python files for perfoming the statistical tests on features importance and performing the CR-probability analysis of the genes. 
  
