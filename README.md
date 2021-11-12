# ML-based-prediction-of-CR-related-genes
Here you will find two folders: 'Data' and 'Code'.

The folder 'Data' contains other 5 folders:

 * 'Ageing and CR genes' : Contains ageing-related genes, CR-related genes (both, orthologs from 'GenDR - Gene Manipulations' and 'GenDR - Gene Expression').
 * 'Datasets'            : Contains most of the datasets created during the work. Not all of them are contained (Wholedatasets, Proteins dataset and Coexpression dataset) due to their big size. They can be, however, sent by request to gusdany@liverpool.ac.uk.
 * 'Features Importance' : Features importance from GO terms and PathDIP pathways. Contains two folders: 'Original' with the features and their corresponding original scores (unresized to the 0-100 rank). The folder 'Scale_abd_Discussed' presents the features with their scaled values, ontology definitions, and statistical tests.
 * 'CR-probabilites'     : Contains two folders: 'Strongest models' has features importances for {PathDIP, CAT} and {GO terms, BRF}. 'NestedCV_outerfolds' contains the outer fold data for each of the ML algorithms in GO terms and PathDIP pathways, the two most predictive datasets.
 * 'External'            : Contains the data from external databases that was used to perform this work, one folder per source: 'GeneFriends', 'Gtex', 'HAGR', 'OMA', 'PathDIP', 'PPI'.


The 'Code' folder is split  in three internal folders: 
  * 'Datasets construction' which contains R files each of which is focused on retrieving different types of features and constructing the datasets.
  * 'Machine learning' which contains python files for conducting the ML predictions and features importance computations.
  * 'Results analysis' which contains R and python files for perfoming the statistical tests on features importance and performing the CR-probability analysis of the genes. 
  
