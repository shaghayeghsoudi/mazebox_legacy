'''
Code pulled from ALL HUMAN SAMPLES.ipynb on June 22, 2021
Sarah Maddox Groves
'''

## Classifier based on high confidence cell labels: Label Spreading
import numpy as np
from sklearn import datasets
from sklearn.semi_supervised import LabelSpreading
def label_spreading_recipe(adata, subsample = True, random_state = 0, fraction = 0.1, genes_to_use = None, kernel = 'knn'):
    '''
    genes_to_use: can be "all","hvg", or a list
    '''
    if subsample:
        adata_sub = sc.pp.subsample(adata, fraction=fraction, copy = True, random_state=random_state)
    else: adata_sub = adata.copy()
    X = pd.DataFrame(adata_sub.X.todense(), index=adata_sub.obs_names, columns = adata_sub.var_names)

    #Depending on which features are being used for label propagation, subset the data 
    if type(genes_to_use) not in [list,str]: 
        print("Warning: Genes_to_use must be str (all, hvg) or list")
        pass
    if type(genes_to_use) == type(list):
        X = X[genes_to_use]
    else:
        if genes_to_use == 'all': pass
        elif genes_to_use == 'hvg':
            sc.pp.highly_variable_genes(adata_sub, inplace=True)
            X = X.T[adata_sub.var['highly_variable']].T
    y = np.array(adata_sub.obs['Phenotype'])
    #remove Generalist, None, and nan labels and use -1 for classifier
    unlabeled = adata_sub.obs['Phenotype'] == 'Generalist'
    y[unlabeled] = -1
    y[[type(i) == type(np.nan) for i in y]] = -1
    y[[i == 'None' for i in y]] = -1
    map_dict = {}
    
    for x, i in enumerate(set(y).difference(set({-1}))):
            map_dict[i] = x
    map_dict[-1] = -1
    
    print("Fitting model...")
    #fit Label Spreading model to data
    LP = LabelSpreading(kernel = kernel)
    LP.fit(X, [map_dict[i] for i in y])
    
    
    inv_dict = {}
    for k in map_dict.keys():
        inv_dict[map_dict[k]] = k
    print(f"Adding attributes 'Phenotype_labelspread' and 'LP_prob_x' for {set(y).difference(set({-1}))}'")
    #add attributes to adata_sub
    adata_sub.obs['Phenotype_labelspread'] = [inv_dict[i] for i in LP.transduction_]
    LP_prob = pd.DataFrame(LP.label_distributions_, index = X.index)
    for c in LP_prob.columns:
        adata_sub.obs[f'LP_prob_{inv_dict[c]}'] = LP_prob[c]
    return adata_sub
adata_sub = label_spreading_recipe(adata,genes_to_use=list(sig_matrix2.index))
scv.pl.scatter(adata_sub, basis = 'umap_sc', color = [f"LP_prob_{x}" for x in ['SCLC-N', 'SCLC-A', 'SCLC-P', 'SCLC-Y', 'SCLC-A2']], figsize = (10,10), cmap = 'viridis')

scv.pl.scatter(adata_sub, basis = 'umap_sc', color = ['Phenotype','Phenotype_labelspread'], figsize = (10,10), palette=pal)
scv.pl.scatter(adata_sub, basis = 'umap_sc', color = 'batch', figsize = (10,10))
import plotly.express as px
df = px.data.election()
fig = px.scatter_ternary(adata_sub.obs, a="LP_prob_SCLC-P", b="LP_prob_SCLC-N", c="LP_prob_SCLC-Y", 
                         color = 'batch')
fig.show()
adata.obs['SCLC-NE_Score_pos'] = adata.obs['SCLC-A_Score_pos']+adata.obs['SCLC-A2_Score_pos']+adata.obs['SCLC-N_Score_pos']
import plotly.express as px
fig = px.scatter_ternary(adata.obs, a="SCLC-P_Score_pos", b="SCLC-NE_Score_pos", c="SCLC-Y_Score_pos", color = 'batch')
fig.show()
import plotly.express as px
fig = px.scatter_ternary(adata.obs, a="SCLC-A2_Score_pos", b="SCLC-N_Score_pos", c="SCLC-A_Score_pos", color = 'batch')
fig.show()


##Self Training Classifier based on SVM
from sklearn.semi_supervised import SelfTrainingClassifier
from sklearn import svm 

def semi_svc_recipe(adata, subsample = True, random_state = 0, fraction = 0.1, genes_to_use = None, kernel = 'knn'):
    '''
    genes_to_use: can be "all","hvg", or a list
    '''
    if subsample:
        adata_sub = sc.pp.subsample(adata, fraction=fraction, copy = True, random_state=random_state)
    else: adata_sub = adata.copy()
    X = pd.DataFrame(adata_sub.X.todense(), index=adata_sub.obs_names, columns = adata_sub.var_names)

    #Depending on which features are being used for label propagation, subset the data 
    if type(genes_to_use) not in [list,str]: 
        print("Warning: Genes_to_use must be str (all, hvg) or list")
        pass
    if type(genes_to_use) == type(list):
        X = X[genes_to_use]
    else:
        if genes_to_use == 'all': pass
        elif genes_to_use == 'hvg':
            if 'highly_variable' not in adata_sub.var.columns:
                sc.pp.highly_variable_genes(adata_sub, inplace=True)
            X = X.T[adata_sub.var['highly_variable']].T
    y = np.array(adata_sub.obs['Phenotype'])
    low_confidence = adata_sub.obs[['SCLC-A_Score', 'SCLC-A2_Score', 'SCLC-N_Score', 'SCLC-P_Score', 'SCLC-Y_Score']].max(axis = 1) < .6
    print("Length of low confidence cells: ", sum(low_confidence))
    print("Length of remaining labeled cells: ", y.shape[0]-sum(low_confidence))
    y[low_confidence] = -1
    #remove Generalist, None, and nan labels and use -1 for classifier
    unlabeled = adata_sub.obs['Phenotype'] == 'Generalist'
    y[unlabeled] = -1
    y[[type(i) == type(np.nan) for i in y]] = -1
    y[[i == 'None' for i in y]] = -1
    
    map_dict = {}
    
    for x, i in enumerate(set(y).difference(set({-1}))):
            map_dict[i] = x
    map_dict[-1] = -1
    
 
    print(np.unique([map_dict[i] for i in y], return_counts = True))
    print("Fitting model...")
    #fit SVM Classifier model to data
    svc = svm.SVC(probability=True)
    self_training_model = SelfTrainingClassifier(svc)
    self_training_model.fit(X, [map_dict[i] for i in y])
    
    inv_dict = {}
    for k in map_dict.keys():
        inv_dict[map_dict[k]] = k
    print(inv_dict)
  
    #add attributes to adata_sub
    print("Adding attributes to adata...")
    print("Phenotype_semi_svc")

    print("semi_svm_proba_x for x in: ", [inv_dict[c] for c in self_training_model.classes_])

    svm_proba = pd.DataFrame(self_training_model.predict_proba(X), columns = self_training_model.classes_, index = X.index)
    adata_sub.obs['Phenotype_semi_svc'] = svm_proba.idxmax(axis = 1)
    for c in svm_proba.columns:
        adata_sub.obs[f'semi_svm_proba_{inv_dict[c]}'] = svm_proba[c]
    return adata_sub, svm_proba
adata_sub,svm_proba = semi_svc_recipe(adata, subsample=True,genes_to_use='hvg')
fig = px.scatter_ternary(adata_sub.obs, a="semi_svm_proba_SCLC-A", b="semi_svm_proba_SCLC-A2", c="semi_svm_proba_SCLC-Y", color = 'Phenotype')
fig.show()
scv.pl.scatter(adata_sub, basis = 'umap_sc', color = [f"semi_svm_proba_{x}" for x in ['SCLC-N', 'SCLC-A', 'SCLC-P', 'SCLC-Y', 'SCLC-A2']], figsize = (10,10), 
               cmap = 'viridis')