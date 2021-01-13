import scipy as sp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
from scvelo.preprocessing.moments import get_connectivities
import seaborn as sns
from sklearn.utils import resample
from scipy.stats import entropy


def lstsq_scoring(adata, sig_matrix, score_name='_Score', col_name='_Col', pheno_name='Phenotype_unfiltered', groupby='cline',
                  log=False,expm1 = False, smooth=True, save=True, save_name='', basis = 'umap', mouse = False):
    # score_name is also used for saving the plots
    # if gene signature is based on un-log-transformed data, use default setting log == False and adata.raw
    # else, use adata.X and log-transformed signature
    if log == True: sig_matrix = np.log1p(sig_matrix)

    print(sig_matrix.head())

    adata_ = adata.copy()
    if mouse == True:
        sig_matrix.index = [i.capitalize() for i in sig_matrix.index]
    sig_matrix = sig_matrix.loc[[i for i in set(sig_matrix.index).intersection(adata_.var_names)]]
    sig_matrix = sig_matrix / np.linalg.norm(sig_matrix, axis=0)

    adata_ = adata_[:, list(sig_matrix.index)]

    if sp.sparse.issparse(adata_.X.transpose()):
        adata_ = adata_.X.transpose().todense()
        print('csc matrix to dense')
    else:
        adata_ = adata_.X.transpose()
    print(sig_matrix.shape, adata_.shape)

    if expm1 == True:
        x, res, rank, s = np.linalg.lstsq(np.expm1(np.array(sig_matrix)), np.expm1(adata_))

    else:
        x, res, rank, s = np.linalg.lstsq(np.array(sig_matrix), adata_)

    pheno_score = pd.DataFrame(x.transpose(), index=adata.obs_names, columns=sig_matrix.columns)
    col_score = pd.DataFrame(index=adata.obs_names)
    for i in sig_matrix:
        adata.obs[str(i) + score_name] = pheno_score[i]
        col = adata.obs[str(i) + score_name]
        #         col /= np.max(np.abs(col))
        ax = scv.pl.scatter(adata, color='lightgrey', show=False, alpha=1, size=200, basis=basis)
        if save == True:
            if smooth == True:
                scv.pl.scatter(adata, color=col, basis=basis, fontsize=24, size=80, smooth=True,
                               color_map='RdBu_r', colorbar=True, title=i, ax=ax, zorder=2,
                               save=f'{save_name}_{i}{score_name}.pdf')
                c = get_connectivities(adata).dot(col)
            else:
                scv.pl.scatter(adata, color=col, basis=basis, fontsize=24, size=80,
                               color_map='RdBu_r', colorbar=True, title=i, ax=ax, zorder=2,
                               save=f'{save_name}_{i}{score_name}.pdf')
                c = col
        else:
            if smooth == True:
                scv.pl.scatter(adata, color=col, basis=basis, fontsize=24, size=80, smooth=True,
                               color_map='RdBu_r', colorbar=True, title=i, ax=ax, zorder=2)
                c = get_connectivities(adata).dot(col)
            else:
                scv.pl.scatter(adata, color=col, basis=basis, fontsize=24, size=80,
                               color_map='RdBu_r', colorbar=True, title=i, ax=ax, zorder=2)
                c = col
        adata.obs[str(i) + col_name] = c
        col_score[str(i)] = c
        if save == True:
            sc.pl.violin(adata, groupby=groupby, keys=str(i) + col_name, save=f'_{i}{score_name}.pdf')
        else:
            sc.pl.violin(adata, groupby=groupby, keys=str(i) + col_name)
            plt.hlines(y=0, xmin=-1, xmax=11)

        plt.close()

    pheno = []
    maxes = []
    for i, r in col_score.iterrows():
        pheno.append(r.idxmax())
        maxes.append(r.max())
    adata.obs[pheno_name] = pheno
    adata.obs['max'] = maxes
    return adata


def lstsq_scoring_background(adata, sig_matrix, log=False, expm1 = False,smooth=True, mouse = False):
    # score_name is also used for saving the plots
    # if gene signature is based on un-log-transformed data, use default setting log == False and adata.raw
    # else, use adata.X and log-transformed signature
    adata_shuffled = adata.copy()
    random = []
    try:
        for r in np.asarray(adata.X.todense()):
            random.append(resample(r))
    except AttributeError:
        for r in np.asarray(adata.X):
            random.append(resample(r))
    adata_shuffled.X = np.asarray(random)

    adata_ = adata_shuffled.copy()
    if mouse == True:
        sig_matrix.index = [i.capitalize() for i in sig_matrix.index]
    sig_matrix = sig_matrix.loc[[i for i in set(sig_matrix.index).intersection(adata_.var_names)]]
    adata_ = adata_[:, list(sig_matrix.index)]
    if sp.sparse.issparse(adata_.X.transpose()):
        adata_ = adata_.X.transpose().todense()
        print('csc matrix to dense')
    else:
        adata_ = adata_.X.transpose()
    if log == True:
        x, res, rank, s = np.linalg.lstsq(np.log1p(np.array(sig_matrix)), adata_)

    elif expm1 == True:
        x, res, rank, s = np.linalg.lstsq(np.expm1(np.array(sig_matrix)), np.expm1(adata_))

    else:
        x, res, rank, s = np.linalg.lstsq(np.array(sig_matrix), adata_)

    pheno_score = pd.DataFrame(x.transpose(), index=adata.obs_names, columns=sig_matrix.columns)
    col_score = pd.DataFrame(index=adata.obs_names)
    if type(get_connectivities(adata)) == type(None):
        sc.pp.neighbors(adata)
    for i in sig_matrix:
        adata.obs[str(i) + 'score'] = pheno_score[i]
        col = adata.obs[str(i) + 'score']
        #         col /= np.max(np.abs(col))
        if smooth == True:
            c = get_connectivities(adata).dot(col)
        else:
            c = col

        adata.obs[str(i) + 'col'] = c
        col_score[str(i)] = c

    pheno = []
    maxes = []
    for i, r in col_score.iterrows():
        pheno.append(r.idxmax())
        maxes.append(r.max())


    del(adata_shuffled)
    del(adata_)
    return col_score

def signature_scoring(adata, sig_matrix, method = 'lstsq',score_name='_Score', col_name='_Col', pheno_name='Phenotype_unfiltered',
                      groupby='cline', smooth=True, save=False, basis = 'umap', num_bg = 20, mouse = False, threshold = 0.05,
                      correction = None, log = False, expm1 = False, dt = 1):
    print("Scoring individual cells...")
    if method == 'lstsq':
        adata_scores = lstsq_scoring(adata, sig_matrix, score_name=score_name, col_name=col_name, mouse=mouse, basis = basis,
                                               groupby=groupby, pheno_name=pheno_name, smooth=smooth, save=save,log = log, expm1 = expm1)

        adata_t1 = adata.copy()

        adata_t1.X = adata.X + np.nan_to_num(adata.layers['velocity'])*dt
        adata_t1_scores = lstsq_scoring(adata_t1, sig_matrix, score_name=f"{score_name}_t1", col_name=col_name, mouse=mouse,
                                     basis=basis, groupby=groupby, pheno_name=pheno_name, smooth=smooth, save=save, log=log,
                                     expm1=expm1)

        for i in sig_matrix:
            adata.obs[str(i) + f"{score_name}_t1"] = adata_t1.obs[str(i) + f"{score_name}_t1"]

        # Plotting
        custom_palette = ['#fc8d62', '#66c2a5', '#FFD43B', '#8da0cb', '#e78ac3', '#a6d854']
        plt.figure(figsize=[12, 8])
        sc.pl.scatter(adata_scores, color=pheno_name, palette=custom_palette, basis=basis)
        plt.show()

        #generate background
        print("Generating background...")
        col_score_shuffled = pd.DataFrame()
        for rand in range(num_bg):
            print(rand)
            col_score = lstsq_scoring_background(adata, sig_matrix,mouse = mouse, smooth=smooth,log = log, expm1 = expm1)
            col_score_shuffled = col_score_shuffled.append(col_score, ignore_index=True)

        thresholds = dict()
        means = dict()
        sds = dict()
        for i in col_score_shuffled:
            means[i] = np.mean(col_score[i])
            sds[i] = np.std(col_score[i])
            if type(correction) == type(None):
                thresholds[i] = np.percentile(col_score[i], 100 - threshold * 100)

            elif correction == 'Bonferroni':
                thresholds[i] = np.percentile(col_score[i],
                                              100 - (threshold * 100 / len(col_score[i])))  # Bonferroni correction
            elif correction == 'Positive':
                thresholds[i] = means[i]

            sns.kdeplot(col_score_shuffled[i])
        plt.show()


    print("Correcting individual cell scores...")
    plt.rcParams["figure.figsize"] = [8, 6]

    col_score_filtered = pd.DataFrame(index=adata.obs_names)

    for i in sig_matrix:
        #sig score is any score above threshold
        #z score is standardized score for any over threshold
        # for generalist/specialist distinction, I want to standardize each score and then choose the positive ones
        adata.obs[f'{i}_Sig{score_name}'] = adata.obs[f'{i}{col_name}'] * (adata.obs[f'{i}{col_name}'] > thresholds[i])
        adata.obs[f'{i}_Z{score_name}'] = [(x - means[i]) / sds[i] for x in adata.obs[f'{i}{col_name}']] * \
                                    (adata.obs[f'{i}{col_name}'] > thresholds[i])

        ax = scv.pl.scatter(adata, color='lightgrey', show=False, alpha=1, size=200, basis=basis)
        scv.pl.scatter(adata, color=f'{i}_Z{score_name}', basis=basis, fontsize=24, size=80, vmin=0, smooth=smooth,
                       color_map='RdBu_r', colorbar=True, title=i, ax=ax, zorder=2)
        col = adata.obs[f'{i}_Z{score_name}']
        #     c = get_connectivities(adata).dot(col)
        col_score_filtered[i] = col
        sc.pl.violin(adata, groupby=groupby, keys=str(i) + f'_Z{score_name}')
        plt.close()
        plt.figure()
        ax = sc.pl.violin(adata, groupby=groupby, keys=str(i) + col_name, show=False)
        ax.hlines(y=thresholds[i], xmin=-.5, xmax=7.5, linestyle='--')
        plt.xticks(rotation=90)
        plt.show()
        plt.close()

def subtype_cells(adata,sig_matrix, basis = 'umap',score_name = 'Z_Score', plot = True, pheno_name = 'Phenotype_filtered',
                  mix_filter = 1, nonnegative = True):
    pheno = []
    maxs = []
    score_entropy = []
    score_df = pd.DataFrame(index=adata.obs_names)
    for i in sig_matrix:
        if nonnegative:
            adata.obs[f'{i}_{score_name}'][adata.obs[f'{i}_{score_name}']<0] = 0
            if f'{i}_{score_name}_t1' in adata.obs.columns:
                adata.obs[f'{i}_{score_name}_t1'][adata.obs[f'{i}_{score_name}_t1']<0] = 0

        score_df[i] = adata.obs[f'{i}_{score_name}']


    for i, r in score_df.iterrows():
        if r.sum() == 0:
            pheno.append('None')
            score_entropy.append(None)
        else:
            if (r / r.sum()).max() < mix_filter:
                pheno.append('Mixed')
            else:
                pheno.append(r.idxmax())
            score_entropy.append(entropy(r.T / r.sum()))

        maxs.append(r.max())

    adata.obs[pheno_name] = pheno
    adata.obs['max_filtered'] = maxs
    adata.obs['score_entropy'] = score_entropy
    custom_palette = ['grey', 'lightgrey', '#fc8d62', '#66c2a5', '#FFD43B', '#8da0cb', '#e78ac3']
    if plot:
        plt.figure(figsize=[10, 12])
        ax = scv.pl.scatter(adata, color='lightgrey', show=False, alpha=1, size=40, basis=basis, figsize=(8, 6))

        scv.pl.scatter(adata, color=pheno_name, palette=custom_palette, basis=basis, zorder=2, ax=ax, size=40,
                   legend_loc='on right')
        plt.show()

