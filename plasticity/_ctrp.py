import os.path as op
import scipy as sp
from scipy.spatial import distance
import pandas as pd
import scanpy as sc
from scvelo import settings
from scvelo import logging as logg
from scvelo.tools.utils import scale, groups_to_bool, strings_to_categoricals
from random import seed
from scipy.sparse import linalg, csr_matrix, issparse
import numpy as np
from matplotlib.colors import ListedColormap, to_rgba
import scvelo as scv
from scvelo.preprocessing.neighbors import get_connectivities
from scipy.sparse import csr_matrix

def write_to_obs(adata, key, vals, cell_subset=None):
    if cell_subset is None:
        adata.obs[key] = vals
    else:
        vals_all = adata.obs[key].copy() if key in adata.obs.keys() else np.zeros(adata.n_obs)
        vals_all[cell_subset] = vals
        adata.obs[key] = vals_all

# Would like to add in a check if the user wants to separate by groups. If some of the groups do not have an endpoint,
# suggest combining groups to find transitions across groups.
def ctrp(data, vkey='velocity', groupby=None, groups=None, self_transitions=False, basis=None,
                    weight_diffusion=0, scale_diffusion=1, eps=1e-2, copy=False, adata_dist = 'X_pca', random_seed = 1):
    '''
    :param data: AnnData object
    :param vkey: velocity key to use for transition matrix
    :param groupby: if string, analysis will be run on each category separately
    :param groups: groups for cell subsetting
    :param self_transitions:
    :param basis: basis for transitions to be restricted to. None means that the transition matrix will contain all transitions,
    'pca' means only transitions within the pca reduced space are found, etc.
    :param weight_diffusion: weight of the diffusion kernel for transition matrix calculation
    :param scale_diffusion: scale of diffusion kernel for transition matrix calculation
    :param eps: tolerance for eigenvalue calculation
    :param copy: if True, return copy of AnnData object with updated parameters
    :param adata_dist: data to be used for distance calculation. Default X_pca means distances between cells will be
    calculated in PCA space
    :param random_seed: Default 1
    :return: adata if copy == True, else None and adata is updated with the obs attributes:
    ctrp: plasticity value for each cell
    absorbing: Boolean value detailing if cell is an end state
    '''

    seed(random_seed)
    adata = data.copy() if copy else data
    logg.info('computing terminal states', r=True)

    strings_to_categoricals(adata)
    if groupby is not None:
        logg.warn(
            "Only set groupby, when you have evident distinct clusters/lineages,"
            " each with an own root and end point."
        )

    categories = [None]
    if groupby is not None and groups is None:
        categories = adata.obs[groupby].cat.categories

    # groupby = 'cell_fate' if groupby is None and 'cell_fate' in adata.obs.keys() else groupby
    # categories = adata.obs[groupby].cat.categories if groupby is not None and groups is None else [None]
    for c, cat in enumerate(categories):
        groups = cat if cat is not None else groups
        cell_subset = groups_to_bool(adata, groups=groups, groupby=groupby)
        _adata = adata if groups is None else adata[cell_subset]
        connectivities = get_connectivities(_adata, 'distances')

        T = scv.tl.transition_matrix(_adata, vkey=vkey, basis=basis, weight_diffusion=weight_diffusion,
                              scale_diffusion=scale_diffusion, self_transitions=self_transitions, backward=False)

        eigvecs_ends = eigs(T, eps=eps, perc=[2, 98])[1]
        ends = csr_matrix.dot(connectivities, eigvecs_ends).sum(1)
        # ends = eigvecs_ends.sum(1)

        ends = scale(np.clip(ends, 0, np.percentile(ends, 98)))
        write_to_obs(adata, 'end_points', ends, cell_subset)

        _adata.obs['col'] = ends

        n_ends = eigvecs_ends.shape[1]
        groups_str = ' (' + groups + ')' if isinstance(groups, str) else ''

        logg.info('    identified ' + str(n_ends) + ' end points' + groups_str)
        print("Dropping absorbing rows for fundamental matrix...")


        absorbing = np.where(ends >=1. - eps)[0]
        T = dropcols_coo(T,np.where(ends >=1. - eps)[0])
        if np.all(T.getnnz(1) == 0) == False:
            dropped = absorbing
        else:
            dropped = np.concatenate((absorbing, np.where(T.getnnz(1) == 0)[0]))
            print("Dropping rows with all zeroes...")
            T = dropcols_coo(T, np.where(T.getnnz(1) == 0)[0])

        I = np.eye(T.shape[0])
        print(T.shape)

        # TODO: still under development. If det(I-T) = 0, inspect transition matrix.
        if np.abs(np.linalg.det(I-T)) == 0:
            print('Determinant equals 0. Solving for fate using pseudoinverse.')
            fate = np.linalg.pinv(I-T)
        else:
            print("Calculating fate...")
            fate = np.linalg.inv(I - T)  # fundamental matrix of markov chain
        if issparse(T): fate = fate.A
        fate = fate / fate.sum(axis=1)  # each row is a different starting cell

        if adata_dist is None:
            try:
                adataX = dropcols_coo(_adata.X, dropped)
                adataX = sp.sparse.csr_matrix.todense(np.expm1(adataX))
            except AttributeError:
                adataX = _adata.X
                mask = np.ones(len(adataX), dtype=bool)
                mask[dropped] = False
                adataX = adataX[mask,:]

        elif adata_dist == 'raw':
            adataX = dropcols_coo(_adata.raw.X, dropped)

            adataX = sp.sparse.csr_matrix.todense(adataX)
        else:
            adataX = _adata.obsm[adata_dist]
            mask = np.ones(len(adataX), dtype=bool)
            mask[dropped] = False
            adataX = adataX[mask,:]
        print("Calculating distances...")
        dist = distance.cdist(adataX, adataX, 'euclidean')
        print("Calculating inner product...")
        ip = np.einsum('ij,ij->i', fate, dist)
        if cell_subset is not None:
            cs = cell_subset.copy()
            new_subset = pd.Series(cs, index = adata.obs_names.values)
            true_subset = np.where(cs == True)[0]
            for i in list(dropped):
                new_subset.iloc[true_subset[i]] = False
        else:
            new_subset = pd.Series([True]*len(adata.obs_names.values), index = adata.obs_names.values)
            for i in list(dropped):
                new_subset.iloc[i] = False
            if 'ctrp' in adata.obs.keys():
                adata.obs['ctrp'] = np.zeros(adata.n_obs)
        write_to_obs(adata, 'ctrp', np.log1p(ip), new_subset)

    adata.obs['absorbing'] = [i for i in adata.obs['end_points'] >=1. - eps]

    return adata if copy else None

def eigs(T, k=100, eps=1e-2, perc=None):
    try:
        eigvals, eigvecs = linalg.eigs(T.T, k=k, which='LR')  # find k eigs with largest real part

        p = np.argsort(eigvals)[::-1]                        # sort in descending order of eigenvalues
        eigvals = eigvals.real[p]
        eigvecs = eigvecs.real[:, p]

        idx = (eigvals >= 1 - eps)                           # select eigenvectors with eigenvalue of 1 within eps tolerance
        eigvals = eigvals[idx]
        eigvecs = np.absolute(eigvecs[:, idx])
        print("Eigenvalues: ", eigvals)
        print(eigvecs.shape)
        if perc is not None:
            lbs, ubs = np.percentile(eigvecs, perc, axis=0)
            eigvecs[eigvecs < lbs] = 0
            eigvecs = np.clip(eigvecs, 0, ubs)
            eigvecs /= eigvecs.max(0)

    except:
        eigvals, eigvecs = np.empty(0), np.zeros(shape=(T.shape[0], 0))
    return eigvals, eigvecs

def dropcols_coo(M, idx_to_drop):
    idx_to_drop = np.unique(idx_to_drop)
    C = M.tocoo()
    keep = ~np.in1d(C.col, idx_to_drop)
    C.data, C.row, C.col = C.data[keep], C.row[keep], C.col[keep]
    C.col -= idx_to_drop.searchsorted(C.col)    # decrement column indices
    C.row -= idx_to_drop.searchsorted(C.row)    # decrement column indices
    C._shape = (C.shape[0]- len(idx_to_drop), C.shape[1] - len(idx_to_drop))
    return C.tocsr()
