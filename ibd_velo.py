import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import igraph
import scvelo as scv
import loompy as lmp
import anndata
import os
# import hdf5plugin
import warnings
import signature_utils as su  #todo : what ?

warnings.filterwarnings('ignore')
os.chdir(r'D:\Michael\velo')


# Read in the clusters
Clusters = pd.read_csv("clusters.csv", delimiter=',', index_col=0) # which cluster every cell belongs to
# Treatment = pd.read_csv("Treatment.csv", delimiter=',', index_col=0) # defined control or treatment for every cell
# Read UMAP exported R
UMAP = pd.read_csv("cell_embeddings_umap.csv", delimiter=',', index_col=0) # the umap value for each cell

# Read filtered feature bc matrix output from cellranger count
counts_a = sc.read_10x_mtx(r"<PATH>\filtered_feature_bc_matrix", var_names='gene_symbols', cache=True)
counts_b = sc.read_10x_mtx(r"<PATH>\filtered_feature_bc_matrix", var_names='gene_symbols', cache=True)
counts_c = sc.read_10x_mtx(r"<PATH>\filtered_feature_bc_matrix", var_names='gene_symbols', cache=True)
counts_d = sc.read_10x_mtx(r"<PATH>\filtered_feature_bc_matrix", var_names='gene_symbols', cache=True)
counts_e = sc.read_10x_mtx(r"<PATH>\filtered_feature_bc_matrix", var_names='gene_symbols', cache=True)


 #M# need this ? : counts.obs.index = ["Control_" + bc for bc in counts_C.obs.index.tolist()]
counts_full = anndata.concat([counts_a, counts_b, counts_c, counts_d, counts_e])

common_ind = list(set(counts_full.obs.index).intersection(Clusters.index.tolist(), UMAP.index.tolist()))
counts_sub = counts_full[common_ind]
UMAP = UMAP.loc[common_ind]
# Transform to Numpy (for formatting)
UMAP = UMAP.to_numpy()

Clusters = Clusters.loc[common_ind]

counts_sub.obs['Clusters'] = Clusters
# Add UMAP to object
counts_sub.obsm["X_umap"] = UMAP

ldata_a = scv.read(r'.loom', cache=True)
ldata_b = scv.read(r'.loom', cache=True)
ldata_c = scv.read(r'.loom', cache=True)
ldata_d = scv.read(r'.loom', cache=True)
ldata_e = scv.read(r'.loom', cache=True)
































