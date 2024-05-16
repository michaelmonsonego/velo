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

ldata = anndata.concat([ldata_a, ldata_b, ldata_c, ldata_d, ldata_e])
ldata.var_names_make_unique()
merged = scv.utils.merge(counts_sub, ldata)

scv.pp.filter_and_normalize(merged, flavor='seurat')
scv.pp.moments(merged, n_pcs=8, n_neighbors=30)
scv.tl.recover_dynamics(merged, n_jobs=8)  # long time
scv.tl.velocity(merged, mode='dynamical')
scv.tl.velocity_graph(merged)
scv.pl.velocity_embedding_stream(merged, color='Clusters', basis="umap", dpi=300) #, save='velocity_embedding_stream.png'

# sc.write(filename='merged.h5ad', adata=merged)
merged = sc.read(filename='merged.h5ad')

#M# velo signatures
#M# 16.5 16:40 - which signature to run? ask shai

#merged = sc.read(filename='merged.h5ad')
#M# hakllmark infg respone geneset
up_sig = pd.read_table(r'D:/user/Downloads/Ifng hallmark sig.txt')
gene_names = up_sig['HALLMARK_INTERFERON_GAMMA_RESPONSE'].tolist()
gene_names.remove('> Genes up-regulated in response to IFNG [GeneID=3458].')
formatted_gene_names = [name.capitalize() for name in gene_names]
#up_sig = up_sig['x'].str.lower().values #M# genes are with capital first letter
su.run_exp_signature_on_obj_val(merged, up_sig=formatted_gene_names)
sc.pl.umap(merged, color=["SigScore"],color_map="magma", save='_exp_sig_hallmark_infgr.png')
su.run_velo_signature_on_obj_val(merged, up_sig=formatted_gene_names)
sc.pl.umap(merged, color=["SigScore"],color_map="magma", save='_velo_sig_hallmark_infgr.png')































