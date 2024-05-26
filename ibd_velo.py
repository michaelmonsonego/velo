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
import signature_utils_velo as su  #todo : what ?

warnings.filterwarnings('ignore')
os.chdir(r'D:\Michael\ibd_velo')

# Read in the clusters
Clusters = pd.read_csv("clusters.csv", delimiter=',', index_col=0) # which cluster every cell belongs to
# Treatment = pd.read_csv("Treatment.csv", delimiter=',', index_col=0) # defined control or treatment for every cell
# Read UMAP exported R
UMAP = pd.read_csv("cell_embeddings_umap.csv", delimiter=',', index_col=0) # the umap value for each cell

# Read filtered feature bc matrix output from cellranger count
counts_a = sc.read_10x_mtx(r"D:\Michael\IBD_Project\A_exp75\outs\filtered_feature_bc_matrix", var_names='gene_symbols', cache=False)
counts_b = sc.read_10x_mtx(r"D:\Michael\IBD_Project\D_exp75\outs\filtered_feature_bc_matrix", var_names='gene_symbols', cache=False)
counts_c = sc.read_10x_mtx(r"D:\Michael\IBD_Project\F_exp75\outs\filtered_feature_bc_matrix", var_names='gene_symbols', cache=False)
counts_d = sc.read_10x_mtx(r"D:\Michael\IBD_Project\G_exp75\outs\filtered_feature_bc_matrix", var_names='gene_symbols', cache=False)
counts_e = sc.read_10x_mtx(r"D:\Michael\IBD_Project\H_exp75\outs\filtered_feature_bc_matrix", var_names='gene_symbols', cache=False)


counts_full = anndata.concat([counts_a, counts_b, counts_c, counts_d, counts_e])
#M# unifying format
Clusters.index = [name.split('_')[-2] for name in Clusters.index.tolist()]
Clusters = Clusters[~Clusters.index.duplicated(keep='first')]
UMAP.index = [name.split('_')[-2] for name in UMAP.index.tolist()]
UMAP = UMAP[~UMAP.index.duplicated(keep='first')]
#M# filter out cells who were filtered out in preprocessing
common_ind = list(set(counts_full.obs.index).intersection(Clusters.index, UMAP.index))
counts_sub = counts_full[counts_full.obs_names.isin(common_ind), :]
counts_sub=counts_sub[~counts_sub.obs.index.duplicated(), :]

#M# debug :
Clusters.index.equals(UMAP.index) # True
c1 = list(set(Clusters.index).intersection(UMAP.index)) # why is this only 6530 and not 6532?

#M# doesnt work : UMAP = UMAP.loc[common_ind]
#M# doesnt work : Clusters = Clusters.loc[common_ind]
# Transform to Numpy (for formatting)
UMAP = UMAP.to_numpy()

counts_sub.obs['Clusters'] = Clusters
# Add UMAP to object
counts_sub.obsm["X_umap"] = UMAP

ldata_a = scv.read(r'D:\Michael\ibd_looms\A_exp75.loom', cache=False)
ldata_b = scv.read(r'D:\Michael\ibd_looms\D_exp75.loom', cache=False)
ldata_c = scv.read(r'D:\Michael\ibd_looms\F_exp75.loom', cache=False)
ldata_d = scv.read(r'D:\Michael\ibd_looms\G_exp75.loom', cache=False)
ldata_e = scv.read(r'D:\Michael\ibd_looms\H_exp75.loom', cache=False)

#M# todo : why do i need to run this twice in order for both replacements to happen ?
ldata_a.obs.index = [bc.replace("A_e-1p75:", "") for bc in ldata_a.obs.index.tolist()]
ldata_a.obs.index = [bc.replace("x", "-1") for bc in ldata_a.obs.index.tolist()]
ldata_b.obs.index = [bc.replace("D_e-1p75:", "") for bc in ldata_b.obs.index.tolist()]
ldata_b.obs.index = [bc.replace("x", "-1") for bc in ldata_b.obs.index.tolist()]
ldata_c.obs.index = [bc.replace("F_e-1p75:", "") for bc in ldata_c.obs.index.tolist()]
ldata_c.obs.index = [bc.replace("x", "-1") for bc in ldata_c.obs.index.tolist()]
ldata_d.obs.index = [bc.replace("G_e-1p75:", "") for bc in ldata_d.obs.index.tolist()]
ldata_d.obs.index = [bc.replace("x", "-1") for bc in ldata_d.obs.index.tolist()]
ldata_e.obs.index = [bc.replace("H_e-1p75:", "") for bc in ldata_e.obs.index.tolist()]
ldata_e.obs.index = [bc.replace("x", "-1") for bc in ldata_e.obs.index.tolist()]

ldata_a = ldata_a[np.isin(ldata_a.obs.index, common_ind)].copy()
ldata_a.var_names_make_unique()
ldata_b = ldata_b[np.isin(ldata_b.obs.index, common_ind)].copy()
ldata_b.var_names_make_unique()
ldata_c = ldata_c[np.isin(ldata_c.obs.index, common_ind)].copy()
ldata_c.var_names_make_unique()
ldata_d = ldata_d[np.isin(ldata_d.obs.index, common_ind)].copy()
ldata_d.var_names_make_unique()
ldata_e = ldata_e[np.isin(ldata_e.obs.index, common_ind)].copy()
ldata_e.var_names_make_unique()


ldata = anndata.concat([ldata_a, ldata_b, ldata_c, ldata_d, ldata_e])
ldata.var_names_make_unique()

ldata=ldata[~ldata.obs.index.duplicated(), :]

merged = scv.utils.merge(counts_sub, ldata)

scv.pp.filter_and_normalize(merged, flavor='seurat')
scv.pp.moments(merged, n_pcs=8, n_neighbors=30)
scv.tl.recover_dynamics(merged, n_jobs=8)  # long time
scv.tl.velocity(merged, mode='dynamical')
scv.tl.velocity_graph(merged)
scv.pl.velocity_embedding_stream(merged, color='Clusters', basis="umap", dpi=300) #, save='velocity_embedding_stream.png'

# sc.write(filename='merged.h5ad', adata=merged)
# merged = sc.read(filename='merged.h5ad')

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































