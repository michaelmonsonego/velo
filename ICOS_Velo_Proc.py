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

warnings.filterwarnings('ignore')
os.chdir(r'D:\Michael\velo')

# Read in the clusters
Clusters = pd.read_csv("clusters.csv", delimiter=',', index_col=0) # which cluster every cell belongs to
Treatment = pd.read_csv("Treatment.csv", delimiter=',', index_col=0) # defined control or treatment for every cell
# Read UMAP exported R
UMAP = pd.read_csv("cell_embeddings_umap.csv", delimiter=',', index_col=0) # the umap value for each cell

# Read filtered feature bc matrix output from cellranger count
counts_C = sc.read_10x_mtx(r"ICOS_shai\antiICOS_C1\filtered_feature_bc_matrix", var_names='gene_symbols', cache=True)
counts_T = sc.read_10x_mtx(r"ICOS_shai\antiICOS_T1\filtered_feature_bc_matrix", var_names='gene_symbols', cache=True)

counts_C.obs.index = ["Control_" + bc for bc in counts_C.obs.index.tolist()]
counts_T.obs.index = ["Treatment_" + bc for bc in counts_T.obs.index.tolist()]

counts_full = anndata.concat([counts_C, counts_T])

# Filter Cells to only clustered ones
common_ind = list(set(counts_full.obs.index).intersection(Clusters.index.tolist(), UMAP.index.tolist()))
counts_sub = counts_full[common_ind]
UMAP = UMAP.loc[common_ind]
# Transform to Numpy (for formatting)
UMAP = UMAP.to_numpy()

Clusters = Clusters.loc[common_ind]
Treatment = Treatment.loc[common_ind]

# Add Clusters from Seurat to object
counts_sub.obs['Clusters'] = Clusters
counts_sub.obs['Treatment'] = Treatment
# Add UMAP to object
counts_sub.obsm["X_umap"] = UMAP

# load loom files for spliced/unspliced matrices for each sample:
ldataC = scv.read(r'ICOS_shai\antiICOS_C1\C1.loom', cache=True)
ldataT = scv.read(r'ICOS_shai\antiICOS_T1\T1.loom', cache=True)

ldataC.obs.index = [bc.replace("C1:", "Control_").replace("x", "-1") for bc in ldataC.obs.index.tolist()]
ldataT.obs.index = [bc.replace("T1:", "Treatment_").replace("x", "-1") for bc in ldataT.obs.index.tolist()]

ldataC = ldataC[np.isin(ldataC.obs.index, common_ind)]
ldataT = ldataT[np.isin(ldataT.obs.index, common_ind)]

ldataC.var_names_make_unique()
ldataT.var_names_make_unique()

ldata = anndata.concat([ldataC, ldataT])

# Merge velocyto and cellranger outputs
merged = scv.utils.merge(counts_sub, ldata)
scv.pp.filter_and_normalize(merged, flavor='seurat')
scv.pp.moments(merged, n_pcs=8, n_neighbors=30)
scv.tl.recover_dynamics(merged, n_jobs=8)  # long time
scv.tl.velocity(merged, mode='dynamical')
scv.tl.velocity_graph(merged)
scv.pl.velocity_embedding_stream(merged, color='Clusters', basis="umap", dpi=300) #, save='velocity_embedding_stream.png'

# scv.write(filename='merged.h5ad', adata=merged)
# merged = sc.read(filename='merged.h5ad')

scv.tl.recover_latent_time(merged)
scv.pl.scatter(merged, color='latent_time', color_map='gnuplot', size=80)

#M#
scv.pl.velocity_embedding(merged, color='Clusters', basis='umap', arrow_size=8, dpi=300)  # M# arrow size wont change
scv.pl.velocity_embedding_grid(merged, color='Clusters', basis='umap', arrow_size=8, dpi=300)  # M# make arrows longer
scv.pl.velocity(merged, ['Ifngr1'], dpi=300)
scv.pl.velocity(merged, ['Ifngr1'], color='Clusters', dpi=300)  # M# same but with clusters color(why does this erase two other graphs?)
scv.pl.velocity(merged, ['Ifngr1'], color='Treatment', dpi=300)  # M# same but with clusters color(why does this erase two other graphs?)

scv.pl.velocity(merged, ['Ifngr1'], color='Clusters', dpi=300, add_outline=True)  # M# looks worse
scv.pl.proportions(merged)
# scv.pl.scatter(merged, 'Cpe', color=['Clusters', 'velocity'],

#M# graph visualisation
scv.tl.paga(merged, groups='Clusters')
scv.pl.paga(merged, basis='umap', size=50, alpha=.1, min_edge_width=2, node_size_scale=1.5)

# M# driver genes :
top_genes = merged.var['fit_likelihood'].sort_values(ascending=False).index  # M# fixme : this works, just using fixme for fun
scv.pl.scatter(merged, basis=top_genes[:15], color='Clusters', ncols=5, frameon=False)
scv.pl.scatter(merged, x='latent_time', y=var_names, frameon=False) # todo: try make work, from dynamical modeling scvelo


scv.tl.rank_dynamical_genes(merged, groupby='Clusters')
df = scv.get_df(merged, 'rank_dynamical_genes/names')
df.head(5)
#M# write dynamical genes to excel
try:
    df.to_excel('dynamical_genes.xlsx')
except Exception as e:
    print("An error occurred:", e)

for cluster in merged.obs['Clusters'].unique().tolist():
    scv.pl.scatter(merged, df[cluster][:5], color='Clusters', ylabel=cluster, frameon=False)

scv.tl.rank_velocity_genes(merged, groupby='Clusters', min_corr=.3) #M# todo: whats the difference between this and rank_dynamical_genes?
df = scv.get_df(merged, 'rank_velocity_genes/names')
df.head()
for cluster in merged.obs['Clusters'].unique().tolist():
    scv.pl.scatter(merged, df[cluster][:5], ylabel=cluster, frameon=False)  # M# why does this not work when i add : color='Clusters'?

for cluster in merged.obs['Clusters'].unique().tolist():
    scv.pl.scatter(merged, df[cluster][:5], ylabel=cluster, frameon=False, color='Clusters')  # M# why does this not work when i add : color='Clusters'?


#M# subsetting from specific genes to see how this afects the velocity umap
genes_to_remove = ['Actb', 'Gm8369', 'Rtp4', 'Rgs1', 'Ifi206', 'Ptpn18', 'S100a6', 'Il2ra', 'Xaf1', 'Sdcbp2', 'Arl5a', 'Snx2']
gene_mask = np.logical_not(np.isin(merged.var_names, genes_to_remove))
sub_merged=merged[:, gene_mask]
scv.pp.filter_and_normalize(sub_merged, flavor='seurat')
scv.pp.moments(sub_merged, n_pcs=8, n_neighbors=30)
scv.tl.recover_dynamics(sub_merged, n_jobs=8)  # long time
scv.tl.velocity(sub_merged, mode='dynamical')
scv.tl.velocity_graph(sub_merged)
scv.pl.velocity_embedding_stream(sub_merged, color='Clusters', basis="umap", dpi=300, save='sub_no_Naive.png')
scv.tl.recover_latent_time(sub_merged)
scv.pl.scatter(sub_merged, color='latent_time', color_map='gnuplot', size=80)


if __name__== '__main__' :
    merged = sc.read(filename='merged.h5ad')

   #M# testing
    scv.pl.velocity(merged, ['Ifngr1'], dpi=300)
    scv.pl.velocity(merged, ['Ifngr1'], color='Clusters', dpi=300)  # M# same but with clusters color(why does this erase two other graphs?)
    scv.pl.velocity(merged, ['Ifngr1'], color='Treatment', dpi=300)  # M# same but with clusters color(why does this erase two other graphs?)
    #M# find most dynamic genes per treatment : genes whos velocities are interesting following treatment
    scv.tl.rank_velocity_genes(merged, groupby='Treatment', min_corr=.3)  # M# todo: whats the difference between this and rank_dynamical_genes?
    df_by_treatment = scv.get_df(merged, 'rank_velocity_genes/names')
    for cluster in df_by_treatment.columns:
        for gene_to_show in df_by_treatment.loc[:4, cluster]:
            scv.pl.velocity(merged, var_names=gene_to_show, color='Treatment', title=f'sub_{cluster}', dpi=300)  # , title=f'sub_{cluster}'
       # genes_to_show= df_by_treatment.loc[:4, cluster]
       # scv.pl.velocity(merged, var_names=genes_to_show, color='Treatment', dpi=300) # , title=f'sub_{cluster}'

    scv.tl.rank_velocity_genes(merged, groupby='Clusters', min_corr=.3)  # M# todo: whats the difference between this and rank_dynamical_genes?
    df_by_clusters = scv.get_df(merged, 'rank_velocity_genes/names')
    for cluster in df_by_clusters.columns:
        genes_to_remove = df_by_clusters.loc[:9, cluster]
        gene_mask = np.logical_not(np.isin(merged.var_names, genes_to_remove))
        sub_merged = merged[:, gene_mask]
        scv.pp.filter_and_normalize(sub_merged, flavor='seurat')
        scv.pp.moments(sub_merged, n_pcs=8, n_neighbors=30)
        scv.tl.recover_dynamics(sub_merged, n_jobs=8)  # long time
        scv.tl.velocity(sub_merged, mode='dynamical')
        scv.tl.velocity_graph(sub_merged)
        scv.pl.velocity_embedding_stream(sub_merged, color='Clusters', basis="umap", title=f'sub_no_{cluster}', dpi=300, save=f'sub_no_{cluster}.png')