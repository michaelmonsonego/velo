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

# sc.write(filename='merged.h5ad', adata=merged)
merged = sc.read(filename='merged.h5ad')

scv.tl.recover_latent_time(merged)
scv.pl.scatter(merged, color='latent_time', color_map='gnuplot', size=80)

#M#
scv.pl.velocity_embedding(merged, color='Clusters', basis='umap', arrow_size=8, dpi=300)  # M# arrow size wont change
scv.pl.velocity_embedding_grid(merged, color='Clusters', basis='umap', arrow_size=8, dpi=300)  # M# make arrows longer
scv.pl.velocity(merged, ['Ifngr1'], dpi=300)
scv.pl.velocity(merged, ['Ifngr1'], color='Clusters', dpi=300, save='Ifngr_clusters.png')  # M# same but with clusters color(why does this erase two other graphs?)
scv.pl.velocity(merged, ['Ifngr1'], color='Treatment', dpi=300, save='Ifngr_treatment.png')  # M# same but with clusters color(why does this erase two other graphs?)
scv.pl.scatter(merged, ['Ifngr1'], color='Clusters', dpi=300, legend_loc='best') # , save='Ifngr_clusters.png'
scv.pl.scatter(merged, ['Ifngr1'], color='Treatment', dpi=300, legend_loc='best')# , save='Ifngr_treatment.png'

scv.pl.velocity(merged, ['Ifngr1'], color='Clusters', dpi=300, add_outline=True)  # M# looks worse
scv.pl.proportions(merged)
# scv.pl.scatter(merged, 'Cpe', color=['Clusters', 'velocity'],

#M# graph visualisation
scv.tl.paga(merged, groups='Clusters')
scv.pl.paga(merged, basis='umap', size=50, alpha=.1, min_edge_width=2, node_size_scale=1.5) # , save='paga_velocity_embedding_stream.png'



# M# driver genes :
top_genes = merged.var['fit_likelihood'].sort_values(ascending=False).index  # M# fixme : this works, just using fixme for fun
scv.pl.scatter(merged, basis=top_genes[:15], color='Clusters', ncols=5, frameon=False) # , save='driver_genes_scatter.png'
scv.pl.scatter(merged, x='latent_time', y=var_names, frameon=False) # todo: try make work, from dynamical modeling scvelo


scv.tl.rank_dynamical_genes(merged, groupby='Clusters')
dynamical_genes = scv.get_df(merged, 'rank_dynamical_genes/names')
dynamical_genes.head(5)
#M# write dynamical genes to excel
try:
    dynamical_genes.to_excel('dynamical_genes.xlsx')
except Exception as e:
    print("An error occurred:", e)

for cluster in merged.obs['Clusters'].unique().tolist():
    scv.pl.scatter(merged, dynamical_genes[cluster][:5], color='Clusters', ylabel=cluster, frameon=False) # , save=f'dynamical_for_{cluster}.png'

scv.tl.rank_velocity_genes(merged, groupby='Clusters', min_corr=.3) #M# todo: whats the difference between this and rank_dynamical_genes?
velocity_genes = scv.get_df(merged, 'rank_velocity_genes/names')
velocity_genes.head()
velocity_genes.to_excel('velocity_genes.xlsx')
for cluster in merged.obs['Clusters'].unique().tolist():
    scv.pl.scatter(merged, velocity_genes[cluster][:5], ylabel=cluster, frameon=False, color='Clusters')  # M# why does this not work when i add : color='Clusters'?
#M# todo : rank_velocity_genes() and rank_dynamical_genes() does same thing? understand difference

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

#todo : make inedxes right after subsetting. getting weird warnings
#M# splitting to 2 objects by treatment
control_merged = merged[merged.obs['Treatment']=='Control'].copy()
treat_merged = merged[merged.obs['Treatment']=='Treatment'].copy()
#M# check
control_merged.obs['Clusters']
control_merged.obsm["X_umap"]
control_merged.obs.index

#M# Control analysis
scv.pp.filter_and_normalize(control_merged, flavor='seurat')
scv.pp.neighbors(control_merged)
scv.pp.moments(control_merged, n_pcs=8, n_neighbors=30)
scv.tl.recover_dynamics(control_merged, n_jobs=8)  # long time
scv.tl.velocity(control_merged, mode='dynamical')
scv.tl.velocity_graph(control_merged)
scv.pl.velocity_embedding_stream(control_merged, color='Clusters', basis="umap", dpi=300) # , save='.png'
scv.tl.paga(control_merged, groups='Clusters')
scv.pl.paga(control_merged, basis='umap', size=50, alpha=.1, min_edge_width=2, node_size_scale=1.5)
#sc.write(filename='control_merged.h5ad', adata=control_merged)

#M# Treatment analysis
treat_merged.obs['Clusters']
scv.pp.filter_and_normalize(treat_merged, flavor='seurat')
scv.pp.moments(treat_merged, n_pcs=8, n_neighbors=30)
scv.tl.recover_dynamics(treat_merged, n_jobs=8)  # long time
scv.tl.velocity(treat_merged, mode='dynamical')
scv.tl.velocity_graph(treat_merged)
scv.pl.velocity_embedding_stream(treat_merged, color='Clusters', basis="umap", dpi=300) # , save='treatment_obj.png'
scv.tl.paga(treat_merged, groups='Clusters')
scv.pl.paga(treat_merged, basis='umap', size=50, alpha=.1, min_edge_width=2, node_size_scale=1.5)
#sc.write(filename='treat_merged.h5ad', adata=treat_merged)




#M# trying functions from scv tutorial
#M# testing for differential kinetics:
merged.obs['clusters'] = merged.obs['Clusters'] #M# difference is c/C
scv.tl.differential_kinetic_test(merged)  # , groupby=None #
differential_kin_df = scv.get_df(merged, ['fit_diff_kinetics', 'fit_pval_kinetics'], precision=2)
differential_kin_df.head(10)
scv.pl.scatter(merged, basis=['Nop58', 'Tuba4a', 'Sp110', 'Sp100'], add_outline='fit_diff_kinetics')# , save='differential_kinetecs_1.png'
scv.pl.scatter(merged, basis=['Bcl2', 'Rgs2', 'Cacna1e', 'Tor3a'], add_outline='fit_diff_kinetics',legend_loc='best') #M# , save='differential_kinetecs_2.png'
#M# recompute velocities considering new knowledge about differential kinetics!
top_genes = merged.var['fit_likelihood'].sort_values(ascending=False).index[:100]
scv.tl.differential_kinetic_test(merged, var_names=top_genes, groupby='clusters')
scv.pl.scatter(merged, basis=top_genes[:15], ncols=5, add_outline='fit_diff_kinetics')
scv.pl.scatter(merged, basis=top_genes[15:30], ncols=5, add_outline='fit_diff_kinetics')
scv.tl.velocity(merged, diff_kinetics=True)
scv.tl.velocity_graph(merged)
scv.pl.velocity_embedding(merged, dpi=120, arrow_size=2, arrow_length=2)
scv.pl.velocity_embedding_stream(merged, color='Clusters', basis="umap", dpi=300, save='recomputed_velocities_stream.png') # , save='.png'

scv.tl.paga(merged, groups='Clusters')
scv.pl.paga(merged, basis='umap', size=50, alpha=.1, min_edge_width=2, node_size_scale=1.5)
scv.tl.recover_latent_time(merged)
scv.pl.scatter(merged, color='latent_time', color_map='gnuplot', size=80)

scv.pl.heatmap(merged, ['Ifngr1', 'Pycard', 'Mrpl20', 'Ly6e'])





#M# removing 100 most dynamic genes from naive like cluster to see influence on velocity stream umap
scv.tl.rank_velocity_genes(merged, groupby='Clusters', min_corr=.3)  # M# todo: whats the difference between this and rank_dynamical_genes?
df_by_clusters = scv.get_df(merged, 'rank_velocity_genes/names')
genes_to_remove = df_by_clusters.loc[:99, '3_Naive_like']
gene_mask = np.logical_not(np.isin(merged.var_names, genes_to_remove))
sub_merged_no_naive = merged[:, gene_mask]
scv.pp.filter_and_normalize(sub_merged_no_naive, flavor='seurat')
scv.pp.moments(sub_merged_no_naive, n_pcs=8, n_neighbors=30)
scv.tl.recover_dynamics(sub_merged_no_naive, n_jobs=8)  # long time
scv.tl.velocity(sub_merged_no_naive, mode='dynamical')
scv.tl.velocity_graph(sub_merged_no_naive)
scv.pl.velocity_embedding_stream(sub_merged_no_naive, color='Clusters', basis="umap", title=f'sub_no_naive', dpi=300)# , save='sub_no_naive.png'
scv.tl.paga(sub_merged_no_naive, groups='Clusters')
scv.pl.paga(sub_merged_no_naive, basis='umap', size=50, alpha=.1, min_edge_width=2, node_size_scale=1.5)#, save='paga_sub_no_naive.png'

#M# code above only calculates new velocity for same umap. doesnt look great

#M# trying to calculate velocity umap stream embedding based on most dynamic genes between treatments
#M# OR trying to calculate velocity based on genes i choose (most dynamic between treatments)
scv.tl.recover_dynamics(merged, var_names='all') #M# maybe edit var_names = [ , , , , , ] or var_names = {list of interesting genes} (same thing)
#M# came from here : https://github.com/theislab/scvelo/issues/200




#M# understanding number of velocity genes
print(merged.var['velocity_genes'].sum(), merged.n_vars) #M# 917 velocity contributing genes
#M# phase plots of top likelyhood genes
top_genes = merged.var_names[merged.var.fit_likelihood.argsort()[::-1]]
scv.pl.scatter(merged, basis=top_genes[:10], ncols=5)

#M# adjusting threshold for velocity genes
modify_velo_params = merged
scv.tl.velocity(modify_velo_params, mode='dynamical', min_r2=1e-5, min_likelihood=0.00001) # original : min_r2=1e-3, min_likelihood=0.001
scv.tl.velocity_graph(modify_velo_params)
scv.pl.velocity_embedding_stream(modify_velo_params, color='Clusters', basis="umap", title=f'sub_no_naive', dpi=300)# , save='sub_no_naive.png'
print(modify_velo_params.var['velocity_genes'].sum(), modify_velo_params.n_vars) #M# 990 velocity contributing genes
scv.tl.paga(modify_velo_params, groups='Clusters')
scv.pl.paga(modify_velo_params, basis='umap', size=50, alpha=.1, min_edge_width=2, node_size_scale=1.5)#, save='paga_sub_no_naive.png'



#M# attempt to make velocity signatures phase plots
scv.tl.rank_velocity_genes(merged, groupby='Clusters', min_corr=.3) #M# todo: whats the difference between this and rank_dynamical_genes?
df = scv.get_df(merged, 'rank_velocity_genes/names')
df.head()
df.to_excel('velocity_genes.xlsx')
trying = merged.var[merged.var['velocity_genes']==True]
velo_genes = trying.sort_values(by='fit_likelihood', ascending=False)

checkout=merged.var["velocity_score"] #M# todo : understand info


















if __name__== '__main__' :
    merged = sc.read(filename='merged.h5ad')

   #M# testing
    scv.pl.velocity(merged, ['Ifngr1'], dpi=300)
    scv.pl.velocity(merged, ['Ifngr1'], color='Clusters', dpi=300)  # M# same but with clusters color(why does this erase two other graphs?)
    scv.pl.velocity(merged, ['Ifngr1'], color='Treatment', dpi=300)  # M# same but with clusters color(why does this erase two other graphs?)
    #M# find most dynamic genes per treatment : genes whos velocities are interesting following treatment
    scv.tl.rank_velocity_genes(merged, groupby='Treatment', min_corr=.3)  # M# todo: whats the difference between this and rank_dynamical_genes?
    df_by_treatment = scv.get_df(merged, 'rank_velocity_genes/names')
    df_by_treatment.to_excel('dynamical_genes_by_treatment.xlsx')
    for cluster in df_by_treatment.columns:
        counter=1
        for gene_to_show in df_by_treatment.loc[:4, cluster]:
            scv.pl.velocity(merged, var_names=gene_to_show, color='Treatment', legend_loc='best', dpi=300)  # , title=f'sub_{cluster}' # , save=f'dynamic_{cluster}_{counter}.png'
            counter=counter+1
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
        
        
        
        
        
        
        
        