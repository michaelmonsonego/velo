import pandas as pd
import numpy as np
from scipy.stats import ranksums, f_oneway
import pingouin as pg 
import scanpy as sc
from sklearn.preprocessing import StandardScaler
import scipy


alpha = 0.9
epsilon = 0.001

# ---------------------------
# post hock analysis for violin plot
# when there are more than 2 groups
# ---------------------------
def post_hock_anaylsis(value_dict):
  dfs = []
  for key,value in value_dict.items():
    df = pd.DataFrame({"score value":value.values.reshape(-1)},index=list(range(len(value))))
    df.insert(loc = 1,column = "ident",value = [key for _ in range(len(value))])
    dfs.append(df)
  
  result = pd.concat(dfs,axis=0)
  posthocs = pg.pairwise_tukey(data = result,dv = "score value",between="ident")
  posthocs = posthocs.loc[posthocs["p-tukey"] <0.05,:]
  return posthocs

def run_f_test(value_dict):
  vecs = value_dict.values()
  vecs = [vec.values.reshape(-1) for vec in vecs]
  return f_oneway(*vecs)[1]
  
# ---------------------------
# wilcoxon test
# gene_exp - reduced only to up signature
# backround_exp
#   if there is down sig, calculates for it
#   else, the rest of the genes in expression not in up_sig
#s
#   returns (-log(rank)) -> if not significant returns 0
# ---------------------------
def wilcoxon_enrcment_test(up_sig, down_sig, exp):
    gene_exp = exp.loc[exp.index.isin(up_sig)]
    if down_sig is None:     
        backround_exp = exp.loc[~exp.index.isin(up_sig)]
    else:
        backround_exp = exp.loc[exp.index.isin(down_sig)]
        
    rank = ranksums(backround_exp,gene_exp,alternative="less")[1] # rank expression of up sig higher than backround
    rank = 1 if rank > 0.05 else rank 
    #print(rank)
    return -1 * np.log(rank)


# ---------------------------
# calculates the signature of the data
#
# returns scores vector of signature calculated per cell
# ---------------------------
def signature_values(exp, up_sig, format_flag, down_sig=None):
    up_sig = pd.DataFrame(up_sig).squeeze()
    # first letter of gene in upper case
 #   if format_flag:
#        up_sig = up_sig.apply(lambda x: x[0].upper() + x[1:].lower())
    # keep genes in sig that appear in exp data
    up_sig = up_sig[up_sig.isin(exp.index)]

    if down_sig is not None:
        down_sig = pd.DataFrame(down_sig).squeeze()
        if format_flag:
            up_sig = up_sig.apply(lambda x: x[0].upper() + x[1:].lower())
        down_sig = down_sig[down_sig.isin(exp.index)]
    
    return exp.apply(lambda cell: wilcoxon_enrcment_test(up_sig, down_sig, cell), axis=0)


# ---------------------------
# Y - scores vector of cells
# W - Adjacency matrix
#
# f_t = alpha * (W * f_(t-1)) + (1-alpha)*Y
#
# returns f/f1
# ---------------------------
def propagation(Y, W):
    W = normW(W)
    f = np.array(Y)
    Y = np.array(Y)
    W = np.array(W.values)
    
    Y1 = np.ones(Y.shape, dtype=np.float64)
    f1 = np.ones(Y.shape, dtype=np.float64)
    flag = True

    while(flag):
        next_f = alpha*(W@f) + (1-alpha)*Y
        next_f1 = alpha*(W@f1) + (1-alpha)*Y1
    
        if np.linalg.norm(next_f - f) <= epsilon and np.linalg.norm(next_f1 - f1) <= epsilon:
            flag = False
        else:
            #print(np.linalg.norm(next_f - f))
            f = next_f
            f1 = next_f1
            
    return np.array(f/f1) 


# ---------------------------
# W - Adjacency matrix
# calculated for W:
# W'_ij = W_ij/sqrt(D_iD_j)
# ---------------------------
def normW(W):
    sum_rows = pd.DataFrame(W.sum(axis=1))
    sum_rows = sum_rows @ sum_rows.T
    sum_rows **= 1/2
    return W / sum_rows


def run_exp_signature_on_obj_val(obj, up_sig, down_sig=None, conn_flag =  True, umap_flag = True):
   # scaler = StandardScaler()
   # exp = scaler.fit_transform(exp.T).T      #
    exp = obj.to_df().T
    sigs_scores = signature_values(exp, up_sig, down_sig) # wilcoxon score
    if conn_flag:

        graph = obj.obsp["connectivities"].toarray()
    else:
        graph = obj.obsp["distances"].toarray()
    #sigs_scores = propagation(sigs_scores,graph)
    obj.obs["SigScore"] = sigs_scores
    if umap_flag:
        sc.pl.umap(obj, color=["SigScore"],color_map="magma")
    else:
        sc.pl.tsne(obj, color=["SigScore"],color_map="magma")

    return obj


def run_velo_signature_on_obj_val(obj, up_sig, down_sig=None, conn_flag =  True, umap_flag = True):
    exp = pd.DataFrame(obj.layers["velocity"], columns=obj.var.index, index=obj.obs.index)
    exp = exp.T.dropna()
    scaler = StandardScaler()
    exp = scaler.fit_transform(exp.T).T      #
    exp = obj.to_df().T
    sigs_scores = signature_values(exp, up_sig, down_sig) # wilcoxon score
    if conn_flag:
        graph = obj.obsp["connectivities"].toarray()
    else:
        graph = obj.obsp["distances"].toarray()
    #sigs_scores = propagation(sigs_scores,graph)
    obj.obs["SigScore"] = sigs_scores
    if umap_flag:
        sc.pl.umap(obj, color=["SigScore"],color_map="magma")
    else:
        sc.pl.tsne(obj, color=["SigScore"],color_map="magma")

    return obj