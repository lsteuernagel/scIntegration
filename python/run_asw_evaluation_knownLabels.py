# Import relevant modules
import pandas as pd
import numpy as np
import scanpy as sc
import louvain
import igraph
import re
import sklearn
import sklearn.metrics
import os
import sys
import json
import gc

# find file helper
def find(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if pattern in name:
                result.append(os.path.join(root, name))
    return result

print("Starting ASW evaluation.")

# Parameters
# get path to parameters and path json
json_file = open(sys.argv[1])
# read json into dictionary
json_str = json_file.read()
parameter_dict = json.loads(json_str)

# general params from dict (could also save these lines, but I think this way it is easier to digest)
integration_folder_path = parameter_dict["integration_folder_path"]
merged_file_h5ad = parameter_dict["merged_file_h5ad"]
integration_results_path = parameter_dict["integration_res_path"]
evaluation_file = parameter_dict["evaluation_file"]
evaluation_file_grouped = parameter_dict["evaluation_file_grouped"] 
global_seed = parameter_dict["global_seed"]
embedding_names = parameter_dict["integration_names"]
label_file_name = parameter_dict["label_file_name"] 
job_id = parameter_dict["job_id"]

print("Read anndata.")

# read adata
adata = sc.read_h5ad(merged_file_h5ad)
#ensure there are no bytestrings 
str_df = adata.obs
str_df = str_df.applymap(lambda x: x.decode() if isinstance(x, bytes) else x)
str_df = str_df.set_index('Cell_ID',drop=False)
adata.obs = str_df
# for features:
str_df = adata.var
str_df = str_df.applymap(lambda x: x.decode() if isinstance(x, bytes) else x)
str_df = str_df.set_index('features',drop=False)
adata.var = str_df

print("Subset anndata.")

# read  label_file_name and subset if wanted
label_df = pd.read_csv(label_file_name,skip_blank_lines=True,header=0)
cell_ids = label_df.iloc[:, 0]
adata = adata[cell_ids, :].copy()

# add labels
cell_labels = label_df.iloc[:, 1]
cell_labels = cell_labels.set_axis(cell_ids)
adata.obs = adata.obs.assign(cell_labels = cell_labels)
print("subsetted to "+str(adata.n_obs))

# clean up
gc.collect()

#### Calculate asw
#init df
asw_all = pd.DataFrame(columns=['reduction','silhouette_score_euclidean','silhouette_score_cosine'])
asw_all_grouped = pd.DataFrame(columns=['reduction','cell_labels','silhouette_score_euclidean','silhouette_score_cosine'])

# read all files
for embedding in embedding_names:
  
    ## load embedding and add (plus subset to relevant cells!)
    print("embedding: "+embedding)
    embedding_path = find(embedding, integration_results_path)[0]
    # add to adata
    df = pd.read_csv(embedding_path,sep="\t",skip_blank_lines=True,index_col=0)
    # subset to relevant cell ids
    df=df.loc[cell_ids] 
    adata.obsm[embedding] = df.to_numpy()
    #print(str(df.shape[0]))
    #print(str(adata.obsm[embedding].shape[0]))
    
    ## Calculate ASW
    # see also: https://scikit-learn.org/stable/modules/clustering.html#silhouette-coefficient
    print(" Calculating asw metrics")
    asw_euclidean = sklearn.metrics.silhouette_samples(adata.obsm[embedding], adata.obs['cell_labels'],metric='euclidean')
    asw_cosine = sklearn.metrics.silhouette_samples(adata.obsm[embedding], adata.obs['cell_labels'],metric='cosine')
    
    ## Mean over all asw per sample
    add = pd.DataFrame({'reduction' : [embedding], 
                        'silhouette_score_euclidean' : [asw_euclidean.mean()], 
                        'silhouette_score_cosine' : [asw_cosine.mean()]
                        })
    asw_all = asw_all.append(add)
    
    ## Mean over each label category asw per sample
    # join in pandas data frame for easy use of groupby 
    asw_res_df = pd.DataFrame({'cell_labels' : adata.obs['cell_labels'], 
                               'silhouette_score_euclidean' : asw_euclidean, 
                               'silhouette_score_cosine' : asw_cosine
                           })
    # take grouped means
    per_label_mean_df = asw_res_df.groupby(['cell_labels']).mean()
    per_label_mean_df['reduction'] = embedding # add reduction
    per_label_mean_df['cell_labels'] = per_label_mean_df.index
    asw_all_grouped = asw_all_grouped.append(per_label_mean_df) ## append
    
    # clean up
    gc.collect()

# set index
asw_all.index = list(range(0, asw_all.shape[0], 1))  
asw_all_grouped.index = list(range(0, asw_all_grouped.shape[0], 1))  

# save results by appending to existing file ---> create file first!
asw_all.to_csv(evaluation_file, sep='\t',index=False, mode='a', header=False)
asw_all_grouped.to_csv(evaluation_file_grouped, sep='\t',index=False, mode='a', header=False)

print("Finalized ASW evaluation job.")




















    
