# Import relevant modules
import pandas as pd
import numpy as np
import scanpy as sc
import louvain
import igraph
import re
import sklearn
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

# Parameters
# get path to parameters and path json
json_file = open(sys.argv[1])
# read json into dictionary
json_str = json_file.read()
parameter_dict = json.loads(json_str)

# general params from dict (could also save these lines, but I think this way it is easier to digest)
batch=parameter_dict["batch_var"]
integration_folder_path = parameter_dict["integration_folder_path"]
merged_file_h5ad = parameter_dict["merged_file_h5ad"]
integration_results_path = parameter_dict["integration_res_path"]
evaluation_file = parameter_dict["evaluation_file"]
global_seed = parameter_dict["global_seed"]
embedding_names = parameter_dict["integration_names"]
target_clusterN = parameter_dict["target_clusterN"]
id_file_name = parameter_dict["id_file_name"] 
subset_cells = parameter_dict["subset_cells"] 
job_id = parameter_dict["job_id"]
start_res = parameter_dict["start_res_asw"]
end_res = parameter_dict["end_res_asw"]
n_neighbors = parameter_dict["k_param"]

#define resolution range, hardcoded atm
# start_res=1
# end_res=12
resolutions = [x*0.5 for x in range(2*start_res, 2*end_res+1)]

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

# read id file name and subset if wanted
if subset_cells:
  #all_cells = pd.read_csv(integration_folder_path+id_file_name,skip_blank_lines=True,header=None)
  json_file = open(integration_folder_path+id_file_name)
  cell_ids = json.load(json_file)
  cell_ids = cell_ids['cell_ids']
  adata = adata[cell_ids, :].copy()
  print("subsetted to "+str(adata.n_obs))

# Calculate asw
#init df
asw_all = pd.DataFrame(columns=['embedding','resolution','n_clusters','silhouette_score_euclidean','silhouette_score_cosine','calinski_harabasz','davies_bouldin'])

# clean up
gc.collect()

# read all files
for embedding in embedding_names:
    print("embedding: "+embedding)
    embedding_path = find(embedding, integration_results_path)[0]
    # add to adata
    df = pd.read_csv(embedding_path,sep="\t",skip_blank_lines=True,index_col=0)
    if subset_cells:
      df=df.loc[cell_ids]
    adata.obsm[embedding] = df.to_numpy()
    print(str(df.shape[0]))
    print(str(adata.obsm[embedding].shape[0]))
    # run pp.neighbors
    print(" Finding neighbor graph")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=df.shape[1],use_rep=embedding,key_added=embedding,random_state=global_seed)
    # run clustering and ASW
    print(" Running clustering at increasing resolutions")
    for res in resolutions:
        key_name = str(embedding)+"_clusters_"+str(res)
        #print("  key: "+key_name)
        sc.tl.leiden(adata,resolution=res,key_added=key_name,random_state=global_seed,neighbors_key=embedding) # could use neighbors_key
        print(" Ran leiden with resolution "+str(res)+" and found "+str(len(set(adata.obs[key_name])))+" clusters")
        if(len(set(adata.obs[key_name])) >= target_clusterN):
            print(" Reached "+str(len(set(adata.obs[key_name])))+" clusters")
            break
    # calculate separation on last result
    # see also: https://scikit-learn.org/stable/modules/clustering.html#silhouette-coefficient
    print(" Calculating asw metrics")
    asw_euclidean = sklearn.metrics.silhouette_samples(adata.obsm[embedding], adata.obs[key_name],metric='euclidean')
    asw_cosine = sklearn.metrics.silhouette_samples(adata.obsm[embedding], adata.obs[key_name],metric='cosine')       
    print(" Calculating calinski_harabasz_score")
    calinski_harabasz_score = sklearn.metrics.calinski_harabasz_score(adata.obsm[embedding], adata.obs[key_name])
    print(" Calculating davies_bouldin_score")
    davies_bouldin_score = sklearn.metrics.davies_bouldin_score(adata.obsm[embedding], adata.obs[key_name])
    # add last result:
    add = pd.DataFrame({'reduction' : [embedding], 
                            'resolution' : [res], 
                            'n_clusters' : [len(set(adata.obs[key_name]))], 
                            'silhouette_score_euclidean' : [asw_euclidean.mean()], 
                            'silhouette_score_cosine' : [asw_cosine.mean()], 
                            'calinski_harabasz' : [calinski_harabasz_score.mean()],
                            'davies_bouldin' : [davies_bouldin_score.mean()]
                       })
    #batch_silhouette_score' : [asw_batch[0]['silhouette_score'].mean()]})
    asw_all = asw_all.append(add)
    # clean up
    gc.collect()
# set index
asw_all.index = list(range(0, asw_all.shape[0], 1))  

# save results by appending to existing file ---> create file first!
asw_all.to_csv(evaluation_file, sep='\t',index=False, mode='a', header=False)

print("Finalized ASW evaluation job.")




















    
