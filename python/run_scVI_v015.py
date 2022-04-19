# Run python scVI
# runs different parameters on one assay+featureset

# Import relevant modules
import sys
import pandas as pd
import numpy as np
import scanpy as sc
import louvain
import igraph
import os
import sys
import json
#import re
import scvi
import gc
from os import listdir
from os.path import isfile, join

print(scvi.__version__)

print("Read parameters")
# get path to parameters and path json
json_file = open(sys.argv[1])
# read json into dictionary
json_str = json_file.read()
parameter_dict = json.loads(json_str)

# general params from dict (could also save these lines, but I think this way it is easier to digest)
# general params from dict (could also save these lines, but I think this way it is easier to digest)
job_id=parameter_dict["job_id"]
batch_var=parameter_dict["batch_var"]
data_filepath_full = parameter_dict["data_filepath_full"]
feature_set_file = parameter_dict["feature_set_file"]
results_path = parameter_dict["output_folder"]
global_seed = parameter_dict["global_seed"]
hvgs_set_name = parameter_dict["hvgs_set_name"]
param_path = parameter_dict["param_path"]
use_cuda = parameter_dict["use_cuda"]
categorical_covariates = parameter_dict["categorical_covariates"]
continuous_covariates = parameter_dict["continuous_covariates"]

scvi.settings.seed = global_seed

# read parameters for scVI
param_df = pd.read_csv(param_path,sep="\t")

# read RNA assay (RNA assay because scVI needs raw data)
print("Read anndata")
adata = sc.read_h5ad(data_filepath_full)
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

# read json into dictionary
json_file = open(feature_set_file)
hvg_dict = json.load(json_file)
hvgs = hvg_dict[hvgs_set_name] # add in hvgs variable

# print("Subset anndata to selected featureset:"+str(hvgs_set_name))
# feature_set_files = [f for f in listdir(filepath_sets) if isfile(join(filepath_sets, f))]
# feature_set_files = [i for i in feature_set_files if hvgs_set_name in i]
# file = feature_set_files[0] # expects only 1 so just use first (should be only) entry
# mydf = pd.read_csv(filepath_sets+file,skip_blank_lines=True)
# hvgs = mydf.iloc[:,0].to_numpy().tolist() # add in hvgs variable

# scVI needs raw data and so we subset X to .raw and then relevant features
print("Copy raw into .X")
adata.X = adata.raw.X.copy()
adata = adata[:, hvgs]

# clean up
gc.collect()

## Run scVI
# https://www.scvi-tools.org/en/stable/user_guide/notebooks/harmonization.html
print("Preparing scVI")
# if only 1 variable --> make list
if not isinstance(categorical_covariates, list):
    categorical_covariates_save=categorical_covariates
    categorical_covariates = list()
    categorical_covariates.insert(0,categorical_covariates_save)
if not isinstance(continuous_covariates, list):
    continuous_covariates_save=continuous_covariates
    continuous_covariates = list()
    continuous_covariates.insert(0,continuous_covariates_save)
length_cov = len(categorical_covariates)+len(continuous_covariates)
if len(categorical_covariates) == 0:
    categorical_covariates = None
if len(continuous_covariates) == 0:
    continuous_covariates = None
print("categorical_covariate_keys: "+str(categorical_covariates))
print("continuous_covariate_keys: "+str(continuous_covariates))
adata_scvi = adata.copy()
# setup for scvi
# Imporant: I am only using the categorical_covariate_keys to describe the batch variable(s) --> unambigous way to specify 1-many cvovariates as batches
scvi.model.SCVI.setup_anndata(adata_scvi, categorical_covariate_keys=categorical_covariates,continuous_covariate_keys=continuous_covariates)

# for all parameter combinations
for index, row in param_df.iterrows():
    print("Running "+str(index+1)+" of "+str(len(param_df.index))+" scVI runs")
    # set up VAE
    vae = scvi.model.SCVI(adata_scvi,
        n_layers=int(row['n_layers']),
        n_latent=int(row['n_latent']),
        n_hidden=int(row['n_hidden']),
        dropout_rate=float(row['dropout_rate']),
        dispersion = str(row['dispersion']),
        gene_likelihood = str(row['gene_likelihood'])
        #use_cuda=use_cuda
    )

    # train
    vae.train(max_epochs=int(row['max_epochs']), early_stopping = bool(row['early_stopping']), use_gpu = use_cuda)
    # get result
    adata_scvi.obsm["X_scVI"] = vae.get_latent_representation()
    output = pd.DataFrame(adata_scvi.obsm["X_scVI"])
    output = output.set_index(adata_scvi.obs_names)
    output2 = output.set_axis(["scVI_" + str(s) for s in output.axes[1].to_list()], axis=1, inplace=False)
    # save
    output2.to_csv(results_path+"scVI_"+str(index)+"_"+str(int(row['max_epochs']))+"_"+
                  str(float(row['dropout_rate']))+"_"+str(int(row['n_layers']))+"_"+
                  str(int(row['n_hidden']))+"_"+str(row['dispersion'])+"_"+str(row['gene_likelihood'])+"_cov"+str(length_cov)+
                  "..scVI.."+str(int(row['n_latent']))+".."+hvgs_set_name+"_"+job_id+".txt", sep='\t',index=True)
    # clean up
    gc.collect()
                  
## End of for
print("Finalized scVI runs")















