# scIntegration
Find optimal integration for a merged single cell dataset

# Overview 

This repository contains a pipeline to integrate a merged single cell dataset using multiple methods and hyperparameter sets and evaluate the resulting integrated low dimensional embeddings with regard to their mixing of Batches / Samples / Datasets and purity of cell types.

This package was designed to integrate the HypoMap data.

After finding good hyperparameters for integrating the data, the: [scHarmonization](https://github.com/lsteuernagel/scHarmonization) pipeline can bu used to run downstream steps such as clustering and marker detection to create a fully harmonized and unified reference dataset such as HypoMap.

Key steps involve:

- Preparation of files, feature sets  via [Seurat](https://satijalab.org/seurat/) 
- Prepare the evaluation basis using short signatures of well-defined marker genes to objectively identify 'ground-truth' cell types across batches
- Running PCA via [Seurat](https://satijalab.org/seurat/) and saving the integrated low dimensional embeddings 
- Running [scVI](https://scvi-tools.org/) and saving the integrated low dimensional embeddings 
- Evaluating the resulting methods using 4 approaches:
  - kNN based mixing: entropy of dataset distribution of neighbors per cell
  - kNN based purity: percentage of neighbors annotated with the same cell type per cell
  - random forest based mixing: out-of-bag prediction performance of a random forest trained to classify batches (i.e. bad performace means good mixing) per cell (but only on a downsampled set of all cells)
  - average silhouette width of clusters based on the integrated embedding to measure general (biology independent) cluster separation [sklearn: Silhouette Coefficient](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.silhouette_score.html)
- Visualizing the overall results to choose a final integration hyperparameter set and obtain a low dimensional embedding of the integrated data 

Additionally the pipeline relies on the following R package: [scUtils](https://github.sf.mpg.de/lsteuernagel/scUtils) and some other minor dependencies.

Most pipeline steps are executed via slurm jobs using a singularity or docker image (with all required dependencies) that can be pulled via dockerhub: [Docker image](https://hub.docker.com/r/lsteuernagel/r_scvi_v2/tags).

The execution and queue management scripts are written in R and require some packages like digest (for hashing), magrittr and dplyr.

# Step-by step 

This section intends to explain the different scripts involved and should allow to reproduce the above steps.

- The pipeline requires a seurat object stored as an rds file, a json with features to exclude from HVG search and a json with multiple entry (one for each cell type), containing genes names.

- The parameters for the pipeline are stored in json format as well (see below).

- The main script is [R/executeIntegration.R] to prepare the integration objects etc. and run the hyperparameter serach with scVI. See the script for detils on the slurm pipeline and which scripts are used for each of the steps outlined above.

- For evaluation [R/executescEvaluation.R] runs the evaluation part of the pipeline. See the script for detils on the slurm pipeline and which scripts are used for each of the steps outlined above.

### Parameters

The parameters are loaded from a json files that can be adjusted (or copied and edited). [R/make_parameter_json.R] can be used to create the json.

For the hypoMap pipeline we used this parameter file: [parameters_integration_v2_2.json] 

It looks like this:

```yaml
{
  "integration_folder_path": "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_v2_integration/",
  "merged_file": "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/hypoMap_merged_filtered.rds",
  "new_name_suffix": "hypoMap_integrated",
  "celltype_signature_file": "data/hypothalamus_celltype_signatures.json",
  "auc_backup_file": "data/hypoMap_celltype_auc_per_cell_result.txt",
  "genes_to_exclude_file": "data/features_exclude_list.json",
  "n_cores": 50,
  "id_column": "Cell_ID",
  "global_seed": 123456,
  "sample_column": "Sample_ID",
  "batch_var": "Batch_ID",
  "feature_set_sizes": [750, 1000, 1250, 1500, 2000, 2500, 3000],
  "assay_name": "RNA",
  "auc_max_rank": 700,
  "block_size": 10000,
  "alpha": 2,
  "thrP": 0.01,
  "smallestPopPercent": 0.01,
  "auc_max_pos_thresh": 0.05,
  "auc_min_pos_thresh": 0.15,
  "detected_cells_filename": "detected_celltypes.json",
  "id_file_name": "downsampled_ids_for_evaluation.json",
  "target_sub_sample": 38700,
  "stepsize": 200,
  "scvi_models_to_test": 100,
  "latent_space_sizes": [50, 65, 80, 95, 110, 140],
  "ntrees_mixing": 20000,
  "sampsize_pct": 0.3333,
  "max_for_norm": 0.01,
  "k_param": 20,
  "dist_type": "cosine",
  "subset_cells": true,
  "target_clusterN": 150,
  "start_res_asw": 5,
  "end_res_asw": 15
}
```

The first parameters manage the input and put put files.

### Key Inputs to define

- **merged_file** points to a merged seurat object of all datasets that has to be build before running this pipeline. 
- **integration_folder_path** is the output directory for all intermediate and final results. This should be different for each project. 
- **celltype_signature_file** requires a json file with named gene lists that represent signatures of ground truth cell types used for evaluation


### Other parameters descriptions (selected):

- **feature_set_sizes** Different numbers of highly variable features to test.
- **sample_column** column in meta.data of merged seurat object with sample level information
- **batch_var** column in meta.data of merged seurat object with batch level information --> This will be the one that is used in the integration
- **auc_max_rank, block_size, alpha, thrP, smallestPopPercent,auc_max_pos_thresh, auc_min_pos_thresh** various parameters for running AUC and then determining thresholds to define which cells are part of a cell type (for purity evaluation)
- **target_sub_sample** some evaluation functions only use a subset of this size of the data
- **scvi_models_to_test** the pipeline uses a "random search"-like approach where out of the grid of all parameter combinations only this many will be tested
- **latent_space_sizes** scvi laten space sie parameters
- **ntrees_mixing, sampsize_pct, max_for_norm** random forest evaluation related paramters.
- **k_param** how many neighbors for nerighbor evaluation
- **subset_cells**
- **target_clusterN, start_res_asw, end_res_asw** ASW is calcualated on a set of preliminary clusters. To determine theses clusters we iterate over different leiden resolutions until target_clusterN is reached.



