#!/bin/bash
#SBATCH --job-name='scHarmonize_integration'
#SBATCH --output='/beegfs/scratch/bruening_scratch/lsteuernagel/data/scHarmonize/hypothalamusMapFull_v4/log_files/%j-scHarmonize_integration.out'
#SBATCH --error='/beegfs/scratch/bruening_scratch/lsteuernagel/data/scHarmonize/hypothalamusMapFull_v4/log_files/%j-scHarmonize_integration.err'
#SBATCH --ntasks=1
#SBATCH --time=240:00:00
#SBATCH --cpus-per-task=54
#SBATCH --partition=blade-b

# need to get 3 variables from call
# singularity image
singularity_image=$1
# $method_script : script that should be executed (depend on method)
method_script=$2
# $param_file : specifies which object, assay etc are used and which parameters of the method
param_file=$3
# $logfile: logfile with clear id if not slurm, else it will use %j
#logfile=$2

# Run script 
echo "singularity exec "$singularity_image" python -u "$method_script" "$param_file
srun singularity exec $singularity_image python -u $method_script $param_file


