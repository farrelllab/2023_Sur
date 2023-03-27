#!/bin/bash -x
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=8:00:00 # (8 hours)
#SBATCH --partition=norm
#SBATCH --mem=96g # 96g RAM
#SBATCH --open-mode=append
#SBATCH --array=1-60
#SBATCH -o /data/CSD/zfext/LTA/logs/11-Modules/FC-tissue_%A-%a.out # Standard out goes to this file
#SBATCH -e /data/CSD/zfext/LTA/logs/11-Modules/FC-tissue_%A-%a.err # Standard err goes to this file

# Paths
script_path="/data/CSD/zfext/LTA/scripts/11-Modules/tissue_fuzzyc_cluster.R"
file_list_path="/data/CSD/zfext/LTA/scripts/11-Modules/SampleSheet.tsv"

# Get sample ID name from file_list
param_c=`sed -n "$SLURM_ARRAY_TASK_ID"p "${file_list_path}" | cut -f 1`
param_m=`sed -n "$SLURM_ARRAY_TASK_ID"p "${file_list_path}" | cut -f 2`
param_smooth=`sed -n "$SLURM_ARRAY_TASK_ID"p "${file_list_path}" | cut -f 3`
param_tissue=`sed -n "$SLURM_ARRAY_TASK_ID"p "${file_list_path}" | cut -f 4`

# Run R script
Rscript ${script_path} ${param_c} ${param_m} ${param_smooth} ${param_tissue}
