#!/bin/bash -x
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=240:00:00 # (7 days)
#SBATCH --partition=norm
#SBATCH --mem=160g # 40g RAM
#SBATCH --open-mode=append
#SBATCH --array=1-4
#SBATCH -o /data/CSD/zfext/LTA/logs/11-Modules/FC-mama_%A-%a.out # Standard out goes to this file
#SBATCH -e /data/CSD/zfext/LTA/logs/11-Modules/FC-mama_%A-%a.err # Standard err goes to this file

# Paths
script_path="/data/CSD/zfext/LTA/scripts/11-Modules/mamasmooth_fuzzyc_cluster.R"
file_list_path="/data/CSD/zfext/LTA/scripts/11-Modules/mamasmoothSample.tsv"

# Get sample ID name from file_list
param_c=`sed -n "$SLURM_ARRAY_TASK_ID"p "${file_list_path}" | cut -f 1`
param_m=`sed -n "$SLURM_ARRAY_TASK_ID"p "${file_list_path}" | cut -f 2`

# Run R script
Rscript ${script_path} ${param_c} ${param_m}
