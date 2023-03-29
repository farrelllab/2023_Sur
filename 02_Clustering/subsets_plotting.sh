#!/bin/bash -x
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -J subsets_plotting
#SBATCH --time=12:00:00 # (12 hours)
#SBATCH --partition=norm
#SBATCH --mem=128g
#SBATCH --open-mode=append
#SBATCH -o /data/CSD/zfext/logs/04d-SubsetsV4/subsets_plotting_%A-%a.out # Standard out goes to this file
#SBATCH -e /data/CSD/zfext/logs/04d-SubsetsV4/subsets_plotting_%A-%a.err # Standard err goes to this file
#SBATCH --array=1-20

# Variables
#output_dir="/data/CSD/zfext/results/03a-cNMF/"
R_script_loc="/data/CSD/zfext/scripts/04d-SubsetsV4/subsets_plotting.r"
sample_list="/data/CSD/zfext/scripts/04d-SubsetsV4/sample_list_new_mama.txt"

# Get sample name out of list
sample_name=`sed -n "$SLURM_ARRAY_TASK_ID"p "${sample_list}" | cut -f 1`

# Run R Script
Rscript ${R_script_loc} ${sample_name}	
