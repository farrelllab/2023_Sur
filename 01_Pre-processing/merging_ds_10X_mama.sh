#!/bin/bash -x
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -J mama_ds_merge
#SBATCH --time=200:00:00 # (200 hours)
#SBATCH --partition=largemem
#SBATCH --mem=800g
#SBATCH --open-mode=append
#SBATCH -o /data/CSD/zfext/logs/04c-SubsetsV3/mama_merge_%A-%a.out # Standard out goes to this file
#SBATCH -e /data/CSD/zfext/logs/04c-SubsetsV3/mama_merge_%A-%a.err # Standard err goes to this file


# Variables
#output_dir="/data/CSD/zfext/results/03a-cNMF/"
R_script_loc="/data/CSD/zfext/scripts/04d-SubsetsV4/merging_mama_dropseq.R"
#sample_list="/data/CSD/zfext/scripts/04c-SubsetsV3/sample_list_classify.txt"

# Get sample name out of list
#sample_name=`sed -n "$SLURM_ARRAY_TASK_ID"p "${sample_list}" | cut -f 1`

# Run R Script
Rscript ${R_script_loc} 
