#!/bin/bash -x
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=5
#SBATCH -J matrix_distances_eps
#SBATCH --time=80:00:00 # (80 hours)
#SBATCH --partition=largemem
#SBATCH --mem=1000g
#SBATCH --open-mode=append
#SBATCH -o /data/CSD/zfext/LTA/logs/13-EPS/matrix_distances_eps_%A-%a.out # Standard out goes to this file
#SBATCH -e /data/CSD/zfext/LTA/logs/13-EPS/matrix_distances_eps_%A-%a.err # Standard err goes to this file
#SBATCH --array=1-20

# Variables
#output_dir="/data/CSD/zfext/results/03a-cNMF/"
R_script_loc="/data/CSD/zfext/LTA/scripts/13-EPS/matrices_distance_threshold.r"
sample_list="/data/CSD/zfext/scripts/04d-SubsetsV4/sample_list_new_mama.txt"

# Get sample name out of list
sample_name=`sed -n "$SLURM_ARRAY_TASK_ID"p "${sample_list}" | cut -f 1`

# Run R Script
Rscript ${R_script_loc} ${sample_name}	
