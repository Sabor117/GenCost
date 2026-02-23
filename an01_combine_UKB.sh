#!/bin/bash

#SBATCH --job-name=an01_combine_UKB_%a
#SBATCH --time=3-00:00:00
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=128G
#SBATCH --account=project_2007428
#SBATCH --array=1-72
#SBATCH --output=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/an01_combine_UKB_ALL_out_%a.log
#SBATCH --error=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/an01_combine_UKB_ALL_out_%a.log

module load r-env

mainDir="/scratch/project_2007428/projects/prj_001_cost_gwas/"

printf "Script starts:\n\n"
date
printf "\n\n==========================\n\n"

printf "First command is:\n"
echo "Rscript ${mainDir}scripts/an_step1_1_UKB_file_combine.R ${SLURM_ARRAY_TASK_ID}"
printf "\n\n"
printf "==========================================\n\n"


Rscript ${mainDir}scripts/an_step1_1_UKB_file_combine.R ${SLURM_ARRAY_TASK_ID}


printf "\n\n"
printf "R SCRIPT COMPLETE.\n\n"
printf "\n\n"
printf "==========================================\n\n"
printf "Script complete. Goodbye.\n\n"
date
printf "==========================================\n\n"


