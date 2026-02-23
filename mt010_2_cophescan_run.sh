#!/bin/bash

#SBATCH --job-name=mt10_2_cophescan_run_%a
#SBATCH --time=3-00:00:00
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=128G
#SBATCH --account=project_2007428
#SBATCH --array=1-3
#SBATCH --output=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/mt10_2_cophescan_run_%a.log
#SBATCH --error=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/mt10_2_cophescan_run_%a.log

##### =========================== #####

### Setup environment

##### =========================== #####

module load r-env

### Clean up .Renviron file in home directory

if test -f ~/.Renviron; then
     sed -i '/TMPDIR/d' ~/.Renviron
fi

### Specify a temp folder path

echo "TMPDIR=/scratch/project_2007428/projects/prj_001_cost_gwas/tmpdir/" >> ~/.Renviron


##### =========================== #####

### Start script

##### =========================== #####

mainDir="/scratch/project_2007428/projects/prj_001_cost_gwas/"

printf "Script starts:\n\n"
printf "NOBODY EXPECTS THE SEB INQUISTION:\n"
date
printf "\n\n==========================\n\n"

printf "First command is:\n"
echo "Rscript ${mainDir}scripts/meta_step10_2_cophescan_run.R ${SLURM_ARRAY_TASK_ID}"
printf "\n\n"
printf "==========================================\n\n"


Rscript ${mainDir}scripts/meta_step10_2_cophescan_run.R ${SLURM_ARRAY_TASK_ID}


printf "\n\n"
printf "RUNNING RUN.\n\n"
printf "Inquisition was expected after all.\n\n"
date
printf "\n\n"
printf "==========================================\n\n"
printf "Script complete. Goodbye.\n\n"
date
printf "==========================================\n\n"

