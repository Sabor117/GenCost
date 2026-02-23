#!/bin/bash

#SBATCH --job-name=qc03_GnoMAD_freq_check_%a
#SBATCH --time=3-00:00:00
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=300G
#SBATCH --account=project_2007428
#SBATCH --array=1-16
#SBATCH --output=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/qc03_GnoMAD_freq_check_%a.log
#SBATCH --error=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/qc03_GnoMAD_freq_check_%a.log

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
date
printf "\n\n==========================\n\n"

printf "First command is:\n"
echo "Rscript ${mainDir}scripts/qc_step3_gnomad_freq_check.R ${SLURM_ARRAY_TASK_ID}"
printf "\n\n"
printf "==========================================\n\n"


Rscript ${mainDir}scripts/qc_step3_gnomad_freq_check.R ${SLURM_ARRAY_TASK_ID}


printf "\n\n"
printf "ALL FREQUENCIES COMPARED.\n\n"
printf "\n\n"
printf "==========================================\n\n"
printf "Script complete. Goodbye.\n\n"
date
printf "==========================================\n\n"
