#!/bin/bash

#SBATCH --job-name=ldpredPRS_01_weights_calc_%a
#SBATCH --time=14-00:00:00
#SBATCH --partition=longrun
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=360G
#SBATCH --account=project_2007428
#SBATCH --array=1-6
#SBATCH --output=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/ldpred2_logs/ldpredPRS_01_weights_calc_%a.log
#SBATCH --error=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/ldpred2_logs/ldpredPRS_01_weights_calc_%a.log

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

printf "Script starts. NOBODY EXPECTS THE SEB INQUISITION:\n\n"
date
printf "\n\n==========================\n\n"

printf "First command is:\n"
echo "Rscript ${mainDir}scripts/ldpred2_pipeline/ldpred2_step1_make_weights.R ${SLURM_ARRAY_TASK_ID}"
printf "\n\n"
printf "==========================================\n\n"


Rscript ${mainDir}scripts/ldpred2_pipeline/ldpred2_step1_make_weights.R ${SLURM_ARRAY_TASK_ID}


printf "\n\n"
printf "Correlations created.\n\n"
date
printf "\n\n"
printf "==========================================\n\n"
printf "Script complete. Okay, somebody expected us. Goodbye.\n\n"
date
printf "==========================================\n\n"
