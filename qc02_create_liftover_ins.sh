#!/bin/bash

#SBATCH --job-name=qc02_Liftover_Inputs
#SBATCH --time=3-00:00:00
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=200G
#SBATCH --account=project_2007428
#SBATCH --output=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/qc02_liftover_ins.log
#SBATCH --error=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/qc02_liftover_ins.log

module load r-env

mainDir="/scratch/project_2007428/projects/prj_001_cost_gwas/"

printf "Script starts:\n\n"
date
printf "\n\n==========================\n\n"

printf "First command is:\n"
echo "Rscript ${mainDir}scripts/qc_step2_make_liftover_ins.R"
printf "\n\n"
printf "==========================================\n\n"


Rscript ${mainDir}scripts/qc_step2_make_liftover_ins.R


printf "\n\n"
printf "ALL INPUTS MADE.\n\n"
printf "\n\n"
printf "==========================================\n\n"
printf "Script complete. Goodbye.\n\n"
date
printf "==========================================\n\n"
