#!/bin/bash

#SBATCH --job-name=mt08_00_LDSC_ins
#SBATCH --time=24:00:00
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=64G
#SBATCH --account=project_2007428
#SBATCH --output=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/mt08_00_LDSC_ins.log
#SBATCH --error=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/mt08_00_LDSC_ins.log

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
echo "Rscript ${mainDir}scripts/meta_step8_meta_LDSC_prep.R"
printf "\n\n"
printf "==========================================\n\n"


Rscript ${mainDir}scripts/meta_step8_meta_LDSC_prep.R


printf "\n\n"
printf "ALL PLOTS MADE.\n\n"
date
printf "\n\n"
printf "==========================================\n\n"
printf "Script complete. Goodbye.\n\n"
date
printf "==========================================\n\n"

