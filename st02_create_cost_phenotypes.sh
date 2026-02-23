#!/bin/bash

#SBATCH --job-name=st02_make_cost_phenos
#SBATCH --time=3-00:00:00
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=128G
#SBATCH --account=project_2007428
#SBATCH --output=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/st02_phenotype_make.log
#SBATCH --error=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/st02_phenotype_make.log

module load r-env

mainDir="/scratch/project_2007428/projects/prj_001_cost_gwas/"

### Clean up .Renviron file in home directory

if test -f ~/.Renviron; then
     sed -i '/TMPDIR/d' ~/.Renviron
fi

### Specify a temp folder path

echo "TMPDIR=/scratch/project_2007428/projects/prj_001_cost_gwas/tmpdir/" >> ~/.Renviron

printf "Script starts:\n\n"
date
printf "\n\n==========================\n\n"

printf "First command is:\n"
echo "Rscript ${mainDir}scripts/st_step3_create_cost_phenotype.R"
printf "\n\n"
printf "==========================================\n\n"


Rscript ${mainDir}scripts/st_step3_create_cost_phenotype.R


printf "\n\n"
printf "Phenotypes made.\n\n"
date
printf "\n\n"
printf "==========================================\n\n"
printf "Script complete. Goodbye.\n\n"
date
printf "==========================================\n\n"

