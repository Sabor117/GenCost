#!/bin/bash

#SBATCH --job-name=mt06_snp_matrices
#SBATCH --time=3-00:00:00
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=200G
#SBATCH --account=project_2007428
#SBATCH --output=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/mt06_snp_matrices.log
#SBATCH --error=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/mt06_snp_matrices.log

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
metaDir=${mainDir}"outputs/METAL_v4/"

printf "Script starts:\n\n"
date
printf "\n\n==========================\n\n"

printf "First command is:\n"
echo "Rscript ${mainDir}scripts/meta_step6_top_snp_comparison.R"
printf "\n\n"
printf "==========================================\n\n"


Rscript ${mainDir}scripts/meta_step6_top_snp_comparison.R


printf "\n\n"
printf "SNP comparison frames created.\n\n"
printf "==========================================\n\n"
printf "Script complete. Goodbye.\n\n"
date
printf "==========================================\n\n"

