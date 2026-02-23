#!/bin/bash

#SBATCH --job-name=mt011_2_coloc_run_%a
#SBATCH --time=3-00:00:00
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=128G
#SBATCH --account=project_2007428
#SBATCH --array=22
#SBATCH --output=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/mt011_2_coloc_run_%a.log
#SBATCH --error=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/mt011_2_coloc_run_%a.log

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
echo "Rscript ${mainDir}scripts/meta_step11_coloc_run.R --run_no ${SLURM_ARRAY_TASK_ID} --pheno IN"
printf "\n\n"
printf "==========================================\n\n"


Rscript ${mainDir}scripts/meta_step11_coloc_run.R  --run_no ${SLURM_ARRAY_TASK_ID} --pheno IN

printf "\n\n"
printf "COLOC IN COMPLETE.\n\n"
echo "Rscript ${mainDir}scripts/meta_step11_coloc_run.R --run_no ${SLURM_ARRAY_TASK_ID} --pheno DRUG"
printf "\n\n"
printf "==========================================\n\n"

Rscript ${mainDir}scripts/meta_step11_coloc_run.R  --run_no ${SLURM_ARRAY_TASK_ID} --pheno DRUG

printf "\n\n"
printf "COLOC DRUG COMPLETE.\n\n"
echo "Rscript ${mainDir}scripts/meta_step11_coloc_run.R --run_no ${SLURM_ARRAY_TASK_ID} --pheno PRIM"
printf "\n\n"
printf "==========================================\n\n"

Rscript ${mainDir}scripts/meta_step11_coloc_run.R  --run_no ${SLURM_ARRAY_TASK_ID} --pheno PRIM

printf "\n\n"
printf "COLOC PRIM COMPLETE.\n\n"
echo "Rscript ${mainDir}scripts/meta_step11_coloc_run.R --run_no ${SLURM_ARRAY_TASK_ID} --pheno INOUT"
printf "\n\n"
printf "==========================================\n\n"

Rscript ${mainDir}scripts/meta_step11_coloc_run.R  --run_no ${SLURM_ARRAY_TASK_ID} --pheno INOUT


printf "\n\n"
printf "COLOC INOUT COMPLETE.\n\n"
printf "\n\n"
printf "==========================================\n\n"
printf "Script complete. Goodbye.\n\n"
date
printf "==========================================\n\n"
