#!/bin/bash

#SBATCH --job-name=mgPRS_01_prepare_sumstats
#SBATCH --time=3-00:00:00
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=64G
#SBATCH --account=project_2007428
#SBATCH --output=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/megaPRS_logs/mgPRS_01_prepare_sumstats.log
#SBATCH --error=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/megaPRS_logs/mgPRS_01_prepare_sumstats.log

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
scriptDir=${mainDir}"scripts/megaPRS_pipeline/"

printf "Script starts:\n\n"
date
printf "\n\n==========================\n\n"

printf "First command is:\n"
echo "Rscript ${scriptDir}/prs_step1_prepare_sumstats_for_MegaPRS.R"
printf "\n\n"
printf "==========================================\n\n"


Rscript ${scriptDir}/prs_step1_prepare_sumstats_for_MegaPRS.R


printf "\n\n"
printf "Sumstats prep complete.\n\n"
date
printf "\n\n"
printf "==========================================\n\n"
printf "Script complete. Goodbye.\n\n"
date
printf "==========================================\n\n"

