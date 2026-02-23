#!/bin/bash

#SBATCH --job-name=mt03_01_ukb_liftover
#SBATCH --time=3-00:00:00
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=200G
#SBATCH --account=project_2007428
#SBATCH --output=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/mt03_01_ukb_liftover.log
#SBATCH --error=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/mt03_01_ukb_liftover.log

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

### Directories and files

mainDir="/scratch/project_2007428/projects/prj_001_cost_gwas/"
outDir="/scratch/project_2007428/projects/prj_001_cost_gwas/processing/misc_data/"

##### =========================== #####

### Start script

##### =========================== #####

printf "Script starts:\n\n"
date
printf "\n\n==========================\n\n"


printf "Command is:\n"
echo "Rscript ${mainDir}scripts/meta_step3_1_make_UKB_liftover_input.R"
printf "\n\n"
printf "==========================================\n\n"

Rscript ${mainDir}scripts/meta_step3_1_make_UKB_liftover_input.R


##### =========================== #####

### Running Liftover

##### =========================== #####

chain_file="/projappl/project_2007428/software/liftOver/chain_files/hg19ToHg38.over.chain.gz"
workingfile="/scratch/project_2007428/projects/prj_001_cost_gwas/processing/misc_data/ukb_ALLSNPs_hg19_liftover_input.txt"

/projappl/project_2007428/software/liftOver/liftOver ${workingfile} ${chain_file} ${outDir}ukb_ALLSNPs_hg38_liftOver_output.out ${outDir}liftOver.err


printf "\n\n==========================================\n\n"
printf "Script complete. Goodbye.\n\n"
date
printf "==========================================\n\n"

