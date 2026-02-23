#!/bin/bash

#SBATCH --job-name=mt07_gzip_mr_mega
#SBATCH --time=3-00:00:00
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=128G
#SBATCH --account=project_2007428
#SBATCH --output=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/mt07_gzip_mr_mega.log
#SBATCH --error=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/mt07_gzip_mr_mega.log

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
megaDir=${mainDir}"outputs/mr_mega/"

printf "Script starts:\n\n"
date
printf "\n\n==========================\n\n"

printf "Gzipping all MR-MEGA results files."
printf "\n\n"
printf "==========================================\n\n"


for nfile in ${megaDir}*.result
    do
        echo "Currently gzipping ${nfile}."
        printf "\n\n----\n\n"
        gzip $nfile
    done

printf "\n\n"
printf "All files zipped.\n\n"
printf "==========================================\n\n"
printf "Script complete. Goodbye.\n\n"
date
printf "==========================================\n\n"

