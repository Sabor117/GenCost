#!/bin/bash

#SBATCH --job-name=mt015_format_files
#SBATCH --time=3-00:00:00
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=128G
#SBATCH --account=project_2007428
#SBATCH --output=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/mt015_format_files.log
#SBATCH --error=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/mt015_format_files.log

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
echo "Rscript ${mainDir}scripts/meta_step15_prep_files_for_upload.R"
printf "\n\n"
printf "==========================================\n\n"

Rscript ${mainDir}scripts/meta_step15_prep_files_for_upload.R

### Run the tar command
cd "$mainDir/outputs/"

tar -czvf gwas_catalogue_upload_MAIN.tar.gz METAL_v4_upload/*_ALL_*.tsv.gz

printf "\nMain files tarred.\n...\n\n"

tar -czvf gwas_catalogue_upload_stratified.tar.gz $(find METAL_v4_upload -type f -name "*.tsv.gz" ! -name "*_ALL_*")

printf "\n\n"
printf "Files formatted.\n\n"
printf "\n\n"
printf "==========================================\n\n"
printf "Script complete. Goodbye.\n\n"
date
printf "==========================================\n\n"
