#!/bin/bash

#SBATCH --job-name=Extr_cost_UKB_fields
#SBATCH --time=24:00:00
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8G
#SBATCH --account=project_2007428
#SBATCH --output=/scratch/project_2007428/projects/prj_001_cost_gwas_pipeline/logs/a01_field_extract.log
#SBATCH --error=/scratch/project_2007428/projects/prj_001_cost_gwas_pipeline/logs/a01_field_extract.log

module load r-env

mainDir="/scratch/project_2007428/projects/prj_001_cost_gwas_pipeline/"
dataDir="/scratch/project_2007428/data/processing/ukbb_78537/phenotypes/"

printf "Script starts:\n\n"
date
printf "\n\n==========================\n\n"

printf "First command is:\n"
echo "Rscript ${mainDir}scripts/step1_1_extract_ukb_fields.R"
printf "\n\n"
printf "==========================================\n\n"


Rscript ${mainDir}scripts/step1_1_extract_ukb_fields.R


printf "\n\n"
printf "Field list created. Extracting phenotypes. command is:\n"
echo "${dataDir}ukbconv ${dataDir}ukb671404.enc_ukb csv -Icosting_fields_list.txt -Oukb671404_cost_phenos"
printf "\n\n"
date
printf "\n\n"
printf "==========================================\n\n"


cd ${mainDir}/processing/

${dataDir}ukbconv ${dataDir}ukb671404.enc_ukb csv -Icosting_fields_list.txt -Oukb671404_cost_phenos


printf "\n==========================================\n\n"
printf "Script complete. Goodbye.\n\n"
date
printf "==========================================\n\n"

