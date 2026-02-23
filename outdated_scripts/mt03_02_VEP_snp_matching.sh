#!/bin/bash

#SBATCH --job-name=mt03_02_VEP_snp_match_%a
#SBATCH --time=3-00:00:00
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=32G
#SBATCH --account=project_2007428
#SBATCH --array=1-56
#SBATCH --output=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/mt03_02_VEP_snp_match_%a.log
#SBATCH --error=/scratch/project_2007428/projects/prj_001_cost_gwas/logs/mt03_02_VEP_snp_match_%a.log

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


##### =========================== #####

### Start script

##### =========================== #####

printf "Script starts:\n\n"
date
printf "\n\n==========================\n\n"

inputDir="/scratch/project_2007428/projects/prj_001_cost_gwas/processing/misc_data/VEP_inputs/"
outputDir="/scratch/project_2007428/projects/prj_001_cost_gwas/processing/misc_data/VEP_outputs/"
cacheDir="/scratch/project_2007428/data/processing/ensembl/vep_cache/"

run_index=$((SLURM_ARRAY_TASK_ID - 1))

file_list=("${inputDir}"gencost_allSNPs_for_VEP_*.txt)

file_input="${file_list[$run_index]}"

### Extracting the base file name

base_file_name=$(basename "$file_input")

### Extracting text between "VEP_" and ".txt"

extracted_name=${base_file_name#*VEP_}
extracted_name=${extracted_name%.txt}

file_output="${outputDir}gencost_allSNPs_VEP_out_${extracted_name}.txt"

printf "Command is:\n"
echo "/projappl/project_2007428/software/vep-111.0/bin/vep -i ${file_input} -o ${file_output} --cache --dir_cache ${cacheDir} --force_overwrite --verbose --check_existing"
printf "\n\n"
printf "==========================================\n\n"

/projappl/project_2007428/software/vep-111.0/bin/vep -i ${file_input} -o ${file_output} --cache --dir_cache ${cacheDir} --force_overwrite --verbose --check_existing

printf "\n\n==========================================\n\n"
printf "Script complete. Goodbye.\n\n"
date
printf "==========================================\n\n"
